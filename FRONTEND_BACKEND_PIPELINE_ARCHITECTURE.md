# Frontend ↔ Backend Pipeline Architecture
## Cross-Program Data Processing Roles & IPC Coordination

**Version:** 3.1  
**Last Updated:** Feb 24, 2026  
**Framework:** UQFF Star-Magic v3.1 (Self-Expanding Physics Backend)  
**Authors:** Daniel T. Murphy, Physics Architecture Team

---

## Executive Summary

Star-Magic uses a **3-tier architecture** separating concerns between user interface, computation orchestration, and backend physics services:

| Tier | Program | Role | Technology Stack |
|------|---------|------|-------------------|
| **Tier 1: GUI Orchestrator** | `source2.cpp` (15,753 lines) | User input → IPC dispatch → Result display | Qt6, VTK, QWebEngine |
| **Tier 2: Physics Backend** | `source2(HEAD PROGRAM).cpp` (4,021 lines) | IPC message receiver → GPU/symbolic compute → response sender | Qt3D, libtorch, SymEngine, pybind11 |
| **Tier 3: IPC Pipeline** | `ipc/` directory (15 files) | Message protocol, channel management, result caching, state synchronization | Named Pipes, SharedMemory, gRPC, protobuf |

**Key Innovation:** Separation of GUI orchestration (source2.cpp) from computation (source2(HEAD PROGRAM).cpp) via inter-process communication, enabling:
- **Parallel computation** - Backend runs simultaneously with UI
- **Scalability** - Multiple backends can serve single GUI
- **VR readiness** - Backend can run in headless VR mode with Qt3D rendering
- **Network access** - Optional gRPC layer for distributed computation

---

## 1. Individual Data Processing Roles

### 1.1 Tier 1: Frontend GUI Orchestrator (source2.cpp)

**File:** `source2.cpp` (15,753 lines)  
**Role:** User-facing interface, computation request dispatcher, result aggregator and displayer

#### Responsibilities:
1. **User Input Collection**
   - Tab 1: Query field (e.g., "Sagittarius A*", "M87", "NGC 3596")
   - Dropdowns for parameter selection (mass, radius, magnetic field, spin period)
   - Time range selectors for simulations

2. **APIFetch Integration**
   - Calls `APIFetch.py` internally
   - Queries 55 APIs (SIMBAD, NASA, VizieR, NED, Gaia, Grok fallback)
   - Returns `bodies_YYYYMMDD_HHMMSS.csv` with observational data
   - Dataset injected into physics calculation pipeline

3. **IPC Message Creation**
   - Converts user query → `CALCULATE_FIELD` or `CALCULATE_GRAVITY` IPC message
   - Populates `CalculateFieldRequest` payload:
     ```cpp
     struct CalculateFieldRequest {
         char system_name[64];    // "SGR 1935+2154"
         double r;               // 1e6 m (radial distance)
         double t;               // Current time (s)
         double tn;              // Negative time factor
         double theta;           // Angular position (rad)
         uint32_t flags;         // Calculation options
     };
     ```
   - Sets message type: `MessageType::CALCULATE_FIELD` (0x0001)
   - Sequence number for ordering
   - Timestamp for performance tracking

4. **Backend Dispatch**
   - Sends IPC message via `physics_service` channel (Named Pipes or SharedMemory)
   - Blocks waiting for `RESPONSE_DATA` message with results
   - Timeout handling for long-running computations

5. **Result Reception & Aggregation**
   - Receives `FieldResponse` payload from backend:
     ```cpp
     struct FieldResponse {
         double F_U;             // Total unified field
         double Ug1, Ug2, Ug3, Ug4;  // Gravity components
         double Um;              // Magnetism
         double Ubi;             // Buoyancy opposition
         double g_compressed;    // MUGE compressed gravity
         double residual;        // |calc - observed|
         double confidence;      // 0-1 confidence
     };
     ```
   - Aggregates multiple field calculations if simulating (parallel requests)
   - Caches results via `result_cache.py` for recall

6. **Result Display & Recall**
   - Tab 2-8: Render results in appropriate format
   - Tab 8 (VTK 3D Simulator): Visualize 3D field profiles
   - Tab 9 (Session Logger): Store results for recall via `CondensedPhysics_OutputData.py`
   - Tab 10+ (Domain-specific): Show equations, plots, predictions

#### Data Flow Example (User Query Path):
```
User Input (Tab 1)
    ↓ "SGR 1935+2154"
APIFetch.py ← calls 55 APIs
    ↓ 
bodies_YYYYMMDD_HHMMSS.csv ← {mass: 1.4e30 kg, radius: 10 km, B: 1e15 T}
    ↓
source2.cpp creates IPC message:
    - type: CALCULATE_FIELD (0x0001)
    - system_name: "sgr1935"
    - r: 1e6, t: 0, theta: 0
    ↓
Named Pipe Channel → physics_service
    ↓ [BACKEND PROCESSES]
    ↓
FieldResponse → Named Pipe Channel
    ↓
source2.cpp receives:
    - F_U: 2.341e-8 N/kg
    - Ug1: 1.023e-8 N/kg
    - confidence: 0.97
    ↓
Tab 8 (VTK): Render 3D gravity field
Tab 9 (Logger): Cache {query, response, timestamp}
```

---

### 1.2 Tier 2: Physics Backend (source2(HEAD PROGRAM).cpp)

**File:** `source2(HEAD PROGRAM).cpp` (4,021 lines)  
**Role:** Headless physics computation engine, dynamic term manager, simulation orchestrator

#### Responsibilities:

1. **IPC Server Initialization** (startup)
   - Reads `ServiceConfig` (from source2.cpp or command line)
   - Creates `physics_service` instance
   - Opens IPC channels (Named Pipes, SharedMemory, optional gRPC)
   - Spawns worker thread pool (default 4 threads)
   - Ready to accept messages

2. **Message Reception & Deserialization**
   ```cpp
   // IPC server loop (lines ~1500-1800)
   MessageHeader header;
   std::vector<uint8_t> payload;
   while (channel->receive(header, payload, timeout_ms)) {
       // Deserialize message type and payload
       switch(header.type) {
           case MessageType::CALCULATE_FIELD:
               handleCalculateField(*(CalculateFieldRequest*)payload.data());
               break;
           case MessageType::REGISTER_TERM:
               handleRegisterTerm(*(RegisterTermRequest*)payload.data());
               break;
           // ... etc
       }
   }
   ```

3. **Field Calculation** (CALCULATE_FIELD message)
   
   **Path 1: Native GPU Compute (libtorch)**
   ```cpp
   // NativeGPU namespace (lines ~400-500)
   namespace NativeGPU {
       // For FFT/matrix-heavy operations:
       auto tensor = torch::from_blob(...).to(getDevice());  // "cuda:0" or "cpu"
       auto fft_result = torch::fft::fft(tensor);
       // Result: 1-10ms for GPU, 100+ms for CPU
   }
   ```
   
   **Path 2: Symbolic Math (SymEngine)**
   ```cpp
   // NativeSymbolic namespace (lines ~500-550)
   namespace NativeSymbolic {
       // For equation solving:
       auto expr = SymEngine::parse("x^2 + 3*x");
       auto deriv = SymEngine::diff(expr, x);  // dF_U/dr
       double result = SymEngine::eval_double(*expr->subs(subs));
   }
   ```

   **Path 3: Integrated MAIN_1_CoAnQi Physics Engine**
   ```cpp
   // If enable_main1_integration = true:
   // Calls compiled MAIN_1_CoAnQi calculator directly
   // Computes:
   // - F_U = SOURCE4::compute_FU_SOURCE4(...)  // Complete unified field
   // - Ug1-4 = gravity components
   // - Um = magnetism
   // - Ubi = buoyancy opposition
   // - g_compressed = MUGE compressed gravity
   ```

   **Path 4: Python Subprocess (Advanced Features)**
   ```cpp
   // PythonBridge namespace (lines ~200-300)
   QProcess advancedFeatures;
   advancedFeatures.start("python", {"advanced_features.py"});
   // Sends JSON request, receives JSON result
   // Used for: AI agents, ONNX export, plugin execution
   ```

4. **Self-Expand: Dynamic Term Registration** (REGISTER_TERM message)
   ```cpp
   // physics_service.h lines ~70-130
   struct DynamicTerm {
       std::string name;           // e.g., "VacuumEnergy_Custom"
       DynamicTermType type;       // GRAVITY_MODIFIER, QUANTUM_COUPLING
       double coefficient;         // Multiplier
       double r_dependence;        // Power law exponent: term ∝ r^exponent
   };
   
   // Backend can accept new terms from VR, external models, or optimization
   // Terms added ADDITIVELY to core UQFF calculations:
   F_U_total = F_U_core + Σ(dynamic_terms[i].evaluate(r,t))
   ```

5. **Self-Update: Parameter Tuning** (UPDATE_PARAMETER message)
   ```cpp
   // physics_service.h lines ~140-180
   struct CalibratedParameters {
       double kappa;        // κ - calibrated to 0.0005/day
       double SSq;          // [SSq] - calibrated to 0.57
       double H_SCm;        // H_SCm - calibrated to 0.99
       double U_UA;         // U_UA - calibrated to 0.0001
       double beta_i;       // β_i - calibrated to 0.603
   };
   
   // At runtime, can adjust these values based on:
   // - VR user feedback
   // - Observational data fit
   // - Auto-optimization (learning_rate = 0.001)
   // Changes apply to next calculation immediately
   ```

6. **Self-Simulate: Time Evolution** (SIM_START/SIM_FRAME/SIM_COMPLETE messages)
   ```cpp
   // physics_service.h lines ~190-230
   struct SimulationConfig {
       double t_start, t_end;  // Time range [s]
       double dt;              // Time step (1 hour default)
       int frames;             // 1000 frames
       bool stream_to_vr;      // Stream results real-time
   };
   
   for (int step = 0; step < total_steps(); step += steps_per_frame()) {
       // Compute field at (r, t)
       FieldResponse frame_result = computeField(r, t);
       
       // Package as SimulationFrame
       struct SimulationFrame {
           double F_U, Ug1-4, Um;  // Field values
           double dF_dt;           // Time derivative
           double orbital_velocity;// v = sqrt(F_U*r)
       };
       
       // Send via IPC (SIM_FRAME message)
       channel->send(MessageHeader(SIM_FRAME), &frame_result);
       
       t += dt;  // Next step
   }
   
   // When complete:
   channel->send(MessageHeader(SIM_COMPLETE), nullptr);
   ```

7. **Pipeline Integration** (PIPELINE_PROCESS message - Feb 24, 2026)
   ```cpp
   // NEW: CondensedPhysics.py integration via IPC
   struct PipelineProcessRequest {
       char object_name[128];      // "SGR 1935+2154"
       uint32_t timeout_ms;        // Computation timeout
       char callback_id[64];       // For async callbacks
   };
   
   // Backend dispatches to CondensedPhysics calculation pipeline:
   // 1. Fetch observational data (APIFetch.py)
   // 2. Run ALL applicable calculators (176 available)
   // 3. Store results in OutputDataStore
   // 4. Return summary via PIPELINE_RESULT
   
   struct PipelineResultPayload {
       double F_U, Ug1-4, Um;      // Field results
       uint32_t calculators_run;   // e.g., 47 calculators executed
       uint32_t calculation_success;// e.g., 45 successful
       double compute_time_ms;     // Total time
       uint32_t json_payload_follows;  // Full result JSON attached
   };
   ```

8. **Response Serialization**
   - Packs `FieldResponse` payload
   - Sets message type: `MessageType::RESPONSE_DATA` (0x1002)
   - Sets sequence number matching original request
   - Sends via IPC channel back to frontendistrator

9. **Logging & Diagnostics**
   - Writes to `physics_service.log`
   - Tracks compute time per message type
   - Records dynamic term registrations
   - Monitors parameter updates for scientific audit trail

#### Data Flow Example (Backend Processing):
```
IPC Message Received (Named Pipe)
    ↓
Deserialize MessageHeader + CalculateFieldRequest payload
    type: CALCULATE_FIELD (0x0001)
    system_name: "sgr1935"
    r: 1e6 m, t: 0 s, theta: 0 rad
    ↓
Physics Computation (Parallel):
    ├─ GPU Path: libtorch FFT (1ms)
    ├─ Symbolic Path: SymEngine derivative (2ms)
    ├─ MAIN_1 Path: Integrated F_U + Ug1-4 (5ms)
    └─ Python Path: Advanced features (if requested)
    ↓
Apply dynamic terms (if registered):
    F_U_final = F_U_core + Σ(registered_terms[i].evaluate(r,t))
    ↓
Package FieldResponse:
    F_U: 2.341e-8, Ug1: 1.023e-8, confidence: 0.97
    ↓
Serialize response + FieldResponse payload
    type: RESPONSE_DATA (0x1002)
    sequence: (same as request)
    ↓
Send via IPC (Named Pipe)
```

---

### 1.3 Tier 3: IPC Pipeline (ipc/ directory)

**Files:** 15 files across 5 layers  
**Role:** Message protocol, channel management, caching, synchronization

#### Layer 1: Protocol Definition (uqff_ipc.h - 644 lines)

**MessageType Enum** (45 distinct message types):
```cpp
// Field calculations
CALCULATE_FIELD         = 0x0001
CALCULATE_GRAVITY       = 0x0002
CALCULATE_BUOYANCY      = 0x0003
CALCULATE_MAGNETISM     = 0x0004

// Dynamic term management
REGISTER_TERM           = 0x0010
UPDATE_PARAMETER        = 0x0012

// Simulation control
SIM_START               = 0x0200
SIM_FRAME               = 0x0210
SIM_COMPLETE            = 0x0212

// Pipeline (new Feb 24, 2026)
PIPELINE_PROCESS        = 0x0300
PIPELINE_RESULT         = 0x0301
PIPELINE_CALLBACK       = 0x0310

// Responses
RESPONSE_SUCCESS        = 0x1000
RESPONSE_ERROR          = 0x1001
RESPONSE_DATA           = 0x1002
```

**MessageHeader** (32 bytes, aligned):
```cpp
struct MessageHeader {  // 32 bytes total
    uint32_t magic;      // 0x55514646 = "UQFF"
    uint32_t version;    // Protocol version (1)
    MessageType type;    // Message type (above enum)
    uint32_t payload_size;  // Size of payload in bytes
    uint64_t timestamp;   // Us since epoch
    uint32_t sequence;    // Ordering number
    uint32_t flags;       // Reserved
};
```

**Payload Structures** (variable-length):
```cpp
// CALCULATE_FIELD
struct CalculateFieldRequest {
    char system_name[64];
    double r, t, tn, theta;
    uint32_t flags;
};

// REGISTER_TERM
struct RegisterTermRequest {
    char term_name[64];
    char description[256];
    double initial_value;
    uint32_t term_type;
};

// VR_FRAME_UPDATE
struct VRFrameUpdate {
    uint64_t frame_number;
    double head_position[3];
    double head_orientation[4];      // Quaternion
    double left_hand[7];             // Pos + quat
    double right_hand[7];
    uint32_t gesture_flags;
};

// PIPELINE_PROCESS (NEW)
struct PipelineProcessRequest {
    char object_name[128];           // e.g., "SGR 1935+2154"
    uint32_t timeout_ms;
    char callback_id[64];            // For async
};

// PIPELINE_RESULT (NEW)
struct PipelineResultPayload {
    char query_id[64];
    char object_type[32];
    uint32_t calculators_run;
    double F_U, Ug1-4, Um;           // Field results
    double g_compressed;
    double g_resonant;
    uint32_t status;                 // 0=success
    uint32_t json_payload_follows;   // Full JSON attached if 1
};
```

#### Layer 2: Channel Abstraction (uqff_ipc.h cont'd)

**IChannel Interface** (abstract):
```cpp
class IChannel {
    virtual bool send(const MessageHeader&, const void* payload) = 0;
    virtual bool receive(MessageHeader&, std::vector<uint8_t>&, int timeout_ms) = 0;
    virtual bool is_connected() const = 0;
};
```

**Implementations:**

1. **NamedPipeChannel** (Windows)
   - Uses `CreateNamedPipe()` / `CreateFile()`
   - Named pipe: `\\.\pipe\StarMagic_UQFF`
   - Advantages: Built-in Windows, synchronous
   - Use case: Single-machine GUI ↔ backend

2. **SharedMemoryChannel** (Cross-platform)
   - Uses `CreateFileMapping()` or `mmap()`
   - Shared buffer: 16 MB default (configurable)
   - Advantages: Zero-copy, fast
   - Use case: Parallel compute tasks

3. **GrpcChannel** (Optional, Network)
   - Full gRPC implementation with protobuf
   - Protocol: `grpc://localhost:50051`
   - Advantages: Network-transparent, structured
   - Use case: Distributed backends, cloud integration

#### Layer 3: Backend Service (physics_service.h - 508 lines)

**ServiceConfig**:
```cpp
struct ServiceConfig {
    std::string ipc_channel_name;    // "uqff_physics"
    size_t shm_buffer_size;          // 16 MB
    std::string grpc_address;        // "localhost:50051"
    bool enable_grpc;                // true in Phase 3
    int worker_threads;              // 4 default
    bool enable_main1_integration;   // Use MAIN_1_CoAnQi
    bool enable_python_integration;  // Use CondensedPhysics
    std::string log_file;            // "physics_service.log"
};
```

**Service Lifecycle**:
```cpp
// Startup
PhysicsService service(config);
service.start();                     // Spawn worker threads, open IPC

// Main loop
while (running) {
    MessageHeader header;
    std::vector<uint8_t> payload;
    
    if (channel.receive(header, payload, timeout_ms)) {
        // Dispatch based on message type
        switch(header.type) {
            case CALCULATE_FIELD:
                auto request = *(CalculateFieldRequest*)payload.data();
                auto response = computeField(request);
                channel.send(MessageHeader(RESPONSE_DATA, sizeof(response)), &response);
                break;
            // ... etc for 44 other message types
        }
    }
}

// Shutdown
service.stop();                      // Graceful shutdown, finish pending
```

#### Layer 4: Python Integration (python_bridge.h/cpp)

Enables Python subprocess coordination:
```cpp
namespace PythonBridge {
    // Launch advanced_features.py subprocess
    QProcess process;
    process.start("python", {"advanced_features.py"});
    
    // Send JSON request
    QJsonObject request;
    request["method"] = "run_ai_agent";
    request["query"] = "Analyze gravity field anomaly";
    
    // Receive JSON response
    QJsonObject response = callAdvancedFeatures("run_ai_agent", request);
    
    // Available methods:
    // - run_ai_agent (Grok 4 AI analysis)
    // - export_to_onnx (Model export)
    // - list_plugins (Available transformations)
    // - run_plugin (Custom computation)
}
```

#### Layer 5: Utility Services (state_sync.py, result_cache.py, etc.)

**state_sync.py**: Synchronize physics engine state across processes
```python
# When UPDATE_PARAMETER received in backend:
# Export state from physics_service
state = service.export_state()  # {kappa: 0.0005, SSq: 0.57, ...}

# Sync to frontend for display
state_sync.push_state(process_id, state)

# Frontend queries if needed
current_params = state_sync.get_state(backend_process_id)
```

**result_cache.py**: Persistent result caching
```python
# Backend computes result
result = {
    'F_U': 2.341e-8,
    'Ug1': 1.023e-8,
    'confidence': 0.97
}

# Cache it (SQLite or Redis)
cache.store(
    key=f"sgr1935_r=1e6_t=0",
    value=result,
    ttl=86400  # 1 day
)

# Frontend can recall
cached = cache.get("sgr1935_r=1e6_t=0")  # instant return
```

**file_watcher.py**: File-based IPC fallback
```python
# If Named Pipes unavailable:
# Frontend writes request to file
with open("ipc_request.json", "w") as f:
    json.dump(CalculateFieldRequest(...), f)

# Backend monitors directory, reads when available
watch_dir("ipc/", 
    on_create=lambda f: handle_request_file(f))

# Backend writes response
with open("ipc_response.json", "w") as f:
    json.dump(FieldResponse(...), f)

# Frontend polls for response
while not os.path.exists("ipc_response.json"):
    time.sleep(0.1)
```

**file_lock.py**: Concurrency control
```python
# Prevent reads while writes in progress
lock = FileLock("ipc_request.json.lock")

with lock.acquire():
    # Safe file operations
    pass
```

---

## 2. Complete Data Flow Diagrams

### 2.1 User Query → Physics Calculation → Display

```
SOURCE2.CPP (GUI Orchestrator)
┌──────────────────────────────────────────────┐
│  User enters: "SGR 1935+2154" in Tab 1       │
│  Clicks "Analyze"                            │
└────────────────────┬────────────────────────┘
                     ↓
┌──────────────────────────────────────────────────────────────────────┐
│  SOURCE2.CPP: Create IPC message                                    │
│  - type: CALCULATE_FIELD (0x0001)                                   │
│  - system_name: "sgr1935"                                           │
│  - r: 1e6 m (distance)                                              │
│  - t: 0 s (initial time)                                            │
│  - tn: 1.0 (negative time factor)                                   │
│  - theta: 0 rad (angle)                                             │
│  - flags: 0x0001 (include magnetism)                                │
└────────────────────┬────────────────────────────────────────────────┘
                     ↓
                [Named Pipe IPC]
                     ↓
┌──────────────────────────────────────────────────────────────────────┐
│  SOURCE2(HEAD PROGRAM).CPP (Physics Backend)                        │
│  Receives IPC message                                               │
└────────────────────┬────────────────────────────────────────────────┘
                     ↓
            [DISPATCH to computation paths]
            ┌───────┬────────┬──────────┐
            ↓       ↓        ↓          ↓
         ┌────┐ ┌────┐  ┌──────┐    ┌────────┐
         │GPU │ │Sym │  │MAIN_1│    │Python  │
         │    │ │Eng │  │CoAnQi│    │Proc    │
         └─┬──┘ └─┬──┘  └───┬──┘    └────┬───┘
           │ 1ms  │ 2ms     │ 5ms        │ 100ms
           └──────┴─────────┴──────────┬─┴────────
                            [MERGE RESULTS]
                                  ↓
         ┌─────────────────────────────────────┐
         │ F_U: 2.341e-8 N/kg                  │
         │ Ug1: 1.023e-8, Ug2: 0.512e-8, ...  │
         │ Um: 0.341e-8 N/kg                   │
         │ Ubi: 1.023e-9 N/kg (buoyancy)       │
         │ g_compressed: 1.512e-8 m/s²         │
         │ confidence: 0.97 (based on data fit) │
         └──────────┬──────────────────────────┘
                    ↓
                [Serialize to FieldResponse]
                    ↓
               [Named Pipe IPC]
                    ↓
         ┌──────────────────────────┐
         │ MESSAGE: RESPONSE_DATA    │
         │ payload_size: 64 bytes    │
         │ sequence: (match request) │
         └──────────┬────────────────┘
                    ↓
              SOURCE2.CPP
         ┌────────────────────────┐
         │ Receives response       │
         │ Unpacks FieldResponse   │
         └────────┬───────────────┘
                  ↓
         ┌──────────────────────────────────────┐
         │ DISPLAY RESULTS:                     │
         │ Tab 3: Show F_U value (2.341e-8)     │
         │ Tab 4: Show field components         │
         │ Tab 8: Render 3D gravity field       │
         │ Tab 9: Cache query + result          │
         │        for later recall              │
         └──────────────────────────────────────┘
```

### 2.2 Simulation Stream (Real-time VR Visualization)

```
SOURCE2.CPP
┌──────────────────────────────┐
│ User clicks "Simulate"       │
│ Config: t=0 to 1 year,       │
│         dt=1 hour,           │
│         1000 frames          │
└────────────┬─────────────────┘
             ↓
    [Create SIM_START message]
             ↓
        [Named Pipe IPC]
             ↓
SOURCE2(HEAD PROGRAM).CPP
┌──────────────────────────────────────────────────┐
│ Receives SIM_START message                       │
│ Allocates simulation state                       │
│ Starts time evolution loop:                      │
│                                                  │
│ for t = 0 to 1_year step 1_hour:                │
│     for r = r_start to r_end step dr:           │
│         compute F_U(r, t)                        │
│         dF/dt ≈ (F(t+dt) - F(t)) / dt           │
│         v_orbit = sqrt(F_U * r)                  │
│                                                  │
│         ← Pack as SimulationFrame                │
│         ← Send via IPC (SIM_FRAME msg)           │
│                                                  │
│     if (t % steps_per_frame == 0):              │
│         yield control (frame complete)          │
│                                                  │
│ When done:                                       │
│     ← Send SIM_COMPLETE message                  │
└──────────────────────────────────────────────────┘
             ↑
         [IPC STREAM - 1000 messages]
             ↓
SOURCE2.CPP (VR Runtime)
┌──────────────────────────────────────────────────┐
│ Real-time visualization as frames arrive:        │
│                                                  │
│ foreach SIM_FRAME message received:              │
│     Update 3D visualization (Tab 8)              │
│     Show spiral field lines ← F_U(r,t)          │
│     Show orbital paths ← v_orbit(r,t)           │
│     Show time slider ← current t                │
│     Refresh 30 FPS (or IPC message rate)        │
│                                                  │
│ When SIM_COMPLETE:                              │
│     Animation finished                          │
│     Full trajectory loaded                      │
│     User can replay, zoom, inspect              │
└──────────────────────────────────────────────────┘
```

### 2.3 Dynamic Term Registration (Self-Expand)

```
SOURCE2.CPP (VR User or AI Agent)
┌──────────────────────────────────────────────────┐
│ User (or AI in advanced_features.py):            │
│ "Add vacuum fluctuation term with coeff 1e-10"  │
└──────────────┬─────────────────────────────────┘
               ↓
    [Create REGISTER_TERM message]
    ┌────────────────────────────────────────┐
    │ term_name: "VacuumFluctuation_Custom"  │
    │ type: VACUUM_ENERGY (0x0002)           │
    │ description: "User-defined Term"       │
    │ initial_value: 1.0                     │
    │ coefficient: 1e-10                     │
    │ r_dependence: -3.0  (∝ 1/r³)           │
    │ t_dependence: 0.0   (time-independent)│
    └────────────────────────────────────────┘
               ↓
          [Named Pipe IPC]
               ↓
SOURCE2(HEAD PROGRAM).CPP
┌──────────────────────────────────────────────────┐
│ physics_service::handleRegisterTerm()            │
│                                                  │
│ Add to DynamicTermRegistry:                      │
│ registry[new_id] = DynamicTerm {                 │
│     name: "VacuumFluctuation_Custom",            │
│     type: VACUUM_ENERGY,                         │
│     coefficient: 1e-10,                          │
│     r_dependence: -3.0,                          │
│     enabled: true,                               │
│     source: "VR"                                 │
│ }                                                │
│                                                  │
│ From this point on, all F_U calculations include:│
│ F_U_total = F_U_core +                           │
│             Σ(registry[i].evaluate(r,t))         │
│                                                  │
│ where evaluate() = coeff * base * r^(-3) * t^0  │
│                  = 1e-10 * 1.0 * (r^-3)          │
└──────────────────────────────────────────────────┘
               ↓
        [RESPONSE_SUCCESS]
               ↓
SOURCE2.CPP
┌──────────────────────────────────────────────────┐
│ Receive confirmation                             │
│ Display in TUI: "New term registered"            │
│ Next calculations automatically use new term     │
└──────────────────────────────────────────────────┘
```

### 2.4 Pipeline Integration (CondensedPhysics - Feb 24, 2026)

```
SOURCE2.CPP
┌─────────────────────────────────────────────┐
│ User clicks "Full Analysis"                 │
│ Input: "SGR 1935+2154"                      │
└────────────┬────────────────────────────────┘
             ↓
    [Create PIPELINE_PROCESS message]
    ┌────────────────────────────────────────┐
    │ type: PIPELINE_PROCESS (0x0300)         │
    │ object_name: "SGR 1935+2154"            │
    │ timeout_ms: 30000 (30 seconds)          │
    │ callback_id: "analysis_01"              │
    └────────────────────────────────────────┘
             ↓
        [Named Pipe IPC]
             ↓
SOURCE2(HEAD PROGRAM).CPP
┌────────────────────────────────────────────────────────┐
│ physics_service::handlePipelineProcess()              │
│                                                        │
│ ① FETCH PHASE (APIFetch.py)                           │
│    └─ Query 55 APIs for observational data            │
│       → bodies_YYYYMMDD_HHMMSS.csv                    │
│       → {mass, radius, distance, magnetic_field, ...} │
│                                                        │
│ ② COMPUTE PHASE (CondensedPhysics.py)                 │
│    └─ Run 176 available calculators in parallel       │
│       Calculator 1: TriadicGravityCalculator          │
│       Calculator 2: MUGECompressed                     │
│       ... (176 total)                                 │
│       └─ Each computes physics equation for this obj  │
│                                                        │
│ ③ STORE PHASE (OutputDataStore)                       │
│    └─ Save all equation results to:                   │
│       CondensedPhysics_OutputData.py                  │
│       {query_id: "analysis_01",                       │
│        equations: [eq1_result, eq2_result, ...],      │
│        timestamp: ...,                                │
│        calculators_success: 174/176}                  │
│                                                        │
│ ④ SUMMARY PHASE (Response generation)                 │
│    └─ Extract key results for response:               │
│       - F_U (from Calculator 5)                       │
│       - g_compressed (from Calculator 87)             │
│       - confidence metrics                            │
│       - compute time: 5.2 seconds                      │
│       - calculators_run: 176, success: 174            │
└────────────────────────────────────────────────────────┘
             ↓
    [Create PIPELINE_RESULT message]
    ┌────────────────────────────────────────────────┐
    │ type: PIPELINE_RESULT (0x0301)                  │
    │ query_id: "analysis_01"                         │
    │ object_type: "Magnetar"                         │
    │ calculators_run: 176                            │
    │ calculation_success: 174                         │
    │ F_U: 2.341e-8                                   │
    │ Ug1-4: {...}                                    │
    │ compute_time_ms: 5200                           │
    │ json_payload_follows: 1                         │
    │ [FULL JSON RESULT ATTACHED]                     │
    └────────────────────────────────────────────────┘
             ↓
        [Named Pipe IPC]
             ↓
SOURCE2.CPP
┌────────────────────────────────────────────┐
│ Receive PIPELINE_RESULT                    │
│ Parse payload:                              │
│  - Display "SGR 1935+2154: Magnetar"        │
│  - Show F_U = 2.341e-8 N/kg                │
│  - Show "174/176 calculators succeeded"     │
│  - Show compute time: 5.2 seconds           │
│                                              │
│ Store full JSON to OutputStore:             │
│  CondensedPhysics_OutputData.py             │
│                                              │
│ User can later:                             │
│  - Recall full equation set (Tab 9 Logger) │
│  - Export JSON/CSV for publication          │
│  - Inspect individual calculator results    │
└────────────────────────────────────────────┘
```

---

## 3. Key Architectural Principles

### 3.1 Separation of Concerns

| Concern | Owner | Why Separated |
|---------|-------|---------------|
| User Interaction | source2.cpp | Responsive GUI, no blocking I/O |
| Physics Computation | source2(HEAD).cpp | CPU/GPU intensive, can timeout, run headless |
| Message Protocol | ipc/ layer | Enable multiple transports (pipes, memory, gRPC) |
| Python Integration | PythonBridge, advanced_features.py | Encapsulate subprocess, JSON IPC |
| Data Persistence | CondensedPhysics_OutputData.py | Decouple compute from storage |

### 3.2 Parallelism Model

**Three Concurrency Levels:**

1. **Frontend-Backend Parallelism** (Asynchronous IPC)
   - source2.cpp sends request, continues UI responsiveness
   - source2(HEAD).cpp computes in background
   - Results streamed back as available (SIM_FRAME messages)
   - No blocking

2. **Computation Parallelism** (Multi-threaded Physics)
   - source2(HEAD).cpp worker thread pool (default 4 threads)
   - Multiple field calculations in parallel
   - One thread per IPC message (if sufficient threads)
   - GPU operations on CUDA cores

3. **Pipeline Parallelism** (Parallel Calculators in CondensedPhysics)
   - 176 calculators run simultaneously for one query
   - Each calculator solves one-to-many equations
   - Results aggregated in OutputDataStore
   - Total pipeline time: 5-30 seconds (vs. sequential 5+ minutes)

### 3.3 Extensibility Points

**Where new code integrates without modifying core:**

1. **New Calculators in CondensedPhysics.py**
   - Add class inheriting from `PhysicsCalculator`
   - Automatically discovered and executed in pipeline
   - Results stored in OutputDataStore

2. **New Message Types in IPC**
   - Add enum to `MessageType` in uqff_ipc.h
   - Add payload struct for new data
   - Add handler in `physics_service::dispatch()`
   - No change to message header (extensible design)

3. **Dynamic Terms via REGISTER_TERM**
   - User/VR can add custom physics terms at runtime
   - Applied additively to F_U calculation
   - No recompilation needed

4. **Parameter Tuning via UPDATE_PARAMETER**
   - Adjust κ, [SSq], β_i, etc. at runtime
   - Changes apply immediately
   - Enables A/B testing, optimization loops

5. **Transport Layer via new IChannel**
   - Currently: NamedPipeChannel, SharedMemoryChannel, GrpcChannel
   - Can add: WebSocketChannel, UnixSocketChannel, MQChannel
   - Polymorphic dispatch in physics_service

---

## 4. Performance Characteristics

### 4.1 Latency Budget

```
User Query "SGR 1935+2154"
│
├─ GUI Paint: 16ms (60 FPS)
│
├─ API Fetch (APIFetch.py): 500-2000ms (depends on network)
│
├─ IPC Send: <1ms (Named Pipe)
│
├─ Backend Parsing: <1ms
│
├─ Physics Computation:
│  ├─ GPU Path (FFT): 1-10ms
│  ├─ Symbolic Path (SymEngine): 2-5ms
│  ├─ MAIN_1 Path: 5-50ms
│  └─ Python Subprocess: 50-500ms
│  └─ Total parallel: max(1ms, 2ms, 5ms, 50ms) = 50ms
│
├─ Response Serialization: <1ms
│
├─ IPC Send: <1ms (Named Pipe)
│
├─ GUI Update (VTK render): 16ms
│
└─ TOTAL: ~550-2000ms to display result
```

### 4.2 Throughput

| Operation | Rate | Limit |
|-----------|------|-------|
| GUI queries per second | 10-30 | UI responsiveness |
| Backend computations parallel | 4 (thread pool) | CPU cores |
| IPC messages per second | 1000+ | Named Pipe bandwidth |
| Simulation frames per second | 30-60 | VR refresh rate |
| Pipeline queries per second | 1 (sequential) | ~5 second compute time |

### 4.3 Memory Usage

| Component | Typical Size |
|-----------|-------------|
| source2.cpp (GUI) | 100-200 MB |
| source2(HEAD).cpp (Backend) | 200-500 MB |
| IPC Shared Buffer | 16 MB (configurable) |
| Result Cache | 100-1000 MB (depends on query count) |
| CondensedPhysics_OutputData.py | 50-500 MB per year of queries |
| **Total (typical)** | **~500 MB - 2 GB** |

---

## 5. Error Handling & Robustness

### 5.1 Frontend Error Paths

```cpp
// source2.cpp
try {
    // Send CALCULATE_FIELD message
    message.send(channel, timeout_ms=5000);
    
    // Wait for response with timeout
    if (!channel.receive(response, timeout_ms=5000)) {
        throw TimeoutException("Backend did not respond within 5s");
    }
    
    // Validate response
    if (!response.is_valid()) {
        throw ValidationException("Response magic number mismatch");
    }
    
    // Display results
    displayResults(response);
    
} catch (const TimeoutException& e) {
    // Backend crashed or hung
    UI.showError("Computation timed out. Check backend logs.");
    UI.suggestAction("Restart backend service");
    
} catch (const ValidationException& e) {
    // Corrupted message
    UI.showError("Data corruption detected over IPC");
    UI.suggestAction("Check physics_service.log for details");
}
```

### 5.2 Backend Error Paths

```cpp
// source2(HEAD PROGRAM).cpp
physics_service::dispatch(MessageHeader header, vector<uint8_t> payload) {
    try {
        switch(header.type) {
            case CALCULATE_FIELD: {
                auto req = *(CalculateFieldRequest*)payload.data();
                auto resp = computeField(req);
                
                if (!resp.success) {
                    // Return RESPONSE_ERROR instead
                    ResponseError err;
                    err.status = 1;
                    err.message = "Computation failed";
                    channel.send(MessageHeader(RESPONSE_ERROR), &err);
                } else {
                    channel.send(MessageHeader(RESPONSE_DATA), &resp);
                }
            } break;
            
            // ... other cases
        }
    } catch (const std::exception& e) {
        // Unexpected error - log and notify frontend
        log.error("Computation error: {}", e.what());
        
        ResponseError err;
        err.status = 2;
        err.message = e.what();
        channel.send(MessageHeader(RESPONSE_ERROR), &err);
    }
}
```

### 5.3 Graceful Degradation

| Failure Scenario | Behavior |
|------------------|----------|
| MAIN_1_CoAnQi not available | Use GPU path + SymEngine instead (slightly lower accuracy) |
| GPU not available | Fall back to CPU (10x slower but still functional) |
| SymEngine not compiled | Use numerical differentiation |
| Python subprocess timeout | Return partial result from GPU path |
| IPC NamedPipe unavailable | Fall back to SharedMemory or FileWatcher |
| Backend not running | Queue messages, retry with exponential backoff |
| Result cache unavailable | Recompute (performance cost, no correctness impact) |

---

## 6. Integration with CondensedPhysics.py

### Key Connection Points

**Phase 1: Request Entry**
```
source2.cpp
    ↓ (APIFetch.py)
bodies_YYYYMMDD_HHMMSS.csv
    ↓ IPData.py
CondensedPhysics_OutputData.py ← ① QUERY STORED
```

**Phase 2: Backend Dispatch**
```
source2(HEAD).cpp receives PIPELINE_PROCESS message
    ↓
physics_service dispatches to:
    → CondensedPhysics.py (176 calculators)
    → MAIN_1_CoAnQi.exe (Optional, for validation)
    → advanced_features.py (AI analysis)
    ↓
Results computed in parallel
```

**Phase 3: Result Storage**
```
All 176 calculator results
    ↓ ② RESULTS STORED
CondensedPhysics_OutputData.py
    ↓
source2.cpp Session Logger (Tab 9)
    ↓ User can RECALL
```

**Session Logger Workflow:**
1. User makes query (stored in OutputData as Step ①)
2. Pipeline runs (results stored as Step ②)
3. User later clicks "Recall" in Tab 9
4. Session Logger retrieves {query, results} from OutputData
5. Equations and plots regenerated without recomputation

---

## 7. Recommended Next Steps

1. **Verify IPC Connectivity**
   - Compile source2(HEAD PROGRAM).cpp with IPC support
   - Launch with `--service` flag
   - source2.cpp should auto-discover on startup

2. **Test Message Round-trip**
   - Add diagnostic logging to uqff_ipc.h
   - Trace: CALCULATE_FIELD request → backend → RESPONSE_DATA response
   - Measure latency for each pathway (GPU/SymEngine/MAIN_1/Python)

3. **Performance Profiling**
   - Use `physics_service.log` to measure compute times
   - Identify bottleneck: API fetch, GPU, SymEngine, or Python?
   - Optimize slowest path first

4. **Scale Testing**
   - Increase worker_threads in ServiceConfig (8, 16, 32)
   - Measure throughput (queries/second)
   - Identify thread pool saturation point

5. **Network Integration**
   - Enable gRPC channel in ServiceConfig
   - Run backend on remote machine
   - Measure: network latency + compute latency

---

## Appendix: Message Flow Reference

### All 45 Message Types Quick Reference

```cpp
// Field Calculations (4 types)
0x0001: CALCULATE_FIELD         ← Response: RESPONSE_DATA + FieldResponse
0x0002: CALCULATE_GRAVITY       ← Response: RESPONSE_DATA + FieldResponse
0x0003: CALCULATE_BUOYANCY      ← Response: RESPONSE_DATA
0x0004: CALCULATE_MAGNETISM     ← Response: RESPONSE_DATA

// Dynamic Term Management (4 types)
0x0010: REGISTER_TERM           ← Response: RESPONSE_SUCCESS
0x0011: UNREGISTER_TERM         ← Response: RESPONSE_SUCCESS
0x0012: UPDATE_PARAMETER        ← Response: RESPONSE_SUCCESS
0x0013: GET_PARAMETER           ← Response: RESPONSE_DATA

// System Management (3 types)
0x0020: GET_SYSTEM_LIST         ← Response: RESPONSE_DATA + list
0x0021: SELECT_SYSTEM           ← Response: RESPONSE_SUCCESS
0x0022: ADD_CUSTOM_SYSTEM       ← Response: RESPONSE_SUCCESS

// State Management (3 types)
0x0030: EXPORT_STATE            ← Response: RESPONSE_DATA + state JSON
0x0031: IMPORT_STATE            ← Response: RESPONSE_SUCCESS
0x0032: SYNC_STATE              ← Response: RESPONSE_SUCCESS

// VR Operations (3 types)
0x0100: VR_FRAME_UPDATE         ← Response: none (telemetry)
0x0101: VR_GESTURE_INPUT        ← Response: RESPONSE_SUCCESS
0x0102: VR_RENDER_REQUEST       ← Response: RESPONSE_DATA + mesh

// Simulation Control (5 types)
0x0200: SIM_START               ← Response: RESPONSE_SUCCESS
0x0201: SIM_STOP                ← Response: RESPONSE_SUCCESS
0x0202: SIM_PAUSE               ← Response: RESPONSE_SUCCESS
0x0203: SIM_RESUME              ← Response: RESPONSE_SUCCESS
0x0210: SIM_FRAME               ← Response: none (streaming output)

// Simulation Metadata (2 types)
0x0211: SIM_PROGRESS            ← Response: none (telemetry)
0x0212: SIM_COMPLETE            ← Response: none (final message)

// Pipeline Operations (3 types - Feb 24, 2026)
0x0300: PIPELINE_PROCESS        ← Response: RESPONSE_DATA + PIPELINE_RESULT
0x0301: PIPELINE_RESULT         ← Response: none (server sends)
0x0310: PIPELINE_CALLBACK       ← Response: none (event stream)

// Responses (3 types)
0x1000: RESPONSE_SUCCESS        ← Operation succeeded
0x1001: RESPONSE_ERROR          ← Operation failed
0x1002: RESPONSE_DATA           ← Data payload follows

// Control (3 types)
0xFF00: PING                    ← Response: PONG
0xFF01: PONG                    ← Response: none
0xFFFF: SHUTDOWN                ← Graceful termination

// TOTAL: 45 message types
```

---

**Document Generated:** February 24, 2026  
**Framework Version:** UQFF Star-Magic v3.1  
**For questions:** See commits 3e66d94 (SOURCE4), 79e73ec (SOURCE115), 59fd4c4 (SOURCE116)
