# UQFF Bridging Architecture

> **Version:** 2.1 CANONICAL (matches ARCHITECTURE_FLOW_DIAGRAM.md v4.1)
> **Updated:** February 21, 2026 (Phase 3 gRPC Complete)
> **CRITICAL:** source2.cpp = PRINCIPAL GUI (FIRST), source2(HEAD PROGRAM).cpp = VR/VM Backend

## IPC Pipeline Status

| Phase | Status | Description |
|-------|--------|-------------|
| Phase 1 | ✅ COMPLETE | IPC Pipeline Connection (SharedMemory + NamedPipe) |
| Phase 2 | ✅ COMPLETE | Physics Backend Service Mode (--service flag) |
| Phase 3 | ✅ COMPLETE | Full gRPC Implementation (port 50051) |
| Phase 4 | ✅ COMPLETE | Astro Graphics IPC Integration |
| Phase 5 | ✅ COMPLETE | Full VR Experience (Gesture → Physics → Render) |

## Complete Data Flow Diagram

This document describes the integration architecture between all UQFF Star-Magic components, including the FTPS bridging layer, JavaScript engine, C++ calculators, Python calculators, and the Qt GUI.

```
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                           EXTERNAL USER INTERFACE                                    │
│                                                                                      │
│   User Query (HTTPS) → Web Browser / External Application / CLI                     │
└─────────────────────────────────────────────────────────────────────────────────────┘
                                        │
                                        │ (HTTPS to external server)
                                        ▼
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                         EXTERNAL FTPS SERVER (Port 990/21)                          │
│                                                                                      │
│   Purpose: Secure file transfer endpoint for remote queries                         │
│   Ports:                                                                             │
│     - Implicit FTPS: Port 990 (TLS from connection start)                           │
│     - Explicit FTPS: Port 21 (upgrades via STARTTLS)                                │
│                                                                                      │
│   Operations:                                                                        │
│     1. Receives user query via HTTPS → converts to req_*.json                       │
│     2. Writes request file to /uqff_data/requests/                                  │
│     3. Polls /uqff_data/responses/ for result                                       │
│     4. Returns result to user via HTTPS                                             │
└─────────────────────────────────────────────────────────────────────────────────────┘
                                        │
                      (File-based RPC via req_*.json files)
                                        ▼
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                    ./uqff_data/requests/req_*.json                                  │
│                                                                                      │
│   Example request:                                                                   │
│   {                                                                                  │
│       "request_id": "20260221_050034_123456",                                       │
│       "timestamp": "2026-02-21T05:00:34.003Z",                                      │
│       "params": {                                                                    │
│           "system": "SgrA*",                                                         │
│           "M": 8.155e36,                                                             │
│           "r": 4.4e19,                                                               │
│           "B0": 1e-4                                                                 │
│       }                                                                              │
│   }                                                                                  │
└─────────────────────────────────────────────────────────────────────────────────────┘
                                        │
                         (Polled by uqff_server.js)
                                        ▼
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                       uqff_server.js (HTTP Port 3141)                               │
│                                                                                      │
│   Two interfaces:                                                                    │
│   ┌─────────────────────────────────┐  ┌─────────────────────────────────┐         │
│   │  1. HTTP REST API (Local)       │  │  2. File-based RPC (FTPS)       │         │
│   │  http://127.0.0.1:3141          │  │  Polls ./uqff_data/requests/    │         │
│   │                                  │  │  Writes to ./uqff_data/responses│         │
│   │  Endpoints:                      │  │                                  │         │
│   │  GET  /api/health               │  │  Enable via:                     │         │
│   │  GET  /api/constants            │  │  UQFF_FILE_RPC=true             │         │
│   │  GET  /api/systems              │  │  UQFF_REQUEST_DIR=./uqff_data/  │         │
│   │  POST /api/compute              │  │                      requests   │         │
│   │  POST /api/batch                │  │                                  │         │
│   └─────────────────────────────────┘  └─────────────────────────────────┘         │
│                                                                                      │
│   Imports: index.js LIBRARY INDEX (23,790 lines, 106 astrophysical systems)         │
│            Usage: const UQFF = require('./index.js'); UQFF.computeSagA(params);     │
└─────────────────────────────────────────────────────────────────────────────────────┘
                                        │
                                        ▼
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                         index.js (JavaScript UQFF LIBRARY INDEX)                     │
│                                                                                      │
│   Lines: 23,790+ | Systems: 106 astrophysical bodies                                │
│   ** NOTE: This is a LIBRARY INDEX, NOT a calculator **                             │
│   ** uqff_server.js imports this library to serve REST API **                       │
│                                                                                      │
│   Core Computations (exported for require()):                                       │
│   ├── Ug1: Magnetic dipole gravity                                                  │
│   ├── Ug2: Charge/reactivity gravity                                                │
│   ├── Ug3: String rotation gravity                                                  │
│   ├── Ug4: Vacuum concentration gravity                                             │
│   ├── F_U_Bi_i: Master buoyancy force                                               │
│   └── compressed_g: 26-layer gravity                                                │
│                                                                                      │
│   Constants (synchronized with shared_constants.h):                                 │
│   ├── RHO_VAC_UA: 7.09e-36 J/m³ (gravitational scale)                               │
│   ├── RHO_VAC_UA_FIELD: 1e-27 J/m³ (field scale)                                    │
│   ├── GRADIENT_RATIO: 7.09e-9 (DPM coupling)                                        │
│   └── ... 40+ physics constants                                                     │
└─────────────────────────────────────────────────────────────────────────────────────┘
                                        │
                                        ▼
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                    ./uqff_data/responses/resp_*.json                                │
│                                                                                      │
│   Example response:                                                                  │
│   {                                                                                  │
│       "request_id": "20260221_050034_123456",                                       │
│       "timestamp": "2026-02-21T05:00:34.103Z",                                      │
│       "status": "success",                                                           │
│       "result": {                                                                    │
│           "system_name": "Sagittarius A*",                                          │
│           "Ug1": 44000.0,                                                            │
│           "Ug2": 2.81e-23,                                                           │
│           "Ug3": 2.81e-19,                                                           │
│           "Ug4": 1.84e-22,                                                           │
│           "F_U_Bi_i": 2.16e+41,                                                      │
│           "g_compressed": 38547.76                                                   │
│       }                                                                              │
│   }                                                                                  │
└─────────────────────────────────────────────────────────────────────────────────────┘
                                        │
                        (Read by uqff_ftps_client.py)
                                        ▼
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                       uqff_ftps_client.py (Python FTPS Client)                      │
│                                                                                      │
│   Lines: 890+ | RFC 4217 Compliant | TLS 1.3 Enforced                               │
│                                                                                      │
│   Features:                                                                          │
│   ├── Explicit FTPS (port 21 + STARTTLS)                                            │
│   ├── Implicit FTPS (port 990)                                                      │
│   ├── Certificate verification                                                      │
│   ├── MD5 checksum validation                                                       │
│   ├── Auto-retry with backoff                                                       │
│   ├── Progress callbacks                                                            │
│   ├── Batch operations                                                              │
│   └── call_uqff_server() - File-based RPC bridge                                    │
└─────────────────────────────────────────────────────────────────────────────────────┘

═══════════════════════════════════════════════════════════════════════════════════════
           source2.cpp - PRINCIPAL GUI APPLICATION (USER STARTS HERE)
═══════════════════════════════════════════════════════════════════════════════════════

┌─────────────────────────────────────────────────────────────────────────────────────┐
│                   source2.cpp (PRINCIPAL GUI - USER STARTS HERE)                    │
│                                                                                      │
│   Lines: 11,058+ | Tabs: 21 | Qt6 + VTK + Chromium                                  │
│                                                                                      │
│   ** THIS IS THE ENTRY POINT - USER LAUNCHES THIS APPLICATION FIRST **              │
│                                                                                      │
│   ┌─────────────────────────────────────────────────────────────────────────────┐   │
│   │                         USER QUERY FIELD                                     │   │
│   │   "Sagittarius A*", "M87", "Betelgeuse", "NGC 3596"...                      │   │
│   └─────────────────────────────────────────────────────────────────────────────┘   │
│                                        │                                             │
│                                        ▼                                             │
│   ┌─────────────────────────────────────────────────────────────────────────────┐   │
│   │                    API FETCH LAYER (APIFetch.py)                             │   │
│   │   55 APIs: SIMBAD → NASA → VizieR → NED → Gaia → Grok fallback              │   │
│   │   Output: bodies_YYYYMMDD_HHMMSS.csv                                         │   │
│   └─────────────────────────────────────────────────────────────────────────────┘   │
│                                                                                      │
│   Tab Layout:                                                                        │
│   ├── Tab 1:  🎛️ MAIN_1 Calculator (PowerShellTerminalWidget)                      │
│   ├── Tab 2:  🐍 QCalc.py (PythonTerminalWidget)                                    │
│   ├── Tab 3:  🧮 Scientific Calculator                                              │
│   ├── Tab 4:  📓 Notebook Editor                                                    │
│   ├── Tab 5:  📚 CondensedPhysics.py                                                │
│   ├── Tab 6:  🤖 CoAnQi_bot (Ollama)                                                │
│   ├── Tab 7:  🧠 SuperGrok4 (xAI)                                                   │
│   ├── Tab 8:  🌌 UQFF Simulator (3D VTK)                                            │
│   ├── Tab 9:  📋 Session Logger (RECALL ACCESS)                                     │
│   ├── Tab 10: ⚖️ Compare C++/Python                                                 │
│   ├── Tab 11: 📐 Equation Display (IPC to QCalc.py)                                 │
│   ├── Tab 12: 🌐 JS Engine (HTTP to uqff_server.js)                                 │
│   └── Tab 13-21: Browser windows for search results                                │
│                                                                                      │
│   Includes:                                                                          │
│   ├── equation_renderer.h - Long-form equation display with IPC                    │
│   ├── source2_widgets_enhanced.h - UQFFJavaScriptWidget                            │
│   ├── source2_mainwindow.h - Qt MOC declarations                                   │
│   ├── source2_event_bus.h - Inter-widget events                                    │
│   └── shared_constants.h - UQFF:: namespace constants                              │
└─────────────────────────────────────────────────────────────────────────────────────┘
         │                    │                    │                    │
         │     SIMULTANEOUS JOINT OPERATION PIPELINE (parallel dispatch)             │
         ▼                    ▼                    ▼                    ▼
┌────────────────┐  ┌────────────────┐  ┌────────────────┐  ┌────────────────┐
│ MAIN_1_CoAnQi  │  │   QCalc.py     │  │ CondensedPhys  │  │ uqff_server.js │
│    .cpp/.exe   │  │  (Python)      │  │    ics.py      │  │  (HTTP 3141)   │
│                │  │                │  │                │  │                │
│ 107,019 lines  │  │  9,100+ lines  │  │  81,626 lines  │  │   504 lines    │
│ 446 modules    │  │  16 assoc.     │  │  12+ assoc.    │  │  imports       │
│ 15+ headers    │  │  programs      │  │  programs      │  │  index.js LIB  │
│ 6,688+ terms   │  │                │  │                │  │                │
└────────────────┘  └────────────────┘  └────────────────┘  └────────────────┘
         │                    │                    │                    │
         │                    └──────┬─────────────┘                    │
         │                           │                                  │
         ▼                           ▼                                  ▼
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                            shared_constants.h / .py / .js                            │
│                                                                                      │
│   Synchronized Constants:                                                           │
│   ├── rho_vac_UA = 7.09e-36 J/m³ (gravitational scale)                              │
│   ├── rho_vac_UA_field = 1e-27 J/m³ (field scale)                                   │
│   ├── GRADIENT_RATIO = 7.09e-9 (DPM field-gravity coupling)                         │
│   └── ... 50+ unified constants                                                     │
└─────────────────────────────────────────────────────────────────────────────────────┘

═══════════════════════════════════════════════════════════════════════════════════════
                    RECIRCULATION LOOP (Shared Files Persistence)
═══════════════════════════════════════════════════════════════════════════════════════

┌─────────────────────────────────────────────────────────────────────────────────────┐
│                     RECIRCULATION PERSISTENCE LOOP                                   │
│                                                                                      │
│   FLOW: APIFetch → IPData.py → CALCULATORS → OPData.py → OutputData → RECALL       │
│                                                                                      │
│   1. bodies_*.csv → IPData.py (input parameters, 431 lines)                         │
│   2. Calculators process → OPData.py (output results, 327 lines)                    │
│   3. OPData.py → uqff_results.json (JSON storage)                                   │
│   4. uqff_results.json → CondensedPhysics_OutputData.py (RECALL STORAGE)            │
│   5. Session Logger (Tab 9) → session_*.json (GUI state persist)                    │
│   6. MAIN_1_CoAnQi → coAnQi_log_*.txt (computation logs)                            │
│   7. Query history → USER RECALL (back to Query Field)                              │
│                                                                                      │
│   ┌─────────────────────────────────────────────────────────────────────────────┐   │
│   │                    PERSISTENCE FILES                                         │   │
│   │   ├── bodies_*.csv (timestamped API fetch results)                          │   │
│   │   ├── IPData.py (input parameter storage)                                   │   │
│   │   ├── OPData.py (output results storage)                                    │   │
│   │   ├── uqff_results.json (JSON computation results)                          │   │
│   │   ├── CondensedPhysics_OutputData.py (RECALL STORAGE)                       │   │
│   │   ├── session_*.json (GUI session state)                                    │   │
│   │   ├── coAnQi_log_*.txt (C++ computation logs)                               │   │
│   │   └── ./uqff_data/*.json (RPC requests/responses)                           │   │
│   └─────────────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────────────┘

═══════════════════════════════════════════════════════════════════════════════════════
                VR/VM BACKEND LAYER (DEVELOPER SIDE - GPU HEAVY)
═══════════════════════════════════════════════════════════════════════════════════════

┌─────────────────────────────────────────────────────────────────────────────────────┐
│                     VR/VM BACKEND (Developer Side - GPU Heavy)                       │
│                                                                                      │
│   ** NOT the entry point - this is for DEVELOPER heavy simulations **               │
│   ** Runs in virtual space with virtual keyboard, virtual goggles **                │
│                                                                                      │
│   ┌────────────────────────────┐    ┌────────────────────────────────────┐          │
│   │     vr_runtime.cpp         │    │     physics_backend.cpp            │          │
│   │     (GPU Heavy VR)         │    │     (CPU Heavy Physics)            │          │
│   │                            │    │                                    │          │
│   │   - Virtual space render   │    │   - Self-Expanding Framework       │          │
│   │   - Virtual keyboard       │    │   - Self-Updating (learning rate)  │          │
│   │   - Virtual goggles (XR)   │    │   - Self-Simulating (time evol)    │          │
│   │   - Astro Graphics (GPU)   │    │   - Distributed Computing pool     │          │
│   └────────────────────────────┘    └────────────────────────────────────┘          │
│                            │                    │                                    │
│                            └────────┬───────────┘                                    │
│                                     ▼                                                │
│   ┌─────────────────────────────────────────────────────────────────────────────┐   │
│   │                    IPC LAYER (Named Pipes + SharedMem + gRPC)                │   │
│   │                                                                              │   │
│   │   Named Pipe: \\.\pipe\StarMagic_UQFF                                       │   │
│   │   SharedMem:  Low-latency VR frame data                                     │   │
│   │   gRPC:       Structured commands (optional deployment)                     │   │
│   │                                                                              │   │
│   │   Message Types (uqff_ipc.h, 403 lines):                                    │   │
│   │   - CALCULATE_FIELD → F_U computation                                       │   │
│   │   - CALCULATE_GRAVITY → Ug1-4 computation                                   │   │
│   │   - VR_FRAME_UPDATE → stream field data to VR                               │   │
│   │   - REGISTER_TERM → add dynamic physics term                                │   │
│   │   - UPDATE_PARAMETER → modify runtime parameter                             │   │
│   │   - SYNC_STATE → synchronize module states                                  │   │
│   └─────────────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────────────┘

═══════════════════════════════════════════════════════════════════════════════════════
                           FILE INVENTORY
═══════════════════════════════════════════════════════════════════════════════════════

| File | Purpose | Lines | Key Functions |
|------|---------|-------|---------------|
| `source2.cpp` | **PRINCIPAL GUI** (user starts here) | 11,058 | MainWindow, 21 tabs |
| `source2_mainwindow.h` | Qt MOC declarations | ~50 | Q_OBJECT macro |
| `source2_widgets_enhanced.h` | Enhanced widgets | 897+ | UQFFJavaScriptWidget |
| `source2_event_bus.h` | Inter-widget events | ~200 | UQFF_LOG_* macros |
| `equation_renderer.h` | Long-form equation display | 790 | IPC to QCalc.py |
| `shared_constants.h` | UQFF C++ constants | 351 | UQFF:: namespace |
| `shared_constants.py` | UQFF Python constants | 250 | UQFFConstants dataclass |
| `MAIN_1_CoAnQi.cpp` | C++ UQFF calculator | 107,019 | 446 modules, 6,688+ terms |
| `QCalc.py` | Python solver | 9,100+ | UnifiedFieldSolver |
| `CondensedPhysics.py` | Pure physics calculator | 81,626 | 176 calculator classes |
| `CondensedPhysics_OutputData.py` | **RECALL STORAGE** | ~500 | Query recall |
| `IPData.py` | Input parameter storage | 431 | Input schema |
| `OPData.py` | Output results storage | 327 | Output schema |
| `index.js` | **LIBRARY INDEX** (NOT calculator) | 23,790 | 106 systems, Ug1-4 |
| `uqff_server.js` | REST API (imports index.js) | 504 | HTTP 3141, polls files |
| `uqff_ftps_client.py` | FTPS client | 890+ | RFC 4217, TLS 1.3 |
| `APIFetch.py` | 55 API client | 1,722 | bodies_*.csv output |
| `vr_runtime.cpp` | VR/VM GPU runtime | ~5K | Developer VR backend |
| `physics_backend.cpp` | CPU physics server | ~12K | Developer physics backend |
| `ipc/uqff_ipc.h` | IPC layer | 403 | Named Pipes, SharedMem |

═══════════════════════════════════════════════════════════════════════════════════════
                           PORT ASSIGNMENTS
═══════════════════════════════════════════════════════════════════════════════════════

| Port | Service | Protocol | Description |
|------|---------|----------|-------------|
| 990 | FTPS Implicit | TLS | External secure (TLS from start) |
| 21 | FTPS Explicit | STARTTLS | External (upgrade to TLS) |
| 3141 | uqff_server.js | HTTP | REST API (π×1000) |
| 8443 | QCalc_API.py | HTTPS | FastAPI (optional) |
| N/A | Named Pipe | IPC | \\.\pipe\StarMagic_UQFF |
| N/A | SharedMem | IPC | Low-latency VR data |

═══════════════════════════════════════════════════════════════════════════════════════
                           ENVIRONMENT VARIABLES
═══════════════════════════════════════════════════════════════════════════════════════

```bash
# uqff_server.js
UQFF_PORT=3141                    # HTTP port (default: 3141 = π×1000)
UQFF_HOST=127.0.0.1               # Bind address
UQFF_FILE_RPC=true                # Enable file-based RPC for FTPS
UQFF_REQUEST_DIR=./uqff_data/requests
UQFF_RESPONSE_DIR=./uqff_data/responses
UQFF_POLL_INTERVAL=1000           # ms

# uqff_ftps_client.py
UQFF_FTPS_HOST=ftps.example.com   # FTPS server hostname
UQFF_FTPS_PORT=990                # FTPS port (990=implicit, 21=explicit)
UQFF_FTPS_USER=username
UQFF_FTPS_PASS=password
UQFF_FTPS_IMPLICIT=true           # Use implicit FTPS
UQFF_FTPS_VERIFY=true             # Verify server certificate
UQFF_FTPS_PASSIVE=true            # Use passive mode
UQFF_FTPS_DATA_DIR=/uqff_data/    # Remote data directory

# Grok AI Integration
XAI_API_KEY=your_grok_api_key
```

═══════════════════════════════════════════════════════════════════════════════════════
                           QUICK START
═══════════════════════════════════════════════════════════════════════════════════════

```bash
# 1. Start the JavaScript REST API server
node uqff_server.js

# 2. Test local computation
curl -X POST http://127.0.0.1:3141/api/compute \
  -H "Content-Type: application/json" \
  -d '{"system":"SgrA*", "M":8.155e36, "r":4.4e19}'

# 3. Run the C++ calculator
./build_msvc/Release/MAIN_1_CoAnQi.exe --batch "Sagittarius A*"

# 4. Run Python solver
python -c "from QCalc import UnifiedFieldSolver; s=UnifiedFieldSolver(); print(s.solve_example())"

# 5. Run FTPS client tests
python test_uqff_ftps_client.py
```

---

**CANONICAL DOCUMENT** - Version 2.0 - Matches ARCHITECTURE_FLOW_DIAGRAM.md v4.0
**Author:** Daniel T. Murphy  
**Framework:** UQFF Star-Magic v4.0 CANONICAL  
**Copyright:** © 2025-2026 Daniel T. Murphy - All Rights Reserved
