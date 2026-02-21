# VR Compositor + Physics Backend Architecture Strategy

> **Version:** 2.1 CANONICAL (matches ARCHITECTURE_FLOW_DIAGRAM.md v4.1)
> **Date:** February 21, 2026  
> **Updated:** Phase 3 gRPC Complete  
> **Author:** GitHub Copilot (Claude Opus 4.5) + Daniel Murphy  
> **Branch:** master  
> **Repository:** Daniel8Murphy0007/Star-Magic  

---

## CANONICAL ARCHITECTURE RULES

1. **USER INPUT goes FIRST** → enters through `source2.cpp` (Principal GUI)
2. **source2.cpp** = Principal GUI application (user-facing, 21 tabs, Qt6)
3. **source2(HEAD PROGRAM).cpp** = VR/VM developer backend (GPU-heavy, virtual space, 2,452 lines)
4. **physics_backend.cpp** = CPU-bound physics server (headless)
5. **index.js** = LIBRARY INDEX (NOT a calculator) - exports 106 systems for require()

## IPC Pipeline Status

| Phase | Status | Description |
|-------|--------|-------------|
| Phase 1 | ✅ COMPLETE | IPC Pipeline Connection (SharedMemory + NamedPipe) |
| Phase 2 | ✅ COMPLETE | Physics Backend Service Mode (--service flag) |
| Phase 3 | ✅ COMPLETE | Full gRPC Implementation (port 50051) |
| Phase 4 | ✅ COMPLETE | Astro Graphics IPC Integration |

---

## Executive Summary

This document formalizes the architectural evolution of Star-Magic UQFF consisting of:

1. **Principal GUI** (`source2.cpp`) - User-facing Qt6 application with 21 tabs
2. **VR Runtime Layer** (`vr_runtime.cpp`) - Developer VR/VM backend for GPU-heavy simulations
3. **Physics Backend Service** (`physics_backend.cpp`) - Headless computation server
4. **Core Computation Engine** (`MAIN_1_CoAnQi.cpp`) - Library mode or interactive calculator

---

## 1. Dual-Program Role Clarification

### 1.1 source2.cpp → Principal GUI (USER STARTS HERE)

| Aspect | Description |
|--------|-------------|
| **Lines** | 11,058 |
| **Purpose** | **PRINCIPAL GUI** - User-facing entry point |
| **Dependencies** | Qt6 + VTK + Chromium |
| **Tabs** | 21 tabs including calculators, AI bots, session logger |
| **Thread Model** | Qt Event Loop + Worker Threads |

**Key Characteristics:**
- **THIS IS THE ENTRY POINT** - User launches this application first
- User query field for astronomical object names
- APIFetch.py integration for 55 API sources
- Dispatches to calculators (MAIN_1, QCalc.py, CondensedPhysics.py, uqff_server.js)
- Session Logger (Tab 9) for RECALL functionality
- Cross-validation between C++ and Python (Tab 10)

### 1.2 source2(HEAD PROGRAM).cpp → VR/VM Developer Backend

| Aspect | Description |
|--------|-------------|
| **Lines** | 2,452 |
| **Purpose** | VR/VM developer backend - GPU-heavy simulations in virtual space |
| **Dependencies** | OpenXR, Vulkan, Astro Graphics |
| **Thread Model** | GPU-bound, real-time priority |
| **IPC** | SharedMemoryChannel, GrpcChannel (via ipc/uqff_ipc.h) |
| **Service Mode** | --service flag for headless operation |

**Key Characteristics:**
- **NOT the user entry point** - Developer backend for heavy GPU work
- Runs in virtual space (VR headset environment)
- VR namespace with VRRuntime class merged from vr_runtime.cpp
- Virtual keyboard input, virtual goggles (OpenXR)
- Manages astronomical graphics program (GPU tasking)
- Task-Oriented Bot for automation (voice/gesture input)
- Communicates with physics_backend.cpp via IPC (Phase 1-3 complete)

### 1.3 physics_backend.cpp → CPU-Bound Physics Server

| Aspect | Description |
|--------|-------------|
| **Lines** | ~12,000 (projected) |
| **Purpose** | Headless computation server for heavy physics calculations |
| **Dependencies** | MAIN_1_CoAnQi.cpp (library mode), pybind11, gRPC |
| **Thread Model** | Distributed Computing Pool |

**Key Characteristics:**
- Self-Expanding Framework (register terms at runtime)
- Self-Updating (parameter optimization via learning rate)
- Self-Simulating (time evolution, validation pipeline)
- Handles IPC requests from vr_runtime.cpp
- Can run headless on server machines

---

## 2. Generational Constant Evolution Philosophy

### 2.1 Understanding Constant Nuances

**CRITICAL:** Similar constants across source files are NOT duplicates. They represent:

1. **Generational Evolution:** Constants evolved during development, with older iterations preserved for:
   - Recognition patterns in C++ development
   - Baseline tracking and validation
   - Method-specific precision requirements

2. **Multi-Solution Physics:** UQFF produces multiple simultaneous solutions:
   - `c` (speed of light) may have quantum corrections in some modules
   - `G` (gravitational constant) has running coupling variations
   - `κ` (decay constant) ranges from 0.0005/day to 0.0015/day per system

3. **Timestamp Hardcoding:** Each `sourceXX.cpp` contains hardcoded:
   - Creation timestamps
   - Baseline measurements
   - Equation evolution history
   - Validation checkpoints

### 2.2 Constant Namespace Strategy

```cpp
// GLOBAL (shared_constants.h)
namespace UQFF::Constants {
    constexpr double c = 2.99792458e8;      // CODATA 2018 exact
    constexpr double G = 6.67430e-11;       // CODATA 2018
}

// MODULE-SPECIFIC (preserved for generational tracking)
namespace InfoParadox {
    constexpr double c = 2.99792458e8;      // Same value, different context
    constexpr double G = 6.67430e-11;       // May have local corrections
}

// SYSTEM-SPECIFIC (observational_systems_config.h)
namespace SGR1745 {
    constexpr double kappa = 0.0005;        // κ calibration for this magnetar
}
```

### 2.3 Files Preserving Generational History

| File Range | Content | Timestamps |
|------------|---------|------------|
| `source13.cpp` - `source35.cpp` | Original UQFF systems | Nov 2025 |
| `source36.cpp` - `source80.cpp` | Expanded physics | Nov-Dec 2025 |
| `source81.cpp` - `source116.cpp` | Advanced modules | Dec 2025 |
| `source117.cpp` - `source173.cpp` | Validation batches | Jan 2026 |
| `source174.cpp` - `source200.cpp` | Wolfram/Grok integration | Jan-Feb 2026 |

**Policy:** All sourceXX.cpp files are preserved. Their hardcoded baselines form the historical record of UQFF evolution.

---

## 3. VR Architecture Data Flow (CANONICAL)

```
┌──────────────────────────────────────────────────────────────────────────────────────┐
│                              USER (starts here)                                       │
│                          [Desktop, CLI, Web Browser]                                  │
└───────────────────────────────────┬──────────────────────────────────────────────────┘
                                    │
                                    ▼
┌──────────────────────────────────────────────────────────────────────────────────────┐
│                 source2.cpp (PRINCIPAL GUI APPLICATION)                              │
│                           11,058 lines | Qt6 + VTK + Chromium                        │
│                                                                                       │
│  ┌─────────────────────────────────────────────────────────────────────────────────┐ │
│  │                           USER QUERY FIELD                                       │ │
│  │   "Sagittarius A*", "M87", "Betelgeuse", "NGC 3596"...                          │ │
│  └─────────────────────────────────────────────────────────────────────────────────┘ │
│                                    │                                                  │
│                                    ▼                                                  │
│  ┌─────────────────────────────────────────────────────────────────────────────────┐ │
│  │                    API FETCH LAYER (APIFetch.py)                                 │ │
│  │   55 APIs: SIMBAD → NASA → VizieR → NED → Gaia → Grok fallback                  │ │
│  │   Output: bodies_YYYYMMDD_HHMMSS.csv → IPData.py                                │ │
│  └─────────────────────────────────────────────────────────────────────────────────┘ │
│                                    │                                                  │
│                                    ▼                                                  │
│  ┌─────────────────────────────────────────────────────────────────────────────────┐ │
│  │             SIMULTANEOUS JOINT OPERATION PIPELINE                                │ │
│  │                    (parallel dispatch to calculators)                            │ │
│  │                                                                                  │ │
│  │  ┌────────────────┐ ┌────────────────┐ ┌────────────────┐ ┌────────────────┐    │ │
│  │  │ MAIN_1_CoAnQi  │ │   QCalc.py     │ │ CondensedPhys  │ │ uqff_server.js │    │ │
│  │  │ (C++ 107K)     │ │ (Python 9K)    │ │ (Python 81K)   │ │ (HTTP:3141)    │    │ │
│  │  └────────────────┘ └────────────────┘ └────────────────┘ └────────────────┘    │ │
│  └─────────────────────────────────────────────────────────────────────────────────┘ │
│                                    │                                                  │
│                                    ▼                                                  │
│  ┌─────────────────────────────────────────────────────────────────────────────────┐ │
│  │                    RECIRCULATION LOOP                                            │ │
│  │   OPData.py → uqff_results.json → CondensedPhysics_OutputData.py → RECALL       │ │
│  └─────────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                       │
│  [PRINCIPAL: User-facing, 21 tabs, query/recall workflow]                            │
└──────────────────────────────────────────────────────────────────────────────────────┘
                   │
                   │ Optional IPC for VR/GPU workloads
                   ▼
┌──────────────────────────────────────────────────────────────────────────────────────┐
│                     VR/VM BACKEND LAYER (Developer Side)                             │
│                                                                                       │
│  ┌─────────────────────────────────────────────────────────────────────────────────┐ │
│  │                    IPC LAYER (Named Pipes + SharedMem + gRPC)                   │ │
│  │   Named Pipe: \\.\pipe\StarMagic_UQFF                                           │ │
│  │   SharedMem: Low-latency VR frame data                                          │ │
│  │   Messages: CALCULATE_FIELD, VR_FRAME_UPDATE, REGISTER_TERM, etc.              │ │
│  └─────────────────────────────────────────────────────────────────────────────────┘ │
│                         │                     │                                       │
│                         ▼                     ▼                                       │
│  ┌─────────────────────────────┐  ┌─────────────────────────────────────────────┐   │
│  │  source2(HEAD PROGRAM).cpp  │  │     physics_backend.cpp                      │   │
│  │  (VR/VM BACKEND - 2,452 ln)│  │     (CPU-bound Physics)                      │   │
│  │                             │  │                                              │   │
│  │  - Virtual space render     │  │  - Self-Expanding (registerDynamicTerm)     │   │
│  │  - Virtual keyboard         │  │  - Self-Updating (learning rate opt)        │   │
│  │  - Virtual goggles (XR)     │  │  - Self-Simulating (time evolution)         │   │
│  │  - Astro Graphics (GPU)     │  │  - Distributed Computing pool               │   │
│  │  - Task Bot (voice/gesture) │  │  - ML Integration (PyTorch)                 │   │
│  │  - --service flag supported │  │                                              │   │
│  │  [2,452 lines, GPU-bound]   │  │  [~12K lines, CPU-bound, async]             │   │
│  └─────────────────────────────┘  └─────────────────────────────────────────────┘   │
│                                                                                       │
│  [DEVELOPER BACKEND: Heavy GPU/CPU simulations in virtual space]                     │
└──────────────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼ Library calls
┌──────────────────────────────────────────────────────────────────────────────────────┐
│                 MAIN_1_CoAnQi.cpp (Core Physics Engine)                              │
│                                                                                       │
│  ┌─────────────────────────────────────────────────────────────────────────────────┐ │
│  │  6,688+ PhysicsTerm classes | 116 SOURCE blocks | 121+ systems                  │ │
│  │  SOURCE1-SOURCE116 integrated modules (446 unique modules)                      │ │
│  │  source1.cpp - source173.cpp (173 files)                                        │ │
│  └─────────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                       │
│  [Core: 107K lines, 16-option menu, standalone or library mode]                      │
└──────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 4. Implementation Phases

### Phase 1: IPC Foundation ✅ COMPLETE
- Created `ipc/uqff_ipc.h` and `ipc/uqff_ipc.cpp`
- Defined message types: `CALCULATE_FIELD`, `REGISTER_TERM`, `UPDATE_PARAMETER`, etc.
- Implemented `SharedMemoryChannel` for low-latency VR data
- Implemented `NamedPipeChannel` for cross-process communication
- Connected all 5 Principal Programs via IPC
- **Commit:** 87168f3

### Phase 2: Physics Backend Service Mode ✅ COMPLETE
- Added `#include "ipc/physics_service.h"` to source2(HEAD PROGRAM).cpp
- Added `--service` flag check in main() for headless operation
- VR namespace with VRRuntime class functional
- **Commit:** 0b1e737

### Phase 3: gRPC Implementation ✅ COMPLETE
- Created `ipc/uqff_service.proto` with message definitions
- Implemented full `GrpcChannel` class in uqff_ipc.h/cpp
- Added `calculateField()`, `getStatus()` methods via gRPC
- Enabled gRPC in `physics_service.h` (enable_grpc = true)
- Port 50051 for gRPC service
- **Commit:** 1e5a722

### Phase 4: Astro Graphics IPC Integration ✅ COMPLETE
- Added IPC integration to `vr/astro_graphics.h` and `vr/astro_graphics.cpp`
- Implemented `calculateFieldViaIPC()` and `calculateAllFieldsViaIPC()` methods
- Added `loadAstroGraphics()` method to VRRuntime class in `source2(HEAD PROGRAM).cpp`
- Connected AstroGraphics to physics_channel_ for real-time field calculations
- Catalog objects can now get F_U values from physics backend via IPC
- **Location:** `source2(HEAD PROGRAM).cpp` (VRRuntime class, lines 240-280)

### Phase 5: Full VR Experience ⏳ PENDING
- Complete OpenXR input handling
- Gesture → Physics → Render loop
- Task Bot automation
- Multi-machine deployment
- **Production deployment**

---

## 5. Current Environmental Conditions

### 5.1 Build Environment
- **OS:** Windows 11 Pro
- **Compiler:** MSVC 2022 (14.44.35219)
- **C++ Standard:** C++20
- **Build System:** CMake + Visual Studio 17 2022
- **Package Manager:** vcpkg (x64-windows)

### 5.2 Available Dependencies
| Dependency | Status | Purpose |
|------------|--------|---------|
| Qt6 | ✅ Found | GUI framework |
| VTK 9.3 | ✅ Found | Scientific visualization |
| CURL | ✅ Found | HTTP/API access |
| SQLite3 | ✅ Found | Local caching |
| OpenCV 4.11 | ✅ Found | Computer vision |
| AWS SDK 1.11 | ✅ Found (disabled) | Cloud sync |
| pybind11 | ⚠️ Conditional | Python integration |
| Wolfram WSTP | ⚠️ Optional | Symbolic math |

---

## 6. Port Assignments

| Port | Service | Protocol | Description |
|------|---------|----------|-------------|
| 990 | FTPS Implicit | TLS | External secure (TLS from start) |
| 21 | FTPS Explicit | STARTTLS | External (upgrade to TLS) |
| 3141 | uqff_server.js | HTTP | REST API (π×1000) |
| 8443 | QCalc_API.py | HTTPS | FastAPI (optional) |
| N/A | Named Pipe | IPC | \\.\pipe\StarMagic_UQFF |
| N/A | SharedMem | IPC | Low-latency VR data |

---

## 7. File Utilization Plan

**Goal:** Tie/utilize every file in the worktree within this framework.

| File Category | Count | Integration Point |
|---------------|-------|-------------------|
| `source1.cpp` - `source173.cpp` | 173 | Physics engine via SOURCE blocks |
| `*_wolfram.cpp` | 50+ | Wolfram export via source174 |
| `Core/Modules/*.cpp` | 20+ | Dynamic term registration |
| Python (`*.py`) | 30+ | pybind11 integration |
| JavaScript (`*.js`) | 3 | index.js LIBRARY, uqff_server.js REST |
| Headers (`*.h`) | 40+ | Shared constants, IPC layer |
| Config (`*.json`, `*.csv`) | 20+ | System parameters, bodies data |

---

## 8. Self-* Capability Preservation

| Capability | Current Implementation | VR Architecture |
|------------|------------------------|-----------------|
| **Self-Expand** | `registerDynamicTerm()` | IPC: `onRegisterTerm()` from VR |
| **Self-Update** | `setDynamicParameter()`, learning rate | IPC: `onUpdateParameter()` loop |
| **Self-Simulate** | `validation_pipeline()`, time evolution | IPC: `onSimulationStep()` streaming |

---

## Appendix: Key Constants Overview

```cpp
// Core UQFF Constants (shared_constants.h)
c       = 2.99792458e8    // Speed of light (m/s)
G       = 6.67430e-11     // Gravitational constant (m³/kg·s²)
hbar    = 1.054571817e-34 // Reduced Planck constant (J·s)
M_sun   = 1.98892e30      // Solar mass (kg)

// Calibrated UQFF Parameters
κ       = 0.0005/day      // Decay constant (system-dependent)
[SSq]   = 0.57            // Spatial scale quotient
H_SCm   ≈ 0.99            // Superconductive modulation
U_UA    ≈ 0.0001          // Universal acceleration

// Generationally-varied (preserved in sourceXX.cpp)
k_η     = 10⁻¹¹³          // String tension (varies by module)
β_i     ≈ 0.603           // Buoyancy scaling (system-specific)
```

---

**CANONICAL DOCUMENT** - Version 2.0 - Matches ARCHITECTURE_FLOW_DIAGRAM.md v4.0  
*Document generated: February 21, 2026*  
*Next action: Git commit, then Phase 0 build fixes*
