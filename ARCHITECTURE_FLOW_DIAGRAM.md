# Star-Magic UQFF Architecture Flow Diagram

> **Version:** 4.2.2 (CANONICAL - DO NOT DEVIATE)
> **Generated:** 2026-02-21
> **Updated:** 2026-02-21 (v4.2.2 Dual Bot Architecture: CoAnQi + Poseidon)
> **Author:** Daniel T. Murphy
> **CRITICAL:** This is the MASTER architecture document. All other docs must match.

---

## CANONICAL ARCHITECTURE RULES (MEMORIZE THESE)

1. **USER INPUT goes FIRST** → enters through `source2.cpp` (Principal GUI)
2. **source2.cpp** = Principal GUI application (user-facing, 21 tabs, Qt6) - USER-FACING
3. **source2(HEAD PROGRAM).cpp** = VR/VM developer backend (GPU-heavy, headless capable) - BACKEND
4. **index.js** = LIBRARY INDEX (NOT a calculator) - exports 106 systems for require()
5. **Recirculation Loop** = bodies_*.csv → IPData → Calculators → OPData → OutputData → RECALL
6. **Simultaneous Joint Pipeline** = All 5 calculators run in parallel via IPC layer (Phase 1-5 complete)
7. **Poseidon TaskBot** = Offline physics maintenance (read/write/compare/validate/update cross-platform)
8. **CoAnQi_bot** = MAIN_1_CoAnQi.cpp EXCLUSIVE specialist (PhysicsTerm mgmt, self-expand, self-update)

## IPC Pipeline Status (Phases 1-5 Complete + v4.2.1)

| Phase | Description | Status | Commit |
|-------|-------------|--------|--------|
| Phase 1 | IPC Pipeline Connection | ✅ Complete | 87168f3 |
| Phase 2 | Physics Backend Service Mode | ✅ Complete | 0b1e737 |
| Phase 3 | Full gRPC Implementation | ✅ Complete | 1e5a722 |
| Phase 4 | Astro Graphics IPC Integration | ✅ Complete | 3351f42 |
| Phase 5 | Full VR Experience | ✅ Complete | e84c434 |
| v3.1a | Cross-Platform IPC (NamedPipeChannel) | ✅ Complete | 8967469 |
| v3.1b | Self-Expanding Physics Backend | ✅ Complete | 81097a8 |
| v4.2.1 | Poseidon TaskBot Integration | ✅ Complete | 277f954 |
| v4.2.2 | Dual Bot Architecture (CoAnQi + Poseidon) | ✅ Complete | 7436b0c |

---

## Complete System Architecture

```
═══════════════════════════════════════════════════════════════════════════════════════════════════════════════
                              STAR-MAGIC UQFF COMPLETE CROSS-PLATFORM ARCHITECTURE
                              (USER → source2.cpp GUI FIRST, VR/VM = Developer Backend)
═══════════════════════════════════════════════════════════════════════════════════════════════════════════════

                                              ┌─────────────────────────┐
                                              │       USER INPUT        │ ◄─── STARTS HERE (FIRST)
                                              │  "Sagittarius A*"       │
                                              │  "M87 Black Hole"       │
                                              │  "NGC 3596"             │
                                              └───────────┬─────────────┘
                                                          │
                                                          ▼
┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                                  source2.cpp (PRINCIPAL GUI APPLICATION)                                     │
│                                     11,058 lines | Qt6 + VTK + Chromium | 21 Tabs                            │
│                                                                                                               │
│   Includes:                     Tabs: 🎛️ MAIN_1 | 🐍 QCalc.py | 🧮 SciCalc | 📓 Notebook | 📚 CondensedPhysics  │
│   source2_mainwindow.h              | 🤖 CoAnQi_bot | 🧠 SuperGrok4 | 🌌 3D Simulator | 📋 Session Logger    │
│   source2_widgets_enhanced.h       | ⚖️ Compare C++/Python | 📐 Equation Display | 🌐 JS Engine (13-21)      │
│   source2_event_bus.h                                                                                         │
│   equation_renderer.h                                                                                         │
│                                          ORCHESTRATES ALL COMPUTATION PATHS                                   │
└────────────────────────────────────────────────────────┬─────────────────────────────────────────────────────┘
                                                         │
              ┌──────────────────────────────────────────┼──────────────────────────────────────────┐
              │                                          │                                          │
              ▼                                          ▼                                          ▼
═══════════════════════════════════════════════════════════════════════════════════════════════════════════════
               API DATA FETCH                      COMPUTATION DISPATCH                   EXTERNAL BRIDGE
═══════════════════════════════════════════════════════════════════════════════════════════════════════════════

┌───────────────────────────────┐     ┌───────────────────────────────────────────────┐     ┌──────────────────┐
│         APIFetch.py           │     │            DISPATCH METHODS                   │     │   FTPS BRIDGE    │
│         1,722 lines           │     │                                               │     │                  │
│                               │     │   pybind11 ─── Python embedded in C++         │     │ uqff_ftps_client │
│   55 ASTRONOMICAL APIs:       │     │   HTTP:3141 ── REST to uqff_server.js         │     │     .py          │
│   ├── SIMBAD (15M objects)    │     │   QProcess ─── Terminal spawn (fallback)      │     │   890+ lines     │
│   ├── NED (300M objects)      │     │   IPC ──────── Named pipes, SharedMem         │     │                  │
│   ├── VizieR (20K+ catalogs)  │     │                                               │     │ Port 990 (FTPS)  │
│   ├── NASA APOD/NeoWs/DONKI   │     │   ┌───────────────────────────────────────┐   │     │ Port 21  (FTP+TLS│
│   ├── Gaia (1.8B stars)       │     │   │   SIMULTANEOUS JOINT OPERATION        │   │     │                  │
│   ├── LIGO/Virgo (GW events)  │     │   │   PIPELINE (all calculators run       │   │     │ RFC 4217, TLS 1.3│
│   └── Grok/OpenAI fallback    │     │   │   in parallel where possible)         │   │     └──────────────────┘
│                               │     │   └───────────────────────────────────────┘   │
│   OUTPUT:                     │     │                                               │
│   bodies_YYYYMMDD_HHMMSS.csv  │────►│   Cross-validation: UQFF vs MUGE vs GR       │
└───────────────────────────────┘     └───────────────────────────────────────────────┘
              │
              ▼
═══════════════════════════════════════════════════════════════════════════════════════════════════════════════
                                     CALCULATOR ECOSYSTEMS (w/ Associated Programs)
═══════════════════════════════════════════════════════════════════════════════════════════════════════════════

┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│  CALCULATOR #1: QCalc.py (Python Unified Field Solver) 9,100+ lines                                          │
│                                                                                                               │
│  CORE:       QCalc.py                  ─── UnifiedFieldSolver, 8 master equations, long_form_equations       │
│                                                                                                               │
│  ASSOCIATED PROGRAMS:                                                                                         │
│  ├── QCalc_API.py              ─── FastAPI REST wrapper (Port 8443 HTTPS)                                    │
│  ├── QCalc_Advanced.py         ─── Advanced physics methods                                                  │
│  ├── QCalc_validation.py       ─── Validation suite, automated testing                                       │
│  ├── QCalc_test.py             ─── Unit tests, integration tests                                             │
│  ├── QCalc_stat.py             ─── Statistical analysis utilities                                            │
│  ├── QCalc_stat_test.py        ─── Statistical test suite                                                    │
│  ├── QCalc_Performance.py      ─── Performance benchmarking                                                  │
│  ├── QCalc_core_uqff.py        ─── Core UQFF equations extraction                                            │
│  ├── QCalc_cpp_equations.py    ─── C++ equation conversions                                                  │
│  ├── QCalc_cpp_extracted.py    ─── Constants extracted from C++                                              │
│  ├── QCalc_js_extracted.py     ─── Constants extracted from JS                                               │
│  ├── QCalc_extracted.py        ─── General extractions                                                       │
│  ├── QCalc_Wolfram_Extensions.py ─ Wolfram symbolic math extensions                                          │
│  ├── QCalc_Phase1_Validation.py ── Phase 1 validation suite                                                  │
│  └── QCalc_test_SOURCE16_50.py ─── SOURCE16-50 test coverage                                                 │
│                                                                                                               │
│  DATA FLOW: APIFetch.py → IPData.py → QCalc.py → OPData.py → uqff_results.json                               │
└──────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│  CALCULATOR #2: CondensedPhysics.py (Pure Physics Calculator) 81,626 lines                                   │
│                                                                                                               │
│  CORE:       CondensedPhysics.py       ─── 176 calculator classes, 111 Model classes, UQFFMasterEquations    │
│                                                                                                               │
│  ASSOCIATED PROGRAMS:                                                                                         │
│  ├── CondensedPhysics_InputData.py  ─── Input data schema and validation                                     │
│  ├── CondensedPhysics_OutputData.py ─── Output storage (RECALL SOURCE)                                       │
│  ├── CondensedPhysics_Validation.py ─── Validation framework                                                 │
│  ├── Phase5_Consolidated.py         ─── 838 lines, 16 systems (Source16-25)                                  │
│  ├── Phase6_Consolidated.py         ─── 2,100+ lines, 11 systems (Source26-36)                               │
│  ├── Phase7_Consolidated.py         ─── 3,645 lines, 14 cosmological systems                                 │
│  ├── InputData.py                   ─── General input schema                                                 │
│  ├── OutputData.py                  ─── General output schema                                                │
│  ├── IPData.py                      ─── Input parameter storage (431 lines)                                  │
│  ├── OPData.py                      ─── Output results storage (327 lines)                                   │
│  ├── CoAnQi_Wrapper.py              ─── C++ to Python bridge wrapper                                         │
│  └── shared_constants.py            ─── Synchronized constants (250 lines)                                   │
│                                                                                                               │
│  DATA FLOW: bodies_*.csv → CondensedPhysics.py → CondensedPhysics_OutputData.py (RECALL STORAGE)             │
└──────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│  CALCULATOR #3: MAIN_1_CoAnQi.cpp (C++ Native Physics Engine) 107,019 lines                                  │
│                                                                                                               │
│  CORE:       MAIN_1_CoAnQi.cpp         ─── 446 modules, 6,688+ physics terms, 16-option interactive menu     │
│                                                                                                               │
│  ASSOCIATED PROGRAMS:                                                                                         │
│  ├── source1.cpp - source173.cpp    ─── 173 physics module files (116 integrated into SOURCE blocks)        │
│  ├── shared_constants.h             ─── UQFF:: namespace constants (351 lines)                              │
│  ├── observational_systems_config.h ─── 35+ astrophysical system parameters                                 │
│  ├── MUGE.h                         ─── Modified Unified Gravity Equations                                  │
│  ├── uqff_self_expanding.h          ─── Self-expanding framework                                            │
│  ├── uqff_dual_physics.h            ─── Dual physics validation                                             │
│  ├── uqff_tracing.h                 ─── Computation tracing                                                 │
│  ├── uqff_cross_platform.h          ─── Cross-platform compatibility                                        │
│  ├── uqff_constants.h               ─── Core constants                                                      │
│  ├── wolfram_wstp_runtime.h         ─── Wolfram WSTP integration                                            │
│  ├── UQFFResultsWidget.h            ─── Results display widget                                               │
│  ├── UQFFSource10.h                 ─── Source10 integration                                                │
│  ├── FluidSolver.h                  ─── Navier-Stokes fluid solver                                          │
│  ├── CelestialBody.h                ─── Celestial body data structures                                      │
│  ├── csv_body_reader.h              ─── bodies.csv reader                                                   │
│  └── UnitTests.h                    ─── Unit test framework                                                  │
│                                                                                                               │
│  OUTPUT: coAnQi_log_*.txt (computation logs), stdout (interactive menu)                                      │
└──────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│  LIBRARY INDEX: index.js (JavaScript UQFF Library) 23,790 lines ◄── NOT A CALCULATOR, IT'S A LIBRARY INDEX   │
│                                                                                                               │
│  CORE:       index.js                  ─── 106 astrophysical systems, CONSTANTS, COUPLING, system functions  │
│                                                                                                               │
│  ASSOCIATED PROGRAMS:                                                                                         │
│  ├── uqff_server.js                 ─── REST API server (Port 3141 HTTP) - CALLS index.js                   │
│  ├── automated_legacy_converter.js  ─── Legacy format conversion                                             │
│  └── (exports module for require())                                                                          │
│                                                                                                               │
│  USAGE: const UQFF = require('./index.js'); UQFF.computeSagA(params);                                        │
│  Server: uqff_server.js imports index.js → /api/compute, /api/batch, /api/systems                            │
└──────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Recirculation Loop (Shared Files: Frontend ↔ Backend ↔ Storage)

```
═══════════════════════════════════════════════════════════════════════════════════════════════════════════════
                           RECIRCULATION THROUGH SHARED FILES
═══════════════════════════════════════════════════════════════════════════════════════════════════════════════

   ┌─────────────────────┐    ┌─────────────────────┐    ┌─────────────────────┐    ┌─────────────────────┐
   │  bodies_*.csv       │    │  uqff_results.json  │    │  session_*.json     │    │  coAnQi_log_*.txt   │
   │  (API fetch output) │    │  (OPData storage)   │    │  (GUI state persist)│    │  (computation logs) │
   └─────────┬───────────┘    └─────────┬───────────┘    └─────────┬───────────┘    └─────────┬───────────┘
             │                          │                          │                          │
   ┌─────────┴───────────┐    ┌─────────┴───────────┐    ┌─────────┴───────────┐    ┌─────────┴───────────┐
   │  IPData.py          │    │  OPData.py          │    │CondensedPhysics_    │    │ shared_constants.h  │
   │  (input params)     │    │  (output results)   │    │ OutputData.py       │    │ shared_constants.py │
   │  431 lines          │    │  327 lines          │    │ (RECALL storage)    │    │ (synchronized)      │
   └─────────┬───────────┘    └─────────┬───────────┘    └─────────┬───────────┘    └─────────┬───────────┘
             │                          │                          │                          │
             └────────────┬─────────────┴──────────────────────────┼──────────────────────────┘
                          │                                        │
                          ▼                                        ▼
              ┌─────────────────────────────────────────────────────────────────────────────────────────┐
              │                         RECIRCULATION LOOP FLOW                                          │
              │                                                                                          │
              │   USER QUERY ──► APIFetch.py ──► bodies_*.csv ──► IPData.py ──► CALCULATORS             │
              │        ▲                                                              │                  │
              │        │                                                              ▼                  │
              │   RECALL ◄─── Session Logger (Tab 9) ◄─── OPData.py ◄─── uqff_results.json              │
              │        │                                                              ▲                  │
              │        │                                                              │                  │
              │        └───────────────── CondensedPhysics_OutputData.py ◄────────────┘                  │
              │                         (Stores: solutions, equation lists,                              │
              │                          simulation sets for USER RECALL)                                │
              └─────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## VR/VM Backend Layer (Developer Side - GPU Heavy)

```
═══════════════════════════════════════════════════════════════════════════════════════════════════════════════
               VR/VM BACKEND LAYER (DEVELOPER SIDE - GPU HEAVY SIMULATIONS IN VIRTUAL SPACE)
                    IPC Communication via Named Pipes + SharedMemory + gRPC
═══════════════════════════════════════════════════════════════════════════════════════════════════════════════

┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                    ipc/ LAYER - SIMULTANEOUS JOINT OPERATION PIPELINE                                        │
│                                                                                                               │
│   ┌───────────────────────┐    ┌───────────────────────┐    ┌───────────────────────┐                        │
│   │     uqff_ipc.h        │    │   python_bridge.h     │    │   physics_service.h   │                        │
│   │     515 lines (v3.1)  │    │   (pybind11 bridge)   │    │   470 lines (v3.1)    │                        │
│   │                       │    │                       │    │                       │                        │
│   │ MessageTypes:         │    │ Embeds Python in C++: │    │ Self-Expand (v3.1):   │                        │
│   │ - CALCULATE_FIELD     │    │ - CondensedPhysics.py │    │ - onRegisterTerm()    │                        │
│   │ - CALCULATE_GRAVITY   │    │ - QCalc.py            │    │ - evaluateDynamicTerms│                        │
│   │ - VR_FRAME_UPDATE     │    │ - Phase5-7.py         │    │ Self-Update (v3.1):   │                        │
│   │ - REGISTER_TERM       │    │                       │    │ - onUpdateParameter() │                        │
│   │ - UPDATE_PARAMETER    │    │                       │    │ - κ, [SSq], β_i tuning│                        │
│   │ - SIM_START/FRAME     │    │                       │    │ Self-Simulate (v3.1): │                        │
│   │ - SIM_COMPLETE        │    │                       │    │ - startSimulation()   │                        │
│   └───────────────────────┘    └───────────────────────┘    └───────────────────────┘                        │
│                                                                                                               │
│   Channels: SharedMemoryChannel + NamedPipeChannel + GrpcChannel (all cross-platform)                        │
│   Windows: Named Pipes (CreateNamedPipe) | Linux/macOS: Unix Domain Sockets                                  │
└──────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
                                              │
                          ┌───────────────────┴───────────────────┐
                          │                                       │
                          ▼                                       ▼
┌─────────────────────────────────────────────┐    ┌─────────────────────────────────────────────┐
│   source2(HEAD PROGRAM).cpp                 │    │   physics_backend.cpp (CPU Heavy Server)    │
│   (VR/VM BACKEND - NOT A GUI)               │    │   (Headless physics computation)            │
│   2,625 lines | VR namespace merged         │    │                                             │
│   PURPOSE: Heavy GPU simulations in VR      │    │   PURPOSE: CPU-bound physics computation    │
│   - Virtual space rendering                 │    │   - Self-Expanding (onRegisterTerm)         │
│   - Virtual keyboard input                  │    │   - Self-Updating (onUpdateParameter)       │
│   - Virtual goggles (OpenXR headset)        │    │   - Self-Simulating (startSimulation)       │
│   - Astro Graphics Program (GPU tasking)    │    │   - Distributed Computing pool              │
│   - --service flag for headless mode        │    │                                             │
│   ASSOCIATED HEADERS (vr/ directory):       │    │   Handles IPC messages from VR backend:     │
│   ├── vr_runtime.h (merged content)         │    │   - CALCULATE_FIELD → F_U computation      │
│   ├── openxr_session.h                      │    │   - VR_FRAME_UPDATE → stream field data    │
│   ├── vulkan_compositor.h                   │    │   - REGISTER_TERM → add physics term       │
│   ├── task_bot.h (voice/gesture bot)        │    │   - SYNC_STATE → synchronize modules       │
│   ├── poseidon_task_bot.h (general maint)   │    │                                             │
│   ├── CoAnQi_bot.h (MAIN_1 specialist)      │    │                                             │
│   └── astro_graphics.h                      │    │                                             │
│                                             │    │                                             │
│   [Lightweight: ~5K lines | GPU-bound]      │    │   [Heavy: ~12K lines | CPU-bound | Async]   │
└─────────────────────────────────────────────┘    └─────────────────────────────────────────────┘

┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                                POSEIDON TASK BOT (v4.2.1 - Offline Physics Maintenance)                      │
│                                                                                                               │
│   Files: vr/poseidon_task_bot.h, vr/poseidon_task_bot.cpp, poseidon_maintenance.py                           │
│                                                                                                               │
│   CAPABILITIES:                                                                                               │
│   ├── ProcessNewPhysics()    : Read/parse/integrate new physics equations                                    │
│   ├── CompareAllCalculators(): Cross-validate C++/Python/JS implementations                                  │
│   ├── ValidatePhysics()      : Run QCalc_validation.py + CondensedPhysics_Validation.py + UnitTests.h        │
│   ├── UpdateAndExpandPhysics(): Register dynamic terms via Self-Expanding Backend (v3.1)                    │
│   ├── SyncConstantsAcrossLanguages(): shared_constants.h ↔ .py ↔ index.js                                   │
│   ├── RegenerateExtractedFiles(): QCalc_cpp_extracted.py, QCalc_js_extracted.py                             │
│   ├── BackupAllPhysicsFiles(): Timestamped backups before any change                                        │
│   ├── FTPSPushMaintenanceBundle(): Secure offline bundle via uqff_ftps_client.py                            │
│   └── ExecuteCommand()       : Voice/script command interface                                                │
│                                                                                                               │
│   INTEGRATION: Uses physics_service.h (v3.1), python_bridge.h (pybind11), NamedPipeChannel                  │
│   OFFLINE-FIRST: All operations work without internet; FTPS only for local/network share                    │
└──────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                       COANQI BOT (v4.2.2 - MAIN_1_CoAnQi.cpp EXCLUSIVE Specialist)                           │
│                                                                                                               │
│   Files: vr/CoAnQi_bot.h, vr/CoAnQi_bot.cpp, task_bot_maintenance.py                                         │
│                                                                                                               │
│   PURPOSE: Dedicated maintenance bot that services MAIN_1_CoAnQi.cpp EXCLUSIVELY.                            │
│            Provides continuity for the 107K+ line C++ physics engine (6,688+ PhysicsTerms).                  │
│                                                                                                               │
│   CAPABILITIES:                                                                                               │
│   ├── RegisterPhysicsTerm()  : Add new PhysicsTerm to MAIN_1_CoAnQi registry                                │
│   ├── UpdatePhysicsTerm()    : Modify existing term parameters dynamically                                   │
│   ├── ValidatePhysicsTerm()  : Validate term against observational data                                     │
│   ├── InjectDynamicTerm()    : Runtime term injection (Self-Expanding v3.1)                                 │
│   ├── OptimizeParameters()   : Self-updating parameter optimization (gradient descent)                      │
│   ├── CloneAndMutate()       : Self-cloning with mutation for parameter sensitivity                         │
│   ├── CalculateSystem()      : Compute UQFF physics for single system                                       │
│   ├── CalculateAllSystems()  : Batch process all 26+ predefined systems                                     │
│   ├── RunSimulation()        : Execute one of 6 simulation modes                                            │
│   ├── PerformStatisticalAnalysis(): Full statistical suite (mean, stddev, correlation)                      │
│   ├── CompareWithQCalc()     : Cross-validate with QCalc.py results                                         │
│   ├── CompareWithCondensedPhysics(): Cross-validate with CondensedPhysics.py                                │
│   └── ExecuteMenuOption()    : Programmatically execute MAIN_1_CoAnQi menu options                          │
│                                                                                                               │
│   DISTINCTION FROM POSEIDON:                                                                                  │
│   ├── CoAnQi_bot = SPECIALIZED for MAIN_1_CoAnQi.cpp ONLY (PhysicsTerm mgmt, simulations)                   │
│   └── Poseidon   = GENERAL CONTRACTOR for entire codebase (all languages, cross-platform)                   │
│                                                                                                               │
│   INTEGRATION: Uses python_bridge.h (pybind11 → task_bot_maintenance.py), physics_service.h, IPC            │
└──────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Port Assignments

| PORT | PROTOCOL | SERVICE | DESCRIPTION |
|------|----------|---------|-------------|
| 990 | FTPS | uqff_ftps_client.py | Implicit FTPS (TLS from connection start) |
| 21 | FTP+TLS | uqff_ftps_client.py | Explicit FTPS (upgrades via STARTTLS) |
| 3141 | HTTP | uqff_server.js | REST API (π×1000) - /api/compute, /api/batch, /api/systems |
| 8443 | HTTPS | QCalc_API.py (FastAPI) | Python REST API - /api/uqff (optional) |
| N/A | IPC | Named Pipe | \\.\pipe\StarMagic_UQFF (VR ↔ Physics Backend) |
| N/A | IPC | SharedMemory | Low-latency field data (VR real-time frames) |
| N/A | IPC | gRPC | Structured commands (optional deployment) |

---

## Complete File Inventory

| CATEGORY | COUNT | FILES |
|----------|-------|-------|
| C++ Physics | 173 | source1.cpp - source173.cpp (116 integrated into MAIN_1_CoAnQi.cpp) |
| C++ Headers | 32 | *.h files (shared_constants, IPC, VR, widgets, etc.) |
| Python Calculators | 16 | QCalc*.py ecosystem |
| Python Support | 30+ | CondensedPhysics*.py, Phase*_Consolidated.py, IP/OPData.py, APIFetch.py |
| JavaScript | 3 | index.js (LIBRARY), uqff_server.js, automated_legacy_converter.js |
| IPC Layer | 3 | ipc/uqff_ipc.h, python_bridge.h, physics_service.h |
| VR Layer | 7 | vr/*.h (runtime, openxr, vulkan, task_bot, poseidon_task_bot, CoAnQi_bot, astro_graphics) |
| Modules System | 10+ | modules/*.py (loader, interface, gaming/*, debug/*) |
| Config/Data | 20+ | *.json, *.csv, observational_systems_config.h |

---

## Quick Reference: Data Flow

```
USER QUERY → source2.cpp GUI → APIFetch.py → bodies_*.csv → IPData.py
                    ↓
            ┌───────┴───────┐
            ▼               ▼
      CALCULATORS     VR/VM Backend
   (parallel dispatch)  (GPU heavy)
            │               │
            ▼               ▼
        OPData.py ← ← ← physics_backend.cpp
            │
            ▼
    uqff_results.json
            │
            ▼
CondensedPhysics_OutputData.py
            │
            ▼
    Session Logger (Tab 9)
            │
            ▼
    USER RECALL (back to top)
```

---

*CANONICAL DOCUMENT - Version 4.2.2 - DO NOT DEVIATE*
*Updated: 2026-02-21 (v4.2.2 Dual Bot Architecture: CoAnQi + Poseidon) by Daniel T. Murphy*
