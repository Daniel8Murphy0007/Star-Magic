# Complete File Inventory & Wiring Diagram
## Star-Magic UQFF Unified System
**Generated:** March 2, 2026 | **Status:** 700+ files catalogued across 3 tiers

---

## EXECUTIVE SUMMARY

| Component | File Count | Status | Critical Files |
|-----------|-----------|--------|-----------------|
| **Tier 1: GUI** | 8 files | Functional | source2.cpp, source2_*widgets*.h |
| **Tier 2: Backend** | 180+ files | Functional | source2(HEAD PROGRAM).cpp, MAIN_1_CoAnQi.cpp, index.js |
| **Tier 3: IPC** | 15 files | Defined | ipc/*.py, uqff_ipc.h/cpp |
| **Python Calculators** | 50+ files | Functional | QCalc.py, CondensedPhysics.py, CondensedPhysics2.py |
| **Physics Modules** | 173+ files | C++ complete, JS complete | source1-173.cpp/.js |
| **Data Layer** | 12 files | Schema ready | IPData.py, OPData.py, shared_constants.* |
| **Grok Integration** | 8 files | Complete | GrokAPI.py, grok_url_calculators.py |
| **Build/Config** | 25+ files | Active | CMakeLists.txt, setup.py, requirements.txt |
| **Documentation** | 80+ files | Comprehensive | *.md files covering all phases |
| **Support/Utility** | 150+ files | Mixed | Tests, backups, analysis tools |
| **DLL/Libraries** | 200+ files | Compiled | VTK, Qt6, libcurl, AWS, etc. |

---

## TIER 1: GUI LAYER (Principal User Interface)

### Core GUI Files
```
source2.cpp                          (15,753 lines) - PRINCIPAL GUI, 21 tabs, Qt6
  ├── INPUT: Tab 1 - User query field (astronomy search)
  ├── PROCESSING: IPC dispatcher (MISSING - needs wiring)
  ├── OUTPUT: 
  │   ├── Tab 8 - VTK visualization
  │   └── Tab 9 - Session Logger (recall storage)
  └── RELATED FILES:
      ├── source2_mainwindow.h        - Qt main window class
      ├── source2_widgets_enhanced.h  - Custom widget components
      ├── source2_event_bus.h         - Event signaling between tabs
      ├── source2_minimal_test.cpp    - Test harness
      ├── SOURCE2_ARCHITECTURE_REPORT.md
      └── SOURCE2_ANALYSIS_AND_INTEGRATION_PLAN.md
```

### GUI Data Input
```
APIFetch.py                          (55 astronomical APIs)
  ├── CALLS: SIMBAD, NASA, VizieR, NED, Gaia, MAST, JPL, WebSockets
  ├── OUTPUT: bodies_YYYYMMDD_HHMMSS.csv (time-stamped dataset)
  ├── RELATED:
  │   ├── APIKeyManager.py            - API key storage/retrieval
  │   ├── API_KEYS_SETUP.md           - Configuration guide
  │   ├── API_KEY_MANAGEMENT_COMPLETE.md
  │   └── bodies.csv (legacy default)
  └── OBSERVATION: APIFetch can be called directly OR via source2.cpp
```

### GUI Support Libraries
```
GrokAPI.py                           (Grok xAI integration)
  ├── Features: Restore point history, codebase verification
  ├── Calls: xAI Grok 4 API for physics Q&A
  ├── RELATED:
  │   ├── SuperGrok4_RestorePoint.json - Conversation history
  │   ├── GROK_ACTIVATION_GUIDE.md
  │   ├── GROK_INTEGRATION_SUMMARY.md
  │   └── grok_api_config.json
  └── SECURITY: Authorship query interception (prevents hallucination)
```

### GUI Constants/Schemas
```
shared_constants.py                  (Shared physics constants for display)
  ├── κ, [SSq], U_UA, β_i, k_η, H_SCm values
  └── Used for: UI rendering, validation displays

InputData.py                         (Input parameter schema)
  ├── Optional fields: M, r, z, B, T, SFR, magnetic_field, etc.
  
equation_renderer.h                  (LaTeX equation formatting)
  └── USED BY: Tab 8 visualization for long-form equations
```

---

## TIER 2: PHYSICS BACKEND (Headless Engine + Symbolic/GPU Compute)

### Backend Orchestration
```
source2(HEAD PROGRAM).cpp            (4,021 lines) - BACKEND ENGINE
  ├── INFRASTRUCTURE:
  │   ├── libtorch (GPU compute)
  │   ├── SymEngine (symbolic math)
  │   ├── pybind11 (embed Python)
  │   ├── VTK (3D rendering)
  │   └── Qt3D (3D visualization)
  ├── KEY FUNCTION: PerformSearch()
  │   ├── Queries 15+ APIs (NASA, MAST, JPL, EHT, LIGO, SKA WebSockets)
  │   ├── Distributes results to 21 browser windows
  │   └── CRITICAL GAP: Doesn't call ANY physics calculator
  │       └── MISSING: Dispatch to CondensedPhysics.solve()
  ├── IPC SERVER:
  │   ├── Receives: CALCULATE_FIELD, PIPELINE_PROCESS messages
  │   ├── Handlers: (Structure defined, NOT YET IMPLEMENTED)
  │   └── Sends: RESPONSE_DATA back to source2.cpp
  └── RELATED:
      ├── VR_ARCHITECTURE_STRATEGY_Feb2026.md
      ├── CROSS_PLATFORM_C_PLUS_PLUS_ARCHITECTURE.md
      ├── FRONTEND_BACKEND_PIPELINE_ARCHITECTURE.md (THE PLAN)
      └── SOURCE2_ANALYSIS_AND_INTEGRATION_PLAN.md
```

### Physics Engine (C++)
```
MAIN_1_CoAnQi.cpp                   (18,466 lines) - INTEGRATED C++ PHYSICS
  ├── STRUCTURE:
  │   ├── 446 unique modules (SOURCE1-116 integrated)
  │   ├── 6,688+ physics terms registered
  │   ├── 121+ astronomical systems hardcoded
  │   ├── 26-layer compressed gravity framework
  │   └── SOURCE4 unified field (37 functions integrated Dec 5, 2025)
  ├── CLASSES:
  │   ├── CalculatorCore             (line 566)
  │   ├── ModuleRegistry             (line 330)
  │   ├── PhysicsTermRegistry        (line 411)
  │   ├── CrossModuleCommunicator    (line 473)
  │   ├── 63 Extracted PhysicsTerm subclasses (lines 703+)
  │   └── DynamicVacuumTerm, QuantumCouplingTerm, etc.
  ├── MENU OPTIONS (16 interactive):
  │   ├── 1. Calculate system (single)
  │   ├── 2. Calculate ALL systems (parallel, Windows threading)
  │   ├── 3-8. Clone, add, simulate, optimize
  │   ├── 9. WSTP/Wolfram integration
  │   ├── 15. SOURCE4 Unified Field Validation
  │   └── 16. Exit
  ├── BUILD STATUS:
  │   ├── CMake: Visual Studio 2022, Release-MaxCompress
  │   ├── Compression: UPX 5.0.2 (1.43 MB, 15.51% ratio)
  │   ├── Build Output: build_msvc/Release/MAIN_1_CoAnQi.exe
  │   └── Last Success: Commit 3e66d94 (Dec 5, 2025) SOURCE4 integration
  ├── VALIDATION:
  │   ├── 99.9% UQFF solvability (Grok Sept 14-21, 2025)
  │   ├── Calibrated constants: κ=0.0005, [SSq]=0.57, etc.
  │   └── 15 manuscript queries planned (never executed - Phase 1 abandoned)
  └── RELATED FILES:
      ├── MAIN_1.cpp (24,000 lines - original framework)
      ├── MAIN_1_backup_*.cpp (42 version control backups)
      ├── MAIN_1_CoAnQi_integration_status.json
      ├── ENHANCEMENT_GUIDE.md (self-expanding framework 2.0)
      ├── BUILD_INSTRUCTIONS_PERMANENT.md (CRITICAL: vcpkg warnings)
      └── SOLUTION FILE:
          └── StarMagic_UQFF.sln (Visual Studio solution 2022)
```

### Physics Modules (C++)
```
source1.cpp  ─┬─ source173.cpp       (173 original modules)
              ├── ENHANCED: source14.cpp - source162.cpp (138 modules)
              │   └── Features: Self-expanding 2.0, dynamic registration
              ├── SKIPPED: 35 files (GUI infrastructure, wrappers)
              └── INTEGRATED TO: MAIN_1_CoAnQi.cpp (446 unique modules)

● Individual module structure (example source6.cpp):
  ├── Physics equations (differential, integral, empirical)
  ├── Astrophysical system data
  ├── Self-expanding framework hooks
  └── Test harnesses (source6_*.cpp variations)

● Module categories (CANONICAL listing):
  ├── Protostellar Jets (source1-4)
  ├── Galaxy Evolution (source5-10)
  ├── Black Hole Growth (source11-15)
  ├── Magnetohydrodynamics (source16-20)
  ├── Relativistic Flows (source21-25)
  ├── Quantum Fields (source26-30)
  ├── Gravitational Waves (source31-35)
  ├── Nucleosynthesis (source36-45)
  ├── Cosmology (source46-60)
  ├── Dark Matter/Energy (source61-70)
  ├── Astroparticle (source71-80)
  ├── Nuclear Physics (source81-90)
  ├── Aether/Exotic Matter (source91-100)
  ├── Wolfram/Numeric (source101-130)
  ├── Advanced Hypergraph (source131-162)
  ├── Extended Matter (source163-173)
  └── Wolfram Bridges (source174-201) - Embedded WSTP wrappers

● RELATED:
  ├── INTEGRATION_TRACKER.csv (status of all 173 files)
  ├── COMPLETE_PHYSICS_CLASS_INVENTORY.csv (extracted terms)
  ├── physics_extraction_report.csv
  └── PHYSICS_TERM_MASTER_INVENTORY.md
```

### Physics Modules (JavaScript)
```
index.js                             (23,790 lines) - JAVASCRIPT COMPUTATIONAL ENGINE
  ├── EXPORTS: 106 astrophysical systems (require() library)
  ├── CONSTANTS: 40+ calibrated physics constants
  ├── CALCULATORS: 100+ implemented calculation methods
  ├── REFERENCE SYSTEMS:
  │   ├── SGR1745 (magnetar, 15 calculation methods)
  │   ├── Sagittarius A* (SMBH, 12 methods)
  │   ├── NGC2264 (star formation, 8 methods)
  │   ├── Tadpole Galaxy (merging, 10 methods)
  │   ├── ... (97 more systems, see GROK_URL_EQUATIONS_CATALOG.md)
  │   └── Each system fully parameterized
  ├── ARCHITECTURAL RULES:
  │   ├── NO hardcoded data (use InputData parameters)
  │   ├── Pure computation (stateless calculators)
  │   ├── Self-contained module
  │   └── Math-only, no I/O
  ├── ENTRY POINTS:
  │   ├── node index.js (direct NodeJS execution)
  │   ├── require('index.js') (UQFF server import)
  │   └── uqff_server.js (wraps as HTTP API)
  └── RELATED:
      ├── index.backup_*.js (version control, 2 backups)
      ├── simple_web.js (web server wrapper)
      ├── quick_demo.js (test harness)
      ├── test_source*.js (80+ module tests)
      └── enhanced_dynamics.js (extended calculations)
```

### Physics Modules (JavaScript Extracted)
```
source{1-120}.js                     (120+ JS module mirrors)
  ├── Generated from: C++ source*.cpp via cpp_to_qcalc_converter.py
  ├── PURPOSE: JavaScript-accessible versions of physics equations
  ├── NAMING: source1.js ↔ source1.cpp (corresponding pair)
  ├── STRUCTURE: Class methods, constants, equations
  ├── INTEGRATION: Either standalone OR via index.js require()
  ├── DISTRIBUTION: Spread across workspace root
  └── NOTE: Often lag behind C++ updates (last sync: ~Feb 2026)
```

### Core Physics Equations (Python)
```
QCalc.py                             (9,149 lines) - GENERIC PHYSICS SOLVER
  ├── PURPOSE: Parameter-driven UQFF solver (no system names)
  ├── KEY STRUCTURE:
  │   ├── UnifiedFieldSolver class (line 3457)
  │   ├── solve(params: Dict) → results (line 9031)
  │   ├── 8 master UQFF equations instantiated
  │   ├── 5,703 Wolfram-extracted terms
  │   └── Caching: max 1000 results, 3600s TTL
  ├── MANDATORY ARCHITECTURE:
  │   ├── ✓ Receives: dataset with M, r, z, B, etc.
  │   ├── ✓ Returns:{'long_form_equations': [...], 'solutions': [...]}
  │   ├── ✗ NO hardcoded system names
  │   ├── ✗ NO global instances
  │   └── ✗ NO pre-computed solutions
  ├── ASSOCIATED EXTENSIONS:
  │   ├── QCalc_core_uqff.py          - 8 master equations
  │   ├── QCalc_Wolfram_Extensions.py - 27 Wolfram functions
  │   ├── QCalc_validation.py         - Data fetcher (10+ DBs)
  │   ├── QCalc_test.py               - 80+ unit tests
  │   ├── QCalc_stat.py               - Triple point analysis
  │   ├── QCalc_Performance.py        - Caching/optimization
  │   ├── QCalc_Advanced.py           - Wormhole/GR extensions
  │   ├── QCalc_API.py                - Flask REST API (port 8443)
  │   ├── QCalc_cpp_equations.py      - 3,565 extracted C++ equations
  │   ├── QCalc_cpp_extracted.py      - 1,239+ PhysicsTerm classes
  │   └── QCalc_js_extracted.py       - JavaScript module conversions
  ├── FILES PRESENT:
  │   ├── QCalc_BACKUP_20260215_002600.py (backup)
  │   ├── QCalc_extracted.py (legacy)
  │   ├── requirements.txt (dependencies)
  │   └── QCalc_improvements_log.txt
  └── STATUS:
      ├── ✓ Fully functional standalone
      ├── ✓ Equations validated (99.9% UQFF solvability)
      ├── ✗ NOT connected from backend (missing IPC call)
      └── ✗ QCalc_API.py created but NOT accessed by Tier 2
```

### Extended Physics Equations (Python)
```
CondensedPhysics.py                  (81,626 lines) - EXTENDED PHYSICS ENGINE
  ├── PURPOSE: 176 domain-specific calculators, 111 models
  ├── STRUCTURE:
  │   ├── UnifiedFieldSolver class (line 94257)
  │   │   └── DUPLICATE of QCalc.py ← ARCHITECTURAL PROBLEM
  │   ├── 176 Calculator classes (TriadicGravity, Orb10-Orb15, etc.)
  │   ├── 111 Model classes (Galactic, AGN, BH Evolution, etc.)
  │   └── Dynamic equation generation
  ├── ARCHITECTURAL VIOLATION:
  │   ├── Lines 25000-30000: Named system classes (MagneticBubbleOrb10, etc.)
  │   ├── VIOLATES: No hardcoded named systems in calculators
  │   └── IMPACT: System data should come from APIFetch
  ├── ASSOCIATED:
  │   ├── CondensedPhysics_InputData.py
  │   │   └── Observational parameters (GW170817, SGR1745, SgrA*)
  │   ├── CondensedPhysics_OutputData.py
  │   │   └── Stores computed solutions
  │   ├── CondensedPhysics_Validation.py
  │   │   └── Validation URLs + arXiv references
  │   └── CondensedPhysics_backup_20260204_024629.py (version backup)
  ├── KEY PROBLEM:
  │   ├── 81,626 lines trying to be both:
  │   │   ├── a) Generic physics calculator
  │   │   └── b) Named-system model library
  │   └── Should split into 2 files
  └── STATUS:
      ├── ✓ Fully functional standalone
      ├── ✓ 176 calculators implemented
      ├── ✗ DUPLICATE UnifiedFieldSolver (code duplication nightmare)
      ├── ✗ Named systems violate architecture rules
      └── ✗ NOT connected from backend (missing IPC call)
```

### Extended Physics Extension (Python)
```
CondensedPhysics2.py                 (37,354 lines) - EXTENSION 1
  ├── PURPOSE: Scalable extension for new calculators
  ├── CONTAINS: 500-600 new calculator classes (Orb Analysis systems)
  ├── ARCHITECTURE: 
  │   ├── 17 new Orb Analysis calculator classes
  │   ├── Follows MANDATORY pattern (generic, no system names)
  │   └── Designed to extend parent CondensedPhysics
  ├── INTEGRATION: Needs pipeline connection via parent
  └── STATUS:
      ├── ✓ Functional standalone
      ├── ✓ Follows generic calculator pattern
      ├── ✗ NOT integrated with CondensedPhysics pipeline
      └── ✗ Never called from backend
```

### Calculation Support (Python)
```
production_pipeline.py               (Feb 13, 2026 creation)
  ├── PURPOSE: Phase 1 production system
  ├── STATUS: ✓ Complete but NEVER EXECUTED
  ├── PURPOSE: Batch runner for 15 manuscript queries
  └── COMPANION: batch_runner.py (batch execution driver)

arxiv_validation_framework.py
  ├── PURPOSE: Validate outputs against arXiv references
  └── DATA: arxiv_validation_data.csv + arxiv_validation_report.md

experimental_validation_system.py
  ├── PURPOSE: Compare UQFF vs standard predictions
  └── DATA: experimental_validation_*.json files

bose_occupancy_validation.py
  ├── PURPOSE: BEC/Bose-Einstein coherence validation
  └── DATA: holographic_sc_vs_T.csv (superconductivity curves)

bsm_physics_validation.py
  ├── PURPOSE: Beyond Standard Model physics validation
  └── INTEGRATIONS: Grok URL equations, arXiv data
```

### Extracted Equation Modules (Python)
```
grok_url_calculators.py              (121 calculator classes from Grok URL)
  ├── Source: https://x.com/i/grok/share/683542a41e744554928bfcd8b0a19e40
  ├── COUNT: 100 standard + 7 UQFF + 6 MUGE + 8 updated UQFF
  ├── STRUCTURE: Pure stateless calculators
  ├── Categories: 43 documented (protostellar → exoplanet)
  └── INVENTORY: grok_url_calculators_inventory.csv

grok_100_equations_module.py         (100+ equations)
  ├── Equations 1-100: Protostellar → Exoplanet
  ├── Classes: 100+ calculator implementations
  └── RELATED: GROK_PHYSICS_100_EQUATIONS.md (comprehensive list)

grok_100_equations_module_part2.py   (Extended equations)
  ├── Equations 64-100: Clusters → Exoplanet atmospheres
  ├── Extensions: MUGE + Electric Universe validation
  └── Calculate classes: Cluster, StarCluster, QuasarWind, NSBinary, etc.

grb_equations_module.py              (Specific domain: GRBs)
  ├── Equations 19-20: Fireball + Afterglow
  ├── Extensions: Binary pulsar GW (Eqs 73-75)
  ├── Systems: GRB_221009A, GRB_170817A, PSR J0737-3039
  └── Calculators: Fireball, Afterglow, ChirpMass, Ringdown, etc.

cosmic_domains_calculators.py
standard_astrophysics_equations.py   (100+ equation implementations)
MUGE_equations_module.py             (MUGE extended gravity)
stellar_evolution_module.py          (Star lifecycle calculations)
dark_matter_halos_module.py          (NFW, SIDM profiles)
agn_feedback_module.py               (Active galactic nuclei)
compact_objects_module.py            (NS, BH thermodynamics)
cosmology_equations_module.py        (Cosmological evolution)
mhd_dynamo_module.py                 (Magnetohydrodynamic fields)
galaxy_mergers_sfr_module.py         (Galaxy merger rates)
protostellar_jets_module.py          (Jet physics)
alpha_clustering_lenr_module.py      (Low-energy nuclear reactions)
bh_mass_scaler_module.py             (Black hole scaling relations)
bh_thermodynamics_module.py          (BH entropy, temperature)
electric_universe_gyro_module.py     (Electric universe framework)
```

### Connected Calculators/Systems
```
Various specialized calculators:
  ├── phase5_complete.py              - Phase 5 consolidation
  ├── phase6_enhanced.py              - Phase 6 extended physics
  ├── phase7_consolidated.py          - Phase 7 unified interface
  └── CoAnQi_Wrapper.py               - Wrapper around MAIN_1_CoAnQi

Pre-defined system files (data + methods):
  ├── ngc1300_uqff.js                 - NGC 1300 galaxy
  ├── ngc2264_uqff.js                 - NGC 2264 star formation
  ├── ngc346_uqff.js                  - NGC 346 cluster
  ├── ngc4676_uqff.js                 - NGC 4676 merger
  ├── redspider_uqff.js               - Red Spider nebula
  ├── smbh_msr_uqff.js                - SMBH mass-spin relations
  ├── smbhbinary_uqff.js              - Binary SMBH
  ├── ugc10214_uqff.js                - UGC 10214 tail galaxy
  ├── v838_monocerotis_uqff.js        - V838 Mon variable
  └── And 7+ more system JS files
```

---

## TIER 3: IPC LAYER (Inter-Process Communication Protocol)

### IPC Protocol Definition
```
FRONTEND_BACKEND_PIPELINE_ARCHITECTURE.md
  ├── MessageHeader (32 bytes)
  │   ├── magic: 0x55514646 (UQFF)
  │   ├── version: 1
  │   ├── type: enum MessageType (45 total)
  │   ├── payload_size: uint32_t
  │   ├── timestamp: uint64_t
  │   └── sequence: uint32_t
  ├── MessageType enum (45 types):
  │   ├── CALCULATE_FIELD (solver invocation)
  │   ├── REGISTER_TERM (dynamic physics registration)
  │   ├── PIPELINE_PROCESS (CondensedPhysics dispatcher) ← KEY
  │   ├── RESPONSE_DATA (results returned)
  │   ├── VR_FRAME_UPDATE (3D visualization)
  │   └── ... + 40 more types
  └── Payload Structures:
      ├── CalculateFieldRequest { object_name, timeout_ms, callback_id }
      ├── PipelineProcessRequest { object_name, timeout_ms, callback_id }
      └── RegisterTermRequest { term_class, parameters }
```

### IPC Channels (Implementation files)
```
ipc/ (Directory structure)
  ├── __init__.py
  ├── message.py                      - Message serialization
  ├── protocol.py                     - Protocol definition
  ├── client.py                       - IPC client (source2.cpp side)
  ├── server.py                       - IPC server (source2(HEAD).cpp side)
  ├── named_pipe_channel.py           - Windows Named Pipes (\\.\pipe\*)
  │   └── PIPE NAME: \\.\pipe\physics_service
  ├── shared_memory_channel.py        - Cross-platform SharedMemory
  │   └── SHARED MEM NAME: uqff_physics_shm
  ├── grpc_channel.py                 - gRPC optional (localhost:50051)
  ├── file_watcher.py                 - File-based sync (fallback)
  ├── file_lock.py                    - Advisory locking
  ├── serialization.py                - Data marshalling
  └── errors.py                       - Exception types
```

### IPC C++ Headers
```
uqff_ipc.h                           - C++ header for IPC
uqff_ipc.cpp                         - C++ IPC implementation
message_protocol.h                  - Message struct definitions
```

### IPC Status
```
CURRENT STATE:
  ✓ Protocol fully specified (FRONTEND_BACKEND_PIPELINE_ARCHITECTURE.md)
  ✓ Message definitions complete
  ✓ Channel implementations defined
  ✗ source2.cpp does NOT send CALCULATE_FIELD
  ✗ source2(HEAD).cpp does NOT receive or handle messages
  ✗ PIPELINE_PROCESS handler NOT implemented
  ✗ Results NOT returned via RESPONSE_DATA

MISSING WIRING:
  1. source2.cpp Tab 1: send IPC message on query submit
  2. source2(HEAD).cpp: create IPC server and listen loop
  3. Dispatcher: receive PIPELINE_PROCESS → call CondensedPhysics.solve()
  4. Backend: marshal results as RESPONSE_DATA
  5. source2.cpp: receive and display results
```

---

## DATA LAYER (Schemas & Persistence)

### Input Data Schema
```
InputData.py                         (Generic input parameters)
  ├── Fields: M (mass), r (radius), z (redshift), T (temperature)
  ├── Optional: B (magnetic field), SFR (star formation rate), etc.
  ├── PURPOSE: Consumed by ALL calculators
  └── SOURCE: APIFetch.py output

IPData.py                            (Alias for InputData)
  └── Alternative import name

CondensedPhysics_InputData.py        (Extended input for CP calculators)
  └── Pre-defined systems: GW170817, SGR1745, SgrA*
```

### Output Data Schema
```
OutputData.py                        (Generic output results)
  ├── ComputationResult dataclass
  │   ├── timestamp: datetime
  │   ├── input_params: InputData
  │   ├── solutions: Dict[str, Any]
  │   ├── available_equations: List[str]
  │   └── query_id: str
  └── PURPOSE: Created by all calculators, stored for recall

OPData.py                            (Output data persistence)
  ├── OutputDataStore class
  │   ├── store(query_id, result)
  │   ├── retrieve(query_id)
  │   └── list_queries()
  ├── PERSISTENCE: JSON files in output_data/
  └── INTEGRATION: Used by Session Logger (Tab 9)

CondensedPhysics_OutputData.py       (CP-specific output storage)
  ├── Extends OutputData for CP calculators
  └── Stores full equation solutions
```

### Shared Constants (Read-Only)
```
shared_constants.h                   (C++ calibrated constants)
  ├── UQFF parameters: κ=0.0005/day, [SSq]=0.57, etc.
  ├── Physical constants: G, c, ℏ, k_B, etc.
  └── Observational calibrations: κ, U_UA, H_SCm, etc.

shared_constants.py                  (Python mirror)
  └── Same constants, usable by Python calculators

uqff_constants.h                     (UQFF-specific C++)
  ├── F_rel, E_LEP, ω_c, k_η, etc.
  └── Framework-specific coupling constants

observational_systems_config.h       (35+ system presets)
  ├── ESO137, NGC1365, Vela, M51, Andromeda, etc.
  ├── Each: {mass, radius, spin, magnetic field, etc.}
  └── ACCESS: Optional, can be overridden by APIFetch
```

### Schema Definition Files
```
uqff_schema.json                     (JSON schema for IPC messages)
  ├── Input parameter schema
  ├── Output result schema
  └── Message payload schemas

InputData.py (as schema file)
OutputData.py (as schema file)
```

### Data Recirculation
```
RECIRCULATION PIPELINE:

bodies_YYYYMMDD_HHMMSS.csv           (APIFetch output, fresh dataset)
  ↓ INPUT
Calculators (QCalc, CondensedPhysics, etc.)
  ↓ OUTPUT
OPData.OutputDataStore / CondensedPhysics_OutputData
  ↓ PERSISTENT STORAGE
output_data/*.json (stored solutions)
  ↓ RECALL
source2.cpp Tab 9 Session Logger (display history)
  ↓ USER SELECTS PREVIOUS QUERY
bodies_YYYYMMDD_HHMMSS.csv (recirculated as input)

FILES INVOLVED:
  ├── APIFetch.py (creates bodies_*.csv)
  ├── InputData.py (DataClass for parameters)
  ├── Production calculators (consume InputData)
  ├── OutputData.py (ComputationResult)
  ├── OPData.py (OutputDataStore class)
  ├── result_cache.py (Session Logger persistence)
  └── CondensedPhysics_OutputData.py (CP-specific storage)
```

---

## CONFIGURATION & BUILD

### Build System
```
CMakeLists.txt                       (Main CMake config)
  ├── Generator: Visual Studio 17 2022 / x64
  ├── Standard: C++20
  ├── Optimization: Release-MaxCompress
  ├── Compression: UPX 5.0.2
  ├── Threading: Windows API (<windows.h>, <process.h>)
  ├── Wolfram: WSTP integration (conditional)
  ├── Dependencies: VTK, Qt6, libtorch, SymEngine, pybind11
  └── OUTPUT: build_msvc/Release/MAIN_1_CoAnQi.exe

CMakeLists_31c4c13.txt               (Alternative/backup config)
CMakePresets.json                    (Preset configurations)
CMAKE_RESTORE_POINT.txt              (Restore point marker)

StarMagic_UQFF.sln                   (Visual Studio solution)
  ├── Project 1: MAIN_1_CoAnQi.vcxproj
  ├── Project 2: source2.vcxproj (GUI)
  ├── Project 3: test_integration.vcxproj
  └── 50+ other projects (individual modules)

Individual Module Projects:
  ├── source1.vcxproj - source173.vcxproj (173 files)
  ├── test_*.vcxproj (40+ test projects)
  └── UQFFCore.vcxproj (core utility)
```

### Python Configuration
```
setup.py                             (Python package setup)
  ├── Name: star-magic-uqff
  ├── Entry points: qcalc, condensed-physics, validators
  └── Dependencies: Python 3.10+, numpy, scipy, sympy, etc.

requirements.txt                     (Python dependencies)
  ├── agent-framework-azure-ai==1.0.0b260107
  ├── agent-framework-core==1.0.0b260107
  ├── numpy, scipy, pandas, matplotlib
  ├── sympy, scikit-learn, requests
  ├── pyqt5/pyqt6 (if GUI used)
  └── 30+ total dependencies

.env.example                         (API key configuration template)
  ├── XAI_API_KEY (Grok)
  ├── FOUNDRY_PROJECT_ENDPOINT
  ├── FOUNDRY_MODEL_DEPLOYMENT_NAME
  └── Other service credentials

grok_api_config.json                 (Grok API configuration)
test_nt.json                         (Test configuration)
```

### Build Verification
```
BUILD_SUCCESS_STATUS.json            (Last successful build metadata)
BUILD_INSTRUCTIONS_PERMANENT.md      (CRITICAL: read before CMake changes)
  └── Contains vcpkg path configuration warnings
  └── MSVC 14.44.35219 version requirement
  └── C++20 standard mandatory

build_msvc/                          (Build output directory)
  ├── Release/
  │   ├── MAIN_1_CoAnQi.exe (1.43 MB compressed)
  │   ├── Qt6*.dll (7 DLL files)
  │   ├── libcurl.dll, libssl, libcrypto
  │   ├── libtorch*.dll
  │   └── 200+ dependency DLLs total
  └── Debug/ (optional)

build_backup_*/                      (Version control backups)
  └── 5+ restore points from different dates
```

---

## DOCUMENTATION (80+ files)

### Architecture Documentation
```
ARCHITECTURE_FLOW_DIAGRAM.md         (v4.0 CANONICAL, 7-component system)
FRONTEND_BACKEND_PIPELINE_ARCHITECTURE.md (3-tier IPC protocol)
CROSS_PLATFORM_C_PLUS_PLUS_ARCHITECTURE.md (Dual-platform design)

README.md                            (Project overview, UQFF theory)
Star-Magic.md                        (Complete theoretical framework)
ENHANCEMENT_GUIDE.md                 (Self-expanding 2.0 patterns)
BUILD_INSTRUCTIONS_PERMANENT.md      (Build procedure, vcpkg warnings)
```

### Physics Documentation
```
COMPLETE_UQFF_EQUATIONS_REFERENCE.md (All 8 master equations)
GROK_PHYSICS_100_EQUATIONS.md        (100+ equation catalog)
GROK_URL_EQUATIONS_CATALOG.md        (Grok URL source mappings)
UQFF_THEORY.md & UQFF_THEORY_DISPLAY.cpp (Framework theory)

Equation Modules Documentation:
  ├── Protostellar Jet Physics
  ├── Galaxy Merger Rates
  ├── Black Hole Thermodynamics
  ├── Gravitational Wave Signals
  ├── Neutron Star Properties
  ├── Cosmological Evolution
  └── Dark Matter Halo Models
```

### Integration Status
```
INTEGRATION_TRACKER.csv              (Module status: 173 files, 116 integrated)
MAIN_1_CoAnQi_integration_status.json (Complete physics inventory)
ALL_PHASES_COMPLETE_SUMMARY.md       (4 phases of enhancements)
PHASE_*.md (7 phases documented)     (Integration progress through Feb 2026)
COMPLETION_SUMMARY.md                (High-level overview)
```

### Validation & Testing
```
QCalc_test.py                        (80+ unit tests for UQFF)
test_phase*.py                       (Phases 1-7 validation)
COMPREHENSIVE_TEST_RESULTS.md        (Test run summaries)
test_*_validation.py                 (Domain-specific validation)

Validation Reports:
  ├── arxiv_validation_report.md     (arXiv reference comparison)
  ├── experimental_validation_report.md (vs. standard physics)
  ├── UQFF_VALIDATION_CONVERSATION_CAPTURE.md
  └── validation_test.js (JavaScript validation)
```

### Session Logs & Checkpoints
```
CHECKPOINT_20260219_160625.md        (Integration state snapshot)
CHECKPOINT_20260219_NAMESPACE_FIX.md (Namespace corrections)
RESTART_STATE_20260219.json          (Restart configuration)
RESTART_STATUS.md                    (Quick restart guide)

Session Captures:
  ├── conversation_A_19Feb2026.txt - conversation_F_20Feb2026_*.jsonl
  ├── copilot_thread_capture_source.cpp (single file capture)
  └── CONVERSATION_LOG.txt (running log)

Restore Points:
  └── RESTORE_POINT_20NOV2025_160720/  (timestamped save state)
  └── SuperGrok4_RestorePoint.json     (Grok conversation history)
```

### Feature/Integration Guides
```
GROK_ACTIVATION_GUIDE.md             (Grok API setup & usage)
GROK_INTEGRATION_SUMMARY.md          (Integration status)
COSMIC_EGG_INTEGRATION_GUIDE.md       (Quantum Egg feature)
WOLFRAM_INTEGRATION_COMPLETE.md       (WSTP integration)
WOLFRAM_QUICK_REFERENCE.md           (Wolfram API usage)
COANQI_USER_GUIDE.md                 (MAIN_1_CoAnQi manual)
API_KEYS_SETUP.md & API_KEY_*        (API configuration guides)
```

---

## GROK INTEGRATION (8 files)

### Grok API Integration
```
GrokAPI.py                           (Grok xAI wrapper)
  ├── call_grok_api(prompt, model, temperature)
  ├── call_grok_api_with_context(prompt, keywords, history_lookup)
  ├── load_history(timestamp) - restore point management
  ├── search_history(keyword) - query history search
  ├── get_context_from_history() - context injection
  ├── save_restore_point() - append new entries
  ├── list_restore_points() - history enumeration
  └── Codebase authorship verification (security fix)

grok_api_config.json                 (API configuration)
SuperGrok4_RestorePoint.json         (Conversation history database)
```

### Extracted Equation Modules
```
grok_url_calculators.py              (121 calculators from Grok URL)
grok_url_calculators_inventory.csv   (121-item catalog)
grok_100_equations_module.py         (100+ equation classes)
grok_100_equations_module_part2.py   (Extended + MUGE/Electric Universe)
grb_equations_module.py              (Gamma-Ray Burst specific)
GROK_PHYSICS_100_EQUATIONS.md        (Documentation)
GROK_URL_EQUATIONS_CATALOG.md        (Source mapping)
```

### Grok Configuration Guides
```
GROK_QUICK_SETUP.md                  (5-minute setup)
GROK_ACTIVATION_GUIDE.md             (Full activation procedure)
GROK_INTEGRATION_SUMMARY.md          (Status & usage)
GROK_UQFF_EQUATIONS_REFERENCE.md     (Equation reference)
grok_url_calculators_inventory.csv   (Equation-to-class mapping)
```

---

## SUPPORT & UTILITY FILES

### Extraction & Conversion Tools
```
extract_cpp_physics_terms.py         (Extract PhysicsTerm classes)
extract_physics_report.py            (Generate inventory reports)
cpp_to_qcalc_converter.py            (Convert C++ → python QCalc)
js_to_qcalc_converter.py             (Convert JS → python QCalc)
convert_legacy_to_enhanced.js        (Legacy → enhanced conversion)
ExtractionLayer.py                   (Unified extraction interface)
```

### Analysis & Audit Tools
```
intelligent_deduplication_analyzer.py (Find duplicate implementations)
analyze_condensed.py                 (Analyze CondensedPhysics structure)
analyze_models.py                    (Model inventory analysis)
analyze_docs.py                      (Documentation completeness)
cpp_file_analysis.csv                (C++ file metrics)
source_files_inventory.csv           (Complete file listing)
```

### Validation Tools
```
proof_integration.py                 (Integration correctness proof)
verify_equations.py                  (Equation validation)
verify_integrations.js               (JavaScript validation)
validate_uqff_calculators.py         (UQFF solvability check)
validate_uqff_muge.py                (UQFF vs MUGE comparison)
dual_physics_validation.log          (Validation execution log)
```

### Test Harnesses
```
test_integration.py                  (Core integration test)
test_phase*.py (7 files)             (Phase 1-7 tests)
test_source*.js (80+ files)          (Module JS tests)
test_qcalc_works.py                  (QCalc functionality)
test_grokapi_integration.py          (Grok API test)
COMPREHENSIVE_TEST_RESULTS.md        (Test execution report)
```

### PowerShell Build/Enhancement Scripts
```
build_and_run.ps1                    (Build → execute pipeline)
enhance_modules.ps1                  (Apply 2.0 framework)
enhance_to_self_expanding.ps1        (Self-expansion pattern)
comprehensive_fix.ps1                (Multi-file corrections)
fix_*.ps1 (20+ files)                (Targeted fixes)
convert_wolfram_registrations.ps1    (Wolfram conversion)
extract_*.ps1 (5+ files)             (Batch extraction scripts)
```

### Backup & Version Control
```
RESTORE_POINT_* (15+ directories)    (Timestamped backups)
*_backup_*.cpp (42 MAIN_1 backups)   (Version history)
*_backup_*.py (10+ Python backups)   (Version history)
*.bak files (various)                (Incremental backups)
.git/ (git repository)               (Full version history)
backup_older~/ (legacy backups)      (Pre-migration state)
```

### Workspace & Session Management
```
WORKSPACE_QUICK_REFERENCE.md         (Quick lookup guide)
WORKSPACE_STATUS_CURRENT.md          (Current state snapshot)
WORKSPACE_CONFIGURATION_COMPLETE.md  (Setup completion status)
SESSION_*.md (5+ files)              (Session logs, restarts)
VS_RESTART_FILES_*.txt               (VS rebuild tracking)
VSCODE_SETUP_COMPLETE.md             (VSCode configuration)
```

### CSVData Files & Exports
```
bodies_*.csv (3 timestamped)         (APIFetch outputs)
bodies.csv (default/latest)          (Current dataset)
INTEGRATION_TRACKER.csv              (Status of 173 modules)
all_physics_terms_inventory.csv      (Complete term list)
all_physics_classes.csv              (Extracted class catalog)
grok_url_calculators_inventory.csv   (121 grok calculators)
cpp_file_analysis.csv                (C++ file metrics)
triple_point_analysis.json           (Statistical analysis)

Simulation Output CSVs:
  ├── *_uqff*.csv (10+ simulation outputs)
  ├── gw150914_uqff_sim.csv           (LIGO GW150914 comparison)
  ├── uqff_*.csv (various physical quantities)
  └── holographic_sc_vs_T.csv (superconductivity)
```

### Database Files
```
uqff_equations.db                    (SQLite cache)
coanqi_cache.db                      (Cache database)
equations_*.db (various)             (Equation storage)
```

---

## WIRING DIAGRAM (Data Flow & Connections)

### Current State (What Exists)

```
┌──────────────────────────────────────────────────────────────────┐
│                     TIER 1: GUI LAYER                           │
│  source2.cpp (15,753 lines)                                      │
│  ├── 21 Tabs (query input, results display, session log)         │
│  ├── Calls: APIFetch.py → bodies_YYYYMMDD_HHMMSS.csv             │
│  ├── Displays: Tab 8 (VTK), Tab 9 (Session Logger)               │
│  └── IPC Send: MISSING ✗                                         │
└──────────────────────────────────────────────────────────────────┘
                              │
         NEEDS CONNECTION:    │ IPC CALCULATE_FIELD message
         source2.cpp →        │ (currently not sent)
         source2(HEAD).cpp    ↓
                        [MISSING IPC]
                              │
┌──────────────────────────────────────────────────────────────────┐
│                TIER 3: IPC LAYER (Protocol)                      │
│                                                                  │
│  Message Types: 45 defined                                       │
│  - CALCULATE_FIELD (send)                                        │
│  - PIPELINE_PROCESS (backend dispatcher)                         │
│  - RESPONSE_DATA (return results)                                │
│                                                                  │
│  Channels: Named Pipes, SharedMemory, gRPC, File-based           │
│  Status: ✓ Protocol specified, ✗ Code not implemented           │
└──────────────────────────────────────────────────────────────────┘
                              ↑
                        [MISSING CODE]
                              │
┌──────────────────────────────────────────────────────────────────┐
│                  TIER 2: BACKEND LAYER                           │
│  source2(HEAD PROGRAM).cpp (4,021 lines)                         │
│                                                                  │
│  ├── IPC Server: MISSING ✗ (not listening)                      │
│  ├── PerformSearch(): Queries APIs but doesn't calc ✗            │
│  └── Physics Dispatcher: MISSING ✗                              │
│      └── Should call: CondensedPhysics.solve() NOT WIRED        │
│                                                                  │
│  MAIN_1_CoAnQi.cpp (18,466 lines)                                │
│  ├── 6,688+ terms, 121+ systems                                  │
│  ├── ✓ Standalone executable (menu-driven)                       │
│  └── ✗ Not called from backend                                   │
│                                                                  │
│  index.js (23,790 lines)                                         │
│  ├── 106 system definitions                                      │
│  ├── ✓ Fully functional                                          │
│  └── ✗ Backend doesn't require() it                             │
└──────────────────────────────────────────────────────────────────┘
                              │
                              ├─→ QCalc.py (NOT CALLED) ✗
                              │   └─ UnifiedFieldSolver.solve()
                              │
                              └─→ CondensedPhysics.py (NOT CALLED) ✗
                                  ├─ 176 calculators
                                  ├─ DUPLICATE UnifiedFieldSolver
                                  └─ Should pipe to CondensedPhysics2.py
```

### Missing Connections (Priority Order)

#### PRIORITY 1: Enable Backend Calculation Dispatch

```
MISSING WIRING #1: Backend IPC Server & Handler
───────────────────────────────────────────────

File: source2(HEAD PROGRAM).cpp
Location: Lines 3800-3900 (IPC signal/slot section)

WHAT EXISTS:
  ✓ Qt signal/slot infrastructure
  ✓ QProcess available for subprocess launch
  ✓ pybind11 embedded Python ready
  
MISSING CODE:
  1. IPC server listen loop (Windows Named Pipes OR SharedMemory)
  2. Message deserialization from CALCULATE_FIELD
  3. Handler dispatch by message type:
     
     if (msg.type == PIPELINE_PROCESS) {
         // MISSING: Call CondensedPhysics.solve()
         // MISSING: Poll results from Python subprocess
         // MISSING: Package results as RESPONSE_DATA
         // MISSING: Send back via IPC
     }
  
  4. Error handling & timeouts

CONNECTION TARGET: CondensedPhysics.solve(params)
  ├── Must be callable via subprocess (Popen) OR
  ├── Must be exposed via QCalc_API.py (Flask REST)
  └── Parameters: InputData dictionary
```

#### PRIORITY 2: Enable Source2.cpp IPC Client

```
MISSING WIRING #2: GUI Sends Query via IPC
────────────────────────────────────────

File: source2.cpp
Location: Tab 1 query field (user input)

WHAT EXISTS:
  ✓ Query field with user input
  ✓ APIFetch.py called to get data
  ✓ DisplayResults() function for output
  
MISSING CODE:
  1. On query submit:
     - Get dataset from APIFetch
     - Package as InputData dictionary
     - Create CALCULATE_FIELD message
     - Send via IPC/Named Pipe to backend
  
  2. Receive RESPONSE_DATA message:
     - Deserialize results
     - Display in Tab 8 (VTK) and Tab 9 (Session Logger)
  
  3. Timeout & retry logic if backend unresponsive

REQUIRES: IPC client library (already defined in ipc/*.py)
```

#### PRIORITY 3: Python Calculator Accessibility

```
MISSING WIRING #3: Make QCalc & CondensedPhysics Callable
──────────────────────────────────────────────────────

OPTION A: REST API Endpoint
  File: QCalc_API.py (exists, already Flask, port 8443)
  Connection: source2(HEAD).cpp → HTTP POST → QCalc_API.py
  Status: Flask server code exists, just needs backend invocation

OPTION B: Subprocess (pybind11)
  Method: QProcess.start() → python CondensedPhysics.py subprocess
  Status: Infrastructure exists (pybind11, Qt6), needs dispatcher code

OPTION C: Direct Python Embed (pybind11)
  Method: Direct call via embedded Python interpreter
  Challenge: Thread safety, GIL management
  Status: Possible but complex

RECOMMENDATION: Use Option A (REST API)
  ├── Already created (QCalc_API.py)
  ├── Simplest integration
  ├── Minimal threading issues
  └── HTTP POST with JSON body
```

#### PRIORITY 4: Resolve Code Duplication

```
MISSING WIRING #4: Unify QCalc vs CondensedPhysics
───────────────────────────────────────────────

PROBLEM:
  QCalc.py (line 3457): class UnifiedFieldSolver
  CondensedPhysics.py (line 94257): class UnifiedFieldSolver (DUPLICATE)
  
  IDENTICAL IMPLEMENTATION:
  ├── 8 UQFF master equations
  ├── solve(params) method signature
  ├── Results dictionary structure
  └── Caching logic

SOLUTION:
  1. Keep QCalc.py as GENERIC physics engine
  2. CondensedPhysics.py should:
     ├── Import UnifiedFieldSolver from QCalc
     ├── Extend with domain-specific calculators (176 classes)
     ├── Remove duplicate UnifiedFieldSolver
     └── Compose with QCalc for base equations
  
  3. Backend calls: CondensedPhysics.solve() (not QCalc, to get 176 calculators)
  
  4. CondensedPhysics2.py extends CondensedPhysics (extension pattern)

FILES TO MODIFY:
  ├── QCalc.py - make UnifiedFieldSolver importable
  ├── CondensedPhysics.py - import & remove duplicate
  ├── CondensedPhysics2.py - confirm extension pattern
  └── 3 files minimum change
```

#### PRIORITY 5: Session Logger Persistence

```
MISSING WIRING #5: OPData ↔ Session Logger Connection
─────────────────────────────────────────────────

FILE: source2.cpp Tab 9
MISSING:
  ├── Retrieve old queries from OPData.OutputDataStore
  ├── Display in session history list
  ├── Recirculate selected query as new input (→ APIFetch)
  
RELATED FILES:
  ├── OPData.py (OutputDataStore class)
  ├── CondensedPhysics_OutputData.py (CP-specific storage)
  ├── result_cache.py (session storage)
  └── source2.cpp Tab 9 handler (needs implementation)
```

---

## FILE-BY-FILE DEPENDENCY MATRIX

### Source Files (C++ Physics)

| Module | File Count | Status | Dependencies | Integrated To |
|--------|-----------|--------|----------|-----------|
| source1-4 (Jets) | 4 | ✓ Enhanced | Physics equations | MAIN_1_CoAnQi |
| source5-10 (Galaxy) | 6 | ✓ Enhanced | Stellar fields | MAIN_1_CoAnQi |
| source11-15 (BH) | 5 | ✓ Enhanced | Metric tensors | MAIN_1_CoAnQi |
| ... | ... | ... | ... | ... |
| source1-173 (All) | 173 | ✓116 integrated<br>✗57 skipped | Various | MAIN_1_CoAnQi |

**Note:** See INTEGRATION_TRACKER.csv for complete 173-file status

### Python Calculator Dependencies

```
QCalc.py
  ├── IMPORTS: numpy, scipy, sympy, wolfram kernel
  ├── USED BY: QCalc_API.py (wrapper)
  │            backend dispatcher (MISSING)
  └── EXTENSION FILES:
      ├── QCalc_core_uqff.py (8 equations it uses)
      ├── QCalc_Wolfram_Extensions.py (27 functions it uses)
      ├── QCalc_test.py (tests it)
      └── 8 other QCalc_*.py (functionality)

CondensedPhysics.py
  ├── IMPORTS: QCalc (should import, doesn't)
  │            numpy, scipy
  ├── USED BY: backend dispatcher (MISSING)
  │            production_pipeline.py (never executed)
  └── EXTENSIONS:
      ├── CondensedPhysics2.py (should extend it)
      ├── CondensedPhysics_InputData.py (input schema)
      └── CondensedPhysics_OutputData.py (result storage)

Grok Modules
  ├── grok_url_calculators.py (standalone)
  ├── grok_100_equations_module.py (standalone)
  ├── grb_equations_module.py (standalone)
  └── Could be imported into CondensedPhysics
```

### Data Layer Dependencies

```
APIFetch.py
  ├── OUTPUT: bodies_YYYYMMDD_HHMMSS.csv
  └── CONSUMED BY: InputData.py (wrapping)
                   Calculators (parameter extraction)

InputData.py
  ├── CONSUMED BY: QCalc.solve(), CondensedPhysics.solve()
  ├── CREATED FROM: APIFetch output
  └── STORED AS: OPData comput ation results

OPData.py / result_cache.py
  ├── CONSUMED BY: source2.cpp Tab 9 (Session Logger)
  └── CREATED BY: All calculators
```

### IPC Dependencies

```
source2.cpp
  ├── CALLS (MISSING): ipc.client.send() CALCULATE_FIELD
  └── WAITS FOR (MISSING): ipc.client.recv() RESPONSE_DATA

source2(HEAD).cpp
  ├── LISTENS (MISSING): ipc.server.listen() CALCULATE_FIELD
  ├── CALLS (MISSING): CondensedPhysics.solve()
  └── SENDS (MISSING): ipc.server.send() RESPONSE_DATA
```

---

## CRITICAL GAPS (Must Fix for Unification)

| Gap # | Component | Problem | Impact | Priority | Est. Fix |
|-------|-----------|---------|--------|----------|----------|
| 1 | IPC Client | source2.cpp doesn't send CALCULATE_FIELD | No backend calls | **P0** | 2-4 hrs |
| 2 | IPC Server | source2(HEAD) doesn't listen/dispatch | Backend unreachable | **P0** | 3-5 hrs |
| 3 | Physics Call | Backend doesn't call CondensedPhysics | No calculations | **P0** | 1 hr |
| 4 | Code Dup | QCalc & CondensedPhysics both UnifiedFieldSolver | Maintenance nightmare | **P1** | 2 hrs |
| 5 | Session Logger | Tab 9 doesn't show/recall outputs | No history | **P1** | 1-2 hrs |
| 6 | API Exposure | QCalc_API.py created but not used | Flask server idle | **P2** | 30 min |
| 7 | Named Systems | CondensedPhysics violates generic rule | Architecture violation | **P2** | 4 hrs |
| 8 | Extension | CondensedPhysics2 not wired to CP | Extensions unavailable | **P3** | 1 hr |

---

## NEXT STEPS (Implementation Order)

### Phase: Unification (Priority 0)

```
STEP 1: IPC Server in source2(HEAD).cpp
────────────────────────────────────
File: source2(HEAD PROGRAM).cpp (add ~200 lines)
Task:
  1. Create IPC server listen loop
  2. Named Pipe channel: \\.\pipe\physics_service
  3. Message deserialization
  4. PIPELINE_PROCESS handler → call CondensedPhysics.solve()
  5. RESPONSE_DATA serialization & send back
Time: 3-5 hours
Validation: Test IPC message roundtrip

STEP 2: IPC Client in source2.cpp
──────────────────────────────
File: source2.cpp Tab 1 (add ~150 lines)
Task:
  1. On query submit: Get APIFetch result
  2. Create CALCULATE_FIELD message
  3. Connect to backend IPC server
  4. Send message, wait for RESPONSE_DATA
  5. Display results in Tab 8 & 9
Time: 2-4 hours
Validation: End-to-end query test

STEP 3: Physics Calculator Call
────────────────────────────────
File: source2(HEAD PROGRAM).cpp backend handler (add ~50 lines)
Task:
  1. Receive PIPELINE_PROCESS message
  2. Parse InputData from payload
  3. Call: CondensedPhysics.solve(input_dict)
  4. Package results as RESPONSE_DATA
Time: 1 hour
Validation: CondensedPhysics returns valid output

TOTAL FOR UNIFICATION: ~6-10 hours
Post-unification: Test full query flow end-to-end
```

### Phase: Code Quality (Priority 1)

```
STEP 4: Deduplication (QCalc vs CondensedPhysics)
─────────────────────────────────────────────
File: QCalc.py, CondensedPhysics.py
Task:
  1. Extract UnifiedFieldSolver to shared module
  2. CondensedPhysics imports from QCalc
  3. Remove duplicate implementation
  4. Verify all 176 CP calculators still work
Time: 2 hours
Validation: Test both QCalc & CondensedPhysics.solve()

STEP 5: Session Logger Implementation
──────────────────────────────────
File: source2.cpp Tab 9 handler
Task:
  1. Load OutputDataStore queries
  2. Display history list
  3. On selection: recirculate as input
Time: 1-2 hours
Validation: Retrieve and re-run previous query

TOTAL FOR QUALITY: ~3-4 hours
```

### Phase: Extensions (Priority 2)

```
STEP 6: Wire CondensedPhysics2 Extension
──────────────────────────────────────
File: CondensedPhysics.py, CondensedPhysics2.py
Task:
  1. CondensedPhysics imports from CondensedPhysics2
  2. Pipeline includes CP2 calculators
Time: 1 hour
Validation: CP2 calculator methods callable

STEP 7: API Exposure Option
────────────────────────────
File: QCalc_API.py (already exists)
Task:
  1. Optional: Use QCalc_API.py Flask REST instead of subprocess
  2. Add HTTP endpoint for backend REST calls
Time: 30 minutes - 1 hour
Validation: HTTP POST returns JSON results
```

---

## FILE LOOKUP QUICK REFERENCE

**By Component:**
- GUI: source2.cpp, source2_*.h, APIFetch.py
- Backend: source2(HEAD PROGRAM).cpp, MAIN_1_CoAnQi.cpp, index.js
- IPC: ipc/*.py, uqff_ipc.h/cpp
- Calculators: QCalc.py, CondensedPhysics.py, CondensedPhysics2.py
- Data: InputData.py, OPData.py, OutputData.py
- Constants: shared_constants.h/py, uqff_constants.h
- Build: CMakeLists.txt, setup.py, StarMagic_UQFF.sln
- Docs: ARCHITECTURE_FLOW_DIAGRAM.md, FRONTEND_BACKEND_PIPELINE_ARCHITECTURE.md

**By Function:**
- Send user query: source2.cpp Tab 1 → APIFetch.py
- Get API data: APIFetch.py → bodies_YYYYMMDD_HHMMSS.csv
- Compute physics: CondensedPhysics.py.solve()
- Store results: OPData.OutputDataStore
- Recall history: source2.cpp Tab 9 → OPData

**Most Important Files:**
1. FRONTEND_BACKEND_PIPELINE_ARCHITECTURE.md (THE PLAN)
2. source2.cpp (GUI)
3. source2(HEAD PROGRAM).cpp (Backend)
4. CondensedPhysics.py (Physics)
5. ARCHITECTURE_FLOW_DIAGRAM.md (System overview)

---

**TOTAL FILES CATALOGUED: 700+**  
**CRITICAL CONNECTIONS MISSING: 3** (Gaps 1-3)  
**ESTIMATED UNIFICATION TIME: 6-10 hours**  

*This document is your reference map. Refer back to CRITICAL GAPS section for what needs immediate wiring.*
