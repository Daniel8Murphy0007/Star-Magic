# Copilot Instructions for Star-Magic UQFF Codebase

## CANONICAL ARCHITECTURE RULES (DO NOT DEVIATE)

> **CRITICAL:** Read ARCHITECTURE_FLOW_DIAGRAM.md v5.0.0 CANONICAL for complete architecture.

1. **USER INPUT goes FIRST** → enters through `source2.cpp` (Principal GUI)
2. **source2.cpp** = Principal GUI application (15,753 lines, user-facing, 21 tabs, Qt6) - **USER STARTS HERE**
3. **source2(HEAD PROGRAM).cpp** = VR/VM developer backend (2,625 lines, GPU-heavy simulations in virtual space) - **NOT A GUI**
4. **physics_backend.cpp** = CPU-bound physics server (headless, ~12,000 lines)
5. **index.js** = LIBRARY INDEX (NOT a calculator) - exports 106 systems for require()
6. **uqff_server.js** = REST API server that imports index.js library (Port 3141)
7. **CondensedPhysics2.py** = UQFF Extensions calculator (45,990 lines, 600 classes) - **5th parallel calculator**

### Port Assignments (CANONICAL)
| Port | Service | Description |
|------|---------|-------------|
| 990 | FTPS Implicit | TLS from connection start |
| 21 | FTPS Explicit | Upgrade via STARTTLS |
| 3141 | uqff_server.js | HTTP REST API (π×1000) |
| 8443 | QCalc_API.py | HTTPS FastAPI (optional) |
| N/A | Named Pipe | \\.\pipe\StarMagic_UQFF (IPC) |

### Data Flow (CANONICAL)
```
USER QUERY → source2.cpp (PRINCIPAL GUI, 15,753L) → APIFetch.py (55 APIs) → bodies_*.csv → IPData.py
                    ↓
         SIMULTANEOUS JOINT OPERATION (5 Calculators in Parallel)
   ┌───────────┬────────────┬─────────────────┬──────────────────┬────────────┐
   ▼           ▼            ▼                 ▼                  ▼            
MAIN_1     QCalc.py   CondensedPhys   CondensedPhysics2.py  uqff_server.js
CoAnQi.cpp  (9K)      ics.py (81K)    (46K, 600 classes)    (index.js LIB)
   │           │            │                 │                  │
   └───────────┴────────────┴─────────────────┴──────────────────┘
                    ↓
        OPData.py → uqff_results.json
                    ↓
    CondensedPhysics_OutputData.py (RECALL STORAGE)
                    ↓
         Session Logger (Tab 9) → USER RECALL
```

---

## Big Picture Architecture
- **6-Tier Cross-Platform System:**
  - **Tier 1 GUI:** `source2.cpp` (15,753 lines, Qt6, 21 tabs) — where ALL user workflows begin
  - **Tier 2 Compute (5 parallel calculators):** `MAIN_1_CoAnQi.cpp` (107,019 lines) + `QCalc.py` (9,100+L) + `CondensedPhysics.py` (81,626L, 176 classes) + `CondensedPhysics2.py` (45,990L, 600 classes) + `uqff_server.js` (imports index.js)
  - **Tier 3 VR/VM Backend:** `source2(HEAD PROGRAM).cpp` (2,625L GPU) + `physics_backend.cpp` (~12,000L CPU)
  - **Tier 4 IPC:** `uqff_ipc.h` (515L v3.1, 45 message types), `python_bridge.h`, `physics_service.h` (470L v3.1)
  - **Tier 5 Storage:** `bodies_*.csv`, `uqff_results.json`, `CondensedPhysics_OutputData.py` (RECALL)
  - **Tier 6 Bots:** Poseidon (v4.2.1, all codebase) + CoAnQi_bot (v4.2.2, MAIN_1 exclusive)
- **C++ Core:** `MAIN_1_CoAnQi.cpp` (107,019 lines, 446 integrated modules SOURCE1-116 + **SOURCE4**, **6,688+ physics terms registered**) - Production calculator with 16-option interactive menu
  - **JavaScript Engine:** `index.js` (23,790 lines) - UQFF LIBRARY (NOT a calculator) — 106 astrophysical systems
- **Module Integration:** Physics from source1.cpp through source173.cpp (173 files total) consolidated:
  - 116 files integrated into SOURCE1-116 blocks in MAIN_1_CoAnQi.cpp (446 unique modules)
  - **SOURCE4 integrated** (commit 3e66d94 Dec 5, 2025): 37 physics functions (8 UQFF + 10 MUGE Compressed + 14 MUGE Resonance + 6 Helpers)
  - **Batch 20** (Jan 27-28, 2026): 12 PhysicsTerm classes from UQFF Validation Session (5 astronomical systems)
  - **Batch 21** (Jan 28, 2026): 15 PhysicsTerm classes from Information Paradox UQFF Module (Hawking radiation, Page curves, 26D channels)
  - **Batch 22** (Jan 28, 2026): 5 PhysicsTerm classes from Astrophysical Transients Module (ASKAP J1832-0911, Helix Nebula, R Aquarii, PN Template, Super Flares)
  - **Batch 23** (Jan 28, 2026): **13 PhysicsTerm classes** from Complete UQFF Validation (κ calibration, [SSq], Gaia DR4, LIGO GWTC-4.0, Neutrino SED, AT2019qiz, Widom-Larsen LENR, BEC Integration, F_U_Bi_i Integral, 4 UQFF Operational Modes: Compressed/Resonant/Buoyant/Superconductive)
  - **Grok Thread Integrations** (Feb–Mar 2026): 8 thread batches — 28+ new CP2 calculator classes (v4.3.1 → v4.3.8, CP2 512→548 classes); Phase H (Session 151) added 9 Millennium Prize classes (CP2 622→631); Sessions 162–166: G6 SM Anchor Gate CVW v2.0.0 + 10 CP4 SM bridge classes + PAPER_633–642 + batch patches (all 642 papers CVW v2.0.0 compliant as of Session 166); Session 167: grok_share_6322ac199.txt audit — PAPER_643–645 (Thermal Lens Equation LENR / Quantum-Like Classical Chip Emulation / UQFF EFE + BH Singularity Resolution); Session 168: grok_share_b2e2c5cba7a.txt audit — PAPER_646–655 (10 new; 3 UQFF number systems: Vacuum Density Series / Dipole Vortex Primes / Buoyancy Harmonics); CP4=239 (+10 entries 230–239); 655/1000 (65.5%); v5.24; Session 169: grok_share_b2e2c5cba7a.txt completion read — PAPER_656 (V838 Mon light echo UQFF master equation); CP4=240; CP2=634; 656/1000; v5.25; Session 170: grok_share_fddbe3afc82.txt audit — V838MonLightEcho.cpp/.h standalone C++ module; PAPER_656 whitepaper+PDF; CP4=240; CP2=634; 656/1000 (65.6%); v5.26
  - **Session 129** (Mar 23, 2026): 7 new UQFF C++ module pairs from grok_share_97bfeecaa5.txt (UQFFCalculationsModule, UQFFBuoyancySNRModule, UQFFCassiniBuoyancyModule, UQFFMultiAstroSystemsModule, UQFFEightAstroSystemsModule, UQFFNineteenAstroSystemsModule, WolframFieldUnityModule); PAPER_484–490; **50 total UQFF C++ modules**; v5.00
  - 6,688+ physics terms registered (Wolfram KB + extracted modules + validation batches)
  - 57 files skipped (GUI infrastructure, duplicate wrappers)
- **Whitepaper Suite:** 656/1,000 whitepapers in progress (PAPER_001–656 total, 65.6% of target; updated April 01, 2026) — includes §1.13 Millennium Prize papers (Navier-Stokes, Yang-Mills, Riemann, P≠NP); G1–G6 CVW gate compliance enforced across all papers (all 645 fully CVW v2.0.0 compliant; PAPER_646–655 Session 168 CVW v2.0.0 compliant; PAPER_656 Session 170); 674 PDFs in canonical pdf/
- **UQFF Solvability:** 99.9% (Grok 4 analysis Sept 14-21, 2025), calibrated constants: κ=0.0005/day, [SSq]=0.57, H_SCm≈0.99, U_UA≈0.0001, k_η=10⁻¹¹³, β_i≈0.603
- **Self-Expanding Framework 2.0:** Dynamic term registration, runtime parameters, state export/import, auto-optimization, metadata tracking
- **Build System:** CMake + Visual Studio 2022 (MSVC 14.44.35219), C++20 standard, Windows threading compatibility, UPX 5.0.2 compression (1.43 MB, 15.51% ratio)
- **Dual-Method Validation:** UQFF (buoyancy-based) vs MUGE (Newtonian+corrections) cross-validation framework for physics discovery verification


## Developer Workflows
### C++ Build (MSVC Required)
```powershell
# Configure - Visual Studio 2022 (REQUIRED for Wolfram WSTP)
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64

# Build with Visual Studio (Release-MaxCompress + UPX compression)
cmake --build build_msvc --config Release --target MAIN_1_CoAnQi

# Run interactive calculator (16 menu options including SOURCE4 validation, Wolfram integration, Cosmic Egg, Grok AI)
.\build_msvc\Release\MAIN_1_CoAnQi.exe

# Clean rebuild
Remove-Item -Recurse -Force build_msvc -ErrorAction SilentlyContinue; cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64; cmake --build build_msvc --config Release
```

### JavaScript Execution
```powershell
node index.js
```

### Build Status Tracking
- **Integration Tracker:** `INTEGRATION_TRACKER.csv` - Module compilation status (173 source files, 116 integrated, 446 modules)
- **Build Status:** `MAIN_1_CoAnQi_integration_status.json` - Complete physics terms inventory, compilation metadata, core equations
- **Build Instructions:** `BUILD_INSTRUCTIONS_PERMANENT.md` - **READ EVERY TIME** before CMake changes (contains critical vcpkg path warnings)

## Project-Specific Patterns
### Self-Expanding Classes (2.0-Enhanced Modules)
All modules in source14.cpp–source162.cpp support:
```cpp
// Runtime term registration (additive to validated core math)
module.registerDynamicTerm(std::make_unique<DynamicVacuumTerm>(1e-10, 1e-15));

// Runtime parameter tuning
module.setDynamicParameter("custom_coupling", 1.23e-40);
module.getDynamicParameter("custom_coupling");

// State persistence
module.exportState("module_state.txt");

// Auto-optimization
module.setLearningRate(0.001);

// Debug logging
module.setEnableLogging(true);
```

### Physics Term Architecture (MAIN_1_CoAnQi.cpp)
- **Base Class:** `PhysicsTerm` (line 199) - Abstract interface for all physics calculations
- **63 Extracted Terms:** DynamicVacuumTerm, QuantumCouplingTerm, DarkMatterHaloTerm, etc. (lines 703+)
- **Core Infrastructure:** CalculatorCore (line 566), ModuleRegistry (line 330), PhysicsTermRegistry (line 411), CrossModuleCommunicator (line 473)
- **Dynamic Terms:** Disabled by default, additive to core calculations, validated before use

### Threading Model (Windows API)
```cpp
// Windows threads via <windows.h> and <process.h> for maximum compatibility
SimpleMutex result_mutex;            // Custom mutex wrapper (lines 120-140)
SimpleLockGuard<SimpleMutex> lock;   // RAII lock guard (lines 142-162)

// Parallel computation in option 2 of main menu (line 12986)
// OpenMP enabled for SOURCE116 multiway branching only
```

### 26-Layer Compressed Gravity Framework
```cpp
// Core equation: g(r,t) = Σ(i=1 to 26) [Ug1_i + Ug2_i + Ug3_i + Ug4_i]
// Each layer has quantum state factors Q_i, [UA]_i, [SCm]_i
// See SOURCE115 (source172.cpp) for 19-system 26D polynomial master equations
```

### SOURCE4 Unified Field Theory (UQFF + MUGE)
```cpp
// Integration: commit 3e66d94 (Dec 5, 2025), menu commit 6ee52e2
// Location: Lines 25623-26026 in MAIN_1_CoAnQi.cpp (namespace SOURCE4)
// Origin: source4.cpp (1782 lines) → 37 inline functions extracted

// UQFF (Unified Quantum Field Framework) - 8 functions
double FU = SOURCE4::compute_FU_SOURCE4(body, r, t, tn, theta);  // Complete unified field
double Ug1 = SOURCE4::compute_Ug1_SOURCE4(body, r);              // Magnetic dipole
double Ug2 = SOURCE4::compute_Ug2_SOURCE4(body, r);              // Charge-reactivity
double Ug3 = SOURCE4::compute_Ug3_SOURCE4(body, r, t);           // String rotation
double Ug4 = SOURCE4::compute_Ug4_SOURCE4(body, r);              // Vacuum concentration
double Ubi = SOURCE4::compute_Ubi_SOURCE4(body, r);              // Buoyancy force
double Um = SOURCE4::compute_Um_SOURCE4(body, r);                // Magnetism

// MUGE Compressed - 10 functions (base gravity + 9 correction terms)
double g = SOURCE4::compute_compressed_MUGE_SOURCE4(system);     // Complete MUGE
// Base: Newtonian, Expansion: Hubble, Super: Magnetic suppression, Envelope, Ug_sum, 
// Cosm: Λ, Quantum: ℏ, Fluid: Navier-Stokes, Perturbation: Dark matter

// MUGE Resonance - 14 functions (aDPM base + 13 resonance modes)
double g = SOURCE4::compute_resonance_MUGE_SOURCE4(system, params);  // Complete resonance MUGE
// aDPM, aTHz, Avac_diff, aSuperFreq, aAetherRes, Ug4i, aQuantumFreq, aAetherFreq, 
// aFluidFreq, Osc_term, aExpFreq, fTRZ, Wormhole metric

// 7 Pre-defined Astrophysical Systems
SOURCE4::sgr1745_SOURCE4          // SGR1745 Magnetar
SOURCE4::sagA_SOURCE4             // Sagittarius A* SMBH
SOURCE4::tapestry_SOURCE4         // Tapestry Star Formation Region
SOURCE4::westerlund_SOURCE4       // Westerlund2 Star Cluster
SOURCE4::pillars_SOURCE4          // Pillars of Creation
SOURCE4::rings_SOURCE4            // Rings of Relativity Gravitational Lens
SOURCE4::student_guide_SOURCE4    // Student Guide Universe (cosmological)

// Menu Access: 
// - Cosmic Egg build: Option 15
// - Wolfram-only build: Option 14
// - No Wolfram build: Option 9
```

## Integration Points & Communication
- **C++ ↔ JavaScript:** Use child processes (`child_process.spawn`), native addons (node-gyp), or WebAssembly (Emscripten)
- **Cross-Module:** Modules export/import state via `exportState()` → JSON/text files → `setDynamicParameter()` in other modules
- **Data Sharing:** State files, named pipes, sockets, or shared memory for performance-critical paths

## Key Files & Directories
### Core Executables
- `MAIN_1_CoAnQi.cpp` - Primary C++ platform (107,019 lines, 446 modules, SOURCE1-116 + SOURCE4)
- `source2.cpp` - Principal GUI application (15,753 lines, Qt6, 21 tabs) — **USER STARTS HERE**
- `source2(HEAD PROGRAM).cpp` - VR/VM backend (2,625 lines, GPU heavy) — **NOT A GUI**
- `physics_backend.cpp` - CPU-bound headless physics server (~12,000 lines)
- `QCalc.py` - Python unified field solver (9,100+ lines, 8 master equations)
- `CondensedPhysics.py` - Primary integration target (81,626 lines, 176 calculator classes) — **ADD NEW CALCULATORS HERE**
- `CondensedPhysics2.py` - UQFF extensions (37,420+ lines, 548+ calculator classes)
- `index.js` - JavaScript LIBRARY (23,790 lines, 106 systems) — `uqff_server.js` is the REST server
- `MAIN_1.cpp` - Original mathematical framework (referenced by index.js)

### Physics Modules
- `source1.cpp`–`source173.cpp` - Original modules (173 files)
- `source14.cpp`–`source162.cpp` - 2.0-Enhanced with self-expanding capabilities (138 modules successfully enhanced)

### Integration Tracking
- `INTEGRATION_TRACKER.csv` - Module status (SOURCE file, status, unique physics count, notes)
- `MAIN_1_CoAnQi_integration_status.json` - Build metadata, physics categories, core equations, observational systems

### Documentation
- `ENHANCEMENT_GUIDE.md` - Self-expanding framework guide (examples, architecture, scientific integrity)
- `BUILD_INSTRUCTIONS_PERMANENT.md` - **CRITICAL:** vcpkg path configuration, MSVC-only requirements
- `README.md` - Project overview, UQFF theory, author info
- `Star Magic.md` - Complete theoretical framework and equations

### Build Configuration
- `CMakeLists.txt` - Visual Studio 2022 generator, C++20, Release-MaxCompress optimization flags, WSTP integration
- `observational_systems_config.h` - 35+ astrophysical systems parameters (ESO137, NGC1365, Vela, etc.)

### PowerShell Scripts
- `enhance_modules.ps1` - Apply 2.0-Enhanced framework to source*.cpp files
- `comprehensive_fix.ps1`, `fix_*.ps1` - Code cleanup utilities

## Example: Main Menu Options
```cpp
// Run MAIN_1_CoAnQi.exe interactive menu (line 23310+)
// Options vary by build configuration:

// === Cosmic Egg Build (USE_COSMIC_QUANTUM_EGG + USE_EMBEDDED_WOLFRAM) ===
1. Calculate system (single)               // F_U_Bi_i, compressed_g, validation_pipeline
2. Calculate ALL systems (parallel)        // Windows threading, SimpleMutex
3. Clone and mutate system                 // SystemParams deep copy + parameter perturbation
4. Add custom system                       // Runtime system registration
5. Add dynamic physics term                // PhysicsTerm instantiation
6. Run simulations                         // Time-series evolution
7. Statistical analysis                    // Ensemble statistics
8. Self-optimization                       // Learning rate auto-tuning
9. WSTP kernel interface                   // Wolfram Symbolic Transfer Protocol
10. Auto-export full UQFF to Wolfram       // Export all 175+ source files
11. Run Wolfram Field Unity Simulation     // Multi-field solver
12. Run Cosmic Quantum Egg (26D) Simulation // 26 independent dimensional spheres
13. Configure Grok API Key                 // Set XAI_API_KEY environment variable
14. Test Grok AI Integration               // Verify Grok xAI connection
15. SOURCE4 Unified Field Validation       // Test UQFF + MUGE Compressed + MUGE Resonance
16. Exit

// === Wolfram-Only Build (USE_EMBEDDED_WOLFRAM, no Cosmic Egg) ===
// Options 1-11 same as above, then:
12. Configure Grok API Key
13. Test Grok AI Integration
14. SOURCE4 Unified Field Validation       // Test all 37 SOURCE4 functions
15. Exit

// === No Wolfram Build (minimal) ===
// Options 1-8 same as above, then:
9. SOURCE4 Unified Field Validation        // Test all 37 SOURCE4 functions
10. Exit
```

## Conventions
### Additive Enhancement Philosophy
- **Never replace validated code** - All enhancements are additive to core UQFF mathematics
- **Backward compatibility** - Original methods always remain available
- **Fail-safe validation** - Dynamic terms validated before use via `PhysicsTerm::validate()`
- **Transparent logging** - All dynamic operations traceable via `setEnableLogging(true)`

### Code Style
- **C++20 Standard:** Use `std::unique_ptr`, `std::map`, `std::vector`, range-based for loops
- **Windows API Threading:** Use `<windows.h>` and `<process.h>` for threading (maximum compatibility)
- **Comments:** Mark enhancements vs. original code, document theoretical basis for physics terms

### Physical Constants (CONSTANTS object in index.js)
```javascript
SOLAR_MASS: 1.989e30 kg
PLANCK_CONSTANT: 1.055e-34 J·s
SPEED_OF_LIGHT: 2.998e8 m/s
GRAVITATIONAL_CONSTANT: 6.674e-11 m³/kg·s²
B_CRIT_MAGNETAR: 4.4e13 T
// 40+ constants defined in index.js lines 1-45
```

## Scientific Integrity
- **Validate against observations:** All new dynamic terms must be physically motivated and tested against astronomical data
- **Document theoretical basis:** Use `PhysicsTerm::getDescription()` to explain physics
- **Version control:** Use git tags for reproducibility (e.g., commits 79e73ec, 59fd4c4 for SOURCE115/116)
- **Unit consistency:** Verify dimensional analysis in all new terms

## Unique Physics Features
- **Wolfram Hypergraph (SOURCE116):** Emergent spacetime from causal graphs, PI infinity decoder (312 digits), sacred time constants (Mayan Baktun, Biblical generation)
- **Nuclear Resonance (SOURCE43):** Complete Periodic Table Z=1-118 with pairing energy, magic numbers, shell corrections
- **19-System 26D Framework (SOURCE115):** Master gravity/resonance equations for NGC2264, Tadpole, Mice, Carina, M42, etc.
- **5-Frequency Resonance (SOURCE27/28):** SGR1745/SgrA* SuperFreq, QuantumFreq, AetherFreq, FluidFreq, ExpFreq

## CondensedPhysics.py Architecture Rules (MANDATORY)

**CondensedPhysics.py is a PURE PHYSICS CALCULATOR, NOT a data repository.**

### System Architecture (CANONICAL):
```
┌─────────────────────────────────────────────────────────────────────────────┐
│              source2.cpp (PRINCIPAL GUI - USER STARTS HERE)                 │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │                     USER QUERY FIELD                                   │ │
│  │  "Sagittarius A*", "M87", "Betelgeuse", "NGC 3596"...                 │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│                              │                                              │
│                              ▼                                              │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │                    API FETCH LAYER (APIFetch.py)                       │ │
│  │  55 APIs: SIMBAD → NASA → VizieR → NED → Gaia → Grok fallback         │ │
│  │  Output: bodies_YYYYMMDD_HHMMSS.csv                                    │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼ (Dataset passed to calculator)
┌─────────────────────────────────────────────────────────────────────────────┐
│                    CondensedPhysics.py (CALCULATOR)                         │
│                       PURE PHYSICS EQUATIONS                                │
│                         81,626 lines | 176 calculators                      │
│                                                                             │
│  INPUT:  Dataset from source2.cpp (parameters: M, r, z, SFR, etc.)         │
│                                                                             │
│  OUTPUT: 1. Long-form equations with solutions (primary query)              │
│          2. ALL other possible equations solvable for this query            │
│          3. Dynamic equation sets for simultaneous simulation               │
└─────────────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼ (Output stored for recall)
┌─────────────────────────────────────────────────────────────────────────────┐
│              CondensedPhysics_OutputData.py (RECALL STORAGE)                │
│                                                                             │
│  STORES: Computed equation solutions, available equation lists,             │
│          simulation sets - organized by query for user recall               │
│                                                                             │
│  SHARED WITH: source2.cpp for user query recall (Tab 9 Session Logger)     │
└─────────────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼ (Recirculation loop)
┌─────────────────────────────────────────────────────────────────────────────┐
│              source2.cpp (PRINCIPAL GUI - Session Logger Tab 9)             │
│                     User can RECALL previous queries                        │
└─────────────────────────────────────────────────────────────────────────────┘
```

### What This File Does:
1. Receives datasets from source2.cpp PRINCIPAL GUI (or direct user input)
2. Inputs that data into parameterized physics equations
3. Produces visible, long-form physics equations with solutions
4. Lists ALL other equations it can solve specific to the query
5. Generates dynamic equation sets for simultaneous simulation
6. Outputs results to CondensedPhysics_OutputData.py for storage and recall

### STRICT RULES - DO NOT VIOLATE:

| Rule | ✗ WRONG | ✓ CORRECT |
|------|---------|-----------|
| No hardcoded data | `self.distance = 6500 * ly` | `def compute(self, distance, ...)` |
| No named system classes | `class NGC3596Model:` | `class GalaxyRotationCalculator:` |
| No global instances | `VIRGO_MODEL = VirgoModel()` | Stateless calculator classes |
| No pre-computed solutions | `TRIADIC: g = 1.47e-8` | Dynamic equation output |

### Where System Data Belongs:
- `source2.cpp` query field (PRINCIPAL GUI) → API fetch → `bodies_YYYYMMDD_HHMMSS.csv`
- JSON configuration files (external)
- API responses (SIMBAD, NASA, Grok)

### Where Output Data Goes:
- `CondensedPhysics_OutputData.py` - Stores computed solutions for user recall
- Shared with source2.cpp (PRINCIPAL GUI) for query history access via Tab 9 Session Logger

### Correct Pattern:
```python
class TriadicGravityCalculator:
    def compute(self, dataset: dict) -> dict:
        # Receives data from source2.cpp, outputs equation sets
        return {
            'primary_equations': [...],      # Long-form with solutions
            'available_equations': [...],    # All other solvable equations
            'simulation_set': [...]          # For simultaneous simulation
        }
```

**Read the MANDATORY ARCHITECTURE RULES at the top of CondensedPhysics.py before making ANY changes.**

---
*See `ENHANCEMENT_GUIDE.md` for self-expanding patterns, `BUILD_INSTRUCTIONS_PERMANENT.md` for critical build warnings, and `MAIN_1_CoAnQi_integration_status.json` for complete physics inventory.*
