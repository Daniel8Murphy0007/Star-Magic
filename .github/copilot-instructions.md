# Copilot Instructions for Star-Magic UQFF Codebase

## Big Picture Architecture
- **Dual-Platform System:**
  - **C++ Core:** `MAIN_1_CoAnQi.cpp` (19,356 lines, 446 integrated physics modules SOURCE1-116, **100 astronomical systems**) - Production calculator with 9-option interactive menu
  - **JavaScript Engine:** `index.js` (23,790 lines) - UQFF computational orchestration layer with 106 astrophysical systems
- **Module Integration:** Physics from source1.cpp through source173.cpp (173 files total) consolidated:
  - 116 files integrated into SOURCE1-116 blocks in MAIN_1_CoAnQi.cpp (446 unique modules, 446 physics terms)
  - 57 files skipped (GUI infrastructure, duplicate wrappers)
- **Self-Expanding Framework 2.0:** Dynamic term registration, runtime parameters, state export/import, auto-optimization, metadata tracking
- **Build System:** CMake + MinGW-w64 GCC 14.2.0, C++17 standard, Windows threading compatibility

## Developer Workflows
### C++ Build (Primary)
```powershell
# Configure - Visual Studio 2022 (Release-MaxCompress optimizations)
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64

# Configure - MinGW (alternative, smaller footprint)
cmake -S . -B build -G "MinGW Makefiles"

# Build with Visual Studio (Release-MaxCompress + UPX compression)
cmake --build build_msvc --config Release --target MAIN_1_CoAnQi

# Build with MinGW
cmake --build build --target MAIN_1_CoAnQi

# Run interactive calculator (9 menu options)
.\build_msvc\Release\MAIN_1_CoAnQi.exe   # Visual Studio optimized
.\build\MAIN_1_CoAnQi.exe                # MinGW

# Clean rebuild (Visual Studio)
Remove-Item -Recurse -Force build_msvc -ErrorAction SilentlyContinue; cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64; cmake --build build_msvc --config Release

# Clean rebuild (MinGW)
Remove-Item -Recurse -Force build -ErrorAction SilentlyContinue; cmake -S . -B build -G "MinGW Makefiles"; cmake --build build
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

### Threading Model (MinGW Compatibility)
```cpp
// Windows threads via <windows.h> and <process.h> (NOT std::thread)
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

## Integration Points & Communication
- **C++ ↔ JavaScript:** Use child processes (`child_process.spawn`), native addons (node-gyp), or WebAssembly (Emscripten)
- **Cross-Module:** Modules export/import state via `exportState()` → JSON/text files → `setDynamicParameter()` in other modules
- **Data Sharing:** State files, named pipes, sockets, or shared memory for performance-critical paths

## Key Files & Directories
### Core Executables
- `MAIN_1_CoAnQi.cpp` - Primary C++ platform (18,466 lines, 446 modules, SOURCE1-116)
- `MAIN_1.cpp` - Original mathematical framework (referenced by index.js)
- `index.js` - JavaScript computational engine (23,790 lines, 106 systems)

### Physics Modules
- `source1.cpp`–`source173.cpp` - Original modules (173 files)
- `source14.cpp`–`source162.cpp` - 2.0-Enhanced with self-expanding capabilities (138 modules successfully enhanced)

### Integration Tracking
- `INTEGRATION_TRACKER.csv` - Module status (SOURCE file, status, unique physics count, notes)
- `MAIN_1_CoAnQi_integration_status.json` - Build metadata, physics categories, core equations, observational systems

### Documentation
- `ENHANCEMENT_GUIDE.md` - Self-expanding framework guide (examples, architecture, scientific integrity)
- `BUILD_INSTRUCTIONS_PERMANENT.md` - **CRITICAL:** vcpkg path warnings, Visual Studio vs. MinGW conflicts
- `README.md` - Project overview, UQFF theory, author info
- `Star Magic.md` - Complete theoretical framework and equations

### Build Configuration
- `CMakeLists.txt` - Visual Studio 2022 + MinGW generators, C++17, Release-MaxCompress optimization flags, WSTP integration
- `observational_systems_config.h` - 35+ astrophysical systems parameters (ESO137, NGC1365, Vela, etc.)
- `Star Magic.md` - Complete theoretical framework and equations

### Build Configuration
- `CMakeLists.txt` - MinGW Makefiles generator, C++17, 155+ source*.cpp targets
- `observational_systems_config.h` - 35+ astrophysical systems parameters (ESO137, NGC1365, Vela, etc.)

### PowerShell Scripts
- `enhance_modules.ps1` - Apply 2.0-Enhanced framework to source*.cpp files
- `comprehensive_fix.ps1`, `fix_*.ps1` - Code cleanup utilities

## Example: Main Menu Options
```cpp
// Run MAIN_1_CoAnQi.exe interactive menu (line 12901-13050+)
1. Calculate system (single)         // F_U_Bi_i, compressed_g, validation_pipeline
2. Calculate ALL systems (parallel)  // Windows threading, SimpleMutex
3. Clone and mutate system          // SystemParams deep copy + parameter perturbation
4. Add custom system                // Runtime system registration
5. Add dynamic physics term         // PhysicsTerm instantiation
6. Run simulations                  // Time-series evolution
7. Statistical analysis             // Ensemble statistics
8. Self-optimization                // Learning rate auto-tuning
9. Exit
```

## Conventions
### Additive Enhancement Philosophy
- **Never replace validated code** - All enhancements are additive to core UQFF mathematics
- **Backward compatibility** - Original methods always remain available
- **Fail-safe validation** - Dynamic terms validated before use via `PhysicsTerm::validate()`
- **Transparent logging** - All dynamic operations traceable via `setEnableLogging(true)`

### Code Style
- **C++17 Standard:** Use `std::unique_ptr`, `std::map`, `std::vector`, range-based for loops
- **Windows Compatibility:** Use `<windows.h>` and `<process.h>` for threading (NOT `<thread>` with MinGW)
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

---
*See `ENHANCEMENT_GUIDE.md` for self-expanding patterns, `BUILD_INSTRUCTIONS_PERMANENT.md` for critical build warnings, and `MAIN_1_CoAnQi_integration_status.json` for complete physics inventory.*
