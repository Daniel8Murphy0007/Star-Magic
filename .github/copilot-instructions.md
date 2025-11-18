# Copilot Instructions for Star-Magic UQFF Codebase

## Big Picture Architecture
- **Primary Platform:** `MAIN_1_CoAnQi.cpp` (18,466 lines) - Conscious Quantum Intelligence UQFF Calculator with 446 integrated physics modules (SOURCE1-116)
- **Module Integration:** All unique physics from source1.cpp through source173.cpp extracted and compiled into single comprehensive framework
- **Self-Expanding System:** 2.0-Enhanced framework with dynamic term registration, runtime parameters, state export/import, and auto-optimization
- **Build System:** CMake + MinGW-w64 GCC 14.2.0, C++17 standard, builds MAIN_1_CoAnQi.exe as primary executable
- **Integration Tracking:** `INTEGRATION_TRACKER.csv` documents all 161 source files, 445 modules, 359+ unique physics terms

## Developer Workflows
- **Build Primary Executable:**
  - Configure: `cmake -S . -B build -G "MinGW Makefiles"`
  - Build: `cmake --build build --target MAIN_1_CoAnQi`
  - Run: `.\build\MAIN_1_CoAnQi.exe` (8-option interactive menu)
- **Quick Build Commands:**
  - Clean: `Remove-Item -Recurse -Force build -ErrorAction SilentlyContinue`
  - Full rebuild: `cmake -S . -B build -G "MinGW Makefiles"; cmake --build build`
- **Integration Status:** Check `INTEGRATION_TRACKER.csv` for module compilation status, physics terms count, and notes

## Project-Specific Patterns
- **Self-Expanding Classes:** All enhanced modules support:
  - `registerDynamicTerm()` to add new physics terms at runtime
  - `setDynamicParameter()`/`getDynamicParameter()` for runtime parameter tuning
  - `exportState()` for state persistence and cross-module communication
  - `setLearningRate()` for auto-optimization
- **Dynamic Terms:** Additive to core calculations; original validated math is always preserved. Dynamic terms are opt-in and disabled by default.
- **Metadata:** Each module tracks metadata for provenance and versioning.
- **Logging:** Enable debug output via `setEnableLogging(true)`.

## Integration Points & Communication
- **JavaScript ↔ C++:** Use child processes, native addons, or WebAssembly for integration. Share data via exported state files or IPC.
- **Cross-Module:** Modules can export/import parameters for collaborative computation (see examples in `ENHANCEMENT_GUIDE.md`).

## Key Files & Directories
- `MAIN_1_CoAnQi.cpp`: Primary platform (18,463 lines, 446 integrated modules SOURCE1-116)
- `INTEGRATION_TRACKER.csv`: Complete module status tracking (173 source files, 116 integrated, 57 skipped)
- `MAIN_1_CoAnQi_integration_status.json`: Build status and physics terms inventory
- `source1.cpp`–`source173.cpp`: Original physics modules extracted and integrated
- `ENHANCEMENT_GUIDE.md`: Self-expanding framework documentation
- `BUILD_INSTRUCTIONS_PERMANENT.md`: Critical build workflow (READ EVERY TIME)
- `CMakeLists.txt`: Build configuration (MinGW Makefiles, C++17)
- `copilot_thread_capture_source.cpp`: Session recovery file (Nov 14 @ 5:03 PM)
- `restore_point_thread_capture_15nov2025_0248.txt`: Session recovery (Nov 15 @ 2:48 AM)

## Example: Using MAIN_1_CoAnQi Interactive Menu
```cpp
// Run the primary executable
.\build\MAIN_1_CoAnQi.exe

// 8-option menu:
// 1. Calculate system (single)
// 2. Calculate ALL systems (parallel)
// 3. Clone and mutate system
// 4. Add custom system
// 5. Add dynamic physics term
// 6. Run simulations
// 7. Statistical analysis
// 8. Self-optimization
// 9. Exit
```

## Conventions
- **Additive Enhancement:** Never replace validated code; new features are always additive.
- **Backward Compatibility:** Original methods remain available.
- **Fail-Safe:** All dynamic terms are validated before use.
- **Transparent Logging:** All dynamic operations can be traced.

## Scientific Integrity
- Validate new dynamic terms against observational data before deployment.
- Document theoretical basis in term descriptions.
- Use version control for reproducibility.

---
*For more details, see `ENHANCEMENT_GUIDE.md` and comments in `index.js` and module files.*
