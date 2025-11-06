# Copilot Instructions for Star-Magic UQFF Codebase

## Big Picture Architecture
- **Core Engine:** `index.js` orchestrates advanced Unified Quantum Field Force (UQFF) calculations, integrating mathematical frameworks from `MAIN_1.cpp` and enhanced C++ modules (`Source14.cpp` through `Source162.cpp`).
- **Module Ecosystem:** Over 138 C++ modules (see `ENHANCEMENT_GUIDE.md`) implement self-expanding physics terms, supporting organic code growth and dynamic runtime extension.
- **Data Flow:** JavaScript (`index.js`) calls C++ modules via child processes, native addons, or WebAssembly. Modules communicate via shared state files, sockets, or memory.

## Developer Workflows
- **Build/Run C++ Modules:**
  - Compile: `g++ -std=c++17 Source14.cpp -o test_source14`
  - Test: `./test_source14`
  - Batch test: `g++ -std=c++17 test_enhanced_modules.cpp -o test_enhanced; ./test_enhanced`
- **Enhancement Script:** Use `enhance_modules.ps1` to upgrade modules with dynamic term support. Backups are stored in `module_backups_YYYYMMDD_HHMMSS/`.
- **Integration:** Use Node.js (`index.js`) to orchestrate and connect modules. For new physics terms, follow the `PhysicsTerm` interface (see `ENHANCEMENT_GUIDE.md`).

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
- `index.js`: Central computational engine
- `Source14.cpp`–`Source162.cpp`: Enhanced physics modules
- `ENHANCEMENT_GUIDE.md`: Full enhancement details and usage patterns
- `MAIN_1.cpp`: Mathematical backbone
- `enhance_modules.ps1`: Enhancement automation script
- `SETUP.md`: Setup and configuration guide

## Example: Registering a Dynamic Term
```cpp
module.registerDynamicTerm(std::make_unique<DarkMatterHaloTerm>(1e12 * M_sun, 20000));
module.setDynamicParameter("custom_coupling", 1.23e-40);
module.exportState("moduleA_shared.txt");
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
