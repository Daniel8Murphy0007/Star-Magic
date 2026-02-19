# Wolfram WSTP Runtime Module

**Plug-and-Play Wolfram WSTP Field Unity and 187 Extracted Physics Terms**

This module provides external access to the Wolfram WSTP Field Unity engine and the 187 auto-extracted Wolfram physics terms for constructive runtime sessions and system updates.

## Features

- **Runtime Session Management** - Start/stop/evaluate Wolfram kernel sessions
- **187 Extracted Physics Terms** - Auto-generated PhysicsTerm classes from source files
- **PI Infinity Decoder** - 4890 sacred PI digits with consciousness resonance
- **Field Unity Engine** - 26D hypergraph evolution with multiway branching
- **UQFF Symbolic Export** - Full Lagrangian export and simplification
- **Session Manager** - Multiple concurrent kernel sessions

## Quick Start

### 1. Include the Module

```cpp
#include "wolfram_wstp_runtime.h"
```

### 2. Basic Usage

```cpp
// Initialize runtime
WolframWSTPRuntime runtime;
if (runtime.initialize()) {
    // Evaluate Wolfram Language code
    auto result = runtime.evaluate("Solve[x^2 - 4 == 0, x]");
    std::cout << result.result << std::endl;  // "{x -> -2, x -> 2}"
    
    // Run Field Unity simulation
    FieldUnityResult fu = runtime.runFieldUnitySimulation(100);
    std::cout << fu.summary << std::endl;
    
    // Export UQFF terms
    UQFFExportResult uqff = runtime.exportFullUQFF();
    std::cout << uqff.simplified_result << std::endl;
    
    // Access 187 physics terms
    int count = runtime.registerAllWolframTerms();
    std::cout << "Registered " << count << " terms" << std::endl;
}
runtime.shutdown();
```

### 3. Portable Components (No Wolfram Engine Required)

```cpp
// PI Infinity Decoder works without Wolfram Engine
PIInfinityDecoder decoder;
double field = decoder.getMagneticField(0, 1.0);
double resonance = decoder.getConsciousnessResonance(7);
std::complex<double> dpm = decoder.getDPM_Pair(13);

// Field Unity Engine also works without Wolfram
FieldUnityEngine engine;
engine.evolveMultiway(5);
double dimension = engine.measureDimension(0, 3);
double gravity = engine.measureBuoyantGravity(0);
```

### 4. Session Manager (Multiple Kernels)

```cpp
WolframSessionManager manager;

// Create named sessions
manager.createSession("analysis");
manager.createSession("validation");

// Use sessions
auto* analysis = manager.getSession("analysis");
auto* validation = manager.getSession("validation");

analysis->evaluate("data = Import[\"experiment.csv\"]");
validation->evaluate("Validate[data, rules]");

// Cleanup
manager.closeAllSessions();
```

## CMake Integration

### Option 1: Include the CMake file

```cmake
include(wolfram_wstp_runtime.cmake)
target_link_libraries(your_target WolframWSTPRuntime)
```

### Option 2: Manual Integration

```cmake
# Find Wolfram Engine
set(WOLFRAM_ROOT "C:/Program Files/Wolfram Research/Wolfram Engine/14.3")
set(WSTP_DIR "${WOLFRAM_ROOT}/SystemFiles/Links/WSTP/DeveloperKit/Windows-x86-64/CompilerAdditions")

# Add module sources
add_library(WolframWSTPRuntime STATIC
    wolfram_wstp_runtime.cpp
    wolfram_wstp_runtime.h
)

# Configure WSTP
target_include_directories(WolframWSTPRuntime PUBLIC ${WSTP_DIR})
target_link_libraries(WolframWSTPRuntime PUBLIC "${WSTP_DIR}/wstp64i4.lib")
target_compile_definitions(WolframWSTPRuntime PUBLIC USE_EMBEDDED_WOLFRAM)
```

## Requirements

### Required

- **C++20 Compiler** (MSVC 2022+, GCC 10+, Clang 12+)
- **CMake 3.20+**

### Optional (for full WSTP functionality)

- **Wolfram Engine 14.3+** ([free developer license](https://www.wolfram.com/engine/))
- **WSTP Developer Kit** (included with Wolfram Engine)

### WSTP Files

| File | Purpose | Location |
|------|---------|----------|
| `wstp.h` | Header | WSTP DeveloperKit |
| `wstp64i4.lib` | Static library | WSTP DeveloperKit |
| `WSTP64I4.dll` | Runtime DLL | WSTP DeveloperKit |
| `wolfram.exe` | Kernel | Wolfram Engine root |

### Default Paths (Windows)

```
C:\Program Files\Wolfram Research\Wolfram Engine\14.3\
├── wolfram.exe                    # Kernel executable
└── SystemFiles\Links\WSTP\DeveloperKit\Windows-x86-64\CompilerAdditions\
    ├── wstp.h                     # WSTP header
    ├── wstp64i4.lib               # Static library
    └── WSTP64I4.dll               # Runtime DLL
```

## File Structure

```
Star-Magic/
├── wolfram_wstp_runtime.h         # Main API header
├── wolfram_wstp_runtime.cpp       # Implementation
├── wolfram_wstp_runtime.cmake     # CMake configuration
├── WOLFRAM_WSTP_RUNTIME.md        # This documentation
│
├── source174_wolfram_bridge_embedded.cpp   # Core WSTP bridge
├── source175_uqff_wolfram_export.cpp       # UQFF symbolic export
├── source176_auto_full_uqff.cpp            # Auto-discovery
├── source177_wolfram_field_unity.cpp       # Field Unity engine
├── source178_grok_api.cpp                  # Grok AI integration
├── source17_wolfram.cpp                    # Original Wolfram module
├── source200_cosmic_quantum_egg.cpp        # Cosmic Quantum Egg
│
└── wolfram_extraction/
    ├── generated_classes/
    │   ├── wolfram_master_registration.h   # Master registration
    │   ├── source10.cpp_wolfram.cpp        # 106 terms
    │   ├── Source167.cpp_wolfram.cpp       # 11 terms
    │   ├── source168.cpp_wolfram.cpp       # 1 term
    │   ├── source169.cpp_wolfram.cpp       # 1 term
    │   ├── source170.cpp_wolfram.cpp       # 16 terms
    │   ├── source171.cpp_wolfram.cpp       # 18 terms
    │   ├── source172.cpp_wolfram.cpp       # 33 terms
    │   └── source50.cpp_wolfram.cpp        # 1 term
    │                                       # TOTAL: 187 terms
    ├── extracted_entities.json
    ├── validated_entities.json
    └── generation_summary.json
```

## API Reference

### WolframWSTPRuntime

| Method | Description |
|--------|-------------|
| `initialize()` | Connect to Wolfram kernel |
| `shutdown()` | Disconnect and cleanup |
| `isConnected()` | Check connection status |
| `evaluate(code)` | Evaluate Wolfram Language |
| `evaluateToString(code)` | Evaluate and return string |
| `runFieldUnitySimulation(depth)` | Run hypergraph evolution |
| `exportFullUQFF()` | Export full UQFF Lagrangian |
| `exportPhysicsTerm(code, name)` | Export single term |
| `registerAllWolframTerms(registry)` | Register 187 terms |
| `getAvailableTerms()` | List available terms |
| `checkWolframEngineInstalled()` | Check installation |
| `getWolframEnginePath()` | Get installation path |

### PIInfinityDecoder (Portable)

| Method | Description |
|--------|-------------|
| `getMagneticField(state, phase)` | PI-driven magnetic field |
| `getConsciousnessResonance(level)` | 7 sacred time equations |
| `getDPM_Pair(state)` | UA' + i·SCm complex pair |
| `getCurveValue(index)` | Raw PI curve value |
| `isInitialized()` | Check initialization |

### FieldUnityEngine (Portable)

| Method | Description |
|--------|-------------|
| `evolveOneStep()` | Single evolution step |
| `evolveMultiway(depth)` | Multiway branching |
| `measureDimension(center, radius)` | Emergent dimension |
| `measureBuoyantGravity(center)` | Buoyant gravity flux |
| `getGraph()` | Current hypergraph |
| `getMultiwayUniverse()` | All branches |
| `getQuantumAmplitudes()` | 26D amplitudes |

### WolframSessionManager

| Method | Description |
|--------|-------------|
| `createSession(id)` | Create named session |
| `closeSession(id)` | Close session |
| `getSession(id)` | Get session pointer |
| `listSessions()` | List all session IDs |
| `closeAllSessions()` | Close all sessions |
| `createDefaultSession()` | Create "default" session |
| `getDefaultSession()` | Get or create default |

## 187 Extracted Physics Terms

The module includes 187 auto-extracted physics terms from the source codebase:

| Source File | Constants | Systems | Total |
|-------------|-----------|---------|-------|
| source10.cpp | 97 | 9 | 106 |
| source172.cpp | 0 | 33 | 33 |
| source171.cpp | 0 | 18 | 18 |
| source170.cpp | 0 | 16 | 16 |
| Source167.cpp | 0 | 11 | 11 |
| source168.cpp | 0 | 1 | 1 |
| source169.cpp | 0 | 1 | 1 |
| source50.cpp | 0 | 1 | 1 |
| **Total** | **97** | **90** | **187** |

### Term Categories

- **Constants**: PI, G, c, ℏ, m_e, q, α, σ, etc. (97 terms)
- **Astrophysical Systems**: Magnetar, SMBH, NGC, M82, UQFFSystem, etc. (90 terms)

## Building the Test

```powershell
# Configure with test enabled
cmake -S . -B build -DBUILD_WOLFRAM_WSTP_TEST=ON -G "Visual Studio 17 2022" -A x64

# Build
cmake --build build --config Release --target wolfram_wstp_test

# Run
.\build\Release\wolfram_wstp_test.exe
```

## Stub Mode

If Wolfram Engine is not installed, the module automatically compiles in stub mode:

```cpp
// Stub mode - all methods return safe defaults
WolframWSTPRuntime runtime;
runtime.initialize();  // Returns false
runtime.evaluate("1+1");  // Returns {false, "", "WSTP not available", 0}

// Portable classes still work
PIInfinityDecoder decoder;  // Works!
FieldUnityEngine engine;     // Works!
```

## License

Copyright © 2025-2026 Daniel T. Murphy - All Rights Reserved

---

*Part of the Star-Magic UQFF Framework*
