# Build Requirements for calculator_advanced

**Module**: Advanced Calculator with UQFF Integration  
**Location**: calculator_advanced/  
**Integration**: Tab 22 in source2.cpp Principal GUI  
**Status**: Code complete (6,300+ lines), build validation pending

---

## Dependencies

### Required Packages

| Package | Version | Purpose | vcpkg Install Command |
|---------|---------|---------|----------------------|
| **Qt6** | 6.10.0+ | GUI framework (Core, Widgets) | `qt:x64-windows` |
| **ANTLR4** | 4.13.2 | Parser runtime for Math.g4 grammar | `antlr4:x64-windows` ✅ |
| **SymEngine** | 0.14.0 | Symbolic mathematics backend | `symengine:x64-windows` |
| **GSL** | 2.8+ | GNU Scientific Library (polynomials ≤ degree 26) | `gsl:x64-windows` |
| **Eigen3** | 3.4.1 | Linear algebra library | `eigen3:x64-windows` ✅ |
| **QCustomPlot** | 2.1.1 | Interactive 2D plotting widget | `qcustomplot:x64-windows` |

✅ = Already installed  
🔲 = Needs installation

---

## Installation via vcpkg

### Prerequisites
- **vcpkg** installed at `C:\vcpkg\`
- **Visual Studio 2022** (MSVC 17.2022+)
- **7zip** version 26.0.0+ (required by vcpkg for extraction)

### Install Missing Dependencies

```powershell
# Update vcpkg to latest
cd C:\vcpkg
git pull
.\bootstrap-vcpkg.bat

# Install required packages
C:\vcpkg\vcpkg.exe install qt:x64-windows symengine:x64-windows gsl:x64-windows qcustomplot:x64-windows --recurse
```

**Note**: Qt6 installation takes 20-40 minutes and requires ~3 GB disk space.

---

## CMake Configuration

### Configure with vcpkg Toolchain

```powershell
cd C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic

cmake -S . -B build_msvc `
  -G "Visual Studio 17 2022" `
  -A x64 `
  -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake
```

### Build calculator_advanced Library

```powershell
# Build library only
cmake --build build_msvc --config Release --target calculator_advanced

# Build Source2 with Tab 22 integration
cmake --build build_msvc --config Release --target Source2
```

### CMakeLists.txt Conditional Logic

The calculator_advanced module is **optional** and only builds when all dependencies are found:

```cmake
if(Qt6_FOUND AND antlr4-runtime_FOUND AND SymEngine_FOUND)
    message(STATUS "calculator_advanced: All dependencies found - building Tab 22")
    add_subdirectory(calculator_advanced)
    
    # Link to Source2
    if(TARGET calculator_advanced)
        target_link_libraries(Source2 PRIVATE calculator_advanced)
    endif()
else()
    message(WARNING "calculator_advanced: Missing dependencies - Tab 22 disabled")
endif()
```

---

## Test Suite

### Run All Tests

```powershell
cd build_msvc
ctest --config Release -R calculator --verbose
```

### Individual Test Files

1. **test_parser.cpp** (10 tests) - ANTLR4 grammar validation
   - Basic expressions, derivatives, integrals, series, parametric, ODE, UQFF patterns
2. **test_symbolic.cpp** (11 tests) - SymEngine operations
   - Differentiation, integration, simplification, expansion, solving, Taylor series
3. **test_polynomial.cpp** (10 tests) - GSL polynomial solver
   - Degrees 1-26, real root filtering, complex roots, multiple roots
4. **test_dimensional.cpp** (11 tests) - SI unit system
   - Base quantities, force, energy, UQFF dimensional analysis
5. **test_uqff_integration.cpp** (15 tests) - End-to-end UQFF validation
   - All 11 equations with realistic physics parameters

**Total**: 57 test cases, ~1,100 lines

---

## Implemented UQFF Equations

### 11 Complete Implementations

1. **F_U_Bi_i** - Universal buoyancy force
   - Parameters: M_p=1.673×10⁻²⁷ kg, r_p=8.4×10⁻¹⁶ m, F_rel=4.30×10³³ N, E_LEP=200 GeV
2. **Um** - Universal magnetism
   - Time evolution: t=115 days, ω_c=1.585×10⁻⁸ rad/s, γ=5×10⁻⁵ day⁻¹
3. **g_MUGE_H** - Hydrogen atom gravity
   - Bohr radius: r=5.29×10⁻¹¹ m
4. **g_Magnetar** - Magnetar at 1 year
   - B_crit=4.4×10¹³ T, M=3.0×10³⁰ kg
5. **g_SgrA** - Sagittarius A* at 1 Gyr
   - M_0=4.3×10⁶ M_☉, t_age=9 Gyr
6. **P_alpha** - Alpha clustering probability
   - E_th=8×10⁻¹³ J (5 MeV threshold)
7. **R_EU** - Electric Universe ratio
   - F_EM/F_g validation (Solar mass test case)
8. **tau_gyro** - Gyroscopic torque
   - I=10⁴⁰ kg·m² (galactic moment of inertia)
9. **g_compressed** - 26-layer compressed gravity
   - Q_i=0.99^i (quantum state decay)
10. **eta_LENR** - Neutron production rate
    - k_η=10⁻¹¹³, [SSq]=0.57 (Widom-Larsen calibration)
11. **Equation Catalog** - Search, retrieval, LaTeX export
    - listEquations(), searchEquations(), getEquation()

---

## Features

### Parser (ANTLR4)
- **Grammar**: Math.g4 (100 lines)
- **Operators**: +, -, *, /, ^, sqrt(), sin(), cos(), log(), exp()
- **Calculus**: ∂/∂ (derivatives), ∫ (integrals), ∑ (summations), ∏ (products)
- **Advanced**: Parametric equations, ODEs, series expansion
- **UQFF Detection**: F_U_Bi_i, Um, g_MUGE, g_Magnetar, g_SgrA patterns

### Symbolic Math (SymEngine)
- **Operations**: Differentiation, integration, simplification, expansion
- **Solving**: Linear, quadratic, polynomial equations
- **Series**: Taylor expansion at arbitrary points
- **LaTeX**: Export equations to LaTeX format

### Polynomial Solver (GSL)
- **Degree Range**: 1 to 26 (maximum)
- **Root Types**: Real and complex
- **Filtering**: Extract real-only roots for physical applications
- **Special Cases**: Zero constant, multiple roots, unity roots

### Dimensional Analysis
- **SI Base Units**: M, L, T, I, Θ, N, J (mass, length, time, current, temperature, amount, luminous intensity)
- **Derived Units**: Force (M·L·T⁻²), Energy (M·L²·T⁻²), Power (M·L²·T⁻³)
- **UQFF Quantities**: F_U_Bi_i, Um, g validation with dimensional consistency

### Plotting (QCustomPlot)
- **2D Graphs**: Line plots, scatter plots
- **Interactions**: Zoom, pan, export to PNG/PDF
- **Multi-plot**: Multiple functions on same axes
- **Styling**: Colors, line widths, axis labels, legends

---

## User Interface (Tab 22)

### CalculatorWidget Features
- **Input Field**: Single-line expression editor with syntax highlighting
- **Result Display**: Multi-line output with scrolling
- **Equation Browser**: List view of available UQFF equations
- **LaTeX Preview**: Rendered equation display
- **Plot Canvas**: Interactive graph area
- **Control Buttons**: Compute, Clear, Export, Help

### Workflow
1. User enters expression (e.g., "∂/∂x (x^2 + 3*x)")
2. Parser validates syntax and detects operation type
3. Symbolic engine computes result
4. Result displayed with LaTeX formatting
5. Optional: Plot function or export to file

---

## Integration Status

### source2.cpp Changes
- **Line 80**: Added `#include "calculator_advanced/include/calculator_widget.h"`
- **Line 210**: Increased `MAX_WINDOWS` from 21 to 22
- **Lines 15440-15455**: Created Tab 22 case (index 21) with CalculatorWidget instantiation
- **Commit**: ac9cdfe (March 3, 2026)

### CMakeLists.txt Changes
- **Line ~136**: Added `add_subdirectory(calculator_advanced)` with dependency checks
- **Line ~500**: Linked calculator_advanced to Source2 target
- **Commit**: ac9cdfe (March 3, 2026)

### Git Status
- **Branch**: master
- **Remote**: https://github.com/Daniel8Murphy0007/Star-Magic.git
- **Files Added**: 17 (14 new sources/tests, 3 modified)
- **Insertions**: 3,190 lines
- **Deletions**: 64 lines

---

## Troubleshooting

### Common Issues

**Issue**: `qt6-base not found` or `qt:x64-windows` fails  
**Solution**: Update vcpkg (`git pull` in C:\vcpkg) and ensure 7zip 26.0.0+ installed

**Issue**: `ANTLR4 parser not found`  
**Solution**: Already installed (4.13.2), verify with `C:\vcpkg\vcpkg.exe list antlr4`

**Issue**: `SymEngine not found`  
**Solution**: Install via `C:\vcpkg\vcpkg.exe install symengine:x64-windows`

**Issue**: Tab 22 doesn't appear in Source2  
**Solution**: Check CMake output for "calculator_advanced: Missing dependencies" warning

**Issue**: Tests fail to compile  
**Solution**: Ensure all 5 test files in `calculator_advanced/tests/` directory

### Verification Commands

```powershell
# Check installed packages
C:\vcpkg\vcpkg.exe list | Select-String "qt|antlr4|symengine|gsl|eigen3|qcustomplot"

# Check CMake configuration
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake 2>&1 | Select-String "calculator"

# Check if target exists
cmake --build build_msvc --target calculator_advanced --config Release --dry-run
```

---

## Next Steps

### For Build Completion
1. Install missing vcpkg dependencies (Qt6, SymEngine, GSL, QCustomPlot)
2. Configure CMake with vcpkg toolchain
3. Build calculator_advanced library
4. Build Source2 with Tab 22 integration
5. Run test suite to validate all 57 test cases
6. Launch Source2.exe and verify Tab 22 appears

### For Production Deployment
1. Verify all tests pass (57/57)
2. Test UQFF equations with real astronomical data
3. Benchmark performance (symbolic operations, plotting)
4. Add user documentation for Tab 22 features
5. Create release build with optimizations
6. Package dependencies for distribution

---

**Author**: Daniel T. Murphy  
**Date**: March 3, 2026  
**Repository**: https://github.com/Daniel8Murphy0007/Star-Magic  
**Thread Reference**: https://x.com/i/grok/share/533da64c6ded4ada90fc83b522d90fe6  
**Analysis**: [GROK_THREAD_B6D9BC22_ANALYSIS.md](../GROK_THREAD_B6D9BC22_ANALYSIS.md)
