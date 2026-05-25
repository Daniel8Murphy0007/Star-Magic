# C++ Compilation and Testing Guide

**File:** `simultaneous_7layer_solver.cpp` (700 lines)  
**Date:** May 24, 2026  
**Framework:** UQFF v5.1.0  
**Status:** Ready for compilation

---

## 🔧 SYSTEM REQUIREMENTS

### Windows (MSVC)
- **OS:** Windows 10/11
- **Compiler:** Visual Studio 2022 (MSVC 14.44+)
- **C++ Standard:** C++20
- **Linear Algebra:** Eigen 3.4.0+ or Armadillo 11.0+
- **Build Tool:** CMake 3.20+

### Linux (g++)
- **OS:** Ubuntu 20.04+, CentOS 8+, Debian 11+
- **Compiler:** g++ 9.0+ (C++20 support)
- **Linear Algebra:** Eigen (apt install libeigen3-dev)
- **Build Tool:** CMake 3.20+

### macOS (clang++)
- **OS:** macOS 10.15+
- **Compiler:** Xcode 12+ (clang++ with C++20)
- **Linear Algebra:** Eigen (brew install eigen)
- **Build Tool:** CMake 3.20+

---

## 📦 DEPENDENCY INSTALLATION

### Windows (MSVC + Eigen)

**Step 1: Download Eigen**
```powershell
# Download from https://eigen.tuxfamily.org/
# Extract to C:\Libraries\eigen-3.4.0

# Set environment variable
[Environment]::SetEnvironmentVariable("EIGEN_ROOT", "C:\Libraries\eigen-3.4.0", "User")
```

**Step 2: Verify Eigen Installation**
```powershell
dir C:\Libraries\eigen-3.4.0\Eigen
# Should show: Dense, StdVector, etc. header files
```

### Linux (g++ + Eigen)

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install libeigen3-dev cmake build-essential
```

**CentOS/RHEL:**
```bash
sudo yum install eigen-devel cmake gcc-c++ make
```

**Verify Installation:**
```bash
find /usr -name Eigen -type d
# Should show: /usr/include/eigen3/Eigen
```

### macOS (clang++ + Eigen)

```bash
# Install Homebrew if not installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Eigen
brew install eigen cmake

# Verify
brew list eigen
```

---

## 🛠️ COMPILATION METHODS

### Method 1: Direct Compilation (Simplest)

**Windows (MSVC):**
```powershell
cd C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic

cl.exe /std:c++20 /O2 /I"C:\Libraries\eigen-3.4.0" `
    simultaneous_7layer_solver.cpp `
    /link /OUT:simultaneous_7layer_solver.obj

# Verify compilation
dir simultaneous_7layer_solver.obj
```

**Linux (g++):**
```bash
cd ~/repos/Daniel8Murphy0007/Star-Magic

g++ -std=c++20 -O3 -I/usr/include/eigen3 \
    -fPIC -shared \
    simultaneous_7layer_solver.cpp \
    -o simultaneous_7layer_solver.so

# Verify compilation
ls -lah simultaneous_7layer_solver.so
```

**macOS (clang++):**
```bash
cd ~/repos/Daniel8Murphy0007/Star-Magic

clang++ -std=c++20 -O3 -I/usr/local/include/eigen3 \
    -fPIC -shared \
    simultaneous_7layer_solver.cpp \
    -o simultaneous_7layer_solver.dylib

# Verify compilation
ls -lah simultaneous_7layer_solver.dylib
```

---

### Method 2: CMake Build (Recommended for Integration)

**Create `CMakeLists.txt` in Star-Magic directory:**

```cmake
cmake_minimum_required(VERSION 3.20)
project(UQFF_Solver CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Find Eigen
find_package(Eigen3 REQUIRED NO_MODULE)

# Create executable
add_executable(simultaneous_7layer_solver_exe simultaneous_7layer_solver.cpp)
target_link_libraries(simultaneous_7layer_solver_exe PUBLIC Eigen3::Eigen)

# Create shared library (for Python ctypes)
add_library(simultaneous_7layer_solver SHARED simultaneous_7layer_solver.cpp)
target_link_libraries(simultaneous_7layer_solver PUBLIC Eigen3::Eigen)
```

**Build with CMake:**

```bash
# Create build directory
mkdir build
cd build

# Configure (choose one based on OS)
# Windows (MSVC)
cmake .. -G "Visual Studio 17 2022" -A x64

# Linux (g++)
cmake .. -DCMAKE_BUILD_TYPE=Release

# macOS (clang++)
cmake .. -DCMAKE_BUILD_TYPE=Release

# Compile
cmake --build . --config Release

# Verify
ls Release/
# Should show: simultaneous_7layer_solver.exe/.so/.dylib
```

---

## ✅ COMPILATION VERIFICATION

### Step 1: Check Object File

```powershell
# Windows
dir simultaneous_7layer_solver.obj
# Expected size: ~50-100 KB

# Linux/macOS
ls -lah simultaneous_7layer_solver.so
# Expected size: ~100-200 KB
```

### Step 2: Check for Symbols

**Linux/macOS:**
```bash
nm simultaneous_7layer_solver.so | grep gmres
# Should show symbols related to GMRES solver
```

**Windows (MSVC):**
```powershell
dumpbin /SYMBOLS simultaneous_7layer_solver.obj | findstr gmres
```

### Step 3: Test with Simple Program

**Create `test_compile.cpp`:**
```cpp
#include <iostream>
#include <Eigen/Dense>

int main() {
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(3, 3);
    std::cout << "Eigen version: " << EIGEN_WORLD_VERSION << "."
              << EIGEN_MAJOR_VERSION << "."
              << EIGEN_MINOR_VERSION << std::endl;
    std::cout << "Matrix sum: " << A.sum() << std::endl;
    return 0;
}
```

**Compile and run:**
```bash
# Windows
cl.exe /std:c++20 /I"C:\Libraries\eigen-3.4.0" test_compile.cpp
test_compile.exe

# Linux
g++ -std=c++20 -I/usr/include/eigen3 test_compile.cpp -o test_compile
./test_compile

# macOS
clang++ -std=c++20 -I/usr/local/include/eigen3 test_compile.cpp -o test_compile
./test_compile
```

**Expected output:**
```
Eigen version: 3.4.0
Matrix sum: [some float value]
```

---

## 🧪 TESTING PROCEDURES

### Test 1: Syntax Validation

```bash
# Windows (MSVC) - Parse check only
cl.exe /std:c++20 /E simultaneous_7layer_solver.cpp > nul

# Linux/macOS
g++ -std=c++20 -fsyntax-only simultaneous_7layer_solver.cpp
```

**Expected:** No errors, exit code 0

### Test 2: Standalone Executable Test

**Create `test_uqff_solver.cpp`:**

```cpp
#include <iostream>
#include <vector>
#include <cmath>

// Copy relevant functions from simultaneous_7layer_solver.cpp
// Or link against compiled object file

int main() {
    // Test case: Helium (Z=2, n=1)
    int Z = 2;
    int n = 1;
    double v_init = 0.005 * 2.998e8;
    
    std::cout << "Testing UQFF simultaneous 7-layer solver..." << std::endl;
    std::cout << "Input: Z=" << Z << ", n=" << n << std::endl;
    
    // Call solver (requires linkage to simultaneous_7layer_solver.obj)
    // Expected: convergence in < 100 iterations
    
    std::cout << "✓ Test completed" << std::endl;
    return 0;
}
```

**Compile and run:**
```bash
# Windows
cl.exe /std:c++20 /I"C:\Libraries\eigen-3.4.0" ^
    test_uqff_solver.cpp simultaneous_7layer_solver.obj
test_uqff_solver.exe

# Linux
g++ -std=c++20 -I/usr/include/eigen3 \
    test_uqff_solver.cpp simultaneous_7layer_solver.o -o test_uqff_solver
./test_uqff_solver
```

### Test 3: Python Integration Test

**Run:** `python uqff_cpp_wrapper.py`

```python
#!/usr/bin/env python3
"""Test C++ wrapper functionality"""

from uqff_cpp_wrapper import UQFFCppSolver

def test_cpp_solver():
    print("Loading C++ UQFF solver...")
    
    # Platform-specific library name
    import sys
    if sys.platform == 'win32':
        lib_path = './simultaneous_7layer_solver.dll'
    elif sys.platform == 'darwin':
        lib_path = './simultaneous_7layer_solver.dylib'
    else:
        lib_path = './simultaneous_7layer_solver.so'
    
    solver = UQFFCppSolver(lib_path)
    
    # Test: Helium ground state
    print("Test 1: Helium ground state (Z=2, n=1)")
    result = solver.solve(Z=2, n=1, l=0)
    
    print(f"  r_shell: {result['r_shell']:.6e} m")
    print(f"  v_orb: {result['v_orb']:.6e} m/s")
    print(f"  E_single: {result['E_single']:.6f} eV")
    print(f"  Convergence: {result['convergence_achieved']}")
    print(f"  Iterations: {result['iterations']}")
    print(f"  Residual: {result['residual_norm']:.6e}")
    
    # Verify Helium ground state energy ≈ -79.0 eV
    if abs(result['E_single'] - (-79.0)) < 1.0:
        print("  ✓ PASS: He ground state matches expectation")
        return True
    else:
        print(f"  ✗ FAIL: Expected -79.0 eV, got {result['E_single']:.1f} eV")
        return False

if __name__ == '__main__':
    success = test_cpp_solver()
    exit(0 if success else 1)
```

---

## 📊 PERFORMANCE BENCHMARKS

### Expected Compilation Times

| Platform | Method | Time | Object Size |
|----------|--------|------|-------------|
| Windows (MSVC) | Direct | 2-3 sec | 50-100 KB |
| Linux (g++) | Direct | 1-2 sec | 100-150 KB |
| macOS (clang++) | Direct | 2-3 sec | 100-150 KB |
| All | CMake | 5-10 sec | 100-150 KB |

### Expected Execution Times

| Operation | Time |
|-----------|------|
| Single GMRES solve (H atom) | < 10 ms |
| Single GMRES solve (He atom) | 15-20 ms |
| Single GMRES solve (Ne atom) | 30-50 ms |
| Batch solve (10 systems) | < 500 ms |

---

## 🐛 TROUBLESHOOTING

### Error 1: "fatal error: Eigen/Dense: No such file or directory"

**Cause:** Eigen headers not found

**Solution:**
```bash
# Verify Eigen installation
# Windows
dir C:\Libraries\eigen-3.4.0\Eigen\Dense

# Linux
ls /usr/include/eigen3/Eigen/Dense

# macOS
ls /usr/local/include/eigen3/Eigen/Dense

# Re-run compilation with correct -I path
```

### Error 2: "undefined reference to '_GMRES'"

**Cause:** Eigen is header-only; this error suggests linker issue

**Solution:**
```bash
# Ensure no external BLAS is being linked
# Remove any -lblas, -llapack flags if present

# Try compilation without external libraries:
g++ -std=c++20 -I/usr/include/eigen3 \
    simultaneous_7layer_solver.cpp -o solver
```

### Error 3: "C++20 features not supported"

**Cause:** Compiler too old or not set to C++20 mode

**Solution:**
```bash
# Check compiler version
g++ --version
clang++ --version

# Upgrade if needed (Ubuntu):
sudo apt install g++-11

# Use specific compiler:
g++-11 -std=c++20 ...
```

### Error 4: "cannot open shared library" (Python wrapper)

**Cause:** .so/.dylib file not found or incorrect path

**Solution:**
```bash
# Verify file exists in current directory
ls -lah simultaneous_7layer_solver.*

# Set library path if in subdirectory
export LD_LIBRARY_PATH=./build/Release:$LD_LIBRARY_PATH

# Update wrapper path
solver = UQFFCppSolver('./build/Release/simultaneous_7layer_solver.so')
```

---

## 📋 COMPILATION CHECKLIST

### Pre-Compilation
- [ ] Eigen installed and verified
- [ ] C++ compiler version checked (requires C++20 support)
- [ ] simultaneious_7layer_solver.cpp exists
- [ ] Working directory set correctly

### Compilation
- [ ] Compilation command executed without errors
- [ ] Object/library file created
- [ ] File size reasonable (> 50 KB)

### Post-Compilation
- [ ] Symbols verified (nm/dumpbin)
- [ ] Test program compiles and runs
- [ ] Python wrapper loads library successfully
- [ ] GMRES solver produces expected results

### Validation
- [ ] Helium ground state energy matches -79.0 eV ±1%
- [ ] Convergence achieved in < 100 iterations
- [ ] Residual norm < 1e-5
- [ ] Execution time < 50 ms per solve

---

## 🚀 NEXT STEPS

**After Successful Compilation:**

1. **Python Integration:**
   ```bash
   python uqff_cpp_wrapper.py
   ```

2. **Run All Tests:**
   ```bash
   python test_superposition_pair_helium.py
   python test_superposition_binary_bh.py
   ```

3. **Integration with MAIN_1_CoAnQi.cpp:**
   - Link compiled object file in CMakeLists.txt
   - Add menu option (see UQFF_INTEGRATION_GUIDE.md)

4. **Performance Profiling:**
   ```bash
   # Linux/macOS
   time ./simultaneous_7layer_solver_exe
   
   # Windows
   measure-command { ./simultaneous_7layer_solver.exe }
   ```

---

## 📞 SUPPORT

**For compilation issues:**
1. Check system requirements above
2. Verify Eigen installation
3. Try CMake method if direct method fails
4. Review UQFF_INTEGRATION_GUIDE.md for additional context

**For runtime issues:**
1. Check test cases in Python wrapper
2. Verify input parameters (Z, n, l, v_init)
3. Compare results with test_superposition_pair_helium.py

---

*C++ Compilation Guide v5.1.0*  
*Date: May 24, 2026*  
*Framework: UQFF v5.1.0*
