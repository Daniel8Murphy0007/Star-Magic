# Star-Magic UQFF: Complete JSON Installation Guide

## 📋 Overview
This guide provides **step-by-step instructions** for installing all JSON-related dependencies required for the Star-Magic UQFF workspace.

---

## 🎯 Phase 4 Workspace State
- **Current Phase**: Phase 4 - Integration & Testing
- **Restore Point**: Git commit `5d0a633`
- **Library**: `build/libUQFFCore.a` (66.52 MB, 118 modules compiled)
- **Primary Platform**: `source2.cpp` (HEAD PROGRAM)

---

## 📦 JSON Dependencies Required

### 1️⃣ **JSON for Modern C++ (nlohmann/json)**
**Purpose**: Primary JSON library for C++ (header-only, no compilation needed)

#### Installation via vcpkg:
```powershell
# Step 1: Verify vcpkg is installed
Test-Path C:\vcpkg

# Step 2: Update vcpkg (if needed)
cd C:\vcpkg
git pull

# Step 3: Install nlohmann-json
.\vcpkg install nlohmann-json:x64-windows

# Step 4: Integrate with Visual Studio/CMake
.\vcpkg integrate install
```

#### Manual Installation (Alternative):
```powershell
# Step 1: Download single header file
Invoke-WebRequest -Uri "https://github.com/nlohmann/json/releases/download/v3.11.3/json.hpp" -OutFile "C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\include\json.hpp"

# Step 2: Create include directory if not exists
New-Item -ItemType Directory -Force -Path "C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\include"

# Step 3: Verify download
Get-Item "C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\include\json.hpp"
```

#### Usage in C++:
```cpp
#include <nlohmann/json.hpp>  // If using vcpkg
// OR
#include "include/json.hpp"   // If using manual installation

using json = nlohmann::json;

int main() {
    json data = {
        {"module", "UQFFCore"},
        {"version", "1.0"},
        {"modules_compiled", 118}
    };
    std::cout << data.dump(4) << std::endl;
    return 0;
}
```

---

### 2️⃣ **RapidJSON** (Alternative high-performance JSON library)
**Purpose**: Ultra-fast JSON parsing/generation (header-only)

#### Installation via vcpkg:
```powershell
# Install RapidJSON
cd C:\vcpkg
.\vcpkg install rapidjson:x64-windows
.\vcpkg integrate install
```

#### Manual Installation:
```powershell
# Step 1: Clone repository
cd C:\Users\tmsjd\source\repos
git clone https://github.com/Tencent/rapidjson.git

# Step 2: Copy headers to workspace
Copy-Item -Path "C:\Users\tmsjd\source\repos\rapidjson\include\rapidjson" -Destination "C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\include\" -Recurse

# Step 3: Verify
Test-Path "C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\include\rapidjson"
```

#### Usage in C++:
```cpp
#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>

using namespace rapidjson;

int main() {
    Document doc;
    doc.SetObject();
    doc.AddMember("module", "UQFFCore", doc.GetAllocator());
    doc.AddMember("modules", 118, doc.GetAllocator());
    
    StringBuffer buffer;
    Writer<StringBuffer> writer(buffer);
    doc.Accept(writer);
    std::cout << buffer.GetString() << std::endl;
    return 0;
}
```

---

### 3️⃣ **JSON Schema Validator**
**Purpose**: Validate JSON against schemas (for configuration files)

#### Installation via vcpkg:
```powershell
cd C:\vcpkg
.\vcpkg install nlohmann-json-schema-validator:x64-windows
.\vcpkg integrate install
```

---

### 4️⃣ **CMake Configuration for JSON Libraries**

#### Update `CMakeLists.txt`:
```cmake
cmake_minimum_required(VERSION 3.28)
project(StarMagicUQFF CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ========== VCPKG INTEGRATION ==========
if(EXISTS "C:/vcpkg/scripts/buildsystems/vcpkg.cmake")
    set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake"
        CACHE STRING "Vcpkg toolchain file")
endif()

# ========== FIND JSON LIBRARIES ==========
find_package(nlohmann_json CONFIG REQUIRED)
# find_package(RapidJSON CONFIG REQUIRED)  # Uncomment if using RapidJSON

# ========== LINK LIBRARIES ==========
add_executable(MAIN_1 MAIN_1.cpp)
target_link_libraries(MAIN_1 PRIVATE nlohmann_json::nlohmann_json)

add_executable(Source2 "source2(HEAD PROGRAM).cpp")
target_link_libraries(Source2 PRIVATE nlohmann_json::nlohmann_json)
```

---

## 🔧 Verification Steps

### Verify vcpkg Installation:
```powershell
# Check vcpkg version
cd C:\vcpkg
.\vcpkg version

# List installed packages
.\vcpkg list
```

### Verify JSON Libraries:
```powershell
# Check nlohmann-json
.\vcpkg list | Select-String "nlohmann"

# Check RapidJSON (if installed)
.\vcpkg list | Select-String "rapidjson"
```

### Test Compilation:
```powershell
# Create test file
@"
#include <iostream>
#include <nlohmann/json.hpp>

int main() {
    nlohmann::json test = {{"status", "success"}, {"modules", 118}};
    std::cout << test.dump(2) << std::endl;
    return 0;
}
"@ | Out-File -FilePath test_json.cpp -Encoding utf8

# Compile
g++ -std=c++17 test_json.cpp -o test_json -I"C:\vcpkg\installed\x64-windows\include"

# Run
.\test_json.exe

# Expected output:
# {
#   "modules": 118,
#   "status": "success"
# }

# Cleanup
Remove-Item test_json.cpp, test_json.exe
```

---

## 🛠️ Troubleshooting

### Issue 1: vcpkg not found
**Solution**:
```powershell
# Install vcpkg
git clone https://github.com/Microsoft/vcpkg.git C:\vcpkg
cd C:\vcpkg
.\bootstrap-vcpkg.bat
```

### Issue 2: PowerShell version too old
**Solution**:
```powershell
# Check version
$PSVersionTable.PSVersion

# If < 7.5.3, install PowerShell 7:
winget install --id Microsoft.Powershell --source winget
```

### Issue 3: CMake cannot find nlohmann_json
**Solution**:
```powershell
# Re-run vcpkg integrate
cd C:\vcpkg
.\vcpkg integrate install

# Clean and reconfigure CMake
Remove-Item -Recurse -Force C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\build
cmake -B build -G "MinGW Makefiles" -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake
```

### Issue 4: Include path errors
**Solution**:
Add to `.vscode/c_cpp_properties.json`:
```json
{
    "configurations": [
        {
            "name": "Win32",
            "includePath": [
                "${workspaceFolder}/**",
                "C:/vcpkg/installed/x64-windows/include"
            ],
            "defines": ["_DEBUG", "UNICODE"],
            "compilerPath": "C:/msys64/mingw64/bin/g++.exe",
            "cStandard": "c17",
            "cppStandard": "c++17",
            "intelliSenseMode": "windows-gcc-x64"
        }
    ],
    "version": 4
}
```

---

## 📚 Additional JSON Tools (Optional)

### JSON Viewer Extension for VS Code:
```powershell
# Install from command line
code --install-extension ZainChen.json
```

### JSON Formatter (command-line):
```powershell
# Install jq for PowerShell JSON formatting
winget install jqlang.jq

# Example usage:
Get-Content config.json | jq .
```

---

## ✅ Installation Checklist

- [ ] PowerShell 7.5.3+ installed
- [ ] vcpkg installed at `C:\vcpkg`
- [ ] vcpkg integrated with `.\vcpkg integrate install`
- [ ] nlohmann-json installed via vcpkg
- [ ] CMakeLists.txt updated with vcpkg toolchain
- [ ] Test compilation successful
- [ ] `.vscode/c_cpp_properties.json` includes vcpkg paths
- [ ] Workspace settings synchronized to restore point `5d0a633`

---

## 🎯 Next Steps After Installation

1. **Reload VS Code** to apply all configuration changes
2. **Verify auto-run tasks** execute on folder open (workspace initialization banner)
3. **Test CMake configuration**:
   ```powershell
   cmake -B build -G "MinGW Makefiles"
   cmake --build build --target UQFFCore
   ```
4. **Proceed to Phase 4 Integration Testing** (see BUILD_STATUS.md)

---

## 📞 Support References

- **nlohmann/json documentation**: https://json.nlohmann.me/
- **RapidJSON documentation**: https://rapidjson.org/
- **vcpkg documentation**: https://vcpkg.io/
- **Session recovery file**: `copilot_thread_capture_source.cpp`
- **Build status**: `BUILD_STATUS.md`

---

**Generated**: 2025-01-14  
**Workspace**: Star-Magic UQFF Phase 4  
**Restore Point**: 5d0a633 (Phase 3 Complete - 118 modules, 66.52 MB)
