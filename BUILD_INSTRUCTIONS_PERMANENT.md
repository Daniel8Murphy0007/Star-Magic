# PERMANENT BUILD INSTRUCTIONS - READ THIS EVERY TIME

**Created: November 16, 2025**
**Last Updated: December 4, 2025 - Phase 32: Qt Networking + Virgo Cluster Integration**
**DO NOT DELETE THIS FILE**

## ✅ SAVE POINT: December 4, 2025 18:58 PM

### FULLY WORKING STATE - WOLFRAM + QT6 + GROK + VIRGO CLUSTER ACTIVE

- **Build System**: Visual Studio 2022 (MSVC 19.44.35207) ✅
- **Wolfram WSTP**: ENABLED (source168-173 integrated) ✅
- **Qt6**: 6.10.0 ENABLED (QCoreApplication + QNetworkAccessManager) ✅
- **Grok AI**: XAI_API_KEY configured ✅
- **Tracing System**: uqff_tracing.h integrated ✅
- **Executable**: 1.35 MB (built Dec 4 @ 18:58) ✅ RUNTIME VERIFIED 17:40:40
- **Menu Options**: 1-12 (includes Wolfram 9-11, Exit 12) ✅
- **Physics Terms**: 6,643+ registered (Virgo Cluster additions) ✅
- **Preprocessor**: USE_EMBEDDED_WOLFRAM=ON, HAVE_QT6 defined ✅
- **Latest Commit**: 5a6346f - Phase 32 workspace update ✅

### Critical Files Modified (Dec 4, 2025)

1. **MAIN_1_CoAnQi.cpp** - Added Qt networking includes (lines 206-210)
   - QCoreApplication for Qt networking (QNetworkAccessManager)
   - Comment formatting cleanup and alignment
   - Status: Compiles and runs successfully
2. **source82_wolfram.cpp** - Virgo Cluster cosmological physics
   - VirgoClusterMassTerm (NFW-like mass profile, ~1.2e15 M_sun)
   - VirgoClusterIntraclusterMediumTerm (hot ICM gas, T~2-4 keV)
   - Total additions: 447+ lines of new physics code
3. **CMakeLists.txt** - Updated to Phase 32 (commit 5a6346f)
4. **Star-Magic.code-workspace** - Updated session restore point

### What YOU Have NOW

1. **Build System**: Visual Studio 2022 (MSVC-ONLY, Wolfram WSTP requires it) ✅
2. **Dependencies**: From vcpkg x64-windows (8 total: Qt6, VTK, CURL, SQLite3, OpenCV, libwebsockets, AWS SDK, Python/pybind11) ✅
3. **Wolfram Integration**: 6 source files (168-173) with 19+11+8 astronomical systems ✅
4. **AI Integration**: Grok API + Wolfram hypergraph working together ✅

---

## 🛑 STOP DOING THIS

### **NEVER** run these commands again

```powershell
# ❌ WRONG - Uses MinGW
cmake . -B .\build -G "MinGW Makefiles"

# ❌ WRONG - Resets everything
git reset --hard

# ❌ WRONG - Commits without understanding
git commit -m "..."
```

---

## ✅ CORRECT BUILD PROCESS (December 1, 2025 - VERIFIED WORKING)

### Step 0: SAVE POINT RESTORE (if needed)

```powershell
# This exact configuration is working as of Dec 1, 2025 3:05 PM
# USE_EMBEDDED_WOLFRAM=ON, Qt6 6.10.0, XAI_API_KEY set
# Executable: build_msvc\Release\MAIN_1_CoAnQi.exe (1.31 MB)

# To restore this exact state:
git log --oneline | Select-Object -First 5  # Find Dec 1, 2025 commit
# git checkout <commit-hash>  # If needed
```

### Step 1: One-Time vcpkg Setup (ALREADY DONE)

```powershell
# ✅ COMPLETED - Dependencies installed for Visual Studio (x64-windows)
cd C:\vcpkg
.\vcpkg install opengl:x64-windows
.\vcpkg install glew:x64-windows
.\vcpkg install glfw3:x64-windows
.\vcpkg install glm:x64-windows
.\vcpkg install assimp:x64-windows
.\vcpkg install yaml-cpp:x64-windows
.\vcpkg install nlohmann-json:x64-windows
.\vcpkg integrate install
```

### Step 2: VERIFIED WORKING BUILD COMMANDS (Dec 1, 2025)

```powershell
# Clean rebuild with Wolfram WSTP enabled
Remove-Item -Recurse -Force build_msvc -ErrorAction SilentlyContinue

# Configure with Wolfram (REQUIRED for menu options 9-11)
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64 -DUSE_EMBEDDED_WOLFRAM=ON

# Build (8 parallel jobs, UPX compression automatic)
cmake --build build_msvc --config Release --target MAIN_1_CoAnQi -j 8

# Run (Wolfram + Grok both active)
.\build_msvc\Release\MAIN_1_CoAnQi.exe
```

### Step 3: Environment Variables (ALREADY SET)

```powershell
# ✅ XAI_API_KEY - Permanent user environment variable
$env:XAI_API_KEY  # Should show your configured API key

# If not set, configure with your xAI API key from https://console.x.ai
# [System.Environment]::SetEnvironmentVariable('XAI_API_KEY', 'your-xai-api-key-here', 'User')
```

---

## 🛑 CRITICAL: DO NOT CHANGE THESE (Dec 1, 2025 Working Config)

### CMakeLists.txt - VERIFIED CORRECT SETTINGS

```cmake
# ✅ Line 32: CORRECT (vcpkg x64-windows)
include_directories(C:/vcpkg/installed/x64-windows/include)

# ✅ Line 38: CORRECT (Qt6 path)
set(CMAKE_PREFIX_PATH "C:/Qt/6.10.0/msvc2022_64" ${CMAKE_PREFIX_PATH})

# ✅ Line 159: CORRECT (Wolfram enabled by default)
option(USE_EMBEDDED_WOLFRAM "Enable embedded Wolfram kernel via WSTP" ON)

# ✅ Line 163: CORRECT (WSTP directory)
set(WSTP_DIR "C:/Program Files/Wolfram Research/Wolfram Engine/14.3/SystemFiles/Links/WSTP/DeveloperKit/Windows-x86-64/CompilerAdditions")

# ✅ Line 191: CORRECT (source178 NOT compiled separately, only included)
add_executable(MAIN_1_CoAnQi MAIN_1_CoAnQi.cpp)
```

### uqff_tracing.h - VERIFIED CORRECT (Dec 1, 2025)

```cpp
// ✅ Lines 40-46: CORRECT (TRACE_ prefix avoids Windows WSTP conflicts)
enum class TraceLevel {
    TRACE_DEBUG = 0,
    TRACE_INFO = 1,
    TRACE_WARN = 2,
    TRACE_ERROR = 3,
    TRACE_FATAL = 4
};
```

### MAIN_1_CoAnQi.cpp - VERIFIED CORRECT

```cpp
// ✅ Line 206: CORRECT (source178 included for Grok)
#include "source178_grok_api.cpp"

// ✅ Line 214: CORRECT (tracing header)
#include "uqff_tracing.h"

// ✅ Lines 21804, 21812, 21907: CORRECT (TraceLevel::TRACE_INFO)
TRACE_EVENT("...", TraceLevel::TRACE_INFO);
```

---

## 🛑 NEVER DO THIS AGAIN

```powershell
# ❌ WRONG - Uses MinGW (Wolfram WSTP requires MSVC)
cmake . -B .\build -G "MinGW Makefiles"

# ❌ WRONG - Disables Wolfram (menu options 9-11 will disappear)
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64 -DUSE_EMBEDDED_WOLFRAM=OFF

# ❌ WRONG - Changes vcpkg path back to mingw
# include_directories(C:/vcpkg/installed/x64-mingw-dynamic/include)

# ❌ WRONG - Recompiles source178 separately (causes duplicate symbols)
# add_executable(MAIN_1_CoAnQi MAIN_1_CoAnQi.cpp source178_grok_api.cpp)

# ❌ WRONG - Resets working state
git reset --hard
```

**To this:**

```cmake
set(CMAKE_PREFIX_PATH "C:/vcpkg/installed/x64-windows")
```

### Step 3: Configure with Visual Studio (NOT MinGW)

```powershell
# Delete old build
Remove-Item -Recurse -Force .\build -ErrorAction SilentlyContinue

# Configure for Visual Studio 2022
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 `
  -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake `
  -DVCPKG_TARGET_TRIPLET=x64-windows

# Build
cmake --build build --config Debug -j 12
```

### Step 4: Verify UQFFCore.hpp Integration

Your CMakeLists.txt **ALREADY HAS** this (lines 233-238):

```cmake
target_sources(UQFFCore PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES Core/UQFFCore.hpp
)
```

This is **CORRECT** - don't change it.

---

## 📋 Git Workflow (Current Situation)

### Your Current Staged Changes

```
Changes to be committed:
  deleted:    glm/ (entire directory - 1539 files)
  deleted:    module_backups_20251104_105304/ (backup files)
  deleted:    nlohmann/json.hpp
  modified:   source7.cpp
  modified:   test_core_integration.vcxproj
```

### Safe Commit Command

```powershell
git commit -m "Remove local GLM and json copies (using vcpkg versions); update source7 and test project"
```

**This is safe because:**

- GLM is installed via vcpkg
- nlohmann/json is installed via vcpkg  
- module_backups are duplicates (originals in Core/Modules/)

---

## 🔄 Breaking the Cycle

### When You Come Back Next Time

1. **Read this file first**: `BUILD_INSTRUCTIONS_PERMANENT.md`
2. **Check git status**: See what's actually changed
3. **DON'T** let AI agent modify CMakeLists.txt unless you explicitly review changes
4. **Use Visual Studio generator**, not MinGW
5. **Verify vcpkg triplet**: x64-windows (not x64-mingw-dynamic)

### Files That Should NEVER Change

- Core/UQFFCore.hpp (8,314 bytes, 155 includes)
- Core/SystemCatalogue.cpp/.hpp
- Core/UQFFModule4.cpp/.hpp
- Core/Modules/*.cpp (128 files, Nov 16 4:21 AM)

### Files That Can Change

- CMakeLists.txt (only to fix vcpkg paths or add new targets)
- Individual source*.cpp files (your code)
- Build artifacts (build/, Debug/, x64/)

---

## 🎯 Your Actual Goal

**Build UQFFCore.lib** (not libUQFFCore.a) with:

- Visual Studio 2022 toolchain
- Core/UQFFCore.hpp as public header
- All 119 Core/Modules compiled
- Dependencies from vcpkg x64-windows

**Then build these executables:**

- MAIN_1.exe
- MAIN_1_CoAnQi.exe
- source4-162.exe (155 programs)
- Source134, 163-167.exe

---

## 📞 Quick Reference Commands

### Check what's changed

```powershell
git status
git diff CMakeLists.txt | Select-Object -First 50
```

### Restore CMakeLists.txt if broken

```powershell
git checkout -- CMakeLists.txt
```

### Clean build

```powershell
Remove-Item -Recurse -Force build, Debug, x64 -ErrorAction SilentlyContinue
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows
cmake --build build --config Debug
```

### List what's built

```powershell
Get-ChildItem -Recurse -Filter "*.exe" | Select-Object Directory, Name, Length
Get-ChildItem -Recurse -Filter "*.lib" | Select-Object Directory, Name, Length
```

---

## ⚡ Emergency Reset (If Everything Breaks)

**ONLY USE IF COMPLETELY BROKEN:**

```powershell
# Save current CMakeLists.txt
Copy-Item CMakeLists.txt CMakeLists_BROKEN_$(Get-Date -Format 'yyyyMMdd_HHmmss').txt

# Restore from RESTORE_POINT
Copy-Item RESTORE_POINT_16NOV2025_651AM\CMakeLists_BACKUP.txt CMakeLists.txt

# Update paths to x64-windows
(Get-Content CMakeLists.txt) -replace 'x64-mingw-dynamic', 'x64-windows' | Set-Content CMakeLists.txt

# Clean build
Remove-Item -Recurse -Force build, Debug, x64 -ErrorAction SilentlyContinue
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows
cmake --build build --config Debug
```

---

## 🔍 Verification Checklist

After building, verify:

- [ ] UQFFCore.lib exists in build/Debug/ (not libUQFFCore.a)
- [ ] MAIN_1.exe runs
- [ ] MAIN_1_CoAnQi.exe runs
- [ ] Can include Core/UQFFCore.hpp in test program
- [ ] All 119 modules accessible via UQFFCore namespace

---

**REMEMBER**: The problem isn't the build - it's using **MinGW instead of Visual Studio**.

**FIX**: Change vcpkg paths from `x64-mingw-dynamic` → `x64-windows` and use VS generator.

**END OF PERMANENT INSTRUCTIONS**
