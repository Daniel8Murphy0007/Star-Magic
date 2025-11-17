# PERMANENT BUILD INSTRUCTIONS - READ THIS EVERY TIME

**Created: November 16, 2025**
**DO NOT DELETE THIS FILE**

## ⚠️ CRITICAL: Your Configuration Status

### Current State (as of Nov 16, 2025 7:30 PM)

- **Git Status**: You have staged DELETIONS of temporary files (glm/, module_backups/, nlohmann/)
- **These deletions are SAFE** - you're using vcpkg versions
- **Only 2 real file changes**: source7.cpp, test_core_integration.vcxproj
- **CMakeLists.txt**: Configured for x64-mingw-dynamic (NOT Visual Studio yet)

### What YOU Want

1. **Build System**: Visual Studio 2022 (not MinGW)
2. **Master Header**: Core/UQFFCore.hpp properly integrated
3. **Dependencies**: From vcpkg x64-windows (NOT x64-mingw-dynamic)
4. **Stop the cycle**: Same steps repeated 14+ times

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

## ✅ CORRECT BUILD PROCESS

### Step 1: One-Time vcpkg Setup (if not done)

```powershell
# Install for Visual Studio (x64-windows), NOT MinGW
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

### Step 2: Update CMakeLists.txt for Visual Studio

**Change this line:**

```cmake
include_directories(C:/vcpkg/installed/x64-mingw-dynamic/include)
```

**To this:**

```cmake
include_directories(C:/vcpkg/installed/x64-windows/include)
```

**And change this:**

```cmake
set(CMAKE_PREFIX_PATH "C:/vcpkg/installed/x64-mingw-dynamic")
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
