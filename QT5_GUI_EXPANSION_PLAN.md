# Qt6 GUI Expansion Plan - Star-Magic UQFF
**Created:** December 1, 2025  
**Updated:** December 3, 2025 @ 13:45 PM  
**Status:** MAIN_1_CoAnQi OPERATIONAL, Source2 BROKEN (AWS SDK iostream errors)  
**Target:** Build source2.cpp (HEAD PROGRAM) with full Qt6 GUI ecosystem

---

## Current Status: Core Computational Engine Complete ✅

- **MAIN_1_CoAnQi.exe**: 2.00 MB (Dec 3 @ 13:21:45), 13 menu options, Wolfram + Grok AI functional
- **Runtime Verified**: Dec 3 @ 13:37:40 (6,643/6,785 terms registered)
- **Physics Terms**: 6,643 active across 116 integrated source files (97.9% registration)
- **Build System**: Visual Studio 2022 (MSVC 19.44.35207), CMake 3.31.0, C++20
- **AI Integrations**: Wolfram WSTP 14.3, Qt6 6.10.0, xAI Grok API
- **Source2 Status**: ❌ BROKEN (8 iostream LNK2001 errors from AWS SDK DLLs)

---

## Phase 1: Qt6 Installation via vcpkg ✅ COMPLETE

### Qt6 Components Installed

```powershell
# Qt6 6.10.0 ALREADY INSTALLED at C:/Qt/6.10.0/msvc2022_64/
# vcpkg integration: qt6-base, qt6-webenginewidgets

# Verification
C:\vcpkg\vcpkg.exe list | Select-String "qt6"
# qt6-base:x64-windows                               6.10.0
# qt6-webenginewidgets:x64-windows                   6.10.0 (used by Grok API)
```

---

## Phase 2: External Dependencies for source2.cpp

### Required Libraries (from source2.cpp analysis)

| Library | Purpose | Installation Command |
|---------|---------|----------------------|
| **VTK** | 3D visualization (scatter plots, charts) | `vcpkg install vtk:x64-windows` |
| **OpenCV** | Computer vision (video/image processing) | `vcpkg install opencv:x64-windows` |
| **libcurl** | HTTP/HTTPS requests (NASA/MAST APIs) | `vcpkg install curl:x64-windows` |
| **SQLite3** | Local database caching | `vcpkg install sqlite3:x64-windows` |
| **AWS SDK** | S3 storage, Cognito authentication | `vcpkg install aws-sdk-cpp[s3,cognito-idp]:x64-windows` |
| **PocketSphinx** | Speech recognition (voice commands) | Manual build or vcpkg (if available) |
| **pybind11** | Python embedding (GPT-4 AI) | `vcpkg install pybind11:x64-windows` |
| **Qalculate** | Symbolic math library | Manual build (not in vcpkg) |
| **WebSocket** | Real-time server communication | `vcpkg install websocketpp:x64-windows` |
| **nlohmann-json** | JSON parsing | `vcpkg install nlohmann-json:x64-windows` ✅ (already installed) |

### Installation Script

```powershell
# Phase 2A: vcpkg-available dependencies
C:\vcpkg\vcpkg.exe install vtk:x64-windows
C:\vcpkg\vcpkg.exe install opencv:x64-windows
C:\vcpkg\vcpkg.exe install curl:x64-windows
C:\vcpkg\vcpkg.exe install sqlite3:x64-windows
C:\vcpkg\vcpkg.exe install aws-sdk-cpp[s3,cognito-idp]:x64-windows
C:\vcpkg\vcpkg.exe install pybind11:x64-windows
C:\vcpkg\vcpkg.exe install websocketpp:x64-windows

# Phase 2B: Manual installs (PocketSphinx, Qalculate)
# TBD - Research best approach for Windows MSVC builds
```

---

## Phase 3: CMakeLists.txt Modifications

### Add source2.cpp Target

```cmake
# Add after MAIN_1_CoAnQi target (around line 230)

if(Qt5_FOUND)
    # source2.cpp - HEAD PROGRAM (Qt5 GUI with full ecosystem)
    add_executable(Source2_HeadProgram source2.cpp)
    
    # Qt5 AutoMOC for Q_OBJECT macros
    set_target_properties(Source2_HeadProgram PROPERTIES AUTOMOC ON AUTOUIC ON AUTORCC ON)
    
    # Link Qt5 modules
    target_link_libraries(Source2_HeadProgram PRIVATE
        Qt5::Core
        Qt5::Widgets
        Qt5::WebEngineWidgets
        Qt5::PrintSupport
        Qt5::Network
        Qt5::Multimedia
        Qt5::Svg
    )
    
    # Link VTK for visualization
    if(VTK_FOUND)
        target_link_libraries(Source2_HeadProgram PRIVATE ${VTK_LIBRARIES})
    endif()
    
    # Link OpenCV for vision processing
    if(OpenCV_FOUND)
        target_link_libraries(Source2_HeadProgram PRIVATE ${OpenCV_LIBS})
    endif()
    
    # Link libcurl for HTTP requests
    if(CURL_FOUND)
        target_link_libraries(Source2_HeadProgram PRIVATE CURL::libcurl)
    endif()
    
    # Link SQLite3 for local database
    if(SQLite3_FOUND)
        target_link_libraries(Source2_HeadProgram PRIVATE SQLite::SQLite3)
    endif()
    
    # Link AWS SDK for cloud services
    if(AWSSDK_FOUND)
        target_link_libraries(Source2_HeadProgram PRIVATE 
            aws-cpp-sdk-s3 
            aws-cpp-sdk-cognito-idp
        )
    endif()
    
    # Link pybind11 for Python embedding
    if(pybind11_FOUND)
        target_link_libraries(Source2_HeadProgram PRIVATE pybind11::embed)
    endif()
    
    # Link WebSocket
    # (Header-only library, no linking needed)
    
    # Link PocketSphinx (if available)
    # TBD - Manual path configuration
    
    # Link Qalculate (if available)
    # TBD - Manual path configuration
    
    message(STATUS "Source2_HeadProgram target added (Qt5 GUI)")
else()
    message(WARNING "Qt5 not found - source2.cpp (HEAD PROGRAM) will not be built")
endif()
```

### Add find_package Calls

```cmake
# Add after Qt5/Qt6 find_package calls (around line 50)

find_package(VTK QUIET)
find_package(OpenCV QUIET)
find_package(CURL QUIET)
find_package(SQLite3 QUIET)
find_package(AWSSDK COMPONENTS s3 cognito-idp QUIET)
find_package(pybind11 QUIET)
```

---

## Phase 4: source2.cpp Code Preparation

### API Keys to Update

Replace placeholder API keys in source2.cpp:

```cpp
// Line 96-104 (approximately)
#define NASA_API_KEY_1 "[PROMPT_FOR_NASA_API_KEY_1]" // Keep existing
#define NASA_API_KEY_2 "[PROMPT_FOR_NASA_API_KEY_2]" // Keep existing
#define MAST_API_KEY "[PROMPT_FOR_MAST_API_KEY]"   // Keep existing
#define OPENAI_API_KEY "your_openai_api_key_here"               // UPDATE THIS
#define COGNITO_CLIENT_ID "your_cognito_client_id"              // UPDATE THIS
```

### Optional Features (Can Be Disabled)

If external dependencies are unavailable, comment out:

- **PocketSphinx**: Voice recognition (lines with `pocketsphinx.h`)
- **Qalculate**: Symbolic math (lines with `qalculate.h`)
- **AWS SDK**: Cloud sync (lines with `aws/`)
- **pybind11**: Python GPT-4 integration (lines with `pybind11/`)

---

## Phase 5: Build and Test

### Build Commands

```powershell
# Clean rebuild with Qt5 GUI target
Remove-Item -Recurse -Force build_msvc -ErrorAction SilentlyContinue

# Configure with all dependencies
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64 `
    -DUSE_EMBEDDED_WOLFRAM=ON `
    -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake

# Build Source2_HeadProgram
cmake --build build_msvc --config Release --target Source2_HeadProgram -j 8

# Run
.\build_msvc\Release\Source2_HeadProgram.exe
```

### Expected Output

- **21 Browser Windows**: ALMA Cycle 12 parallel search interface
- **VTK Visualizations**: 3D scatter plots, charts, graphs
- **NASA/MAST Integration**: Real-time astronomical data fetching
- **Voice Commands**: PocketSphinx speech recognition
- **Cloud Sync**: AWS S3 data synchronization
- **Python AI**: GPT-4 summarization via pybind11

---

## Phase 6: Enable 8 Excluded Modules

Once Qt5 is working, re-enable these modules in CMakeLists.txt:

### Qt5-Dependent Modules

1. **ScientificCalculatorDialog.cpp** - Full GUI scientific calculator
2. **SymEngine.cpp** - Symbolic math engine GUI
3. **SymEngineAllocator.cpp** - Memory management for SymEngine
4. **ForceModule.cpp** - Force calculation GUI (formerly for.cpp)
5. **InputModule.cpp** - Input handling GUI (formerly in.cpp)

### OpenGL-Dependent Modules

6. **FluidSolver.cpp** - Fluid dynamics visualization (requires GLEW)
7. **SIMPlugin.cpp** - Simulation plugin with 3D rendering (requires GLEW)

### Self-Healing Module

8. **HydrogenResonanceUQFFModule.cpp** - Requires internal class repairs

### Build After Enabling

```cmake
# Add to CMakeLists.txt after Source2_HeadProgram
add_executable(ScientificCalculator ScientificCalculatorDialog.cpp)
target_link_libraries(ScientificCalculator PRIVATE Qt5::Widgets Qt5::Core)

# ... (similar for other 7 modules)
```

---

## Phase 7: Integration Testing

### Test Matrix

| Component | Test | Expected Result |
|-----------|------|-----------------|
| **MAIN_1_CoAnQi** | Menu option 1-13 | All options functional ✅ |
| **Source2_HeadProgram** | Launch GUI | 21 browser windows open |
| **VTK Charts** | Plot test data | 3D scatter plot renders |
| **NASA API** | Fetch APOD | Daily image displayed |
| **Wolfram WSTP** | Option 9-11 | Kernel communication works |
| **Grok API** | Option 12 | AI diagnostic response |
| **Voice Input** | Speak command | PocketSphinx recognizes |
| **AWS S3** | Sync data | Cloud upload successful |
| **Python AI** | GPT-4 summary | pybind11 executes Python |

---

## Timeline Estimate

| Phase | Task | Duration | Dependencies |
|-------|------|----------|--------------|
| **1** | Install Qt5 via vcpkg | 30-60 min | Internet, disk space |
| **2** | Install VTK/OpenCV/AWS SDK | 60-120 min | Internet, disk space |
| **3** | Modify CMakeLists.txt | 15-30 min | None |
| **4** | Update source2.cpp API keys | 5-10 min | API key acquisition |
| **5** | Build Source2_HeadProgram | 10-20 min | All dependencies |
| **6** | Enable 8 excluded modules | 30-60 min | Qt5, OpenGL/GLEW |
| **7** | Integration testing | 60-120 min | All components built |

**Total Estimate**: 3-6 hours (depending on download speeds and complexity)

---

## Rollback Plan

If Qt5 GUI expansion causes issues:

```powershell
# Return to working MAIN_1_CoAnQi state
git checkout savepoint-wolfram-qt6-grok-dec1-2025

# Rebuild core calculator only
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64 -DUSE_EMBEDDED_WOLFRAM=ON
cmake --build build_msvc --config Release --target MAIN_1_CoAnQi -j 8
```

---

## Success Criteria

- [ ] source2.cpp compiles without errors
- [ ] Source2_HeadProgram.exe launches GUI
- [ ] 21 browser windows open simultaneously
- [ ] VTK visualizations render correctly
- [ ] NASA/MAST APIs fetch data
- [ ] Voice commands work (if PocketSphinx installed)
- [ ] AWS cloud sync functional (if configured)
- [ ] Python GPT-4 integration active (if pybind11 configured)
- [ ] All 8 excluded modules compile and run
- [ ] MAIN_1_CoAnQi.exe still functional (no regression)

---

**Next Steps:**
1. Test MAIN_1_CoAnQi.exe with current Qt6 setup
2. Begin Phase 1 (Qt5 installation) if test passes
3. Proceed sequentially through phases 2-7
4. Document any issues encountered for future reference

**Author:** GitHub Copilot  
**Date:** December 1, 2025  
**Version:** 1.0-Initial
