# Star-Magic UQFF Project
# Unified Quantum Field Force Calculations
# Author: Daniel T. Murphy
cmake_minimum_required(VERSION 3.10)

# ========== VCPKG TOOLCHAIN INTEGRATION ==========
# CRITICAL: Must be set BEFORE project() command for find_package() to work
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND EXISTS "C:/vcpkg/scripts/buildsystems/vcpkg.cmake")
    set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")
endif()

# Project name and programming language
# Build Status: February 21, 2026 - v4.2.1 Poseidon TaskBot Integration
# Session: Added Poseidon offline physics maintenance bot to VR/VM Backend
# IPC Layer: Cross-platform (Windows + Linux/macOS) with NamedPipeChannel + SharedMemoryChannel + gRPC
# v4.2 Features:
#   Poseidon TaskBot: Offline physics maintenance (read/write/compare/validate/update)
#   Self-Expand: onRegisterTerm() - dynamic physics terms from VR/external
#   Self-Update: onUpdateParameter() - runtime κ, [SSq], β_i tuning
#   Self-Simulate: startSimulation() - time evolution with VR streaming
# Previous: v3.1 Self-Expanding Physics Backend (Feb 21, 2026)
# Framework Status: 99.9% solvability (Grok 4 analysis Sept 14-21, 2025)
# Calibrated Variables: κ=0.0005/day, [SSq]=0.57, H_SCm=0.99, U_UA=0.0001, k_η=10⁻¹¹³, β_i=0.603
# PhysicsTerm Classes: 6,688+ registered
# Compiler: MSVC 19.44.35219.0 (Visual Studio 2022)
# Standard: C++20 strict compliance
# Executable: 1.43 MB (9,217,024 -> 1,429,504 bytes, 15.51% UPX compression)
# Last Update: February 21, 2026 @ v4.2.1 Poseidon TaskBot
project(StarMagic_UQFF LANGUAGES CXX)

# ========== BUILD OPTIONS ==========
# CMake options for profiling, benchmarking, and WebAssembly builds
option(ENABLE_PROFILING "Enable performance profiling with gprof" OFF)
option(ENABLE_BENCHMARK "Enable Google Benchmark tests" OFF)
option(ENABLE_WASM "Build for WebAssembly using Emscripten" OFF)
option(ENABLE_TESTS "Enable unit testing framework" ON)
option(ENABLE_PLUGINS "Enable plugin system compilation" ON)

# ========== FORCE MSVC ON WINDOWS ==========
# WSTP libraries (wstp64i4.lib, wstp64i4.dll) are MSVC-compiled ONLY
# MinGW CANNOT link with these libraries - binary incompatibility
if(WIN32 AND NOT MSVC)
    message(FATAL_ERROR "\n"
        "============================================================\n"
        "ERROR: This project REQUIRES MSVC on Windows!\n"
        "\n"
        "REASON: Wolfram WSTP libraries are MSVC-compiled only.\n"
        "        MinGW cannot link with MSVC binaries.\n"
        "\n"
        "SOLUTION: Use Visual Studio generator:\n"
        "  cmake -G \"Visual Studio 17 2022\" -A x64 -S . -B build\n"
        "============================================================\n")
endif()

# Global include directories for all targets
include_directories(C:/vcpkg/installed/x64-windows/include)

# ========== OPTION A FULL RECOVERY: Qt6 + ANTLR4 + SymEngine ==========
# Using existing Qt 6.10.0 MSVC 2022 installation
set(CMAKE_PREFIX_PATH "C:/Qt/6.10.0/msvc2022_64" ${CMAKE_PREFIX_PATH})
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Find Qt6 (upgraded from Qt5)
find_package(Qt6 COMPONENTS Core Widgets Network WebEngineWidgets QUIET)
if(Qt6_FOUND)
    message(STATUS "Qt6 found - GUI components enabled (Option A Full Recovery)")
    add_definitions(-DHAVE_QT6)
    if(TARGET Qt6::WebEngineWidgets)
        message(STATUS "Qt6::WebEngineWidgets found - source2 web browser enabled")
    else()
        message(STATUS "Qt6::WebEngineWidgets NOT found - source2 will have limited functionality")
    endif()
else()
    # Fallback to Qt5 if Qt6 not found
    find_package(Qt5 COMPONENTS Core Widgets WebEngineWidgets PrintSupport QUIET)
    if(Qt5_FOUND)
        message(STATUS "Qt5 found - GUI components enabled")
        add_definitions(-DHAVE_QT5)
    else()
        message(STATUS "Qt not found - GUI components will be commented out")
    endif()
endif()

# Find VTK for scientific visualization (source2)
find_package(VTK CONFIG QUIET)
if(VTK_FOUND)
    message(STATUS "VTK found at version ${VTK_VERSION} - source2 visualization enabled")
else()
    message(STATUS "VTK NOT found - source2 will compile with -DNO_VTK")
endif()

# Find CURL for HTTP requests (source2)
find_package(CURL CONFIG QUIET)
if(CURL_FOUND)
    message(STATUS "CURL found - source2 web API access enabled")
else()
    message(STATUS "CURL NOT found - source2 will compile with -DNO_CURL")
endif()

# Find SQLite3 for local database (source2)
find_package(unofficial-sqlite3 CONFIG QUIET)
if(unofficial-sqlite3_FOUND)
    message(STATUS "SQLite3 found - source2 database caching enabled")
else()
    message(STATUS "SQLite3 NOT found - source2 will compile with -DNO_SQLITE")
endif()

# Find OpenCV for computer vision (source2)
find_package(OpenCV CONFIG QUIET)
if(OpenCV_FOUND)
    message(STATUS "OpenCV found at version ${OpenCV_VERSION} - source2 vision processing enabled")
else()
    message(STATUS "OpenCV NOT found - source2 will compile with -DNO_OPENCV")
endif()

# Find AWS SDK for cloud services (source2)
find_package(AWSSDK CONFIG COMPONENTS cognito-idp s3 QUIET)
if(AWSSDK_FOUND)
    message(STATUS "AWS SDK found - source2 cloud sync enabled")
else()
    message(STATUS "AWS SDK NOT found - source2 will compile with -DNO_AWS")
endif()

# Find ANTLR4 runtime (for parser code patterns)
find_package(antlr4-runtime CONFIG QUIET)
if(antlr4-runtime_FOUND)
    message(STATUS "ANTLR4 runtime found - parser components enabled")
    add_definitions(-DHAVE_ANTLR4)
endif()

# Find SymEngine (for symbolic math patterns)
find_package(SymEngine CONFIG QUIET)
if(SymEngine_FOUND)
    message(STATUS "SymEngine found - symbolic math components enabled")
    add_definitions(-DHAVE_SYMENGINE)
endif()

# Add calculator_advanced subdirectory (Grok Thread Iterations #30-32)
# Advanced scientific calculator with ANTLR4 parser, SymEngine backend, 11 UQFF equations
if(Qt6_FOUND AND antlr4-runtime_FOUND AND SymEngine_FOUND)
    message(STATUS "calculator_advanced: All dependencies found - building Tab 22")
    add_subdirectory(calculator_advanced)
else()
    if(NOT Qt6_FOUND)
        message(STATUS "calculator_advanced: Qt6 not found - skipping")
    endif()
    if(NOT antlr4-runtime_FOUND)
        message(STATUS "calculator_advanced: ANTLR4 not found - skipping")
    endif()
    if(NOT SymEngine_FOUND)
        message(STATUS "calculator_advanced: SymEngine not found - skipping")
    endif()
endif()

# Set C++ standard to C++20 for maximum VR/VM functionality
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler-specific flags for MSVC
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /std:c++20 /permissive- /Zc:__cplusplus /GR- /arch:AVX2 /wd4100 /wd4244 /wd4267 /wd4305 /wd4996 /wd4005")
    # Warning suppressions: /wd4100=unreferenced params, /wd4244=conversion loss, /wd4267=size_t conversion, /wd4305=truncation, /wd4996=deprecated POSIX, /wd4005=macro redefinition
    set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od")
    set(CMAKE_CXX_FLAGS_RELEASE "/Os /GL /Gw /Gy /Zc:inline /GF /Oi")  # /Os = Favor small code, /GL = Whole Program Optimization, /Gw = Global data optimization, /Gy = Function-level linking, /Zc:inline = Remove unreferenced inline, /GF = String pooling, /Oi = Intrinsics
    
    # RECOMMENDED FINAL RELEASE CONFIGURATION "Release-MaxCompress"
    # Configuration verified against VS2022 Property Pages best practices
    
    # STEP 1 of 8: Whole Program Optimization Control Points
    # /GL (Whole Program Optimization) + Linker → /LTCG + /OPT:REF + /OPT:ICF
    # Applies to: All 175 .cpp files (source*.cpp, MAIN_1*.cpp, Core/Modules/*.cpp)
    # Effect: Code + identical COMDAT folding across entire codebase
    
    # STEP 2 of 8: Favor Size Optimization
    # /Os instead of /O2 or /Ox
    # Effect: Prioritizes code size over speed (smaller binaries, better cache utilization)
    # VS2022 Property Pages: C/C++ → Optimization → Favor Size or Speed = Favor small code (/Os)
    
    # STEP 3 of 8: Global Data Optimization + Inline Function Elimination
    # /Gw (Global Data Optimization) - separate COMDAT sections for each global variable
    # /Zc:inline (Remove Unreferenced Inline) - eliminates unused inline/template functions
    # Effect: CRITICAL for this project - removes duplicate constexpr/static data across 175+ translation units
    # Expected savings: 500KB-2MB (massive static const duplication in physics modules)
    # VS2022 Property Pages: C/C++ → Optimization → Enable Global Data Optimization = Yes (/Gw)
    
    # STEP 7 of 8: Standards Conformance for Better LTCG Optimization
    # /permissive- - strict C++ standards conformance mode
    # /Zc:__cplusplus - report correct __cplusplus macro value (201703L for C++17)
    # Effect: Cleaner codegen enables LTCG to identify more merging opportunities
    # Benefit: More aggressive template instantiation folding, better type deduction
    # Minor size savings: 10-50KB (better template deduplication)
    # VS2022 Property Pages: C/C++ → Language → Conformance mode = Yes (/permissive-)
    
    # ADDITIONAL RELEASE-MAXCOMPRESS FLAGS:
    # /GF - String Pooling: merges identical string literals across compilation units
    # /GR- - Disable RTTI: removes Run-Time Type Information (typeid, dynamic_cast)
    # Effect: /GF saves 20-100KB (duplicate strings), /GR- saves 50-200KB (RTTI tables)
    # VS2022 Property Pages:
    #   C/C++ → Code Generation → Enable String Pooling = Yes (/GF)
    #   C/C++ → Language → Enable Run-Time Type Information = No (/GR-)
    
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG /OPT:REF /OPT:ICF /MERGE:.rdata=.text /SUBSYSTEM:CONSOLE /MANIFESTUAC:level='asInvoker' /DEBUG:NONE /PDBALTPATH:none /INCREMENTAL:NO")
    # Explicit LTCG ensures link-time code generation with /GL compiler flag
    set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /LTCG /OPT:REF /OPT:ICF /MERGE:.rdata=.text /DEBUG:NONE")
    set(CMAKE_STATIC_LINKER_FLAGS_RELEASE "${CMAKE_STATIC_LINKER_FLAGS_RELEASE} /LTCG")
    
    # STEP 4 of 8: Section Merging for Reduced PE Overhead
    # /MERGE:.rdata=.text - merges read-only data section into code section
    # Effect: Reduces PE header overhead, fewer sections = smaller binary
    # Benefit: .rdata contains string literals, vtables, const data - merging with .text saves section alignment padding
    # Expected savings: 50-200KB (section alignment is typically 4KB, each merged section eliminates padding)
    # VS2022 Property Pages: Linker → Advanced → Merge Sections = .rdata=.text
    
    # STEP 5 of 8: Subsystem Optimization + Minimal Manifest
    # /SUBSYSTEM:CONSOLE - required for MAIN_1_CoAnQi.cpp (uses std::cout)
    # /MANIFESTUAC:level='asInvoker' - minimal UAC manifest (no elevation required)
    # Effect: Smaller embedded manifest, no admin prompt overhead
    # Note: /SUBSYSTEM:WINDOWS would save ~10KB but breaks console I/O (std::cout becomes silent)
    # VS2022 Property Pages: 
    #   Linker → System → SubSystem = Console (/SUBSYSTEM:CONSOLE)
    #   Linker → Manifest File → UAC Execution Level = asInvoker
    
    # STEP 6 of 8: Strip Debug Information for Maximum Size Reduction
    # /DEBUG:NONE - completely disables PDB generation (no debug symbols)
    # /PDBALTPATH:none - removes PDB path from PE header
    # Effect: Eliminates multi-MB PDB overhead, removes debug directory from executable
    # Expected savings: 100-500KB in .exe (CodeView debug info removed), 5-50MB .pdb file not created
    # Trade-off: No stack traces on crash, no source-level debugging in production builds
    # VS2022 Property Pages: Linker → Debugging → Generate Debug Info = No (/DEBUG:NONE)
    
    # VS2022 Property Pages equivalents:
    # C/C++ → Optimization → Whole Program Optimization = Yes (/GL)
    # Linker → Optimization → Enable Whole-Program Optimization = Yes (/LTCG)
    # Linker → Optimization → References = Yes (/OPT:REF)
    # Linker → Optimization → Enable COMDAT Folding = Yes (/OPT:ICF)
endif()

# Force synchronized PDB writes for parallel compilation (prevents C1041 error)
if(MSVC)
    add_compile_options(/FS)
endif()

# Export compile commands for IDE integration
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ========== PROFILING CONFIGURATION ==========
# Enable performance profiling with gprof on Linux/Unix systems
if(ENABLE_PROFILING)
    if(UNIX AND NOT APPLE)
        message(STATUS "Profiling ENABLED - Adding -pg flags for gprof")
        
        # Add -pg flag for all targets (generates gmon.out)
        add_compile_options(-pg)
        add_link_options(-pg)
        
        # Apply to specific targets
        if(TARGET MAIN_1_CoAnQi)
            target_compile_options(MAIN_1_CoAnQi PRIVATE -pg)
            target_link_options(MAIN_1_CoAnQi PRIVATE -pg)
        endif()
        
        if(TARGET source2)
            target_compile_options(source2 PRIVATE -pg)
            target_link_options(source2 PRIVATE -pg)
        endif()
        
        message(STATUS "Profiling: Run executable, then analyze with 'gprof ./MAIN_1_CoAnQi gmon.out > analysis.txt'")
    else()
        message(STATUS "Profiling: Only supported on Linux (gprof). Use Visual Studio Profiler on Windows.")
    endif()
endif()

# ========== BENCHMARK CONFIGURATION ==========
# Enable Google Benchmark for performance testing
if(ENABLE_BENCHMARK)
    find_package(benchmark QUIET)
    if(benchmark_FOUND)
        message(STATUS "Google Benchmark FOUND - Enabling benchmark tests")
        
        # Define benchmark flag
        add_definitions(-DENABLE_BENCHMARK)
        
        # Link benchmark library to targets
        if(TARGET MAIN_1_CoAnQi)
            target_link_libraries(MAIN_1_CoAnQi PRIVATE benchmark::benchmark)
        endif()
        
        if(TARGET source2)
            target_link_libraries(source2 PRIVATE benchmark::benchmark)
        endif()
        
        message(STATUS "Benchmark: Link with benchmark::benchmark library")
    else()
        message(STATUS "Google Benchmark NOT found - Install with: vcpkg install benchmark")
        message(STATUS "Benchmark tests will be DISABLED")
    endif()
endif()

# ========== WEBASSEMBLY CONFIGURATION ==========
# Configure build for WebAssembly using Emscripten compiler
if(ENABLE_WASM)
    if(EMSCRIPTEN)
        message(STATUS "WebAssembly build ENABLED (Emscripten detected)")
        
        # Change executable suffix to .html for browser loading
        set(CMAKE_EXECUTABLE_SUFFIX ".html")
        
        # Add Emscripten-specific flags
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s WASM=1 -s ALLOW_MEMORY_GROWTH=1")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s EXPORTED_RUNTIME_METHODS=['ccall','cwrap']")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --bind")  # Enable Embind for C++/JS interop
        
        # Optimize for size in WASM builds
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Os")
        
        message(STATUS "WebAssembly: Executable suffix = .html")
        message(STATUS "WebAssembly: Memory growth enabled")
        message(STATUS "WebAssembly: ccall/cwrap runtime methods exported")
    else()
        message(WARNING "ENABLE_WASM is ON but Emscripten not detected!")
        message(WARNING "Use: emcmake cmake -B build-wasm -DENABLE_WASM=ON")
    endif()
endif()

# ========== WOLFRAM ENGINE INTEGRATION (OPTIONAL) ==========
option(USE_EMBEDDED_WOLFRAM "Enable embedded Wolfram kernel via WSTP" ON)

if(USE_EMBEDDED_WOLFRAM)
    # WSTP-Specific Settings: C:\Program Files\Wolfram Research\Wolfram Engine\14.3\SystemFiles\Links\WSTP\DeveloperKit\Windows-x86-64\CompilerAdditions\wstp
    set(WSTP_DIR "C:/Program Files/Wolfram Research/Wolfram Engine/14.3/SystemFiles/Links/WSTP/DeveloperKit/Windows-x86-64/CompilerAdditions")
    
    if(EXISTS "${WSTP_DIR}")
        # C/C++ → General → Additional Include Directories
        include_directories("${WSTP_DIR}")
        
        # Linker → General → Additional Library Directories
        link_directories("${WSTP_DIR}")
        
        # Preprocessor → Preprocessor Definitions: USE_EMBEDDED_WOLFRAM, _CRT_SECURE_NO_WARNINGS, WIN64
        add_definitions(-DUSE_EMBEDDED_WOLFRAM -D_CRT_SECURE_NO_WARNINGS -DWIN64)
        
        message(STATUS "Wolfram Engine integration ENABLED (WSTP found at ${WSTP_DIR})")
    else()
        message(WARNING "USE_EMBEDDED_WOLFRAM is ON but WSTP_DIR not found: ${WSTP_DIR}")
        message(WARNING "Wolfram integration will be DISABLED. Install Wolfram Engine or adjust WSTP_DIR path.")
    endif()
else()
    message(STATUS "Wolfram Engine integration DISABLED (use -DUSE_EMBEDDED_WOLFRAM=ON to enable)")
endif()

# ========== MAIN EXECUTABLES ==========
# Core mathematical framework executables
add_executable(MAIN_1 MAIN_1.cpp)

# MAIN_1_CoAnQi with optional Wolfram integration + Grok AI API
if(USE_EMBEDDED_WOLFRAM)
    # source174 is included via source176, don't add it separately to avoid duplicate definitions
    # source178_grok_api.cpp is included directly in MAIN_1_CoAnQi.cpp (line 206), NOT compiled separately
    # Phase 4 Wolfram extraction: Files included directly in MAIN_1_CoAnQi.cpp, not compiled separately
    add_executable(MAIN_1_CoAnQi 
        MAIN_1_CoAnQi.cpp
        ipc/uqff_ipc.cpp
        Core/Modules/WolframFieldUnityModule.cpp
    )
    
    # Disable AUTOMOC for MAIN_1_CoAnQi - it's a pure C++ physics calculator, not a Qt GUI app
    # Q_OBJECT macros in commented sections should not trigger AutoMoc
    set_target_properties(MAIN_1_CoAnQi PROPERTIES AUTOMOC OFF AUTOUIC OFF AUTORCC OFF)
    
    # Linker → Input → Additional Dependencies: wstp64i4.lib + Windows system libraries
    target_link_libraries(MAIN_1_CoAnQi 
        wstp64i4
        kernel32
        user32
        advapi32
        shell32
    )
    
    # Add Qt5 libraries if available (required for Grok API QNetworkAccessManager)
    if(Qt5_FOUND)
        target_link_libraries(MAIN_1_CoAnQi
            Qt5::Core
            Qt5::Widgets
            Qt5::WebEngineWidgets
            Qt5::PrintSupport
            Qt5::Network
        )
    endif()
    
    # Add Qt6 libraries if available (required for Grok API QNetworkAccessManager)
    if(Qt6_FOUND)
        target_link_libraries(MAIN_1_CoAnQi
            Qt6::Core
            Qt6::Widgets
            Qt6::Network
        )
    endif()
    
    # Physics term extractor (Wolfram-powered)
    add_executable(PhysicsTermExtractor physics_term_extractor_main.cpp)
    target_link_libraries(PhysicsTermExtractor 
        wstp64i4
        kernel32
        user32
        advapi32
        shell32
    )
    
    message(STATUS "MAIN_1_CoAnQi compiled WITH Wolfram WSTP support")
    message(STATUS "PhysicsTermExtractor compiled WITH Wolfram WSTP support")
else()
    add_executable(MAIN_1_CoAnQi 
        MAIN_1_CoAnQi.cpp
        ipc/uqff_ipc.cpp
        Core/Modules/WolframFieldUnityModule.cpp
    )
    
    # Disable AUTOMOC for MAIN_1_CoAnQi - pure C++ physics calculator
    set_target_properties(MAIN_1_CoAnQi PROPERTIES AUTOMOC OFF AUTOUIC OFF AUTORCC OFF)
    
    message(STATUS "MAIN_1_CoAnQi compiled WITHOUT Wolfram support")
endif()

# STEP 8 of 8: UPX Ultimate Packer (Post-Build Compression)
# UPX 5.0.2 (Latest as of Nov 22, 2025) - Download from https://github.com/upx/upx/releases
# Compresses final executable using UPX with maximum compression settings
# Effect: Dramatic size reduction (50-70% compression ratio typical)
# Target: Get MAIN_1_CoAnQi.exe under 3 MB (from ~5-8 MB uncompressed)
# Requirements: UPX 5.0.2+ installed at C:/Tools/upx/upx.exe
# Trade-offs: Slower startup (decompression overhead ~50-200ms), incompatible with some antivirus
# VS2022 Property Pages: Build Events → Post-Build Event → Command Line
# NOTE: UPX compression DISABLED to avoid "already packed" errors during rebuilds
#find_program(UPX_EXECUTABLE upx PATHS "C:/Tools/upx" "C:/Tools/UPX" "C:/Program Files/upx" ENV PATH DOC "UPX 5.0.2 executable packer")
#if(UPX_EXECUTABLE)
#    add_custom_command(TARGET MAIN_1_CoAnQi POST_BUILD
#        COMMAND "${UPX_EXECUTABLE}" --best --lzma "$<TARGET_FILE:MAIN_1_CoAnQi>"
#        COMMENT "STEP 8: Compressing MAIN_1_CoAnQi.exe with UPX 5.0.2 (--best --lzma)"
#        VERBATIM
#    )
#    message(STATUS "UPX compression ENABLED for MAIN_1_CoAnQi (found at ${UPX_EXECUTABLE})")
#else()
#    message(WARNING "UPX not found - executable will not be compressed.")
#    message(WARNING "Download UPX 5.0.2: https://github.com/upx/upx/releases/download/v5.0.2/upx-5.0.2-win64.zip")
#    message(WARNING "Extract to C:/Tools/upx/ and add to PATH")
#endif()
message(STATUS "UPX compression DISABLED for MAIN_1_CoAnQi during development")

# OpenSSL DLL Auto-Copy for Qt6 TLS/HTTPS Support (Grok API Integration)
# Automatically copies OpenSSL DLLs after building MAIN_1_CoAnQi
# Required for: Grok AI integration (Option 13), Qt6 QNetworkAccessManager HTTPS requests
# Installation: winget install ShiningLight.OpenSSL.Light
if(EXISTS "C:/Program Files/OpenSSL-Win64/bin/libssl-3-x64.dll")
    add_custom_command(TARGET MAIN_1_CoAnQi POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "C:/Program Files/OpenSSL-Win64/bin/libssl-3-x64.dll"
            "C:/Program Files/OpenSSL-Win64/bin/libcrypto-3-x64.dll"
            $<TARGET_FILE_DIR:MAIN_1_CoAnQi>
        COMMENT "Copying OpenSSL DLLs for Qt6 TLS/HTTPS support (Grok API)"
    )
    message(STATUS "OpenSSL DLL auto-copy ENABLED for MAIN_1_CoAnQi")
else()
    message(STATUS "OpenSSL not found at 'C:/Program Files/OpenSSL-Win64/' - Grok API will require manual DLL setup")
    message(STATUS "Install via: winget install ShiningLight.OpenSSL.Light")
endif()

# ========== IPC LIBRARY (Phase 1 - VR Architecture) ==========
# Shared memory and gRPC channels for UQFF inter-process communication
# Phase 2: Added PhysicsService for backend service mode
# Phase 2 Fix: Added PythonBridge for CondensedPhysics.py/QCalc.py integration
set(IPC_SOURCES
    ipc/uqff_ipc.cpp
    ipc/physics_service.cpp
    ipc/python_bridge.cpp
)
set(IPC_HEADERS
    ipc/uqff_ipc.h
    ipc/physics_service.h
    ipc/python_bridge.h
)

# Source2: HEAD PROGRAM - Primary User Base Platform (Qt6 + optional VTK, AWS SDK, OpenCV, etc.)
# Status: ENABLED with Qt6 - optional dependencies compile conditionally
# To install missing dependencies via vcpkg (x64-windows for MSVC):
#   vcpkg install qt6-webengine:x64-windows
#   vcpkg install vtk:x64-windows opencv:x64-windows
#   vcpkg install curl:x64-windows sqlite3:x64-windows
#   vcpkg install aws-sdk-cpp[cognito-idp,s3]:x64-windows
if(Qt6_FOUND)
    add_executable(Source2 
        source2.cpp 
        source2_mainwindow.h 
        UQFFResultsWidget.h 
        equation_renderer.h 
        UQFF_ScientificCalculator.cpp
        ${IPC_SOURCES} 
        ${IPC_HEADERS}
        ${CMAKE_SOURCE_DIR}/resources/Star-Magic.rc
    )
    
    # Enable Qt AUTOMOC/AUTOUIC/AUTORCC for Source2 (Q_OBJECT in header file)
    set_target_properties(Source2 PROPERTIES 
        AUTOMOC ON 
        AUTOUIC ON 
        AUTORCC ON
        # Use dynamic CRT (required for vcpkg AWS SDK compatibility)
        MSVC_RUNTIME_LIBRARY "MultiThreadedDLL$<$<CONFIG:Debug>:Debug>"
    )
    
    # IPC include directory (Phase 1 - VR Architecture)
    target_include_directories(Source2 PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/ipc
        ${CMAKE_CURRENT_SOURCE_DIR}/resources
    )
    
    # Base Qt6 dependencies (always required)
    target_link_libraries(Source2 PRIVATE 
        Qt6::Core 
        Qt6::Widgets 
        Qt6::Network
    )
    
    # Link calculator_advanced library (Tab 22)
    if(TARGET calculator_advanced)
        target_link_libraries(Source2 PRIVATE calculator_advanced)
        target_include_directories(Source2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
        message(STATUS "Source2: Linked with calculator_advanced (Tab 22)")
    endif()
    
    # Qt6 WebEngineWidgets (conditional)
    if(TARGET Qt6::WebEngineWidgets)
        target_link_libraries(Source2 PRIVATE Qt6::WebEngineWidgets)
    else()
        target_compile_definitions(Source2 PRIVATE NO_WEBENGINE)
        message(STATUS "Source2: Compiling without web browser (install qt6-webengine)")
    endif()
    
    # VTK (conditional)
    if(VTK_FOUND)
        target_link_libraries(Source2 PRIVATE ${VTK_LIBRARIES})
        target_include_directories(Source2 PRIVATE ${VTK_INCLUDE_DIRS})
    else()
        target_compile_definitions(Source2 PRIVATE NO_VTK)
    endif()
    
    # CURL (conditional)
    if(CURL_FOUND)
        target_link_libraries(Source2 PRIVATE CURL::libcurl)
    else()
        target_compile_definitions(Source2 PRIVATE NO_CURL)
    endif()
    
    # SQLite3 (conditional)
    if(unofficial-sqlite3_FOUND)
        target_link_libraries(Source2 PRIVATE unofficial::sqlite3::sqlite3)
    else()
        target_compile_definitions(Source2 PRIVATE NO_SQLITE)
    endif()
    
    # OpenCV (conditional)
    if(OpenCV_FOUND)
        target_link_libraries(Source2 PRIVATE ${OpenCV_LIBS})
    else()
        target_compile_definitions(Source2 PRIVATE NO_OPENCV)
    endif()
    
    # CRITICAL: Ensure consistent MSVC runtime library for all dependencies
    # This prevents LNK2001 iostream/STL symbol conflicts
    set_target_properties(Source2 PROPERTIES 
        MSVC_RUNTIME_LIBRARY "MultiThreadedDLL$<$<CONFIG:Debug>:Debug>"
    )
    
    # Explicit linker flags for MSVC runtime (prevents CRT mismatch)
    if(MSVC)
        target_link_options(Source2 PRIVATE 
            $<$<CONFIG:Release>:/NODEFAULTLIB:LIBCMT>     # Exclude static CRT
            $<$<CONFIG:Release>:/NODEFAULTLIB:LIBCMTD>    # Exclude static debug CRT
            $<$<CONFIG:Debug>:/NODEFAULTLIB:LIBCMT>
            $<$<CONFIG:Debug>:/NODEFAULTLIB:MSVCRT>       # Exclude release CRT in debug
        )
        # Explicitly link MSVC dynamic runtime
        target_link_libraries(Source2 PRIVATE 
            $<$<CONFIG:Release>:msvcrt.lib>
            $<$<CONFIG:Debug>:msvcrtd.lib>
        )
    endif()
    
    # AWS SDK - TEMPORARILY DISABLED for CRT debugging
    # TODO: Re-enable after resolving iostream symbol conflicts
    target_compile_definitions(Source2 PRIVATE NO_AWS)
    message(STATUS "Source2: AWS SDK temporarily disabled (CRT debugging)")
    
    # # AWS SDK (conditional) - RE-ENABLED with proper CRT configuration
    # # Solution: Ensure consistent /MD runtime and proper link order
    # if(AWSSDK_FOUND)
    #     # Link AWS SDK libraries in correct order (dependencies first)
    #     target_link_libraries(Source2 PRIVATE 
    #         aws-cpp-sdk-s3
    #         aws-cpp-sdk-cognito-idp
    #         aws-cpp-sdk-core
    #         ${AWSSDK_LINK_LIBRARIES}
    #     )
    #     # Add AWS SDK include directories
    #     if(AWSSDK_INCLUDE_DIRS)
    #         target_include_directories(Source2 PRIVATE ${AWSSDK_INCLUDE_DIRS})
    #     endif()
    #     message(STATUS "Source2: AWS SDK enabled (S3 + Cognito)")
    # else()
    #     target_compile_definitions(Source2 PRIVATE NO_AWS)
    #     message(STATUS "Source2: AWS SDK disabled (not found)")
    # endif()
    
    # Python (pybind11) - Phase 2 Fix: HYBRID APPROACH
    # - Source2 GUI code uses NO_PYTHON to avoid Qt+pybind conflicts
    # - IPC layer (python_bridge.cpp) uses PYTHON_BRIDGE_ENABLED for --service mode
    # - CondensedPhysics.py (81K lines) and QCalc.py (9K lines) work via IPC
    find_package(Python3 COMPONENTS Interpreter Development QUIET)
    find_package(pybind11 CONFIG QUIET)
    
    if(Python3_FOUND AND pybind11_FOUND)
        # Link pybind11 for python_bridge.cpp in IPC layer
        target_link_libraries(Source2 PRIVATE 
            Python3::Python
            pybind11::embed
        )
        target_include_directories(Source2 PRIVATE ${Python3_INCLUDE_DIRS})
        # NO_PYTHON: Skips Qt-conflicting code in source2.cpp
        # PYTHON_BRIDGE_ENABLED: Enables python_bridge.cpp for --service mode
        target_compile_definitions(Source2 PRIVATE NO_PYTHON PYTHON_BRIDGE_ENABLED)
        message(STATUS "Source2: Python via IPC (--service mode uses CondensedPhysics.py/QCalc.py)")
    else()
        target_compile_definitions(Source2 PRIVATE NO_PYTHON)
        message(STATUS "Source2: Python disabled (pybind11 not found)")
    endif()
    
    # Qalculate, PocketSphinx, MPI, LLVM - Not installed, always disabled
    target_compile_definitions(Source2 PRIVATE NO_QALCULATE NO_POCKETSPHINX NO_MPI NO_LLVM)
    
    # Post-build: Deploy Python wrappers and OpenSSL DLLs
    if(WIN32)
        add_custom_command(TARGET Source2 POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E echo "Deploying Python wrappers and dependencies for Source2..."
            # TEMPORARILY DISABLED: MAIN_1_CoAnQi.exe copy (has build errors)
            # Copy MAIN_1_CoAnQi.exe (Tab 1 - Embedded Calculator Terminal)
            # COMMAND ${CMAKE_COMMAND} -E copy_if_different
            #     "$<TARGET_FILE:MAIN_1_CoAnQi>"
            #     "$<TARGET_FILE_DIR:Source2>/MAIN_1_CoAnQi.exe"
            # Copy Python wrappers (UQFF, Symbolic Math, API Fetch)
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_SOURCE_DIR}/CoAnQi_Wrapper.py"
                "$<TARGET_FILE_DIR:Source2>/CoAnQi_Wrapper.py"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_SOURCE_DIR}/SymbolicMath_Wrapper.py"
                "$<TARGET_FILE_DIR:Source2>/SymbolicMath_Wrapper.py"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_SOURCE_DIR}/APIFetch.py"
                "$<TARGET_FILE_DIR:Source2>/APIFetch.py"
            # Copy QCalc.py (Tab 2 - 8 UQFF Master Equations)
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_SOURCE_DIR}/QCalc.py"
                "$<TARGET_FILE_DIR:Source2>/QCalc.py"
            # Copy CondensedPhysics.py (Tab 5 - General Model Solver)
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_SOURCE_DIR}/CondensedPhysics.py"
                "$<TARGET_FILE_DIR:Source2>/CondensedPhysics.py"
            # Copy ValidationData.py (Tab 7 - Validation wrapper for 6000+ URLs)
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_SOURCE_DIR}/ValidationData.py"
                "$<TARGET_FILE_DIR:Source2>/ValidationData.py"
            # Copy GrokAPI.py (Tab 7 - xAI Grok API wrapper for SuperGrok4)
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_SOURCE_DIR}/GrokAPI.py"
                "$<TARGET_FILE_DIR:Source2>/GrokAPI.py"
            # Copy OllamaAPI.py (Tab 6 - Ollama API wrapper for CoAnQi_bot)
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_SOURCE_DIR}/OllamaAPI.py"
                "$<TARGET_FILE_DIR:Source2>/OllamaAPI.py"
            # Copy CondensedPhysics dependency files (data modules)
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_SOURCE_DIR}/CondensedPhysics_InputData.py"
                "$<TARGET_FILE_DIR:Source2>/CondensedPhysics_InputData.py"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_SOURCE_DIR}/CondensedPhysics_OutputData.py"
                "$<TARGET_FILE_DIR:Source2>/CondensedPhysics_OutputData.py"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_SOURCE_DIR}/CondensedPhysics_Validation.py"
                "$<TARGET_FILE_DIR:Source2>/CondensedPhysics_Validation.py"
            # Copy OpenSSL DLLs (required for MAIN_1_CoAnQi.exe Grok integration)
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "C:/Program Files/OpenSSL-Win64/bin/libssl-3-x64.dll"
                "C:/Program Files/OpenSSL-Win64/bin/libcrypto-3-x64.dll"
                "$<TARGET_FILE_DIR:Source2>/"
            # Copy Visual C++ Runtime DLLs (required for MAIN_1_CoAnQi.exe)
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "C:/Windows/System32/vcruntime140.dll"
                "C:/Windows/System32/vcruntime140_1.dll"
                "C:/Windows/System32/msvcp140.dll"
                "C:/Windows/System32/msvcp140_1.dll"
                "C:/Windows/System32/msvcp140_2.dll"
                "C:/Windows/System32/concrt140.dll"
                "$<TARGET_FILE_DIR:Source2>/"
            # Copy Wolfram WSTP DLL (required for MAIN_1_CoAnQi.exe physics calculations)
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "C:/Program Files/Wolfram Research/Wolfram Engine/14.3/SystemFiles/Links/WSTP/DeveloperKit/Windows-x86-64/SystemAdditions/wstp64i4.dll"
                "$<TARGET_FILE_DIR:Source2>/"
            # Deploy Qt DLLs (including WebEngine) using windeployqt
            COMMAND "C:/Qt/6.10.0/msvc2022_64/bin/windeployqt.exe" --release --no-translations "$<TARGET_FILE:Source2>"
            COMMAND ${CMAKE_COMMAND} -E echo "✅ MAIN_1_CoAnQi.exe + Python wrappers (UQFF, QCalc, CondensedPhysics, ValidationData, GrokAPI, OllamaAPI) + OpenSSL + VC++ Runtime + WSTP DLLs + Qt6 deployed"
            COMMENT "Deploying MAIN_1_CoAnQi.exe, Python wrappers, OpenSSL, VC++ Runtime, and Wolfram WSTP for cross-platform integration"
        )
    endif()
    
    message(STATUS "Source2 (HEAD PROGRAM) build target ENABLED")
    
    # ========== Source2_HEAD_PROGRAM (VR/VM GPU Backend) ==========
    # HEAD PROGRAM: GPU-heavy VR/VM backend with QWebEngineView, threaded executor
    # Features: DraggableButton, jqMathRenderer, ThreadedSymPyExecutor
    # Separated from Source2.cpp to avoid runtime conflicts with invasive dependencies
    option(BUILD_SOURCE2_HEAD_PROGRAM "Build Source2 HEAD PROGRAM (VR/VM backend)" OFF)
    if(BUILD_SOURCE2_HEAD_PROGRAM)
        # Note: filename has spaces/parentheses - must skip AUTOMOC for this file
        # Q_OBJECT macros handled manually without MOC (no Qt signals/slots used)
        set_source_files_properties("source2(HEAD PROGRAM).cpp" PROPERTIES SKIP_AUTOMOC ON)
        
        # VR/Astro graphics sources
        # NOTE: External VR files (task_bot, poseidon_task_bot, etc.) have incomplete implementations
        # VR runtime is embedded in source2(HEAD PROGRAM).cpp via namespace VR {}
        # Enable these when APIs are complete:
        # set(VR_SOURCES
        #     ${CMAKE_SOURCE_DIR}/vr/astro_graphics.cpp
        #     ${CMAKE_SOURCE_DIR}/vr/vulkan_compositor.cpp
        #     ${CMAKE_SOURCE_DIR}/vr/vr_runtime.cpp
        #     ${CMAKE_SOURCE_DIR}/vr/openxr_session.cpp
        #     ${CMAKE_SOURCE_DIR}/vr/task_bot.cpp
        #     ${CMAKE_SOURCE_DIR}/vr/poseidon_task_bot.cpp
        #     ${CMAKE_SOURCE_DIR}/vr/CoAnQi_bot.cpp
        # )
        
        add_executable(Source2_HEAD_PROGRAM 
            "source2(HEAD PROGRAM).cpp"
            ${IPC_SOURCES}
            ${IPC_HEADERS}
            # ${VR_SOURCES}  # Disabled until VR files are complete
            ${CMAKE_SOURCE_DIR}/resources/Star-Magic.rc  # Windows icon/version resources
        )
        
        set_target_properties(Source2_HEAD_PROGRAM PROPERTIES 
            AUTOMOC OFF
            AUTOUIC ON 
            AUTORCC ON
            OUTPUT_NAME "Source2_HEAD_PROGRAM"
            MSVC_RUNTIME_LIBRARY "MultiThreadedDLL$<$<CONFIG:Debug>:Debug>"
        )
        
        target_include_directories(Source2_HEAD_PROGRAM PRIVATE 
            ${CMAKE_CURRENT_SOURCE_DIR}/ipc
            ${CMAKE_CURRENT_SOURCE_DIR}/vr
            ${CMAKE_CURRENT_SOURCE_DIR}/resources
        )
        
        # Link Qt6 (including WebEngineWidgets for jqMath)
        target_link_libraries(Source2_HEAD_PROGRAM PRIVATE 
            Qt6::Core 
            Qt6::Widgets 
            Qt6::Network
        )
        
        # QtWebEngine for jqMath GPU rendering
        if(TARGET Qt6::WebEngineWidgets)
            target_link_libraries(Source2_HEAD_PROGRAM PRIVATE Qt6::WebEngineWidgets)
            target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE HAVE_WEBENGINE)
            message(STATUS "  Source2_HEAD_PROGRAM: Qt6::WebEngineWidgets linked (jqMath GPU rendering)")
        else()
            message(STATUS "  Source2_HEAD_PROGRAM: Qt6::WebEngineWidgets NOT found - jqMath disabled")
        endif()
        
        # pybind11 for ThreadedSymPyExecutor
        if(pybind11_FOUND)
            target_link_libraries(Source2_HEAD_PROGRAM PRIVATE pybind11::embed)
        else()
            target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE NO_PYTHON)
        endif()
        
        # VTK for 3D visualization (optional)
        if(VTK_FOUND)
            target_link_libraries(Source2_HEAD_PROGRAM PRIVATE ${VTK_LIBRARIES})
        else()
            target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE NO_VTK)
        endif()
        
        # CURL for API calls
        if(CURL_FOUND)
            target_link_libraries(Source2_HEAD_PROGRAM PRIVATE CURL::libcurl)
        else()
            target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE NO_CURL)
        endif()
        
        # AWS SDK (DISABLED for HEAD PROGRAM - causes runtime library conflicts)
        # Enable when AWS SDK is rebuilt with matching runtime library
        target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE NO_AWS)
        # AWS linking disabled - causes MSVC runtime library conflicts
        # if(aws-cpp-sdk-core_FOUND)
        #     target_link_libraries(Source2_HEAD_PROGRAM PRIVATE 
        #         aws-cpp-sdk-core
        #         aws-cpp-sdk-s3
        #         aws-cpp-sdk-cognito-idp
        #     )
        # else()
        #     target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE NO_AWS)
        # endif()
        
        # Qt3D (DISABLED - not installed, causes build errors)
        # Phase 0D: Quick fix for IPC testing - disable Qt3D includes
        target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE NO_QT3D)
        
        # LibTorch / PyTorch C++ (DISABLED - not installed, Phase 0 uses QCalc instead)
        # Phase 0D: Disable libtorch for IPC testing
        target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE NO_LIBTORCH)
        
        # SymEngine (DISABLED - not installed)
        # Phase 0D: Disable SymEngine for IPC testing
        target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE NO_SYMENGINE)
        
        # OpenCV (optional - for vision processing)
        if(OpenCV_FOUND)
            target_link_libraries(Source2_HEAD_PROGRAM PRIVATE ${OpenCV_LIBS})
        else()
            target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE NO_OPENCV)
        endif()
        
        # SQLite (optional - for local caching)
        if(SQLite3_FOUND)
            target_link_libraries(Source2_HEAD_PROGRAM PRIVATE SQLite::SQLite3)
        else()
            target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE NO_SQLITE)
        endif()
        
        # libwebsockets (optional - for real-time data feeds)
        # Not configured yet - disable WebSocket features
        target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE NO_WEBSOCKETS)
        
        # Qalculate (optional - for symbolic math)
        if(Qalculate_FOUND)
            target_link_libraries(Source2_HEAD_PROGRAM PRIVATE ${Qalculate_LIBRARIES})
        else()
            target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE NO_QALCULATE)
        endif()
        
        # PocketSphinx (optional - for voice commands)
        if(PocketSphinx_FOUND)
            target_link_libraries(Source2_HEAD_PROGRAM PRIVATE ${PocketSphinx_LIBRARIES})
        else()
            target_compile_definitions(Source2_HEAD_PROGRAM PRIVATE NO_POCKETSPHINX)
        endif()
        
        message(STATUS "Source2_HEAD_PROGRAM (VR/VM GPU backend) build target ENABLED")
        message(STATUS "  Build with: cmake --build . --target Source2_HEAD_PROGRAM")
    endif()

    # ========== VR RUNTIME (Phase 3 - SCAFFOLDED) ==========
    # Phase 3: VR Runtime scaffold with OpenXR + Vulkan
    # Uses dedicated vr/ directory with modular components:
    #   - vr_runtime.cpp: Main VR compositor
    #   - openxr_session.cpp: OpenXR session management
    #   - vulkan_compositor.cpp: Vulkan GPU rendering
    #   - task_bot.cpp: Voice/gesture command processing
    #   - astro_graphics.cpp: Astronomical visualization
    # Note: Phase 2 Python bridge is complete in Source2.exe --service mode
    
    option(BUILD_VR_RUNTIME "Build VRRuntime target (Phase 3)" OFF)
    if(BUILD_VR_RUNTIME)
        # Find OpenXR (optional - stubs used if not found)
        find_package(OpenXR CONFIG QUIET)
        if(OpenXR_FOUND)
            message(STATUS "OpenXR found - VR headset support enabled")
            add_definitions(-DHAVE_OPENXR)
        else()
            message(STATUS "OpenXR NOT found - using stub implementation")
            message(STATUS "  To enable: vcpkg install openxr-loader[vulkan]:x64-windows")
        endif()
        
        # Find Vulkan (optional - stubs used if not found)
        find_package(Vulkan QUIET)
        if(Vulkan_FOUND)
            message(STATUS "Vulkan found at ${Vulkan_VERSION} - GPU compositor enabled")
            add_definitions(-DHAVE_VULKAN)
        else()
            message(STATUS "Vulkan NOT found - using stub implementation")
            message(STATUS "  To enable: Install Vulkan SDK from https://vulkan.lunarg.com/")
        endif()
        
        # VR Runtime source files
        set(VR_RUNTIME_SOURCES
            ${CMAKE_SOURCE_DIR}/vr/vr_runtime.cpp
            ${CMAKE_SOURCE_DIR}/vr/openxr_session.cpp
            ${CMAKE_SOURCE_DIR}/vr/vulkan_compositor.cpp
            ${CMAKE_SOURCE_DIR}/vr/task_bot.cpp
            ${CMAKE_SOURCE_DIR}/vr/poseidon_task_bot.cpp
            ${CMAKE_SOURCE_DIR}/vr/CoAnQi_bot.cpp
            ${CMAKE_SOURCE_DIR}/vr/astro_graphics.cpp
            ${CMAKE_SOURCE_DIR}/ipc/uqff_ipc.cpp
        )
        
        set(VR_RUNTIME_HEADERS
            ${CMAKE_SOURCE_DIR}/vr/vr_runtime.h
            ${CMAKE_SOURCE_DIR}/vr/openxr_session.h
            ${CMAKE_SOURCE_DIR}/vr/vulkan_compositor.h
            ${CMAKE_SOURCE_DIR}/vr/task_bot.h
            ${CMAKE_SOURCE_DIR}/vr/poseidon_task_bot.h
            ${CMAKE_SOURCE_DIR}/vr/CoAnQi_bot.h
            ${CMAKE_SOURCE_DIR}/vr/astro_graphics.h
            ${CMAKE_SOURCE_DIR}/ipc/uqff_ipc.h
        )
        
        add_executable(VRRuntime ${VR_RUNTIME_SOURCES} ${VR_RUNTIME_HEADERS})
        
        target_compile_definitions(VRRuntime PRIVATE VR_RUNTIME_STANDALONE)
        
        set_target_properties(VRRuntime PROPERTIES 
            MSVC_RUNTIME_LIBRARY "MultiThreadedDLL$<$<CONFIG:Debug>:Debug>"
        )
        
        target_include_directories(VRRuntime PRIVATE
            ${CMAKE_SOURCE_DIR}
            ${CMAKE_SOURCE_DIR}/vr
            ${CMAKE_SOURCE_DIR}/ipc
        )
        
        # OpenXR linking
        if(OpenXR_FOUND)
            target_link_libraries(VRRuntime PRIVATE OpenXR::openxr_loader)
        endif()
        
        # Vulkan linking
        if(Vulkan_FOUND)
            target_link_libraries(VRRuntime PRIVATE Vulkan::Vulkan)
        endif()
        
        # VTK for scientific visualization (optional)
        if(VTK_FOUND)
            target_link_libraries(VRRuntime PRIVATE ${VTK_LIBRARIES})
            target_include_directories(VRRuntime PRIVATE ${VTK_INCLUDE_DIRS})
            target_compile_definitions(VRRuntime PRIVATE HAVE_VTK)
        endif()
        
        # CURL for web APIs (optional)
        if(CURL_FOUND)
            target_link_libraries(VRRuntime PRIVATE CURL::libcurl)
            target_compile_definitions(VRRuntime PRIVATE HAVE_CURL)
        endif()
        
        # OpenCV for gesture recognition (optional)
        if(OpenCV_FOUND)
            target_link_libraries(VRRuntime PRIVATE ${OpenCV_LIBS})
            target_compile_definitions(VRRuntime PRIVATE HAVE_OPENCV)
        endif()
        
        # JSON library for catalog parsing
        find_package(nlohmann_json CONFIG QUIET)
        if(nlohmann_json_FOUND)
            target_link_libraries(VRRuntime PRIVATE nlohmann_json::nlohmann_json)
            target_compile_definitions(VRRuntime PRIVATE HAVE_NLOHMANN_JSON)
        endif()
        
        message(STATUS "VRRuntime build target ENABLED (Phase 3 scaffold + v4.2.1 Dual Bot)")
        message(STATUS "  Sources: vr/vr_runtime.cpp, vr/openxr_session.cpp, vr/vulkan_compositor.cpp")
        message(STATUS "           vr/task_bot.cpp, vr/poseidon_task_bot.cpp, vr/CoAnQi_bot.cpp, vr/astro_graphics.cpp")
        message(STATUS "  Bots: CoAnQi_bot (MAIN_1_CoAnQi.cpp specialist), Poseidon (general contractor)")
    else()
        message(STATUS "VRRuntime DISABLED - set BUILD_VR_RUNTIME=ON to enable (Phase 3)")
    endif()
else()
    message(STATUS "Source2 DISABLED - Qt6 required (install qt6-base:x64-windows)")
endif()

# ========== BUILD SUMMARY (November 30, 2025 - Commit 92c9908) ==========
# PhysicsTerm Classes:
#   - MAIN_1_CoAnQi.cpp: 1,172 classes (Batches 1-17, SOURCE1-117)
#   - wolfram_physics_classes.cpp: 5,703 classes (Batch 18, Nov 22-23)
#   - wolfram_extraction/: 187 classes (Batch 19, Nov 23-30)
#   - TOTAL DEFINED: 7,062 classes
#   - TOTAL REGISTERED: 6,784 classes (96.1% registration rate)
#   - UNREGISTERED: 278 classes (helpers, experimental, commented)
#
# Mathematical Computation Framework:
#   - 10 distinct calculation patterns identified
#   - Computational efficiency: 81% O(1) constant-time operations
#   - 15% O(n) linear complexity, 4% complex algorithms
#   - Verification: WOLFRAM_LAST_EVALUATION_REPORT.md (496 lines)
#
# Modularized Sources: 4 modules (source4, source5, source6, source7)
#   - 46 classes (source4: UQFF + MUGE)
#   - 45 classes (source5: UQFF Enhanced)
#   - 29 classes (source6: Graphics + Physics)
#   - 24 classes (source7: MUGE 7-System)
#   - Total Modularized: 144 classes
#
# Main Platforms:
#   - MAIN_1.cpp: Original mathematical framework
#   - MAIN_1_CoAnQi.cpp: Production calculator (12-option menu, SOURCE1-117, 6,784 registered classes)
#   - UQFFCore Library: Consolidated architecture (Core/Modules/*.cpp)
#
# WSTP Integration:
#   - Status: Dormant (USE_EMBEDDED_WOLFRAM=ON available)
#   - Potential: 5,703 live Wolfram queries when activated
#   - Dependencies: wstp64i4.lib, WSTP64I4.dll (MSVC only)
#
# Executable Size:
#   - Uncompressed: ~9.14 MB (MSVC Release-MaxCompress)
#   - UPX compressed: 1.34 MB (85.3% reduction)
#
# Build Configuration:
#   - Compiler: MSVC 19.44.35219.0 (Visual Studio 2022)
#   - Standard: C++20 strict compliance
#   - Optimization: /Os /GL /LTCG /Gw /Gy /Zc:inline /GF /Oi
#   - Remote: origin/master synced (commit 775c7c9 tagged v1.0-wolfram-5890-complete)

# ========== ADDITIONAL SOURCE MODULES ==========

# ========== MODULARIZED WOLFRAM SOURCE MODULES ==========
# Source4, Source5, Source6, Source7 follow 4-file modularization pattern
# Each has: main infrastructure, category_A, category_B, simulation_harness

# Source4: UQFF + MUGE (46 PhysicsTerm classes across 4 files)
# - source4_wolfram.cpp: 24 core classes (helpers, unified field, astrophysical systems)
# - source4_wolfram_compressed.cpp: 9 MUGE compressed components
# - source4_wolfram_resonance.cpp: 13 MUGE resonance components
# - source4_simulation_harness.cpp: Interactive simulation engine
add_executable(source4_wolfram source4_wolfram.cpp source4.cpp)
add_executable(source4_compressed source4_wolfram_compressed.cpp source4.cpp)
add_executable(source4_resonance source4_wolfram_resonance.cpp source4.cpp)
add_executable(source4_simulator 
    source4_simulation_harness.cpp
    source4_wolfram.cpp
    source4_wolfram_compressed.cpp
    source4_wolfram_resonance.cpp
    source4.cpp
)

# Source5: UQFF Enhanced (45 PhysicsTerm classes across 4 files)
# - source5_wolfram.cpp: 21 core classes (Ug1-3, Um, helpers, wrappers)
# - source5_wolfram_compressed.cpp: 10 compressed components + wrapper
# - source5_wolfram_resonance.cpp: 14 resonance components + wrapper
# - source5_simulation_harness.cpp: 6-option interactive menu
add_executable(source5_wolfram source5_wolfram.cpp Source5.cpp)
add_executable(source5_compressed source5_wolfram_compressed.cpp Source5.cpp)
add_executable(source5_resonance source5_wolfram_resonance.cpp Source5.cpp)
add_executable(source5_simulator 
    source5_simulation_harness.cpp
    source5_wolfram.cpp
    source5_wolfram_compressed.cpp
    source5_wolfram_resonance.cpp
    Source5.cpp
)

# Source6: Graphics + Physics (29 PhysicsTerm classes across 4 files)
# - source6_wolfram.cpp: Infrastructure (CelestialBody, PhysicsTermRegistry)
# - source6_wolfram_graphics.cpp: 14 graphics classes (OpenGL, Vulkan, mesh, texture, etc.)
# - source6_wolfram_physics.cpp: 15 UQFF physics classes (helpers + Ug1-4 + Ubi + Um + metrics)
# - source6_register.cpp: Registration implementation
# - source6_simulation_harness.cpp: Interactive menu
add_executable(Source6 Source6.cpp)  # Refactored: 2.0-Enhanced with UQFFModule6
add_executable(source6_wolfram source6_wolfram.cpp)
add_executable(source6_simulator 
    source6_simulation_harness.cpp
    source6_wolfram.cpp
)

# Source7: MUGE 7-System (24 PhysicsTerm classes across 4 files)
# - source7_wolfram.cpp: 2 infrastructure classes (YAMLConfigLoader, ArchiveMediaManager)
# - source7_wolfram_compressed.cpp: 9 compressed MUGE components
# - source7_wolfram_resonance.cpp: 13 resonance MUGE components
# - source7_simulation_harness.cpp: 7-option menu (7 MUGE systems)
add_executable(Source7 source7.cpp)  # Refactored: 2.0-Enhanced with UQFFModule7, FluidSolver, MUGE
add_executable(source7_wolfram source7_wolfram.cpp)
add_executable(source7_compressed source7_wolfram_compressed.cpp)
add_executable(source7_resonance source7_wolfram_resonance.cpp)
add_executable(source7_simulator 
    source7_simulation_harness.cpp
    source7_wolfram.cpp
    source7_wolfram_compressed.cpp
    source7_wolfram_resonance.cpp
)

# Source9: UQFF Unified Field Calculator (CoAnQi Namespace Organization)
add_executable(Source9 source9.cpp)  # Refactored: 2.0-Enhanced with UQFFModule9, CoAnQi, MUGE

# Source10: UQFF Catalogue & Calculator (F_U_Bi_i, compressed_g, 26-layer)
add_executable(Source10 source10.cpp)  # Refactored: 2.0-Enhanced with UQFFModule10, 15 systems

# Source11: UQFF Unified Field Calculator (compute_FU, Ug1-4, Ubi, Um)
add_executable(Source11 source11.cpp)  # Refactored: 2.0-Enhanced with UQFFModule11, full physics

# Source13: Magnetar SGR 1745-2900 MUGE Calculator (12 physics terms)
add_executable(Source13 source13.cpp)  # Refactored: 2.0-Enhanced with UQFFModule13, magnetar physics

# ========== STANDARD SOURCE MODULES ==========
# Single-file executables (not yet modularized)
# add_executable(source2 source2.cpp)  # Requires Qt5 - Reserved for HEAD PROGRAM
add_executable(source4 source4.cpp)  # Original monolithic version (keep for compatibility)
add_executable(Source5 Source5.cpp)  # Original monolithic version (keep for compatibility)
# add_executable(source9 source9.cpp)  # Requires VTK library
add_executable(source10 source10.cpp)
# Source13 already defined above

# ========== ALIGNED UQFF MODULES (Jan 2026) ==========
# All modules below have: UQFFConfigN singleton, UQFFModuleN wrapper, 8 standard methods, uqff_constants.h
add_executable(source14 source14.cpp)  # SGR 0501+4516 Magnetar (12-term MUGE, B-field decay, time-reversal)
add_executable(Source15 source15.cpp)  # Sagittarius A* SMBH (12-term MUGE, M(t) growth, Kerr spin)
add_executable(source16 source16.cpp)  # NGC 2014/2020 Tapestry (11-term MUGE, star formation, stellar wind)
add_executable(source17 source17.cpp)  # Westerlund 2 Cluster (11-term MUGE, 30000 M_sun, τ_SF=2 Myr)
add_executable(source18 source18.cpp)  # Pillars of Creation M16 (12-term MUGE, E(t) erosion photoevaporation)
add_executable(source19 source19.cpp)  # Rings of Relativity (11-term MUGE, L(t) lensing, z=0.5 Einstein ring)
add_executable(source20 source20.cpp)  # Galaxy NGC 2525 (12-term MUGE, g_BH SMBH, M_SN(t) SN 2018gv)
add_executable(source21 source21.cpp)  # NGC 3603 Extreme Cluster (13-term MUGE, P(t) cavity pressure, wind-blown bubble)
add_executable(source22 source22.cpp)  # Bubble Nebula NGC 7635 (11-term MUGE, g_wind stellar wind feedback, O-type star)
add_executable(source23 source23.cpp)  # Antennae Galaxies NGC 4038/4039 (12-term MUGE, I(t) merger interaction, M(t) star formation)
add_executable(source24 source24.cpp)  # Horsehead Nebula Barnard 33 (11-term MUGE, E(t) erosion photoevaporation, sigma Orionis UV radiation)
add_executable(source25 source25.cpp)  # NGC 1275 Perseus A (13-term MUGE, C(t) cooling flow, AGN feedback, M_BH=8e8 M_sun)
add_executable(source26 source26.cpp)  # HUDF Galaxies Galore (12-term MUGE, M(t) star formation, I(t) inter-galaxy interaction, z=3.5)
add_executable(source27 source27.cpp)  # NGC 1792 "The Stellar Forge" Starburst Galaxy (9-term MUGE, M(t) star formation, supernova feedback, z=0.0095)
add_executable(source28 source28.cpp)  # Andromeda Galaxy M31 (11-term MUGE, Ug1-4, Lambda, quantum, EM, fluid, resonant, DM 80%, z=-0.001 blueshift)
add_executable(source29 source29.cpp)  # Sombrero Galaxy M104 (12-term MUGE, SMBH 1e9 M_sun, dust lane, SC correction, z=0.0063 Virgo Cluster)
add_executable(source30 source30.cpp)  # Saturn Evolution (11-term MUGE, g_sun+g_saturn+T_ring+Ug_sum+Lambda+quantum+EM+fluid+resonant+DM+wind, SC correction, z=0)
add_executable(source31 source31.cpp)  # M16 Eagle Nebula Pillars of Creation (8-term MUGE, M_sf star formation, E_rad radiation erosion, SFR=1 M_sun/yr, z=0.0015)
add_executable(source32 source32.cpp)  # Crab Nebula SNR M1 (12-term MUGE, r(t) expansion, pulsar wind P=5e31W, magnetic, superconductivity, z=0.0015)
add_executable(source33 source33.cpp)  # SGR 1745-2900 Magnetar (10-term MUGE, B=2e10 T high-B superconductivity, P=3.76s spin, Galactic Center z=0)
add_executable(source34 source34.cpp)  # SGR 1745-2900 Frequency/Resonance (11-term MUGE, DPM+THz+VacDiff+SuperFreq+AetherRes+Ug4i+QuantumFreq+AetherFreq+FluidFreq+Osc+ExpFreq, all frequency-based)
add_executable(source35 source35.cpp)  # Sagittarius A* SMBH Frequency/Resonance (11-term MUGE, M=4.3e6 M_sun, r_s=1.27e10 m, f_DPM=1e9 Hz scaled for SMBH, accretion dynamics)
add_executable(source36 source36.cpp)
add_executable(source37 source37.cpp)
add_executable(source38 source38.cpp)
add_executable(source39 source39.cpp)
add_executable(source40 source40.cpp)
add_executable(source41 source41.cpp)
add_executable(source42 source42.cpp)
add_executable(source43 source43.cpp)
add_executable(source44 source44.cpp)
add_executable(source45 source45.cpp)
add_executable(source46 source46.cpp)
add_executable(source47 source47.cpp)
add_executable(source48 source48.cpp)
add_executable(Source48-5 Source48-5.cpp)
add_executable(source49 source49.cpp)
add_executable(source50 source50.cpp)
add_executable(source52 source52.cpp)
add_executable(source54 source54.cpp)
add_executable(source56 source56.cpp)
add_executable(source57 source57.cpp)
add_executable(source60 source60.cpp)
add_executable(source64 source64.cpp)
add_executable(source65 source65.cpp)
add_executable(source66 source66.cpp)
add_executable(source67 source67.cpp)
add_executable(source68 source68.cpp)
add_executable(source69 source69.cpp)
add_executable(source70 source70.cpp)
add_executable(source71 source71.cpp)
add_executable(source72 source72.cpp)
add_executable(source73 source73.cpp)
add_executable(source74 source74.cpp)
add_executable(source76 source76.cpp)
add_executable(source77 source77.cpp)
add_executable(source78 source78.cpp)
add_executable(source79 source79.cpp)
add_executable(source80 source80.cpp)
add_executable(source81 source81.cpp)
add_executable(source82 source82.cpp)
add_executable(source83 source83.cpp)
add_executable(source84 source84.cpp)
add_executable(source85 source85.cpp)
add_executable(source86 source86.cpp)
add_executable(source87 source87.cpp)
add_executable(source88 source88.cpp)
add_executable(source89 source89.cpp)
add_executable(source90 source90.cpp)
add_executable(source91 source91.cpp)
add_executable(source92 source92.cpp)
add_executable(source93 source93.cpp)
add_executable(source94 source94.cpp)
add_executable(source95 source95.cpp)
add_executable(source96 source96.cpp)
add_executable(source97 source97.cpp)
add_executable(source98 source98.cpp)
add_executable(source100 source100.cpp)
add_executable(source101 source101.cpp)
add_executable(source102 source102.cpp)
add_executable(Source103 Source103.cpp)
add_executable(source104 source104.cpp)
add_executable(source105 source105.cpp)
add_executable(source106 source106.cpp)
add_executable(source107 source107.cpp)
add_executable(source108 source108.cpp)
add_executable(source109 source109.cpp)
add_executable(source110 source110.cpp)
add_executable(source111 source111.cpp)
add_executable(source112 source112.cpp)
add_executable(source113 source113.cpp)
add_executable(source114 source114.cpp)
add_executable(Source115 Source115.cpp)
add_executable(source116 source116.cpp)
add_executable(source117 source117.cpp)
add_executable(source118 source118.cpp)
add_executable(Source119 Source119.cpp)
add_executable(source120 source120.cpp)
add_executable(source121 source121.cpp)
add_executable(Source122 Source122.cpp)
add_executable(source123 source123.cpp)
add_executable(source124 source124.cpp)
add_executable(source125 source125.cpp)
add_executable(source126 source126.cpp)
add_executable(source127 source127.cpp)
add_executable(source128 source128.cpp)
add_executable(source129 source129.cpp)
add_executable(source130 source130.cpp)
add_executable(source131 source131.cpp)
add_executable(source132 source132.cpp)
add_executable(source133 source133.cpp)
add_executable(Source134 Source134.cpp)
add_executable(source135 source135.cpp)
add_executable(source136 source136.cpp)
add_executable(source137 source137.cpp)
add_executable(source138 source138.cpp)
add_executable(Source139 Source139.cpp)
add_executable(source140 source140.cpp)
add_executable(source141 source141.cpp)
add_executable(source142 source142.cpp)
add_executable(Source143 Source143.cpp)
add_executable(source144 source144.cpp)
add_executable(source145 source145.cpp)
add_executable(source146 source146.cpp)
add_executable(source147 source147.cpp)
add_executable(source148 source148.cpp)
add_executable(source149 source149.cpp)
add_executable(source150 source150.cpp)
add_executable(source151 source151.cpp)
add_executable(source152 source152.cpp)
add_executable(source153 source153.cpp)
add_executable(Source154 Source154.cpp)
add_executable(Source155 Source155.cpp)
add_executable(Source156 Source156.cpp)
add_executable(Source157 Source157.cpp)
add_executable(Source158 Source158.cpp)
add_executable(Source159 Source159.cpp)
add_executable(Source160 Source160.cpp)
add_executable(Source161 Source161.cpp)
add_executable(Source162 Source162.cpp)
add_executable(Source163 Source163.cpp)
add_executable(Source164 Source164.cpp)
add_executable(Source165 Source165.cpp)
add_executable(Source166 Source166.cpp)
add_executable(Source167 Source167.cpp)

# ========== SOURCE200 COSMIC QUANTUM EGG (WSTP Required) ==========
# Added: January 22, 2026
# Purpose: 26D Chaotic Dimensional Structure - Integrates into MAIN_1_CoAnQi.cpp
# Requires: Wolfram WSTP (includes source174_wolfram_bridge_embedded.cpp)
if(USE_EMBEDDED_WOLFRAM AND EXISTS "${WSTP_DIR}")
    add_executable(source200 source200_cosmic_quantum_egg.cpp)
    target_compile_definitions(source200 PRIVATE USE_EMBEDDED_WOLFRAM USE_COSMIC_QUANTUM_EGG)
    target_link_libraries(source200 
        wstp64i4
        kernel32
        user32
        advapi32
    )
    message(STATUS "source200 (Cosmic Quantum Egg) enabled with WSTP integration")
else()
    message(STATUS "source200 (Cosmic Quantum Egg) SKIPPED - requires USE_EMBEDDED_WOLFRAM=ON")
endif()

# ========== UQFF CORE LIBRARY ==========
# Added: November 14, 2025
# Purpose: Consolidated Core architecture from extracted modules
# Updated: Matches successful build_backup_20251116_072753 configuration

# Foundation modules
set(CORE_FOUNDATION
    Core/SystemCatalogue.cpp
    Core/UQFFModule4.cpp
)

# Collect all Core module source files
file(GLOB CORE_MODULES "Core/Modules/*.cpp")

# Remove backup files and files with unmet dependencies
list(FILTER CORE_MODULES EXCLUDE REGEX ".*_CORRUPTED_BACKUP\\.cpp$")
list(FILTER CORE_MODULES EXCLUDE REGEX ".*_HEALED_20251115_0140\\.cpp$")

# Conditionally exclude modules based on available dependencies
if(NOT Qt5_FOUND)
    list(FILTER CORE_MODULES EXCLUDE REGEX ".*/ForceModule\\.cpp$")
    list(FILTER CORE_MODULES EXCLUDE REGEX ".*/ScientificCalculatorDialog\\.cpp$")
    list(FILTER CORE_MODULES EXCLUDE REGEX ".*/InputModule\\.cpp$")
    list(FILTER CORE_MODULES EXCLUDE REGEX ".*/SIMPlugin\\.cpp$")
endif()

if(NOT GMP_FOUND)
    list(FILTER CORE_MODULES EXCLUDE REGEX ".*/SymEngine\\.cpp$")
    list(FILTER CORE_MODULES EXCLUDE REGEX ".*/SymEngineAllocator\\.cpp$")
endif()

# FluidSolver requires VTK - exclude it
list(FILTER CORE_MODULES EXCLUDE REGEX ".*/FluidSolver\\.cpp$")

# Build UQFFCore as static library
add_library(UQFFCore STATIC
    ${CORE_FOUNDATION}
    ${CORE_MODULES}
)

# Install UQFFCore.hpp master header
target_sources(UQFFCore PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES Core/UQFFCore.hpp
)

target_include_directories(UQFFCore PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Core
    ${CMAKE_CURRENT_SOURCE_DIR}/Core/Modules
    C:/vcpkg/installed/x64-windows/include
)

# Configure OpenGL/GLEW, GLFW, GLM, Assimp, yaml-cpp, nlohmann-json from vcpkg installations
set(CMAKE_PREFIX_PATH "C:/vcpkg/installed/x64-windows")
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(glm CONFIG REQUIRED)
find_package(assimp CONFIG REQUIRED)
find_package(yaml-cpp CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

# Optional: Qt5 and GMP
find_package(Qt5 COMPONENTS Core Widgets QUIET)
find_package(GMP QUIET)

# Add GLFW_DLL definition to match successful build
add_definitions(-DGLFW_DLL)

# Link all required libraries to UQFFCore (matches build_backup configuration)
target_link_libraries(UQFFCore PUBLIC 
    OpenGL::GL
    GLEW::GLEW 
    glfw
    glm::glm
    assimp::assimp
    yaml-cpp::yaml-cpp
    nlohmann_json::nlohmann_json
)

# If Qt5 was found, link it
if(Qt5_FOUND)
    message(STATUS "Qt5 found: ${Qt5_VERSION} - enabling GUI modules")
    target_link_libraries(UQFFCore PUBLIC Qt5::Core Qt5::Widgets)
else()
    message(STATUS "Qt5 not found - GUI modules will be excluded")
endif()

# Build summary
list(LENGTH CORE_MODULES CORE_MODULE_COUNT)
message(STATUS "==========================================")
message(STATUS "UQFFCore Library Build Summary")
message(STATUS "==========================================")
message(STATUS "Core Modules: ${CORE_MODULE_COUNT} modules")
message(STATUS "Qt5 Found: ${Qt5_FOUND}")
message(STATUS "GMP Found: ${GMP_FOUND}")
message(STATUS "==========================================")

# add_executable(Source13 source13.cpp)
# add_executable(Source14 source14.cpp)
# add_executable(Source147 Source147.cpp)
# add_executable(Source148 Source148.cpp)

# ========== COMPILER OPTIONS ==========

# Compiler-specific options
if(MSVC)
    target_compile_options(Source134 PRIVATE /Zi)
    target_compile_options(MAIN_1 PRIVATE /Zi /Od)
    target_compile_options(MAIN_1_CoAnQi PRIVATE /Zi /Od)
    target_compile_options(Source163 PRIVATE /Zi)
    target_compile_options(Source164 PRIVATE /Zi)
    target_compile_options(Source165 PRIVATE /Zi)
    target_compile_options(Source166 PRIVATE /Zi)
    target_compile_options(Source167 PRIVATE /Zi)
else()
    target_compile_options(Source134 PRIVATE -g)
    target_compile_options(MAIN_1 PRIVATE -g -O0)
    target_compile_options(MAIN_1_CoAnQi PRIVATE -g -O0)
    target_compile_options(Source163 PRIVATE -g -O0)
    target_compile_options(Source164 PRIVATE -g -O0)
    target_compile_options(Source165 PRIVATE -g -O0)
    target_compile_options(Source166 PRIVATE -g -O0)
    target_compile_options(Source167 PRIVATE -g -O0)
endif()

# ========== OPTIONAL: Add include directories ==========
# include_directories(${CMAKE_SOURCE_DIR}/include)

# ========== OPTIONAL: Link libraries ==========
# If you need to link math library or others
# target_link_libraries(MAIN_1 m)
# target_link_libraries(Source134 m)

# ========== INSTALL TARGETS ==========
install(TARGETS Source134 MAIN_1 MAIN_1_CoAnQi Source163 Source164 Source165 Source166 Source167
    RUNTIME DESTINATION bin
)

# When Source2 is fixed, add it back:
# install(TARGETS MAIN_1 Source134 Source2
#     RUNTIME DESTINATION bin
# )

# ============================================================================
# COANQI ENHANCEMENTS - Plugin System, WebAssembly, ML, Testing
# ============================================================================

# ========== PLUGIN COMPILATION ==========
# Build shared library plugins for CoAnQi extensibility
if(EXISTS "${CMAKE_SOURCE_DIR}/src/plugins/ramanujan_validator.cpp")
    add_library(ramanujan_validator SHARED src/plugins/ramanujan_validator.cpp)
    
    # Link against UQFFCore if it's the core library
    if(TARGET UQFFCore)
        target_link_libraries(ramanujan_validator PUBLIC UQFFCore)
    endif()
    
    # If we create a CoAnQiCore library in the future, link it here
    # target_link_libraries(ramanujan_validator PUBLIC CoAnQiCore)
    
    message(STATUS "Plugin: ramanujan_validator.dll will be built")
    
    install(TARGETS ramanujan_validator
        LIBRARY DESTINATION lib/plugins
        RUNTIME DESTINATION bin/plugins
    )
else()
    message(STATUS "Plugin source not found: src/plugins/ramanujan_validator.cpp - skipping")
endif()

# ========== WEBASSEMBLY BUILD ==========
# Build CoAnQi for web deployment using Emscripten
if(EMSCRIPTEN)
    message(STATUS "Emscripten detected - building WebAssembly version")
    
    set(CMAKE_EXECUTABLE_SUFFIX ".html")
    
    if(EXISTS "${CMAKE_SOURCE_DIR}/src/webassembly/main.cpp")
        add_executable(coanqi_wasm src/webassembly/main.cpp)
        
        # Link WebAssembly-specific options
        target_link_options(coanqi_wasm PUBLIC 
            -s WASM=1 
            -s ALLOW_MEMORY_GROWTH=1
            -s EXPORTED_FUNCTIONS=['_main','_CoAnQi_calculate','_CoAnQi_sync']
            -s EXPORTED_RUNTIME_METHODS=['ccall','cwrap']
            --bind
        )
        
        # Enable IndexedDB filesystem persistence
        target_link_options(coanqi_wasm PUBLIC -lidbfs.js)
        
        message(STATUS "WebAssembly: coanqi_wasm.html will be built")
    else()
        message(STATUS "WebAssembly source not found: src/webassembly/main.cpp - skipping")
    endif()
endif()

# ========== MACHINE LEARNING SUPPORT ==========
# Integrate PyTorch for AI-enhanced validation
find_package(Torch QUIET)
if(Torch_FOUND)
    message(STATUS "PyTorch found at ${Torch_DIR} - enabling ML features")
    
    # Enable ML compilation flag
    add_definitions(-DCOANQI_ML_ENABLED)
    
    # Link PyTorch to core library
    if(TARGET UQFFCore)
        target_compile_definitions(UQFFCore PUBLIC COANQI_ML_ENABLED)
        target_link_libraries(UQFFCore PUBLIC ${TORCH_LIBRARIES})
    endif()
    
    # If source2 becomes a target, link it too
    if(TARGET source2)
        target_compile_definitions(source2 PUBLIC COANQI_ML_ENABLED)
        target_link_libraries(source2 PUBLIC ${TORCH_LIBRARIES})
    endif()
    
    message(STATUS "ML: PyTorch ${Torch_VERSION} linked successfully")
else()
    message(STATUS "PyTorch NOT found - ML features disabled (compile without COANQI_ML_ENABLED)")
    message(STATUS "To enable: Install PyTorch C++ (libtorch) from pytorch.org")
endif()

# ========== UNIT TESTING ==========
# Enable CTest framework for automated testing
enable_testing()

# Add tests subdirectory if it exists
if(EXISTS "${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt")
    add_subdirectory(tests)
    message(STATUS "Testing: Unit tests enabled in tests/ directory")
else()
    message(STATUS "Testing: tests/CMakeLists.txt not found - create it to enable unit tests")
    message(STATUS "Testing: Example structure:")
    message(STATUS "  tests/")
    message(STATUS "    CMakeLists.txt")
    message(STATUS "    test_mathematical_validator.cpp")
    message(STATUS "    test_plugins.cpp")
    message(STATUS "    test_ml_integration.cpp")
endif()

# ========== PYBIND11 PYTHON BINDINGS ==========
option(BUILD_PYTHON_BINDINGS "Build Python bindings via PyBind11" OFF)

if(BUILD_PYTHON_BINDINGS)
    message(STATUS "PyBind11: Looking for Python and pybind11...")
    
    find_package(Python3 COMPONENTS Interpreter Development QUIET)
    find_package(pybind11 CONFIG QUIET)
    
    if(Python3_FOUND AND pybind11_FOUND)
        message(STATUS "PyBind11: Python ${Python3_VERSION} found at ${Python3_EXECUTABLE}")
        message(STATUS "PyBind11: pybind11 found at ${pybind11_DIR}")
        
        pybind11_add_module(uqff_core uqff_pybind.cpp)
        target_compile_features(uqff_core PRIVATE cxx_std_17)
        target_compile_definitions(uqff_core PRIVATE VERSION_INFO="3.0.0")
        
        # Install to Python site-packages
        install(TARGETS uqff_core
            LIBRARY DESTINATION ${Python3_SITELIB}
            RUNTIME DESTINATION ${Python3_SITELIB}
        )
        
        message(STATUS "PyBind11: uqff_core module target created")
    elseif(Python3_FOUND AND NOT pybind11_FOUND)
        message(STATUS "PyBind11: pybind11 not found - install via 'pip install pybind11' and rerun cmake")
        message(STATUS "PyBind11: Or use 'pip install .' from the project root")
    else()
        message(STATUS "PyBind11: Python3 not found - skipping Python bindings")
    endif()
else()
    message(STATUS "PyBind11: Disabled (use -DBUILD_PYTHON_BINDINGS=ON to enable)")
endif()

# ========== ENHANCEMENT SUMMARY ==========
message(STATUS "==========================================")
message(STATUS "CoAnQi Enhancement Features")
message(STATUS "==========================================")
message(STATUS "Plugins:      ramanujan_validator (if available)")
message(STATUS "WebAssembly:  ${EMSCRIPTEN}")
message(STATUS "ML (PyTorch): ${Torch_FOUND}")
message(STATUS "Testing:      ${BUILD_TESTING}")
message(STATUS "==========================================")
message(STATUS "To build all enhancements:")
message(STATUS "  cmake --build build --config Release")
message(STATUS "  ctest --test-dir build --config Release")
message(STATUS "==========================================")
message(STATUS "For WebAssembly build:")
message(STATUS "  emcmake cmake -S . -B build_wasm")
message(STATUS "  cmake --build build_wasm")
message(STATUS "  python -m http.server 8000 --directory build_wasm")
message(STATUS "  Open: http://localhost:8000/coanqi_wasm.html")
message(STATUS "==========================================")