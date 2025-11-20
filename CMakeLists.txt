# Star-Magic UQFF Project
# Unified Quantum Field Force Calculations
# Author: Daniel T. Murphy
cmake_minimum_required(VERSION 3.10)

# Project name and programming language
project(StarMagic_UQFF LANGUAGES CXX)

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

# Set C++ standard to C++20 for maximum VR/VM functionality
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler-specific flags for MSVC
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /std:c++20 /permissive- /Zc:__cplusplus /GR-")
    set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od")
    set(CMAKE_CXX_FLAGS_RELEASE "/Os /GL /Gw /Zc:inline /GF")  # /Os = Favor small code, /GL = Whole Program Optimization, /Gw = Global data optimization, /Zc:inline = Remove unreferenced inline, /GF = String pooling
    
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

# Export compile commands for IDE integration
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ========== WOLFRAM ENGINE INTEGRATION (OPTIONAL) ==========
option(USE_EMBEDDED_WOLFRAM "Enable embedded Wolfram kernel via WSTP" OFF)

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

# MAIN_1_CoAnQi with optional Wolfram integration
if(USE_EMBEDDED_WOLFRAM)
    # source174 is included via source176, don't add it separately to avoid duplicate definitions
    add_executable(MAIN_1_CoAnQi MAIN_1_CoAnQi.cpp)
    
    # Linker → Input → Additional Dependencies: wstp64i4.lib + Windows system libraries
    target_link_libraries(MAIN_1_CoAnQi 
        wstp64i4
        kernel32
        user32
        advapi32
        shell32
    )
    
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
    add_executable(MAIN_1_CoAnQi MAIN_1_CoAnQi.cpp)
    message(STATUS "MAIN_1_CoAnQi compiled WITHOUT Wolfram support")
endif()

# STEP 8 of 8: UPX Ultimate Packer (Post-Build Compression)
# Compresses final executable using UPX with maximum compression settings
# Effect: Dramatic size reduction (50-70% compression ratio typical)
# Target: Get MAIN_1_CoAnQi.exe under 3 MB (from ~5-8 MB uncompressed)
# Requirements: UPX installed at C:/Tools/UPX/upx.exe (download from https://upx.github.io/)
# Trade-offs: Slower startup (decompression overhead ~50-200ms), incompatible with some antivirus
# VS2022 Property Pages: Build Events → Post-Build Event → Command Line
find_program(UPX_EXECUTABLE upx PATHS "C:/Tools/UPX" "C:/Tools/upx" "C:/Program Files/upx" DOC "UPX executable packer")
if(UPX_EXECUTABLE)
    add_custom_command(TARGET MAIN_1_CoAnQi POST_BUILD
        COMMAND ${UPX_EXECUTABLE} --best --ultra-brute "$<TARGET_FILE:MAIN_1_CoAnQi>"
        COMMENT "STEP 8: Compressing MAIN_1_CoAnQi.exe with UPX (--best --ultra-brute)"
        VERBATIM
    )
    message(STATUS "UPX compression ENABLED for MAIN_1_CoAnQi (found at ${UPX_EXECUTABLE})")
else()
    message(WARNING "UPX not found - executable will not be compressed. Install UPX to C:/Tools/upx/ for Step 8 compression.")
    message(WARNING "Download from: https://github.com/upx/upx/releases")
endif()

# Source2: HEAD PROGRAM - Primary User Base Platform (REQUIRES Qt5, VTK, AWS SDK, OpenCV, etc.)
# Uncomment after installing dependencies via vcpkg:
#   vcpkg install qt5-base:x64-mingw-dynamic qt5-webengine:x64-mingw-dynamic
#   vcpkg install vtk:x64-mingw-dynamic opencv:x64-mingw-dynamic
#   vcpkg install curl:x64-mingw-dynamic sqlite3:x64-mingw-dynamic
#   vcpkg install aws-sdk-cpp[cognito-idp,s3]:x64-mingw-dynamic
# add_executable(Source2 source2.cpp)
# target_link_libraries(Source2 PRIVATE UQFFCore Qt5::Widgets Qt5::WebEngineWidgets ${VTK_LIBRARIES} ${OpenCV_LIBS} CURL::libcurl SQLite::SQLite3 aws-cpp-sdk-cognito-idp aws-cpp-sdk-s3)

# ========== ADDITIONAL SOURCE MODULES ==========

# ========== ALL SOURCE MODULES (EXCLUDING 1, 3, 8, 11, 12) ==========
# Complete build configuration - all modules ready
# add_executable(source2 source2.cpp)  # Requires Qt5
add_executable(source4 source4.cpp)
add_executable(Source5 Source5.cpp)
# add_executable(Source6 Source6.cpp)  # Missing 3DGraphics.h
# add_executable(source7 source7.cpp)  # Missing Animation.h
# add_executable(source9 source9.cpp)  # Requires VTK library
add_executable(source10 source10.cpp)
# add_executable(source13 source13.cpp)  # Missing main()
add_executable(source14 source14.cpp)
# add_executable(source15 source15.cpp)  # Missing main()
add_executable(source16 source16.cpp)
add_executable(source17 source17.cpp)
add_executable(source18 source18.cpp)
add_executable(source19 source19.cpp)
add_executable(source20 source20.cpp)
add_executable(source21 source21.cpp)
add_executable(source22 source22.cpp)
add_executable(source23 source23.cpp)
add_executable(source24 source24.cpp)
add_executable(source25 source25.cpp)
add_executable(source26 source26.cpp)
add_executable(source27 source27.cpp)
add_executable(source28 source28.cpp)
add_executable(source29 source29.cpp)
add_executable(source30 source30.cpp)
add_executable(source31 source31.cpp)
add_executable(source32 source32.cpp)
add_executable(source33 source33.cpp)
add_executable(source34 source34.cpp)
add_executable(source35 source35.cpp)
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