# Star-Magic UQFF Project
# Unified Quantum Field Force Calculations
# Author: Daniel T. Murphy
cmake_minimum_required(VERSION 3.10)

# Project name and programming language
project(StarMagic_UQFF LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable debugging symbols
set(CMAKE_BUILD_TYPE Debug)

# Compiler-specific flags
if(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od /W4")
else()
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
endif()

# Export compile commands for IDE integration
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ========== MAIN EXECUTABLES ==========

# Source134: Abell 2256 Galaxy Cluster (fully functional)
add_executable(Source134 Source134.cpp)

# MAIN_1: Core mathematical framework (ENABLED - repaired Nov 9, 2025)
add_executable(MAIN_1 MAIN_1.cpp)

# MAIN_1_CoAnQi: Conscious Quantum Intelligence UQFF Calculator (273 physics terms)
add_executable(MAIN_1_CoAnQi MAIN_1_CoAnQi.cpp)

# Source2: HEAD PROGRAM - Primary User Base Platform (REQUIRES Qt5, VTK, AWS SDK, OpenCV, etc.)
# Uncomment after installing dependencies via vcpkg:
#   vcpkg install qt5-base:x64-mingw-dynamic qt5-webengine:x64-mingw-dynamic
#   vcpkg install vtk:x64-mingw-dynamic opencv:x64-mingw-dynamic
#   vcpkg install curl:x64-mingw-dynamic sqlite3:x64-mingw-dynamic
#   vcpkg install aws-sdk-cpp[cognito-idp,s3]:x64-mingw-dynamic
# add_executable(Source2 source2.cpp)
# target_link_libraries(Source2 PRIVATE UQFFCore Qt5::Widgets Qt5::WebEngineWidgets ${VTK_LIBRARIES} ${OpenCV_LIBS} CURL::libcurl SQLite::SQLite3 aws-cpp-sdk-cognito-idp aws-cpp-sdk-s3)

# ========== ADDITIONAL SOURCE MODULES ==========

# Source163-167: Latest modules (Nov 2025)
add_executable(Source163 Source163.cpp)
add_executable(Source164 Source164.cpp)
add_executable(Source165 Source165.cpp)
add_executable(Source166 Source166.cpp)
add_executable(Source167 Source167.cpp)

# Source modules 1-162 (add as needed)
# Example: Uncomment to build specific modules

# ========== UQFF CORE LIBRARY ==========
# Added: November 14, 2025
# Purpose: Consolidated Core architecture from extracted modules

# Foundation modules (Week 1)
set(CORE_FOUNDATION
    Core/SystemCatalogue.cpp
)

# Adaptive framework (Week 2)
set(CORE_ADAPTIVE
    Core/UQFFModule4.cpp
    Core/FluidSolver.cpp
)

# Extracted modules (Weeks 3-4): All extracted modules
file(GLOB CORE_MODULES "Core/Modules/*.cpp")
# Remove false keywords and modules with external GUI dependencies (OpenGL/Qt) or needing class updates
list(FILTER CORE_MODULES EXCLUDE REGEX ".*/(FluidSolver|SIMPlugin|ScientificCalculatorDialog|SymEngine|SymEngineAllocator|ForceModule|InputModule|HydrogenResonanceUQFFModule)\\.cpp$")

# Build UQFFCore as static library
add_library(UQFFCore STATIC
    ${CORE_FOUNDATION}
    ${CORE_ADAPTIVE}
    ${CORE_MODULES}
)

target_include_directories(UQFFCore PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Core
)

# Optional: Build integration test
option(BUILD_CORE_TESTS "Build UQFFCore integration tests" ON)
if(BUILD_CORE_TESTS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test_core_integration.cpp")
    add_executable(test_core_integration test_core_integration.cpp)
    target_link_libraries(test_core_integration PRIVATE UQFFCore)
    enable_testing()
    add_test(NAME CoreIntegrationTest COMMAND test_core_integration)
endif()

# Build summary
list(LENGTH CORE_MODULES CORE_MODULE_COUNT)
message(STATUS "UQFFCore Library: ${CORE_MODULE_COUNT} extracted modules + foundation + adaptive")
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