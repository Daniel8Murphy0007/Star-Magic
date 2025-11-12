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

# MAIN_1_CoAnQi: Conscious Quantum Intelligence UQFF Calculator (50+ physics terms)
add_executable(MAIN_1_CoAnQi MAIN_1_CoAnQi.cpp)

# Source2: HEAD PROGRAM (DISABLED - requires Qt framework)
# add_executable(Source2 "source2(HEAD PROGRAM).cpp")

# ========== ADDITIONAL SOURCE MODULES ==========

# Source modules 1-162 (add as needed)
# Example: Uncomment to build specific modules
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
else()
    target_compile_options(Source134 PRIVATE -g)
    target_compile_options(MAIN_1 PRIVATE -g -O0)
    target_compile_options(MAIN_1_CoAnQi PRIVATE -g -O0)
endif()

# ========== OPTIONAL: Add include directories ==========
# include_directories(${CMAKE_SOURCE_DIR}/include)

# ========== OPTIONAL: Link libraries ==========
# If you need to link math library or others
# target_link_libraries(MAIN_1 m)
# target_link_libraries(Source134 m)

# ========== INSTALL TARGETS ==========
install(TARGETS Source134 MAIN_1 MAIN_1_CoAnQi
    RUNTIME DESTINATION bin
)

# When Source2 is fixed, add it back:
# install(TARGETS MAIN_1 Source134 Source2
#     RUNTIME DESTINATION bin
# )