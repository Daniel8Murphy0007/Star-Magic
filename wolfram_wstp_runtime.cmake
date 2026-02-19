# wolfram_wstp_runtime.cmake
# CMake configuration for Wolfram WSTP Runtime Module
#
# USAGE:
# include(wolfram_wstp_runtime.cmake)
# target_link_libraries(your_target WolframWSTPRuntime)
#
# Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
# Date: February 19, 2026

cmake_minimum_required(VERSION 3.20)

# ═══════════════════════════════════════════════════════════════════════════════
# WOLFRAM ENGINE DETECTION
# ═══════════════════════════════════════════════════════════════════════════════

set(WOLFRAM_VERSIONS "14.3;14.2;14.1;14.0;13.3;13.2")
set(WOLFRAM_ROOT "" CACHE PATH "Path to Wolfram Engine installation")
set(WOLFRAM_WSTP_FOUND FALSE)

# Auto-detect Wolfram Engine on Windows
if(WIN32 AND NOT WOLFRAM_ROOT)
    foreach(VER ${WOLFRAM_VERSIONS})
        set(CANDIDATE "C:/Program Files/Wolfram Research/Wolfram Engine/${VER}")
        if(EXISTS "${CANDIDATE}/wolfram.exe")
            set(WOLFRAM_ROOT "${CANDIDATE}")
            message(STATUS "Found Wolfram Engine ${VER} at ${WOLFRAM_ROOT}")
            break()
        endif()
    endforeach()
endif()

# Set WSTP paths
if(WOLFRAM_ROOT AND EXISTS "${WOLFRAM_ROOT}/wolfram.exe")
    set(WSTP_DIR "${WOLFRAM_ROOT}/SystemFiles/Links/WSTP/DeveloperKit/Windows-x86-64/CompilerAdditions")
    
    if(EXISTS "${WSTP_DIR}/wstp.h")
        set(WOLFRAM_WSTP_FOUND TRUE)
        set(WSTP_INCLUDE_DIR "${WSTP_DIR}")
        set(WSTP_LIBRARY "${WSTP_DIR}/wstp64i4.lib")
        set(WSTP_RUNTIME_DLL "${WSTP_DIR}/WSTP64I4.dll")
        message(STATUS "WSTP SDK found at ${WSTP_DIR}")
    else()
        message(WARNING "WSTP SDK not found at ${WSTP_DIR}")
    endif()
else()
    message(STATUS "Wolfram Engine not found - WSTP runtime will use stub implementation")
endif()

# ═══════════════════════════════════════════════════════════════════════════════
# WOLFRAM WSTP RUNTIME MODULE
# ═══════════════════════════════════════════════════════════════════════════════

# Define the module as an OBJECT library for easy integration
add_library(WolframWSTPRuntime STATIC
    ${CMAKE_CURRENT_LIST_DIR}/wolfram_wstp_runtime.cpp
    ${CMAKE_CURRENT_LIST_DIR}/wolfram_wstp_runtime.h
)

target_include_directories(WolframWSTPRuntime PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)

# Enable C++20 features
target_compile_features(WolframWSTPRuntime PUBLIC cxx_std_20)

# Platform-specific settings
if(MSVC)
    target_compile_options(WolframWSTPRuntime PRIVATE /W4 /permissive-)
endif()

# ═══════════════════════════════════════════════════════════════════════════════
# WSTP INTEGRATION (conditional)
# ═══════════════════════════════════════════════════════════════════════════════

if(WOLFRAM_WSTP_FOUND)
    # Add WSTP include directory
    target_include_directories(WolframWSTPRuntime PUBLIC ${WSTP_INCLUDE_DIR})
    
    # Link WSTP library
    target_link_libraries(WolframWSTPRuntime PUBLIC ${WSTP_LIBRARY})
    
    # Enable full WSTP functionality
    target_compile_definitions(WolframWSTPRuntime PUBLIC USE_EMBEDDED_WOLFRAM)
    
    # Copy WSTP DLL to output directory (for runtime)
    if(EXISTS ${WSTP_RUNTIME_DLL})
        add_custom_command(TARGET WolframWSTPRuntime POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${WSTP_RUNTIME_DLL}"
            "$<TARGET_FILE_DIR:WolframWSTPRuntime>"
            COMMENT "Copying WSTP runtime DLL"
        )
    endif()
    
    message(STATUS "WolframWSTPRuntime: Full WSTP integration enabled")
else()
    message(STATUS "WolframWSTPRuntime: Using stub implementation (no Wolfram Engine)")
endif()

# ═══════════════════════════════════════════════════════════════════════════════
# GENERATED TERM CLASSES (187 physics terms)
# ═══════════════════════════════════════════════════════════════════════════════

set(WOLFRAM_GENERATED_CLASSES_DIR "${CMAKE_CURRENT_LIST_DIR}/wolfram_extraction/generated_classes")

if(EXISTS "${WOLFRAM_GENERATED_CLASSES_DIR}/wolfram_master_registration.h")
    # Add generated classes as additional sources
    file(GLOB WOLFRAM_GENERATED_SOURCES "${WOLFRAM_GENERATED_CLASSES_DIR}/*.cpp")
    
    add_library(WolframPhysicsTerms STATIC ${WOLFRAM_GENERATED_SOURCES})
    
    target_include_directories(WolframPhysicsTerms PUBLIC
        ${WOLFRAM_GENERATED_CLASSES_DIR}
        ${CMAKE_CURRENT_LIST_DIR}
    )
    
    target_compile_features(WolframPhysicsTerms PUBLIC cxx_std_20)
    
    # Link to main runtime
    target_link_libraries(WolframWSTPRuntime PUBLIC WolframPhysicsTerms)
    
    message(STATUS "WolframPhysicsTerms: 187 extracted physics term classes included")
else()
    message(STATUS "WolframPhysicsTerms: Generated classes not found (optional)")
endif()

# ═══════════════════════════════════════════════════════════════════════════════
# INTERFACE TARGET FOR HEADER-ONLY USAGE
# ═══════════════════════════════════════════════════════════════════════════════

add_library(WolframWSTPRuntime_HeaderOnly INTERFACE)
target_include_directories(WolframWSTPRuntime_HeaderOnly INTERFACE
    ${CMAKE_CURRENT_LIST_DIR}
)

if(WOLFRAM_WSTP_FOUND)
    target_include_directories(WolframWSTPRuntime_HeaderOnly INTERFACE ${WSTP_INCLUDE_DIR})
    target_compile_definitions(WolframWSTPRuntime_HeaderOnly INTERFACE USE_EMBEDDED_WOLFRAM)
endif()

# ═══════════════════════════════════════════════════════════════════════════════
# STANDALONE TEST TARGET
# ═══════════════════════════════════════════════════════════════════════════════

option(BUILD_WOLFRAM_WSTP_TEST "Build standalone test for Wolfram WSTP Runtime" OFF)

if(BUILD_WOLFRAM_WSTP_TEST)
    add_executable(wolfram_wstp_test
        ${CMAKE_CURRENT_LIST_DIR}/wolfram_wstp_runtime.cpp
    )
    
    target_compile_definitions(wolfram_wstp_test PRIVATE WOLFRAM_WSTP_RUNTIME_STANDALONE_TEST)
    target_include_directories(wolfram_wstp_test PRIVATE ${CMAKE_CURRENT_LIST_DIR})
    target_compile_features(wolfram_wstp_test PUBLIC cxx_std_20)
    
    if(WOLFRAM_WSTP_FOUND)
        target_include_directories(wolfram_wstp_test PRIVATE ${WSTP_INCLUDE_DIR})
        target_link_libraries(wolfram_wstp_test ${WSTP_LIBRARY})
        target_compile_definitions(wolfram_wstp_test PRIVATE USE_EMBEDDED_WOLFRAM)
        
        # Copy DLL
        if(EXISTS ${WSTP_RUNTIME_DLL})
            add_custom_command(TARGET wolfram_wstp_test POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${WSTP_RUNTIME_DLL}"
                "$<TARGET_FILE_DIR:wolfram_wstp_test>"
            )
        endif()
    endif()
    
    message(STATUS "WolframWSTPRuntime: Test executable will be built")
endif()

# ═══════════════════════════════════════════════════════════════════════════════
# EXPORT VARIABLES
# ═══════════════════════════════════════════════════════════════════════════════

set(WOLFRAM_WSTP_RUNTIME_FOUND TRUE CACHE INTERNAL "")
set(WOLFRAM_WSTP_RUNTIME_LIBRARIES WolframWSTPRuntime CACHE INTERNAL "")
set(WOLFRAM_WSTP_RUNTIME_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR} CACHE INTERNAL "")
set(WOLFRAM_WSTP_TERMS_COUNT 187 CACHE INTERNAL "Number of extracted Wolfram physics terms")
