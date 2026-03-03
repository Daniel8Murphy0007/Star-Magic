cmake_minimum_required(VERSION 3.20)
project(CalculatorAdvanced VERSION 1.0.0 LANGUAGES CXX)

# C++20 Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Qt6 Settings
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Find packages via vcpkg
find_package(Qt6 COMPONENTS Core Widgets WebEngineWidgets REQUIRED)
find_package(antlr4-runtime CONFIG REQUIRED)
find_package(symengine CONFIG REQUIRED)
find_package(GSL REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(QCustomPlot CONFIG REQUIRED)

# Source files
set(CALCULATOR_SOURCES
    src/antlr4_parser.cpp
    src/symengine_wrapper.cpp
    src/equation_solver.cpp
    src/dimensional_analysis.cpp
    src/polynomial_solver.cpp
    src/uqff_equations.cpp
    src/plotter_widget.cpp
    src/calculator_widget.cpp
)

set(CALCULATOR_HEADERS
    include/antlr4_parser.h
    include/symengine_wrapper.h
    include/equation_solver.h
    include/dimensional_analysis.h
    include/polynomial_solver.h
    include/uqff_equations.h
    include/plotter_widget.h
    include/calculator_widget.h
)

# ANTLR4 Grammar generation (to be added after grammar file extraction)
# antlr4_generate(Equation.g4 EquationLexer.g4 EquationParser.g4)

# Static library for calculator
add_library(calculator_advanced STATIC
    ${CALCULATOR_SOURCES}
    ${CALCULATOR_HEADERS}
)

target_include_directories(calculator_advanced PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}  # For generated ANTLR4 files
)

target_link_libraries(calculator_advanced PUBLIC
    Qt6::Core
    Qt6::Widgets
    Qt6::WebEngineWidgets
    antlr4_static
    symengine
    GSL::gsl
    GSL::gslcblas
    Eigen3::Eigen
    qcustomplot
)

# Optional: Test executable
option(BUILD_CALCULATOR_TESTS "Build calculator tests" ON)
if(BUILD_CALCULATOR_TESTS)
    add_executable(test_calculator
        tests/test_parser.cpp
        tests/test_symbolic.cpp
        tests/test_polynomial.cpp
        tests/test_dimensional.cpp
        tests/test_uqff_integration.cpp
    )
    
    target_link_libraries(test_calculator PRIVATE
        calculator_advanced
    )
    
    enable_testing()
    add_test(NAME CalculatorTests COMMAND test_calculator)
endif()

# Installation
install(TARGETS calculator_advanced
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/
    DESTINATION include/calculator_advanced
    FILES_MATCHING PATTERN "*.h"
)

# Print configuration summary
message(STATUS "==========================================================")
message(STATUS "Calculator Advanced Configuration")
message(STATUS "==========================================================")
message(STATUS "C++ Standard:          ${CMAKE_CXX_STANDARD}")
message(STATUS "Build Type:            ${CMAKE_BUILD_TYPE}")
message(STATUS "Qt6 Version:           ${Qt6_VERSION}")
message(STATUS "ANTLR4:                Found")
message(STATUS "SymEngine:             Found")
message(STATUS "GSL:                   Found")
message(STATUS "Eigen3:                Found")
message(STATUS "QCustomPlot:           Found")
message(STATUS "Build Tests:           ${BUILD_CALCULATOR_TESTS}")
message(STATUS "==========================================================")
