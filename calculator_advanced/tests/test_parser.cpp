// Test suite for ANTLR4 Parser
// Tests expression parsing, derivative/integral detection, error handling
// Part of calculator_advanced framework extracted from Grok thread Iterations #30-32

#include "../include/antlr4_parser.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <string>

// Test helper macro
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        std::cerr << "TEST FAILED: " << message << std::endl; \
        return false; \
    }

// Test Case 1: Basic Expression Parsing
bool test_basic_expressions() {
    std::cout << "Running test_basic_expressions..." << std::endl;
    ANTLR4Parser parser;
    
    // Test simple arithmetic
    auto result1 = parser.parse("2 + 3 * 4");
    TEST_ASSERT(result1.isValid, "Failed to parse '2 + 3 * 4'");
    TEST_ASSERT(result1.type == EquationType::FUNCTIONAL, "Wrong type for arithmetic expression");
    
    // Test with variables
    auto result2 = parser.parse("x^2 + 2*x + 1");
    TEST_ASSERT(result2.isValid, "Failed to parse 'x^2 + 2*x + 1'");
    TEST_ASSERT(result2.variables.size() >= 1, "Failed to detect variable 'x'");
    
    // Test with functions
    auto result3 = parser.parse("sin(x) + cos(x)");
    TEST_ASSERT(result3.isValid, "Failed to parse 'sin(x) + cos(x)'");
    TEST_ASSERT(result3.functions.size() >= 2, "Failed to detect sin/cos functions");
    
    std::cout << "test_basic_expressions PASSED" << std::endl;
    return true;
}

// Test Case 2: Derivative Parsing
bool test_derivative_parsing() {
    std::cout << "Running test_derivative_parsing..." << std::endl;
    ANTLR4Parser parser;
    
    // Test basic derivative d/dx
    auto result1 = parser.parse("d/dx(x^2)");
    TEST_ASSERT(result1.isValid, "Failed to parse 'd/dx(x^2)'");
    TEST_ASSERT(result1.type == EquationType::DERIVATIVE, "Wrong type for derivative");
    
    // Test partial derivative ∂/∂
    auto result2 = parser.parse("∂/∂x(x*y)");
    TEST_ASSERT(result2.isValid, "Failed to parse partial derivative");
    TEST_ASSERT(result2.type == EquationType::DERIVATIVE, "Wrong type for partial derivative");
    
    // Test higher order derivative
    auto result3 = parser.parse("d/dx(d/dx(sin(x)))");
    TEST_ASSERT(result3.isValid, "Failed to parse second derivative");
    
    std::cout << "test_derivative_parsing PASSED" << std::endl;
    return true;
}

// Test Case 3: Integral Parsing
bool test_integral_parsing() {
    std::cout << "Running test_integral_parsing..." << std::endl;
    ANTLR4Parser parser;
    
    // Test indefinite integral
    auto result1 = parser.parse("∫ x^2 dx");
    TEST_ASSERT(result1.isValid, "Failed to parse indefinite integral");
    TEST_ASSERT(result1.type == EquationType::INTEGRAL, "Wrong type for integral");
    
    // Test definite integral with limits
    auto result2 = parser.parse("∫[0,1] x dx");
    TEST_ASSERT(result2.isValid, "Failed to parse definite integral");
    
    std::cout << "test_integral_parsing PASSED" << std::endl;
    return true;
}

// Test Case 4: Series and Summation Parsing
bool test_series_parsing() {
    std::cout << "Running test_series_parsing..." << std::endl;
    ANTLR4Parser parser;
    
    // Test summation
    auto result1 = parser.parse("∑[i=1,10] i^2");
    TEST_ASSERT(result1.isValid, "Failed to parse summation");
    TEST_ASSERT(result1.type == EquationType::SERIES, "Wrong type for summation");
    
    // Test product
    auto result2 = parser.parse("∏[i=1,5] i");
    TEST_ASSERT(result2.isValid, "Failed to parse product");
    
    // Test series expansion
    auto result3 = parser.parse("series(exp(x), x, 5)");
    TEST_ASSERT(result3.isValid, "Failed to parse series expansion");
    TEST_ASSERT(result3.type == EquationType::SERIES, "Wrong type for series");
    
    std::cout << "test_series_parsing PASSED" << std::endl;
    return true;
}

// Test Case 5: Parametric Equation Parsing
bool test_parametric_parsing() {
    std::cout << "Running test_parametric_parsing..." << std::endl;
    ANTLR4Parser parser;
    
    // Test parametric equation
    auto result1 = parser.parse("x(t) = cos(t)");
    TEST_ASSERT(result1.isValid, "Failed to parse parametric equation");
    TEST_ASSERT(result1.type == EquationType::PARAMETRIC, "Wrong type for parametric");
    
    std::cout << "test_parametric_parsing PASSED" << std::endl;
    return true;
}

// Test Case 6: ODE Parsing
bool test_ode_parsing() {
    std::cout << "Running test_ode_parsing..." << std::endl;
    ANTLR4Parser parser;
    
    // Test ODE dy/dt = ...
    auto result1 = parser.parse("dy/dt = -k*y");
    TEST_ASSERT(result1.isValid, "Failed to parse ODE");
    TEST_ASSERT(result1.type == EquationType::ODE, "Wrong type for ODE");
    
    std::cout << "test_ode_parsing PASSED" << std::endl;
    return true;
}

// Test Case 7: UQFF Equation Detection
bool test_uqff_detection() {
    std::cout << "Running test_uqff_detection..." << std::endl;
    ANTLR4Parser parser;
    
    // Test F_U_Bi_i detection
    auto result1 = parser.parseUQFF("F_U_Bi_i(r, t, M)");
    TEST_ASSERT(result1.isValid, "Failed to parse F_U_Bi_i");
    TEST_ASSERT(result1.isUQFF, "Failed to detect UQFF pattern");
    
    // Test Um detection
    auto result2 = parser.parseUQFF("Um(t, r, n)");
    TEST_ASSERT(result2.isValid, "Failed to parse Um");
    TEST_ASSERT(result2.isUQFF, "Failed to detect Um UQFF pattern");
    
    // Test g_MUGE detection
    auto result3 = parser.parseUQFF("g_MUGE_H(r, t)");
    TEST_ASSERT(result3.isValid, "Failed to parse g_MUGE_H");
    TEST_ASSERT(result3.isUQFF, "Failed to detect g_MUGE UQFF pattern");
    
    std::cout << "test_uqff_detection PASSED" << std::endl;
    return true;
}

// Test Case 8: Error Handling
bool test_error_handling() {
    std::cout << "Running test_error_handling..." << std::endl;
    ANTLR4Parser parser;
    
    // Test invalid syntax
    auto result1 = parser.parse("2 +* 3");
    TEST_ASSERT(!result1.isValid, "Should reject invalid syntax '2 +* 3'");
    TEST_ASSERT(!result1.errorMessage.empty(), "Should provide error message");
    
    // Test mismatched parentheses
    auto result2 = parser.parse("sin(x");
    TEST_ASSERT(!result2.isValid, "Should reject mismatched parentheses");
    
    // Test empty input
    auto result3 = parser.parse("");
    TEST_ASSERT(!result3.isValid, "Should reject empty input");
    
    std::cout << "test_error_handling PASSED" << std::endl;
    return true;
}

// Test Case 9: Complex Expressions
bool test_complex_expressions() {
    std::cout << "Running test_complex_expressions..." << std::endl;
    ANTLR4Parser parser;
    
    // Test nested functions
    auto result1 = parser.parse("exp(sin(x^2 + 1))");
    TEST_ASSERT(result1.isValid, "Failed to parse nested functions");
    
    // Test multiple operations
    auto result2 = parser.parse("(2*x + 3)/(x^2 - 1)");
    TEST_ASSERT(result2.isValid, "Failed to parse rational expression");
    
    // Test scientific notation
    auto result3 = parser.parse("1.23e-45 * x");
    TEST_ASSERT(result3.isValid, "Failed to parse scientific notation");
    
    std::cout << "test_complex_expressions PASSED" << std::endl;
    return true;
}

// Test Case 10: Validation Method
bool test_validation() {
    std::cout << "Running test_validation..." << std::endl;
    ANTLR4Parser parser;
    
    // Test valid expression validation
    bool valid1 = parser.validate("x^2 + 1");
    TEST_ASSERT(valid1, "Should validate 'x^2 + 1'");
    
    // Test invalid expression validation
    bool valid2 = parser.validate("2 +* 3");
    TEST_ASSERT(!valid2, "Should invalidate '2 +* 3'");
    
    // Test getLastError
    parser.validate("invalid**syntax");
    std::string error = parser.getLastError();
    TEST_ASSERT(!error.empty(), "Should return error message");
    
    std::cout << "test_validation PASSED" << std::endl;
    return true;
}

// Main test runner
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "ANTLR4 Parser Test Suite" << std::endl;
    std::cout << "calculator_advanced framework" << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::vector<bool (*)(void)> tests = {
        test_basic_expressions,
        test_derivative_parsing,
        test_integral_parsing,
        test_series_parsing,
        test_parametric_parsing,
        test_ode_parsing,
        test_uqff_detection,
        test_error_handling,
        test_complex_expressions,
        test_validation
    };
    
    int passed = 0;
    int failed = 0;
    
    for (auto test : tests) {
        try {
            if (test()) {
                passed++;
            } else {
                failed++;
            }
        } catch (const std::exception& e) {
            std::cerr << "EXCEPTION: " << e.what() << std::endl;
            failed++;
        }
        std::cout << std::endl;
    }
    
    std::cout << "========================================" << std::endl;
    std::cout << "Results: " << passed << " passed, " << failed << " failed" << std::endl;
    std::cout << "========================================" << std::endl;
    
    return (failed == 0) ? 0 : 1;
}
