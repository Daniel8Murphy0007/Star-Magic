// Test suite for Polynomial Solver (GSL wrapper)
// Tests polynomial root finding from degree 1 to 26
// Part of calculator_advanced framework extracted from Grok thread Iterations #30-32

#include "../include/polynomial_solver.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

// Test helper macro
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        std::cerr << "TEST FAILED: " << message << std::endl; \
        return false; \
    }

// Helper function to check if value is approximately zero
bool isNearZero(double value, double epsilon = 1e-8) {
    return std::abs(value) < epsilon;
}

// Test Case 1: Linear Equation (degree 1)
bool test_linear_equation() {
    std::cout << "Running test_linear_equation..." << std::endl;
    
    // Test: x + 2 = 0, solution: x = -2
    Polynomial poly;
    poly.coefficients = {2.0, 1.0};  // 2 + 1*x
    poly.degree = 1;
    
    auto roots = solve(poly);
    std::cout << "  Roots of x+2=0: ";
    for (const auto& root : roots) {
        std::cout << root << " ";
    }
    std::cout << std::endl;
    
    TEST_ASSERT(roots.size() == 1, "Linear equation should have 1 root");
    TEST_ASSERT(std::abs(roots[0].real() + 2.0) < 1e-8, "Root should be -2");
    TEST_ASSERT(isNearZero(roots[0].imag()), "Root should be real");
    
    std::cout << "test_linear_equation PASSED" << std::endl;
    return true;
}

// Test Case 2: Quadratic Equation (degree 2)
bool test_quadratic_equation() {
    std::cout << "Running test_quadratic_equation..." << std::endl;
    
    // Test: x^2 - 5*x + 6 = 0, solutions: x = 2, 3
    Polynomial poly;
    poly.coefficients = {6.0, -5.0, 1.0};  // 6 - 5*x + 1*x^2
    poly.degree = 2;
    
    auto roots = solve(poly);
    std::cout << "  Roots of x^2-5*x+6=0: ";
    for (const auto& root : roots) {
        std::cout << root << " ";
    }
    std::cout << std::endl;
    
    TEST_ASSERT(roots.size() == 2, "Quadratic equation should have 2 roots");
    
    std::cout << "test_quadratic_equation PASSED" << std::endl;
    return true;
}

// Test Case 3: Cubic Equation (degree 3)
bool test_cubic_equation() {
    std::cout << "Running test_cubic_equation..." << std::endl;
    
    // Test: x^3 - 6*x^2 + 11*x - 6 = 0, solutions: x = 1, 2, 3
    Polynomial poly;
    poly.coefficients = {-6.0, 11.0, -6.0, 1.0};  // -6 + 11*x - 6*x^2 + 1*x^3
    poly.degree = 3;
    
    auto roots = solve(poly);
    std::cout << "  Roots of x^3-6*x^2+11*x-6=0: ";
    for (const auto& root : roots) {
        std::cout << root << " ";
    }
    std::cout << std::endl;
    
    TEST_ASSERT(roots.size() == 3, "Cubic equation should have 3 roots");
    
    std::cout << "test_cubic_equation PASSED" << std::endl;
    return true;
}

// Test Case 4: Quartic Equation (degree 4)
bool test_quartic_equation() {
    std::cout << "Running test_quartic_equation..." << std::endl;
    
    // Test: x^4 - 1 = 0, solutions: x = ±1, ±i
    Polynomial poly;
    poly.coefficients = {-1.0, 0.0, 0.0, 0.0, 1.0};  // -1 + 0*x + 0*x^2 + 0*x^3 + 1*x^4
    poly.degree = 4;
    
    auto roots = solve(poly);
    std::cout << "  Roots of x^4-1=0: ";
    for (const auto& root : roots) {
        std::cout << root << " ";
    }
    std::cout << std::endl;
    
    TEST_ASSERT(roots.size() == 4, "Quartic equation should have 4 roots");
    
    std::cout << "test_quartic_equation PASSED" << std::endl;
    return true;
}

// Test Case 5: Real Root Filtering
bool test_real_root_filtering() {
    std::cout << "Running test_real_root_filtering..." << std::endl;
    
    // Test: x^4 - 1 = 0, real roots: x = ±1
    Polynomial poly;
    poly.coefficients = {-1.0, 0.0, 0.0, 0.0, 1.0};
    poly.degree = 4;
    
    auto realRoots = solveReal(poly);
    std::cout << "  Real roots of x^4-1=0: ";
    for (const auto& root : realRoots) {
        std::cout << root << " ";
    }
    std::cout << std::endl;
    
    TEST_ASSERT(realRoots.size() == 2, "Should have 2 real roots (±1)");
    
    std::cout << "test_real_root_filtering PASSED" << std::endl;
    return true;
}

// Test Case 6: High Degree Polynomial (degree 10)
bool test_high_degree_polynomial() {
    std::cout << "Running test_high_degree_polynomial..." << std::endl;
    
    // Test: x^10 - 1 = 0 (10th roots of unity)
    Polynomial poly;
    poly.coefficients.resize(11, 0.0);
    poly.coefficients[0] = -1.0;  // constant term
    poly.coefficients[10] = 1.0;  // x^10 term
    poly.degree = 10;
    
    auto roots = solve(poly);
    std::cout << "  Roots of x^10-1=0: " << roots.size() << " roots found" << std::endl;
    
    TEST_ASSERT(roots.size() == 10, "Should have 10 roots (10th roots of unity)");
    
    std::cout << "test_high_degree_polynomial PASSED" << std::endl;
    return true;
}

// Test Case 7: Maximum Degree Support (degree 26)
bool test_maximum_degree() {
    std::cout << "Running test_maximum_degree..." << std::endl;
    
    // Test: x^26 - 1 = 0 (26th roots of unity)
    Polynomial poly;
    poly.coefficients.resize(27, 0.0);
    poly.coefficients[0] = -1.0;  // constant term
    poly.coefficients[26] = 1.0;  // x^26 term
    poly.degree = 26;
    
    auto roots = solve(poly);
    std::cout << "  Roots of x^26-1=0: " << roots.size() << " roots found" << std::endl;
    
    TEST_ASSERT(roots.size() == 26, "Should have 26 roots");
    
    std::cout << "test_maximum_degree PASSED" << std::endl;
    return true;
}

// Test Case 8: Zero Constant Term
bool test_zero_constant() {
    std::cout << "Running test_zero_constant..." << std::endl;
    
    // Test: x^2 + x = 0, solutions: x = 0, -1
    Polynomial poly;
    poly.coefficients = {0.0, 1.0, 1.0};  // 0 + 1*x + 1*x^2
    poly.degree = 2;
    
    auto roots = solve(poly);
    std::cout << "  Roots of x^2+x=0: ";
    for (const auto& root : roots) {
        std::cout << root << " ";
    }
    std::cout << std::endl;
    
    TEST_ASSERT(roots.size() == 2, "Should have 2 roots");
    
    std::cout << "test_zero_constant PASSED" << std::endl;
    return true;
}

// Test Case 9: Complex Roots
bool test_complex_roots() {
    std::cout << "Running test_complex_roots..." << std::endl;
    
    // Test: x^2 + 1 = 0, solutions: x = ±i
    Polynomial poly;
    poly.coefficients = {1.0, 0.0, 1.0};  // 1 + 0*x + 1*x^2
    poly.degree = 2;
    
    auto roots = solve(poly);
    std::cout << "  Roots of x^2+1=0: ";
    for (const auto& root : roots) {
        std::cout << root << " ";
    }
    std::cout << std::endl;
    
    TEST_ASSERT(roots.size() == 2, "Should have 2 complex roots");
    
    // Check that real roots filtering returns empty
    auto realRoots = solveReal(poly);
    TEST_ASSERT(realRoots.size() == 0, "Should have no real roots");
    
    std::cout << "test_complex_roots PASSED" << std::endl;
    return true;
}

// Test Case 10: Multiple Real Roots
bool test_multiple_real_roots() {
    std::cout << "Running test_multiple_real_roots..." << std::endl;
    
    // Test: (x-1)^2 = x^2 - 2*x + 1 = 0, double root at x = 1
    Polynomial poly;
    poly.coefficients = {1.0, -2.0, 1.0};  // 1 - 2*x + 1*x^2
    poly.degree = 2;
    
    auto roots = solve(poly);
    std::cout << "  Roots of (x-1)^2=0: ";
    for (const auto& root : roots) {
        std::cout << root << " ";
    }
    std::cout << std::endl;
    
    TEST_ASSERT(roots.size() == 2, "Should detect 2 roots (double root)");
    
    std::cout << "test_multiple_real_roots PASSED" << std::endl;
    return true;
}

// Main test runner
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Polynomial Solver Test Suite (GSL)" << std::endl;
    std::cout << "calculator_advanced framework" << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::vector<bool (*)(void)> tests = {
        test_linear_equation,
        test_quadratic_equation,
        test_cubic_equation,
        test_quartic_equation,
        test_real_root_filtering,
        test_high_degree_polynomial,
        test_maximum_degree,
        test_zero_constant,
        test_complex_roots,
        test_multiple_real_roots
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
