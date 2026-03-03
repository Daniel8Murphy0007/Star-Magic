// Test suite for SymEngine Wrapper
// Tests symbolic differentiation, integration, simplification, solving
// Part of calculator_advanced framework extracted from Grok thread Iterations #30-32

#include "../include/symengine_wrapper.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

// Test helper macro
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        std::cerr << "TEST FAILED: " << message << std::endl; \
        return false; \
    }

// Test Case 1: Basic Expression Construction
bool test_expression_construction() {
    std::cout << "Running test_expression_construction..." << std::endl;
    
    // Test simple polynomial
    SymbolicExpression expr1("x^2 + 2*x + 1");
    TEST_ASSERT(!expr1.toString().empty(), "Failed to construct 'x^2 + 2*x + 1'");
    
    // Test trigonometric expression
    SymbolicExpression expr2("sin(x) + cos(x)");
    TEST_ASSERT(!expr2.toString().empty(), "Failed to construct 'sin(x) + cos(x)'");
    
    // Test exponential
    SymbolicExpression expr3("exp(x)");
    TEST_ASSERT(!expr3.toString().empty(), "Failed to construct 'exp(x)'");
    
    std::cout << "test_expression_construction PASSED" << std::endl;
    return true;
}

// Test Case 2: Differentiation
bool test_differentiation() {
    std::cout << "Running test_differentiation..." << std::endl;
    
    // Test polynomial derivative
    SymbolicExpression expr1("x^2");
    auto deriv1 = expr1.differentiate("x");
    std::string result1 = deriv1.toString();
    std::cout << "  d/dx(x^2) = " << result1 << std::endl;
    TEST_ASSERT(!result1.empty(), "Failed to differentiate x^2");
    
    // Test sin derivative
    SymbolicExpression expr2("sin(x)");
    auto deriv2 = expr2.differentiate("x");
    std::string result2 = deriv2.toString();
    std::cout << "  d/dx(sin(x)) = " << result2 << std::endl;
    TEST_ASSERT(!result2.empty(), "Failed to differentiate sin(x)");
    
    // Test product rule
    SymbolicExpression expr3("x*exp(x)");
    auto deriv3 = expr3.differentiate("x");
    std::string result3 = deriv3.toString();
    std::cout << "  d/dx(x*exp(x)) = " << result3 << std::endl;
    TEST_ASSERT(!result3.empty(), "Failed to differentiate x*exp(x)");
    
    std::cout << "test_differentiation PASSED" << std::endl;
    return true;
}

// Test Case 3: Integration
bool test_integration() {
    std::cout << "Running test_integration..." << std::endl;
    
    // Test polynomial integration
    SymbolicExpression expr1("x");
    auto integ1 = expr1.integrate("x");
    std::string result1 = integ1.toString();
    std::cout << "  ∫ x dx = " << result1 << std::endl;
    TEST_ASSERT(!result1.empty(), "Failed to integrate x");
    
    // Test power integration
    SymbolicExpression expr2("x^2");
    auto integ2 = expr2.integrate("x");
    std::string result2 = integ2.toString();
    std::cout << "  ∫ x^2 dx = " << result2 << std::endl;
    TEST_ASSERT(!result2.empty(), "Failed to integrate x^2");
    
    // Test constant integration
    SymbolicExpression expr3("5");
    auto integ3 = expr3.integrate("x");
    std::string result3 = integ3.toString();
    std::cout << "  ∫ 5 dx = " << result3 << std::endl;
    TEST_ASSERT(!result3.empty(), "Failed to integrate constant");
    
    std::cout << "test_integration PASSED" << std::endl;
    return true;
}

// Test Case 4: Simplification
bool test_simplification() {
    std::cout << "Running test_simplification..." << std::endl;
    
    // Test algebraic simplification
    SymbolicExpression expr1("x + x");
    auto simplified1 = expr1.simplify();
    std::string result1 = simplified1.toString();
    std::cout << "  simplify(x + x) = " << result1 << std::endl;
    TEST_ASSERT(!result1.empty(), "Failed to simplify x + x");
    
    // Test trigonometric simplification
    SymbolicExpression expr2("sin(x)^2 + cos(x)^2");
    auto simplified2 = expr2.simplify();
    std::string result2 = simplified2.toString();
    std::cout << "  simplify(sin^2 + cos^2) = " << result2 << std::endl;
    TEST_ASSERT(!result2.empty(), "Failed to simplify trig identity");
    
    std::cout << "test_simplification PASSED" << std::endl;
    return true;
}

// Test Case 5: Expansion
bool test_expansion() {
    std::cout << "Running test_expansion..." << std::endl;
    
    // Test polynomial expansion
    SymbolicExpression expr1("(x + 1)^2");
    auto expanded1 = expr1.expand();
    std::string result1 = expanded1.toString();
    std::cout << "  expand((x+1)^2) = " << result1 << std::endl;
    TEST_ASSERT(!result1.empty(), "Failed to expand (x+1)^2");
    
    // Test product expansion
    SymbolicExpression expr2("(x + 1)*(x - 1)");
    auto expanded2 = expr2.expand();
    std::string result2 = expanded2.toString();
    std::cout << "  expand((x+1)*(x-1)) = " << result2 << std::endl;
    TEST_ASSERT(!result2.empty(), "Failed to expand product");
    
    std::cout << "test_expansion PASSED" << std::endl;
    return true;
}

// Test Case 6: Substitution
bool test_substitution() {
    std::cout << "Running test_substitution..." << std::endl;
    
    // Test numeric substitution
    SymbolicExpression expr1("x^2 + 1");
    SymbolicExpression value("2");
    auto substituted = expr1.substitute("x", value);
    std::string result = substituted.toString();
    std::cout << "  x^2+1 at x=2: " << result << std::endl;
    TEST_ASSERT(!result.empty(), "Failed to substitute x=2");
    
    std::cout << "test_substitution PASSED" << std::endl;
    return true;
}

// Test Case 7: Evaluation
bool test_evaluation() {
    std::cout << "Running test_evaluation..." << std::endl;
    
    // Test numeric evaluation
    SymbolicExpression expr1("x^2 + 2*x + 1");
    std::map<std::string, double> values1 = {{"x", 2.0}};
    double result1 = expr1.evaluate(values1);
    std::cout << "  Evaluate x^2+2*x+1 at x=2: " << result1 << std::endl;
    TEST_ASSERT(std::abs(result1 - 9.0) < 1e-10, "Wrong evaluation of x^2+2*x+1 at x=2");
    
    // Test trigonometric evaluation
    SymbolicExpression expr2("sin(x)");
    std::map<std::string, double> values2 = {{"x", 0.0}};
    double result2 = expr2.evaluate(values2);
    std::cout << "  Evaluate sin(x) at x=0: " << result2 << std::endl;
    TEST_ASSERT(std::abs(result2) < 1e-10, "Wrong evaluation of sin(0)");
    
    std::cout << "test_evaluation PASSED" << std::endl;
    return true;
}

// Test Case 8: LaTeX Conversion
bool test_latex_conversion() {
    std::cout << "Running test_latex_conversion..." << std::endl;
    
    // Test polynomial LaTeX
    SymbolicExpression expr1("x^2 + 2*x + 1");
    std::string latex1 = expr1.toLatex();
    std::cout << "  LaTeX(x^2+2*x+1): " << latex1 << std::endl;
    TEST_ASSERT(!latex1.empty(), "Failed to generate LaTeX for polynomial");
    
    // Test fraction LaTeX
    SymbolicExpression expr2("1/x");
    std::string latex2 = expr2.toLatex();
    std::cout << "  LaTeX(1/x): " << latex2 << std::endl;
    TEST_ASSERT(!latex2.empty(), "Failed to generate LaTeX for fraction");
    
    std::cout << "test_latex_conversion PASSED" << std::endl;
    return true;
}

// Test Case 9: Symbolic Solver - Single Equation
bool test_solving_single_equation() {
    std::cout << "Running test_solving_single_equation..." << std::endl;
    
    SymbolicSolver solver;
    
    // Test linear equation: x + 1 = 0
    SymbolicExpression eq1("x + 1");
    auto solutions1 = solver.solve(eq1, "x");
    std::cout << "  Solutions to x+1=0: " << solutions1.size() << " found" << std::endl;
    TEST_ASSERT(solutions1.size() > 0, "Failed to solve x+1=0");
    
    // Test quadratic equation: x^2 - 1 = 0
    SymbolicExpression eq2("x^2 - 1");
    auto solutions2 = solver.solve(eq2, "x");
    std::cout << "  Solutions to x^2-1=0: " << solutions2.size() << " found" << std::endl;
    TEST_ASSERT(solutions2.size() > 0, "Failed to solve x^2-1=0");
    
    std::cout << "test_solving_single_equation PASSED" << std::endl;
    return true;
}

// Test Case 10: Taylor Series Expansion
bool test_taylor_series() {
    std::cout << "Running test_taylor_series..." << std::endl;
    
    SymbolicSolver solver;
    
    // Test exp(x) series
    SymbolicExpression func1("exp(x)");
    auto series1 = solver.taylorSeries(func1, "x", 3, 0.0);
    std::string result1 = series1.toString();
    std::cout << "  Taylor(exp(x), x=0, order=3): " << result1 << std::endl;
    TEST_ASSERT(!result1.empty(), "Failed to compute Taylor series of exp(x)");
    
    // Test sin(x) series
    SymbolicExpression func2("sin(x)");
    auto series2 = solver.taylorSeries(func2, "x", 3, 0.0);
    std::string result2 = series2.toString();
    std::cout << "  Taylor(sin(x), x=0, order=3): " << result2 << std::endl;
    TEST_ASSERT(!result2.empty(), "Failed to compute Taylor series of sin(x)");
    
    std::cout << "test_taylor_series PASSED" << std::endl;
    return true;
}

// Test Case 11: Critical Points
bool test_critical_points() {
    std::cout << "Running test_critical_points..." << std::endl;
    
    SymbolicSolver solver;
    
    // Test parabola critical point: f(x) = x^2
    SymbolicExpression func1("x^2");
    auto criticals1 = solver.findCriticalPoints(func1, "x");
    std::cout << "  Critical points of x^2: " << criticals1.size() << " found" << std::endl;
    TEST_ASSERT(criticals1.size() > 0, "Failed to find critical point of x^2");
    
    std::cout << "test_critical_points PASSED" << std::endl;
    return true;
}

// Main test runner
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "SymEngine Wrapper Test Suite" << std::endl;
    std::cout << "calculator_advanced framework" << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::vector<bool (*)(void)> tests = {
        test_expression_construction,
        test_differentiation,
        test_integration,
        test_simplification,
        test_expansion,
        test_substitution,
        test_evaluation,
        test_latex_conversion,
        test_solving_single_equation,
        test_taylor_series,
        test_critical_points
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
