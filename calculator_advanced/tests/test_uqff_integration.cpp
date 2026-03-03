// Test suite for UQFF Equation Integration
// Tests all 11 UQFF equations with sample parameters
// Part of calculator_advanced framework extracted from Grok thread Iterations #30-32

#include "../include/uqff_equations.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>

// Test helper macro
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        std::cerr << "TEST FAILED: " << message << std::endl; \
        return false; \
    }

// Test Case 1: F_U_Bi_i Calculation
bool test_F_U_Bi_i() {
    std::cout << "Running test_F_U_Bi_i..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    // Test parameters for universal buoyancy force
    std::map<std::string, double> params;
    params["F_rel"] = 4.30e33;      // N (LEP baseline)
    params["E_cm"] = 13000.0;       // GeV (LHC energy)
    params["E_LEP"] = 200.0;        // GeV
    params["Q_wave"] = 1e12;        // Wave factor
    params["M"] = 1.67e-27;         // kg (proton mass)
    params["r"] = 1e-15;            // m (nuclear radius)
    params["t"] = 1e-21;            // s
    
    double result = catalog.calculate("F_U_Bi_i", params);
    std::cout << "  F_U_Bi_i = " << result << " N" << std::endl;
    
    TEST_ASSERT(!std::isnan(result), "F_U_Bi_i should not be NaN");
    TEST_ASSERT(!std::isinf(result), "F_U_Bi_i should not be infinite");
    TEST_ASSERT(result > 0, "F_U_Bi_i should be positive for these params");
    
    std::cout << "test_F_U_Bi_i PASSED" << std::endl;
    return true;
}

// Test Case 2: Um Time-Dependent Magnetism
bool test_Um() {
    std::cout << "Running test_Um..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    // Test parameters for universal magnetism
    std::map<std::string, double> params;
    params["t"] = 1e7;              // s (115 days)
    params["r"] = 1e10;             // m
    params["n"] = 1.0;              // Harmonic order
    params["gamma"] = 5e-5;         // day^-1
    params["phi"] = 1.0;            // Phase factor
    params["P_SCm"] = 1.0;          // SCm pressure
    params["f_Heav"] = 0.01;        // Heaviside factor
    params["f_quasi"] = 0.01;       // Quasi factor
    
    double result = catalog.calculate("Um", params);
    std::cout << "  Um = " << result << " T·m (at t=115 days)" << std::endl;
    
    TEST_ASSERT(!std::isnan(result), "Um should not be NaN");
    TEST_ASSERT(!std::isinf(result), "Um should not be infinite");
    
    std::cout << "test_Um PASSED" << std::endl;
    return true;
}

// Test Case 3: g_MUGE_H (Hydrogen Atom)
bool test_g_MUGE_H() {
    std::cout << "Running test_g_MUGE_H..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    // Test parameters for hydrogen atom MUGE
    std::map<std::string, double> params;
    params["r"] = 5.29e-11;         // m (Bohr radius)
    params["t"] = 0.0;              // s
    params["m_eff"] = 1.67e-27;     // kg (proton mass)
    
    double result = catalog.calculate("g_MUGE_H", params);
    std::cout << "  g_MUGE_H = " << result << " m/s^2 (at Bohr radius)" << std::endl;
    
    TEST_ASSERT(!std::isnan(result), "g_MUGE_H should not be NaN");
    TEST_ASSERT(!std::isinf(result), "g_MUGE_H should not be infinite");
    TEST_ASSERT(result > 0, "g_MUGE_H should be positive");
    
    std::cout << "test_g_MUGE_H PASSED" << std::endl;
    return true;
}

// Test Case 4: g_Magnetar (Field Decay)
bool test_g_Magnetar() {
    std::cout << "Running test_g_Magnetar..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    // Test parameters for magnetar gravity
    std::map<std::string, double> params;
    params["M"] = 3.0e30;           // kg (1.5 solar masses)
    params["r"] = 1e4;              // m (10 km)
    params["t"] = 3.15e7;           // s (1 year)
    params["B"] = 1e13;             // T (field strength)
    
    double result = catalog.calculate("g_Magnetar", params);
    std::cout << "  g_Magnetar = " << result << " m/s^2" << std::endl;
    
    TEST_ASSERT(!std::isnan(result), "g_Magnetar should not be NaN");
    TEST_ASSERT(!std::isinf(result), "g_Magnetar should not be infinite");
    TEST_ASSERT(result > 0, "g_Magnetar should be positive");
    
    std::cout << "test_g_Magnetar PASSED" << std::endl;
    return true;
}

// Test Case 5: g_SgrA (Sgr A* Evolution)
bool test_g_SgrA() {
    std::cout << "Running test_g_SgrA..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    // Test parameters for Sgr A*
    std::map<std::string, double> params;
    params["r"] = 1e11;             // m (~1 AU)
    params["t"] = 1e9 * 3.15e7;     // s (1 Gyr)
    params["B"] = 0.1;              // T (weak field at 1 AU)
    
    double result = catalog.calculate("g_SgrA", params);
    std::cout << "  g_SgrA = " << result << " m/s^2" << std::endl;
    
    TEST_ASSERT(!std::isnan(result), "g_SgrA should not be NaN");
    TEST_ASSERT(!std::isinf(result), "g_SgrA should not be infinite");
    
    std::cout << "test_g_SgrA PASSED" << std::endl;
    return true;
}

// Test Case 6: P_alpha (Clustering Probability)
bool test_P_alpha() {
    std::cout << "Running test_P_alpha..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    // Test parameters for alpha clustering
    std::map<std::string, double> params;
    params["F_U_Bi_i"] = 1e-12;     // J (above threshold)
    
    double result = catalog.calculate("P_alpha", params);
    std::cout << "  P_alpha = " << result << " (probability)" << std::endl;
    
    TEST_ASSERT(!std::isnan(result), "P_alpha should not be NaN");
    TEST_ASSERT(result >= 0.0 && result <= 1.0, "P_alpha should be between 0 and 1");
    
    std::cout << "test_P_alpha PASSED" << std::endl;
    return true;
}

// Test Case 7: R_EU (Electric Universe Ratio)
bool test_R_EU() {
    std::cout << "Running test_R_EU..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    // Test parameters for EU ratio
    std::map<std::string, double> params;
    params["q"] = 1.6e-19;          // C (elementary charge)
    params["Um"] = 1e20;            // T·m
    params["rho_vac"] = 7.09e-36;   // J/m^3
    params["v"] = 1e5;              // m/s
    params["r"] = 1e11;             // m
    params["M"] = 1.989e30;         // kg (solar mass)
    params["m"] = 1.67e-27;         // kg (proton mass)
    
    double result = catalog.calculate("R_EU", params);
    std::cout << "  R_EU = " << result << " (F_EM/F_g)" << std::endl;
    
    TEST_ASSERT(!std::isnan(result), "R_EU should not be NaN");
    TEST_ASSERT(!std::isinf(result), "R_EU should not be infinite");
    
    std::cout << "test_R_EU PASSED" << std::endl;
    return true;
}

// Test Case 8: tau_gyro (Gyroscopic Torque)
bool test_tau_gyro() {
    std::cout << "Running test_tau_gyro..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    // Test parameters for gyroscopic torque
    std::map<std::string, double> params;
    params["I"] = 1e40;             // kg·m^2 (moment of inertia)
    params["omega"] = 1e-5;         // rad/s (angular velocity)
    params["alpha"] = 1e-10;        // rad/s^2 (angular acceleration)
    
    double result = catalog.calculate("tau_gyro", params);
    std::cout << "  tau_gyro = " << result << " N·m" << std::endl;
    
    TEST_ASSERT(!std::isnan(result), "tau_gyro should not be NaN");
    TEST_ASSERT(!std::isinf(result), "tau_gyro should not be infinite");
    
    std::cout << "test_tau_gyro PASSED" << std::endl;
    return true;
}

// Test Case 9: g_compressed (26-Layer)
bool test_g_compressed() {
    std::cout << "Running test_g_compressed..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    // Test parameters for 26-layer compressed gravity
    std::map<std::string, double> params;
    params["M"] = 1.989e30;         // kg (solar mass)
    params["r"] = 1e11;             // m (1 AU)
    params["t"] = 0.0;              // s
    
    double result = catalog.calculate("g_compressed", params);
    std::cout << "  g_compressed (26 layers) = " << result << " m/s^2" << std::endl;
    
    TEST_ASSERT(!std::isnan(result), "g_compressed should not be NaN");
    TEST_ASSERT(!std::isinf(result), "g_compressed should not be infinite");
    
    std::cout << "test_g_compressed PASSED" << std::endl;
    return true;
}

// Test Case 10: eta_LENR (Neutron Production)
bool test_eta_LENR() {
    std::cout << "Running test_eta_LENR..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    // Test parameters for LENR neutron rate
    std::map<std::string, double> params;
    params["SSq"] = 0.57;           // Calibration constant
    params["n"] = 1.0;              // Harmonic order
    params["t"] = 1e4;              // s
    params["Um"] = 1e20;            // T·m
    
    double result = catalog.calculate("eta_LENR", params);
    std::cout << "  eta_LENR = " << result << " s^-1" << std::endl;
    
    TEST_ASSERT(!std::isnan(result), "eta_LENR should not be NaN");
    TEST_ASSERT(!std::isinf(result), "eta_LENR should not be infinite");
    TEST_ASSERT(result >= 0, "eta_LENR should be non-negative");
    
    std::cout << "test_eta_LENR PASSED" << std::endl;
    return true;
}

// Test Case 11: Equation Listing
bool test_equation_listing() {
    std::cout << "Running test_equation_listing..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    auto equations = catalog.listEquations();
    std::cout << "  Available equations: " << equations.size() << std::endl;
    for (const auto& eq : equations) {
        std::cout << "    - " << eq << std::endl;
    }
    
    TEST_ASSERT(equations.size() >= 10, "Should have at least 10 equations");
    
    std::cout << "test_equation_listing PASSED" << std::endl;
    return true;
}

// Test Case 12: Equation Search
bool test_equation_search() {
    std::cout << "Running test_equation_search..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    // Search for gravity equations
    auto gravityEqs = catalog.searchEquations("gravity");
    std::cout << "  Gravity equations found: " << gravityEqs.size() << std::endl;
    for (const auto& eq : gravityEqs) {
        std::cout << "    - " << eq << std::endl;
    }
    
    TEST_ASSERT(gravityEqs.size() > 0, "Should find gravity equations");
    
    // Search for magnetism equations
    auto magnetismEqs = catalog.searchEquations("magnetism");
    std::cout << "  Magnetism equations found: " << magnetismEqs.size() << std::endl;
    
    TEST_ASSERT(magnetismEqs.size() > 0, "Should find magnetism equations");
    
    std::cout << "test_equation_search PASSED" << std::endl;
    return true;
}

// Test Case 13: Invalid Equation Name
bool test_invalid_equation() {
    std::cout << "Running test_invalid_equation..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    std::map<std::string, double> params;
    
    bool threw_exception = false;
    try {
        catalog.calculate("InvalidEquation", params);
    } catch (const std::exception& e) {
        threw_exception = true;
        std::cout << "  Expected exception: " << e.what() << std::endl;
    }
    
    TEST_ASSERT(threw_exception, "Should throw for invalid equation name");
    
    std::cout << "test_invalid_equation PASSED" << std::endl;
    return true;
}

// Test Case 14: Missing Parameters
bool test_missing_parameters() {
    std::cout << "Running test_missing_parameters..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    // Try to calculate F_U_Bi_i with missing parameters
    std::map<std::string, double> params;
    params["F_rel"] = 4.30e33;  // Only one parameter
    
    bool threw_exception = false;
    try {
        catalog.calculate("F_U_Bi_i", params);
    } catch (const std::exception& e) {
        threw_exception = true;
        std::cout << "  Expected exception: " << e.what() << std::endl;
    }
    
    TEST_ASSERT(threw_exception, "Should throw for missing parameters");
    
    std::cout << "test_missing_parameters PASSED" << std::endl;
    return true;
}

// Test Case 15: Equation Retrieval
bool test_equation_retrieval() {
    std::cout << "Running test_equation_retrieval..." << std::endl;
    
    UQFFEquationCatalog catalog;
    
    auto eq = catalog.getEquation("F_U_Bi_i");
    std::cout << "  F_U_Bi_i equation: " << eq.name << std::endl;
    std::cout << "  Description: " << eq.description << std::endl;
    std::cout << "  LaTeX: " << eq.latexFormula << std::endl;
    
    TEST_ASSERT(eq.name == "F_U_Bi_i", "Should retrieve correct equation");
    TEST_ASSERT(!eq.description.empty(), "Should have description");
    TEST_ASSERT(!eq.latexFormula.empty(), "Should have LaTeX formula");
    
    std::cout << "test_equation_retrieval PASSED" << std::endl;
    return true;
}

// Main test runner
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "UQFF Integration Test Suite" << std::endl;
    std::cout << "calculator_advanced framework" << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::vector<bool (*)(void)> tests = {
        test_F_U_Bi_i,
        test_Um,
        test_g_MUGE_H,
        test_g_Magnetar,
        test_g_SgrA,
        test_P_alpha,
        test_R_EU,
        test_tau_gyro,
        test_g_compressed,
        test_eta_LENR,
        test_equation_listing,
        test_equation_search,
        test_invalid_equation,
        test_missing_parameters,
        test_equation_retrieval
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
