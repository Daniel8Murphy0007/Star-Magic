// StepFunctionModule.h
// Modular C++ implementation of the Step Function (S) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes S(r - R_b) = 1 for r > R_b, 0 otherwise; activates U_g2 outside outer field bubble.
// Pluggable: #include "StepFunctionModule.h"
// StepFunctionModule mod; mod.computeU_g2(1.5e13); mod.updateVariable("R_b", new_value);
// Variables in std::map; example for Sun at r=1.496e13 m (R_b: S=1, U_g2≈1.18e53 J/m³); r=1e11 m: S=0, U_g2=0.
// Approximations: S=1 at r=R_b; (1 + δ_sw v_sw)=5001; H_SCm=1; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef STEP_FUNCTION_MODULE_H
#define STEP_FUNCTION_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class StepFunctionModule {
private:
    std::map<std::string, double> variables;
    double computeS_r_Rb(double r);
    double computeU_g2(double r);

public:
    // Constructor: Initialize with framework defaults (Sun)
    StepFunctionModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeS_r_Rb(double r);  // Step: 1 if r > R_b, 0 otherwise
    double computeU_g2(double r);  // U_g2 with S(r - R_b) (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // STEP_FUNCTION_MODULE_H
