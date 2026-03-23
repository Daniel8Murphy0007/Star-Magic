// OuterFieldBubbleModule.h
// Modular C++ implementation of the Radius of the Outer Field Bubble (R_b) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes R_b=1.496e13 m (100 AU); defines S(r - R_b) step function in Universal Gravity U_g2 term.
// Pluggable: #include "OuterFieldBubbleModule.h"
// OuterFieldBubbleModule mod; mod.computeU_g2(1.5e13); mod.updateVariable("R_b", new_value);
// Variables in std::map; example for Sun at t=0; S=1 for r >= R_b, 0 otherwise.
// Approximations: S step=1 at r=R_b; δ_sw v_sw=5001; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef OUTER_FIELD_BUBBLE_MODULE_H
#define OUTER_FIELD_BUBBLE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class OuterFieldBubbleModule {
private:
    std::map<std::string, double> variables;
    double computeS_r_Rb(double r);
    double computeU_g2(double r);

public:
    // Constructor: Initialize with framework defaults (Sun)
    OuterFieldBubbleModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeR_b();  // 1.496e13 m (100 AU)
    double computeR_bInAU();  // 100 AU
    double computeS_r_Rb(double r);  // Step function
    double computeU_g2(double r);  // U_g2 with S(r - R_b) (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // OUTER_FIELD_BUBBLE_MODULE_H
