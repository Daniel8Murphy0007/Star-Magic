// HeliosphereThicknessModule.h
// Modular C++ implementation of the Heliosphere Thickness Factor (H_SCm) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes H_SCm ≈1 (unitless) and its scaling in Universal Gravity U_g2 term.
// Pluggable: #include "HeliosphereThicknessModule.h"
// HeliosphereThicknessModule mod; mod.computeU_g2(0.0, 0.0); mod.updateVariable("H_SCm", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0, r=R_b=1.496e13 m.
// Approximations: S(r - R_b)=1; δ_sw v_sw=5001; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef HELIOSPHERE_THICKNESS_MODULE_H
#define HELIOSPHERE_THICKNESS_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class HeliosphereThicknessModule {
private:
    std::map<std::string, double> variables;
    double computeH_SCm();
    double computeU_g2(double t, double t_n);

public:
    // Constructor: Initialize with framework defaults
    HeliosphereThicknessModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeH_SCm();  // ≈1 (unitless)
    double computeU_g2(double t, double t_n);  // U_g2 with H_SCm (J/m^3)
    double computeU_g2_no_H(double t, double t_n);  // Without H_SCm variation

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // HELIOSPHERE_THICKNESS_MODULE_H
