// CorePenetrationModule.h
// Modular C++ implementation of the Planetary Core Penetration Factor (P_core) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes P_core ≈1 (unitless for Sun, ~1e-3 for planets); scales P_core in Universal Gravity U_g3 term.
// Pluggable: #include "CorePenetrationModule.h"
// CorePenetrationModule mod; mod.computeU_g3(0.0); mod.updateVariable("P_core", new_value);
// Variables in std::map; example for Sun at t=0; planet mode with P_core=1e-3.
// Approximations: cos(ω_s t π)=1 at t=0; E_react=1e46; B_j=1e3 T.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef CORE_PENETRATION_MODULE_H
#define CORE_PENETRATION_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class CorePenetrationModule {
private:
    std::map<std::string, double> variables;
    double computeU_g3(double t);

public:
    // Constructor: Initialize with framework defaults (Sun)
    CorePenetrationModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeP_core();  // ≈1 for Sun (unitless)
    double computeU_g3(double t);  // U_g3 with P_core (J/m^3)
    double computeU_g3_planet(double t);  // For planet P_core=1e-3

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // CORE_PENETRATION_MODULE_H
