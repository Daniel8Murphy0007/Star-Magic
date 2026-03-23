// GalacticDistanceModule.h
// Modular C++ implementation of the Distance from Galactic Center (d_g) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes d_g=2.55e20 m (~27,000 ly) and conversions; scales M_bh / d_g in U_bi and Ug4.
// Pluggable: #include "GalacticDistanceModule.h"
// GalacticDistanceModule mod; mod.computeMbhOverDg(); mod.updateVariable("d_g", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0.
// Approximations: cos(π t_n)=1; ε_sw * ρ_vac,sw ≈0; α=0.001 s^-1; f_feedback=0.1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef GALACTIC_DISTANCE_MODULE_H
#define GALACTIC_DISTANCE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class GalacticDistanceModule {
private:
    std::map<std::string, double> variables;
    double computeDgInLy();
    double computeDgInPc();
    double computeMbhOverDg();
    double computeU_b1();
    double computeU_g4();

public:
    // Constructor: Initialize with framework defaults
    GalacticDistanceModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeDg();  // d_g in m (2.55e20)
    double computeDgInLy();
    double computeDgInPc();
    double computeMbhOverDg();  // M_bh / d_g (kg/m)
    double computeU_b1();  // Universal Buoyancy example (J/m^3)
    double computeU_g4();  // Ug4 example (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // GALACTIC_DISTANCE_MODULE_H
