// GalacticSpinModule.h
// Modular C++ implementation of the Galactic Spin Rate (Ω_g) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes Ω_g=7.3e-16 rad/s and its role in Universal Buoyancy U_bi; scales Ω_g * (M_bh / d_g) in F_U.
// Pluggable: #include "GalacticSpinModule.h"
// GalacticSpinModule mod; mod.computeU_b1(); mod.updateVariable("Omega_g", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; integrates β_i, U_gi, etc.
// Approximations: cos(π t_n)=1; (1 + ε_sw ρ_vac,sw)≈1; U_UA=1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef GALACTIC_SPIN_MODULE_H
#define GALACTIC_SPIN_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class GalacticSpinModule {
private:
    std::map<std::string, double> variables;
    double computeOmega_g();
    double computeMbhOverDg();
    double computeU_bi(int i);

public:
    // Constructor: Initialize with framework defaults
    GalacticSpinModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeOmega_g();  // 7.3e-16 rad/s
    double computeMbhOverDg();  // M_bh / d_g (kg/m)
    double computeU_bi(int i);  // U_bi for i=1-4 (J/m^3, example i=1)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // GALACTIC_SPIN_MODULE_H
