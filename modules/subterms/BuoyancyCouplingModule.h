// BuoyancyCouplingModule.h
// Modular C++ implementation of the Buoyancy Coupling Constants (β_i) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes the Universal Buoyancy terms U_bi = -β_i * U_gi * Ω_g * (M_bh / d_g) * E_react for i=1 to 4 (Ug1-Ug4).
// Pluggable: #include "BuoyancyCouplingModule.h"
// BuoyancyCouplingModule mod; mod.computeU_bi(1); mod.updateVariable("beta", new_value);
// Variables in std::map; β_i=0.6 uniform (unitless); opposes gravity with 60% scaling.
// Approximations: cos(π t_n)=1 at t_n=0; ε_sw * ρ_vac,sw ≈0; U_UA=1; computes per i or sum.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef BUOYANCY_COUPLING_MODULE_H
#define BUOYANCY_COUPLING_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class BuoyancyCouplingModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> computeAllU_bi();

public:
    // Constructor: Initialize with framework defaults
    BuoyancyCouplingModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeBeta(int i);  // β_i = 0.6 for all i
    double computeU_bi(int i);  // U_bi for specific i (Ug1-4)
    std::vector<double> computeAllU_bi();  // All four U_bi
    double computeF_U_contribution();  // Sum β_i terms in F_U

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print all U_bi
    void printU_bi();
};

#endif // BUOYANCY_COUPLING_MODULE_H
