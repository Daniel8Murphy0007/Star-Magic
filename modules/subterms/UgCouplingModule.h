// UgCouplingModule.h
// Modular C++ implementation of the Coupling Constants for Ug Ranges (k_i) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes scaled Universal Gravity terms k_i * U_gi for i=1-4 (Ug1-Ug4), with k1=1.5, k2=1.2, k3=1.8, k4=1.0 (unitless).
// Pluggable: #include "UgCouplingModule.h"
// UgCouplingModule mod; mod.computeSumK_Ugi(); mod.updateVariable("U_g1", new_value);
// Variables in std::map; sum contributes to F_U; placeholders for full U_gi equations.
// Approximations: t_n=0, cos(π t_n)=1; δ_def=0, etc.; example values from Sun at t=0.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UG_COUPLING_MODULE_H
#define UG_COUPLING_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class UgCouplingModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> k_values;  // [k1, k2, k3, k4]
    std::vector<double> computeAllK_Ugi();

public:
    // Constructor: Initialize with framework defaults
    UgCouplingModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeK_i(int i);  // k_i for specific i (1-4)
    double computeU_gi(int i);  // Placeholder U_gi (J/m^3)
    double computeK_Ugi(int i);  // k_i * U_gi
    std::vector<double> computeAllK_Ugi();  // All four k_i * U_gi
    double computeSumK_Ugi();  // Sum for F_U contribution

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print all k_i * U_gi
    void printK_Ugi();
};

#endif // UG_COUPLING_MODULE_H
