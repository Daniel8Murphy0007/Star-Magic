// UgIndexModule.h
// Modular C++ implementation of the Index for Discrete Universal Gravity Ranges (i) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module uses i=1 to 4 to label Ug1-Ug4; computes sum_{i=1}^4 k_i * U_gi for F_U contribution.
// Pluggable: #include "UgIndexModule.h"
// UgIndexModule mod; mod.computeSumKUgi(); mod.updateVariable("U_g1", new_value);
// Variables in std::map; defaults for Sun at t=0; i labels: 1=Internal Dipole, 2=Outer Bubble, 3=Magnetic Disk, 4=Star-BH.
// Approximations: k_i from coupling; sum ≈1.42e53 J/m³ (Ug2 dominant).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UG_INDEX_MODULE_H
#define UG_INDEX_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class UgIndexModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> k_values;  // [k1=1.5, k2=1.2, k3=1.8, k4=1.0]
    std::vector<double> computeAllKUgi();

public:
    // Constructor: Initialize with framework defaults
    UgIndexModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    int getIndexRange();  // i=1 to 4
    double computeU_gi(int i);  // U_gi for i=1-4 (J/m^3)
    double computeK_i(int i);   // k_i for i
    double computeKUgi(int i);  // k_i * U_gi
    double computeSumKUgi(int i_min=1, int i_max=4);  // Sum for F_U

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print breakdown by i
    void printIndexBreakdown();
};

#endif // UG_INDEX_MODULE_H
