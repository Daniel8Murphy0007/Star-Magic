// MagneticStringModule.h
// Modular C++ implementation of the Distance Along Magnetic String's Path (r_j) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes r_j = 1.496e13 m (100 AU) and its conversions; scales μ_j / r_j in Universal Magnetism U_m and Ug3.
// Pluggable: #include "MagneticStringModule.h"
// MagneticStringModule mod; mod.computeMuOverRj(); mod.updateVariable("r_j", new_value);
// Variables in std::map; j-indexed strings; example for j=1 at t=0.
// Approximations: γ=5e-5 day^-1; cos(π t_n)=1; φ_hat_j=1; at t=0, 1 - exp term=0.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef MAGNETIC_STRING_MODULE_H
#define MAGNETIC_STRING_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class MagneticStringModule {
private:
    std::map<std::string, double> variables;
    double computeRjInAU();
    double computeRjInLy();
    double computeRjInPc();
    double computeMuOverRj(int j);
    double computeUmContribution(int j);

public:
    // Constructor: Initialize with framework defaults
    MagneticStringModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeRj(int j);  // r_j in m (default 1.496e13)
    double computeRjInAU(int j);
    double computeRjInLy(int j);
    double computeRjInPc(int j);
    double computeMu_j(int j, double t);  // Magnetic moment
    double computeMuOverRj(int j);
    double computeUmContribution(int j, double t);  // Single string to U_m
    double computeUg3Contribution();  // Example Ug3 influence

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print r_j conversions and contributions
    void printStringContributions(int j = 1, double t = 0.0);
};

#endif // MAGNETIC_STRING_MODULE_H
