// MagneticMomentModule.h
// Modular C++ implementation of the Magnetic Moment of the j-th String (μ_j) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes μ_j = (10^3 + 0.4 sin(ω_c t)) * 3.38e20 T·m^3; scales μ_j / r_j in Universal Magnetism U_m and Ug3.
// Pluggable: #include "MagneticMomentModule.h"
// MagneticMomentModule mod; mod.computeMu_j(0.0); mod.updateVariable("base_mu", new_value);
// Variables in std::map; j-indexed; example for j=1 at t=0.
// Approximations: ω_c=2.5e-6 rad/s; at t=0, sin=0, μ_j≈3.38e23 T·m^3 (adjusted for example consistency).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef MAGNETIC_MOMENT_MODULE_H
#define MAGNETIC_MOMENT_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class MagneticMomentModule {
private:
    std::map<std::string, double> variables;
    double computeMu_j(int j, double t);
    double computeUmContrib(int j, double t);

public:
    // Constructor: Initialize with framework defaults
    MagneticMomentModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeMu_j(int j, double t);  // T·m^3
    double computeB_j(double t);  // Base field 10^3 + 0.4 sin(ω_c t) T
    double computeUmContrib(int j, double t);  // Example U_m single string (J/m^3)
    double computeUg3Contrib(double t);  // Example Ug3 (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print μ_j and contributions
    void printMomentContributions(int j = 1, double t = 0.0);
};

#endif // MAGNETIC_MOMENT_MODULE_H
