// PiConstantModule.h
// Modular C++ implementation of the Mathematical Constant Pi (π) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes π ≈3.14159 (unitless) and its role in oscillatory terms like cos(π t_n), sin(ω_c t), with ω_c=2π / period.
// Pluggable: #include "PiConstantModule.h"
// PiConstantModule mod; mod.computeCosPiTn(0.0); mod.updateVariable("t_n", new_value);
// Variables in std::map; examples for U_m μ_j and U_g1 cos(π t_n) at t=0, t_n=0.
// Approximations: π=3.141592653589793; sin(ω_c * 0)=0; cos(π * 0)=1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef PI_CONSTANT_MODULE_H
#define PI_CONSTANT_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class PiConstantModule {
private:
    std::map<std::string, double> variables;
    double computeCosPiTn(double t_n);
    double computeSinOmegaCT(double t);
    double computeMuJExample(double t);
    double computeUg1CosTerm(double t_n);

public:
    // Constructor: Initialize with framework defaults
    PiConstantModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computePi();  // ≈3.141592653589793 (unitless)
    double computeCosPiTn(double t_n);  // cos(π t_n)
    double computeSinOmegaCT(double t);  // sin(ω_c t), ω_c=2π / period
    double computeMuJExample(double t);  // Example μ_j with sin(ω_c t)
    double computeUg1CosTerm(double t_n);  // cos(π t_n) in U_g1

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // PI_CONSTANT_MODULE_H
