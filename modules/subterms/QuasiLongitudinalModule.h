// QuasiLongitudinalModule.h
// Modular C++ implementation of the Quasi-Longitudinal Wave Factor (f_quasi) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes f_quasi=0.01 (unitless) and its scaling (1 + f_quasi) in Universal Magnetism U_m term.
// Pluggable: #include "QuasiLongitudinalModule.h"
// QuasiLongitudinalModule mod; mod.computeUmContribution(0.0); mod.updateVariable("f_quasi", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; minor 1% increase in U_m.
// Approximations: 1 - e^{-γ t cos(π t_n)}=0 at t=0; φ_hat_j=1; P_SCm=1; f_Heaviside=0.01 (1 + 10^13 f=1e11+1).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef QUASI_LONGITUDINAL_MODULE_H
#define QUASI_LONGITUDINAL_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class QuasiLongitudinalModule {
private:
    std::map<std::string, double> variables;
    double computeQuasiFactor();
    double computeUmBase(int j, double t);
    double computeUmContribution(int j, double t);

public:
    // Constructor: Initialize with framework defaults
    QuasiLongitudinalModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeF_quasi();  // 0.01 (unitless)
    double computeQuasiFactor();  // 1 + f_quasi = 1.01
    double computeUmContribution(int j, double t);  // U_m single string (J/m^3)
    double computeUmWithNoQuasi(int j, double t);  // Without quasi

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print U_m comparison (with/without quasi)
    void printUmComparison(int j = 1, double t = 0.0);
};

#endif // QUASI_LONGITUDINAL_MODULE_H
