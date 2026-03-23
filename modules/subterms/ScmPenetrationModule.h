// ScmPenetrationModule.h
// Modular C++ implementation of the [SCm] Penetration Factor (P_SCm) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes P_SCm ≈1 (unitless for Sun, ~1e-3 for planets); scales P_SCm in Universal Magnetism U_m term.
// Pluggable: #include "ScmPenetrationModule.h"
// ScmPenetrationModule mod; mod.computeUmContribution(0.0); mod.updateVariable("P_SCm", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; full penetration for plasma cores.
// Approximations: 1 - e^{-γ t cos(π t_n)}=0 at t=0; φ_hat_j=1; μ_j / r_j=2.26e10 T m².
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SCM_PENETRATION_MODULE_H
#define SCM_PENETRATION_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class ScmPenetrationModule {
private:
    std::map<std::string, double> variables;
    double computeUmBase(double t);
    double computeUmContribution(double t);

public:
    // Constructor: Initialize with framework defaults (Sun)
    ScmPenetrationModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeP_SCm();  // ≈1 for Sun (unitless)
    double computeUmContribution(double t);  // U_m with P_SCm (J/m^3)
    double computeUmPlanet(double t);  // For planet P_SCm=1e-3

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // SCM_PENETRATION_MODULE_H
