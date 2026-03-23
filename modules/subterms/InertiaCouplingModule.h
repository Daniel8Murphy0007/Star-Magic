// InertiaCouplingModule.h
// Modular C++ implementation of the Inertia Coupling Constants (λ_i) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes λ_i=1.0 (unitless, uniform for i=1-4) and scales U_i in F_U: -∑_i [λ_i U_i E_react].
// Pluggable: #include "InertiaCouplingModule.h"
// InertiaCouplingModule mod; mod.computeSumInertiaTerms(0.0); mod.updateVariable("rho_vac_SCm", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; U_i ≈1.38e-47 J/m³, contrib ≈ -0.138 J/m³.
// Approximations: Uniform λ_i=1.0; cos(π t_n)=1; f_TRZ=0.1; ω_s=2.5e-6 rad/s; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef INERTIA_COUPLING_MODULE_H
#define INERTIA_COUPLING_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class InertiaCouplingModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> computeAllInertiaTerms(double t);

public:
    // Constructor: Initialize with framework defaults
    InertiaCouplingModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeLambda_i(int i);  // λ_i=1.0 (unitless)
    double computeU_i(int i, double t);  // U_i for i=1-4 (J/m^3)
    double computeInertiaTerm(int i, double t);  // -λ_i U_i E_react
    double computeSumInertiaTerms(double t);  // Sum for F_U contribution (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print inertia breakdown
    void printInertiaBreakdown(double t = 0.0);
};

#endif // INERTIA_COUPLING_MODULE_H
