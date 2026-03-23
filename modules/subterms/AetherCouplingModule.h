// AetherCouplingModule.h
// Modular C++ implementation of the Aether Coupling Constant (η) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes the Aether metric perturbation A_μν = g_μν + η * T_s^{μν}, where η is the dimensionless Aether coupling constant.
// Pluggable: #include "AetherCouplingModule.h"
// AetherCouplingModule mod; mod.computePerturbation(); mod.updateVariable("eta", new_value);
// Variables in std::map; supports diagonal metric components [1, -1, -1, -1] for flat Minkowski.
// Includes stress-energy tensor T_s from ρ_vac_UA, ρ_vac_SCm, ρ_vac_A; computes perturbed metric.
// Approximations: Diagonal T_s ≈ 1.123e7 J/m^3; perturbation ~1e-15; weak coupling preserves flatness.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef AETHER_COUPLING_MODULE_H
#define AETHER_COUPLING_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class AetherCouplingModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> g_mu_nu;  // Background metric [1, -1, -1, -1]
    std::vector<double> computePerturbedMetric();

public:
    // Constructor: Initialize with framework defaults
    AetherCouplingModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeT_s();  // Stress-energy tensor scalar approx (J/m^3)
    double computePerturbation();  // η * T_s
    std::vector<double> computeA_mu_nu();  // Perturbed metric A_μν

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print perturbed metric
    void printPerturbedMetric();
};

#endif // AETHER_COUPLING_MODULE_H
