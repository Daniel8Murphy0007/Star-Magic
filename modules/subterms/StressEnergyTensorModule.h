// StressEnergyTensorModule.h
// Modular C++ implementation of the Stress-Energy Tensor (T_s^{μν}) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes T_s^{μν} ≈1.123e7 J/m³ (diagonal scalar); perturbs A_μν = g_μν + η T_s^{μν} (~1.123e-15).
// Pluggable: #include "StressEnergyTensorModule.h"
// StressEnergyTensorModule mod; mod.computeA_mu_nu(); mod.updateVariable("rho_vac_A", new_value);
// Variables in std::map; diagonal [tt, xx, yy, zz]; example for Sun at t_n=0.
// Approximations: Diagonal T_s = T_s_base + ρ_vac_A; η=1e-22; g_μν=[1,-1,-1,-1].
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef STRESS_ENERGY_TENSOR_MODULE_H
#define STRESS_ENERGY_TENSOR_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class StressEnergyTensorModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> g_mu_nu;  // Background [1, -1, -1, -1]
    double computeT_s();  // Scalar approx J/m³
    std::vector<double> computeA_mu_nu();

public:
    // Constructor: Initialize with framework defaults
    StressEnergyTensorModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeT_s();  // 1.123e7 J/m³
    double computePerturbation();  // η * T_s ≈1.123e-15
    std::vector<double> computeA_mu_nu();  // Perturbed metric

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print tensor and metric
    void printTensorAndMetric();
};

#endif // STRESS_ENERGY_TENSOR_MODULE_H
