// SolarCycleFrequencyModule.h
// Modular C++ implementation of the Solar Cycle Frequency (ω_c) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ω_c = 2π / 3.96e8 s⁻¹ (~1.59e-8 rad/s, period ~12.55 years); used in sin(ω_c t) for μ_j in U_m.
// Pluggable: #include "SolarCycleFrequencyModule.h"
// SolarCycleFrequencyModule mod; mod.computeMuJExample(0.0); mod.updateVariable("period", new_value);
// Variables in std::map; example for Sun at t=0 (sin=0, μ_j=3.38e23 T·m³); t~1 year: slight increase.
// Approximations: Period=3.96e8 s (~12.55 yr); base B_j=1e3 T.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SOLAR_CYCLE_FREQUENCY_MODULE_H
#define SOLAR_CYCLE_FREQUENCY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class SolarCycleFrequencyModule {
private:
    std::map<std::string, double> variables;
    double computeOmega_c();
    double computeSinOmegaCT(double t);

public:
    // Constructor: Initialize with framework defaults
    SolarCycleFrequencyModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeOmega_c();  // 2π / period s⁻¹
    double computeSinOmegaCT(double t);  // sin(ω_c t)
    double computeMuJExample(double t);  // (10^3 + 0.4 sin(ω_c t)) * 3.38e20 T·m³

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // SOLAR_CYCLE_FREQUENCY_MODULE_H
