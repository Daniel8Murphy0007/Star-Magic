// ScmReactivityDecayModule.h
// Modular C++ implementation of the [SCm] Reactivity Decay Rate (κ) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes κ=0.0005 day⁻¹ (~5.8e-6 s⁻¹); used in E_react = 10^46 * exp(-κ t) for decay in U_m, U_bi, etc.
// Pluggable: #include "ScmReactivityDecayModule.h"
// ScmReactivityDecayModule mod; mod.computeE_react(0.0); mod.updateVariable("kappa_day", new_value);
// Variables in std::map; example for Sun at t=0 (E_react=1e46); t=2000 days: ~3.68e45.
// Approximations: t in days; timescale ~5.5 years; integrates into U_m example.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SCM_REACTIVITY_DECAY_MODULE_H
#define SCM_REACTIVITY_DECAY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class ScmReactivityDecayModule {
private:
    std::map<std::string, double> variables;
    double computeKappa_s();  // κ in s⁻¹
    double computeE_react(double t_day);
    double computeUmExample(double t_day);

public:
    // Constructor: Initialize with framework defaults
    ScmReactivityDecayModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeKappa_day();  // 0.0005 day⁻¹
    double computeKappa_s();    // ~5.8e-6 s⁻¹
    double computeE_react(double t_day);  // 1e46 * exp(-κ t)
    double computeUmExample(double t_day);  // Simplified U_m with E_react (J/m³)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print decay effects
    void printDecayEffects(double t_day = 2000.0);
};

#endif // SCM_REACTIVITY_DECAY_MODULE_H
