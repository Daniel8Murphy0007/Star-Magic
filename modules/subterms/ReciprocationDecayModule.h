// ReciprocationDecayModule.h
// Modular C++ implementation of the Reciprocation Decay Rate (γ) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes γ=0.00005 day⁻¹ (~5.8e-10 s⁻¹); used in exp(-γ t cos(π t_n)) for U_m decay.
// Pluggable: #include "ReciprocationDecayModule.h"
// ReciprocationDecayModule mod; mod.computeOneMinusExp(1000.0, 0.0); mod.updateVariable("gamma_day", new_value);
// Variables in std::map; example for t=1000 days, t_n=0; 1-exp ≈0.049.
// Approximations: cos(π t_n)=1; timescale ~55 years; μ_j / r_j=2.26e10 T m².
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef RECIPROCATION_DECAY_MODULE_H
#define RECIPROCATION_DECAY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class ReciprocationDecayModule {
private:
    std::map<std::string, double> variables;
    double computeGamma_s();  // γ in s⁻¹
    double computeCosPiTn(double t_n);
    double computeExpTerm(double t_day, double t_n);
    double computeOneMinusExp(double t_day, double t_n);

public:
    // Constructor: Initialize with framework defaults
    ReciprocationDecayModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeGamma_day();  // 0.00005 day⁻¹
    double computeGamma_s();    // ~5.8e-10 s⁻¹
    double computeCosPiTn(double t_n);  // cos(π t_n)
    double computeExpTerm(double t_day, double t_n);  // exp(-γ t cos(π t_n))
    double computeOneMinusExp(double t_day, double t_n);  // 1 - exp(...)
    double computeUmExample(double t_day, double t_n, double mu_over_rj = 2.26e10);  // Simplified U_m contrib

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print decay effects
    void printDecayEffects(double t_day = 1000.0, double t_n = 0.0);
};

#endif // RECIPROCATION_DECAY_MODULE_H
