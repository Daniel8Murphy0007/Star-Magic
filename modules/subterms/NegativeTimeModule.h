// NegativeTimeModule.h
// Modular C++ implementation of the Negative Time Factor (t_n) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes t_n = t - t_0 (s or days, allows t_n < 0); used in cos(π t_n) for oscillations and exp(-γ t cos(π t_n)) for growth/decay.
// Pluggable: #include "NegativeTimeModule.h"
// NegativeTimeModule mod; mod.computeCosPiTn(1000.0); mod.updateVariable("t_0", new_value);
// Variables in std::map; defaults t_0=0, t=0 (t_n=0); example for U_m term with t_n negative.
// Approximations: cos even function; γ=5e-5 day^-1; at t_n=-1, exp term negative (growth).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef NEGATIVE_TIME_MODULE_H
#define NEGATIVE_TIME_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class NegativeTimeModule {
private:
    std::map<std::string, double> variables;
    double computeCosPiTn(double t_n);
    double computeExpTerm(double gamma, double t, double t_n);

public:
    // Constructor: Initialize with framework defaults
    NegativeTimeModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeT_n(double t);  // t_n = t - t_0 (s/days)
    double computeCosPiTn(double t);  // cos(π t_n)
    double computeExpTerm(double gamma, double t);  // exp(-γ t cos(π t_n))
    double computeOneMinusExp(double gamma, double t);  // 1 - exp(-γ t cos(π t_n))
    double computeUmExample(double t, double mu_over_rj = 2.26e10);  // Simplified U_m contrib

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print t_n effects (positive/negative)
    void printTnEffects(double t, double gamma = 5e-5);
};

#endif // NEGATIVE_TIME_MODULE_H
