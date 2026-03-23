// FeedbackFactorModule.h
// Modular C++ implementation of the Feedback Factor (f_feedback) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes f_feedback=0.1 for ΔM_BH=1 dex (10x mass increase); scales (1 + f_feedback) in U_g4 term.
// Pluggable: #include "FeedbackFactorModule.h"
// FeedbackFactorModule mod; mod.computeU_g4(0.0, 0.0); mod.updateVariable("f_feedback", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; ΔM_BH=1 dex → M_bh_final=10*M_bh_initial.
// Approximations: cos(π t_n)=1; e^{-α t}=1 at t=0; α=0.001 day^-1 (scaled to s^-1 if needed).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef FEEDBACK_FACTOR_MODULE_H
#define FEEDBACK_FACTOR_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class FeedbackFactorModule {
private:
    std::map<std::string, double> variables;
    double computeDeltaM_BH();  // 1 dex = log10(10) = factor of 10
    double computeM_bh_final();
    double computeU_g4(double t, double t_n);
    double computeU_g4_no_feedback(double t, double t_n);

public:
    // Constructor: Initialize with framework defaults
    FeedbackFactorModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeF_feedback();  // f_feedback=0.1 (unitless)
    double computeDeltaM_BH();   // 1 dex
    double computeM_bh_final();  // 10 * M_bh_initial
    double computeU_g4(double t, double t_n);  // With feedback (J/m^3)
    double computeU_g4_no_feedback(double t, double t_n);  // Without (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print U_g4 comparison (with/without feedback)
    void printU_g4_comparison(double t = 0.0, double t_n = 0.0);
};

#endif // FEEDBACK_FACTOR_MODULE_H
