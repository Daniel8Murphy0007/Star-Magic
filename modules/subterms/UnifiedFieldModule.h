// UnifiedFieldModule.h
// Modular C++ implementation of the Unified Field Strength (F_U) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes F_U as normalized vacuum energy density (J/m³) from Ug, Um, Ub, Ui, and Aether terms across 26 quantum levels.
// Pluggable: #include "UnifiedFieldModule.h"
// UnifiedFieldModule mod; mod.computeFU(double t); mod.updateVariable("U_g1", new_value);
// Variables in std::map; defaults for Sun at t=0 (level 13); normalization via coupling constants.
// Approximations: Dominant Um ~2.28e65 J/m³; Aether small; cos(π t_n)=1 at t_n=0.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UNIFIED_FIELD_MODULE_H
#define UNIFIED_FIELD_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class UnifiedFieldModule {
private:
    std::map<std::string, double> variables;
    double computeUgSum();
    double computeUm();
    double computeUbSum();
    double computeUi();
    double computeAether();

public:
    // Constructor: Initialize with framework defaults (Sun at t=0)
    UnifiedFieldModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: F_U(t) in J/m³
    double computeFU(double t);

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print component breakdown
    void printComponentBreakdown(double t);
};

#endif // UNIFIED_FIELD_MODULE_H
