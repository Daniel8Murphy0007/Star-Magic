// TimeReversalZoneModule.h
// Modular C++ implementation of the Time-Reversal Zone Factor (f_TRZ) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes f_TRZ=0.1 (unitless); scales (1 + f_TRZ) in Universal Inertia U_i term for TRZ enhancement.
// Pluggable: #include "TimeReversalZoneModule.h"
// TimeReversalZoneModule mod; mod.computeU_i(0.0, 0.0); mod.updateVariable("f_TRZ", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; U_i ≈1.38e-47 J/m³ (with, +10%); without: ≈1.25e-47 J/m³.
// Approximations: λ_i=1.0; cos(π t_n)=1; ω_s=2.5e-6 rad/s; ρ_sum=7.80e-36 J/m³.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef TIME_REVERSAL_ZONE_MODULE_H
#define TIME_REVERSAL_ZONE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class TimeReversalZoneModule {
private:
    std::map<std::string, double> variables;
    double computeTRZFactor();
    double computeU_i_base(double t, double t_n);
    double computeU_i(double t, double t_n);

public:
    // Constructor: Initialize with framework defaults
    TimeReversalZoneModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeF_TRZ();  // 0.1 (unitless)
    double computeTRZFactor();  // 1 + f_TRZ = 1.1
    double computeU_i(double t, double t_n);  // U_i with TRZ (J/m^3)
    double computeU_i_no_TRZ(double t, double t_n);  // Without TRZ

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print U_i comparison (with/without TRZ)
    void printUiComparison(double t = 0.0, double t_n = 0.0);
};

#endif // TIME_REVERSAL_ZONE_MODULE_H
