// StellarRotationModule.h
// Modular C++ implementation of the Stellar/Planetary Rotation Rate (ω_s) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ω_s=2.5e-6 rad/s (~29-day Sun period); scales ω_s(t) in U_g3 cos(ω_s t π) and U_i ω_s cos(π t_n).
// Pluggable: #include "StellarRotationModule.h"
// StellarRotationModule mod; mod.computeU_g3(0.0); mod.updateVariable("omega_s", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; U_g3 ≈1.8e49 J/m³, U_i ≈1.38e-47 J/m³.
// Approximations: cos(π t_n)=1; f_TRZ=0.1; λ_i=1.0; ρ_vac sum=7.80e-36 J/m³.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef STELLAR_ROTATION_MODULE_H
#define STELLAR_ROTATION_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class StellarRotationModule {
private:
    std::map<std::string, double> variables;
    double computeOmega_s_t(double t);  // ω_s(t), simplified constant
    double computeU_g3(double t);
    double computeU_i(double t, double t_n);

public:
    // Constructor: Initialize with framework defaults (Sun)
    StellarRotationModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeOmega_s();  // 2.5e-6 rad/s
    double computeOmega_s_t(double t);  // ω_s(t) (rad/s)
    double computePeriod_days();  // ~29 days
    double computeU_g3(double t);  // U_g3 example (J/m^3)
    double computeU_i(double t, double t_n);  // U_i example (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // STELLAR_ROTATION_MODULE_H
