// StellarMassModule.h
// Modular C++ implementation of the Stellar/Planetary Mass (M_s) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes M_s=1.989e30 kg (1 M_sun for Sun); scales M_s / r^2 in Universal Gravity U_g1 and U_g2 terms.
// Pluggable: #include "StellarMassModule.h"
// StellarMassModule mod; mod.computeU_g2(1.496e13); mod.updateVariable("M_s", new_value);
// Variables in std::map; example for Sun at r=1.496e13 m; U_g2 ≈1.18e53 J/m³.
// Approximations: S(r - R_b)=1; (1 + δ_sw v_sw)=5001; H_SCm=1; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef STELLAR_MASS_MODULE_H
#define STELLAR_MASS_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class StellarMassModule {
private:
    std::map<std::string, double> variables;
    double computeM_sOverR2(double r);
    double computeU_g1(double r);
    double computeU_g2(double r);

public:
    // Constructor: Initialize with framework defaults (Sun)
    StellarMassModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeM_s();  // 1.989e30 kg
    double computeM_sInMsun();  // 1 M_sun
    double computeM_sOverR2(double r);  // M_s / r^2 (kg/m²)
    double computeU_g1(double r);  // U_g1 example (J/m^3)
    double computeU_g2(double r);  // U_g2 example (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // STELLAR_MASS_MODULE_H
