// SurfaceTemperatureModule.h
// Modular C++ implementation of the Surface Temperature (T_s) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes T_s=5778 K (Sun effective); potential scaling T_s / T_s_ref in B_j for U_g3 magnetic strings.
// Pluggable: #include "SurfaceTemperatureModule.h"
// SurfaceTemperatureModule mod; mod.computeU_g3_example(0.0, 5778.0); mod.updateVariable("T_s", new_value);
// Variables in std::map; example for Sun at t=0; T_s=5778 K → U_g3≈1.8e49 J/m³ (full); T_s=10000 K: ~3.11e49 J/m³.
// Approximations: T_s_ref=5778 K (Sun); cos(ω_s t π)=1; P_core=1; E_react=1e46; hypothetical B_j scaling.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SURFACE_TEMPERATURE_MODULE_H
#define SURFACE_TEMPERATURE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class SurfaceTemperatureModule {
private:
    std::map<std::string, double> variables;
    double computeB_j_hypothetical(double t, double T_s);
    double computeU_g3_example(double t, double T_s);

public:
    // Constructor: Initialize with framework defaults (Sun)
    SurfaceTemperatureModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeT_s();  // 5778 K (Sun)
    double computeB_j_hypothetical(double t, double T_s);  // Scaled B_j (T)
    double computeU_g3_example(double t, double T_s);  // U_g3 with scaling (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // SURFACE_TEMPERATURE_MODULE_H
