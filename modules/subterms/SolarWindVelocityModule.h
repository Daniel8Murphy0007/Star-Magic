// SolarWindVelocityModule.h
// Modular C++ implementation of the Solar Wind Velocity (v_sw) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes v_sw=5e5 m/s (500 km/s); scales (1 + δ_sw v_sw) in Universal Gravity U_g2 term.
// Pluggable: #include "SolarWindVelocityModule.h"
// SolarWindVelocityModule mod; mod.computeU_g2(1.496e13); mod.updateVariable("v_sw", new_value);
// Variables in std::map; example for Sun at r=R_b=1.496e13 m; amplification ~5001x.
// Approximations: S(r - R_b)=1; H_SCm=1; E_react=1e46; ρ_sum=7.80e-36 J/m³.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SOLAR_WIND_VELOCITY_MODULE_H
#define SOLAR_WIND_VELOCITY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class SolarWindVelocityModule {
private:
    std::map<std::string, double> variables;
    double computeModulationFactor();
    double computeU_g2(double r);

public:
    // Constructor: Initialize with framework defaults (Sun)
    SolarWindVelocityModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeV_sw();  // 5e5 m/s
    double computeV_swKmS();  // 500 km/s
    double computeModulationFactor();  // 1 + δ_sw v_sw
    double computeU_g2(double r);  // U_g2 with modulation (J/m^3)
    double computeU_g2_no_sw(double r);  // Without v_sw (set=0)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // SOLAR_WIND_VELOCITY_MODULE_H
