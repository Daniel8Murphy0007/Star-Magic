// UaVacuumDensityModule.h
// Modular C++ implementation of the Vacuum Energy Density of [UA] (ρ_vac,[UA]) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ρ_vac,[UA] = 7.09e-36 J/m³ (Sun, level 13); scales in U_g2, U_i, T_s terms.
// Pluggable: #include "UaVacuumDensityModule.h"
// UaVacuumDensityModule mod; mod.computeU_g2_example(1.496e13); mod.updateVariable("rho_vac_UA", new_value);
// Variables in std::map; example for Sun at r=1.496e13 m; U_g2 ≈1.18e53 J/m³, U_i ≈1.38e-47 J/m³.
// Approximations: S(r - R_b)=1; (1 + δ_sw v_sw)=5001; λ_i=1.0; f_TRZ=0.1; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UA_VACUUM_DENSITY_MODULE_H
#define UA_VACUUM_DENSITY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class UaVacuumDensityModule {
private:
    std::map<std::string, double> variables;
    double computeU_g2_base(double r);
    double computeU_i_base(double t, double t_n);

public:
    // Constructor: Initialize with framework defaults (Sun, level 13)
    UaVacuumDensityModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeRho_vac_UA();  // 7.09e-36 J/m³
    double computeU_g2_example(double r);  // U_g2 with ρ_vac,[UA] (J/m³)
    double computeU_i_example(double t, double t_n);  // U_i with ρ_vac,[UA] (J/m³)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // UA_VACUUM_DENSITY_MODULE_H
