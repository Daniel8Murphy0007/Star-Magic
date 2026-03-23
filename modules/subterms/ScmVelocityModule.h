// ScmVelocityModule.h
// Modular C++ implementation of the [SCm] Velocity (v_SCm) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes v_SCm = 1e8 m/s (~c/3); scales in E_react = ρ_vac,[SCm] v_SCm² / ρ_vac,A * exp(-κ t) for U_m, U_bi, etc.
// Pluggable: #include "ScmVelocityModule.h"
// ScmVelocityModule mod; mod.computeE_react(0.0); mod.updateVariable("v_sc m", new_value);
// Variables in std::map; example for Sun at t=0 (E_react=1e46 J); t=2000 days: scales down via exp.
// Approximations: κ=0.0005 day⁻¹; ρ_vac,[SCm]=7.09e-37 J/m³; ρ_vac,A=1e-23 J/m³; U_m base=2.28e65 J/m³.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SCM_VELOCITY_MODULE_H
#define SCM_VELOCITY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class ScmVelocityModule {
private:
    std::map<std::string, double> variables;
    double computeE_react_base();
    double computeE_react(double t_day);
    double computeUmExample(double t_day);

public:
    // Constructor: Initialize with framework defaults
    ScmVelocityModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeV_sc m();  // 1e8 m/s
    double computeE_react(double t_day);  // ρ_[SCm] v_SCm² / ρ_A * exp(-κ t)
    double computeUmExample(double t_day);  // Simplified U_m with E_react (J/m³)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print velocity effects
    void printVelocityEffects(double t_day = 2000.0);
};

#endif // SCM_VELOCITY_MODULE_H
