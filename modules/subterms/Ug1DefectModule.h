// Ug1DefectModule.h
// Modular C++ implementation of the Ug1 Defect Factor (δ_def) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes δ_def = 0.01 * sin(0.001 t) (unitless); scales (1 + δ_def) in Universal Gravity U_g1 term.
// Pluggable: #include "Ug1DefectModule.h"
// Ug1DefectModule mod; mod.computeU_g1(0.0, 1.496e11); mod.updateVariable("amplitude", new_value);
// Variables in std::map; example for Sun at t=0 (δ_def=0, U_g1≈4.51e31 J/m³); t=1570.8 days: +1%.
// Approximations: α=0.001 day⁻¹; cos(π t_n)=1 at t_n=0; μ_s=3.38e23 T·m³; ∇(M_s/r)≈M_s/r²=8.89e7 m/s².
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UG1_DEFECT_MODULE_H
#define UG1_DEFECT_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class Ug1DefectModule {
private:
    std::map<std::string, double> variables;
    double computeDelta_def(double t_day);
    double computeU_g1(double t_day, double r);

public:
    // Constructor: Initialize with framework defaults
    Ug1DefectModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeDelta_def(double t_day);  // 0.01 * sin(0.001 t)
    double computeU_g1(double t_day, double r);  // U_g1 with defect (J/m^3)
    double computePeriod_years();  // ~17.22 years

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // UG1_DEFECT_MODULE_H
