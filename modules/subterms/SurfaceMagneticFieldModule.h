// SurfaceMagneticFieldModule.h
// Modular C++ implementation of the Surface Magnetic Field (B_s) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes B_s range [1e-4, 0.4] T for Sun; influences B_j in U_g3 magnetic strings (scaled by B_s / B_ref).
// Pluggable: #include "SurfaceMagneticFieldModule.h"
// SurfaceMagneticFieldModule mod; mod.computeU_g3_example(0.0); mod.updateVariable("B_s_min", new_value);
// Variables in std::map; example for Sun at t=0; quiet Sun B_s=1e-4 T → U_g3≈4.5e45 J/m³.
// Approximations: B_ref=0.4 T (max sunspot); cos(ω_s t π)=1; P_core=1; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SURFACE_MAGNETIC_FIELD_MODULE_H
#define SURFACE_MAGNETIC_FIELD_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class SurfaceMagneticFieldModule {
private:
    std::map<std::string, double> variables;
    double computeB_j(double t, double B_s);
    double computeU_g3_example(double t, double B_s);

public:
    // Constructor: Initialize with framework defaults (Sun)
    SurfaceMagneticFieldModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeB_s_min();  // 1e-4 T (quiet Sun)
    double computeB_s_max();  // 0.4 T (sunspot max)
    double computeB_j(double t, double B_s);  // Scaled B_j (T)
    double computeU_g3_example(double t, double B_s);  // U_g3 with B_j (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // SURFACE_MAGNETIC_FIELD_MODULE_H
