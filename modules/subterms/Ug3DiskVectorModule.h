// Ug3DiskVectorModule.h
// Modular C++ implementation of the Unit Vector in the Ug3 Disk Plane (φ̂_j) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes φ̂_j (unit vector, magnitude=1; e.g., [cos θ_j, sin θ_j, 0]); scales in Universal Magnetism U_m term.
// Pluggable: #include "Ug3DiskVectorModule.h"
// Ug3DiskVectorModule mod; mod.computeUmContribution(0.0, 1); mod.updateVariable("theta_j", new_value);
// Variables in std::map; example for j=1 at t=0, θ_j=0 (φ̂_j=[1,0,0], U_m≈2.28e65 J/m³).
// Approximations: φ̂_j magnitude=1; 1 - exp=0 at t=0; μ_j / r_j=2.26e10 T m².
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UG3_DISK_VECTOR_MODULE_H
#define UG3_DISK_VECTOR_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class Ug3DiskVectorModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> computePhiHat_j(int j);
    double computeUmBase(double t);
    double computeUmContribution(double t, int j);

public:
    // Constructor: Initialize with framework defaults
    Ug3DiskVectorModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    std::vector<double> computePhiHat_j(int j);  // Unit vector [cos θ_j, sin θ_j, 0]
    double computePhiHatMagnitude(int j);  // 1.0 (normalized)
    double computeUmContribution(double t, int j);  // U_m single string (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print φ̂_j and U_m
    void printVectorAndUm(int j = 1, double t = 0.0);
};

#endif // UG3_DISK_VECTOR_MODULE_H
