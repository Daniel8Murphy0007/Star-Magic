// GalacticBlackHoleModule.h
// Modular C++ implementation of the Mass of the Galactic Black Hole (M_bh) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes M_bh=8.15e36 kg ≈4.1e6 M_sun; scales M_bh / d_g in Universal Buoyancy U_bi and Ug4.
// Pluggable: #include "GalacticBlackHoleModule.h"
// GalacticBlackHoleModule mod; mod.computeU_b1(); mod.updateVariable("M_bh", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0.
// Approximations: cos(π t_n)=1; (1 + ε_sw ρ_vac,sw)≈1; α=0.001 s^-1; f_feedback=0.1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef GALACTIC_BLACK_HOLE_MODULE_H
#define GALACTIC_BLACK_HOLE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class GalacticBlackHoleModule {
private:
    std::map<std::string, double> variables;
    double computeM_bhInMsun();
    double computeMbhOverDg();
    double computeU_b1();
    double computeU_g4();

public:
    // Constructor: Initialize with framework defaults
    GalacticBlackHoleModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeM_bh();  // 8.15e36 kg
    double computeM_bhInMsun();  // ≈4.1e6 M_sun
    double computeMbhOverDg();  // M_bh / d_g (kg/m)
    double computeU_b1();  // Universal Buoyancy example (J/m^3)
    double computeU_g4();  // Ug4 example (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // GALACTIC_BLACK_HOLE_MODULE_H
