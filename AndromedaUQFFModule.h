// AndromedaUQFFModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (UQFF) for Andromeda Galaxy Evolution.
// This module can be plugged into a base program by including this header and linking the .cpp.
// Usage: #include "AndromedaUQFFModule.h"
// AndromedaUQFFModule mod; mod.computeG(t); mod.updateVariable("M", new_value);
// Variables stored in std::map for dynamic updates.
// Includes base gravity with expansion and TRZ, BH term, dust friction a_dust, EM/Aether term.
// Approximations: z=-0.001 (blueshift); dust scaled by 1e-12; EM normalized to proton mass.
// Andromeda params: M=1e12 Msun, r=1.04e21 m, M_BH=1.4e8 Msun, v_orbit=2.5e5 m/s, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef ANDROMEDA_UQFF_MODULE_H
#define ANDROMEDA_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class AndromedaUQFFModule {
private:
    std::map<std::string, double> variables;
    double computeHz();
    double computeADust();
    double computeEMBase();
    double computeEMTerm();

public:
    // Constructor: Initialize with Andromeda defaults
    AndromedaUQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: g_Andromeda(r, t)
    double computeG(double t);

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print evolution table (0-10 Gyr, 2 Gyr steps)
    void printEvolutionTable();
};

#endif // ANDROMEDA_UQFF_MODULE_H
