// MUGEModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (MUGE) for multiple astronomical systems.
// This module integrates compressed UQFF (from documents) and resonance-based UQFF models.
// Can be plugged into a base program by including this header and linking the .cpp.
// Usage: #include "MUGEModule.h"
// MUGEModule mod("Magnetar"); mod.computeG(t); mod.updateVariable("M", new_value);
// Supports systems: Magnetar SGR 1745-2900, Sagittarius A*, Tapestry of Blazing Starbirth, Westerlund 2,
// Pillars of Creation, Rings of Relativity, Students Guide to the Universe.
// Variables in std::map for dynamic updates; both models computable via computeG_compressed(double t) and computeG_resonance(double t).
// All terms conserved: base gravity, expansion, superconductivity, Ug terms, Lambda, quantum integral, fluid, DM perturbations, system-specific (e.g., stellar wind, lensing).
// Approximations: Ug1/Ug2/Ug4 negligible in compressed; integral normalized; DM fraction variable; resonance freqs tuned per system.
// Associated text: getEquationText() for both models.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef MUGE_MODULE_H
#define MUGE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

enum class SystemType {
    MAGNETAR_SGR_1745_2900,
    SAGITTARIUS_A,
    TAPESTRY_BLAZING_STARBIRTH,
    WESTERLUND_2,
    PILLARS_CREATION,
    RINGS_RELATIVITY,
    STUDENTS_GUIDE_UNIVERSE
};

class MUGEModule {
private:
    std::map<std::string, double> variables;
    SystemType current_system;
    double computeH(double t, double z);
    double computeQuantumTerm();
    double computeFluidTerm(double g_base);
    double computeDMTerm();
    double computeUgSum();
    double computeLambdaTerm();
    double computeResonantTerm(double t);
    double computeEMTerm();
    double computeSystemSpecificTerm(double t);
    // Resonance-specific helpers
    double computeADPM();
    double computeATHz();
    double computeAvacDiff();
    double computeASuperFreq();
    double double computeAAetherRes();
    double computeUg4i();
    double computeAQuantumFreq();
    double double computeAAetherFreq();
    double computeAFluidFreq();
    double computeOscTerm(double t);
    double computeAExpFreq();
    double computeFTRZ();

public:
    // Constructor: Initialize with system-specific defaults
    MUGEModule(SystemType sys = SystemType::MAGNETAR_SGR_1745_2900);

    // Set system
    void setSystem(SystemType sys);

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations: Compressed and Resonance MUGE g(r, t)
    double computeG_compressed(double t);
    double computeG_resonance(double t);

    // Output descriptive texts
    std::string getEquationText_compressed();
    std::string getEquationText_resonance();

    // Print all current variables
    void printVariables();
};

#endif // MUGE_MODULE_H
