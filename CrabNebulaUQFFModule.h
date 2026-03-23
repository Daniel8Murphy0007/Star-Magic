
// CrabNebulaUQFFModule.h
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration) for Crab Nebula Supernova Remnant Evolution.
// This module can be plugged into a base program (e.g., 'crab_sim.cpp') by including this header and linking the .cpp.
// Usage in base: #include "CrabNebulaUQFFModule.h"
// CrabNebulaUQFFModule mod; mod.computeF(t); mod.updateVariable("M", {new_real, new_imag});
// All variables are stored in a std::map for dynamic addition/subtraction/update, using complex<double> for real/imaginary components.
// Nothing is negligible: Includes all terms - base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Integral approximated as integrand * x2 (quadratic root); imag parts small and not fully scaled; LENR dominant due to low ω_0; x2 from quadratic solver approx.
// Crab Nebula params: M=1e31 kg, r=4.73e16 m, L_X=1e27 W, B0=3e-8 T, t=3.06e10 s, ω_0=1e-12 s^-1, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 14, 2025.

#ifndef CRAB_NEBULA_UQFF_MODULE_H
#define CRAB_NEBULA_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

using cdouble = std::complex<double>;

class CrabNebulaUQFFModule {
private:
    std::map<std::string, cdouble> variables;
    cdouble computeIntegrand(double t);
    cdouble computeDPM_resonance();
    cdouble computeX2();
    cdouble computeQuadraticRoot(cdouble a, cdouble b, cdouble c);
    cdouble computeLENRTerm();
    double computeG(double t);
    cdouble computeQ_wave(double t);
    cdouble computeUb1();
    cdouble computeUi(double t);

public:
    // Constructor: Initialize all variables with Crab Nebula defaults
    CrabNebulaUQFFModule();

    // Dynamic variable operations (complex)
    void updateVariable(const std::string& name, cdouble value);
    void addToVariable(const std::string& name, cdouble delta);
    void subtractFromVariable(const std::string& name, cdouble delta);

    // Core computation: Full F_U_Bi_i(r, t) for Crab Nebula (approx integral)
    cdouble computeF(double t);

    // Sub-equations
    cdouble computeCompressed(double t);  // Integrand
    cdouble computeResonant();
    cdouble computeBuoyancy();
    cdouble computeSuperconductive(double t);
    double computeCompressedG(double t);  // g(r,t)

    // Output descriptive text of the equation
    std::string getEquationText();

    // Print all current variables (for debugging/updates)
    void printVariables();
};

#endif // CRAB_NEBULA_UQFF_MODULE_H

