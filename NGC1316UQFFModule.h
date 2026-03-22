// NGC1316UQFFModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (MUGE & UQFF Integration)
// for NGC 1316 "Hubble Spies Cosmic Dust Bunnies" MUGE Evolution.
// This module models NGC 1316's gravitational dynamics, incorporating merger-induced mass growth,
// AGN jet feedback, dust lane density modulation, and star cluster disruption.
// Usage: #include "NGC1316UQFFModule.h"; NGC1316UQFFModule mod; mod.computeG(t, r);
// Variables in std::map for dynamic updates; M(t) growths via M_merge(t) = 1e10 Msun exp(-t/tau).
// Approximations: psi_integral normalized; H(t,z) standard Omega_m=0.3, Omega_Lambda=0.7.
// NGC 1316 params: M=5e11 Msun, r=46 kpc, M_BH=1e8 Msun, rho_dust=1e-21 kg/m^3, z=0.005.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef NGC1316_UQFF_MODULE_H
#define NGC1316_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class NGC1316UQFFModule {
private:
    std::map<std::string, double> variables;
    double computeHtz(double z_val);
    double computeMmerge(double t);
    double computeRt(double t);
    double computeFenv(double t);
    double computeUg1(double t);
    double computeUg2(double t);
    double computeUg3prime(double t);
    double computeUg4(double t);
    double computeUi(double t);
    double computePsiIntegral(double r, double t);
    double computeQuantumTerm(double t_Hubble_val, double r);
    double computeFluidTerm(double g_base);
    double computeDMTerm(double r);
    double computeUgSum(double r);

public:
    // Constructor: Initialize with NGC 1316 defaults
    NGC1316UQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: g_NGC1316(r, t)
    double computeG(double t, double r);

    // Output descriptive text of the equation
    std::string getEquationText();

    // Print all current variables (for debugging)
    void printVariables();
};

#endif // NGC1316_UQFF_MODULE_H
