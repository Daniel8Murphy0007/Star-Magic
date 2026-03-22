// NGC4676UQFFModule.h
// Modular C++ implementation of MUGE & UQFF Integration for NGC 4676 (The Mice) Colliding Galaxies.
// Novel physics: THz Aether-modulated H_eff(z), time-growing Ug2_THz term.
// Params: M_A=M_B=5e10 Msun, r=50 kpc (effective), SFR=5 Msun/yr, d=10 kpc, v_rel=400 km/s, z=0.022.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef NGC4676_UQFF_MODULE_H
#define NGC4676_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class NGC4676UQFFModule {
private:
    std::map<std::string, double> variables;
    double computeHtz(double z_val);
    double computeHeffz(double z_val);      // NEW: Aether-modulated H_eff(z)
    double computeMmerge(double t);
    double computeFenv(double t);
    double computeUg1(double t);
    double computeUg2(double t);
    double computeUg2THz(double t);         // NEW: time-growing THz Ug2 term
    double computeUg3prime(double t);
    double computeUg4(double t);
    double computeUi(double t);
    double computePsiIntegral(double r, double t);
    double computeQuantumTerm(double t_Hubble_val, double r);
    double computeFluidTerm(double g_base);
    double computeDMTerm(double r);
    double computeUgSum(double r);
    double computeMsfFactor(double t);
    double computeRt(double t);

public:
    NGC4676UQFFModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    double computeG(double t, double r);
    std::string getEquationText();
    void printVariables();
};

#endif // NGC4676_UQFF_MODULE_H
