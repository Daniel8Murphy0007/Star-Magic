// M51UQFFModule.h
// Modular C++ implementation of MUGE & UQFF Integration for M51 Whirlpool Galaxy.
// Novel physics: Ug4 BH reaction energy decay, r(t) expansion, spiral density wave psi.
// Tidal companion NGC 5195 interactions via F_tidal and computeMmerge.
// Params: M_vis=1.2e11 Msun, M_DM=4e10 Msun, r=23.58 kpc, z=0.002, SFR=1 Msun/yr.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef M51_UQFF_MODULE_H
#define M51_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class M51UQFFModule {
private:
    std::map<std::string, double> variables;
    double computeHtz(double z_val);
    double computeMmerge(double t);     // NGC 5195 merger mass evolution
    double computeRt(double t);         // radial expansion r(t) = r0 + v_r*t
    double computeFenv(double t);
    double computeUg1(double t);        // magnetic dipole I*A*omega*B
    double computeUg2(double t);        // superconductor field B^2/(2*mu_0)
    double computeUg3prime(double t);   // tidal: G*M_NGC5195/d^2
    double computeUg4(double t);        // BH reaction: k_4*E_react*exp(-0.0005*t)
    double computeUi(double t);         // inertia with vacuum density ratio + F_RZ
    double computePsiSpiral(double r, double theta, double t);  // spiral wave amplitude
    double computePsiIntegral(double r, double t);
    double computeQuantumTerm(double t_Hubble_val, double r);
    double computeFluidTerm(double g_base);
    double computeDMTerm(double r);
    double computeUgSum(double r);
    double computeMsfFactor(double t);

public:
    M51UQFFModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    double computeG(double t, double r);
    std::string getEquationText();
    void printVariables();
};

#endif // M51_UQFF_MODULE_H
