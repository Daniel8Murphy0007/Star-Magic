// NGC1300UQFFModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (MUGE & UQFF Integration)
// for Barred Spiral Galaxy NGC 1300 Evolution.
// Models NGC 1300's gravitational dynamics: bar-driven gas funneling, spiral arm density waves,
// star formation, dust lanes, and dark matter.
// Usage: #include "NGC1300UQFFModule.h"; NGC1300UQFFModule mod; mod.computeG(t, r);
// NGC 1300 params: M=1e11 Msun, r=11.79 kpc, SFR=1 Msun/yr, v_arm=200 km/s, B=1e-5 T, z=0.005.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef NGC1300_UQFF_MODULE_H
#define NGC1300_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class NGC1300UQFFModule {
private:
    std::map<std::string, double> variables;
    double computeHtz(double z_val);
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
    double computeMsfFactor(double t);
    double computeRt(double t);

public:
    NGC1300UQFFModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    double computeG(double t, double r);
    std::string getEquationText();
    void printVariables();
};

#endif // NGC1300_UQFF_MODULE_H
