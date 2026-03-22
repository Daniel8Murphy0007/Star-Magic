// UGC10214UQFFModule.h
// Modular C++ implementation of MUGE & UQFF Integration for UGC 10214 (Tadpole Galaxy).
// Models tidal interaction with companion VV29c, streaming tail, and dark matter.
// Params: M_total=1e11 Msun, r=55 kpc, SFR=4.67 Msun/yr, M_dwarf=3.5e9 Msun, v_tail=400 km/s.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UGC10214_UQFF_MODULE_H
#define UGC10214_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class UGC10214UQFFModule {
private:
    std::map<std::string, double> variables;
    double computeHtz(double z_val);
    double computeMmerge(double t);
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
    UGC10214UQFFModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    double computeG(double t, double r);
    std::string getEquationText();
    void printVariables();
};

#endif // UGC10214_UQFF_MODULE_H
