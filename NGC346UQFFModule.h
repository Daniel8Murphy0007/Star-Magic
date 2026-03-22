// NGC346UQFFModule.h
// Modular C++ implementation of UQFF Integration for NGC 346 (SMC Star-Forming Nebula).
// Novel physics: Um = q*v_rad*B (universal magnetism from radial velocity blueshift),
//               Ug3 magnetic strings = G*M/r^2 * (rho_gas/rho_vac_UA).
// Params: M=1000 Msun, r=5 pc, SFR=0.1 Msun/yr, v_rad=-10 km/s (blueshift), z=0.0006.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef NGC346_UQFF_MODULE_H
#define NGC346_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class NGC346UQFFModule {
private:
    std::map<std::string, double> variables;
    double computeHtz(double z_val);
    double computeFenv(double t);
    double computeUm(double t);            // NOVEL: Um = q*v_rad*B
    double computeUg1(double t);
    double computeUg2(double t);
    double computeUg3strings(double r);   // NOVEL: Ug3 magnetic strings disk
    double computeUg4(double t);
    double computeUi(double t);            // NOVEL: inertia with vacuum density ratio
    double computeEcore(double r);         // core energy = Ug3 + Ui*rho
    double computeTempCore(double r);      // core temperature from Ug3
    double computePsiIntegral(double r, double t);
    double computeQuantumTerm(double t_Hubble_val, double r);
    double computeFluidTerm(double g_base);
    double computeDMTerm(double r);
    double computeUgSum(double r);
    double computeMsfFactor(double t);
    double computeRt(double t);

public:
    NGC346UQFFModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    double computeG(double t, double r);
    std::string getEquationText();
    void printVariables();
};

#endif // NGC346_UQFF_MODULE_H
