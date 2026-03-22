// NGC2264UQFFModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (MUGE & UQFF Integration)
// for Cone Nebula (NGC 2264) Evolution.
// Models stellar winds, pillar erosion, protostar formation, dust/gas densities, and dark matter.
// NGC 2264 params: M=100 Msun, r=3.31e16 m, SFR=0.01 Msun/yr, v_wind=20 km/s, rho=1e-20 kg/m^3, z=0.0008.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef NGC2264_UQFF_MODULE_H
#define NGC2264_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class NGC2264UQFFModule {
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
    NGC2264UQFFModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    double computeG(double t, double r);
    std::string getEquationText();
    void printVariables();
};

#endif // NGC2264_UQFF_MODULE_H
