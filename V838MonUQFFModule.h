// V838MonUQFFModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (MUGE & UQFF Integration)
// for V838 Monocerotis Light Echo Evolution.
// Models light echo intensity, outburst luminosity, dust scattering, gravitational modulation via Ug1,
// time-reversal (f_TRZ), and Aether ([UA]) effects.
// Usage: #include "V838MonUQFFModule.h"; V838MonUQFFModule mod; mod.computeIecho(t, r);
// V838 Mon params: M_s=8 Msun, L_outburst=2.3e38 W, rho_0=1e-22 kg/m^3, d=6.1 kpc, B=1e-5 T.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef V838MON_UQFF_MODULE_H
#define V838MON_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class V838MonUQFFModule {
private:
    std::map<std::string, double> variables;
    double computeUg1(double t, double r);
    double computeRhodust(double r, double t);
    double computeIechoBase(double r);
    double computeTRZCorrection();
    double computeUAscCorrection();

public:
    V838MonUQFFModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    double computeIecho(double t, double r);
    std::string getEquationText();
    void printVariables();
};

#endif // V838MON_UQFF_MODULE_H
