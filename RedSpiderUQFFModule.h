// RedSpiderUQFFModule.h
// Modular C++ implementation of UQFF frequency-domain gravity for Red Spider Nebula (NGC 6537).
// Novel physics: g = Sigma f_i * lambda_Planck / (2*pi) — all UQFF terms as frequency contributions.
// Params: r=7.1e15 m, v_exp=3e5 m/s, z=0.0015, t_age=1900 yr.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef REDSPIDER_UQFF_MODULE_H
#define REDSPIDER_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class RedSpiderUQFFModule {
private:
    std::map<std::string, double> variables;
    double computeHtz(double z_val);
    double computeFtotal();        // sum of all frequency contributions
    double computeFdpm();          // f_DPM: di-pseudo-monopole frequency
    double computeFreactivity();   // f_react: charge reactivity
    double computeFsuper();        // f_super: SuperFreq constant
    double computeFfluid();        // f_fluid: Navier-Stokes frequency
    double computeFquantum();      // f_quantum: quantum frequency
    double computeFAether();       // f_Aether: aether frequency
    double computeFTHz();          // f_THz contribution
    double computeFenv(double t);

public:
    RedSpiderUQFFModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    // Returns acceleration a = f_total * lambda_Planck / (2*pi) [m/s^2]
    double computeG(double t, double r);
    std::string getEquationText();
    void printVariables();
};

#endif // REDSPIDER_UQFF_MODULE_H
