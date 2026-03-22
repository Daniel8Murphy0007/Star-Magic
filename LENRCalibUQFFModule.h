// LENRCalibUQFFModule.h
// Modular C++ implementation of UQFF calibration for LENR K_n neutron production.
// Novel physics: eta(t,n) = k_eta * exp(-[SS_q]^n * 2^6 * exp(-pi - t/yr)) * Um / rho_vac_UA
//               Non-local exponential: exp(-[SS_q]^n * 2^6 * exp(-pi - t/yr))
//               delta_n = (2*pi)^(n/6); rho_vac_state(n,t) with quantum state coupling.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef LENRCALIB_UQFF_MODULE_H_ROOT
#define LENRCALIB_UQFF_MODULE_H_ROOT

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class LENRCalibUQFFModule {
private:
    std::map<std::string, double> variables;
    std::string m_scenario;

    // NOVEL: non-local exponential coupling
    double computeNonLocalExp(int n, double t);
    // delta_n = (2*pi)^(n/6)
    double computeDeltaN(int n);
    // rho_vac state: 1e-23 * (0.1)^n * exp(-[SS_q]^n * 2^6 * exp(-pi - t/yr))
    double computeRhoVacState(int n, double t);
    // Um for LENR context
    double computeUm();

public:
    LENRCalibUQFFModule();
    void setScenario(const std::string& scenario);  // "hydride", "wires", "corona"
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    // eta(t,n) = k_eta * nonLocalExp(n,t) * Um / rho_vac_UA
    double computeEta(double t, int n = 1);
    // rho_vac_UA':SCm(n,t) for state-dependent vacuum density
    double computeRhoVacUASCm(int n, double t);
    std::string getEquationText();
    void printVariables();
};

#endif // LENRCALIB_UQFF_MODULE_H_ROOT
