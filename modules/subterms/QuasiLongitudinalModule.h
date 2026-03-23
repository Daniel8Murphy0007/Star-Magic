// QuasiLongitudinalModule.h
// Quasi-longitudinal EM wave energy density: E_QL = ε_0 E² / 2 — longitudinal component
// of vacuum EM field perpendicular to UQFF propagation axis. Modulates Ug2 reactivity term.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L5330

#ifndef QUASI_LONGITUDINAL_MODULE_H
#define QUASI_LONGITUDINAL_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class QuasiLongitudinalModule {
private:
    std::map<std::string, double> variables;

public:
    QuasiLongitudinalModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeEnergyDensity();           // E_QL = ε_0 E² / 2 [J/m³]
    double computeLongitudinalAmplitude();   // E_L = E × cos(θ) [V/m]
    double computeUg2Modulation();           // Ug2 × (1 + E_QL/E_vac) [m/s²]
    double computeWaveVector(double omega);  // k = ω/c [rad/m]
    std::string getEquationText();
    void printVariables();
};

#endif // QUASI_LONGITUDINAL_MODULE_H
