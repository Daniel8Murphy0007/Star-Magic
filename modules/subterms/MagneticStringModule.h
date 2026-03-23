// MagneticStringModule.h
// Models magnetic string resonance in UQFF: T_s = (μ_0 I²)/(4π) ln(L/a)
// String tension T_s ≈ 1e12 N; I = string current (A); L = string length; a = core radius.
// Magnetic string topology links DPM spheres, mediating field transmission across quantum levels.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L2661

#ifndef MAGNETIC_STRING_MODULE_H
#define MAGNETIC_STRING_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class MagneticStringModule {
private:
    std::map<std::string, double> variables;

public:
    MagneticStringModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeStringTension();           // T_s = (μ_0 I²)/(4π) ln(L/a) [N]
    double computeStringEnergy(double l);    // E = T_s × l [J]
    double computeResonanceFreq();           // f_string = (1/2L) sqrt(T_s / μ_linear [Hz]
    std::string getEquationText();
    void printVariables();
};

#endif // MAGNETIC_STRING_MODULE_H
