// NegativeTimeModule.h
// Models time-reversal physics: g(t<0) = g(0) exp(iωt) — oscillatory complex gravity in
// pre-causal zones. f_TRZ correction in Andromeda UQFF. Provides TRZ boundary radius term.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L4800

#ifndef NEGATIVE_TIME_MODULE_H
#define NEGATIVE_TIME_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class NegativeTimeModule {
private:
    std::map<std::string, double> variables;

public:
    NegativeTimeModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeTRZFraction();         // f_TRZ = r_TRZ / r [unitless, ~0.1]
    double computeOscillation(double t); // Re[exp(iω t)] = cos(ω t) phase factor
    double computeGravityNegT(double g0, double t); // g(t<0) = g0 cos(ω t) [m/s²]
    double computeTRZRadius();           // r_TRZ = r/f_TRZ [m]
    std::string getEquationText();
    void printVariables();
};

#endif // NEGATIVE_TIME_MODULE_H
