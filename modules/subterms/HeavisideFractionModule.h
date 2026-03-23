// HeavisideFractionModule.h
// Models the Heaviside fraction H(f) = θ(f - f_c) step function: 0 below threshold, 1 above.
// Used as phase selector in UQFF quantum state transitions and vacuum energy plateaus.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L3669

#ifndef HEAVISIDE_FRACTION_MODULE_H
#define HEAVISIDE_FRACTION_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class HeavisideFractionModule {
private:
    std::map<std::string, double> variables;

public:
    HeavisideFractionModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeHeaviside(double f);          // H(f): 0 if f < f_c, 1 if f >= f_c
    double computeSmoothHeaviside(double f);    // Sigmoid approximation: 1/(1+exp(-k(f-f_c)))
    double computeFraction(double f1, double f2);  // Fraction of [f1,f2] above f_c
    std::string getEquationText();
    void printVariables();
};

#endif // HEAVISIDE_FRACTION_MODULE_H
