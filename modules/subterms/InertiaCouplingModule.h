// InertiaCouplingModule.h
// Computes inertia coupling I_c = M r² Ω² / c² linking rotational inertia to relativistic
// correction. Modulates Ug3 string-rotation terms in UQFF. At Sun: I_c ≈ 7.8e-10 (unitless).
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L4246

#ifndef INERTIA_COUPLING_MODULE_H
#define INERTIA_COUPLING_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class InertiaCouplingModule {
private:
    std::map<std::string, double> variables;

public:
    InertiaCouplingModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeInertiaCoupling();             // I_c = M r² Ω² / c² [unitless]
    double computeRelativisticCorrection();      // γ_I = 1 + I_c [unitless]
    double computeUg3InertiaTerm();             // Ug3 × (1 + I_c)
    std::string getEquationText();
    void printVariables();
};

#endif // INERTIA_COUPLING_MODULE_H
