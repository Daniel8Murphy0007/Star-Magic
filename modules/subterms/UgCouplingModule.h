// UgCouplingModule.h
// Computes scaled Universal Gravity terms k_i × U_gi for i=1-4 (Ug1-Ug4) in UQFF.
// Coupling constants: k1=1.5, k2=1.2, k3=1.8, k4=1.0 (all unitless)
// Sum ∑ k_i U_gi contributes to F_U; k3 largest (string/rotation dominant).
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L2448

#ifndef UG_COUPLING_MODULE_H
#define UG_COUPLING_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class UgCouplingModule {
private:
    std::map<std::string, double> variables;

public:
    UgCouplingModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // k_i coupling constant for index i (1=1.5, 2=1.2, 3=1.8, 4=1.0)
    double computeK(int i);

    // k_i × U_gi for specific i
    double computeScaledUg(int i);

    // Full sum ∑ k_i U_gi for i=1-4
    double computeSumK_Ugi();

    std::string getEquationText();
    void printVariables();
    void printScaledTerms();
};

#endif // UG_COUPLING_MODULE_H
