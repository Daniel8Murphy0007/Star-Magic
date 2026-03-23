// BuoyancyCouplingModule.h
// Computes Universal Buoyancy terms U_bi in UQFF for i=1-4 (Ug1-Ug4).
// U_bi = -β_i × U_gi × Ω_g × (M_bh / d_g) × E_react × (1 + ε_sw ρ_vac,sw) × U_UA × cos(π t_n)
// β_i = 0.6 uniform (unitless) for all i; 60% gravitational counterforce.
// At Sun params t_n=0: U_b1 ≈ -1.94e27 J/m³ (Ug1-dominated, opposes collapse).
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L2082

#ifndef BUOYANCY_COUPLING_MODULE_H
#define BUOYANCY_COUPLING_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class BuoyancyCouplingModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> computeAllU_biInternal();

public:
    BuoyancyCouplingModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // β_i = 0.6 for all i
    double computeBeta(int i);

    // U_bi for specific i (Ug1-4)
    double computeU_bi(int i);

    // All four U_bi values
    std::vector<double> computeAllU_bi();

    // Sum of β_i terms for F_U
    double computeF_U_contribution();

    std::string getEquationText();
    void printVariables();
    void printU_bi();
};

#endif // BUOYANCY_COUPLING_MODULE_H
