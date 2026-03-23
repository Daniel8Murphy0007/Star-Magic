// Abell2256UQFFModule.h
// Modular C++ implementation of the UQFF Force for Abell 2256 (Galaxy Cluster, merging ICM) in the UQFF framework.
// Computes F_U_Bi_i,enhanced = ∫_{x1}^{x2} [-F0 + DPM_mom + DPM_grav + DPM_stab + LENR + Activation + DE + EM + Neutron + Rel + Sweet_vac + Kozima] dx
// Abell 2256: M=1.5e15 M_sun, r=1.42e25 m, level=13; F ≈ -1.23e218 N (largest cluster UQFF output)
// Cluster physics: merger shocks, radio relics, ICM plasma, Kozima coupling dominant
// Pluggable: #include "Abell2256UQFFModule.h"
// Abell2256UQFFModule mod; mod.computeF_U_Bi(0.0, 1.42e25, 0.0);
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025. Source: grok_share_b0a3dc1d.txt

#ifndef ABELL_2256_UQFF_MODULE_H
#define ABELL_2256_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class Abell2256UQFFModule {
private:
    std::map<std::string, double> variables;
    double computeDPM_momentum_term(double r);
    double computeDPM_gravity_term(double r);
    double computeDPM_stability_term();
    double computeLENR_term();
    double computeActivation_term(double t);
    double computeDE_term(double L_x);
    double computeEM_term();
    double computeNeutron_term();
    double computeRel_term(double E_cm_eff);
    double computeSweet_vac_term();
    double computeKozima_term();
    double computeIntegrand(double x, double t);
    double computeIntegral(double x1, double x2, double t, int n_points = 1000);

public:
    // Constructor: Initialize with Abell 2256 defaults
    Abell2256UQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: F_U_Bi_i,enhanced (N)
    double computeF_U_Bi(double x1, double x2, double t);

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // ABELL_2256_UQFF_MODULE_H
