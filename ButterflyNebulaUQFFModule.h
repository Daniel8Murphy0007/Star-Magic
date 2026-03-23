// ButterflyNebulaUQFFModule.h
// Modular C++ implementation of the UQFF Force for Butterfly Nebula (NGC 6302) in the UQFF framework.
// Computes F_U_Bi_i,enhanced = ∫_{x1}^{x2} [-F0 + DPM_mom + DPM_grav + DPM_stab + LENR + Activation + DE + EM + Neutron + Rel + Sweet_vac + Kozima] dx
// NGC 6302: M=0.64 M_sun, r=3.22e19 m, level=13; F ≈ -2.09e212 N (repulsive stabilization)
// Kozima: k_n σ_n (ω_LENR/ω_0) ≈ 7.85e30 N; Sweet-Patterson: ρ_[UA] DPM_stab V ≈ 7.09e-39 N
// Pluggable: #include "ButterflyNebulaUQFFModule.h"
// ButterflyNebulaUQFFModule mod; mod.computeF_U_Bi(0.0, 3.22e19, 0.0);
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025. Source: grok_share_b0a3dc1d.txt

#ifndef BUTTERFLY_NEBULA_UQFF_MODULE_H
#define BUTTERFLY_NEBULA_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class ButterflyNebulaUQFFModule {
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
    // Constructor: Initialize with NGC 6302 defaults
    ButterflyNebulaUQFFModule();

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

#endif // BUTTERFLY_NEBULA_UQFF_MODULE_H
