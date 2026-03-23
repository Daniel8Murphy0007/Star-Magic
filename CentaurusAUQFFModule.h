// CentaurusAUQFFModule.h
// Modular C++ implementation of the UQFF Force for NGC 5128 (Centaurus A, Radio Galaxy) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes F_U_Bi_i,enhanced as integral from x1 to x2 of [-F0 + DPM terms + LENR + activation + DE + EM + neutron + rel + Sweet + Kozima].
// Pluggable: #include "CentaurusAUQFFModule.h"
// CentaurusAUQFFModule mod; mod.computeF_U_Bi(0.0, 1.17e23, 0.0); mod.updateVariable("M", new_value);
// Variables in std::map; defaults for NGC 5128 (M=5.5e9 M_sun, r=1.17e23 m, level=13); ~ -8.32e217 N at t=0.
// Approximations: Integral approx via average * Δx; cos(θ)=1; ω_LENR / ω_0 tuned; Sweet/Kozima small/negligible.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef CENTAURUS_A_UQFF_MODULE_H
#define CENTAURUS_A_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class CentaurusAUQFFModule {
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
    // Constructor: Initialize with NGC 5128 defaults
    CentaurusAUQFFModule();

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

#endif // CENTAURUS_A_UQFF_MODULE_H
