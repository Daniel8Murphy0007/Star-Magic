// ScmVelocityModule.h
// [SCm] superconducting medium propagation velocity: v_SCm = c / n_SCm.
// Refractive index n_SCm ≈ sqrt(1 + ρ_SCm/ρ_vac). Near-vacuum: n ≈ 1.0000001.
// v_SCm ≈ c - δv where δv tiny; models UQFF signal propagation inside [SCm] medium.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L9357

#ifndef SCM_VELOCITY_MODULE_H
#define SCM_VELOCITY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class ScmVelocityModule {
private:
    std::map<std::string, double> variables;

public:
    ScmVelocityModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeRefractiveIndex();         // n_SCm = sqrt(1 + ρ_SCm/ρ_vac)
    double computeSCmVelocity();             // v_SCm = c / n_SCm [m/s]
    double computeSlowdown();               // δv = c - v_SCm [m/s]
    double computePropagationDelay(double d); // Δt = d/v_SCm - d/c [s]
    std::string getEquationText();
    void printVariables();
};

#endif // SCM_VELOCITY_MODULE_H
