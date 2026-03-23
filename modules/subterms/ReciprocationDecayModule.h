// ReciprocationDecayModule.h
// Reciprocation decay rate: γ = γ_0 × exp(-t / τ_rec) — models the decay of UQFF
// reciprocal field coupling over cosmic time. At τ_rec = 1 Gyr and t=10 Gyr: γ ≈ 4.5e-5 γ_0.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L5712

#ifndef RECIPROCATION_DECAY_MODULE_H
#define RECIPROCATION_DECAY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class ReciprocationDecayModule {
private:
    std::map<std::string, double> variables;

public:
    ReciprocationDecayModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeDecayRate(double t);       // γ = γ_0 exp(-t/τ_rec) [s⁻¹]
    double computeHalfLife();               // t_½ = τ_rec ln(2) [s]
    double computeResidual(double t0, double t1); // ∫ γ dt over [t0,t1]
    double computeUgRecModulation(double Ug, double t); // Ug × γ/γ_0
    std::string getEquationText();
    void printVariables();
};

#endif // RECIPROCATION_DECAY_MODULE_H
