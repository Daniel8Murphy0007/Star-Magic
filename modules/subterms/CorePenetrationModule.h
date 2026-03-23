// CorePenetrationModule.h
// Computes UQFF field penetration depth into compact stellar cores:
// δ = (ρ_core / ρ_avg)^n × r_core — exponential density-driven field attenuation.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L5171

#ifndef CORE_PENETRATION_MODULE_H
#define CORE_PENETRATION_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class CorePenetrationModule {
private:
    std::map<std::string, double> variables;

public:
    CorePenetrationModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computePenetrationDepth();        // δ = (ρ_core/ρ_avg)^n × r_core [m]
    double computeFieldAttenuation();        // α = exp(-δ/λ) [unitless]
    double computeEffectiveUg(double Ug);    // Ug_eff = Ug × α [m/s²]
    double computeCoreDensity();             // ρ_core from M, r [kg/m³]
    std::string getEquationText();
    void printVariables();
};

#endif // CORE_PENETRATION_MODULE_H
