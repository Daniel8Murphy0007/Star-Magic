// StressEnergyTensorModule.h
// Perfect-fluid stress-energy tensor: T_μν = (ρ+p)u_μ u_ν + p g_μν.
// Provides trace T = -ρ + 3p and energy conditions for UQFF vacuum state validation.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L7383

#ifndef STRESS_ENERGY_TENSOR_MODULE_H
#define STRESS_ENERGY_TENSOR_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class StressEnergyTensorModule {
private:
    std::map<std::string, double> variables;

public:
    StressEnergyTensorModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeTrace();                   // T = g^μν T_μν = -ρ + 3p
    double computePressure();                // p = ρ c² / 3 (radiation)
    double computeEnergyDensity();           // ρ = E / V [J/m³]
    double computeEinsteinTerm();            // G_μν = 8πG/c⁴ T_μν contribution
    double computeVacuumTerm();             // T_vac = -ρ_vac g_μν
    std::string getEquationText();
    void printVariables();
};

#endif // STRESS_ENERGY_TENSOR_MODULE_H
