// SolarWindBuoyancyModule.h
// Computes buoyancy modulation by solar wind density: (1 + ε_sw × ρ_vac,sw)
// ε_sw = 0.001 (unitless); ρ_vac,sw = 8e-21 J/m³; correction ≈ 8e-24 (negligible).
// In U_bi: the modulation factor stabilizes heliosphere/nebulae via minor SW density effect.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L2277

#ifndef SOLAR_WIND_BUOYANCY_MODULE_H
#define SOLAR_WIND_BUOYANCY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class SolarWindBuoyancyModule {
private:
    std::map<std::string, double> variables;
    double computeModulationFactorInternal();

public:
    SolarWindBuoyancyModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // ε_sw = 0.001 (fixed)
    double computeEpsilon_sw();

    // (1 + ε_sw × ρ_vac,sw) modulation factor
    double computeModulationFactor();

    // Example U_b1 with full modulation
    double computeU_b1();

    std::string getEquationText();
    void printVariables();
};

#endif // SOLAR_WIND_BUOYANCY_MODULE_H
