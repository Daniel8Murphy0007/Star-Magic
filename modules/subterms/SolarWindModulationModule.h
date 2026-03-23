// SolarWindModulationModule.h
// Solar wind dynamic modulation: v_sw(t) = v_0 × (1 + A sin(ω t)) where A≈0.2, v_0=400 km/s.
// ω = 2π f_sc. Feeds ε_sw = (v_sw - v_0)/v_0 into BuoyancyCouplingModule.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L6445

#ifndef SOLAR_WIND_MODULATION_MODULE_H
#define SOLAR_WIND_MODULATION_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class SolarWindModulationModule {
private:
    std::map<std::string, double> variables;

public:
    SolarWindModulationModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeVelocity(double t);        // v_sw(t) = v_0(1 + A sin(ω t)) [m/s]
    double computeEpsilonSW(double t);       // ε_sw = (v-v_0)/v_0 [unitless, ≈0.2]
    double computePressure(double t);        // P_sw = ½ ρ_sw v² [Pa]
    double computeModulationEnvelope();      // A = peak amplitude [~0.2]
    std::string getEquationText();
    void printVariables();
};

#endif // SOLAR_WIND_MODULATION_MODULE_H
