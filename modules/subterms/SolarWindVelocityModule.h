// SolarWindVelocityModule.h
// Baseline solar wind velocity range: v_sw ∈ [400, 800] km/s depending on solar activity.
// Slow wind: 400 km/s (equatorial); fast wind: 700–800 km/s (coronal holes). CME: up to 2000 km/s.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L6631

#ifndef SOLAR_WIND_VELOCITY_MODULE_H
#define SOLAR_WIND_VELOCITY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class SolarWindVelocityModule {
private:
    std::map<std::string, double> variables;

public:
    SolarWindVelocityModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeSlowWind();                // v_slow = 4e5 m/s
    double computeFastWind();                // v_fast = 7.5e5 m/s
    double computeAlfvenMachNumber();        // M_A = v_sw / v_Alfven [unitless]
    double computeKineticPressure();         // P_kin = ½ ρ v² [Pa]
    double computePDLFromV(double v);        // Pressure-driven lens parameter
    std::string getEquationText();
    void printVariables();
};

#endif // SOLAR_WIND_VELOCITY_MODULE_H
