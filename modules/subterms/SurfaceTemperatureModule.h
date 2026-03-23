// SurfaceTemperatureModule.h
// Effective stellar surface temperature: T_s = (L / (4π σ r²))^(1/4) [K].
// From Stefan-Boltzmann. Sun: T_eff ≈ 5778 K. Magnetar crust: ~1e9 K. Used in SFR modulation.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L7730

#ifndef SURFACE_TEMPERATURE_MODULE_H
#define SURFACE_TEMPERATURE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class SurfaceTemperatureModule {
private:
    std::map<std::string, double> variables;

public:
    SurfaceTemperatureModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeSurfaceTemp();             // T_s = (L/4πσr²)^0.25 [K]
    double computeLuminosity();             // L = 4πr²σT⁴ [W]
    double computePeakWavelength();          // λ_peak = b/T [m] (Wien)
    double computeThermalVelocity();         // v_th = sqrt(3kT/m) [m/s]
    double computeBlackbodyPower();          // P = σ T⁴ [W/m²]
    std::string getEquationText();
    void printVariables();
};

#endif // SURFACE_TEMPERATURE_MODULE_H
