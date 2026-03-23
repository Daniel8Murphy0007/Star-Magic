// HeliosphereThicknessModule.h
// Models heliosphere thickness L_helio = (P_SW / P_ISMF)^(1/6) × r_termshock
// ~100 AU; modulates solar wind buoyancy terms and inner heliosphere boundary conditions.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L3866

#ifndef HELIOSPHERE_THICKNESS_MODULE_H
#define HELIOSPHERE_THICKNESS_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class HeliosphereThicknessModule {
private:
    std::map<std::string, double> variables;

public:
    HeliosphereThicknessModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeTerminationShock();    // r_TS ≈ 90 AU [m]
    double computeHeliopause();          // r_HP ≈ 120 AU [m]
    double computeHelioThickness();      // L = r_HP - r_TS [m]
    double computePressureBalance();     // P_SW = P_ISMF balance factor
    std::string getEquationText();
    void printVariables();
};

#endif // HELIOSPHERE_THICKNESS_MODULE_H
