// GalacticDistanceModule.h
// Computes galactic scale factor d_g (virial radius r_200) used in buoyancy: M_bh / d_g
// d_g ≈ 2.55e20 m (default, ~8.3 kpc for Milky Way virial scale)
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L2894

#ifndef GALACTIC_DISTANCE_MODULE_H
#define GALACTIC_DISTANCE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class GalacticDistanceModule {
private:
    std::map<std::string, double> variables;

public:
    GalacticDistanceModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeD_g();              // d_g = r_200 virial radius [m]
    double computeM_BH_over_d_g();    // M_bh / d_g for buoyancy term [kg/m]
    std::string getEquationText();
    void printVariables();
};

#endif // GALACTIC_DISTANCE_MODULE_H
