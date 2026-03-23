// GalacticSpinModule.h
// Models galactic spin Ω_g = v_disk / r_disk in UQFF buoyancy equation.
// Ω_g = 7.3e-16 rad/s (Milky Way default); appears in U_bi = -β_i U_gi Ω_g (M_bh/d_g) E_react
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L3504

#ifndef GALACTIC_SPIN_MODULE_H
#define GALACTIC_SPIN_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class GalacticSpinModule {
private:
    std::map<std::string, double> variables;

public:
    GalacticSpinModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeOmega_g();           // Ω_g = v_disk / r_disk [rad/s]
    double computeAngularMomentum();   // L = M v r [kg m² / s]
    double computePeriod();            // T = 2π / Ω_g [s]
    std::string getEquationText();
    void printVariables();
};

#endif // GALACTIC_SPIN_MODULE_H
