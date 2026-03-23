// StellarMassModule.h
// Stellar mass evolution: M_s(t) = M_0 × (1 - γ_ML × t) with mass-loss rate γ_ML [s⁻¹].
// For solar: γ_ML ≈ 1.4e-14 s⁻¹; M Sun loses ~1e-14 M☉/yr. Feeds MUGE mass corrections.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L6823

#ifndef STELLAR_MASS_MODULE_H
#define STELLAR_MASS_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class StellarMassModule {
private:
    std::map<std::string, double> variables;

public:
    StellarMassModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeMassAtTime(double t);      // M_s = M_0(1 - γ_ML t) [kg]
    double computeMassLossRate();            // dM/dt = -γ_ML M_0 [kg/s]
    double computeLifetime();               // t_life = 1/γ_ML [s]
    double computeLuminosity();             // L ∝ M^4 main sequence [W]
    std::string getEquationText();
    void printVariables();
};

#endif // STELLAR_MASS_MODULE_H
