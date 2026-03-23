// MagneticMomentModule.h
// Computes stellar/galactic magnetic moment μ = I A_vort where A_vort is the vortex area.
// For solar magnetar class: μ ≈ 1e21 A·m². Provides Ug1 dipole contribution normalization.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L4427

#ifndef MAGNETIC_MOMENT_MODULE_H
#define MAGNETIC_MOMENT_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class MagneticMomentModule {
private:
    std::map<std::string, double> variables;

public:
    MagneticMomentModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeVortexArea();          // A = π r_vort² [m²]
    double computeMagneticMoment();      // μ = I A_vort [A·m²]
    double computeDipolePotential(double r); // V = μ/(4π μ_0 r³)
    double computeUg1Contribution();     // Ug1 dipole term [m/s²]
    std::string getEquationText();
    void printVariables();
};

#endif // MAGNETIC_MOMENT_MODULE_H
