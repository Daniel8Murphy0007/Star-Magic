// SurfaceMagneticFieldModule.h
// Stellar surface magnetic field: B_s = μ_0 M_mag / (4π r_star³) — dipole approximation.
// Magnetar: B ≈ 4.4e13 T = B_crit. Modulates Ug1 dipole and MUGE suppression factor.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L7563

#ifndef SURFACE_MAGNETIC_FIELD_MODULE_H
#define SURFACE_MAGNETIC_FIELD_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class SurfaceMagneticFieldModule {
private:
    std::map<std::string, double> variables;

public:
    SurfaceMagneticFieldModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeSurfaceField();            // B_s = μ_0 M_mag / (4π r³) [T]
    double computeFieldAtR(double r);        // B(r) = B_s (r_star/r)³ [T]
    double computeCyclotronFrequency();      // f_c = q B / (2π m_e) [Hz]
    double computeSynchrotronPower();        // P = σ_T c B² / (6π) × γ² [W]
    double computeMUGESuppression();         // (1 - B/B_crit) factor
    std::string getEquationText();
    void printVariables();
};

#endif // SURFACE_MAGNETIC_FIELD_MODULE_H
