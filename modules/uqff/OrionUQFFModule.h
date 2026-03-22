// OrionUQFFModule.h
// Doc 34: Orion Nebula MUGE — H-alpha resonance, SFR M_sf(t), P_rad Trapezium, W_stellar.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef ORION_UQFF_MODULE_H
#define ORION_UQFF_MODULE_H

#include "UQFFModuleBase.h"
#include <string>

class OrionUQFFModule : public UQFFModuleBase {
private:
    double computeHtz() const;
    double computeMsfFactor(double t) const;
    double computeUgSum() const;
    double computeQuantumPsi() const;
    double computeFluidTerm(double g_base) const;
    double computeWstellar(double t) const;
    double computePrad() const;
    double computeResonantPsi(double t) const;
    double computeDMPert() const;

public:
    OrionUQFFModule();

    // Core API
    double computeG(double t) override;
    std::string getEquationText() const override;

protected:
    void onVariableUpdated(const std::string& name, double value) override;
};

#endif // ORION_UQFF_MODULE_H
