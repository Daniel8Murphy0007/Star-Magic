// EagleUQFFModule.h
// Doc 36: Eagle Nebula / Pillars of Creation UQFF.
// New terms: W_stellar = ρ v_wind², P_rad = L ρ/(4π r² c m_H) (density-scaled).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef EAGLE_UQFF_MODULE_H
#define EAGLE_UQFF_MODULE_H

#include "UQFFModuleBase.h"
#include <string>

class EagleUQFFModule : public UQFFModuleBase {
private:
    double computeHtz() const;
    double computeWstellar() const;
    double computePrad() const;
    double computeMsfFactor(double t) const;
    double computeUgSum() const;
    double computeQuantumPsi() const;
    double computeFluidTerm(double g_base) const;
    double computeDMPert() const;

public:
    EagleUQFFModule();
    double computeG(double t) override;
    std::string getEquationText() const override;

protected:
    void onVariableUpdated(const std::string& name, double value) override;
};

#endif // EAGLE_UQFF_MODULE_H
