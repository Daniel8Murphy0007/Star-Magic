// MultiUQFFModule.h
// Docs 34a + 34b: Multi-system UQFF dispatcher for 15 systems.
// Systems: UniverseDiameter, HydrogenAtom, HydrogenResonancePToE, LagoonNebula,
//          SpiralsSupernovae, NGC6302, OrionNebula, UniverseGuide (batch 1)
//          GalaxiesGalore, StellarForge, SombreroGalaxy, Saturn, CrabNebula,
//          NewStars (batch 2)
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef MULTI_UQFF_MODULE_H
#define MULTI_UQFF_MODULE_H

#include "UQFFModuleBase.h"
#include <string>

class MultiUQFFModule : public UQFFModuleBase {
private:
    std::string current_system;

    void loadSystem(const std::string& system);
    double computeHtz() const;
    double computeUgSum() const;
    double computeQuantumPsi() const;
    double computeFluidTerm(double g_base) const;
    double computeMsfFactor(double t) const;
    double computeDMPert() const;
    double computeHRes(double t) const;  // H-res = A_res sin(2πf_res t)

public:
    // Constructor: provide system name
    explicit MultiUQFFModule(const std::string& system = "OrionNebula");

    // Switch system at runtime
    void setSystem(const std::string& system);
    std::string getSystem() const { return current_system; }

    // Core API
    double computeG(double t) override;
    std::string getEquationText() const override;

protected:
    void onVariableUpdated(const std::string& name, double value) override;
};

#endif // MULTI_UQFF_MODULE_H
