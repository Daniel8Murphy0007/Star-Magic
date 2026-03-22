// MultiUQFFCompressionModule.h
// Docs 39, 40, 41: Multi-system UQFF Compression Cycle 2 — 29 systems.
// Modular F_env(t) = ΣF_i(t), Ug3' generalized, ψ_total consolidated.
// Supports: MagnetarSGR1745, SagittariusA, TapestryStarbirth, Westerlund2,
//   PillarsCreation, RingsRelativity, NGC2525, NGC3603, BubbleNebula,
//   AntennaeGalaxies, HorseheadNebula, NGC1275, NGC1792, HubbleUltraDeepField,
//   StudentsGuideUniverse, SombreroGalaxy, Saturn, EagleNebula, CrabNebula,
//   HydrogenAtom, HydrogenResonance, OrionNebula, GalaxiesGalore, NewStars,
//   StellarForge, LagoonNebula, SpiralsSupernovae, NGC6302, UniverseDiameter.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef MULTI_UQFF_COMPRESSION_MODULE_H
#define MULTI_UQFF_COMPRESSION_MODULE_H

#include "UQFFModuleBase.h"
#include <string>

class MultiUQFFCompressionModule : public UQFFModuleBase {
private:
    std::string current_system;

    void loadSystem(const std::string& system);
    double computeHtz() const;
    double computeF_env(double t) const;
    double computeMsfFactor(double t) const;
    double computeUgSum() const;
    double computeQuantumPsi() const;
    double computeFluidTerm(double g_base) const;
    double computeDMPert() const;
    double computeH_res(double t) const;  // for atomic/quantum systems

public:
    explicit MultiUQFFCompressionModule(const std::string& system = "MagnetarSGR1745");

    void setSystem(const std::string& system);
    std::string getSystem() const { return current_system; }

    double computeG(double t) override;
    std::string getEquationText() const override;

protected:
    void onVariableUpdated(const std::string& name, double value) override;
};

#endif // MULTI_UQFF_COMPRESSION_MODULE_H
