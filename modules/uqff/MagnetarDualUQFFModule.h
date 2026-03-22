// MagnetarDualUQFFModule.h
// Doc 39b: SGR 1745-2900 dual-mode UQFF (compressed vs frequency).
// Mode "compressed": unified H(t,z)/F_env(t) + BH Ug3', ~1.782e39 m/s²
// Mode "frequency":  full resonance sum a_DPM+a_THz+..., ~1.773e-9 m/s²
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef MAGNETAR_DUAL_UQFF_MODULE_H
#define MAGNETAR_DUAL_UQFF_MODULE_H

#include "UQFFModuleBase.h"
#include <string>

class MagnetarDualUQFFModule : public UQFFModuleBase {
private:
    std::string current_mode;  // "compressed" or "frequency"

    // --- Compressed mode helpers ---
    double computeHtz(double z) const;
    double computeF_env(double t) const;
    double computeQuantumTerm() const;
    double computeFluidTerm(double g_base) const;
    double computeUgSum() const;
    double computeDMPertTerm() const;

    // --- Frequency mode helpers ---
    double computeADPM() const;
    double computeATHz() const;
    double computeAvac_diff() const;
    double computeAsuper_freq() const;
    double computeAaether_res(double t) const;
    double computeUg4i(double t) const;
    double computeAquantum_freq() const;
    double computeAAether_freq() const;
    double computeAfluid_freq() const;
    double computeOsc_term(double t) const;
    double computeAexp_freq(double t) const;
    double computeFTRZ() const;

    double computeCompressed(double t);
    double computeFrequency(double t);

public:
    explicit MagnetarDualUQFFModule(const std::string& mode = "compressed");

    void setMode(const std::string& mode);
    std::string getMode() const { return current_mode; }

    double computeG(double t) override;
    std::string getEquationText() const override;

protected:
    void onVariableUpdated(const std::string& name, double value) override;
};

#endif // MAGNETAR_DUAL_UQFF_MODULE_H
