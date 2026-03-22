// MUGEUQFFModule38.h
// Extends MUGEUQFFModule29 to 38 source documents, 14 system types.
// Adds: F_torque (spirals), F_shock (wind shocks), QG_term/DM_term/GW_term (cosmological composite F_cosmo).
// Auto-cascade: updating QG_term, DM_term, or GW_term recalculates F_cosmo.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef MUGE_UQFF_MODULE38_H
#define MUGE_UQFF_MODULE38_H

#include "MUGEUQFFModule29.h"

enum class SystemType38 {
    // All 8 from Doc41 (29)
    SOMBRERO_GALAXY,
    SATURN,
    M16_EAGLE,
    CRAB_NEBULA,
    HYDROGEN_ATOM,
    HYDROGEN_RESONANCE,
    UNIVERSE_DIAMETER,
    // 6 new from Doc42 (38)
    LAGOON_NEBULA,
    SPIRALS_SN,
    NGC6302,
    ORION_NEBULA,
    YOUNG_STARS_OUTFLOW,
    EAGLE_NEBULA,
    GRAVITY_BIGBANG,
    GENERIC
};

class MUGEUQFFModule38 {
private:
    std::map<std::string, double> variables;
    SystemType38 current_system;

    double computeHtz(double z);
    double computeFenv(double t);
    std::complex<double> computePsiTotal(double t);
    double computeQuantumTerm(double t_Hubble_val);
    double computeFluidTerm(double g_base);
    double computeDMTerm();
    double computeUgSum();
    void   recomputeFcosmo();

public:
    MUGEUQFFModule38(SystemType38 sys = SystemType38::GENERIC);

    void setSystem(SystemType38 sys);

    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeUQFF(double t);
    double computeHres(double t);
    double computeDuniverse();

    std::string getEquationText();
    std::string getSolutions(double t);

    void printVariables();
};

#endif // MUGE_UQFF_MODULE38_H
