// BigBangGravityUQFFModule.h
// Doc 38: Cosmic evolution / Big Bang UQFF.
// New terms: QG_term (Planck-scale), DM_term (Ω_DM fraction), GW_term (sinusoidal),
//            M(t)=M_total*(t/t_H), r(t)=c*t, z(t)=t_H/t−1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef BIG_BANG_GRAVITY_UQFF_MODULE_H
#define BIG_BANG_GRAVITY_UQFF_MODULE_H

#include "UQFFModuleBase.h"
#include <string>

class BigBangGravityUQFFModule : public UQFFModuleBase {
private:
    // Cosmic radial and mass evolution
    double computeR(double t) const;
    double computeM(double t) const;
    double computeZ(double t) const;
    double computeHtz(double z) const;

    // New physics terms
    double computeQGTerm(double t) const;
    double computeGWTerm(double t) const;
    double computeDMTerm(double g_base) const;

    double computeUgSum(double r) const;
    double computeQuantumPsi() const;
    double computeFluidTerm(double g_base) const;

public:
    BigBangGravityUQFFModule();
    double computeG(double t) override;
    std::string getEquationText() const override;
};

#endif // BIG_BANG_GRAVITY_UQFF_MODULE_H
