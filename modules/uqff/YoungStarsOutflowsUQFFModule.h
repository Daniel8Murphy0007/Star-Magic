// YoungStarsOutflowsUQFFModule.h
// Doc 35: NGC 346-analogue — time-growing outflow pressure P_outflow(t).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef YOUNG_STARS_OUTFLOWS_UQFF_MODULE_H
#define YOUNG_STARS_OUTFLOWS_UQFF_MODULE_H

#include "UQFFModuleBase.h"
#include <string>

class YoungStarsOutflowsUQFFModule : public UQFFModuleBase {
private:
    double computeHtz() const;
    double computePoutflow(double t) const;
    double computeMsfFactor(double t) const;
    double computeUgSum() const;
    double computeQuantumPsi() const;
    double computeFluidTerm(double g_base) const;
    double computeDMPert() const;

public:
    YoungStarsOutflowsUQFFModule();
    double computeG(double t) override;
    std::string getEquationText() const override;

protected:
    void onVariableUpdated(const std::string& name, double value) override;
};

#endif // YOUNG_STARS_OUTFLOWS_UQFF_MODULE_H
