// TimeReversalZoneModule.h
// Defines the time-reversal zone (TRZ) via: f_TRZ = r_TRZ / r — fraction of orbit radius
// supporting CPT-reversed dynamics. f_TRZ = 0.1 (canonical). |g_TRZ| = g (1 + f_TRZ).
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L7891

#ifndef TIME_REVERSAL_ZONE_MODULE_H
#define TIME_REVERSAL_ZONE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class TimeReversalZoneModule {
private:
    std::map<std::string, double> variables;

public:
    TimeReversalZoneModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeTRZFraction();             // f_TRZ = r_TRZ/r [unitless, ~0.1]
    double computeTRZRadius();               // r_TRZ = f_TRZ × r [m]
    double computeGravityCorrection(double g); // g_eff = g × (1 + f_TRZ) [m/s²]
    double computeCPTSymmetry();             // CPT-factor signature ±1
    double computeTRZEnergy();               // E_TRZ = m g_eff r_TRZ [J]
    std::string getEquationText();
    void printVariables();
};

#endif // TIME_REVERSAL_ZONE_MODULE_H
