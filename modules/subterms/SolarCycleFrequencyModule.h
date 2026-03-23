// SolarCycleFrequencyModule.h
// Solar 11-year activity cycle: f_sc = 1/(11 yr) ≈ 2.88e-9 Hz.
// Provides periodic solar wind modulation, magnetic polarity reversal factor, K_p index scaling.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L6284

#ifndef SOLAR_CYCLE_FREQUENCY_MODULE_H
#define SOLAR_CYCLE_FREQUENCY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class SolarCycleFrequencyModule {
private:
    std::map<std::string, double> variables;

public:
    SolarCycleFrequencyModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeCycleFrequency();          // f_sc = 1/(11×3.156e7) ≈ 2.88e-9 Hz
    double computeAmplitude(double t);       // A cos(2π f_sc t) [unitless modulator]
    double computePolarityPhase(double t);   // 22-year Hale cycle sign bit
    double computeSSNModulation(double ssn); // F_10.7 proxy for SFR→v_sw
    std::string getEquationText();
    void printVariables();
};

#endif // SOLAR_CYCLE_FREQUENCY_MODULE_H
