// UQFFCompressionModule.h
// Modular C++ implementation of UQFF Compression Cycle with 12-channel F_env.
// Dispatcher for multiple astrophysical systems: Magnetar, SagittariusA, Pillars.
// Novel physics: computePsiTotal = qvB + 2A*cos(kx)*cos(omega*t) + (2pi/13.8)*A*Re[exp(i(kx-omega*t))]
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UQFF_COMPRESSION_MODULE_H
#define UQFF_COMPRESSION_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class UQFFCompressionModule {
private:
    std::map<std::string, double> variables;
    std::string m_system;

    double computeHtz(double z_val);
    double computeMsfFactor(double t);
    // NEW: 12-channel environmental factor
    double computeFenv12(double t);
    double computeUg3prime();          // external mass contribution
    // NOVEL: psi_total = qvB + compression waves + (2pi/13.8) resonance
    double computePsiTotal(double t, double x);
    double computeUg1(double t);
    double computeUg2(double t);
    double computeUg4(double t);
    double computeUi(double t);
    double computeQuantumTerm(double t_Hubble_val, double r);
    double computeFluidTerm(double g_base);
    double computeDMTerm(double r);
    double computeUgSum(double r);

public:
    UQFFCompressionModule();
    // Set active system: "Magnetar", "SagittariusA", "Pillars", "Generic"
    void setSystem(const std::string& system_name);
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    double computeG(double t, double r);
    std::string getEquationText();
    void printVariables();
};

#endif // UQFF_COMPRESSION_MODULE_H
