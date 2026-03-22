// MUGEUQFFModuleFinal.h
// Final MUGE/UQFF module for 7 canonical SOURCE4 systems (Doc 42.a).
// Features: computeResonanceAcc() with 10 THz/aether/quantum/oscillatory terms;
//           getSolutions() with side-by-side compressed UQFF vs resonance H_res.
// Canonical systems: MAGNETAR_SGR1745, SGR_A, TAPESTRY_STARBIRTH, WESTERLUND2,
//                    PILLARS_CREATION, RINGS_RELATIVITY, STUDENTS_GUIDE.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef MUGE_UQFF_MODULE_FINAL_H
#define MUGE_UQFF_MODULE_FINAL_H

#include <map>
#include <string>
#include <cmath>
#include <complex>
#include <iostream>
#include <sstream>
#include <iomanip>

enum class SystemTypeFinal {
    MAGNETAR_SGR1745,
    SGR_A,
    TAPESTRY_STARBIRTH,
    WESTERLUND2,
    PILLARS_CREATION,
    RINGS_RELATIVITY,
    STUDENTS_GUIDE,
    GENERIC
};

class MUGEUQFFModuleFinal {
private:
    std::map<std::string, double> variables;
    SystemTypeFinal current_system;

    double computeHtz(double z);
    double computeFenv(double t);
    std::complex<double> computePsiTotal(double t);
    double computeQuantumTerm(double t_Hubble_val);
    double computeFluidTerm(double g_base);
    double computeDMTerm();
    double computeUgSum();

public:
    MUGEUQFFModuleFinal(SystemTypeFinal sys = SystemTypeFinal::GENERIC);

    void setSystem(SystemTypeFinal sys);

    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeUQFF(double t);          // Compressed g_UQFF
    double computeHres(double t);          // Resonance: base + 10-term acc
    double computeResonanceAcc(double t);  // 10-term resonance accelerations
    double computeDuniverse();

    std::string getEquationText();
    std::string getSolutions(double t);    // Side-by-side compressed + resonance

    void printVariables();
};

#endif // MUGE_UQFF_MODULE_FINAL_H
