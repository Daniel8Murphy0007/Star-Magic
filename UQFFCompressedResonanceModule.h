// UQFFCompressedResonanceModule.h
// Modular C++ implementation of Compressed and Resonance UQFF Equations for Multi-System Evolution.
// Supports compressed g_UQFF(r,t) unified form; resonance mode adds oscillatory terms.
// Systems: YoungStars, Eagle, BigBang, M51, NGC1316, V838Mon, NGC1300, Guide.
// Usage: #include "UQFFCompressedResonanceModule.h"; mod.setSystem("Eagle"); mod.setMode("resonance"); mod.computeG(t);
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UQFF_COMPRESSED_RESONANCE_MODULE_H
#define UQFF_COMPRESSED_RESONANCE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class UQFFCompressedResonanceModule {
private:
    std::map<std::string, double> variables;
    std::string current_system;
    std::string mode;  // "compressed" or "resonance"
    double computeHtz(double z_val);
    double computeFenv(double t);
    double computeUgSum();
    double computePsiTotal(double t);
    double computeResonanceTerm(double t);
    double computeQuantumTerm(double t_Hubble_val);
    double computeFluidTerm(double g_base);
    double computeDMTerm();
    double computeMsfFactor(double t);

public:
    UQFFCompressedResonanceModule();
    void setSystem(const std::string& sys_name);
    void setMode(const std::string& m);
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    double computeG(double t, double r = 0.0);
    std::string getEquationText();
    void printVariables();
};

#endif // UQFF_COMPRESSED_RESONANCE_MODULE_H
