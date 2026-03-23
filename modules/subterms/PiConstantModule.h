// PiConstantModule.h
// Encodes π through UQFF Leibniz/arctan series: π = 4 Σ_{k=0}^∞ (-1)^k / (2k+1).
// Used in rotational DPM sphere area, orbital period factors, and F_U_Bi normalization.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L4999

#ifndef PI_CONSTANT_MODULE_H
#define PI_CONSTANT_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class PiConstantModule {
private:
    std::map<std::string, double> variables;

public:
    PiConstantModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computePiLeibniz(int terms);      // π/4 = Σ (-1)^k/(2k+1)
    double computeSphereArea();              // A = 4π r²
    double computeSphereVolume();            // V = (4/3)π r³
    double computeOrbitalPeriod();           // T = 2π r / v [s]
    double getPi();                          // returns M_PI constant
    std::string getEquationText();
    void printVariables();
};

#endif // PI_CONSTANT_MODULE_H
