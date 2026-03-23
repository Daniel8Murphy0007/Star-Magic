// StellarRotationModule.h
// Stellar angular velocity: ω_s = 2π / P_rot — tracks rotation period P_rot [s] and
// rotational velocity v_rot = ω r_star. Feeds Ug3 string-rotation term in UQFF.
// Watermark: Copyright - Daniel T. Murphy. Source: grok_share_b0a3dc1d.txt L7020

#ifndef STELLAR_ROTATION_MODULE_H
#define STELLAR_ROTATION_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class StellarRotationModule {
private:
    std::map<std::string, double> variables;

public:
    StellarRotationModule();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeAngularVelocity();         // ω_s = 2π / P_rot [rad/s]
    double computeEquatorialVelocity();      // v_eq = ω r_star [m/s]
    double computeBreakupVelocity();         // v_break = sqrt(G M/r) [m/s]
    double computeUg3StringTerm();           // Ug3 ∝ ω × T_string contribution
    double computeGyroscopicFactor();        // Ω_gyro = I ω / (M r²) [rad/s]
    std::string getEquationText();
    void printVariables();
};

#endif // STELLAR_ROTATION_MODULE_H
