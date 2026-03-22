// NebularUQFFModule.h
// UQFF for Nebular Cloud Analysis (Drawing 32), LENR, Higgs, NGC 346 (Doc 43.b).
// Equations: Electric field (14-18), neutron rate (15-17,19), transmutation (20),
//            Higgs (24), Ug3 star formation (28), blueshift (29), neutrinos (30),
//            decay (31), DNA (32), buoyancy (33), geometric star positions.
// Non-local term: [SSq]^{n26} * exp(-(pi + t))
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef NEBULAR_UQFF_MODULE_H
#define NEBULAR_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <complex>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <utility>

enum class NebularSystemType {
    NEBULA_CLOUD,
    NGC346,
    LENR_CELL,
    HIGGS_PHYSICS,
    GENERIC
};

class NebularUQFFModule {
private:
    std::map<std::string, double> variables;
    NebularSystemType current_system;

    double computeNonLocalTerm(double t, int n26);
    double computeUg3internal(double t, double r, double theta, int n);

public:
    NebularUQFFModule(NebularSystemType sys = NebularSystemType::GENERIC);

    void setSystem(NebularSystemType sys);

    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core equations (Drawing 32)
    double computeElectricField();               // Eq14-18 LENR
    double computeNeutronRate();                 // η eq15-17,19
    double computeTransmutationEnergy();         // Eq20
    double computeHiggsMass();                   // Eq24
    double computeStarFormationTemp(double t, double r);    // Eq28 Ug3
    double computeRadialVelocity(double dlambda_over_lambda); // Eq29 blueshift
    double computeNeutrinoEnergy(double t);      // Eq30
    double computeUniversalDecay(double t);      // Eq31
    double computeDNAEnergy(double t);           // Eq32
    double computeBuoyancyRatio(double V_little, double V_big); // Eq33
    double computeGeometricAngle(const std::vector<std::pair<double,double>>& positions);
    double computeAccuracy(const std::string& scenario);

    double computeUQFF(double t);                // Weighted sum

    std::string getEquationText();
    std::string getSolutions(double t);

    void printVariables();
};

#endif // NEBULAR_UQFF_MODULE_H
