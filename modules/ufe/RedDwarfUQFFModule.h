// RedDwarfUQFFModule.h
// UQFF for Red Dwarf LENR/Exploding Wire/Solar Corona/Higgs/NGC346/Pi-Calcs (Doc 43.c).
// Equations: W_mag (magnetic work), Um (vacuum magnetism), UH (Higgs-hydrogen),
//            Ug3 (star formation), neutron rate (k_η=2.75e8), delta_N series,
//            Pi/Basel series (S(2)=pi^2/6), buoyancy odd-n series (~-0.8887),
//            transmutation Q (~0.78 MeV), Higgs mass, branching ratio.
// Copyright - Daniel T. Murphy.

#ifndef RED_DWARF_UQFF_MODULE_H
#define RED_DWARF_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <utility>

enum class RedDwarfSystemType {
    LENR_CELL,
    EXPLODING_WIRE,
    SOLAR_CORONA,
    COLLIDER_HIGGS,
    NGC346,
    PI_CALCS
};

class RedDwarfUQFFModule {
private:
    std::map<std::string, double> variables;
    RedDwarfSystemType current_system;

    double computeNonLocalTerm(double t);
    double computePiSeriesInternal(double s, int terms);
    double computeBuoyancySeriesInternal(double x, int terms);

public:
    RedDwarfUQFFModule(RedDwarfSystemType sys = RedDwarfSystemType::LENR_CELL);

    void setSystem(RedDwarfSystemType sys);

    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core physics equations
    double computeWmag();                         // Magnetic work [eV]: 15e9*B_kG*R_km*(v/c)
    double computeUm(double t);                   // Vacuum magnetism energy
    double computeUH(double t, int n);            // Higgs-hydrogen coupling
    double computeUg3(double t, double r, double theta, int n);  // Star formation (k3 variant)
    double computeNeutronRate(double t);          // k_eta * exp(-NL) * Um / rho_UA
    double computeDeltaN(int n);                  // (2*pi)^n / 6
    double computePiSeriesS(double s);            // Basel: sum 1/n^s; S(2)=pi^2/6
    double computeBuoyancySeries(double x);       // sum_{n=odd} 1/x^{(pi+1)^n}; x=3 ~ -0.8887
    double computeTransmutationQ();               // (Mn - Mp - me)*c^2 ~ 0.78 MeV
    double computeHiggsMass();                    // k_Higgs * 125 * mu * kappa_F
    double computeBranchingRatio(const std::string& channel);

    double computeUQFF(double t);

    std::string getEquationText();
    std::string getSolutions(double t);

    void printVariables();
};

#endif // RED_DWARF_UQFF_MODULE_H
