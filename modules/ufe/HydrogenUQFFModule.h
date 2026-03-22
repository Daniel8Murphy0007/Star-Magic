// HydrogenUQFFModule.h
// UQFF for Compressed Space / Hydrogen Energy Levels (Doc 43.e).
// Equations: E_space = E0·SCF·CF·LF·HFF·PTF·QSF ~ 5.52e-104 J (page 85, layers=5)
//            Three-leg proofset: conservation~1, vac ratio 1.683e-97, q energy 4.136e-14 eV
// Copyright - Daniel T. Murphy.

#ifndef HYDROGEN_UQFF_MODULE_H
#define HYDROGEN_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>

enum class HydrogenSystemType {
    COMPRESSED_SPACE_85,   // page 85 — 5 layers → ~5.52e-104 J
    COMPRESSED_SPACE_86,   // page 86 — extended analysis
    HYDROGEN_LEVELS,       // standard hydrogen spectral levels
    GENERIC
};

class HydrogenUQFFModule {
private:
    std::map<std::string, double> variables;
    HydrogenSystemType current_system;

    double computeSCF();   // space compression factor
    double computeCF();    // coupling factor
    double computeLF(int layers);  // layer factor
    double computeHFF();   // Higgs frequency factor
    double computePTF();   // precession time factor
    double computeQSF();   // quantum scaling factor

public:
    HydrogenUQFFModule(HydrogenSystemType sys = HydrogenSystemType::COMPRESSED_SPACE_85);

    void setSystem(HydrogenSystemType sys);

    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core equations
    double computeEspace(int layers = 5);           // compressed space energy
    double computeHydrogenLevel(int n);             // E_n = -13.6 / n^2 eV
    double computeHydrogenTransition(int n1, int n2);  // delta E [eV]

    // Three-leg proofset
    double computeConservation(double E_in, double E_out);   // E_out/E_in ~ 1
    double computeVacDensityRatio();                          // ~ 1.683e-97
    double computeQuantumEnergy();                            // ~ 4.136e-14 eV

    double computeThreeLegProofset(double E_space);

    double computeUQFF(double t);

    std::string getEquationText();
    std::string getSolutions(double t);

    void printVariables();
};

#endif // HYDROGEN_UQFF_MODULE_H
