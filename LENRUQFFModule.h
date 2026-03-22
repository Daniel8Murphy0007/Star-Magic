// LENRUQFFModule.h
// Modular C++ implementation of UQFF for LENR (Low Energy Nuclear Reactions).
// Models Widom-Larsen mechanism: electron acceleration to 0.78 MeV threshold, neutron production.
// Three scenarios: metallic hydride cells, exploding wires, solar corona.
// Novel physics: Omega plasma frequency, modified Fermi decay rate with UQFF fields.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef LENR_UQFF_MODULE_H_ROOT
#define LENR_UQFF_MODULE_H_ROOT

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class LENRUQFFModule {
private:
    std::map<std::string, double> variables;
    std::string m_scenario;

    double computeOmegaPlasma();      // Omega = sqrt(4*pi*rho_e*e^2/m_e)
    double computeEfield();            // E = (m_e*c^2/e)*(Omega/c)
    double computeUm();                // Um = q*v_drift*B
    double computeUg1();
    double computeUg3();
    double computeUg4();
    // Modified neutron production rate via Fermi decay (Heaviside-corrected)
    double computeNeutronRate(double W_val, double Delta);
    double computeRhoE();              // electron density from field

public:
    LENRUQFFModule();
    void setScenario(const std::string& scenario);  // "hydride", "wires", "corona"
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    // Returns neutron production rate eta [cm^-2 s^-1]
    double computeEta(double t);
    // Returns UQFF gravity contribution
    double computeG(double t, double r);
    std::string getEquationText();
    void printVariables();
};

#endif // LENR_UQFF_MODULE_H_ROOT
