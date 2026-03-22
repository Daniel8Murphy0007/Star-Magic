// LENRUQFFModule.cpp
// Implementation of UQFF for LENR — electron acceleration, neutron production, 3 scenarios.
// Bug fix: uses ternary conditional for Heaviside step (std::theta does not exist in C++).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "LENRUQFFModule.h"

LENRUQFFModule::LENRUQFFModule() : m_scenario("hydride") {
    variables["G"] = 6.6743e-11;
    variables["c"] = 3e8;
    variables["hbar"] = 1.0546e-34;
    variables["pi"] = 3.141592653589793;
    variables["q"] = 1.602e-19;          // electron charge
    variables["m_e"] = 9.109e-31;        // electron mass
    variables["m_p"] = 1.673e-27;        // proton mass
    variables["m_n"] = 1.675e-27;        // neutron mass
    variables["year_to_s"] = 3.156e7;

    // Threshold energy: W > Delta = 0.78 MeV (electron capture threshold)
    variables["Delta"] = 0.78e6 * variables["q"];          // J
    variables["G_F"] = 1.166e-5;           // Fermi constant (GeV^-2, natural units)
    // m_tilde: effective mass parameter ~1 MeV/c^2
    variables["m_tilde"] = 1e6 * variables["q"] / (variables["c"] * variables["c"]);

    // Vacuum densities
    variables["rho_vac_UA"] = 7.09e-36;
    variables["rho_vac_SCm"] = 7.09e-37;

    // Default: hydride scenario
    variables["E_field"] = 2e11;         // V/m
    variables["eta_target"] = 1e13;      // cm^-2 s^-1
    variables["B"] = 0.1;               // T (1 kG)
    variables["v_drift"] = 1e6;          // m/s (electron drift)
    variables["rho_e_num"] = 8e28;       // #/m^3 (electron number density in Pd hydride)
    variables["I_Alfven"] = 17e3;        // A (Alfven current)

    variables["t"] = 0.0;
    variables["r"] = 1e-3;              // characteristic scale
    variables["k_4"] = 1.0;
    variables["omega_spin"] = 1e10;
    variables["I_dipole"] = 1e8;
    variables["A_dipole"] = 1e3;
    variables["H_aether"] = 1e-4;
    variables["mu_0"] = 4 * variables["pi"] * 1e-7;
}

void LENRUQFFModule::setScenario(const std::string& scenario) {
    m_scenario = scenario;
    if (scenario == "hydride") {
        variables["E_field"] = 2e11;
        variables["eta_target"] = 1e13;
        variables["B"] = 0.1;
    } else if (scenario == "wires") {
        variables["E_field"] = 28.8e11;
        variables["eta_target"] = 1e8;
        variables["B"] = 1.0;       // ~1 T in wire plasma
    } else if (scenario == "corona") {
        variables["E_field"] = 1.2e-3;
        variables["eta_target"] = 7e-3;
        variables["B"] = 1e-4;     // ~1 G solar corona
    }
}

void LENRUQFFModule::updateVariable(const std::string& name, double value) { variables[name] = value; }
void LENRUQFFModule::addToVariable(const std::string& name, double delta) { variables[name] += delta; }
void LENRUQFFModule::subtractFromVariable(const std::string& name, double delta) { variables[name] -= delta; }

// Omega = sqrt(4*pi*rho_e*e^2/m_e) — plasma angular frequency
double LENRUQFFModule::computeOmegaPlasma() {
    return std::sqrt(4.0 * variables["pi"] * variables["rho_e_num"] * variables["q"] * variables["q"] / variables["m_e"]);
}
// Self-consistent field: E = (m_e*c^2/e)*(Omega/c)
double LENRUQFFModule::computeEfield() {
    double Omega = computeOmegaPlasma();
    return (variables["m_e"] * variables["c"] * variables["c"] / variables["q"]) * (Omega / variables["c"]);
}
double LENRUQFFModule::computeRhoE() {
    return variables["rho_e_num"] * variables["m_e"];
}
double LENRUQFFModule::computeUm() {
    return variables["q"] * variables["v_drift"] * variables["B"];
}
double LENRUQFFModule::computeUg1() {
    return variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"] * variables["B"];
}
double LENRUQFFModule::computeUg3() {
    return (variables["G"] * 1.989e30 / (variables["r"] * variables["r"])) * (variables["rho_e_num"] * variables["m_e"] / variables["rho_vac_UA"]);
}
double LENRUQFFModule::computeUg4() {
    return variables["k_4"] * 1e10 * std::exp(-0.0005 * variables["t"]);
}

// Modified Fermi decay rate — Heaviside via ternary: (W > Delta) ? rate : 0
double LENRUQFFModule::computeNeutronRate(double W_val, double Delta) {
    if (W_val <= Delta) return 0.0;
    double G_F_SI = variables["G_F"] / (1.0e18 * variables["q"]);  // approximate SI conversion
    double m_c2 = variables["m_tilde"] * variables["c"] * variables["c"];
    double W_diff = W_val - Delta;
    // Fermi decay: Gamma ~ G_F^2*(m_tilde*c^2)^4/(2*pi*hbar^3) * (W-Delta)^2
    double rate = (G_F_SI * G_F_SI * std::pow(m_c2, 4) / (2.0 * variables["pi"] * std::pow(variables["hbar"], 3))) * W_diff * W_diff;
    return rate;
}

double LENRUQFFModule::computeEta(double t) {
    variables["t"] = t;
    // Electron energy from external field
    double W_electron = variables["q"] * variables["E_field"] * 1e-10;  // energy over ~10 pm path
    double Um_val = computeUm();
    double Ug3_val = computeUg3();
    // Total electron energy including UQFF
    double W_total = W_electron + std::abs(Um_val) + std::abs(Ug3_val);

    double rate = computeNeutronRate(W_total, variables["Delta"]);
    // Scale to neutron flux (phenomenological calibration)
    double eta = variables["rho_e_num"] * rate / (variables["m_n"] * variables["c"]);
    return eta;
}

double LENRUQFFModule::computeG(double t, double r) {
    variables["t"] = t;
    variables["r"] = r;
    double Ug_sum = computeUm() + computeUg1() + computeUg3() + computeUg4();
    double lambda_term = 1.1e-52 * variables["c"] * variables["c"] / 3.0;
    return Ug_sum + lambda_term;
}

std::string LENRUQFFModule::getEquationText() {
    return "eta_LENR(t) = rho_e * Gamma(W_total, Delta) / (m_n * c)\n"
           "Gamma = G_F^2*(m_tilde*c^2)^4/(2*pi*hbar^3)*(W-Delta)^2 * Theta(W-Delta)\n"
           "W_total = q*E*d + |Um| + |Ug3|; Omega = sqrt(4*pi*rho_e*e^2/m_e)\n"
           "3 scenarios: hydride (E=2e11 V/m, eta~1e13 cm^-2/s), wires, solar corona\n"
           "NOTE: Theta = ternary (W>Delta)?rate:0 — std::theta does not exist in C++";
}
void LENRUQFFModule::printVariables() {
    std::cout << "LENR UQFF Variables [scenario: " << m_scenario << "]:\n";
    for (const auto& pair : variables) std::cout << pair.first << " = " << std::scientific << pair.second << "\n";
}
