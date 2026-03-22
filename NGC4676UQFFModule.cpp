// NGC4676UQFFModule.cpp
// Implementation of MUGE for NGC 4676 (The Mice) — novel THz H_eff(z) and Ug2_THz terms.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "NGC4676UQFFModule.h"
#include <complex>

NGC4676UQFFModule::NGC4676UQFFModule() {
    variables["G"] = 6.6743e-11;
    variables["c"] = 3e8;
    variables["hbar"] = 1.0546e-34;
    variables["Lambda"] = 1.1e-52;
    variables["q"] = 1.602e-19;
    variables["pi"] = 3.141592653589793;
    variables["t_Hubble"] = 13.8e9 * 3.156e7;
    variables["year_to_s"] = 3.156e7;
    variables["H0"] = 70.0;
    variables["Mpc_to_m"] = 3.086e22;
    variables["kpc_to_m"] = 3.086e19;
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    double Msun = 1.989e30;

    // NGC 4676 (The Mice A+B) parameters
    variables["M_A"] = 5e10 * Msun;
    variables["M_B"] = 5e10 * Msun;
    variables["M_visible"] = variables["M_A"] + variables["M_B"];
    variables["M_DM"] = 2e10 * Msun;     // dark matter halo
    variables["M"] = variables["M_visible"] + variables["M_DM"];
    variables["M0"] = variables["M"];
    variables["SFR"] = 5.0 * Msun / variables["year_to_s"];
    variables["r"] = 50 * variables["kpc_to_m"];
    variables["d"] = 10 * variables["kpc_to_m"];     // separation A-B
    variables["v_rel"] = 400e3;
    variables["z"] = 0.022;
    variables["tau_merge"] = 1.7e8 * variables["year_to_s"];
    variables["rho_fluid"] = 1e-21;
    variables["V"] = 1e52;
    variables["B"] = 1e-5;
    variables["B_crit"] = 1e11;

    // THz Aether coupling (novel)
    variables["f_THz"] = 0.05;

    variables["Delta_x"] = 1e-10;
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;
    variables["A"] = 1e-10;
    variables["k"] = 1e18;
    variables["omega"] = 1e-15;
    variables["sigma"] = 1e19;
    variables["v_r"] = 1e3;
    variables["t"] = 0.0;

    variables["mu_0"] = 4 * variables["pi"] * 1e-7;
    variables["rho_vac_SCm"] = 7.09e-37;
    variables["rho_vac_UA"] = 7.09e-36;
    variables["lambda_I"] = 1.0;
    variables["omega_i"] = 1e-8;
    variables["t_n"] = 0.0;
    variables["F_RZ"] = 0.01;
    variables["k_SF"] = 1e-10;
    variables["omega_spin"] = 1e-15;
    variables["I_dipole"] = 1e21;
    variables["A_dipole"] = 1e15;
    variables["H_aether"] = 1e-6;
    variables["k_4"] = 1.0;
    variables["delta_rho_over_rho"] = 1e-5;
    variables["f_TRZ"] = 0.01;
}

void NGC4676UQFFModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
    if (name == "Delta_x") variables["Delta_p"] = variables["hbar"] / value;
    else if (name == "M_A" || name == "M_B") {
        variables["M_visible"] = variables["M_A"] + variables["M_B"];
        variables["M"] = variables["M_visible"] + variables["M_DM"];
        variables["M0"] = variables["M"];
    }
}
void NGC4676UQFFModule::addToVariable(const std::string& name, double delta) { variables[name] += delta; }
void NGC4676UQFFModule::subtractFromVariable(const std::string& name, double delta) { variables[name] -= delta; }

double NGC4676UQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// NOVEL: THz Aether-modulated Hubble parameter H_eff(z) = H(z)*(1 + f_THz*log(1+z))
double NGC4676UQFFModule::computeHeffz(double z_val) {
    return computeHtz(z_val) * (1.0 + variables["f_THz"] * std::log(1.0 + z_val));
}

double NGC4676UQFFModule::computeMmerge(double t) {
    return variables["M_B"] * std::exp(-t / variables["tau_merge"]);
}
double NGC4676UQFFModule::computeFenv(double t) {
    double F_tidal = variables["G"] * computeMmerge(t) / (variables["d"] * variables["d"]);
    double F_SF = variables["k_SF"] * variables["SFR"] / 1.989e30;
    double F_bridge = variables["rho_fluid"] * variables["v_rel"] * variables["v_rel"];
    return F_tidal + F_SF + F_bridge;
}
double NGC4676UQFFModule::computeUg1(double t) {
    double mu_dipole = variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"];
    return mu_dipole * variables["B"];
}
double NGC4676UQFFModule::computeUg2(double t) {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}
// NOVEL: time-growing THz Ug2 term — Ug2*(1 + f_THz*H_eff*t/t_Hubble)
double NGC4676UQFFModule::computeUg2THz(double t) {
    double H_eff = computeHeffz(variables["z"]);
    return computeUg2(t) * (1.0 + variables["f_THz"] * H_eff * t / variables["t_Hubble"]);
}
double NGC4676UQFFModule::computeUg3prime(double t) {
    return variables["G"] * computeMmerge(t) / (variables["d"] * variables["d"]);
}
double NGC4676UQFFModule::computeUg4(double t) {
    return variables["k_4"] * 1e43 * std::exp(-0.0005 * t);
}
double NGC4676UQFFModule::computeUi(double t) {
    return variables["lambda_I"] * (variables["rho_vac_SCm"] / variables["rho_vac_UA"]) * variables["omega_i"] * std::cos(variables["pi"] * variables["t_n"]) * (1 + variables["F_RZ"]);
}
double NGC4676UQFFModule::computeMsfFactor(double t) { return variables["SFR"] * t / variables["M0"]; }
double NGC4676UQFFModule::computeRt(double t) { return variables["r"] + variables["v_r"] * t; }

double NGC4676UQFFModule::computePsiIntegral(double r, double t) {
    std::complex<double> psi(variables["A"] * std::exp(-r * r / (2 * variables["sigma"] * variables["sigma"])) * std::exp(std::complex<double>(0, -variables["omega"] * t)));
    return std::norm(psi);
}
double NGC4676UQFFModule::computeQuantumTerm(double t_Hubble_val, double r) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * computePsiIntegral(r, variables["t"]);
}
double NGC4676UQFFModule::computeFluidTerm(double g_base) { return variables["rho_fluid"] * variables["V"] * g_base; }
double NGC4676UQFFModule::computeDMTerm(double r) {
    return (variables["M_visible"] + variables["M_DM"]) * (variables["delta_rho_over_rho"] + 3 * variables["G"] * variables["M"] / (r * r * r));
}
double NGC4676UQFFModule::computeUgSum(double r) {
    double Ug_base = variables["G"] * variables["M"] / (r * r);
    variables["Ug1"] = computeUg1(variables["t"]);
    variables["Ug2THz"] = computeUg2THz(variables["t"]);   // use THz-augmented version
    variables["Ug3"] = computeUg3prime(variables["t"]);
    variables["Ug4"] = computeUg4(variables["t"]);
    return Ug_base + variables["Ug1"] + variables["Ug2THz"] + variables["Ug3"] + variables["Ug4"];
}

double NGC4676UQFFModule::computeG(double t, double r) {
    variables["t"] = t;
    double m_factor = 1.0 + computeMsfFactor(t);
    double expansion = 1.0 + computeHeffz(variables["z"]) * t;   // use H_eff
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double f_env = computeFenv(t);
    double tr_factor = 1.0 + variables["f_TRZ"];
    double g_base = (variables["G"] * variables["M"] * m_factor / (r * r)) * expansion * sc_correction * (1.0 + f_env) * tr_factor;
    double ug_sum = computeUgSum(r) - variables["G"] * variables["M"] / (r * r);
    double lambda_term = variables["Lambda"] * variables["c"] * variables["c"] / 3.0;
    double ui_term = computeUi(t);
    double quantum_term = computeQuantumTerm(variables["t_Hubble"], r);
    double fluid_term = computeFluidTerm(g_base);
    double dm_term = computeDMTerm(r);
    return g_base + ug_sum + lambda_term + ui_term + quantum_term + fluid_term + dm_term;
}

std::string NGC4676UQFFModule::getEquationText() {
    return "g_NGC4676(r,t) = (G M(t)/r^2)(1+H_eff(z))(1-B/B_crit)(1+F_env) + Sum Ugi + Ui + Lambda c^2/3\n"
           "H_eff(z) = H(z)*(1 + f_THz*log(1+z))  [THz Aether-modulated expansion]\n"
           "Ug2_THz = Ug2*(1 + f_THz*H_eff*t/t_Hubble)  [time-growing magnetic term]\n"
           "F_env = F_tidal + F_SF + F_bridge; f_THz=0.05; M_A=M_B=5e10 Msun, d=10 kpc, v_rel=400 km/s";
}
void NGC4676UQFFModule::printVariables() {
    std::cout << "NGC 4676 (The Mice) Variables:\n";
    for (const auto& pair : variables) std::cout << pair.first << " = " << std::scientific << pair.second << "\n";
}
