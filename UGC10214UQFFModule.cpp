// UGC10214UQFFModule.cpp
// Implementation of MUGE for UGC 10214 (Tadpole Galaxy) — VV29c tidal dynamics.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "UGC10214UQFFModule.h"
#include <complex>

UGC10214UQFFModule::UGC10214UQFFModule() {
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

    // UGC 10214 Tadpole Galaxy parameters
    variables["M_visible"] = 7e10 * Msun;
    variables["M_DM"] = 3e10 * Msun;
    variables["M"] = variables["M_visible"] + variables["M_DM"];
    variables["M0"] = variables["M"];
    variables["SFR"] = 4.67 * Msun / variables["year_to_s"];
    variables["r"] = 55 * variables["kpc_to_m"];
    variables["z"] = 0.032;
    variables["rho_fluid"] = 1e-21;
    variables["V"] = 1e52;
    variables["B"] = 1e-5;
    variables["B_crit"] = 1e11;

    // Tidal companion VV29c
    variables["M_dwarf"] = 3.5e9 * Msun;
    variables["d_dwarf"] = 110 * variables["kpc_to_m"];
    variables["v_tail"] = 400e3;
    variables["tau_merge"] = 2.5e8 * variables["year_to_s"];

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
    variables["scale_macro"] = 1e-12;
}

void UGC10214UQFFModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
    if (name == "Delta_x") variables["Delta_p"] = variables["hbar"] / value;
    else if (name == "M_visible" || name == "M_DM") {
        variables["M"] = variables["M_visible"] + variables["M_DM"];
        variables["M0"] = variables["M"];
    }
}
void UGC10214UQFFModule::addToVariable(const std::string& name, double delta) { variables[name] += delta; }
void UGC10214UQFFModule::subtractFromVariable(const std::string& name, double delta) { variables[name] -= delta; }

double UGC10214UQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}
double UGC10214UQFFModule::computeMmerge(double t) {
    return variables["M_dwarf"] * std::exp(-t / variables["tau_merge"]);
}
double UGC10214UQFFModule::computeFenv(double t) {
    double F_tidal = variables["G"] * computeMmerge(t) / (variables["d_dwarf"] * variables["d_dwarf"]);
    double F_SF = variables["k_SF"] * variables["SFR"] / 1.989e30;
    double F_tail = variables["rho_fluid"] * variables["v_tail"] * variables["v_tail"];
    return F_tidal + F_SF + F_tail;
}
double UGC10214UQFFModule::computeUg1(double t) {
    double mu_dipole = variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"];
    return mu_dipole * variables["B"];
}
double UGC10214UQFFModule::computeUg2(double t) {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}
double UGC10214UQFFModule::computeUg3prime(double t) {
    return variables["G"] * computeMmerge(t) / (variables["d_dwarf"] * variables["d_dwarf"]);
}
double UGC10214UQFFModule::computeUg4(double t) {
    return variables["k_4"] * 1e42 * std::exp(-0.0005 * t);
}
double UGC10214UQFFModule::computeUi(double t) {
    return variables["lambda_I"] * (variables["rho_vac_SCm"] / variables["rho_vac_UA"]) * variables["omega_i"] * std::cos(variables["pi"] * variables["t_n"]) * (1 + variables["F_RZ"]);
}
double UGC10214UQFFModule::computeMsfFactor(double t) { return variables["SFR"] * t / variables["M0"]; }
double UGC10214UQFFModule::computeRt(double t) { return variables["r"] + variables["v_r"] * t; }

double UGC10214UQFFModule::computePsiIntegral(double r, double t) {
    std::complex<double> psi(variables["A"] * std::exp(-r * r / (2 * variables["sigma"] * variables["sigma"])) * std::exp(std::complex<double>(0, -variables["omega"] * t)));
    return std::norm(psi);
}
double UGC10214UQFFModule::computeQuantumTerm(double t_Hubble_val, double r) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * computePsiIntegral(r, variables["t"]);
}
double UGC10214UQFFModule::computeFluidTerm(double g_base) { return variables["rho_fluid"] * variables["V"] * g_base; }
double UGC10214UQFFModule::computeDMTerm(double r) {
    return (variables["M_visible"] + variables["M_DM"]) * (variables["delta_rho_over_rho"] + 3 * variables["G"] * variables["M"] / (r * r * r));
}
double UGC10214UQFFModule::computeUgSum(double r) {
    double Ug_base = variables["G"] * variables["M"] / (r * r);
    variables["Ug1"] = computeUg1(variables["t"]);
    variables["Ug2"] = computeUg2(variables["t"]);
    variables["Ug3"] = computeUg3prime(variables["t"]);
    variables["Ug4"] = computeUg4(variables["t"]);
    return Ug_base + variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

double UGC10214UQFFModule::computeG(double t, double r) {
    variables["t"] = t;
    double m_factor = 1.0 + computeMsfFactor(t);
    double expansion = 1.0 + computeHtz(variables["z"]) * t;
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

std::string UGC10214UQFFModule::getEquationText() {
    return "g_UGC10214(r,t) = (G M(t)/r^2)(1+H(t,z))(1-B/B_crit)(1+F_env) + Sum Ugi + Ui + Lambda c^2/3 + quantum + fluid + DM\n"
           "F_env = F_tidal(Mmerge) + F_SF + F_tail; F_tidal = G*M_dwarf*exp(-t/tau_merge)/d^2; F_tail = rho*v_tail^2\n"
           "M_dwarf=3.5e9 Msun, v_tail=400 km/s, tau_merge=2.5e8 yr. HST ACS Tadpole Galaxy.";
}
void UGC10214UQFFModule::printVariables() {
    std::cout << "UGC 10214 Tadpole Galaxy Variables:\n";
    for (const auto& pair : variables) std::cout << pair.first << " = " << std::scientific << pair.second << "\n";
}
