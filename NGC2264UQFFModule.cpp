// NGC2264UQFFModule.cpp
// Implementation of MUGE for Cone Nebula NGC 2264 stellar formation region.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "NGC2264UQFFModule.h"
#include <complex>

NGC2264UQFFModule::NGC2264UQFFModule() {
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
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    double M_sun_val = 1.989e30;

    // NGC 2264 parameters
    variables["M_visible"] = 80 * M_sun_val;
    variables["M_DM"] = 20 * M_sun_val;
    variables["M"] = variables["M_visible"] + variables["M_DM"];
    variables["M0"] = variables["M"];
    variables["SFR"] = 0.01 * M_sun_val / variables["year_to_s"];
    variables["r"] = 3.31e16;
    variables["z"] = 0.0008;
    variables["v_wind"] = 20e3;
    variables["t"] = 3e6 * variables["year_to_s"];

    variables["rho_fluid"] = 1e-20;
    variables["V"] = 1e48;
    variables["B"] = 1e-5;
    variables["B_crit"] = 1e11;
    variables["Delta_x"] = 1e-10;
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;

    variables["A"] = 1e-10;
    variables["k"] = 1e20;
    variables["omega"] = 1e-14;
    variables["x"] = 0.0;
    variables["v"] = variables["v_wind"];
    variables["sigma"] = 1e15;

    variables["Ug1"] = 0.0; variables["Ug2"] = 0.0;
    variables["Ug3"] = 0.0; variables["Ug4"] = 0.0; variables["Ui"] = 0.0;
    variables["mu_0"] = 4 * variables["pi"] * 1e-7;
    variables["rho_vac_SCm"] = 7.09e-37;
    variables["rho_vac_UA"] = 7.09e-36;
    variables["lambda_I"] = 1.0;
    variables["omega_i"] = 1e-8;
    variables["t_n"] = 0.0;
    variables["F_RZ"] = 0.01;
    variables["k_4"] = 1.0;
    variables["k_SF"] = 1e-10;
    variables["omega_spin"] = 1e-5;
    variables["I_dipole"] = 1e18;
    variables["A_dipole"] = 1e12;
    variables["H_aether"] = 1e-6;
    variables["delta_rho_over_rho"] = 1e-5;

    variables["scale_macro"] = 1e-12;
    variables["f_TRZ"] = 0.1;
    variables["f_sc"] = 1.0;
    variables["v_r"] = 1e3;
    variables["rho"] = variables["rho_fluid"];
}

void NGC2264UQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) variables[name] = value;
    else variables[name] = value;
    if (name == "Delta_x") variables["Delta_p"] = variables["hbar"] / value;
    else if (name == "M_visible" || name == "M_DM") {
        variables["M"] = variables["M_visible"] + variables["M_DM"];
        variables["M0"] = variables["M"];
    }
}
void NGC2264UQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) variables[name] += delta;
    else variables[name] = delta;
}
void NGC2264UQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

double NGC2264UQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}
double NGC2264UQFFModule::computeMsfFactor(double t) { return variables["SFR"] * t / variables["M0"]; }
double NGC2264UQFFModule::computeRt(double t) { return variables["r"] + variables["v_r"] * t; }

double NGC2264UQFFModule::computeFenv(double t) {
    double F_wind = variables["rho_fluid"] * std::pow(variables["v_wind"], 2);
    double F_SF = variables["k_SF"] * variables["SFR"] / 1.989e30;
    double F_erode = 0.05 * (t / (3e6 * variables["year_to_s"]));
    return F_wind + F_SF + F_erode;
}

double NGC2264UQFFModule::computeUg1(double t) {
    double mu_dipole = variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"];
    return mu_dipole * variables["B"];
}
double NGC2264UQFFModule::computeUg2(double t) {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}
double NGC2264UQFFModule::computeUg3prime(double t) {
    double M_star = 20 * 1.989e30;
    double r_star = 1e10;
    return (variables["G"] * M_star) / (r_star * r_star);
}
double NGC2264UQFFModule::computeUg4(double t) {
    return variables["k_4"] * 1e40 * std::exp(-0.0005 * t);
}
double NGC2264UQFFModule::computeUi(double t) {
    return variables["lambda_I"] * (variables["rho_vac_SCm"] / variables["rho_vac_UA"]) * variables["omega_i"] * std::cos(variables["pi"] * variables["t_n"]) * (1 + variables["F_RZ"]);
}

double NGC2264UQFFModule::computePsiIntegral(double r, double t) {
    double A = variables["A"], omega = variables["omega"], sigma = variables["sigma"];
    std::complex<double> psi(A * std::exp(-r*r / (2 * sigma * sigma)) * std::exp(std::complex<double>(0, 0 - omega * t)));
    return std::norm(psi);
}
double NGC2264UQFFModule::computeQuantumTerm(double t_Hubble_val, double r) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * computePsiIntegral(r, variables["t"]);
}
double NGC2264UQFFModule::computeFluidTerm(double g_base) { return variables["rho_fluid"] * variables["V"] * g_base; }
double NGC2264UQFFModule::computeDMTerm(double r) {
    double pert = variables["delta_rho_over_rho"];
    double curv = 3 * variables["G"] * variables["M"] / (r * r * r);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

double NGC2264UQFFModule::computeUgSum(double r) {
    double Ug_base = (variables["G"] * variables["M"]) / (r * r);
    variables["Ug1"] = computeUg1(variables["t"]);
    variables["Ug2"] = computeUg2(variables["t"]);
    variables["Ug3"] = computeUg3prime(variables["t"]);
    variables["Ug4"] = computeUg4(variables["t"]);
    return Ug_base + variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

double NGC2264UQFFModule::computeG(double t, double r) {
    variables["t"] = t;
    double m_factor = 1.0 + computeMsfFactor(t);
    double expansion = 1.0 + computeHtz(variables["z"]) * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double f_env = computeFenv(t);
    double tr_factor = 1.0 + variables["f_TRZ"];

    double g_base = (variables["G"] * variables["M"] * m_factor / (r * r)) * expansion * sc_correction * (1.0 + f_env) * tr_factor;
    double ug_sum = computeUgSum(r) - g_base;
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"]) / 3.0;
    double ui_term = computeUi(t);
    double quantum_term = computeQuantumTerm(variables["t_Hubble"], r);
    double fluid_term = computeFluidTerm(g_base);
    double dm_term = computeDMTerm(r);

    return g_base + ug_sum + lambda_term + ui_term + quantum_term + fluid_term + dm_term;
}

std::string NGC2264UQFFModule::getEquationText() {
    return "g_NGC2264(r,t) = (G M(t)/r^2)(1+H(t,z))(1-B/B_crit)(1+F_env) + Sum Ugi + Ui + Lambda c^2/3 + quantum + fluid + DM\n"
           "F_env = F_wind + F_SF + F_erode; F_wind = rho v_wind^2; Ug3'=G M_star/r_star^2\n"
           "Adaptations: Hubble ACS 2002; SFR=0.01 Msun/yr; M=100 Msun. g~1e-9 m/s^2 at t=3 Myr.";
}
void NGC2264UQFFModule::printVariables() {
    std::cout << "NGC 2264 Variables:\n";
    for (const auto& pair : variables) std::cout << pair.first << " = " << std::scientific << pair.second << "\n";
}
