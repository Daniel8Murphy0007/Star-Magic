// NGC1300UQFFModule.cpp
// Implementation of MUGE for barred spiral galaxy NGC 1300.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "NGC1300UQFFModule.h"
#include <complex>

NGC1300UQFFModule::NGC1300UQFFModule() {
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
    double kpc_val = 3.086e19;

    // NGC 1300 parameters
    variables["M_visible"] = 7e10 * M_sun_val;
    variables["M_DM"] = 3e10 * M_sun_val;
    variables["M"] = variables["M_visible"] + variables["M_DM"];
    variables["M0"] = variables["M"];
    variables["SFR"] = 1 * M_sun_val / variables["year_to_s"];
    variables["r"] = 11.79e3 * kpc_val;
    variables["z"] = 0.005;
    variables["v_arm"] = 200e3;
    variables["t"] = 1e9 * variables["year_to_s"];

    variables["rho_fluid"] = 1e-21;
    variables["V"] = 1e50;
    variables["B"] = 1e-5;
    variables["B_crit"] = 1e11;
    variables["Delta_x"] = 1e-10;
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;

    variables["A"] = 1e-10;
    variables["k"] = 1e20;
    variables["omega"] = 1e-15;
    variables["x"] = 0.0;
    variables["v"] = variables["v_arm"];
    variables["sigma"] = 1e3 * kpc_val;

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
    variables["omega_spin"] = 1e-4;
    variables["I_dipole"] = 1e20;
    variables["A_dipole"] = 1e15;
    variables["H_aether"] = 1e-6;
    variables["delta_rho_over_rho"] = 1e-5;

    variables["scale_macro"] = 1e-12;
    variables["f_TRZ"] = 0.1;
    variables["f_sc"] = 1.0;
    variables["v_r"] = 1e3;
    variables["rho"] = variables["rho_fluid"];
}

void NGC1300UQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) variables[name] = value;
    else variables[name] = value;
    if (name == "Delta_x") variables["Delta_p"] = variables["hbar"] / value;
    else if (name == "M") {
        variables["M_visible"] = 0.7 * value;
        variables["M_DM"] = 0.3 * value;
        variables["M0"] = value;
    }
}

void NGC1300UQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) variables[name] += delta;
    else variables[name] = delta;
}
void NGC1300UQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

double NGC1300UQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

double NGC1300UQFFModule::computeMsfFactor(double t) {
    return variables["SFR"] * t / variables["M0"];
}

double NGC1300UQFFModule::computeRt(double t) {
    return variables["r"] + variables["v_r"] * t;
}

double NGC1300UQFFModule::computeFenv(double t) {
    double F_bar = 0.1 * (variables["G"] * variables["M"]) / (variables["r"] * variables["r"]);
    double F_SF = variables["k_SF"] * variables["SFR"] / 1.989e30;
    double F_wave = variables["rho_fluid"] * std::pow(variables["v_arm"], 2);
    return F_bar + F_SF + F_wave;
}

double NGC1300UQFFModule::computeUg1(double t) {
    double mu_dipole = variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"];
    return mu_dipole * variables["B"];
}

double NGC1300UQFFModule::computeUg2(double t) {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}

double NGC1300UQFFModule::computeUg3prime(double t) {
    double M_bar = 0.2 * variables["M"];
    double r_bar = 0.3 * variables["r"];
    return (variables["G"] * M_bar) / (r_bar * r_bar);
}

double NGC1300UQFFModule::computeUg4(double t) {
    double E_react = 1e46 * std::exp(-0.0005 * t);
    return variables["k_4"] * E_react;
}

double NGC1300UQFFModule::computeUi(double t) {
    return variables["lambda_I"] * (variables["rho_vac_SCm"] / variables["rho_vac_UA"]) * variables["omega_i"] * std::cos(variables["pi"] * variables["t_n"]) * (1 + variables["F_RZ"]);
}

double NGC1300UQFFModule::computePsiIntegral(double r, double t) {
    double A = variables["A"];
    double m = 2.0;
    double omega = variables["omega"];
    double sigma = variables["sigma"];
    std::complex<double> psi_spiral(A * std::exp(-r*r / (2 * sigma * sigma)) * std::exp(std::complex<double>(0, m * 0 - omega * t)));
    return std::norm(psi_spiral);
}

double NGC1300UQFFModule::computeQuantumTerm(double t_Hubble_val, double r) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double psi_int = computePsiIntegral(r, variables["t"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * psi_int;
}

double NGC1300UQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

double NGC1300UQFFModule::computeDMTerm(double r) {
    double pert = variables["delta_rho_over_rho"];
    double curv = 3 * variables["G"] * variables["M"] / (r * r * r);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

double NGC1300UQFFModule::computeUgSum(double r) {
    double Ug_base = (variables["G"] * variables["M"]) / (r * r);
    variables["Ug1"] = computeUg1(variables["t"]);
    variables["Ug2"] = computeUg2(variables["t"]);
    variables["Ug3"] = computeUg3prime(variables["t"]);
    variables["Ug4"] = computeUg4(variables["t"]);
    return Ug_base + variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

double NGC1300UQFFModule::computeG(double t, double r) {
    variables["t"] = t;
    double msf_factor = computeMsfFactor(t);
    double m_factor = 1.0 + msf_factor;
    double Hz = computeHtz(variables["z"]);
    double expansion = 1.0 + Hz * t;
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

std::string NGC1300UQFFModule::getEquationText() {
    return "g_NGC1300(r, t) = (G * M(t) / r(t)^2) * (1 + H(t,z)) * (1 - B/B_crit) * (1 + F_env(t)) + "
           "(Ug1 + Ug2 + Ug3' + Ug4) + Ui + Lambda c^2/3 + quantum + fluid + DM\n"
           "F_env(t) = F_bar + F_SF + F_wave; F_bar = 0.1 G M / r^2; F_wave = rho v_arm^2;\n"
           "Ug3' = G M_bar / r_bar^2 (bar as external mass);\n"
           "Adaptations: Hubble ACS 2004; SFR=1 Msun/yr; M=1e11 Msun. g~2e36 m/s^2 at t=1 Gyr.";
}

void NGC1300UQFFModule::printVariables() {
    std::cout << "NGC 1300 Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}
