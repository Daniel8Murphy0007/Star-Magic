// NGC346UQFFModule.cpp
// Implementation of UQFF for NGC 346 — novel Um magnetism and Ug3 magnetic string disk.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "NGC346UQFFModule.h"
#include <complex>

NGC346UQFFModule::NGC346UQFFModule() {
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
    variables["pc_to_m"] = 3.086e16;
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    double Msun = 1.989e30;

    // NGC 346 (SMC) parameters
    variables["M_visible"] = 800 * Msun;
    variables["M_DM"] = 200 * Msun;
    variables["M"] = variables["M_visible"] + variables["M_DM"];
    variables["M0"] = variables["M"];
    variables["SFR"] = 0.1 * Msun / variables["year_to_s"];
    variables["r"] = 5 * variables["pc_to_m"];   // ~5 pc
    variables["z"] = 0.0006;
    variables["rho_gas"] = 1e-20;        // molecular cloud density
    variables["v_rad"] = -10e3;          // radial velocity: blueshift (m/s)
    variables["B"] = 1e-8;              // ~10 nT field in SMC
    variables["B_crit"] = 1e11;
    variables["V"] = 1e45;
    variables["rho_fluid"] = variables["rho_gas"];

    // Vacuum densities
    variables["rho_vac_SCm"] = 7.09e-37;
    variables["rho_vac_UA"] = 7.09e-36;
    variables["rho_plasm"] = 1e-25;       // local plasma vacuum density

    variables["Delta_x"] = 1e-10;
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;
    variables["A"] = 1e-10;
    variables["omega"] = 1e-14;
    variables["sigma"] = 1e16;
    variables["v_r"] = 1e3;
    variables["t"] = 0.0;

    variables["mu_0"] = 4 * variables["pi"] * 1e-7;
    variables["lambda_I"] = 1.0;
    variables["omega_i"] = 1e-8;
    variables["t_n"] = 0.0;
    variables["F_RZ"] = 0.01;
    variables["k_SF"] = 1e-10;
    variables["omega_spin"] = 1e-12;
    variables["I_dipole"] = 1e17;
    variables["A_dipole"] = 1e10;
    variables["H_aether"] = 1e-7;
    variables["k_4"] = 1.0;
    variables["delta_rho_over_rho"] = 1e-5;
    variables["f_TRZ"] = 0.01;
}

void NGC346UQFFModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
    if (name == "Delta_x") variables["Delta_p"] = variables["hbar"] / value;
    else if (name == "M_visible" || name == "M_DM") {
        variables["M"] = variables["M_visible"] + variables["M_DM"];
        variables["M0"] = variables["M"];
    }
}
void NGC346UQFFModule::addToVariable(const std::string& name, double delta) { variables[name] += delta; }
void NGC346UQFFModule::subtractFromVariable(const std::string& name, double delta) { variables[name] -= delta; }

double NGC346UQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// NOVEL: Um = q * v_rad * B — universal magnetism from charged particle blueshift
double NGC346UQFFModule::computeUm(double t) {
    return variables["q"] * variables["v_rad"] * variables["B"];
}

double NGC346UQFFModule::computeUg1(double t) {
    double mu_dipole = variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"];
    return mu_dipole * variables["B"];
}
double NGC346UQFFModule::computeUg2(double t) {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}

// NOVEL: Ug3 magnetic strings = G*M/r^2 * (rho_gas/rho_vac_UA)
double NGC346UQFFModule::computeUg3strings(double r) {
    return (variables["G"] * variables["M"] / (r * r)) * (variables["rho_gas"] / variables["rho_vac_UA"]);
}

double NGC346UQFFModule::computeUg4(double t) {
    return variables["k_4"] * 1e35 * std::exp(-0.0005 * t);
}

// NOVEL: Ui = lambda_I * (rho_vac_UA/rho_plasm) * omega_i * cos(pi*t_n)
double NGC346UQFFModule::computeUi(double t) {
    return variables["lambda_I"] * (variables["rho_vac_UA"] / variables["rho_plasm"]) * variables["omega_i"] * std::cos(variables["pi"] * variables["t_n"]) * (1 + variables["F_RZ"]);
}

double NGC346UQFFModule::computeEcore(double r) {
    return computeUg3strings(r) + computeUi(variables["t"]) * variables["rho_gas"];
}
double NGC346UQFFModule::computeTempCore(double r) {
    return 1.424e7 * computeUg3strings(r) * variables["rho_vac_SCm"];
}

double NGC346UQFFModule::computeFenv(double t) {
    double F_SF = variables["k_SF"] * variables["SFR"] / 1.989e30;
    double F_collapse = 0.05 * (t / (1e6 * variables["year_to_s"]));
    return F_SF + F_collapse;
}
double NGC346UQFFModule::computeMsfFactor(double t) { return variables["SFR"] * t / variables["M0"]; }
double NGC346UQFFModule::computeRt(double t) { return variables["r"] + variables["v_r"] * t; }

double NGC346UQFFModule::computePsiIntegral(double r, double t) {
    std::complex<double> psi(variables["A"] * std::exp(-r * r / (2 * variables["sigma"] * variables["sigma"])) * std::exp(std::complex<double>(0, -variables["omega"] * t)));
    return std::norm(psi);
}
double NGC346UQFFModule::computeQuantumTerm(double t_Hubble_val, double r) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * computePsiIntegral(r, variables["t"]);
}
double NGC346UQFFModule::computeFluidTerm(double g_base) { return variables["rho_fluid"] * variables["V"] * g_base; }
double NGC346UQFFModule::computeDMTerm(double r) {
    return (variables["M_visible"] + variables["M_DM"]) * (variables["delta_rho_over_rho"] + 3 * variables["G"] * variables["M"] / (r * r * r));
}

double NGC346UQFFModule::computeUgSum(double r) {
    double Ug_base = variables["G"] * variables["M"] / (r * r);
    variables["Um"] = computeUm(variables["t"]);
    variables["Ug1"] = computeUg1(variables["t"]);
    variables["Ug2"] = computeUg2(variables["t"]);
    variables["Ug3str"] = computeUg3strings(r);
    variables["Ug4"] = computeUg4(variables["t"]);
    return Ug_base + variables["Um"] + variables["Ug1"] + variables["Ug2"] + variables["Ug3str"] + variables["Ug4"];
}

double NGC346UQFFModule::computeG(double t, double r) {
    variables["t"] = t;
    double m_factor = 1.0 + computeMsfFactor(t);
    double expansion = 1.0 + computeHtz(variables["z"]) * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double f_env = computeFenv(t);
    double g_base = (variables["G"] * variables["M"] * m_factor / (r * r)) * expansion * sc_correction * (1.0 + f_env);
    double ug_sum = computeUgSum(r) - variables["G"] * variables["M"] / (r * r);
    double lambda_term = variables["Lambda"] * variables["c"] * variables["c"] / 3.0;
    double ui_term = computeUi(t);
    double quantum_term = computeQuantumTerm(variables["t_Hubble"], r);
    double fluid_term = computeFluidTerm(g_base);
    double dm_term = computeDMTerm(r);
    return g_base + ug_sum + lambda_term + ui_term + quantum_term + fluid_term + dm_term;
}

std::string NGC346UQFFModule::getEquationText() {
    return "g_NGC346(r,t) = (G M(t)/r^2)(1+H(t,z))(1-B/B_crit)(1+F_env) + Sum Ugi + Ui + Lambda c^2/3\n"
           "NOVEL Um = q*v_rad*B (blueshift magnetism); Ug3_strings = G*M/r^2*(rho_gas/rho_vac_UA)\n"
           "Ui = lambda_I*(rho_UA/rho_plasm)*omega_i*cos(pi*t_n); M=1000 Msun (SMC HII region)\n"
           "T_core = 1.424e7 * Ug3_strings * rho_vac; E_core = Ug3 + Ui*rho_gas";
}
void NGC346UQFFModule::printVariables() {
    std::cout << "NGC 346 (SMC Nebula) Variables:\n";
    for (const auto& pair : variables) std::cout << pair.first << " = " << std::scientific << pair.second << "\n";
}
