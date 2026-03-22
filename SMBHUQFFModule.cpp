// SMBHUQFFModule.cpp
// UQFF M-sigma relation validation with oscillating magnetic moment and quantum state scaling.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "SMBHUQFFModule.h"
#include <complex>

SMBHUQFFModule::SMBHUQFFModule() {
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

    // SMBH M-sigma parameters
    variables["M_bh"] = 1e9 * Msun;         // default SMBH mass ~1e9 Msun
    variables["sigma"] = 200e3;              // bulge velocity dispersion (m/s)
    variables["R_bulge"] = 1 * variables["kpc_to_m"];  // bulge radius
    variables["t"] = 4.543e9 * variables["year_to_s"]; // default ~Milky Way age
    variables["z"] = 0.0;
    variables["B"] = 1e-9;                  // intergalactic field
    variables["B_crit"] = 1e11;
    variables["mu_0"] = 4 * variables["pi"] * 1e-7;

    // Oscillating magnetic moment: mu_j(t) = (1e3 + 0.4*sin(omega_c*t)) * mu_base_factor
    variables["omega_c"] = 2 * variables["pi"] / 3.96e8;  // ~12.5 years period
    variables["mu_base_factor"] = 3.38e20;      // base magnetic moment scaling
    variables["gamma"] = 0.00005 / variables["year_to_s"];  // decay rate (day^-1 converted)

    // Golden ratio and quantum state scaling
    variables["phi"] = (1.0 + std::sqrt(5.0)) / 2.0;   // golden ratio ~1.618
    variables["n_states"] = 6.0;                          // number of quantum states
    variables["f_feedback"] = 0.063;                     // metal retention calibration

    // Galactic coupling
    variables["omega_s"] = 2 * variables["pi"] / (200e6 * variables["year_to_s"]);  // galactic rotation
    variables["k_galactic"] = 1e-15;                   // galactic coupling constant

    // Vacuum densities
    variables["rho_vac_SCm"] = 7.09e-37;
    variables["rho_vac_UA"] = 7.09e-36;

    variables["f_TRZ"] = 0.01;
    variables["rho"] = 1e-23;          // SMBH environment density
    variables["V"] = 1e55;
    variables["k_SF"] = 1e-10;
}

void SMBHUQFFModule::updateVariable(const std::string& name, double value) { variables[name] = value; }
void SMBHUQFFModule::addToVariable(const std::string& name, double delta) { variables[name] += delta; }
void SMBHUQFFModule::subtractFromVariable(const std::string& name, double delta) { variables[name] -= delta; }

double SMBHUQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// NOVEL: mu_j(t) = (1e3 + 0.4*sin(omega_c*t)) * mu_base_factor
double SMBHUQFFModule::computeMuJ(double t) {
    return (1e3 + 0.4 * std::sin(variables["omega_c"] * t)) * variables["mu_base_factor"];
}

// NOVEL: delta_n = phi*(2*pi)^(n/6)
double SMBHUQFFModule::computeDeltaN(int n) {
    return variables["phi"] * std::pow(2.0 * variables["pi"], n / 6.0);
}

// NOVEL: rho_vac state coupling with double-exponential
double SMBHUQFFModule::computeRhoVacState(int n, double t) {
    double t_yr = t / variables["year_to_s"];
    double SS_q = 0.57;  // [SSq] calibrated constant
    double inner_exp = std::exp(-variables["pi"] - t_yr);
    return variables["rho_vac_UA"] * std::pow(variables["rho_vac_SCm"] / variables["rho_vac_UA"], n) * std::exp(-std::exp(-variables["pi"] - t_yr));
}

// Um = mu_j(t) * B (magnetic moment coupling)
double SMBHUQFFModule::computeUm(double t) {
    return computeMuJ(t) * variables["B"];
}

double SMBHUQFFModule::computeUg1(double t) {
    double mu_j = computeMuJ(t);
    double I_eff = mu_j / (variables["sigma"] * variables["R_bulge"]);
    return I_eff * variables["B"];
}

double SMBHUQFFModule::computeOmegaSKgal() {
    return variables["omega_s"] * variables["k_galactic"];
}

double SMBHUQFFModule::computeRhoVacSCm(double t) {
    double delta_n_sum = 0.0;
    for (int n = 1; n <= 6; ++n) delta_n_sum += computeDeltaN(n);
    return variables["rho_vac_SCm"] * delta_n_sum * std::exp(-variables["gamma"] * t);
}

double SMBHUQFFModule::computeFenv(double t) {
    return variables["f_feedback"] * computeRhoVacSCm(t) / variables["rho_vac_UA"];
}

double SMBHUQFFModule::computeUgSum(double r, double t) {
    return computeUm(t) + computeUg1(t) + computeOmegaSKgal();
}

double SMBHUQFFModule::computeG(double t, double r) {
    variables["t"] = t;
    double g_newton = variables["G"] * variables["M_bh"] / (r * r);
    double f_env = computeFenv(t);
    double expansion = 1.0 + computeHtz(variables["z"]) * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double g_base = g_newton * expansion * sc_correction * (1.0 + f_env);
    double ug_total = computeUgSum(r, t);
    double lambda_term = variables["Lambda"] * variables["c"] * variables["c"] / 3.0;
    return g_base + ug_total + lambda_term;
}

// M-sigma relation: log(M_bh/Msun) ~ 8.45 + 4.38*log(sigma/200 km/s) — Tremaine+2002
double SMBHUQFFModule::computeMsigmaPredict(double sigma_km_s) {
    return 1.989e30 * std::pow(10.0, 8.45 + 4.38 * std::log10(sigma_km_s / 200.0));
}

std::string SMBHUQFFModule::getEquationText() {
    return "g_SMBH(t,r) = (G*M_bh/r^2)(1+H(t,z))(1-B/B_crit)(1+F_env) + Um + Ug1 + omega_s*k_gal + Lambda c^2/3\n"
           "NOVEL: mu_j(t) = (1e3 + 0.4*sin(omega_c*t)) * 3.38e20; omega_c = 2pi/3.96e8 rad/s\n"
           "delta_n = phi*(2*pi)^(n/6); rho_vac_state: double-exp coupling exp(-exp(-pi - t/yr))\n"
           "M-sigma: log(M_bh/Msun) = 8.45 + 4.38*log(sigma/200 km/s); f_feedback=0.063";
}
void SMBHUQFFModule::printVariables() {
    std::cout << "SMBH M-sigma Variables:\n";
    for (const auto& pair : variables) std::cout << pair.first << " = " << std::scientific << pair.second << "\n";
}
