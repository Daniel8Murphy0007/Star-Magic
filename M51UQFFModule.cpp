// M51UQFFModule.cpp
// Implementation of MUGE for M51 (Whirlpool Galaxy) — tidal NGC 5195, spiral density wave.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "M51UQFFModule.h"
#include <complex>

M51UQFFModule::M51UQFFModule() {
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

    // M51 Whirlpool Galaxy parameters
    variables["M_visible"] = 1.2e11 * Msun;
    variables["M_DM"] = 4e10 * Msun;
    variables["M"] = variables["M_visible"] + variables["M_DM"];
    variables["M0"] = variables["M"];
    variables["SFR"] = 1.0 * Msun / variables["year_to_s"];
    variables["r"] = 23.58 * variables["kpc_to_m"];
    variables["z"] = 0.002;
    variables["rho_fluid"] = 1e-21;
    variables["V"] = 1e52;
    variables["B"] = 1e-5;
    variables["B_crit"] = 1e11;

    // Central SMBH and NGC 5195 companion
    variables["M_BH"] = 1e6 * Msun;
    variables["M_NGC5195"] = 1e10 * Msun;
    variables["d_NGC5195"] = 50 * variables["kpc_to_m"];     // companion separation
    variables["tau_merge"] = 5e8 * variables["year_to_s"];   // merger timescale

    // Spiral density wave parameters (m=2 arm mode)
    variables["A_spiral"] = 1e-10;          // wave amplitude
    variables["sigma_spiral"] = 10 * variables["kpc_to_m"];  // spiral scale
    variables["m_arm"] = 2.0;               // 2-arm spiral
    variables["omega_spiral"] = 1e-16;      // pattern speed
    variables["theta"] = 0.0;              // azimuthal angle

    variables["Delta_x"] = 1e-10;
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;
    variables["v_r"] = 1e3;               // radial drift velocity
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
    variables["E_react"] = 1e42;          // reaction energy from M_BH (J)
    variables["delta_rho_over_rho"] = 1e-5;
    variables["f_TRZ"] = 0.01;
}

void M51UQFFModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
    if (name == "Delta_x") variables["Delta_p"] = variables["hbar"] / value;
    else if (name == "M_visible" || name == "M_DM") {
        variables["M"] = variables["M_visible"] + variables["M_DM"];
        variables["M0"] = variables["M"];
    }
}
void M51UQFFModule::addToVariable(const std::string& name, double delta) { variables[name] += delta; }
void M51UQFFModule::subtractFromVariable(const std::string& name, double delta) { variables[name] -= delta; }

double M51UQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// NGC 5195 merger mass evolution  
double M51UQFFModule::computeMmerge(double t) {
    return variables["M_NGC5195"] * std::exp(-t / variables["tau_merge"]);
}

// Radial expansion
double M51UQFFModule::computeRt(double t) { return variables["r"] + variables["v_r"] * t; }

double M51UQFFModule::computeFenv(double t) {
    double F_tidal = variables["G"] * computeMmerge(t) / (variables["d_NGC5195"] * variables["d_NGC5195"]);
    double F_SF = variables["k_SF"] * variables["SFR"] / 1.989e30;
    return F_tidal + F_SF;
}

// Ug1: magnetic dipole I*A*omega*B
double M51UQFFModule::computeUg1(double t) {
    return variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"] * variables["B"];
}
// Ug2: superconductor field B^2/(2*mu_0)
double M51UQFFModule::computeUg2(double t) {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}
// Ug3': tidal companion G*M_NGC5195(t)/d^2
double M51UQFFModule::computeUg3prime(double t) {
    return variables["G"] * computeMmerge(t) / (variables["d_NGC5195"] * variables["d_NGC5195"]);
}
// Ug4: BH reaction energy decay k_4*E_react*exp(-0.0005*t)
double M51UQFFModule::computeUg4(double t) {
    return variables["k_4"] * variables["E_react"] * std::exp(-0.0005 * t);
}
// Ui: inertia with vacuum density ratio + F_RZ
double M51UQFFModule::computeUi(double t) {
    return variables["lambda_I"] * (variables["rho_vac_SCm"] / variables["rho_vac_UA"]) * variables["omega_i"] * std::cos(variables["pi"] * variables["t_n"]) * (1 + variables["F_RZ"]);
}

// Spiral density wave: psi = A*exp(-r^2/2*sigma^2)*exp(i*(m*theta - omega*t))
double M51UQFFModule::computePsiSpiral(double r, double theta, double t) {
    double envelope = variables["A_spiral"] * std::exp(-r * r / (2 * variables["sigma_spiral"] * variables["sigma_spiral"]));
    std::complex<double> phase = std::exp(std::complex<double>(0, variables["m_arm"] * theta - variables["omega_spiral"] * t));
    std::complex<double> psi = envelope * phase;
    return std::norm(psi);
}
double M51UQFFModule::computePsiIntegral(double r, double t) {
    return computePsiSpiral(r, variables["theta"], t);
}
double M51UQFFModule::computeQuantumTerm(double t_Hubble_val, double r) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * computePsiIntegral(r, variables["t"]);
}
double M51UQFFModule::computeFluidTerm(double g_base) { return variables["rho_fluid"] * variables["V"] * g_base; }
double M51UQFFModule::computeDMTerm(double r) {
    return (variables["M_visible"] + variables["M_DM"]) * (variables["delta_rho_over_rho"] + 3 * variables["G"] * variables["M"] / (r * r * r));
}
double M51UQFFModule::computeMsfFactor(double t) { return variables["SFR"] * t / variables["M0"]; }

double M51UQFFModule::computeUgSum(double r) {
    double Ug_base = variables["G"] * variables["M"] / (r * r);
    variables["Ug1"] = computeUg1(variables["t"]);
    variables["Ug2"] = computeUg2(variables["t"]);
    variables["Ug3p"] = computeUg3prime(variables["t"]);
    variables["Ug4"] = computeUg4(variables["t"]);
    return Ug_base + variables["Ug1"] + variables["Ug2"] + variables["Ug3p"] + variables["Ug4"];
}

double M51UQFFModule::computeG(double t, double r) {
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

std::string M51UQFFModule::getEquationText() {
    return "g_M51(r,t) = (G M(t)/r^2)(1+H(t,z))(1-B/B_crit)(1+F_env) + Sum Ugi + Ui + Lambda c^2/3\n"
           "Ug1 = I*A*omega*B (magnetic dipole); Ug3' = G*M_NGC5195(t)/d^2 (tidal)\n"
           "Ug4 = k_4*E_react*exp(-0.0005*t) [BH reaction decay]; r(t) = r0 + v_r*t\n"
           "psi_spiral = A*exp(-r^2/2sigma^2)*exp(i*(m*theta - omega*t));  M_NGC5195=1e10 Msun";
}
void M51UQFFModule::printVariables() {
    std::cout << "M51 Whirlpool Galaxy Variables:\n";
    for (const auto& pair : variables) std::cout << pair.first << " = " << std::scientific << pair.second << "\n";
}
