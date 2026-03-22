// NGC1316UQFFModule.cpp
// Implementation of MUGE for NGC 1316 (Dust Bunnies) merger dynamics.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "NGC1316UQFFModule.h"
#include <complex>

// Constructor: NGC 1316-specific values
NGC1316UQFFModule::NGC1316UQFFModule() {
    // Universal constants
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

    // NGC 1316 parameters
    variables["M_visible"] = 3.5e11 * M_sun_val;
    variables["M_DM"] = 1.5e11 * M_sun_val;
    variables["M"] = variables["M_visible"] + variables["M_DM"];
    variables["M0"] = variables["M"];
    variables["M_spiral"] = 1e10 * M_sun_val;
    variables["d_spiral"] = 50e3 * kpc_val;
    variables["M_BH"] = 1e8 * M_sun_val;
    variables["M_cluster"] = 1e6 * M_sun_val;
    variables["r"] = 46e3 * kpc_val;
    variables["z"] = 0.005;
    variables["tau_merge"] = 1e9 * variables["year_to_s"];
    variables["t"] = 2e9 * variables["year_to_s"];

    // Dynamics
    variables["rho_dust"] = 1e-21;
    variables["V"] = 1e51;
    variables["B"] = 1e-4;
    variables["B_crit"] = 1e11;
    variables["Delta_x"] = 1e-10;
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;

    // Wave/oscillatory for dust lanes
    variables["A"] = 1e-10;
    variables["k"] = 1e20;
    variables["omega"] = 1e-16;
    variables["x"] = 0.0;
    variables["v"] = 1e3;
    variables["sigma"] = 2e3 * kpc_val;

    // Ug subterms & Ui
    variables["Ug1"] = 0.0;
    variables["Ug2"] = 0.0;
    variables["Ug3"] = 0.0;
    variables["Ug4"] = 0.0;
    variables["Ui"] = 0.0;
    variables["mu_0"] = 4 * variables["pi"] * 1e-7;
    variables["rho_vac_SCm"] = 7.09e-37;
    variables["rho_vac_UA"] = 7.09e-36;
    variables["lambda_I"] = 1.0;
    variables["omega_i"] = 1e-8;
    variables["t_n"] = 0.0;
    variables["F_RZ"] = 0.01;
    variables["k_4"] = 1.0;
    variables["k_cluster"] = 1e-12;
    variables["omega_spin"] = 1e-3;
    variables["I_dipole"] = 1e20;
    variables["A_dipole"] = 1e15;
    variables["H_aether"] = 1e-5;
    variables["delta_rho_over_rho"] = 1e-5;

    // Scales
    variables["scale_macro"] = 1e-12;
    variables["f_TRZ"] = 0.1;
    variables["f_sc"] = 1.0;
    variables["v_r"] = 1e3;
    variables["rho"] = variables["rho_dust"];
}

void NGC1316UQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding." << std::endl;
        variables[name] = value;
    }
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "M_visible" || name == "M_DM") {
        variables["M"] = variables["M_visible"] + variables["M_DM"];
        variables["M0"] = variables["M"];
    }
}

void NGC1316UQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) variables[name] += delta;
    else variables[name] = delta;
}
void NGC1316UQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

double NGC1316UQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

double NGC1316UQFFModule::computeMmerge(double t) {
    return 1e10 * 1.989e30 * std::exp(-t / variables["tau_merge"]);
}

double NGC1316UQFFModule::computeRt(double t) {
    return variables["r"] + variables["v_r"] * t;
}

double NGC1316UQFFModule::computeFenv(double t) {
    double F_tidal = (variables["G"] * variables["M_spiral"]) / (variables["d_spiral"] * variables["d_spiral"]);
    double F_cluster = variables["k_cluster"] * (variables["M_cluster"] / 1.989e30);
    return F_tidal + F_cluster;
}

double NGC1316UQFFModule::computeUg1(double t) {
    double mu_dipole = variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"];
    return mu_dipole * variables["B"];
}

double NGC1316UQFFModule::computeUg2(double t) {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}

double NGC1316UQFFModule::computeUg3prime(double t) {
    return (variables["G"] * variables["M_spiral"]) / (variables["d_spiral"] * variables["d_spiral"]);
}

double NGC1316UQFFModule::computeUg4(double t) {
    double E_react = 1e46 * std::exp(-0.0005 * t);
    return variables["k_4"] * E_react;
}

double NGC1316UQFFModule::computeUi(double t) {
    return variables["lambda_I"] * (variables["rho_vac_SCm"] / variables["rho_vac_UA"]) * variables["omega_i"] * std::cos(variables["pi"] * variables["t_n"]) * (1 + variables["F_RZ"]);
}

double NGC1316UQFFModule::computePsiIntegral(double r, double t) {
    double A = variables["A"];
    double m = 2.0;
    double omega = variables["omega"];
    double sigma = variables["sigma"];
    std::complex<double> psi_dust(A * std::exp(-r*r / (2 * sigma * sigma)) * std::exp(std::complex<double>(0, m * 0 - omega * t)));
    return std::norm(psi_dust);
}

double NGC1316UQFFModule::computeQuantumTerm(double t_Hubble_val, double r) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double psi_int = computePsiIntegral(r, variables["t"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * psi_int;
}

double NGC1316UQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_dust"] * variables["V"] * g_base;
}

double NGC1316UQFFModule::computeDMTerm(double r) {
    double pert = variables["delta_rho_over_rho"];
    double curv = 3 * variables["G"] * variables["M"] / (r * r * r);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

double NGC1316UQFFModule::computeUgSum(double r) {
    double Ug_base = (variables["G"] * variables["M"]) / (r * r);
    variables["Ug1"] = computeUg1(variables["t"]);
    variables["Ug2"] = computeUg2(variables["t"]);
    variables["Ug3"] = computeUg3prime(variables["t"]);
    variables["Ug4"] = computeUg4(variables["t"]);
    return Ug_base + variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

double NGC1316UQFFModule::computeG(double t, double r) {
    variables["t"] = t;
    double m_merge = computeMmerge(t);
    double m_factor = 1.0 + m_merge / variables["M0"];
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

std::string NGC1316UQFFModule::getEquationText() {
    return "g_NGC1316(r, t) = (G * M(t) / r(t)^2) * (1 + H(t, z)) * (1 - B(t) / B_crit) * (1 + F_env(t)) + "
           "(U_g1 + U_g2 + U_g3' + U_g4) + U_i + (Lambda * c^2 / 3) + "
           "(hbar / sqrt(Dx * Dp)) * integral(psi H psi dV) * (2pi / t_Hubble) + "
           "rho_dust * V * g + (M_visible + M_DM) * (drho/rho + 3 G M / r^3)\n"
           "Where: M(t) = M * (1 + M_merge(t)); M_merge(t) = 1e10 Msun * exp(-t/tau); r(t) = r0 + v_r t;\n"
           "H(t, z) = H0 * sqrt(Om (1+z)^3 + OL); F_env(t) = F_tidal + F_cluster;\n"
           "F_tidal = G M_spiral / d^2; F_cluster = k_cluster * M_cluster; U_g1 = mu_dipole * B;\n"
           "U_g2 = B_super^2 / (2 mu0); U_g3' = G M_spiral / d^2; U_g4 = k4 * E_react(t);\n"
           "U_i = lambda_I * (rho_SCm/rho_UA) * omega_i * cos(pi t_n) * (1 + F_RZ);\n"
           "Adaptations: Hubble ACS 2003; M=5e11 Msun; rho_dust=1e-21 kg/m^3. g ~2e37 m/s^2 at t=2 Gyr.";
}

void NGC1316UQFFModule::printVariables() {
    std::cout << "NGC 1316 Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}
