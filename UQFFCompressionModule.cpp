// UQFFCompressionModule.cpp
// Implementation of UQFF Compression Cycle — 12-channel F_env, 3 pre-defined systems.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "UQFFCompressionModule.h"
#include <complex>

UQFFCompressionModule::UQFFCompressionModule() : m_system("Generic") {
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
    double Msun = 1.989e30;

    // Default: Generic system
    variables["M"] = 1e11 * Msun;
    variables["M0"] = variables["M"];
    variables["SFR"] = 1.0 * Msun / variables["year_to_s"];
    variables["r"] = 10 * 3.086e19;
    variables["z"] = 0.0;
    variables["rho_fluid"] = 1e-21;
    variables["V"] = 1e51;
    variables["B"] = 1e-5;
    variables["B_crit"] = 1e11;

    // F_env 12 channels (default weights)
    variables["F_wind"] = 1e-3;
    variables["F_rad"] = 1e-4;
    variables["F_SN"] = 1e-4;
    variables["F_BH"] = 0.0;
    variables["F_erode"] = 1e-5;
    variables["F_lensing"] = 1e-6;
    variables["F_mag"] = 1e-5;
    variables["F_decay"] = 0.0;
    variables["F_coll"] = 0.0;
    variables["F_evo"] = 1e-6;
    variables["F_merge"] = 0.0;
    variables["F_sf"] = 1e-4;

    // External mass for Ug3'
    variables["M_ext"] = 4e6 * Msun;    // default: companion mass
    variables["r_ext"] = 50 * 3.086e19;  // companion distance

    // Compression wave parameters
    variables["A_wave"] = 1e-10;
    variables["k_wave"] = 1e20;
    variables["omega_wave"] = 1e-14;
    variables["x_wave"] = 0.0;
    variables["v"] = 1e5;            // ion velocity for qvB
    variables["A_comp"] = 1e-12;     // compression amplitude
    variables["T_baktun"] = 13.8 * variables["year_to_s"];  // 2pi/13.8 resonance

    variables["sigma_disp"] = 1e19;
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
    variables["t"] = 0.0;
    variables["Delta_x"] = 1e-10;
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;
}

void UQFFCompressionModule::setSystem(const std::string& system_name) {
    m_system = system_name;
    double Msun = 1.989e30;
    if (system_name == "Magnetar") {
        variables["M"] = 2 * Msun;
        variables["r"] = 1e4;
        variables["B"] = 1e11;
        variables["M_ext"] = 4e6 * Msun;   // Sgr A* gravitational influence
        variables["r_ext"] = 8e3 * 3.086e19;
        variables["F_BH"] = 0.0;
        variables["F_erode"] = 0.0;
    } else if (system_name == "SagittariusA") {
        variables["M"] = 4e6 * Msun;
        variables["r"] = 1e10;
        variables["B"] = 1e-3;
        variables["F_BH"] = 1e42 / (variables["G"] * variables["M"]);  // normalized
        variables["M_ext"] = 1e10 * Msun;   // GC stellar mass
        variables["r_ext"] = 1e18;
    } else if (system_name == "Pillars") {
        variables["M"] = 1e3 * Msun;
        variables["r"] = 3e16;        // ~1 pc
        variables["B"] = 1e-8;
        variables["F_erode"] = 0.05;  // initially, grows with t
        variables["rho_fluid"] = 1e-20;
        variables["F_wind"] = variables["rho_fluid"] * (20e3) * (20e3);  // v_wind=20 km/s
    }
}

void UQFFCompressionModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
    if (name == "Delta_x") variables["Delta_p"] = variables["hbar"] / value;
}
void UQFFCompressionModule::addToVariable(const std::string& name, double delta) { variables[name] += delta; }
void UQFFCompressionModule::subtractFromVariable(const std::string& name, double delta) { variables[name] -= delta; }

double UQFFCompressionModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}
double UQFFCompressionModule::computeMsfFactor(double t) { return variables["SFR"] * t / variables["M0"]; }

// 12-channel F_env: includes time-dependent erosion for Pillars
double UQFFCompressionModule::computeFenv12(double t) {
    double t_yr = t / variables["year_to_s"];
    double F_erode_t = variables["F_erode"] * (1.0 + t / (1e6 * variables["year_to_s"]));  // grows with time
    return variables["F_wind"] + variables["F_rad"] + variables["F_SN"]
         + variables["F_BH"] + F_erode_t + variables["F_lensing"]
         + variables["F_mag"] + variables["F_decay"] + variables["F_coll"]
         + variables["F_evo"] + variables["F_merge"] + variables["F_sf"];
}

double UQFFCompressionModule::computeUg3prime() {
    return variables["G"] * variables["M_ext"] / (variables["r_ext"] * variables["r_ext"]);
}

// NOVEL: psi_total = qvB + 2A*cos(kx)*cos(omega*t) + (2pi/13.8)*A*Re[exp(i(kx-omega*t))]
double UQFFCompressionModule::computePsiTotal(double t, double x) {
    double qvB = variables["q"] * variables["v"] * variables["B"];
    double comp_wave = 2.0 * variables["A_wave"] * std::cos(variables["k_wave"] * x) * std::cos(variables["omega_wave"] * t);
    double kx_ot = variables["k_wave"] * x - variables["omega_wave"] * t;
    double resonance = (2.0 * variables["pi"] / 13.8) * variables["A_comp"] * std::cos(kx_ot);
    return qvB + comp_wave + resonance;
}

double UQFFCompressionModule::computeUg1(double t) {
    return variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"] * variables["B"];
}
double UQFFCompressionModule::computeUg2(double t) {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}
double UQFFCompressionModule::computeUg4(double t) {
    return variables["k_4"] * 1e40 * std::exp(-0.0005 * t);
}
double UQFFCompressionModule::computeUi(double t) {
    return variables["lambda_I"] * (variables["rho_vac_SCm"] / variables["rho_vac_UA"]) * variables["omega_i"] * std::cos(variables["pi"] * variables["t_n"]) * (1 + variables["F_RZ"]);
}
double UQFFCompressionModule::computeQuantumTerm(double t_Hubble_val, double r) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double psi2 = variables["A_wave"] * variables["A_wave"] * std::exp(-r * r / (2 * variables["sigma_disp"] * variables["sigma_disp"]));
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * psi2;
}
double UQFFCompressionModule::computeFluidTerm(double g_base) { return variables["rho_fluid"] * variables["V"] * g_base; }
double UQFFCompressionModule::computeDMTerm(double r) {
    return variables["M"] * (variables["delta_rho_over_rho"] + 3 * variables["G"] * variables["M"] / (r * r * r));
}
double UQFFCompressionModule::computeUgSum(double r) {
    double Ug_base = variables["G"] * variables["M"] / (r * r);
    double psi = computePsiTotal(variables["t"], variables["x_wave"]);
    double u3p = computeUg3prime();
    variables["Ug1"] = computeUg1(variables["t"]);
    variables["Ug2"] = computeUg2(variables["t"]);
    variables["Ug4"] = computeUg4(variables["t"]);
    return Ug_base + psi + u3p + variables["Ug1"] + variables["Ug2"] + variables["Ug4"];
}

double UQFFCompressionModule::computeG(double t, double r) {
    variables["t"] = t;
    double m_factor = 1.0 + computeMsfFactor(t);
    double expansion = 1.0 + computeHtz(variables["z"]) * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double f_env = computeFenv12(t);
    double g_base = (variables["G"] * variables["M"] * m_factor / (r * r)) * expansion * sc_correction * (1.0 + f_env);
    double ug_sum = computeUgSum(r) - variables["G"] * variables["M"] / (r * r);
    double lambda_term = variables["Lambda"] * variables["c"] * variables["c"] / 3.0;
    double ui_term = computeUi(t);
    double qt = computeQuantumTerm(variables["t_Hubble"], r);
    double ft = computeFluidTerm(g_base);
    double dt = computeDMTerm(r);
    return g_base + ug_sum + lambda_term + ui_term + qt + ft + dt;
}

std::string UQFFCompressionModule::getEquationText() {
    return "g_compression(r,t) = (G M(t)/r^2)(1+H_z)(1-B/B_crit)(1+F_env12) + Psi + Ug3' + Sum Ugi\n"
           "Psi_total = qvB + 2A*cos(kx)*cos(omega*t) + (2pi/13.8)*A*Re[exp(i(kx-omega*t))]\n"
           "F_env12 = F_wind + F_rad + F_SN + F_BH + F_erode(t) + F_lensing + F_mag + F_decay + F_coll + F_evo + F_merge + F_sf\n"
           "Systems: Magnetar (M=2Msun,r=1e4,B=1e11), SagittariusA (M=4e6Msun), Pillars (erosion)";
}
void UQFFCompressionModule::printVariables() {
    std::cout << "UQFF Compression [" << m_system << "] Variables:\n";
    for (const auto& pair : variables) std::cout << pair.first << " = " << std::scientific << pair.second << "\n";
}
