// SMBHBinaryUQFFModule.cpp
// Frequency-domain UQFF for SMBH Binary gravitational wave source (LISA target).
// NOVEL: a = f_total * lambda_Planck / (2*pi), same method as RedSpider but for binary BH.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "SMBHBinaryUQFFModule.h"

SMBHBinaryUQFFModule::SMBHBinaryUQFFModule() {
    variables["G"] = 6.6743e-11;
    variables["c"] = 3e8;
    variables["hbar"] = 1.0546e-34;
    variables["pi"] = 3.141592653589793;
    variables["lambda_Planck"] = 1.616e-35;
    variables["H0"] = 70.0;
    variables["Mpc_to_m"] = 3.086e22;
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["year_to_s"] = 3.156e7;
    double Msun = 1.989e30;

    // SMBH Binary parameters
    variables["M1"] = 4e6 * Msun;
    variables["M2"] = 2e6 * Msun;
    variables["M"] = variables["M1"] + variables["M2"];
    variables["mu"] = (variables["M1"] * variables["M2"]) / variables["M"];   // reduced mass
    variables["r_init"] = 0.1 * 9.461e15;   // 0.1 light-year in meters
    variables["t_coal"] = 1.555e7;           // coalescence time (s)
    variables["z"] = 0.1;
    variables["rho"] = 1e-20;
    variables["B"] = 1e-5;
    variables["B_crit"] = 1e11;
    variables["mu_0"] = 4 * variables["pi"] * 1e-7;

    // Frequency channels (Hz)
    variables["f_super"] = 1.411e16;     // SuperFreq
    variables["f_fluid"] = 5.070e-8;     // Navier-Stokes (scaled for SMBH)
    variables["f_quantum"] = 1.445e-17;  // quantum coupling
    variables["f_Aether"] = 1.576e-35;   // Aether coupling
    variables["f_react"] = 1e8;          // charge reactivity (SMBH-scale)
    variables["f_DPM"] = 1e10;           // DPM frequency
    variables["f_THz"] = 1e10;           // THz coupling
    variables["rho_vac_plasm"] = 1e-30;

    variables["H_aether"] = 1e-6;
    variables["rho_vac_SCm"] = 7.09e-37;
    variables["rho_vac_UA"] = 7.09e-36;
    variables["omega_spin"] = 1e-12;
    variables["I_dipole"] = 1e25;
    variables["A_dipole"] = 1e18;
    variables["k_4"] = 1.0;
    variables["t"] = 0.0;
}

void SMBHBinaryUQFFModule::updateVariable(const std::string& name, double value) { variables[name] = value; }
void SMBHBinaryUQFFModule::addToVariable(const std::string& name, double delta) { variables[name] += delta; }
void SMBHBinaryUQFFModule::subtractFromVariable(const std::string& name, double delta) { variables[name] -= delta; }

double SMBHBinaryUQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// Orbital separation decreasing toward coalescence (Peters formula simplified)
double SMBHBinaryUQFFModule::computeRsep(double t) {
    double r0 = variables["r_init"];
    double t_c = variables["t_coal"];
    if (t >= t_c) return 1e8;  // Schwarzschild radius scale at merger
    return r0 * std::pow(1.0 - t / t_c, 0.25);
}

double SMBHBinaryUQFFModule::computeFsuper() { return variables["f_super"]; }
double SMBHBinaryUQFFModule::computeFfluid() { return variables["f_fluid"]; }
double SMBHBinaryUQFFModule::computeFquantum() { return variables["f_quantum"]; }
double SMBHBinaryUQFFModule::computeFAether() { return variables["f_Aether"]; }

double SMBHBinaryUQFFModule::computeFUg1(double r) {
    double mu_dipole = variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"];
    return mu_dipole * variables["B"] / r;
}
double SMBHBinaryUQFFModule::computeFUg2() {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}
double SMBHBinaryUQFFModule::computeFUg3(double r) {
    return variables["G"] * variables["M2"] / (r * r);
}
double SMBHBinaryUQFFModule::computeFUg4(double t) {
    return variables["k_4"] * 1e40 * std::exp(-0.0005 * t);
}

double SMBHBinaryUQFFModule::computeFtotal(double t, double r) {
    return computeFsuper() + computeFfluid() + computeFquantum() + computeFAether()
         + computeFUg1(r) + computeFUg2() + computeFUg3(r) + computeFUg4(t)
         + variables["f_DPM"] * variables["rho_vac_plasm"] / variables["c"]
         + variables["f_THz"];
}

double SMBHBinaryUQFFModule::computeG(double t, double r) {
    variables["t"] = t;
    double r_sep = computeRsep(t);
    double f_total = computeFtotal(t, r_sep);
    double a_freq = f_total * variables["lambda_Planck"] / (2.0 * variables["pi"]);
    double expansion = 1.0 + computeHtz(variables["z"]) * t;
    return a_freq * expansion;
}

double SMBHBinaryUQFFModule::computeCoalescenceTime() { return variables["t_coal"]; }

std::string SMBHBinaryUQFFModule::getEquationText() {
    return "a_SMBHBinary(t) = f_total * lambda_P / (2*pi); r_sep(t) = r0*(1 - t/t_coal)^(1/4)\n"
           "f_total = f_super + f_fluid + f_quantum + f_Aether + f_Ug1 + f_Ug2 + f_Ug3 + f_Ug4 + f_DPM + f_THz\n"
           "M1=4e6 Msun, M2=2e6 Msun, t_coal=1.555e7 s; LISA gravitational wave source target";
}
void SMBHBinaryUQFFModule::printVariables() {
    std::cout << "SMBH Binary Variables:\n";
    for (const auto& pair : variables) std::cout << pair.first << " = " << std::scientific << pair.second << "\n";
}
