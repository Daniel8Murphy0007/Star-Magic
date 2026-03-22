// RedSpiderUQFFModule.cpp
// Frequency-domain UQFF for Red Spider Nebula (NGC 6537).
// NOVEL: a = (f_DPM + f_react + f_super + f_fluid + f_quantum + f_Aether + f_THz) * lambda_Planck / (2*pi)
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "RedSpiderUQFFModule.h"

RedSpiderUQFFModule::RedSpiderUQFFModule() {
    variables["G"] = 6.6743e-11;
    variables["c"] = 3e8;
    variables["hbar"] = 1.0546e-34;
    variables["pi"] = 3.141592653589793;
    variables["lambda_Planck"] = 1.616e-35;   // Planck length (m)
    variables["H0"] = 70.0;
    variables["Mpc_to_m"] = 3.086e22;
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["year_to_s"] = 3.156e7;

    // Red Spider Nebula (NGC 6537) parameters
    variables["r"] = 7.1e15;        // effective radius (m, ~0.23 pc)
    variables["M"] = 0.6 * 1.989e30;  // central WD mass ~0.6 Msun
    variables["rho_lobe"] = 1e-22;  // lobe density (kg/m^3)
    variables["rho_fil"] = 1e-20;   // filament density (kg/m^3)
    variables["v_exp"] = 3e5;        // bipolar outflow velocity (m/s, ~300 km/s)
    variables["z"] = 0.0015;
    variables["t_age"] = 1900 * variables["year_to_s"];   // nebula age

    // Frequency channels (Hz)
    variables["f_super"] = 1.411e16;    // SuperFreq (SGR1745/SgrA* derived)
    variables["f_fluid"] = 1.269e-14;   // Navier-Stokes frequency
    variables["f_quantum"] = 1.445e-17; // quantum gravity coupling frequency
    variables["f_Aether"] = 1.576e-35;  // Aether coupling frequency
    variables["f_react"] = 1e10;        // charge-reactivity base frequency
    variables["f_DPM"] = 1e12;          // di-pseudo-monopole frequency
    variables["f_THz"] = 1e12;          // THz contribution

    variables["B"] = 1e-3;          // field ~1 mT in nebula
    variables["q"] = 1.602e-19;
    variables["v_ion"] = 1e5;        // ion velocity for DPM calc
    variables["rho_vac_plasm"] = 1e-30;  // plasma vacuum density
    variables["rho_fluid"] = variables["rho_lobe"];
    variables["t"] = 0.0;
    variables["mu_0"] = 4 * variables["pi"] * 1e-7;
}

void RedSpiderUQFFModule::updateVariable(const std::string& name, double value) { variables[name] = value; }
void RedSpiderUQFFModule::addToVariable(const std::string& name, double delta) { variables[name] += delta; }
void RedSpiderUQFFModule::subtractFromVariable(const std::string& name, double delta) { variables[name] -= delta; }

double RedSpiderUQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// DPM frequency: f_DPM = f_DPM_base * rho_vac_plasm / c
double RedSpiderUQFFModule::computeFdpm() {
    return variables["f_DPM"] * variables["rho_vac_plasm"] / variables["c"];
}
// Reactivity frequency
double RedSpiderUQFFModule::computeFreactivity() { return variables["f_react"]; }
// SuperFreq (calibrated from SGR1745/SgrA*)
double RedSpiderUQFFModule::computeFsuper() { return variables["f_super"]; }
// Fluid frequency (Navier-Stokes coupling)
double RedSpiderUQFFModule::computeFfluid() { return variables["f_fluid"]; }
// Quantum frequency
double RedSpiderUQFFModule::computeFquantum() { return variables["f_quantum"]; }
// Aether frequency
double RedSpiderUQFFModule::computeFAether() { return variables["f_Aether"]; }
// THz contribution
double RedSpiderUQFFModule::computeFTHz() { return variables["f_THz"]; }

// Environmental factor (for expansion correction)
double RedSpiderUQFFModule::computeFenv(double t) {
    double F_wind = variables["rho_lobe"] * variables["v_exp"] * variables["v_exp"];
    double F_age = 0.01 * (t / variables["t_age"]);
    return F_wind + F_age;
}

double RedSpiderUQFFModule::computeFtotal() {
    return computeFdpm() + computeFreactivity() + computeFsuper() + computeFfluid() + computeFquantum() + computeFAether() + computeFTHz();
}

// NOVEL frequency-domain acceleration: a = f_total * lambda_Planck / (2*pi)
double RedSpiderUQFFModule::computeG(double t, double r) {
    variables["t"] = t;
    double f_total = computeFtotal();
    double a_freq = f_total * variables["lambda_Planck"] / (2.0 * variables["pi"]);

    // Expansion correction via H(z)
    double H_z = computeHtz(variables["z"]);
    double expansion = 1.0 + H_z * t;

    // Environmental scaling
    double f_env_scale = 1.0 + computeFenv(t);

    return a_freq * expansion * f_env_scale;
}

std::string RedSpiderUQFFModule::getEquationText() {
    return "a_RedSpider(t) = (f_DPM + f_react + f_super + f_fluid + f_quantum + f_Aether + f_THz) * lambda_P / (2*pi)\n"
           "NOVEL: all UQFF terms as frequency contributions to acceleration [m/s^2]\n"
           "f_DPM = f_DPM_base * rho_vac_plasm/c; f_super=1.411e16 Hz (SGR1745 calibrated)\n"
           "NGC 6537 Red Spider Nebula; r=7.1e15 m, v_exp=300 km/s, t_age=1900 yr";
}
void RedSpiderUQFFModule::printVariables() {
    std::cout << "Red Spider Nebula (NGC 6537) Variables:\n";
    for (const auto& pair : variables) std::cout << pair.first << " = " << std::scientific << pair.second << "\n";
}
