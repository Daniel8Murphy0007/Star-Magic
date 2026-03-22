// UQFFCompressedResonanceModule.cpp
// Multi-system compressed and resonance UQFF equation dispatcher.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "UQFFCompressedResonanceModule.h"
#include <complex>

UQFFCompressedResonanceModule::UQFFCompressedResonanceModule()
    : current_system("Guide"), mode("compressed") {
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
    variables["M_sun"] = 1.989e30;
    variables["kpc"] = 3.086e19;

    // General defaults
    variables["M"] = 1e41;
    variables["M0"] = variables["M"];
    variables["SFR"] = 6e19;
    variables["r"] = 1e20;
    variables["z"] = 0.005;
    variables["M_visible"] = 0.7 * variables["M"];
    variables["M_DM"] = 0.3 * variables["M"];
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
    variables["omega"] = 1e15;
    variables["x"] = 0.0;
    variables["v"] = 1e3;
    variables["Ug1"] = 0.0; variables["Ug2"] = 0.0;
    variables["Ug3"] = 0.0; variables["Ug4"] = 0.0;
    variables["scale_macro"] = 1e-12;
    variables["f_TRZ"] = 0.1;
    variables["f_sc"] = 1.0;
    variables["delta_rho"] = 1e-5 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];
    variables["F_wind"] = 0.0;
}

void UQFFCompressedResonanceModule::setSystem(const std::string& sys_name) {
    current_system = sys_name;
    double Msun = variables["M_sun"];
    double kpc = variables["kpc"];
    double yr_s = variables["year_to_s"];
    if (sys_name == "YoungStars") {
        variables["M"] = 1000 * Msun; variables["r"] = 3e17;
        variables["SFR"] = 0.1 * Msun / yr_s;
        variables["rho_fluid"] = 1e-20; variables["B"] = 1e-6; variables["z"] = 0.0006;
    } else if (sys_name == "Eagle") {
        variables["M"] = 1e4 * Msun; variables["r"] = 2e17;
        variables["SFR"] = 0.5 * Msun / yr_s;
        variables["rho_fluid"] = 1e-21; variables["B"] = 3e-5; variables["z"] = 0.002;
    } else if (sys_name == "BigBang") {
        variables["rho_fluid"] = 8e-27; variables["r"] = 1e26;
        variables["z"] = 1100; variables["SFR"] = 0;
        variables["M"] = 1e53; variables["B"] = 1e-10;
        variables["t"] = 13.8e9 * yr_s;
    } else if (sys_name == "M51") {
        variables["M"] = 1.6e11 * Msun; variables["r"] = 23e3 * kpc;
        variables["SFR"] = 2 * Msun / yr_s;
        variables["rho_fluid"] = 1e-21; variables["B"] = 1e-5; variables["z"] = 0.005;
    } else if (sys_name == "NGC1316") {
        variables["M"] = 5e11 * Msun; variables["r"] = 23e3 * kpc;
        variables["SFR"] = 0.1 * Msun / yr_s;
        variables["rho_fluid"] = 1e-22; variables["B"] = 1e-5; variables["z"] = 0.006;
    } else if (sys_name == "V838Mon") {
        variables["M"] = 8 * Msun; variables["r"] = 2e13;
        variables["SFR"] = 0;
        variables["rho_fluid"] = 1e-22; variables["B"] = 1e-6; variables["z"] = 0.005;
    } else if (sys_name == "NGC1300") {
        variables["M"] = 1e11 * Msun; variables["r"] = 12e3 * kpc;
        variables["SFR"] = 1 * Msun / yr_s;
        variables["rho_fluid"] = 1e-21; variables["B"] = 1e-5; variables["z"] = 0.005;
    } else {  // Guide
        variables["M"] = Msun; variables["r"] = 1e11;
        variables["SFR"] = 1e-10 * Msun / yr_s;
        variables["rho_fluid"] = 1e-20; variables["B"] = 1e-5; variables["z"] = 0;
    }
    variables["M_visible"] = 0.7 * variables["M"];
    variables["M_DM"] = 0.3 * variables["M"];
    variables["M0"] = variables["M"];
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["rho"] = variables["rho_fluid"];
    variables["delta_rho"] = 1e-5 * variables["rho_fluid"];
}

void UQFFCompressedResonanceModule::setMode(const std::string& m) { mode = m; }

void UQFFCompressedResonanceModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
    if (name == "M") {
        variables["M_visible"] = 0.7 * value;
        variables["M_DM"] = 0.3 * value;
        variables["M0"] = value;
    }
    if (name == "Delta_x") variables["Delta_p"] = variables["hbar"] / value;
}
void UQFFCompressedResonanceModule::addToVariable(const std::string& name, double delta) {
    updateVariable(name, variables[name] + delta);
}
void UQFFCompressedResonanceModule::subtractFromVariable(const std::string& name, double delta) {
    updateVariable(name, variables[name] - delta);
}

double UQFFCompressedResonanceModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}
double UQFFCompressedResonanceModule::computeFenv(double t) { return 0.1; }
double UQFFCompressedResonanceModule::computeUgSum() { return 1e-10; }
double UQFFCompressedResonanceModule::computePsiTotal(double t) {
    return variables["q"] * variables["v"] * variables["B"]
         + 2 * variables["A"] * std::cos(variables["k"] * variables["x"] + variables["omega"] * t);
}
double UQFFCompressedResonanceModule::computeResonanceTerm(double t) {
    if (mode != "resonance") return 0.0;
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega"] * t)));
    double g_base = (variables["G"] * variables["M"]) / (variables["r"] * variables["r"]);
    return (2 * variables["pi"] / 13.8) * exp_term.real() * g_base;
}
double UQFFCompressedResonanceModule::computeQuantumTerm(double t_Hubble_val) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double psi = computePsiTotal(variables["t"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * psi;
}
double UQFFCompressedResonanceModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}
double UQFFCompressedResonanceModule::computeDMTerm() {
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3 * variables["G"] * variables["M"] / std::pow(variables["r"], 3);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}
double UQFFCompressedResonanceModule::computeMsfFactor(double t) {
    return variables["SFR"] * t / variables["M0"];
}

double UQFFCompressedResonanceModule::computeG(double t, double r_in) {
    if (r_in > 0) variables["r"] = r_in;
    variables["t"] = t;
    if (current_system == "BigBang") variables["r"] = variables["c"] * t;
    if (current_system == "V838Mon") {
        double rho_d = variables["rho_fluid"] * std::exp(-1.0 * (variables["G"] * variables["M"] / (variables["r"] * variables["r"])));
        return (600000 * 3.826e26) / (4 * variables["pi"] * variables["r"] * variables["r"]) * 1e-12 * rho_d;
    }
    double Hz = computeHtz(variables["z"]);
    double expansion = 1 + Hz * t;
    double sc = 1 - variables["B"] / variables["B_crit"];
    double msf = computeMsfFactor(t);
    double mfact = 1 + msf;
    double fenv = computeFenv(t);
    double g_base = (variables["G"] * variables["M"] * mfact / (variables["r"] * variables["r"])) * expansion * sc * (1 + fenv);
    double ugsum = computeUgSum();
    double lambda_t = variables["Lambda"] * variables["c"] * variables["c"] / 3;
    double qterm = computeQuantumTerm(variables["t_Hubble"]);
    double fterm = computeFluidTerm(g_base);
    double dmterm = computeDMTerm();
    double res_term = computeResonanceTerm(t);
    return g_base + ugsum + lambda_t + qterm + fterm + dmterm + res_term;
}

std::string UQFFCompressedResonanceModule::getEquationText() {
    std::string eq = "g_UQFF(r,t) = (G M(t)/r^2) (1 + H(t,z)) (1 - B/B_crit) (1 + F_env) + "
                     "Sum Ug_i + Lambda c^2/3 + (hbar/sqrt(Dx Dp)) integral(psi H psi dV) (2pi/t_Hubble) + "
                     "rho V g + (M_vis + M_DM)(drho/rho + 3GM/r^3)";
    if (mode == "resonance")
        eq += " + 2 A cos(kx + omega t) g_base + (2pi/13.8) Re[A exp(i(kx - omega t))] g_base";
    eq += "\nM(t)=M(1 + SFR t / M0); System: " + current_system + " Mode: " + mode;
    return eq;
}

void UQFFCompressedResonanceModule::printVariables() {
    std::cout << "System: " << current_system << " Mode: " << mode << "\nVariables:\n";
    for (auto& p : variables) std::cout << p.first << " = " << std::scientific << p.second << "\n";
}
