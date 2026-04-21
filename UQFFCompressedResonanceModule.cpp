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
    variables["F_rad"] = 0.0;
    variables["F_SN"] = 0.0;
    variables["F_BH"] = 0.0;
    variables["F_eta"] = 0.0;
    variables["F_lensing"] = 0.0;
    variables["F_merge"] = 0.0;
    variables["F_decay"] = 0.0;
    variables["F_coll"] = 0.0;
    variables["F_evo"] = 0.0;
    variables["M_ext"] = 1e40;
    variables["r_ext"] = 5e20;
    variables["SC_m"] = 1.0;
    variables["A_res"] = 1e-3;
    variables["f_res"] = 1.0 / variables["t_Hubble"];
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
        variables["F_wind"] = 2e-3; variables["F_rad"] = 5e-4; variables["M_ext"] = 300 * Msun; variables["r_ext"] = 1e18;
    } else if (sys_name == "Eagle") {
        variables["M"] = 1e4 * Msun; variables["r"] = 2e17;
        variables["SFR"] = 0.5 * Msun / yr_s;
        variables["rho_fluid"] = 1e-21; variables["B"] = 3e-5; variables["z"] = 0.002;
        variables["F_wind"] = 5e-3; variables["F_rad"] = 1e-3; variables["F_evo"] = 2e-4; variables["M_ext"] = 1500 * Msun; variables["r_ext"] = 6e17;
    } else if (sys_name == "BigBang") {
        variables["rho_fluid"] = 8e-27; variables["r"] = 1e26;
        variables["z"] = 1100; variables["SFR"] = 0;
        variables["M"] = 1e53; variables["B"] = 1e-10;
        variables["t"] = 13.8e9 * yr_s;
        variables["F_rad"] = 1e-2; variables["F_decay"] = 2e-3; variables["M_ext"] = 0.0; variables["r_ext"] = 1.0;
    } else if (sys_name == "M51") {
        variables["M"] = 1.6e11 * Msun; variables["r"] = 23e3 * kpc;
        variables["SFR"] = 2 * Msun / yr_s;
        variables["rho_fluid"] = 1e-21; variables["B"] = 1e-5; variables["z"] = 0.005;
        variables["F_wind"] = 1e-3; variables["F_SN"] = 8e-4; variables["F_merge"] = 4e-4; variables["M_ext"] = 2e10 * Msun; variables["r_ext"] = 5e20;
    } else if (sys_name == "NGC1316") {
        variables["M"] = 5e11 * Msun; variables["r"] = 23e3 * kpc;
        variables["SFR"] = 0.1 * Msun / yr_s;
        variables["rho_fluid"] = 1e-22; variables["B"] = 1e-5; variables["z"] = 0.006;
        variables["F_SN"] = 5e-4; variables["F_merge"] = 7e-4; variables["M_ext"] = 6e10 * Msun; variables["r_ext"] = 4e20;
    } else if (sys_name == "V838Mon") {
        variables["M"] = 8 * Msun; variables["r"] = 2e13;
        variables["SFR"] = 0;
        variables["rho_fluid"] = 1e-22; variables["B"] = 1e-6; variables["z"] = 0.005;
        variables["F_rad"] = 3e-3; variables["F_eta"] = 2e-4; variables["M_ext"] = 2 * Msun; variables["r_ext"] = 1e14;
    } else if (sys_name == "NGC1300") {
        variables["M"] = 1e11 * Msun; variables["r"] = 12e3 * kpc;
        variables["SFR"] = 1 * Msun / yr_s;
        variables["rho_fluid"] = 1e-21; variables["B"] = 1e-5; variables["z"] = 0.005;
        variables["F_wind"] = 8e-4; variables["F_SN"] = 4e-4; variables["F_lensing"] = 1e-4; variables["M_ext"] = 8e9 * Msun; variables["r_ext"] = 3e20;
    } else {  // Guide
        variables["M"] = Msun; variables["r"] = 1e11;
        variables["SFR"] = 1e-10 * Msun / yr_s;
        variables["rho_fluid"] = 1e-20; variables["B"] = 1e-5; variables["z"] = 0;
        variables["F_wind"] = 1e-5; variables["F_rad"] = 1e-6; variables["M_ext"] = 0.0; variables["r_ext"] = 1.0;
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
double UQFFCompressedResonanceModule::computeFenv(double t) {
    double time_fraction = t / std::max(variables["t_Hubble"], 1.0);
    double wind_term = variables["F_wind"] * (1.0 + 0.25 * computeMsfFactor(t));
    double decay_term = variables["F_decay"] * std::exp(-time_fraction);
    double evo_term = variables["F_evo"] * (1.0 + time_fraction);
    return wind_term + variables["F_rad"] + variables["F_SN"] + variables["F_BH"]
         + variables["F_eta"] + variables["F_lensing"] + variables["F_merge"]
         + decay_term + variables["F_coll"] + evo_term;
}
double UQFFCompressedResonanceModule::computeUgSum() {
    double safe_r_ext = std::max(variables["r_ext"], 1.0);
    double ug3prime = variables["G"] * variables["M_ext"] / (safe_r_ext * safe_r_ext);
    double ug1 = variables["Ug1"] != 0.0 ? variables["Ug1"] : variables["q"] * variables["v"] * variables["B"];
    double ug2 = variables["Ug2"] != 0.0 ? variables["Ug2"] : 0.5 * variables["rho_fluid"] * variables["v"] * variables["v"];
    double ug4 = variables["Ug4"] != 0.0 ? variables["Ug4"] : variables["Lambda"] * variables["c"] * variables["c"];
    variables["Ug3"] = ug3prime;
    return ug1 + ug2 + ug3prime + ug4;
}
double UQFFCompressedResonanceModule::computePsiTotal(double t) {
    return variables["q"] * variables["v"] * variables["B"]
         + 2 * variables["A"] * std::cos(variables["k"] * variables["x"] + variables["omega"] * t);
}
double UQFFCompressedResonanceModule::computeResonanceTerm(double t) {
    if (mode != "resonance") return 0.0;
    double oscillatory = variables["A_res"] * std::sin(2.0 * variables["pi"] * variables["f_res"] * t);
    return oscillatory + computeFenv(t) * variables["SC_m"];
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
    std::string eq = "g_UQFF(r,t) = (G M(t)/r^2) * (1 + H(t,z) * t) * (1 - B/B_crit) * (1 + F_env(t)) + "
                     "(Ug1 + Ug2 + Ug3' + Ug4) + Lambda*c^2/3 + "
                     "(hbar/sqrt(Delta_x*Delta_p)) * integral(psi_total * H_op * psi_total dV) * (2*pi/t_Hubble) + "
                     "rho_fluid*V*g + (M_visible + M_DM) * (delta_rho/rho + 3*G*M/r^3)";
    if (mode == "resonance")
        eq += " + H_res, where H_res = A_res*sin(2*pi*f_res*t) + F_env(t)*SC_m";
    eq += "\nM(t)=M*(1 + SFR*t/M0); Ug3' = G*M_ext/r_ext^2; System: " + current_system + "; Mode: " + mode;
    return eq;
}

void UQFFCompressedResonanceModule::printVariables() {
    std::cout << "System: " << current_system << " Mode: " << mode << "\nVariables:\n";
    for (auto& p : variables) std::cout << p.first << " = " << std::scientific << p.second << "\n";
}
