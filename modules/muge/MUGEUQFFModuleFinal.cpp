// MUGEUQFFModuleFinal.cpp
// Final MUGE/UQFF for 7 canonical SOURCE4 systems with 10-term resonance acceleration.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#include "MUGEUQFFModuleFinal.h"

MUGEUQFFModuleFinal::MUGEUQFFModuleFinal(SystemTypeFinal sys) : current_system(sys) {
    double Msun = 1.989e30;

    // Universal constants
    variables["G"]            = 6.6743e-11;
    variables["c"]            = 3e8;
    variables["hbar"]         = 1.0546e-34;
    variables["Lambda"]       = 1.1e-52;
    variables["q"]            = 1.602e-19;
    variables["pi"]           = 3.141592653589793;
    variables["t_Hubble"]     = 13.8e9 * 3.156e7;
    variables["year_to_s"]    = 3.156e7;
    variables["H0"]           = 70.0;
    variables["Mpc_to_m"]     = 3.086e22;
    variables["Omega_m"]      = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["z"]            = 0.0;
    variables["B_crit"]       = 1e11;
    variables["D_p"]          = 4.4e26;

    // Generic defaults
    variables["M"]            = Msun;
    variables["r"]            = 1e10;
    variables["M_visible"]    = 0.15 * variables["M"];
    variables["M_DM"]         = 0.85 * variables["M"];
    variables["rho_fluid"]    = 1e-20;
    variables["V"]            = 1e3;
    variables["B"]            = 1e-5;
    variables["Delta_x"]      = 1e-10;
    variables["Delta_p"]      = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;
    variables["A_wave"]       = 1e-10;
    variables["k_wave"]       = 1e20;
    variables["omega_wave"]   = 1e15;
    variables["x_wave"]       = 0.0;
    variables["delta_rho"]    = 1e-21;
    variables["rho"]          = variables["rho_fluid"];

    // Ug components
    variables["Ug1"]       = 0.0;
    variables["Ug2"]       = 0.0;
    variables["Ug3_prime"] = 0.0;
    variables["Ug4"]       = 0.0;

    // F_env
    const char* fkeys[] = {
        "F_wind","F_erode","F_merge","F_SN","F_rad","F_fil",
        "F_BH","F_dust","F_ring","F_mag","F_tech","F_shell","F_cosmo",
        "F_torque","F_shock","QG_term","DM_term","GW_term"
    };
    for (const auto& k : fkeys) variables[k] = 0.0;
    variables["SC_m"] = 1.0;

    // Resonance acceleration params (Doc 42.a)
    variables["fTHz"]      = 1e12;        // Hz
    variables["Evac_neb"]  = 7.09e-36;    // J/m^3
    variables["Evac_ISM"]  = 7.09e-37;    // J/m^3
    variables["vexp"]      = 1e3;         // m/s
    variables["aDPM"]      = 1e-10;       // m/s^2 base
    variables["DeltaEvac"] = 6.381e-36;   // J/m^3
    variables["Fsuper"]    = 6.287e-19;
    variables["UA_SCm"]    = 10.0;        // [UA']:[SCm] ratio
    variables["omega_i"]   = 1e-8;        // rad/s
    variables["fTRZ"]      = 0.1;
    variables["k4"]        = 1.0;
    variables["Ereact"]    = 1e-20;       // J
    variables["freact"]    = 1e10;        // Hz
    variables["fquantum"]  = 1.445e-17;
    variables["fAether"]   = 1.576e-35;
    variables["ffluid"]    = 1.269e-14;
    variables["fosc"]      = 4.57e14;     // Hz (optical)
    variables["fexp"]      = 1e-15;

    // Resonance base vars
    variables["A_res"]     = 1.0;
    variables["f_res"]     = 1e15;
    variables["k_A"]       = 1.0;
    variables["Z_nuc"]     = 1.0;
    variables["A_nuc"]     = 1.0;
    variables["A_H"]       = 1.0;
    variables["delta_pair"]= 0.0;
    variables["E_bind"]    = 13.6;
    variables["h_planck"]  = 4.1356677e-15;
    variables["k_res"]     = 1.0;
    variables["A1"]        = 1.0;
    variables["A2"]        = 1.0;
    variables["f_dp"]      = 1e3;
    variables["phi_dp"]    = 0.0;
    variables["k_0"]       = 1.0;
    variables["N_nuc"]     = 1.0;
    variables["Z_magic"]   = 0.0;
    variables["N_magic"]   = 0.0;
    variables["S_shell"]   = 0.0;

    variables["M_total"]   = 1e53;
    variables["r_c"]       = 1e27;
    variables["k_cosm"]    = 0.0;
    variables["t"]         = 1e6 * variables["year_to_s"];

    setSystem(sys);
}

void MUGEUQFFModuleFinal::setSystem(SystemTypeFinal sys) {
    current_system = sys;
    double Msun = 1.989e30;

    switch (sys) {
        case SystemTypeFinal::MAGNETAR_SGR1745:
            variables["M"]       = 1.5 * Msun;
            variables["r"]       = 1e4;          // 10 km
            variables["z"]       = 0.0009;
            variables["M_ext"]   = 4.1e6 * Msun;
            variables["r_ext"]   = 2.84e15;      // ~0.3 ly
            variables["B"]       = 1e10;          // 10^14 G = 10^10 T
            variables["rho_fluid"]= 1e-15;
            variables["V"]       = 4.189e12;
            variables["t"]       = 3.799e10;      // ~1203 yr
            variables["Evac_neb"]= 7.09e-36;
            variables["Evac_ISM"]= 7.09e-37;
            variables["DeltaEvac"]= 6.381e-36;
            variables["vexp"]    = 1e3;
            break;
        case SystemTypeFinal::SGR_A:
            variables["M"]       = 4.1e6 * Msun;
            variables["r"]       = 1.2e10;
            variables["z"]       = 0.0009;
            variables["B"]       = 1e2;
            variables["t"]       = 1e10 * variables["year_to_s"];
            break;
        case SystemTypeFinal::TAPESTRY_STARBIRTH:
            variables["M"]       = 1e5 * Msun;
            variables["r"]       = 1e20;
            variables["z"]       = 0.01;
            variables["SFR"]     = 100.0 * Msun / variables["year_to_s"];
            variables["F_wind"]  = 1e-10;
            break;
        case SystemTypeFinal::WESTERLUND2:
            variables["M"]       = 1e4 * Msun;
            variables["r"]       = 2e19;
            variables["z"]       = 0.0036;
            variables["F_rad"]   = -1e-11;
            break;
        case SystemTypeFinal::PILLARS_CREATION:
            variables["M"]       = 1e3 * Msun;
            variables["r"]       = 1e19;
            variables["z"]       = 0.002;
            variables["F_erode"] = -1e-12;
            break;
        case SystemTypeFinal::RINGS_RELATIVITY:
            variables["M"]       = 1e11 * Msun;
            variables["r"]       = 1e21;
            variables["z"]       = 0.05;
            variables["F_ring"]  = 1e-9;
            break;
        case SystemTypeFinal::STUDENTS_GUIDE:
            variables["M"]       = 1e12 * Msun;
            variables["r"]       = 1e21;
            variables["z"]       = 0.0;
            break;
        case SystemTypeFinal::GENERIC:
        default:
            break;
    }

    // Recompute deps
    variables["Delta_p"]   = variables["hbar"] / variables["Delta_x"];
    variables["M_visible"] = 0.15 * variables["M"];
    variables["M_DM"]      = 0.85 * variables["M"];
    if (variables.count("M_ext") && variables.count("r_ext")) {
        variables["Ug3_prime"] = (variables["G"] * variables["M_ext"]) /
                                  (variables["r_ext"] * variables["r_ext"]);
    }
    variables["S_shell"] = 0.1 * (variables["Z_magic"] + variables["N_magic"]);
    variables["f_res"]   = (variables["E_bind"] / variables["h_planck"]) *
                            (variables["A_H"] / variables["A_nuc"]) *
                            (1.0 + variables["S_shell"]);
}

void MUGEUQFFModuleFinal::updateVariable(const std::string& name, double value) {
    variables[name] = value;
    if (name == "Delta_x") variables["Delta_p"] = variables["hbar"] / value;
    else if (name == "M") {
        variables["M_visible"] = 0.15 * value;
        variables["M_DM"]      = 0.85 * value;
    } else if (name == "QG_term" || name == "DM_term" || name == "GW_term") {
        variables["F_cosmo"] = variables["QG_term"] + variables["DM_term"] + variables["GW_term"];
    }
}

void MUGEUQFFModuleFinal::addToVariable(const std::string& name, double delta) {
    if (variables.count(name)) variables[name] += delta;
}

void MUGEUQFFModuleFinal::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

double MUGEUQFFModuleFinal::computeHtz(double z) {
    double Hz_kms = variables["H0"] * std::sqrt(
        variables["Omega_m"] * std::pow(1.0 + z, 3.0) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

double MUGEUQFFModuleFinal::computeFenv(double t) {
    double fenv = 0.0;
    const std::string fkeys[] = {
        "F_wind","F_erode","F_merge","F_SN","F_rad","F_fil",
        "F_BH","F_dust","F_ring","F_mag","F_tech","F_shell","F_cosmo",
        "F_torque","F_shock"
    };
    for (const auto& k : fkeys) {
        if (variables.count(k)) fenv += variables.at(k);
    }
    if (variables.count("SFR")) {
        double t_yr = t / variables["year_to_s"];
        fenv += (variables["SFR"] * t_yr) / variables["M"];
    }
    return fenv;
}

std::complex<double> MUGEUQFFModuleFinal::computePsiTotal(double t) {
    double cos_t = 2.0 * variables["A_wave"] *
                   std::cos(variables["k_wave"] * variables["x_wave"]) *
                   std::cos(variables["omega_wave"] * t);
    double exp_r = variables["A_wave"] *
                   std::cos(variables["k_wave"] * variables["x_wave"] -
                             variables["omega_wave"] * t);
    return std::complex<double>(cos_t + exp_r, 0.0);
}

double MUGEUQFFModuleFinal::computeQuantumTerm(double t_Hubble_val) {
    double unc   = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    auto   psi   = computePsiTotal(variables["t"]);
    double integ = variables["integral_psi"] * std::norm(psi);
    return (variables["hbar"] / unc) * integ * (2.0 * variables["pi"] / t_Hubble_val);
}

double MUGEUQFFModuleFinal::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

double MUGEUQFFModuleFinal::computeDMTerm() {
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3.0 * variables["G"] * variables["M"] /
                   std::pow(variables["r"], 3.0);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

double MUGEUQFFModuleFinal::computeUgSum() {
    variables["Ug1"] = (variables["G"] * variables["M"]) /
                        (variables["r"] * variables["r"]);
    variables["Ug4"] = variables["Ug1"];
    return variables["Ug1"] + variables["Ug2"] +
           variables["Ug3_prime"] + variables["Ug4"];
}

double MUGEUQFFModuleFinal::computeUQFF(double t) {
    variables["t"] = t;
    double Htz     = computeHtz(variables["z"]);
    double sc_corr = 1.0 - (variables["B"] / variables["B_crit"]);
    double fenv    = computeFenv(t);

    double g_base  = (variables["G"] * variables["M"] * (1.0 + fenv) /
                       (variables["r"] * variables["r"])) *
                      (1.0 + Htz * t) * sc_corr;

    return g_base + computeUgSum() +
           variables["Lambda"] * std::pow(variables["c"], 2.0) / 3.0 +
           computeQuantumTerm(variables["t_Hubble"]) +
           computeFluidTerm(g_base) +
           computeDMTerm();
}

// 10-term resonance acceleration (Doc 42.a core contribution)
double MUGEUQFFModuleFinal::computeResonanceAcc(double t) {
    double fT   = variables["fTHz"];
    double En   = variables["Evac_neb"];
    double Ei   = variables["Evac_ISM"];
    double ve   = variables["vexp"];
    double aD   = variables["aDPM"];
    double c    = variables["c"];
    double DEv  = variables["DeltaEvac"];
    double Fs   = variables["Fsuper"];
    double UASm = variables["UA_SCm"];
    double oi   = variables["omega_i"];
    double fRZ  = variables["fTRZ"];
    double k4   = variables["k4"];
    double Er   = variables["Ereact"];
    double fr   = variables["freact"];
    double fq   = variables["fquantum"];
    double fA   = variables["fAether"];
    double ff   = variables["ffluid"];
    double fo   = variables["fosc"];
    double fex  = variables["fexp"];
    double V    = variables["V"];
    double pi   = variables["pi"];

    double aTHz      = fT * En * ve * aD / (Ei * c);
    double avac_diff = DEv * ve * ve * aD / (En * c * c);
    double aSuperF   = Fs * fT * aD / (En * c);
    double aAethR    = UASm * oi * fT * aD * (1.0 + fRZ);
    double Ug4i      = k4 * Er * fr * aD / (En * c);
    double aQuantF   = fq * En * aD / (Ei * c);
    double aAethF    = fA * En * aD / (Ei * c);
    double aFluidF   = ff * En * V / (Ei * c);
    double OscTerm   = fo * std::sin(2.0 * pi * fo * t) * 1e-20;
    double aExpF     = fex * En * aD / (Ei * c);

    return aTHz + avac_diff + aSuperF + aAethR + Ug4i +
           aQuantF + aAethF + aFluidF + OscTerm + aExpF;
}

double MUGEUQFFModuleFinal::computeHres(double t) {
    // Resonance base
    double A_res = variables["k_A"] * variables["Z_nuc"] *
                   (variables["A_nuc"] / variables["A_H"]) *
                   (1.0 + variables["delta_pair"]);
    double sin_t = A_res * std::sin(2.0 * variables["pi"] * variables["f_res"] * t);
    double U_dp  = variables["k_res"] *
                   (variables["A1"] * variables["A2"] /
                    std::pow(variables["f_dp"], 2.0)) *
                   std::cos(variables["phi_dp"]);
    double k_nuc = variables["k_0"] *
                   (variables["N_nuc"] / variables["Z_nuc"]) *
                   (1.0 + variables["delta_pair"]);
    double res_base = sin_t + U_dp * variables["SC_m"] * k_nuc + variables["S_shell"];

    double fenv     = computeFenv(t);
    double res_acc  = computeResonanceAcc(t);
    return res_base + fenv * variables["SC_m"] + res_acc;
}

double MUGEUQFFModuleFinal::computeDuniverse() {
    double Hz    = computeHtz(variables["z"]);
    double H0_si = (variables["H0"] * 1e3) / variables["Mpc_to_m"];
    double hbar  = variables["hbar"];
    double Dx    = variables["Delta_x"];
    double Dp    = variables["Delta_p"];

    return 2.0 * variables["D_p"] *
           (1.0 + Hz * variables["t"]) *
           (1.0 + variables["Lambda"] * variables["c"] * variables["c"] / (3.0 * H0_si * H0_si)) *
           (1.0 + hbar / (std::sqrt(Dx * Dp) * variables["G"] * variables["M_total"])) *
           (1.0 + variables["k_cosm"] * variables["r_c"] * variables["r_c"]);
}

std::string MUGEUQFFModuleFinal::getEquationText() {
    return
        "MUGE Final (Doc 42.a) — 7 Canonical SOURCE4 Systems:\n\n"
        "Compressed UQFF:\n"
        "  g_UQFF = (GM/r^2)(1+H(t,z))(1-B/B_crit)(1+F_env)\n"
        "         + Ug_sum + Lambda*c^2/3 + quantum + fluid + DM\n\n"
        "Resonance H_res = res_base + F_env*SC_m + res_acc\n"
        "  res_acc = aTHz + avac_diff + aSuperFreq + aAetherRes + Ug4i\n"
        "          + aQuantumFreq + aAetherFreq + aFluidFreq + OscTerm + aExpFreq\n"
        "  aTHz    = fTHz * Evac_neb * vexp * aDPM / (Evac_ISM * c)\n"
        "  avac_diff = DeltaEvac * vexp^2 * aDPM / (Evac_neb * c^2)\n"
        "  OscTerm = fosc * sin(2pi*fosc*t) * 1e-20";
}

std::string MUGEUQFFModuleFinal::getSolutions(double t) {
    double g_uqff   = computeUQFF(t);
    double res_acc  = computeResonanceAcc(t);
    double h_res    = computeHres(t);
    double D_uni    = computeDuniverse();
    double fenv     = computeFenv(t);
    double ug       = computeUgSum();
    double lam_t    = variables["Lambda"] * std::pow(variables["c"], 2.0) / 3.0;
    double g_base   = (variables["G"] * variables["M"] /
                        (variables["r"] * variables["r"]));
    double Htz_t    = computeHtz(variables["z"]) * t;
    double sc_adj   = 1.0 - variables["B"] / variables["B_crit"];

    std::ostringstream ss;
    ss << std::scientific << std::setprecision(4);
    ss << "MUGEUQFFModuleFinal Solutions t=" << t << " s (system=" << static_cast<int>(current_system) << "):\n\n";
    ss << "=== Compressed UQFF ===\n";
    ss << "  g_base    = " << g_base   << " m/s^2\n";
    ss << "  H(t,z)    = " << Htz_t    << "\n";
    ss << "  SC adj    = " << sc_adj   << "\n";
    ss << "  F_env(t)  = " << fenv     << "\n";
    ss << "  Ug_sum    = " << ug       << " m/s^2\n";
    ss << "  Lambda_t  = " << lam_t    << "\n";
    ss << "  g_UQFF    = " << g_uqff   << " m/s^2\n\n";
    ss << "=== Resonance H_res ===\n";
    ss << "  res_acc   = " << res_acc  << " m/s^2 (10-term sum)\n";
    ss << "  H_res(t)  = " << h_res    << "\n\n";
    ss << "D_universe  = " << D_uni    << " m\n";
    return ss.str();
}

void MUGEUQFFModuleFinal::printVariables() {
    std::cout << "MUGEUQFFModuleFinal variables (system=" << static_cast<int>(current_system) << "):\n";
    for (const auto& p : variables) {
        std::cout << "  " << p.first << " = " << p.second << "\n";
    }
}
