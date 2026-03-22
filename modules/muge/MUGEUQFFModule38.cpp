// MUGEUQFFModule38.cpp
// MUGE/UQFF compressed gravity for 38 source documents, 14 system types.
// Extends the 29-system module with 6 new systems and new F_env terms.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#include "MUGEUQFFModule38.h"

MUGEUQFFModule38::MUGEUQFFModule38(SystemType38 sys) : current_system(sys) {
    double Msun = 1.989e30;

    // Universal constants
    variables["G"]            = 6.6743e-11;
    variables["c"]            = 3e8;
    variables["hbar"]         = 1.0546e-34;
    variables["Lambda"]       = 1.1e-52;
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
    variables["Ug1"]         = 0.0;
    variables["Ug2"]         = 0.0;
    variables["Ug3_prime"]   = 0.0;
    variables["Ug4"]         = 0.0;

    // F_env (13 original + 5 new: F_torque, F_shock, QG_term, DM_term, GW_term)
    const char* fkeys[] = {
        "F_wind","F_erode","F_merge","F_SN","F_rad","F_fil",
        "F_BH","F_dust","F_ring","F_mag","F_tech","F_shell","F_cosmo",
        "F_torque","F_shock","QG_term","DM_term","GW_term"
    };
    for (const auto& k : fkeys) variables[k] = 0.0;
    variables["SC_m"] = 1.0;

    // Resonance
    variables["A_res"]     = 1.0;
    variables["f_res"]     = 1e15;
    variables["k_A"]       = 1.0;
    variables["Z_nuc"]     = 1.0;
    variables["A_nuc"]     = 1.0;
    variables["A_H"]       = 1.0;
    variables["delta_pair"]= 0.0;
    variables["E_bind"]    = 13.6;
    variables["h_planck"]  = 4.1356677e-15;
    variables["U_dp"]      = 0.0;
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

    // Universe
    variables["M_total"]   = 1e53;
    variables["r_c"]       = 1e27;
    variables["k_cosm"]    = 0.0;

    variables["t"]         = 1e6 * variables["year_to_s"];

    setSystem(sys);
}

void MUGEUQFFModule38::recomputeFcosmo() {
    variables["F_cosmo"] = variables["QG_term"] + variables["DM_term"] + variables["GW_term"];
}

void MUGEUQFFModule38::setSystem(SystemType38 sys) {
    current_system = sys;
    double Msun = 1.989e30;

    switch (sys) {
        // --- Original 8 from Doc 41 ---
        case SystemType38::SOMBRERO_GALAXY:
            variables["M"]      = 1e11 * Msun;
            variables["r"]      = 1.5e20;
            variables["z"]      = 0.0023;
            variables["M_ext"]  = 1e9 * Msun;
            variables["r_ext"]  = 1e19;
            variables["F_dust"] = 1e-11;
            variables["F_BH"]   = -1e-10;
            break;
        case SystemType38::SATURN:
            variables["M"]      = 5.68e26;
            variables["r"]      = 6e7;
            variables["z"]      = 0.0;
            variables["F_ring"] = 1e-10;
            variables["F_wind"] = 1e-11;
            break;
        case SystemType38::M16_EAGLE:
            variables["M"]      = 8e3 * Msun;
            variables["r"]      = 3.3e19;
            variables["z"]      = 0.002;
            variables["SFR"]    = 1.0 * Msun / variables["year_to_s"];
            variables["F_erode"]= -1e-12;
            break;
        case SystemType38::CRAB_NEBULA:
            variables["M"]      = 1.4 * Msun;
            variables["r"]      = 2.1e19;
            variables["z"]      = 0.004;
            variables["F_wind"] = 1e-10;
            variables["F_mag"]  = 1e-9;
            variables["B"]      = 100.0;
            break;
        case SystemType38::HYDROGEN_ATOM:
            variables["M"]      = 1.6726e-27 + 9.109e-31;
            variables["r"]      = 5.29e-11;
            variables["z"]      = 0.0;
            break;
        case SystemType38::HYDROGEN_RESONANCE:
            variables["z"]      = 0.0;
            variables["f_res"]  = (variables["E_bind"] / variables["h_planck"]) *
                                   (variables["A_H"] / variables["A_nuc"]) *
                                   (1.0 + variables["S_shell"]);
            break;
        case SystemType38::UNIVERSE_DIAMETER:
            variables["M"]      = 1e53;
            variables["r"]      = 1e27;
            variables["z"]      = 1100.0;
            variables["M_total"]= 1e53;
            break;

        // --- 6 new from Doc 42 ---
        case SystemType38::LAGOON_NEBULA:
            variables["M"]      = 1e4 * Msun;
            variables["r"]      = 6.6e19;
            variables["z"]      = 0.002;
            variables["SFR"]    = 2.0 * Msun / variables["year_to_s"];
            variables["F_rad"]  = -1e-12;
            break;
        case SystemType38::SPIRALS_SN:
            variables["M"]       = 1e10 * Msun;
            variables["r"]       = 3e20;
            variables["z"]       = 0.01;
            variables["F_torque"]= 1e-11;
            variables["F_SN"]    = 1e-10;
            break;
        case SystemType38::NGC6302:
            variables["M"]      = 2e3 * Msun;
            variables["r"]      = 3.3e18;
            variables["z"]      = 0.003;
            variables["F_shock"]= 1e-11;
            break;
        case SystemType38::ORION_NEBULA:
            variables["M"]      = 2000.0 * Msun;
            variables["r"]      = 1.18e17;
            variables["z"]      = 0.00034;
            variables["F_wind"] = 1e-11;
            variables["F_rad"]  = -1e-12;
            break;
        case SystemType38::YOUNG_STARS_OUTFLOW:
            variables["M"]      = 5e3 * Msun;
            variables["r"]      = 2e19;
            variables["z"]      = 0.001;
            variables["F_wind"] = 1e-10;
            break;
        case SystemType38::EAGLE_NEBULA:
            variables["M"]      = 8e3 * Msun;
            variables["r"]      = 3.3e19;
            variables["z"]      = 0.002;
            variables["F_wind"] = 1e-11;
            variables["F_rad"]  = -1e-12;
            break;
        case SystemType38::GRAVITY_BIGBANG:
            variables["M"]      = 1e53;
            variables["r"]      = 1e27;
            variables["z"]      = 1100.0;
            variables["QG_term"]= 1e-50;
            variables["DM_term"]= 1e-10;
            variables["GW_term"]= 1e-20;
            recomputeFcosmo();
            break;
        case SystemType38::GENERIC:
        default:
            break;
    }

    // Recompute standard dependencies
    variables["Delta_p"]    = variables["hbar"] / variables["Delta_x"];
    variables["M_visible"]  = 0.15 * variables["M"];
    variables["M_DM"]       = 0.85 * variables["M"];
    if (variables.count("M_ext") && variables.count("r_ext")) {
        variables["Ug3_prime"] = (variables["G"] * variables["M_ext"]) /
                                  (variables["r_ext"] * variables["r_ext"]);
    }
    variables["S_shell"] = 0.1 * (variables["Z_magic"] + variables["N_magic"]);
    variables["f_res"]   = (variables["E_bind"] / variables["h_planck"]) *
                            (variables["A_H"] / variables["A_nuc"]) *
                            (1.0 + variables["S_shell"]);
}

void MUGEUQFFModule38::updateVariable(const std::string& name, double value) {
    variables[name] = value;
    // Auto-cascade F_cosmo when any of its components change
    if (name == "QG_term" || name == "DM_term" || name == "GW_term") {
        recomputeFcosmo();
    }
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "M") {
        variables["M_visible"] = 0.15 * value;
        variables["M_DM"]      = 0.85 * value;
    }
}

void MUGEUQFFModule38::addToVariable(const std::string& name, double delta) {
    if (variables.count(name)) variables[name] += delta;
    if (name == "QG_term" || name == "DM_term" || name == "GW_term") {
        recomputeFcosmo();
    }
}

void MUGEUQFFModule38::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

double MUGEUQFFModule38::computeHtz(double z) {
    double Hz_kms = variables["H0"] * std::sqrt(
        variables["Omega_m"] * std::pow(1.0 + z, 3.0) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

double MUGEUQFFModule38::computeFenv(double t) {
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

std::complex<double> MUGEUQFFModule38::computePsiTotal(double t) {
    double cos_t = 2.0 * variables["A_wave"] *
                   std::cos(variables["k_wave"] * variables["x_wave"]) *
                   std::cos(variables["omega_wave"] * t);
    double exp_r = variables["A_wave"] *
                   std::cos(variables["k_wave"] * variables["x_wave"] -
                             variables["omega_wave"] * t);
    return std::complex<double>(cos_t + exp_r, 0.0);
}

double MUGEUQFFModule38::computeQuantumTerm(double t_Hubble_val) {
    double unc   = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    auto   psi   = computePsiTotal(variables["t"]);
    double integ = variables["integral_psi"] * std::norm(psi);
    return (variables["hbar"] / unc) * integ * (2.0 * variables["pi"] / t_Hubble_val);
}

double MUGEUQFFModule38::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

double MUGEUQFFModule38::computeDMTerm() {
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3.0 * variables["G"] * variables["M"] /
                   std::pow(variables["r"], 3.0);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

double MUGEUQFFModule38::computeUgSum() {
    variables["Ug1"] = (variables["G"] * variables["M"]) /
                        (variables["r"] * variables["r"]);
    variables["Ug4"] = variables["Ug1"];
    return variables["Ug1"] + variables["Ug2"] +
           variables["Ug3_prime"] + variables["Ug4"];
}

double MUGEUQFFModule38::computeUQFF(double t) {
    variables["t"] = t;
    double Htz     = computeHtz(variables["z"]);
    double sc_corr = 1.0 - (variables["B"] / variables["B_crit"]);
    double fenv    = computeFenv(t);
    double m_fac   = 1.0 + fenv;

    double g_base  = (variables["G"] * variables["M"] * m_fac /
                       (variables["r"] * variables["r"])) *
                      (1.0 + Htz * t) * sc_corr;

    double ug_sum  = computeUgSum();
    double lam_t   = variables["Lambda"] * std::pow(variables["c"], 2.0) / 3.0;
    double quant   = computeQuantumTerm(variables["t_Hubble"]);
    double fluid   = computeFluidTerm(g_base);
    double dm      = computeDMTerm();

    return g_base + ug_sum + lam_t + quant + fluid + dm;
}

double MUGEUQFFModule38::computeHres(double t) {
    double A_res  = variables["k_A"] * variables["Z_nuc"] *
                    (variables["A_nuc"] / variables["A_H"]) *
                    (1.0 + variables["delta_pair"]);
    double sin_t  = A_res * std::sin(2.0 * variables["pi"] * variables["f_res"] * t);
    double U_dp   = variables["k_res"] *
                    (variables["A1"] * variables["A2"] /
                     std::pow(variables["f_dp"], 2.0)) *
                    std::cos(variables["phi_dp"]);
    double k_nuc  = variables["k_0"] *
                    (variables["N_nuc"] / variables["Z_nuc"]) *
                    (1.0 + variables["delta_pair"]);
    double res_b  = sin_t + U_dp * variables["SC_m"] * k_nuc + variables["S_shell"];

    double fenv   = computeFenv(t);
    return res_b + fenv * variables["SC_m"];
}

double MUGEUQFFModule38::computeDuniverse() {
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

std::string MUGEUQFFModule38::getEquationText() {
    return
        "MUGE UQFF 38-System (Doc 42):\n"
        "Same core as 29-system plus:\n"
        "  F_env += F_torque (spiral rotation) + F_shock (wind collision)\n"
        "  F_cosmo = QG_term + DM_term + GW_term  [auto-cascade on update]\n"
        "  SFR-dependent F_env: fenv += (SFR * t_yr) / M  [LAGOON, ORION, EAGLE]\n"
        "New systems: LAGOON_NEBULA, SPIRALS_SN, NGC6302, ORION_NEBULA,\n"
        "             YOUNG_STARS_OUTFLOW, EAGLE_NEBULA, GRAVITY_BIGBANG";
}

std::string MUGEUQFFModule38::getSolutions(double t) {
    double g    = computeUQFF(t);
    double hres = computeHres(t);
    double D    = computeDuniverse();
    double fenv = computeFenv(t);

    std::ostringstream ss;
    ss << std::scientific << std::setprecision(4);
    ss << "MUGEUQFFModule38 Solutions t=" << t << " s (system=" << static_cast<int>(current_system) << "):\n";
    ss << "  g_UQFF     = " << g    << " m/s^2\n";
    ss << "  H_res(t)   = " << hres << "\n";
    ss << "  D_universe = " << D    << " m\n";
    ss << "  F_env(t)   = " << fenv << "\n";
    ss << "  F_cosmo    = " << variables["F_cosmo"] << " [QG+DM+GW]\n";
    return ss.str();
}

void MUGEUQFFModule38::printVariables() {
    std::cout << "MUGEUQFFModule38 variables (system=" << static_cast<int>(current_system) << "):\n";
    for (const auto& p : variables) {
        std::cout << "  " << p.first << " = " << p.second << "\n";
    }
}
