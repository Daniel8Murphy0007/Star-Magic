// MUGEUQFFModule29.cpp
// MUGE/UQFF compressed gravity for 29 source documents, 8 system types.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#include "MUGEUQFFModule29.h"

MUGEUQFFModule29::MUGEUQFFModule29(SystemType29 sys) : current_system(sys) {
    // Universal constants
    variables["G"]           = 6.6743e-11;
    variables["c"]           = 3e8;
    variables["hbar"]        = 1.0546e-34;
    variables["Lambda"]      = 1.1e-52;
    variables["q"]           = 1.602e-19;
    variables["pi"]          = 3.141592653589793;
    variables["t_Hubble"]    = 13.8e9 * 3.156e7;
    variables["year_to_s"]   = 3.156e7;
    variables["H0"]          = 70.0;            // km/s/Mpc
    variables["Mpc_to_m"]    = 3.086e22;
    variables["Omega_m"]     = 0.3;
    variables["Omega_Lambda"]= 0.7;
    variables["z"]           = 0.0;
    variables["B_crit"]      = 1e11;            // T
    variables["D_p"]         = 4.4e26;          // m (Hubble radius)

    // Generic defaults
    double Msun = 1.989e30;
    variables["M"]           = Msun;
    variables["r"]           = 1e10;
    variables["M_visible"]   = 0.15 * variables["M"];
    variables["M_DM"]        = 0.85 * variables["M"];
    variables["rho_fluid"]   = 1e-20;
    variables["V"]           = 1e3;
    variables["B"]           = 1e-5;
    variables["Delta_x"]     = 1e-10;
    variables["Delta_p"]     = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"]= 1.0;
    variables["A_wave"]      = 1e-10;
    variables["k_wave"]      = 1e20;
    variables["omega_wave"]  = 1e15;
    variables["x_wave"]      = 0.0;
    variables["delta_rho"]   = 1e-21;
    variables["rho"]         = variables["rho_fluid"];

    // Ug components
    variables["Ug1"]         = 0.0;
    variables["Ug2"]         = 0.0;
    variables["Ug3_prime"]   = 0.0;
    variables["Ug4"]         = 0.0;

    // F_env components (13)
    variables["F_wind"]  = 0.0; variables["F_erode"]  = 0.0;
    variables["F_merge"] = 0.0; variables["F_SN"]     = 0.0;
    variables["F_rad"]   = 0.0; variables["F_fil"]    = 0.0;
    variables["F_BH"]    = 0.0; variables["F_dust"]   = 0.0;
    variables["F_ring"]  = 0.0; variables["F_mag"]    = 0.0;
    variables["F_tech"]  = 0.0; variables["F_shell"]  = 0.0;
    variables["F_cosmo"] = 0.0;
    variables["SC_m"]    = 1.0;

    // Resonance vars
    variables["A_res"]    = 1.0;
    variables["f_res"]    = 1e15;
    variables["k_A"]      = 1.0;
    variables["Z_nuc"]    = 1.0;
    variables["A_nuc"]    = 1.0;
    variables["A_H"]      = 1.0;
    variables["delta_pair"]= 0.0;
    variables["E_bind"]   = 13.6;             // eV
    variables["h_planck"] = 4.1356677e-15;    // eV·s
    variables["U_dp"]     = 0.0;
    variables["k_res"]    = 1.0;
    variables["A1"]       = 1.0;
    variables["A2"]       = 1.0;
    variables["f_dp"]     = 1e3;
    variables["phi_dp"]   = 0.0;
    variables["k_0"]      = 1.0;
    variables["N_nuc"]    = 1.0;
    variables["Z_magic"]  = 0.0;
    variables["N_magic"]  = 0.0;
    variables["S_shell"]  = 0.0;

    // Universe
    variables["M_total"]  = 1e53;
    variables["r_c"]      = 1e27;
    variables["k_cosm"]   = 0.0;

    variables["t"]        = 1e6 * variables["year_to_s"];

    setSystem(sys);
}

void MUGEUQFFModule29::setSystem(SystemType29 sys) {
    current_system = sys;
    double Msun = 1.989e30;
    double mp   = 1.6726e-27;
    double me   = 9.109e-31;

    switch (sys) {
        case SystemType29::SOMBRERO_GALAXY:
            variables["M"]     = 1e11 * Msun;
            variables["r"]     = 1.5e20;
            variables["z"]     = 0.0023;
            variables["M_ext"] = 1e9 * Msun;
            variables["r_ext"] = 1e19;
            variables["F_dust"]= 1e-11;
            variables["F_BH"]  = -1e-10;
            break;
        case SystemType29::SATURN:
            variables["M"]     = 5.68e26;
            variables["r"]     = 6e7;
            variables["z"]     = 0.0;
            variables["M_Sun"] = 1.989e30;
            variables["r_orbit"]= 1.43e12;
            variables["F_ring"]= 1e-10;
            variables["F_wind"]= 1e-11;
            break;
        case SystemType29::M16_EAGLE:
            variables["M"]     = 8e3 * Msun;
            variables["r"]     = 3.3e19;
            variables["z"]     = 0.002;
            variables["SFR"]   = 1.0 * Msun / variables["year_to_s"];
            variables["F_erode"]= -1e-12;
            break;
        case SystemType29::CRAB_NEBULA:
            variables["M"]     = 1.4 * Msun;
            variables["r"]     = 2.1e19;
            variables["z"]     = 0.004;
            variables["F_wind"]= 1e-10;
            variables["F_mag"] = 1e-9;
            variables["B"]     = 100.0;           // ~100 muG
            break;
        case SystemType29::HYDROGEN_ATOM:
            variables["M"]     = mp + me;
            variables["r"]     = 5.29e-11;        // Bohr radius
            variables["z"]     = 0.0;
            variables["E_n"]   = -13.6;            // eV
            variables["F_tech"]= 1e-20;
            break;
        case SystemType29::HYDROGEN_RESONANCE:
            variables["z"]     = 0.0;
            variables["f_res"] = (variables["E_bind"] / variables["h_planck"]) *
                                 (variables["A_H"] / variables["A_nuc"]) *
                                 (1 + variables["S_shell"]);
            variables["SC_m"]  = 1.0;
            break;
        case SystemType29::UNIVERSE_DIAMETER:
            variables["M"]     = 1e53;
            variables["r"]     = 1e27;
            variables["z"]     = 1100.0;
            variables["M_total"]= 1e53;
            variables["D_p"]   = 4.4e26;
            variables["k_cosm"]= 0.0;
            break;
        case SystemType29::GENERIC:
        default:
            break;
    }

    // Recompute dependencies
    variables["Delta_p"]  = variables["hbar"] / variables["Delta_x"];
    variables["M_visible"]= 0.15 * variables["M"];
    variables["M_DM"]     = 0.85 * variables["M"];
    if (variables.count("M_ext") && variables.count("r_ext")) {
        variables["Ug3_prime"] = (variables["G"] * variables["M_ext"]) /
                                  (variables["r_ext"] * variables["r_ext"]);
    }
    variables["S_shell"]  = 0.1 * (variables["Z_magic"] + variables["N_magic"]);
    variables["f_res"]    = (variables["E_bind"] / variables["h_planck"]) *
                             (variables["A_H"] / variables["A_nuc"]) *
                             (1.0 + variables["S_shell"]);
}

void MUGEUQFFModule29::updateVariable(const std::string& name, double value) {
    variables[name] = value;
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "M") {
        variables["M_visible"] = 0.15 * value;
        variables["M_DM"]      = 0.85 * value;
    }
}

void MUGEUQFFModule29::addToVariable(const std::string& name, double delta) {
    if (variables.count(name)) variables[name] += delta;
}

void MUGEUQFFModule29::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

double MUGEUQFFModule29::computeHtz(double z) {
    double Hz_kms = variables["H0"] * std::sqrt(
        variables["Omega_m"] * std::pow(1.0 + z, 3.0) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

double MUGEUQFFModule29::computeFenv(double t) {
    double fenv = 0.0;
    const std::string fkeys[] = {
        "F_wind","F_erode","F_merge","F_SN","F_rad","F_fil",
        "F_BH","F_dust","F_ring","F_mag","F_tech","F_shell","F_cosmo"
    };
    for (const auto& k : fkeys) fenv += variables.count(k) ? variables.at(k) : 0.0;
    if (variables.count("SFR")) {
        double t_yr = t / variables["year_to_s"];
        fenv += (variables["SFR"] * t_yr) / variables["M"];
    }
    return fenv;
}

std::complex<double> MUGEUQFFModule29::computePsiTotal(double t) {
    double cos_term = 2.0 * variables["A_wave"] *
                      std::cos(variables["k_wave"] * variables["x_wave"]) *
                      std::cos(variables["omega_wave"] * t);
    std::complex<double> exp_term(
        variables["A_wave"] * std::cos(variables["k_wave"] * variables["x_wave"] -
                                        variables["omega_wave"] * t), 0.0);
    return std::complex<double>(cos_term + exp_term.real(), 0.0);
}

double MUGEUQFFModule29::computeQuantumTerm(double t_Hubble_val) {
    double unc   = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    auto   psi   = computePsiTotal(variables["t"]);
    double integ = variables["integral_psi"] * std::norm(psi);
    return (variables["hbar"] / unc) * integ * (2.0 * variables["pi"] / t_Hubble_val);
}

double MUGEUQFFModule29::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

double MUGEUQFFModule29::computeDMTerm() {
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3.0 * variables["G"] * variables["M"] /
                   std::pow(variables["r"], 3.0);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

double MUGEUQFFModule29::computeUgSum() {
    variables["Ug1"] = (variables["G"] * variables["M"]) /
                        (variables["r"] * variables["r"]);
    variables["Ug4"] = variables["Ug1"] * 1.0;
    return variables["Ug1"] + variables["Ug2"] +
           variables["Ug3_prime"] + variables["Ug4"];
}

double MUGEUQFFModule29::computeUQFF(double t) {
    variables["t"] = t;
    double Htz      = computeHtz(variables["z"]);
    double expansion= 1.0 + Htz * t;
    double sc_corr  = 1.0 - (variables["B"] / variables["B_crit"]);
    double fenv     = computeFenv(t);
    double m_factor = 1.0 + fenv;

    double g_base   = (variables["G"] * variables["M"] * m_factor /
                        (variables["r"] * variables["r"])) * expansion * sc_corr;

    double ug_sum   = computeUgSum();
    double lambda_t = variables["Lambda"] * std::pow(variables["c"], 2.0) / 3.0;
    double quant    = computeQuantumTerm(variables["t_Hubble"]);
    double fluid    = computeFluidTerm(g_base);
    double dm       = computeDMTerm();

    return g_base + ug_sum + lambda_t + quant + fluid + dm;
}

double MUGEUQFFModule29::computeHres(double t) {
    // Resonance base: A_res sin(2πf_res t) + U_dp*SC_m*k_nuc + S_shell
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
    double S_sh   = variables["S_shell"];
    double res_base = sin_t + U_dp * variables["SC_m"] * k_nuc + S_sh;

    double fenv   = computeFenv(t);
    return res_base + fenv * variables["SC_m"];
}

double MUGEUQFFModule29::computeDuniverse() {
    double Hz     = computeHtz(variables["z"]);
    double t      = variables["t"];
    double G      = variables["G"];
    double hbar   = variables["hbar"];
    double Dx     = variables["Delta_x"];
    double Dp     = variables["Delta_p"];
    double M_tot  = variables["M_total"];
    double D_p    = variables["D_p"];
    double k      = variables["k_cosm"];
    double r_c    = variables["r_c"];
    double Lambda = variables["Lambda"];
    double c      = variables["c"];
    double H0_si  = (variables["H0"] * 1e3) / variables["Mpc_to_m"];

    double exp_fac     = 1.0 + Hz * t;
    double lambda_fac  = 1.0 + Lambda * c * c / (3.0 * H0_si * H0_si);
    double quant_fac   = 1.0 + hbar / (std::sqrt(Dx * Dp) * G * M_tot);
    double curv_fac    = 1.0 + k * r_c * r_c;

    return 2.0 * D_p * exp_fac * lambda_fac * quant_fac * curv_fac;
}

std::string MUGEUQFFModule29::getEquationText() {
    return
        "MUGE UQFF 29-System (Doc 41):\n"
        "g_UQFF(r,t) = [G M(t)/r^2] * (1+H(t,z)) * (1-B/B_crit) * (1+F_env(t))\n"
        "            + (Ug1+Ug2+Ug3'+Ug4)\n"
        "            + Lambda*c^2/3\n"
        "            + (hbar/sqrt(Delta_x*Delta_p)) * |psi_total|^2 * (2pi/t_Hubble)\n"
        "            + rho_fluid*V*g_base\n"
        "            + (M_vis+M_DM)*(delta_rho/rho + 3GM/r^3)\n\n"
        "H_res(t) = A_res*sin(2pi*f_res*t) + U_dp*SC_m*k_nuc + S_shell + F_env(t)*SC_m\n\n"
        "D_universe = 2*D_p*(1+Hz*t)*(1+Lambda*c^2/(3*H0^2))*(1+hbar/(sqrt(Dx*Dp)*G*M))*(1+k*r_c^2)\n\n"
        "F_env = F_wind+F_erode+F_merge+F_SN+F_rad+F_fil+F_BH+F_dust+F_ring+F_mag+F_tech+F_shell+F_cosmo";
}

std::string MUGEUQFFModule29::getSolutions(double t) {
    double g    = computeUQFF(t);
    double hres = computeHres(t);
    double D    = computeDuniverse();
    double fenv = computeFenv(t);
    double ug   = computeUgSum();
    double lam  = variables["Lambda"] * std::pow(variables["c"], 2.0) / 3.0;

    std::ostringstream ss;
    ss << std::scientific << std::setprecision(4);
    ss << "MUGEUQFFModule29 Solutions t=" << t << " s (system=" << static_cast<int>(current_system) << "):\n";
    ss << "  g_UQFF    = " << g    << " m/s^2\n";
    ss << "  H_res(t)  = " << hres << " (mixed units)\n";
    ss << "  D_universe= " << D    << " m\n";
    ss << "  F_env(t)  = " << fenv << "\n";
    ss << "  Ug_sum    = " << ug   << " m/s^2\n";
    ss << "  Lambda_t  = " << lam  << "\n";
    return ss.str();
}

void MUGEUQFFModule29::printVariables() {
    std::cout << "MUGEUQFFModule29 variables (system=" << static_cast<int>(current_system) << "):\n";
    for (const auto& p : variables) {
        std::cout << "  " << p.first << " = " << p.second << "\n";
    }
}
