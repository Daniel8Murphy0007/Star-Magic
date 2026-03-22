// MultiUQFFCompressionModule.cpp
// Docs 39, 40, 41: 29-system UQFF Compression Cycle 2.
//
// Core equation (all systems):
//   g(r,t) = [G M(t)/r²] * (1+H(t,z)) * (1-B/B_crit) * (1+F_env(t))
//          + (Ug1 + Ug3' + Ug4) + (Λc²/3)
//          + quantum_psi + ρ_fluid V g_base + DM_pert
//          [+ H_res(t) for atomic/quantum systems]
//
// F_env(t) = 1 + ΣF_i(t)  — system-specific:
//   * winds:   ρ v_wind²
//   * erosion: 1-exp(-t/τ)
//   * decay:   exp(-t/τ)
//   * lensing: 1+0.1 sin(2π t/t_H)
//   * merger:  0.1 M (1-exp(-t/τ_merge)) / M
//   * BH:      G M_ext/r_ext²
//   * SN:     -M_SN (1-exp(-t/τ_SN)) / M
//   * fil:     1e-10 B r
//
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#include "MultiUQFFCompressionModule.h"
#include "UQFFConstants.h"
#include <cmath>
#include <sstream>

MultiUQFFCompressionModule::MultiUQFFCompressionModule(const std::string& system) {
    // Universal constants
    variables["G"]           = UQFF::G;
    variables["c"]           = UQFF::c;
    variables["hbar"]        = UQFF::hbar;
    variables["Lambda"]      = UQFF::Lambda;
    variables["pi"]          = UQFF::pi;
    variables["t_Hubble"]    = UQFF::t_Hubble;
    variables["year_to_s"]   = UQFF::year_to_s;
    variables["H0"]          = UQFF::H0;
    variables["Mpc_to_m"]    = UQFF::Mpc_to_m;
    variables["Omega_m"]     = UQFF::Omega_m;
    variables["Omega_Lambda"] = UQFF::Omega_Lambda;
    variables["B_crit"]      = UQFF::B_crit_magnetar;
    variables["m_H"]         = UQFF::m_H;
    // Quantum defaults
    variables["Delta_x"]     = 1.0e-10;
    variables["Delta_p"]     = UQFF::hbar / 1.0e-10;
    variables["integral_psi"] = 1.0;
    variables["f_sc"]        = 10.0;
    variables["delta_rho_over_rho"] = 1.0e-5;
    // H_res params
    variables["A_res"]       = 1.0e-10;
    variables["f_res"]       = 1.0e15;
    variables["SC_m"]        = 1.0;
    variables["is_resonance"] = 0.0;

    loadSystem(system);
}

void MultiUQFFCompressionModule::loadSystem(const std::string& sys) {
    current_system = sys;
    double Ms = UQFF::M_sun;

    // Reset per-system fields
    auto resetField = [&](const std::string& k, double v = 0.0) { variables[k] = v; };
    resetField("M"); resetField("r", 1.0); resetField("z"); resetField("t_default");
    resetField("SFR"); resetField("M0"); resetField("M_visible"); resetField("M_DM");
    resetField("M_ext"); resetField("r_ext", 1.0); resetField("v_wind");
    resetField("M_SN"); resetField("tau_merge", 1.0e17); resetField("tau_SN", 1.0e16);
    resetField("tau_erode", 1.0e14); resetField("tau_decay", 1.0e10);
    resetField("B", 1.0e-5); resetField("rho_fluid", 1.0e-20);
    resetField("is_resonance");

    // ---- System parameters ----
    if (sys == "MagnetarSGR1745") {
        variables["M"]    = 2.8 * Ms;    variables["r"] = 1.0e4;
        variables["z"]    = 0.026;        variables["t_default"] = 1.0e3 * UQFF::year_to_s;
        variables["M_ext"]= 4.0e6 * Ms;  variables["r_ext"] = 8.0e9;
        variables["v_wind"]= 1.0e5;       variables["B"] = 1.0e11;
        variables["tau_decay"]= 1.0e3 * UQFF::year_to_s;
        variables["rho_fluid"]= 1.0e17;
    } else if (sys == "SagittariusA") {
        variables["M"]    = 4.0e6 * Ms;  variables["r"] = 1.0e10;
        variables["z"]    = 0.0;          variables["t_default"] = 1.0e6 * UQFF::year_to_s;
        variables["v_wind"]= 1.0e8;
    } else if (sys == "TapestryStarbirth" || sys == "Westerlund2") {
        variables["M"]    = 1.0e4 * Ms;  variables["r"] = 1.0e18;
        variables["z"]    = 0.001;        variables["t_default"] = 5.0e6 * UQFF::year_to_s;
        variables["SFR"]  = 0.1 * Ms;    variables["v_wind"] = 1.0e3;
    } else if (sys == "PillarsCreation") {
        variables["M"]    = 800.0 * Ms;  variables["r"] = 3.0e17;
        variables["z"]    = 0.0018;       variables["t_default"] = 2.0e6 * UQFF::year_to_s;
        variables["SFR"]  = 0.1 * Ms;    variables["v_wind"] = 1.0e4;
        variables["tau_erode"] = 2.0e6 * UQFF::year_to_s;
    } else if (sys == "RingsRelativity") {
        variables["M"]    = 1.0e11 * Ms; variables["r"] = 1.0e21;
        variables["z"]    = 0.5;          variables["t_default"] = 10.0e9 * UQFF::year_to_s;
    } else if (sys == "NGC2525") {
        variables["M"]    = 1.0e10 * Ms; variables["r"] = 1.0e20;
        variables["z"]    = 0.01;         variables["t_default"] = 1.0e9 * UQFF::year_to_s;
        variables["SFR"]  = 1.0 * Ms;    variables["v_wind"] = 1.0e3;
        variables["M_ext"]= 1.0e9 * Ms;  variables["r_ext"] = 1.0e19;
        variables["M_SN"] = 10.0 * Ms;   variables["tau_SN"] = 1.0e8 * UQFF::year_to_s;
    } else if (sys == "NGC3603") {
        variables["M"]    = 2.0e4 * Ms;  variables["r"] = 2.0e18;
        variables["z"]    = 0.001;        variables["t_default"] = 3.0e6 * UQFF::year_to_s;
        variables["SFR"]  = 0.2 * Ms;    variables["v_wind"] = 2.0e3;
        variables["tau_erode"] = 3.0e6 * UQFF::year_to_s;
    } else if (sys == "BubbleNebula") {
        variables["M"]    = 5.0e3 * Ms;  variables["r"] = 5.0e17;
        variables["z"]    = 0.001;        variables["t_default"] = 4.0e6 * UQFF::year_to_s;
        variables["SFR"]  = 0.05 * Ms;   variables["v_wind"] = 5.0e3;
        variables["tau_erode"] = 4.0e6 * UQFF::year_to_s;
    } else if (sys == "AntennaeGalaxies") {
        variables["M"]    = 1.0e11 * Ms; variables["r"] = 5.0e20;
        variables["z"]    = 0.025;        variables["t_default"] = 5.0e8 * UQFF::year_to_s;
        variables["SFR"]  = 10.0 * Ms;   variables["v_wind"] = 1.0e4;
        variables["M_ext"]= 5.0e10 * Ms; variables["r_ext"] = 1.0e20;
        variables["tau_merge"] = 5.0e8 * UQFF::year_to_s;
    } else if (sys == "HorseheadNebula") {
        variables["M"]    = 1.0e3 * Ms;  variables["r"] = 1.0e17;
        variables["z"]    = 0.0;          variables["t_default"] = 1.0e6 * UQFF::year_to_s;
        variables["SFR"]  = 0.01 * Ms;   variables["v_wind"] = 1.0e3;
        variables["tau_erode"] = 1.0e6 * UQFF::year_to_s;
    } else if (sys == "NGC1275") {
        variables["M"]    = 1.0e11 * Ms; variables["r"] = 1.0e21;
        variables["z"]    = 0.017;        variables["t_default"] = 1.0e9 * UQFF::year_to_s;
        variables["SFR"]  = 0.5 * Ms;    variables["v_wind"] = 1.0e4;
        variables["M_ext"]= 8.0e9 * Ms;  variables["r_ext"] = 1.0e19;
        variables["B"]    = 1.0e-4;       // filament field
    } else if (sys == "NGC1792") {
        variables["M"]    = 5.0e10 * Ms; variables["r"] = 5.0e20;
        variables["z"]    = 0.012;        variables["t_default"] = 8.0e8 * UQFF::year_to_s;
        variables["SFR"]  = 2.0 * Ms;    variables["v_wind"] = 2.0e3;
        variables["M_SN"] = 20.0 * Ms;   variables["tau_SN"] = 8.0e8 * UQFF::year_to_s;
    } else if (sys == "HubbleUltraDeepField") {
        variables["M"]    = 1.0e12 * Ms; variables["r"] = 1.0e23;
        variables["z"]    = 10.0;         variables["t_default"] = 10.0e9 * UQFF::year_to_s;
    } else if (sys == "StudentsGuideUniverse") {
        variables["M"]    = 1.0 * Ms;    variables["r"] = UQFF::AU;
        variables["z"]    = 0.0;          variables["t_default"] = 4.35e17;
    } else if (sys == "SombreroGalaxy") {
        variables["M"]    = 8.0e11 * Ms; variables["r"] = 4.73e20;
        variables["z"]    = 0.002;        variables["t_default"] = 4.35e17;
        variables["M_DM"] = 0.8 * 8.0e11 * Ms;
    } else if (sys == "Saturn") {
        variables["M"]    = 5.68e26;     variables["r"] = 6.027e7;
        variables["z"]    = 0.0;          variables["t_default"] = 4.35e17;
    } else if (sys == "EagleNebula") {
        variables["M"]    = 5.0e3 * Ms;  variables["r"] = 3.31e17;
        variables["z"]    = 0.0018;       variables["t_default"] = 10.0e6 * UQFF::year_to_s;
        variables["SFR"]  = 0.2 * Ms;    variables["v_wind"] = 2.0e6;
    } else if (sys == "CrabNebula") {
        variables["M"]    = 5.0 * Ms;    variables["r"] = 5.203e16;
        variables["z"]    = 0.00002;      variables["t_default"] = 971.0 * UQFF::year_to_s;
        variables["v_wind"]= 1.34e6;     // expansion velocity
    } else if (sys == "HydrogenAtom") {
        variables["M"]    = UQFF::m_H;   variables["r"] = UQFF::a0;
        variables["z"]    = 0.0;          variables["t_default"] = 4.35e17;
        variables["is_resonance"] = 1.0;
    } else if (sys == "HydrogenResonance") {
        variables["M"]    = UQFF::m_H;   variables["r"] = UQFF::a0;
        variables["z"]    = 0.0;          variables["t_default"] = 4.35e17;
        variables["f_res"] = UQFF::c / (4.0 * UQFF::a0);
        variables["is_resonance"] = 1.0;
    } else if (sys == "OrionNebula") {
        variables["M"]    = 2.0e3 * Ms;  variables["r"] = 1.18e17;
        variables["z"]    = 0.0004;       variables["t_default"] = 1.0e6 * UQFF::year_to_s;
        variables["SFR"]  = 0.1 * Ms;    variables["v_wind"] = 2.5e5;
    } else if (sys == "GalaxiesGalore") {
        variables["M"]    = 1.0e11 * Ms; variables["r"] = 1.543e21;
        variables["z"]    = 1.0;          variables["t_default"] = 4.35e17;
    } else if (sys == "NewStars" || sys == "StellarForge") {
        variables["M"]    = 5.0e3 * Ms;  variables["r"] = 5.0e17;
        variables["z"]    = 0.001;        variables["t_default"] = 5.0e6 * UQFF::year_to_s;
        variables["SFR"]  = 0.05 * Ms;   variables["v_wind"] = 5.0e2;
    } else if (sys == "LagoonNebula") {
        variables["M"]    = 1.0e4 * Ms;  variables["r"] = 5.203e17;
        variables["z"]    = 0.0001;       variables["t_default"] = 2.0e6 * UQFF::year_to_s;
        variables["SFR"]  = 0.5 * Ms;
    } else if (sys == "SpiralsSupernovae") {
        variables["M"]    = 1.0e11 * Ms; variables["r"] = 1.543e21;
        variables["z"]    = 0.002;        variables["t_default"] = 4.35e17;
    } else if (sys == "NGC6302") {
        variables["M"]    = 1.0 * Ms;    variables["r"] = 1.514e16;
        variables["z"]    = 0.00001;      variables["t_default"] = 1.0e4 * UQFF::year_to_s;
    } else if (sys == "UniverseDiameter") {
        variables["M"]    = 1.5e53;      variables["r"] = 4.4e26;
        variables["z"]    = 1100.0;       variables["t_default"] = 4.35e17;
        variables["M_DM"] = 0.85 * 1.5e53;
    }

    if (variables["M0"] == 0.0) variables["M0"] = variables["M"];
    if (variables["M_visible"] == 0.0) variables["M_visible"] = variables["M"];
    variables["V"] = 1.0 / std::max(variables["rho_fluid"], 1.0e-50);
}

void MultiUQFFCompressionModule::setSystem(const std::string& system) {
    loadSystem(system);
}

void MultiUQFFCompressionModule::onVariableUpdated(const std::string& name, double value) {
    if (name == "Delta_x" && value > 0.0)
        variables["Delta_p"] = UQFF::hbar / value;
    else if (name == "M") {
        variables["M0"]       = value;
        variables["M_visible"]= value;
    } else if (name == "rho_fluid" && value > 0.0)
        variables["V"] = 1.0 / value;
}

double MultiUQFFCompressionModule::computeHtz() const {
    double z = variables.at("z");
    return variables.at("H0") * std::sqrt(
        variables.at("Omega_m") * std::pow(1.0 + z, 3) +
        variables.at("Omega_Lambda")) * 1.0e3 / variables.at("Mpc_to_m");
}

// F_env(t) = 1 + system-specific environmental sum
double MultiUQFFCompressionModule::computeF_env(double t) const {
    double f_env = 1.0;
    double t_yr  = t / variables.at("year_to_s");
    double rho   = variables.at("rho_fluid");
    double vw    = variables.at("v_wind");
    std::string s = current_system;

    // Wind term (most systems)
    if (vw > 0.0) f_env += rho * vw * vw;

    // System-specific additions
    if (s == "MagnetarSGR1745") {
        double M_mag = 1.0e40;
        double Mc2   = variables.at("M") * variables.at("c") * variables.at("c");
        double D_t   = std::exp(-t / variables.at("tau_decay"));
        double BH    = variables.at("G") * variables.at("M_ext") /
                       (variables.at("r_ext") * variables.at("r_ext"));
        f_env += M_mag / Mc2 + D_t + BH;
    } else if (s == "PillarsCreation" || s == "BubbleNebula" ||
               s == "HorseheadNebula" || s == "NGC3603") {
        double E_t = 1.0 - std::exp(-t / variables.at("tau_erode"));
        f_env += E_t;
    } else if (s == "RingsRelativity") {
        double L_t = 1.0 + 0.1 * std::sin(2.0 * variables.at("pi") * t / variables.at("t_Hubble"));
        f_env += L_t;
    } else if (s == "AntennaeGalaxies") {
        double M_merge = 0.1 * variables.at("M") *
                         (1.0 - std::exp(-t / variables.at("tau_merge")));
        f_env += M_merge / std::max(variables.at("M"), 1.0e-30);
        // BH influence
        if (variables.at("M_ext") > 0.0)
            f_env += variables.at("G") * variables.at("M_ext") /
                     (variables.at("r_ext") * variables.at("r_ext"));
    } else if (s == "NGC2525" || s == "NGC1792") {
        double M_loss = variables.at("M_SN") *
                        (1.0 - std::exp(-t_yr / (variables.at("tau_SN") / variables.at("year_to_s"))));
        f_env -= M_loss / std::max(variables.at("M"), 1.0e-30);
        if (s == "NGC2525" && variables.at("M_ext") > 0.0)
            f_env += variables.at("G") * variables.at("M_ext") /
                     (variables.at("r_ext") * variables.at("r_ext"));
    } else if (s == "NGC1275") {
        double F_fil = 1.0e-10 * variables.at("B") * variables.at("r");
        f_env += F_fil;
        if (variables.at("M_ext") > 0.0)
            f_env += variables.at("G") * variables.at("M_ext") /
                     (variables.at("r_ext") * variables.at("r_ext"));
    } else if (s == "HubbleUltraDeepField") {
        f_env += 0.01 * (t / variables.at("t_Hubble"));
    } else if (s == "SagittariusA") {
        // GW correction: (GM)²/(c⁴ r) * ω_dot²  (schematic)
        double omega_dot = 1.0e-3;
        double GM = variables.at("G") * variables.at("M");
        double r  = variables.at("r");
        double c4 = std::pow(variables.at("c"), 4);
        f_env += GM * GM / (c4 * r) * omega_dot * omega_dot;
    }

    return f_env;
}

double MultiUQFFCompressionModule::computeMsfFactor(double t) const {
    double SFR = variables.at("SFR");
    double M0  = variables.at("M0");
    if (SFR == 0.0 || M0 == 0.0) return 0.0;
    return SFR * (t / variables.at("year_to_s")) / M0;
}

double MultiUQFFCompressionModule::computeUgSum() const {
    double G   = variables.at("G");
    double M   = variables.at("M");
    double r   = variables.at("r");
    double Ug1 = G * M / (r * r);
    double Ug3 = (variables.at("M_ext") > 0.0) ?
        G * variables.at("M_ext") / (variables.at("r_ext") * variables.at("r_ext")) : 0.0;
    double Ug4 = Ug1 * variables.at("f_sc");
    return Ug1 + Ug3 + Ug4;
}

double MultiUQFFCompressionModule::computeQuantumPsi() const {
    double unc = std::sqrt(variables.at("Delta_x") * variables.at("Delta_p"));
    return (variables.at("hbar") / unc) *
           variables.at("integral_psi") *
           (2.0 * variables.at("pi") / variables.at("t_Hubble"));
}

double MultiUQFFCompressionModule::computeFluidTerm(double g_base) const {
    return variables.at("rho_fluid") * variables.at("V") * g_base;
}

double MultiUQFFCompressionModule::computeDMPert() const {
    double r   = variables.at("r");
    double drr = variables.at("delta_rho_over_rho");
    double pert = drr + 3.0 * variables.at("G") * variables.at("M") / std::pow(r, 3);
    return (variables.at("M_visible") + variables.at("M_DM")) * pert;
}

// H_res = A_res * sin(2π f_res t) + F_env * SC_m  (for atomic/quantum systems)
double MultiUQFFCompressionModule::computeH_res(double t) const {
    double A = variables.at("A_res");
    double f = variables.at("f_res");
    double pi = variables.at("pi");
    double F_env = computeF_env(t);
    double SC_m  = variables.at("SC_m");
    return A * std::sin(2.0 * pi * f * t) + F_env * SC_m;
}

double MultiUQFFCompressionModule::computeG(double t) {
    variables["t"] = t;

    double Hz        = computeHtz();
    double expansion = 1.0 + Hz * t;
    double sc_corr   = 1.0 - variables["B"] / variables["B_crit"];
    double f_env     = computeF_env(t);
    double msf       = computeMsfFactor(t);
    double m_factor  = 1.0 + msf;

    double r       = variables["r"];
    double g_base  = (variables["G"] * variables["M"] * m_factor / (r * r))
                     * expansion * sc_corr * f_env;

    double ug_sum  = computeUgSum();
    double lam     = variables["Lambda"] * variables["c"] * variables["c"] / 3.0;
    double quantum = computeQuantumPsi();
    double fluid   = computeFluidTerm(g_base);
    double dm_pert = computeDMPert();

    double h_res = (variables["is_resonance"] != 0.0) ? computeH_res(t) : 0.0;

    return g_base + ug_sum + lam + quantum + fluid + dm_pert + h_res;
}

std::string MultiUQFFCompressionModule::getEquationText() const {
    std::ostringstream oss;
    oss << "MultiUQFFCompressionModule — system: " << current_system << "\n"
        << "g(r,t) = [G M(t)/r²] * (1+H(t,z)) * (1-B/B_crit) * (1+F_env(t))\n"
        << "       + (Ug1 + Ug3' + Ug4) + (Λc²/3)\n"
        << "       + quantum_psi + ρ_fluid V g_base + DM_pert\n"
        << "       [+ H_res(t) for atomic/quantum systems]\n"
        << "\n"
        << "F_env(t) = 1 + ΣF_i(t)  [modular sum per system]:\n"
        << "  F_wind    = ρ v_wind²\n"
        << "  F_erode   = 1-exp(-t/τ_erode)         [Pillars, Bubble, Horsehead, NGC3603]\n"
        << "  F_lensing = 1+0.1 sin(2πt/t_H)        [RingsRelativity]\n"
        << "  F_merge   = 0.1M(1-exp(-t/τ_merge))/M [AntennaeGalaxies]\n"
        << "  F_SN      = -M_SN(1-exp(-t/τ_SN))/M   [NGC2525, NGC1792]\n"
        << "  F_fil     = 1e-10 B r                  [NGC1275 filaments]\n"
        << "  F_BH      = GM_ext/r_ext²              [NGC1275, NGC2525, Magnetar]\n"
        << "  F_decay   = exp(-t/τ_decay)            [MagnetarSGR1745]\n"
        << "  F_evo     = 0.01*(t/t_H)               [HubbleUltraDeepField]\n"
        << "\n"
        << "H_res(t) = A_res sin(2πf_res t) + F_env(t)*SC_m  [HydrogenAtom/Resonance]\n"
        << "Supports 29 systems — see loadSystem() for full parameter table.";
    return oss.str();
}
