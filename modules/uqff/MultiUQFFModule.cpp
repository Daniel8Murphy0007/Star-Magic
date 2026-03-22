// MultiUQFFModule.cpp
// Docs 34a + 34b: 15-system UQFF dispatcher.
//
// General equation:
//   g(r,t) = [G M(t)/r²] * (1 + H(t,z)) * (1 - B/B_crit) * (1 + F_env(t))
//          + (Ug1 + Ug3' + Ug4)
//          + (Λc²/3)
//          + quantum ψ_total term
//          + ρ_fluid V g_base
//          + DM perturbation
//   [+ H_res(t) for atomic/quantum systems]
//
// Resonance solutions per system are in UQFFResonanceValues.h.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#include "MultiUQFFModule.h"
#include "UQFFConstants.h"
#include <cmath>
#include <complex>
#include <sstream>

MultiUQFFModule::MultiUQFFModule(const std::string& system) {
    // Universal constants - loaded once
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
    // Quantum uncertainty defaults
    variables["Delta_x"]     = 1.0e-10;
    variables["Delta_p"]     = UQFF::hbar / 1.0e-10;
    variables["integral_psi"] = 1.0;
    variables["f_sc"]        = 10.0;
    // Resonance defaults
    variables["A_res"]       = 1.0e-10;
    variables["f_res"]       = 1.0e15;   // Hz
    loadSystem(system);
}

void MultiUQFFModule::loadSystem(const std::string& sys) {
    current_system = sys;
    double M_sun   = UQFF::M_sun;

    // Reset per-system fields
    variables["M"]        = 0.0;
    variables["r"]        = 0.0;
    variables["z"]        = 0.0;
    variables["t_default"]= 0.0;
    variables["SFR"]      = 0.0;
    variables["M0"]       = 0.0;
    variables["M_visible"]= 0.0;
    variables["M_DM"]     = 0.0;
    variables["M_ext"]    = 0.0;
    variables["r_ext"]    = 1.0;
    variables["v_wind"]   = 0.0;
    variables["v_exp"]    = 0.0;
    variables["rho_fluid"]= 1.0e-20;
    variables["B"]        = 1.0e-5;
    variables["is_resonance_sys"] = 0.0;

    if (sys == "UniverseDiameter") {
        variables["M"]        = 1.5e53;
        variables["r"]        = 4.4e26;
        variables["z"]        = 1100.0;
        variables["t_default"]= 4.35e17;
        variables["M_visible"]= 0.15 * variables["M"];
        variables["M_DM"]     = 0.85 * variables["M"];
    } else if (sys == "HydrogenAtom") {
        variables["M"]        = UQFF::m_H;
        variables["r"]        = UQFF::a0;
        variables["z"]        = 0.0;
        variables["t_default"]= 4.35e17;
        variables["M_visible"]= variables["M"];
        variables["is_resonance_sys"] = 1.0;
    } else if (sys == "HydrogenResonancePToE") {
        variables["M"]        = UQFF::m_H;
        variables["r"]        = UQFF::a0;
        variables["z"]        = 0.0;
        variables["t_default"]= 4.35e17;
        variables["M_visible"]= variables["M"];
        variables["is_resonance_sys"] = 1.0;
        variables["f_res"]    = UQFF::c / (4.0 * UQFF::a0);  // proton-electron Lyman limit
    } else if (sys == "LagoonNebula") {
        variables["M"]        = 1.0e4 * M_sun;
        variables["r"]        = 5.203e17;
        variables["z"]        = 0.0001;
        variables["t_default"]= 2.0e6 * UQFF::year_to_s;
        variables["SFR"]      = 0.5 * M_sun;
        variables["M0"]       = variables["M"];
        variables["M_visible"]= variables["M"];
    } else if (sys == "SpiralsSupernovae") {
        variables["M"]        = 1.0e11 * M_sun;
        variables["r"]        = 1.543e21;
        variables["z"]        = 0.002;
        variables["t_default"]= 4.35e17;
        variables["M_visible"]= variables["M"];
        variables["M_DM"]     = variables["M"] * 5.0;  // DM dominated
    } else if (sys == "NGC6302") {
        variables["M"]        = 1.0 * M_sun;
        variables["r"]        = 1.514e16;
        variables["z"]        = 0.00001;
        variables["t_default"]= 10.0e3 * UQFF::year_to_s;
        variables["M_visible"]= variables["M"];
    } else if (sys == "OrionNebula") {
        variables["M"]        = 2000.0 * M_sun;
        variables["r"]        = 1.18e17;
        variables["z"]        = 0.0004;
        variables["t_default"]= 1.0e6 * UQFF::year_to_s;
        variables["SFR"]      = 0.1 * M_sun;
        variables["M0"]       = variables["M"];
        variables["M_visible"]= variables["M"];
        variables["v_wind"]   = 2.5e5;
    } else if (sys == "UniverseGuide") {
        variables["M"]        = 1.0 * M_sun;
        variables["r"]        = UQFF::AU;
        variables["z"]        = 0.0;
        variables["t_default"]= 4.35e17;
        variables["M_visible"]= variables["M"];
    } else if (sys == "GalaxiesGalore") {
        variables["M"]        = 1.0e11 * M_sun;
        variables["r"]        = 1.543e21;
        variables["z"]        = 1.0;
        variables["t_default"]= 4.35e17;
        variables["M_visible"]= 0.1 * variables["M"];
        variables["M_DM"]     = 0.9 * variables["M"];
    } else if (sys == "StellarForge") {
        variables["M"]        = 5.0e3 * M_sun;
        variables["r"]        = 5.0e17;
        variables["z"]        = 0.001;
        variables["t_default"]= 5.0e6 * UQFF::year_to_s;
        variables["SFR"]      = 0.05 * M_sun;
        variables["M0"]       = variables["M"];
        variables["M_visible"]= variables["M"];
    } else if (sys == "SombreroGalaxy") {
        variables["M"]        = 8.0e11 * M_sun;
        variables["r"]        = 4.73e20;
        variables["z"]        = 0.002;
        variables["t_default"]= 4.35e17;
        variables["M_visible"]= 0.2 * variables["M"];
        variables["M_DM"]     = 0.8 * variables["M"];
    } else if (sys == "Saturn") {
        variables["M"]        = 5.68e26;
        variables["r"]        = 6.027e7;
        variables["z"]        = 0.0;
        variables["t_default"]= 4.35e17;
        variables["M_visible"]= variables["M"];
    } else if (sys == "CrabNebula") {
        variables["M"]        = 5.0 * M_sun;
        variables["r"]        = 5.203e16;
        variables["z"]        = 0.00002;
        variables["t_default"]= 971.0 * UQFF::year_to_s;
        variables["v_exp"]    = 1.34e6;   // m/s
        variables["M_visible"]= variables["M"];
    } else if (sys == "NewStars") {
        variables["M"]        = 5.0e3 * M_sun;
        variables["r"]        = 5.0e17;
        variables["z"]        = 0.001;
        variables["t_default"]= 5.0e6 * UQFF::year_to_s;
        variables["SFR"]      = 0.1 * M_sun;
        variables["M0"]       = variables["M"];
        variables["M_visible"]= variables["M"];
    }

    if (variables["M0"] == 0.0) variables["M0"] = variables["M"];
    if (variables["M_visible"] == 0.0) variables["M_visible"] = variables["M"];
    variables["V"] = 1.0 / ((variables["rho_fluid"] > 0.0) ? variables["rho_fluid"] : 1.0e-20);
}

void MultiUQFFModule::setSystem(const std::string& system) {
    loadSystem(system);
}

void MultiUQFFModule::onVariableUpdated(const std::string& name, double value) {
    if (name == "Delta_x" && value > 0.0)
        variables["Delta_p"] = UQFF::hbar / value;
    else if (name == "M") {
        variables["M0"]       = value;
        variables["M_visible"] = value;
    } else if (name == "rho_fluid" && value > 0.0) {
        variables["V"] = 1.0 / value;
    }
}

double MultiUQFFModule::computeHtz() const {
    double z  = variables.at("z");
    return variables.at("H0") * std::sqrt(
        variables.at("Omega_m") * std::pow(1.0 + z, 3) +
        variables.at("Omega_Lambda")) * 1.0e3 / variables.at("Mpc_to_m");
}

double MultiUQFFModule::computeUgSum() const {
    double G   = variables.at("G");
    double M   = variables.at("M");
    double r   = variables.at("r");
    double Ug1 = G * M / (r * r);
    double Ug3 = (variables.at("M_ext") > 0.0) ?
        G * variables.at("M_ext") / (variables.at("r_ext") * variables.at("r_ext")) : 0.0;
    double Ug4 = Ug1 * variables.at("f_sc");
    return Ug1 + Ug3 + Ug4;
}

double MultiUQFFModule::computeQuantumPsi() const {
    double unc = std::sqrt(variables.at("Delta_x") * variables.at("Delta_p"));
    return (variables.at("hbar") / unc) *
           variables.at("integral_psi") *
           (2.0 * variables.at("pi") / variables.at("t_Hubble"));
}

double MultiUQFFModule::computeFluidTerm(double g_base) const {
    return variables.at("rho_fluid") * variables.at("V") * g_base;
}

double MultiUQFFModule::computeMsfFactor(double t) const {
    double SFR = variables.at("SFR");
    if (SFR == 0.0) return 0.0;
    double M0  = variables.at("M0");
    if (M0 == 0.0) return 0.0;
    return SFR * (t / variables.at("year_to_s")) / M0;
}

double MultiUQFFModule::computeDMPert() const {
    double r   = variables.at("r");
    double drr = variables.at("delta_rho_over_rho");
    if (!hasVariable("delta_rho_over_rho")) return 0.0;
    double pert = drr + 3.0 * variables.at("G") * variables.at("M") / std::pow(r, 3);
    return (variables.at("M_visible") + variables.at("M_DM")) * pert;
}

// H_res = A_res * sin(2π f_res t)  [resonance term for atomic/quantum systems]
double MultiUQFFModule::computeHRes(double t) const {
    double A_res = variables.at("A_res");
    double f_res = variables.at("f_res");
    return A_res * std::sin(2.0 * variables.at("pi") * f_res * t);
}

double MultiUQFFModule::computeG(double t) {
    variables["t"] = t;

    double Hz        = computeHtz();
    double expansion = 1.0 + Hz * t;
    double sc_corr   = 1.0 - variables["B"] / variables["B_crit"];
    double msf       = computeMsfFactor(t);
    double m_factor  = 1.0 + msf;

    // Expansion of CrabNebula adds Ug2-like kinetic term
    double v_exp = variables.at("v_exp");
    double env_factor = 1.0 + msf;
    if (v_exp > 0.0) {
        double r = variables["r"];
        env_factor += v_exp * v_exp / r;
    }

    double r       = variables["r"];
    double g_base  = (variables["G"] * variables["M"] * m_factor / (r * r))
                     * expansion * sc_corr * env_factor;

    double ug_sum    = computeUgSum();
    double lam_term  = variables["Lambda"] * variables["c"] * variables["c"] / 3.0;
    double quantum   = computeQuantumPsi();
    double fluid     = computeFluidTerm(g_base);

    double dm_pert = 0.0;
    if (hasVariable("delta_rho_over_rho"))
        dm_pert = computeDMPert();

    double h_res = (variables["is_resonance_sys"] != 0.0) ? computeHRes(t) : 0.0;

    return g_base + ug_sum + lam_term + quantum + fluid + dm_pert + h_res;
}

std::string MultiUQFFModule::getEquationText() const {
    std::ostringstream oss;
    oss << "MultiUQFFModule — system: " << current_system << "\n"
        << "g(r,t) = [G M(t)/r²] * (1+H(t,z)) * (1-B/B_crit) * (1+F_env(t))\n"
        << "       + (Ug1 + Ug3' + Ug4) + (Λc²/3)\n"
        << "       + [ℏ/√(ΔxΔp)] ∫ψ_total dV (2π/t_H) + ρ_fluid V g_base\n"
        << "       + DM_pert [ + H_res(t) for atomic systems ]\n"
        << "M(t) = M*(1+SFR*t_yr/M0); F_env uses v_exp (CrabNebula) or wind.\n"
        << "H_res(t) = A_res*sin(2π f_res t)  [HydrogenAtom, HydrogenResonancePToE]\n"
        << "Resonance solutions: See UQFFResonanceValues.h";
    return oss.str();
}
