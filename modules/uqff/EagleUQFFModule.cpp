// EagleUQFFModule.cpp
// Doc 36: Eagle Nebula / Pillars of Creation UQFF.
//
// Key new terms:
//   W_stellar = rho * v_wind^2                          (wind pressure, Eagle form)
//   P_rad     = L * rho / (4*pi*r^2 * c * m_H)         (density-scaled radiation pressure)
//
// General equation:
//   g_Eagle(r,t) = [G M(t)/r²] * (1+H(t,z)) * (1-B/B_crit) * (1+W_stellar+P_rad+M_sf(t))
//               + (Ug1 + Ug3' + Ug4) + (Λc²/3)
//               + quantum ψ_total + ρ_fluid V g_base + DM_pert
//
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#include "EagleUQFFModule.h"
#include "UQFFConstants.h"
#include <cmath>
#include <sstream>

EagleUQFFModule::EagleUQFFModule() {
    variables["G"]          = UQFF::G;
    variables["c"]          = UQFF::c;
    variables["hbar"]       = UQFF::hbar;
    variables["Lambda"]     = UQFF::Lambda;
    variables["pi"]         = UQFF::pi;
    variables["t_Hubble"]   = UQFF::t_Hubble;
    variables["year_to_s"]  = UQFF::year_to_s;
    variables["H0"]         = UQFF::H0;      // 67.15 km/s/Mpc
    variables["Mpc_to_m"]   = UQFF::Mpc_to_m;
    variables["Omega_m"]    = UQFF::Omega_m;
    variables["Omega_Lambda"]= UQFF::Omega_Lambda;
    variables["B_crit"]     = UQFF::B_crit_magnetar;
    variables["m_H"]        = UQFF::m_H;

    // Eagle Nebula / Pillars of Creation (Doc 36)
    variables["M"]          = 5.0e3 * UQFF::M_sun;   // 5000 M☉ (M16 cluster)
    variables["M0"]         = variables["M"];
    variables["r"]          = 3.31e17;                // m (~10.7 pc)
    variables["z"]          = 0.0018;
    variables["t_default"]  = 10.0e6 * UQFF::year_to_s;

    // SFR
    variables["SFR"]        = 0.2 * UQFF::M_sun;     // 0.2 M☉/yr (active region)

    // Stellar wind (NGC6611 = Trapezium of Eagle)
    variables["v_wind"]     = 2.0e6;                  // m/s (2000 km/s OB stars)
    variables["rho_fluid"]  = 1.0e-20;                // kg/m^3

    // Radiation source NGC6611
    variables["L_source"]   = 3.83e33;                // W  (10 L☉ for whole cluster)

    // EM
    variables["B"]          = 2.0e-4;
    variables["f_sc"]       = 10.0;

    // Quantum
    variables["Delta_x"]    = 1.0e-10;
    variables["Delta_p"]    = UQFF::hbar / variables["Delta_x"];
    variables["integral_psi"] = 1.0;
    variables["delta_rho_over_rho"] = 1.0e-5;

    variables["V"]          = 1.0 / variables["rho_fluid"];
    variables["M_visible"]  = variables["M"];
    variables["M_DM"]       = 0.0;
    variables["M_ext"]      = 0.0;
    variables["r_ext"]      = 1.0;
}

void EagleUQFFModule::onVariableUpdated(const std::string& name, double value) {
    if (name == "Delta_x" && value > 0.0)
        variables["Delta_p"] = UQFF::hbar / value;
    else if (name == "M") {
        variables["M0"]       = value;
        variables["M_visible"]= value;
    } else if (name == "rho_fluid" && value > 0.0)
        variables["V"] = 1.0 / value;
}

double EagleUQFFModule::computeHtz() const {
    double z = variables.at("z");
    return variables.at("H0") * std::sqrt(
        variables.at("Omega_m") * std::pow(1.0 + z, 3) +
        variables.at("Omega_Lambda")) * 1.0e3 / variables.at("Mpc_to_m");
}

// W_stellar = rho * v_wind^2  (Eagle form: includes gas density)
double EagleUQFFModule::computeWstellar() const {
    return variables.at("rho_fluid") * variables.at("v_wind") * variables.at("v_wind");
}

// P_rad = L * rho / (4*pi*r^2*c*m_H)  (density-scaled radiation pressure)
double EagleUQFFModule::computePrad() const {
    double L   = variables.at("L_source");
    double rho = variables.at("rho_fluid");
    double r   = variables.at("r");
    return L * rho / (4.0 * variables.at("pi") * r * r * variables.at("c") * variables.at("m_H"));
}

double EagleUQFFModule::computeMsfFactor(double t) const {
    double SFR = variables.at("SFR");
    double M0  = variables.at("M0");
    if (SFR == 0.0 || M0 == 0.0) return 0.0;
    return SFR * (t / variables.at("year_to_s")) / M0;
}

double EagleUQFFModule::computeUgSum() const {
    double G   = variables.at("G");
    double M   = variables.at("M");
    double r   = variables.at("r");
    double Ug1 = G * M / (r * r);
    double Ug3 = (variables.at("M_ext") > 0.0) ?
        G * variables.at("M_ext") / (variables.at("r_ext") * variables.at("r_ext")) : 0.0;
    double Ug4 = Ug1 * variables.at("f_sc");
    return Ug1 + Ug3 + Ug4;
}

double EagleUQFFModule::computeQuantumPsi() const {
    double unc = std::sqrt(variables.at("Delta_x") * variables.at("Delta_p"));
    return (variables.at("hbar") / unc) *
           variables.at("integral_psi") *
           (2.0 * variables.at("pi") / variables.at("t_Hubble"));
}

double EagleUQFFModule::computeFluidTerm(double g_base) const {
    return variables.at("rho_fluid") * variables.at("V") * g_base;
}

double EagleUQFFModule::computeDMPert() const {
    double r   = variables.at("r");
    double drr = variables.at("delta_rho_over_rho");
    double pert = drr + 3.0 * variables.at("G") * variables.at("M") / std::pow(r, 3);
    return (variables.at("M_visible") + variables.at("M_DM")) * pert;
}

double EagleUQFFModule::computeG(double t) {
    variables["t"] = t;

    double Hz        = computeHtz();
    double expansion = 1.0 + Hz * t;
    double sc_corr   = 1.0 - variables["B"] / variables["B_crit"];
    double msf       = computeMsfFactor(t);
    double m_factor  = 1.0 + msf;
    double W_st      = computeWstellar();
    double P_rad     = computePrad();
    double env_factor = 1.0 + W_st + P_rad + msf;

    double r       = variables["r"];
    double g_base  = (variables["G"] * variables["M"] * m_factor / (r * r))
                     * expansion * sc_corr * env_factor;

    double ug_sum  = computeUgSum();
    double lam     = variables["Lambda"] * variables["c"] * variables["c"] / 3.0;
    double quantum = computeQuantumPsi();
    double fluid   = computeFluidTerm(g_base);
    double dm_pert = computeDMPert();

    return g_base + ug_sum + lam + quantum + fluid + dm_pert;
}

std::string EagleUQFFModule::getEquationText() const {
    return
        "g_Eagle(r,t) = [G M(t)/r²] * (1+H(t,z)) * (1-B/B_crit) * (1+W_stellar+P_rad+M_sf(t))\n"
        "             + (Ug1 + Ug3' + Ug4) + (Λc²/3)\n"
        "             + quantum_psi + ρ_fluid V g_base + DM_pert\n"
        "\n"
        "KEY TERMS (Eagle form, Doc 36):\n"
        "  W_stellar = ρ_fluid * v_wind²                     [wind pressure density-coupled]\n"
        "  P_rad     = L * ρ_fluid / (4π r² c m_H)          [density-scaled radiation pressure]\n"
        "  (Eagle distinguishes from Orion: density weighting on both wind and radiation)\n"
        "\n"
        "M(t) = M*(1+SFR*t_yr/M0); H0=67.15 km/s/Mpc; L_NGC6611=3.83e33 W";
}
