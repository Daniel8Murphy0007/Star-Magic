// BigBangGravityUQFFModule.cpp
// Doc 38: Cosmic evolution / Big Bang UQFF.
//
// Equation:
//   g_cosmic(t) = g_base(t) + QG_term(t) + DM_term + GW_term(t)
//               + Ug_sum + Λc²/3 + quantum_psi + fluid_term
//
// Unique new physics:
//   QG_term = (ℏc/l_p²) * (t/t_p)              [Planck-scale quantum gravity]
//   DM_term = 0.268 * g_base                    [DM fraction Ω_DM = 0.268]
//   GW_term = h_strain * c²/λ_gw * sin(2πt/t_gw) [GW periodic acceleration]
//   M(t) = M_total * (t/t_Hubble)               [linear cosmic mass assembly]
//   r(t) = c * t                                [Hubble horizon radius]
//   z(t) = t_Hubble/t − 1                       [approximate redshift]
//
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#include "BigBangGravityUQFFModule.h"
#include "UQFFConstants.h"
#include <cmath>
#include <sstream>

BigBangGravityUQFFModule::BigBangGravityUQFFModule() {
    variables["G"]          = UQFF::G;
    variables["c"]          = UQFF::c;
    variables["hbar"]       = UQFF::hbar;
    variables["Lambda"]     = UQFF::Lambda;
    variables["pi"]         = UQFF::pi;
    variables["t_Hubble"]   = UQFF::t_Hubble;
    variables["year_to_s"]  = UQFF::year_to_s;

    // Planck scale
    variables["l_planck"]   = UQFF::l_planck;    // 1.616e-35 m
    variables["t_planck"]   = UQFF::t_planck;    // 5.391e-44 s

    // Total mass of observable universe at Hubble time
    variables["M_total"]    = 1.0e53;             // kg
    variables["H0"]         = UQFF::H0;
    variables["Mpc_to_m"]   = UQFF::Mpc_to_m;
    variables["Omega_m"]    = UQFF::Omega_m;
    variables["Omega_Lambda"]= UQFF::Omega_Lambda;
    variables["Omega_DM"]   = UQFF::Omega_DM;    // 0.268

    // Gravitational wave parameters
    variables["h_strain"]   = UQFF::h_strain_default;  // 1e-21
    variables["lambda_gw"]  = UQFF::lambda_gw_default; // 1e13 m
    variables["t_gw"]       = 1.0e7;              // s (GW period)

    // Quantum ψ_total
    variables["Delta_x"]    = 1.0e-10;
    variables["Delta_p"]    = UQFF::hbar / variables["Delta_x"];
    variables["integral_psi"] = 1.0;

    // Fluid / gas density (CMB era)
    variables["rho_fluid"]  = 1.0e-26;            // kg/m^3 (critical density)
    variables["V"]          = 1.0 / variables["rho_fluid"];

    // EM
    variables["B"]          = 1.0e-9;             // T (primordial field)
    variables["B_crit"]     = UQFF::B_crit_magnetar;
    variables["f_sc"]       = 10.0;
}

// r(t) = c * t  (naive Hubble horizon)
double BigBangGravityUQFFModule::computeR(double t) const {
    return variables.at("c") * t;
}

// M(t) = M_total * (t/t_Hubble)  (linear assembly)
double BigBangGravityUQFFModule::computeM(double t) const {
    return variables.at("M_total") * (t / variables.at("t_Hubble"));
}

// z(t) = t_Hubble / t - 1  (approximate)
double BigBangGravityUQFFModule::computeZ(double t) const {
    if (t <= 0.0) return 1.0e6;
    return variables.at("t_Hubble") / t - 1.0;
}

double BigBangGravityUQFFModule::computeHtz(double z) const {
    return variables.at("H0") * std::sqrt(
        variables.at("Omega_m") * std::pow(1.0 + z, 3) +
        variables.at("Omega_Lambda")) * 1.0e3 / variables.at("Mpc_to_m");
}

// QG_term = (ℏc / l_p²) * (t / t_p)
double BigBangGravityUQFFModule::computeQGTerm(double t) const {
    double lp  = variables.at("l_planck");
    double tp  = variables.at("t_planck");
    return (variables.at("hbar") * variables.at("c") / (lp * lp)) * (t / tp);
}

// GW_term = h_strain * c² / lambda_gw * sin(2π t / t_gw)
double BigBangGravityUQFFModule::computeGWTerm(double t) const {
    double h   = variables.at("h_strain");
    double lam = variables.at("lambda_gw");
    double tgw = variables.at("t_gw");
    return h * variables.at("c") * variables.at("c") / lam *
           std::sin(2.0 * variables.at("pi") * t / tgw);
}

// DM_term = Omega_DM * g_base
double BigBangGravityUQFFModule::computeDMTerm(double g_base) const {
    return variables.at("Omega_DM") * g_base;
}

double BigBangGravityUQFFModule::computeUgSum(double r) const {
    double G = variables.at("G");
    double M = computeM(variables.count("t") ? variables.at("t") : variables.at("t_Hubble"));
    double Ug1 = (r > 0.0) ? G * M / (r * r) : 0.0;
    double Ug4 = Ug1 * variables.at("f_sc");
    return Ug1 + Ug4;
}

double BigBangGravityUQFFModule::computeQuantumPsi() const {
    double unc = std::sqrt(variables.at("Delta_x") * variables.at("Delta_p"));
    return (variables.at("hbar") / unc) *
           variables.at("integral_psi") *
           (2.0 * variables.at("pi") / variables.at("t_Hubble"));
}

double BigBangGravityUQFFModule::computeFluidTerm(double g_base) const {
    return variables.at("rho_fluid") * variables.at("V") * g_base;
}

double BigBangGravityUQFFModule::computeG(double t) {
    variables["t"] = t;

    double r = computeR(t);
    double M = computeM(t);
    double z = computeZ(t);

    double Hz        = computeHtz(z);
    double expansion = 1.0 + Hz * t;
    double sc_corr   = 1.0 - variables["B"] / variables["B_crit"];

    // Base gravity (Newtonian + expansion + SC)
    double g_base = (r > 0.0) ?
        variables["G"] * M / (r * r) * expansion * sc_corr : 0.0;

    // Unique Big Bang terms
    double qg_term  = computeQGTerm(t);
    double gw_term  = computeGWTerm(t);
    double dm_term  = computeDMTerm(g_base);

    double ug_sum   = computeUgSum(r);
    double lam_term = variables["Lambda"] * variables["c"] * variables["c"] / 3.0;
    double quantum  = computeQuantumPsi();
    double fluid    = computeFluidTerm(g_base);

    return g_base + qg_term + gw_term + dm_term + ug_sum + lam_term + quantum + fluid;
}

std::string BigBangGravityUQFFModule::getEquationText() const {
    return
        "g_cosmic(t) = g_base(t) + QG_term(t) + DM_term + GW_term(t)\n"
        "            + (Ug1 + Ug4) + (Λc²/3)\n"
        "            + [ℏ/√(ΔxΔp)] ∫ψ_total dV (2π/t_H) + ρ_fluid V g_base\n"
        "\n"
        "UNIQUE BigBang TERMS (Doc 38):\n"
        "  QG_term = (ℏc/l_p²) * (t/t_p)                   [Planck quantum gravity]\n"
        "  DM_term = Ω_DM * g_base = 0.268 * g_base         [DM fractional gravity]\n"
        "  GW_term = h_strain * c²/λ_gw * sin(2πt/t_gw)    [GW periodic acceleration]\n"
        "\n"
        "COSMIC EVOLUTION:\n"
        "  M(t) = M_total * (t/t_Hubble)    [linear mass assembly,  M_total~10⁵³ kg]\n"
        "  r(t) = c * t                     [Hubble horizon radius]\n"
        "  z(t) = t_Hubble/t − 1            [approximate redshift-time relation]\n"
        "  H(t,z) = H0 sqrt(Ω_m(1+z)³+Ω_Λ)  (Friedmann)";
}
