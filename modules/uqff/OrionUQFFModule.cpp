// OrionUQFFModule.cpp
// Doc 34: Orion Nebula MUGE — full compressed UQFF implementation.
//
// Equation (compressed):
//   g_Orion(r, t) = [G M(t)/r²] * (1 + H(t,z)) * (1 + W_stellar(t) + P_rad + M_sf(t)) * sc_corr
//                + (Ug1 + Ug3' + Ug4) + (Λ c² / 3)
//                + [ℏ/√(Δx Δp)] * ∫ψ_total dV * (2π/t_H)
//                + ρ_fluid V g_base
//                + (M_vis + M_DM)(δρ/ρ + 3GM/r³)
//
// Approximations:
//   M(t) = M * (1 + SFR * t_yr / M0)                    [star formation growth]
//   W_stellar(t) = v_wind² * (1 + t/t_age)              [time-growing stellar wind]
//   P_rad = L_Trap / (4π r² c m_H)                      [Trapezium radiation pressure per H]
//   Resonant ψ: 2A cos(kx)cos(ωt) + (2π/13.8)A Re[exp(i(kx-ωt))]  [H-alpha wave]
//   sc_corr = 1 - B/B_crit
//
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#include "OrionUQFFModule.h"
#include "UQFFConstants.h"
#include <cmath>
#include <complex>
#include <iomanip>
#include <sstream>

OrionUQFFModule::OrionUQFFModule() {
    // Universal constants
    variables["G"]          = UQFF::G;
    variables["c"]          = UQFF::c;
    variables["hbar"]       = UQFF::hbar;
    variables["Lambda"]     = UQFF::Lambda;
    variables["pi"]         = UQFF::pi;
    variables["t_Hubble"]   = UQFF::t_Hubble;
    variables["year_to_s"]  = UQFF::year_to_s;
    variables["H0"]         = UQFF::H0;
    variables["Mpc_to_m"]   = UQFF::Mpc_to_m;
    variables["Omega_m"]    = UQFF::Omega_m;
    variables["Omega_Lambda"] = UQFF::Omega_Lambda;
    variables["m_H"]        = UQFF::m_H;

    // Orion Nebula system parameters (Doc 34)
    variables["M"]          = 2000.0 * UQFF::M_sun;  // 2000 M☉
    variables["M0"]         = variables["M"];
    variables["r"]          = 1.18e17;                // m  (~3.8 pc radius)
    variables["z"]          = 0.0004;
    variables["t_age"]      = 1.0e6 * UQFF::year_to_s;  // 1 Myr

    // SFR / M_sf
    variables["SFR"]        = 0.1 * UQFF::M_sun;     // kg/yr (0.1 M☉/yr)
    variables["M_visible"]  = variables["M"];
    variables["M_DM"]       = 0.0;

    // Stellar wind
    variables["v_wind"]     = 2.5e5;                  // m/s (250 km/s Trapezium)

    // Trapezium radiation (L_Trap ≈ 1.53e32 W)
    variables["L_Trap"]     = 1.53e32;                // W

    // EM / B field
    variables["B"]          = 1.0e-4;                 // T  (nebula B ~100 μG)
    variables["B_crit"]     = UQFF::B_crit_magnetar;

    // H-alpha wave parameters
    variables["A_wave"]     = 1.0e-10;                // m
    variables["k_wave"]     = 2.0 * UQFF::pi / UQFF::lambda_Halpha;  // m^-1
    variables["omega_wave"] = 2.0 * UQFF::pi * UQFF::c / UQFF::lambda_Halpha;  // rad/s
    variables["x_wave"]     = 0.0;                    // m

    // Fluid / gas
    variables["rho_fluid"]  = 1.67e-20;               // kg/m^3 (H II region ~10 cm^-3)
    variables["V"]          = 1.0 / variables["rho_fluid"];
    variables["delta_rho_over_rho"] = 1.0e-5;

    // Quantum uncertainty
    variables["Delta_x"]    = 1.0e-10;                // m
    variables["Delta_p"]    = UQFF::hbar / variables["Delta_x"];
    variables["integral_psi"] = 1.0;

    // Ug sub-terms
    variables["f_sc"]       = 10.0;                   // Ug4 = Ug1 * f_sc
    variables["M_ext"]      = 0.0;                    // No dominant external body
    variables["r_ext"]      = 1.0;                    // prevent /0
}

void OrionUQFFModule::onVariableUpdated(const std::string& name, double value) {
    if (name == "Delta_x" && value > 0.0)
        variables["Delta_p"] = UQFF::hbar / value;
    else if (name == "M") {
        variables["M0"]       = value;
        variables["M_visible"] = value;
    } else if (name == "rho_fluid" && value > 0.0) {
        variables["V"] = 1.0 / value;
    }
}

// --- private helpers ---

double OrionUQFFModule::computeHtz() const {
    double z   = variables.at("z");
    double Hz  = variables.at("H0") * std::sqrt(
                    variables.at("Omega_m") * std::pow(1.0 + z, 3) +
                    variables.at("Omega_Lambda")) * 1.0e3 / variables.at("Mpc_to_m");
    return Hz;
}

double OrionUQFFModule::computeMsfFactor(double t) const {
    double t_yr = t / variables.at("year_to_s");
    return variables.at("SFR") * t_yr / variables.at("M0");
}

double OrionUQFFModule::computeUgSum() const {
    double G   = variables.at("G");
    double M   = variables.at("M");
    double r   = variables.at("r");
    double Ug1 = G * M / (r * r);
    double Ug3_prime = (variables.at("M_ext") > 0.0) ?
        G * variables.at("M_ext") / (variables.at("r_ext") * variables.at("r_ext")) : 0.0;
    double Ug4 = Ug1 * variables.at("f_sc");
    return Ug1 + Ug3_prime + Ug4;
}

double OrionUQFFModule::computeQuantumPsi() const {
    double dx  = variables.at("Delta_x");
    double dp  = variables.at("Delta_p");
    double unc = std::sqrt(dx * dp);
    return (variables.at("hbar") / unc) *
           variables.at("integral_psi") *
           (2.0 * variables.at("pi") / variables.at("t_Hubble"));
}

double OrionUQFFModule::computeFluidTerm(double g_base) const {
    return variables.at("rho_fluid") * variables.at("V") * g_base;
}

// W_stellar(t) = v_wind² * (1 + t/t_age)
double OrionUQFFModule::computeWstellar(double t) const {
    double vw    = variables.at("v_wind");
    double t_age = variables.at("t_age");
    return vw * vw * (1.0 + t / t_age);
}

// P_rad = L_Trap / (4π r² c m_H)
double OrionUQFFModule::computePrad() const {
    double r   = variables.at("r");
    double L   = variables.at("L_Trap");
    return L / (4.0 * variables.at("pi") * r * r * variables.at("c") * variables.at("m_H"));
}

// Resonant ψ: 2A cos(kx)cos(ωt) + (2π/13.8)A Re[exp(i(kx-ωt))]
double OrionUQFFModule::computeResonantPsi(double t) const {
    double A   = variables.at("A_wave");
    double k   = variables.at("k_wave");
    double om  = variables.at("omega_wave");
    double x   = variables.at("x_wave");
    double standing  = 2.0 * A * std::cos(k * x) * std::cos(om * t);
    std::complex<double> phase = std::exp(std::complex<double>(0.0, k * x - om * t));
    double traveling = (2.0 * variables.at("pi") / 13.8) * A * phase.real();
    return standing + traveling;
}

double OrionUQFFModule::computeDMPert() const {
    double r   = variables.at("r");
    double drr = variables.at("delta_rho_over_rho");
    double pert = drr + 3.0 * variables.at("G") * variables.at("M") / std::pow(r, 3);
    return (variables.at("M_visible") + variables.at("M_DM")) * pert;
}

// --- public compute ---

double OrionUQFFModule::computeG(double t) {
    variables["t"] = t;

    double Hz          = computeHtz();
    double expansion   = 1.0 + Hz * t;
    double sc_corr     = 1.0 - variables["B"] / variables["B_crit"];
    double msf         = computeMsfFactor(t);
    double m_factor    = 1.0 + msf;
    double W_st        = computeWstellar(t);
    double P_rad       = computePrad();
    double env_factor  = 1.0 + W_st + P_rad + msf;

    double r       = variables["r"];
    double g_base  = (variables["G"] * variables["M"] * m_factor / (r * r))
                     * expansion * sc_corr * env_factor;

    double ug_sum      = computeUgSum();
    double lambda_term = variables["Lambda"] * variables["c"] * variables["c"] / 3.0;
    double quantum     = computeQuantumPsi();
    double fluid       = computeFluidTerm(g_base);
    double dm_pert     = computeDMPert();
    double psi_res     = computeResonantPsi(t);   // additive oscillatory contribution

    return g_base + ug_sum + lambda_term + quantum + fluid + dm_pert + psi_res;
}

std::string OrionUQFFModule::getEquationText() const {
    return
        "g_Orion(r, t) = [G M(t)/r²] * (1 + H(t,z)) * (1 + W_st(t) + P_rad + M_sf(t)) * (1 - B/B_crit)\n"
        "              + (Ug1 + Ug3' + Ug4)\n"
        "              + (Λ c² / 3)\n"
        "              + [ℏ/√(Δx Δp)] * ∫ψ_total H ψ_total dV * (2π/t_Hubble)\n"
        "              + ρ_fluid V g_base\n"
        "              + (M_vis + M_DM)(δρ/ρ + 3GM/r³)\n"
        "              + ψ_resonant(t)\n"
        "\n"
        "Where:\n"
        "  M(t)   = M*(1 + SFR*t_yr/M0)                       [star formation]\n"
        "  W_stellar(t) = v_wind²*(1 + t/t_age)               [wind growth]\n"
        "  P_rad  = L_Trap / (4π r² c m_H)                   [Trapezium radiation]\n"
        "  ψ_res  = 2A cos(kx)cos(ωt) + (2π/13.8)A Re[e^{i(kx-ωt)}] [H-alpha]\n"
        "  Ug1 = GM/r²; Ug4 = f_sc * Ug1; Ug3' = GM_ext/r_ext²\n"
        "  H(t,z) = H0 sqrt(Ω_m(1+z)³ + Ω_Λ)  (Friedmann)";
}
