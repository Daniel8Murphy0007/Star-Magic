/**
 * UQFF_LA002_Constants.h
 * ─────────────────────────────────────────────────────────────────────────────
 * PhD-Edition Constants for UQFFLearningAssessment_002.cpp
 *
 * Provides the extended physical and UQFF constants required by the PhD
 * Research Edition assessment.  All values sourced from:
 *   - grok_share_c020496d9e.txt (UQFF equations across astrophysical systems)
 *   - Daniel T. Murphy UQFF manuscript (v4.80+)
 *   - 2024 LEP calibration for F_rel, CERN 2025 anyon data, JWST/Gaia DR4
 *
 * Usage:
 *   #include "UQFF_LA002_Constants.h"
 *   using namespace UQFF_LA002;
 *
 * Compatibility: C++17 or later.  No runtime dependencies.
 * Author  : Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date    : 2026
 * Version : 1.0.0
 */

#pragma once
#ifndef UQFF_LA002_CONSTANTS_H
#define UQFF_LA002_CONSTANTS_H

#include <array>
#include <cmath>

namespace UQFF_LA002 {

// ─────────────────────────────────────────────────────────────────────────────
// §1  Mathematical constants
// ─────────────────────────────────────────────────────────────────────────────
static constexpr double PI          = 3.14159265358979323846;
static constexpr double TWO_PI      = 2.0 * PI;
static constexpr double E_EULER     = 2.71828182845904523536;

// ─────────────────────────────────────────────────────────────────────────────
// §2  Standard physical constants (SI)
// ─────────────────────────────────────────────────────────────────────────────
static constexpr double G            = 6.6743e-11;    // N m² kg⁻²  (gravitational)
static constexpr double c_light      = 2.998e8;        // m s⁻¹      (speed of light)
static constexpr double hbar         = 1.0546e-34;     // J s         (reduced Planck)
static constexpr double h_planck     = 6.626e-34;      // J s         (Planck)
static constexpr double Lambda       = 1.1e-52;        // m⁻²         (cosmological Λ)
static constexpr double k_B          = 1.381e-23;      // J K⁻¹       (Boltzmann)
static constexpr double mu_B         = 9.274e-24;      // J T⁻¹       (Bohr magneton)
static constexpr double m_e          = 9.109e-31;      // kg           (electron mass)
static constexpr double m_p          = 1.6726e-27;     // kg           (proton mass)
static constexpr double m_u          = 1.66054e-27;    // kg           (atomic mass unit)

// ─────────────────────────────────────────────────────────────────────────────
// §3  Astrophysical reference scales
// ─────────────────────────────────────────────────────────────────────────────
static constexpr double M_sun        = 1.989e30;       // kg    (solar mass)
static constexpr double R_sun        = 6.957e8;        // m     (solar radius)
static constexpr double pc           = 3.086e16;       // m     (parsec)
static constexpr double kpc          = 3.086e19;       // m     (kiloparsec)
static constexpr double Mpc          = 3.086e22;       // m     (megaparsec)
static constexpr double yr_s         = 3.156e7;        // s     (Julian year)
static constexpr double Myr_s        = 3.156e13;       // s     (megayear)
static constexpr double t_Hubble     = 4.355e17;       // s     (13.8 Gyr)
static constexpr double H_0          = 2.268e-18;      // s⁻¹   (Hubble constant, 70 km/s/Mpc)

// ─────────────────────────────────────────────────────────────────────────────
// §4  UQFF core pipeline constants (standard, shared with _001.cpp)
// ─────────────────────────────────────────────────────────────────────────────
static constexpr double beta_i       = 0.61;           // dimensionless  (buoyancy coupling)
static constexpr double omega_g      = 7.3e-16;        // rad s⁻¹        (UQFF angular freq)
static constexpr double U_UA         = 1.0e-11;        // dimensionless  (vacuum coupling)
static constexpr double f_TRZ        = 0.1;            // dimensionless  (TRZ fraction)
static constexpr double rho_vac_UA   = 7.09e-36;       // J m⁻³          (UA vacuum energy density)
static constexpr double rho_vac_SCm  = 7.09e-37;       // J m⁻³          (SCm vacuum energy density)
static constexpr double B_crit       = 4.4e13;         // T              (QED critical field)

// ─────────────────────────────────────────────────────────────────────────────
// §5  PhD Edition — New UQFF constants (grok_share_c020496d9e.txt)
// ─────────────────────────────────────────────────────────────────────────────

// §5.1  Relativistic calibration (2024 LEP data)
// Source: "F_rel updated to 4.31 × 10^33 N (2024 LEP)" [grok l.267]
static constexpr double F_rel        = 4.31e33;        // N   (relativistic buoyancy calibration)
static constexpr double E_LEP        = 3.204e-8;       // J   (200 GeV = 200 × 1.602e-10 J)
static constexpr double Q_wave       = 1.0e12;         // —   (quantum coupling amplifier)

// §5.2  Superconductive Shell Quotient [SSq]
// Source: "[SSq] = 0.57 (post-calibration)" [grok l.820, session 118]
// Definition: [SSq] = log(ρ_SCm/ρ_UA′) · n · e^{–(π–t)}  → static value = 0.57
static constexpr double SSq          = 0.57;           // —   (vacuum shell entanglement index)

// §5.3  Vacuum decay/magnetic constants
// Source: "γ (5e-5 day^{-1})" [grok, Um discussion]
static constexpr double gamma_decay  = 5.787e-10;      // s⁻¹  (= 5e-5 / 86400, vacuum decay)
static constexpr double delta_k_eta  = 7.25e8;         // —    (nuclear binding energy differential)
static constexpr double k_Ub         = 0.1;            // —    (buoyancy calibration constant)
static constexpr double f_Ub_std     = 2.20e7;         // —    (calibrated f_Ub for star-forming regions)

// §5.4  DPM 4-Component Integral constants
// Source: F_U_Bi_i = ∫₀^{x₂} [...] dx, x₂ ≈ -1.35×10^172 m [grok l.288-296]
static constexpr double F_LENR_ref   = 1.56e36;        // N    (LENR dominant integrand term)
static constexpr double x2_abs_log   = 172.0;          // —    (log₁₀|x₂|; x₂ ≈ −1.35×10^172 m)
static constexpr double DPM_res_W2   = 1.67e3;         // —    (DPM resonance factor, Westerlund 2)
static constexpr double DPM_res_SgrA = 1.67e7;         // —    (DPM resonance factor, Sgr A*)
static constexpr double DPM_res_Pil  = 1.67e7;         // —    (DPM resonance factor, Pillars)
static constexpr double F_UBii_W2    = 2.11e208;       // N    (calibrated DPM integral, W2)
static constexpr double F_UBii_SgrA  = -8.31e211;      // N    (calibrated DPM integral, SgrA*)
static constexpr double F_UBii_Pil   = 2.11e212;       // N    (calibrated DPM integral, Pillars×DPM)

// §5.5  UTe2 Topological Superconductor
// Source: "B_threshold=16T; f_topo=0.1–0.3 topological factor (Andreev STM)" [grok l.4152,4159]
static constexpr double B_thr_UTe2   = 16.0;           // T    (UTe2 superconducting threshold)
static constexpr double f_topo_ref   = 0.20;           // —    (center of 0.1–0.3 range)

// §5.6  Anyon Condensate (CERN 2025)
// Source: "F_UBii,anyons = -F_rel × (E_anyons/E_LEP) × Q_wave × g × exp(-δ_c²/(2σ²))" [grok l.4154]
static constexpr double E_anyons_ref = 1.602e-13;      // J    (1 MeV Ising braiding energy)
static constexpr double delta_c_any  = 1.686;          // —    (anyon collapse threshold, Press–Schechter)
static constexpr double sigma_any    = 1.0;            // —    (non-semisimple TQFT variance)

// §5.7  Hydrogen Resonance — Nuclear Shell Model
// Source: "H_res = A_res sin(2πf_res t) + U_dp × SC_m × k_nuc + S_shell" [grok l.1090]
static constexpr double k_A_hres     = 1.0e-3;         // —    (H_res amplitude coefficient)
static constexpr double A_H          = 1.0;            // —    (hydrogen reference mass number)
static constexpr double k0_nuc       = 0.1;            // —    (nuclear coupling baseline)
// Magic numbers (nuclear shell closures, Mayer–Jensen)
static constexpr std::array<int, 7> MAGIC_NUMBERS = {2, 8, 20, 28, 50, 82, 126};
// Fe-56 peak binding energy reference
static constexpr double E_bind_Fe56  = 8.79e6 * 1.602e-19 * 56.0; // J  (total ~492 MeV)

// §5.8  Buoyancy Harmonics and Vacuum Density Series
// Source: "H_m = Σ_{k=1}^m (1/k)·f_Ub"; "U_g2 = Σ H_m (1–e^{–[SSq]m}) cos(ω t_n)" [grok l.827,214]
// Source: "Vacuum Density Series: Σ_{n=1}^∞ (1/n^26)·[SSq]^n" [grok l.824]
// Partial sum Li₂₆([SSq]) ≈ 0.570 for [SSq]=0.57
static constexpr double Li26_SSq     = 0.570;          // —    (26th polylogarithm partial sum)
static constexpr double V_ratio_std  = 1.0 / 33.0;    // —    (V_little/V_big proto-shell ratio)

// ─────────────────────────────────────────────────────────────────────────────
// §6  Inline utility functions
// ─────────────────────────────────────────────────────────────────────────────

/** Compute dynamic [SSq](n, t_n) at layer n and normalised cosmic time t_n.
 *  [SSq](n,t) = log(ρ_SCm/ρ_UA) · n · e^{–(π–t_n)} */
inline double SSq_dynamic(int n, double t_n) {
    double log_ratio = std::log(rho_vac_SCm / rho_vac_UA);   // negative: log(0.1) ≈ -2.303
    return log_ratio * n * std::exp(-(PI - t_n));
}

/** Compute f_Ub with standard proto-shell volume ratio. */
inline double f_Ub_value(double V_ratio = V_ratio_std) {
    return k_Ub * delta_k_eta * (rho_vac_UA / rho_vac_SCm) * V_ratio;
}

/** Compute UTe2 δ_n series coefficient at layer n.
 *  δ_{n,UTe2} = (2π)^{n/6} · e^{–[SSq]n/26} · (1+f_topo) · e^{–π}  */
inline double delta_n_UTe2(int n, double f_topo = f_topo_ref) {
    double phi_term   = std::pow(TWO_PI, static_cast<double>(n) / 6.0);
    double decay_term = std::exp(-SSq * n / 26.0);
    double topo_boost = 1.0 + f_topo;
    double cosmic_exp = std::exp(-PI);         // e^{–π} baseline (t_n → 0)
    return phi_term * decay_term * topo_boost * cosmic_exp;
}

/** Compute anyon buoyancy force at gravitational field g.
 *  F_UBii,anyons = –F_rel · (E_anyons/E_LEP) · Q_wave · g · exp(–δ_c²/(2σ²)) */
inline double F_UBii_anyons(double g_field, double E_a = E_anyons_ref,
                             double d_c = delta_c_any, double sig = sigma_any) {
    double gauss = std::exp(-d_c * d_c / (2.0 * sig * sig));
    return -F_rel * (E_a / E_LEP) * Q_wave * g_field * gauss;
}

/** Partial sum of Vacuum Density Series Σ_{n=1}^{N} (1/n^26) · [SSq]^n. */
inline double vacuum_density_series(int N_terms, double ssq = SSq) {
    double result = 0.0;
    double ssq_n  = ssq;  // [SSq]^n accumulated
    for (int n = 1; n <= N_terms; ++n) {
        result += ssq_n / std::pow(static_cast<double>(n), 26.0);
        ssq_n  *= ssq;
    }
    return result;
}

/** H_m harmonic coefficient at order m: H_m = Σ_{k=1}^m (1/k) · f_Ub */
inline double H_m_harmonic(int m, double f_ub = f_Ub_std) {
    double H = 0.0;
    for (int k = 1; k <= m; ++k) H += f_ub / static_cast<double>(k);
    return H;
}

} // namespace UQFF_LA002

#endif // UQFF_LA002_CONSTANTS_H
