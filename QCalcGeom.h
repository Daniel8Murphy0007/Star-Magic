/**
 * @file QCalcGeom.h
 * @brief BSFG Geometric Physics Calculator — Public API (v1.0.0)
 *
 * Provides the complete public interface for the Blinking-Star Field Geometry
 * (BSFG) calculator module, along with the three companion number-systems:
 *   - VDS  — Vacuum Density Series   Σ SSq^n / n^26
 *   - DVP  — Dipole Vortex Primes    26!/113, r_q quantization
 *   - BSH  — Buoyancy Series Harmonics U_g2, H_m saturation
 *   - BH26 — BH26 eigenvalue ladder λ_k = k(k+25)
 *
 * Mathematical origin:
 *   BSFG metric:  g_μν = diag(1+ε, −1+ε, r², r²sin²θ)
 *       where ε = η · T_s00(r) · cos(π t_n)
 *       and   T_s00(r) = M_⊙ c² / ((4/3)π r³)  (stellar aether stress energy)
 *
 * Canonical constants are derived from CP4 session tags:
 *   _S147_  — Sessions 147 (DVP/BSH)
 *   _S148_  — Session 148 (BSFG metric, ETA_BSFG, C_FIELD)
 *   _S149_  — Session 149 (BSFG geodesics, horizon, field equations)
 *
 * Requirements boundary:
 *   All public functions must satisfy the 40 tests in qcalcgeom_tests.cpp
 *   (Session 150, commit bdb9244). Run QCALCGEOM::runQCalcGeomTests() to verify.
 *
 * Wolfram WSTP interface:
 *   Symbolic/exact computations are declared in namespace geom_w and
 *   implemented in qcalcgeom_wolfram.h (Phase C).  Guarded by USE_EMBEDDED_WOLFRAM.
 *
 * Integration into MAIN_1_CoAnQi.cpp:
 *   Cosmic-Egg build  : menu case 20 calls runQCalcGeomTests()
 *   Wolfram-only build: menu case 19 calls runQCalcGeomTests()
 *
 * Build:
 *   Requires C++17. No external dependencies beyond standard library.
 *   cl /std:c++17 /O2 /EHsc /D_USE_MATH_DEFINES /I. QCalcGeom.cpp ...
 *
 * Author   : Daniel T. Murphy
 * Created  : Session 150 — March 27, 2026
 * Version  : 1.3.0
 *
 * Revision history:
 *   v1.1.0 — Session 150  : Original 40 BSFG/VDS/DVP/BSH/BH26 tests
 *   v1.2.0 — Session 151  : +neg-buoy +poly26 +UQFF-comp (60 C++ tests)
 *   v1.3.0 — Session 202  : +VDS/DVP/DH26 variant branches + coupling (70 C++ tests)
 */

#ifndef QCALCGEOM_H
#define QCALCGEOM_H

#include <cstdint>
#include <cmath>
#include <string>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

namespace QCALCGEOM {

// ============================================================================
// SECTION 1 — MODULE VERSION
// ============================================================================

constexpr int    QCALCGEOM_VERSION_MAJOR = 1;
constexpr int    QCALCGEOM_VERSION_MINOR = 3;
constexpr int    QCALCGEOM_VERSION_PATCH = 0;
constexpr const char* QCALCGEOM_VERSION_STR = "1.3.0-S202";

// C++ standard gate — matches CMakeLists.txt /std:c++20 project setting
// MSVC reports correct value only with /Zc:__cplusplus; fall back to _MSVC_LANG
#if defined(_MSC_VER)
  static_assert(_MSVC_LANG >= 201703L, "QCalcGeom requires C++17 or later (C++20 preferred)");
#else
  static_assert(__cplusplus >= 201703L, "QCalcGeom requires C++17 or later (C++20 preferred)");
#endif

// ============================================================================
// SECTION 2 — CANONICAL CONSTANTS
// All values must match qcalcgeom_tests.cpp Section 1 exactly.
// Do NOT alter these without a corresponding update to the test file.
// ============================================================================

// Aether-metric coupling η  (_S148_ETA)
constexpr double ETA_BSFG      = 1.0e-22;

// Speed of light c  (_S148_C_LIGHT)
constexpr double C_LIGHT       = 3.0e8;               // m/s

// Solar mass M_⊙  (_S148_MS)
constexpr double M_SUN         = 1.989e30;             // kg

// Solar luminosity L_⊙
constexpr double L_SUN         = 3.828e26;             // W

// Solar radius R_⊙  (_S148_RS)
constexpr double R_SUN         = 6.96e8;               // m

// Gravitational constant G  (_S149_G_N)
constexpr double G_NEWTON      = 6.674e-11;            // m³/(kg·s²)

// Reduced Planck constant ℏ  (_S149_HBAR)
constexpr double HBAR          = 1.055e-34;            // J·s

// Planck constant h  (_S149_H_PL)
constexpr double H_PLANCK      = 6.626e-34;            // J·s

// Boltzmann constant k_B  (_S149_KB)
constexpr double K_BOLTZ       = 1.381e-23;            // J/K

// Planck length l_P  (_S149_LP)
constexpr double L_PLANCK      = 1.616e-35;            // m

// Astronomical unit  (_S149_AU)
constexpr double AU_METERS     = 1.496e11;             // m

// Observed cosmological constant Λ_obs  (_S149_LAM_OBS)
constexpr double LAMBDA_OBS    = 1.1e-52;              // m^{-2}

// DVP prime P_special = 30th prime  (_S147_DVP_PRIME)
constexpr long long DVP_PRIME  = 113LL;

// 26! floating-point approximation  (_S147_FAC26)
constexpr double FAC26_APPROX  = 4.0329146113e26;

// C_num dominant numerator: c²·M_⊙/((4/3)π) ≈ 4.273e46  (_S148_C_FIELD)
constexpr double C_NUM_SOLAR   = 4.273e46;             // m³·kg/m³ · c²

// Default [SSq] = 0.57 — triple-convergence constant (CMB/Kepler/ALMA)
constexpr double SSQ_DEFAULT   = 0.57;

// Einstein coupling κ_E = 8πG/c⁴
// Derived — do NOT hardcode the numerical value in implementations
inline constexpr double kappa_E_value() noexcept {
    return 8.0 * M_PI * G_NEWTON / (C_LIGHT * C_LIGHT * C_LIGHT * C_LIGHT);
}

// ─── Reference values for test comparison (from CP4 #149–#157 docstrings) ──

constexpr double EPS_PRIME_REF  = 5.47e-11;            // |ε′| at R_SUN, t_n=0  [m^{-1}]
constexpr double R_R0R0_REF     = 1.57e-19;            // R^r_{0r0} at R_SUN, t_n=0  [m^{-2}]
constexpr double R_H_BSFG_REF   = 1.62e8;              // Blinking horizon r_h at t_n=1  [m]
constexpr double T_H_BSFG_REF   = 3.37e-12;            // Hawking temperature at r_h  [K]
constexpr double R_CROSS_AU_REF = 0.360;               // Aether–Newton crossover  [AU]
constexpr double H_ETA_REF      = 6.626e-56;           // η · h_Planck  [J·s]
constexpr double AMP_FACTOR_REF = 1.2e4;               // G_00 / (κ_E · T_s00) at R_SUN
constexpr double LAMBDA_EFF_REF = 1.312e-45;           // Λ_eff at R_SUN  [m^{-2}]
constexpr double R_Q_AU_REF     = 0.0973;              // Proplyd quantization radius  [AU]

// BH26 ReRing oscillation centre frequency (_S146_RERING_BB)
constexpr double RERING_BB_HZ   = 1.15e14;            // Hz

// ─── Buoyancy canonical parameters (SOURCE4 compute_Ubi_SOURCE4 formula) ───
// Source: MAIN_1_CoAnQi.cpp SOURCE4 namespace, session 151 neg-buoy integration.
// Ubi = −β_i · Ug_field · Ω_g · M_bh / d_g · wind_mod · U_UA · cos(π·t_n)

constexpr double BETA_I_BSFG    = 0.6;                // β_i — buoyancy coupling constant
constexpr double OMEGA_G_BSFG   = 7.3e-16;            // Ω_g  — galactic spin rate  [rad/s]
constexpr double M_BH_BSFG      = 8.15e36;            // M_bh — canonical BH mass   [kg]
constexpr double D_G_BSFG       = 2.55e20;            // d_g  — GC distance          [m]
constexpr double EPS_SW_BSFG    = 0.001;              // ε_sw — solar wind modulation
constexpr double RHO_SW_BSFG    = 8.0e-21;            // ρ_sw — solar wind density  [kg/m³]
constexpr double U_UA_BSFG      = 1.0;                // U_UA — Aether buoyancy factor

/// Reference Ubi at r=R_SUN, t_n=0 (canonical BSFG, SOURCE4 parameters)
/// = −β_i · (G·M_⊙²/R_⊙²) · (Ω_g·M_bh/d_g) · wind_mod · U_UA · 1.0 ≈ −7.63e33
constexpr double UBI_BSFG_REF   = -7.63e33;

// ─── Session 202: VDS/DVP/DH26 variant-branch reference values ──────────────

/// VDS_prime = ∂/∂z Li_{26}(z)|_{z=0.57} = Li_{25}(0.57)/0.57 ≈ 1.0
/// Sensitivity of the Vacuum Density Series to the [SSq] calibration constant.
constexpr double VDS_PRIME_REF    = 1.0;       // dimensionless

/// BH26 spectral ladder sum Σ_{k=1}^{10} k(k+25) = 1760  (closed form: N(N+1)(2N+1)/6 + 25N(N+1)/2)
constexpr double BH26_SPECTRAL_N10 = 1760.0;  // dimensionless eigenvalue sum

/// Multiplicity of degree-1 spherical harmonics on S^{25}: C(26,25) = 26
constexpr int    BH26_DEG_K1       = 26;       // integer

// ============================================================================
// SECTION 3 — RESULT STRUCTS
// Each struct is the return type of one public API function.
// Field names must remain stable — qcalcgeom_tests.cpp accesses them directly
// by name when the test file is refactored to call QCalcGeom functions.
// ============================================================================

/**
 * @brief BSFG metric and curvature at a given (r, t_n).
 *
 * Covers tests T01–T06 (group BSFG-METRIC).
 *
 * Metric ansatz:  g_00 = 1+ε,  g_rr = −1+ε,  g_ΩΩ = r²
 * where ε = η · T_s00(r) · cos(π t_n).
 */
struct BSFGMetricResult {
    double eps;       ///< ε(r, t_n) — Aether density perturbation                [dimensionless]
    double eps_p;     ///< ε′ = dε/dr                                              [m^{-1}]
    double eps_pp;    ///< ε″ = d²ε/dr²                                            [m^{-2}]
    double A00;       ///< g_00 component = 1 + ε
    double Arr;       ///< g_rr component = −1 + ε
    double R_r0r0;    ///< Riemann component R^r_{0r0} = ε″/2 − (ε′)²/2          [m^{-2}]
    double R_00;      ///< Ricci tensor R_{00} = 3·R^r_{0r0}                       [m^{-2}]
    double R_rr;      ///< Ricci tensor R_{rr}                                     [m^{-2}]
    double R_scalar;  ///< Ricci scalar R = R_00/g^00 + R_rr/g^rr                 [m^{-2}]
    double Kretschner;///< Kretschner invariant K = 12·(R^r_{0r0})²               [m^{-4}]
};

/**
 * @brief BSFG blinking horizon geometry.
 *
 * Covers tests T07–T08 (group BSFG-GEOM).
 *
 * The BSFG horizon exists when A00(r_h)=0, i.e. r_h = (η · C_num)^{1/3}.
 * It is a dynamical (blinking) horizon parametric in t_n.
 */
struct BSFGHorizonResult {
    bool   exists;      ///< true if C_num·η > 0 (horizon physically present)
    double r_h;         ///< Horizon radius r_h = (η·C_num)^{1/3}                  [m]
    double T_H;         ///< Hawking temperature T_H = ℏκ/(2πk_Bc)                [K]
    double kappa_surf;  ///< Surface gravity κ = c² |dA00/dr|_{r_h}/2              [s^{-2}]
    double r_h_over_Rs; ///< Dimensionless r_h / R_⊙
};

/**
 * @brief BSFG Einstein field equation deviations at (r, t_n).
 *
 * Covers tests T28–T29 (group COSMO).
 *
 * G_00 / (κ_E · T_s00) >> 1 signals non-Einstein (Aether-dominated) regime.
 * Λ_eff = κ_E · η · T_s00 / 2 is the effective cosmological constant at radius r.
 */
struct BSFGFieldEqResult {
    double amp_factor;   ///< G_00 / (κ_E · T_s00): GR deviation factor           [dimensionless]
    double Lambda_eff;   ///< Effective Λ at r: κ_E·η·T_s00/2                     [m^{-2}]
    double Lambda_ratio; ///< Λ_eff / Λ_obs
    double rho_vac_eff;  ///< Effective vacuum density: Λ_eff·c²/(8πG)            [kg/m³]
};

/**
 * @brief BSFG geodesic and orbital-quantization results at (r, t_n).
 *
 * Covers tests T09–T10, T34, T39 (group BSFG-GEOM / CHALLENGE).
 *
 * r_cross is where Aether velocity equals Newtonian circular velocity.
 * δJ/J is the fractional correction to orbital angular momentum from ε.
 */
struct BSFGGeodesicResult {
    double v2_newton;       ///< Newtonian circular v² = GM/r                      [m²/s²]
    double v2_aether;       ///< Aether-corrected circular v² at r                 [m²/s²]
    double r_cross_m;       ///< Crossover radius r_cross                           [m]
    double r_cross_AU;      ///< r_cross in AU
    double h_eta;           ///< BSFG Aether action quantum η·h                    [J·s]
    double delta_J_over_J;  ///< |δJ/J| fractional angular-momentum Aether correction
};

/**
 * @brief BSFG holonomy and extra-dimension topology at (r, t_n, loop_area).
 *
 * Covers tests T11–T12, T36 (group BSFG-GEOM / CHALLENGE).
 *
 * The BSFG holonomy group is SO⁺(3,1)×U(1)^{22}, consistent with 26 total
 * dimensions (4 physical + 22 compact U(1) fibres).
 * G2 and Spin(7) holonomies are excluded (both require Ricci-flat manifolds).
 */
struct BSFGHolonomyResult {
    double delta_phi;       ///< Phase accumulation over loop area A: δφ = ω_{0r}·A  [rad]
    double omega_0r;        ///< Off-diagonal metric connection ω_{0r}               [m^{-1}]
    int    n_extra_flat;    ///< Extra flat (U(1)) dimensions = 22 for BSFG
    bool   G2_excluded;     ///< true: holonomy ≠ G₂ (Ricci non-flat)
    bool   Spin7_excluded;  ///< true: holonomy ≠ Spin(7) (Ricci non-flat)
};

/**
 * @brief Vacuum Density Series (VDS) result.
 *
 * Covers tests T13–T16 (group VDS).
 *
 * VDS(SSq, N) = Σ_{n=1}^{N} SSq^n / n^{26}
 * Converges absolutely for |SSq| < 1.  At SSq=0.57, the n=1 term
 * contributes 0.57 / 1 = 0.57 and subsequent terms are < 5×10^{-9}.
 */
struct VDSResult {
    double value;       ///< Σ_{n=1}^{N} SSq^n / n^{26}                           [dimensionless]
    bool   converged;   ///< true if tail bound < 1×10^{-12}
    double tail_bound;  ///< Geometric-series upper bound on truncation error
    int    n_terms_used;///< Terms actually summed (may be < n_terms if converged early)
};

/**
 * @brief Dipole Vortex Primes (DVP) arithmetic result.
 *
 * Covers tests T17–T21 (group DVP).
 *
 * 26! mod 113 = 12 (non-zero by Wilson's theorem, since 113 is prime and 26 < 113).
 * r_q = (2/26!)^{1/26} AU is the proplyd orbital quantization radius.
 */
struct DVPResult {
    long long fac26_mod_113; ///< 26! mod 113 = 12 (computed via exact integer arithmetic)
    bool      non_repeating; ///< true: 26! mod p ≠ 0 (Wilson: p prime, n < p → n! mod p ≠ 0)
    double    r_q_AU;        ///< Proplyd quantization radius (2/26!)^{1/26}        [AU]
    double    r_q_m;         ///< r_q in metres
};

/**
 * @brief Buoyancy Series Harmonics (BSH) result.
 *
 * Covers tests T22–T24 (group BSH).
 *
 * U_g2  = f_Ub · (1 + SSq · Σ_{m=1}^{m_max} (1 − exp(−SSq·m)))
 * H_max = f_Ub · (m_max + 1/2) at saturation (BSH fully developed)
 */
struct BSHResult {
    double U_g2;       ///< Buoyancy harmonic amplitude                            [Hz]
    double H_m_max;    ///< Maximum harmonic H at m = m_max                        [Hz]
    bool   saturated;  ///< true if (1 − exp(−SSq·m_max)) > 1 − 1×10^{−6}
};

/**
 * @brief BH26 Kaluza–Klein eigenvalue and spectral bin for mode k.
 *
 * Covers tests T25–T27 (group BH26).
 *
 * BSFG compactifies 22 extra dimensions on T^{22}.  The KK spectrum on B^{26}
 * gives eigenvalues λ_k = k(k+25) for the 26-sphere Laplacian.
 * The first ladder rung λ_1=26, final checked rung λ_26=1326.
 */
struct BH26Result {
    double lambda_k;    ///< Eigenvalue λ_k = k(k+25)                              [dimensionless]
    double freq_bin_hz; ///< Corresponding spectral bin frequency                  [Hz]
    bool   finite;      ///< true for all k ≥ 1 (all eigenvalues are finite)
};

/**
 * @brief BSFG buoyancy (Ubi) force coupling at (r, t_n).
 *
 * Covers tests T41–T50 (groups NEG-BUOY + NEG-TIME, Session 151).
 *
 * Follows SOURCE4 compute_Ubi_SOURCE4 formula exactly:
 *   Ubi = −β_i · Ug_field · (Ω_g · M_bh / d_g) · wind_mod · U_UA · cos(π·t_n)
 * where:
 *   Ug_field = G · M_⊙² / r²   (BSFG Newtonian self-energy coupling)
 *   wind_mod = 1 + ε_sw · ρ_sw  (solar wind buoyancy modulation ≈ 1)
 *
 * Sign physics (SOURCE4 canonical, _S151_UBI):
 *   t_n = 0  → cos = +1 → Ubi < 0  : buoyancy OPPOSES gravity (normal stabilisation)
 *   t_n = ±1 → cos = −1 → Ubi > 0  : buoyancy AIDS collapse   (negentropic infall phase)
 *   t_n = ±½ → cos =  0 → Ubi ≈ 0  : zero-buoyancy crossover
 *
 * Negative time (t_n < 0): cosine is even → Ubi(t_n) = Ubi(−t_n) exactly.
 * This is the BSFG realisation of source106 NegativeTimeModule time-reversal symmetry.
 */
struct BSFGBuoyancyResult {
    double Ubi;           ///< Buoyancy coupling: −β_i·Ug·orbit·cos(π·t_n)
    double Ug_field;      ///< BSFG gravitational coupling G·M_⊙²/r²             [N-equiv]
    double orbit_factor;  ///< Ω_g·M_bh/d_g · wind_mod · U_UA
    double cos_tn;        ///< cos(π·t_n) — time phase modulator
    bool   negative;      ///< true when Ubi < 0: buoyancy OPPOSES gravity (normal)
    bool   inverted;      ///< true when Ubi > 0: buoyancy AIDS collapse  (negentropic)
    bool   zero_crossing; ///< true when |cos(π·t_n)| < 1e-10             (half-phase)
};

/**
 * @brief Result of the 26th-order polynomial derivative expansion.
 *
 * Covers tests T51–T56 (group POLY26, Session 151 Phase G).
 *
 * Applies the master formula from the 26D projection framework:
 *   d^{26}/dr^{26} (c/r^k) = (k+25)!//(k-1)! · c / r^{k+26}
 *
 * The Pochhammer rising-factorial coefficient k*(k+1)*...*(k+25) is
 * computed exactly as a double product to avoid overflow.
 * At cosmic radii (r ≳ 1 AU) with physical c, the term is < 10^{-280}
 * (negligible) — confirming no singularity at non-zero r.
 */
struct Poly26Result {
    double value;           ///< (k+25)!/(k-1)! * c / r^{k+26}               [SI]
    double factorial_ratio; ///< Pochhammer (k)_{26} = k*(k+1)*...*(k+25)     [dimensionless]
    double r_power;         ///< r^{k+26}                                      [m^{k+26}]
    bool   negligible;      ///< true when |value| < POLY26_NEGLIGIBILITY_THR
};

/// Threshold below which a poly26 term is considered physically negligible
constexpr double POLY26_NEGLIGIBILITY_THR = 1.0e-100;

/**
 * @brief Result of the UQFF compressed-field matrix evaluation.
 *
 * Covers tests T57–T60 (group UQFF-COMP, Session 151 Phase G).
 *
 * The 3×3 UQFF_comp tensor diagonal entries are:
 *   m00 = poly26_derivative(1, G·M_⊙, r).value  (26th r-deriv of U_g base)
 *   m11 = poly26_derivative(26, κ=1, r).value    (26th r-deriv of U_m base)
 *   m22 = 26! * G_N / rho^{27}                   (26th ρ-deriv of U_b = G_N/ρ)
 *
 * Off-diagonal entries are 13th-order cross-coupling estimates:
 *   cross_d13 = sqrt(|poly_d13_ug| * |poly_d13_um|)  (geometric mean coupling)
 *
 * Positive definiteness (eigenvalue_min > 0) confirms that the 26th-order
 * expansion does not introduce negative-energy modes.
 */
struct UQFFCompResult {
    double m00;           ///< Diagonal: 26th r-deriv of gravitational U_g term
    double m11;           ///< Diagonal: 26th r-deriv of magnetic U_m term
    double m22;           ///< Diagonal: 26th ρ-deriv of buoyancy U_b term
    double cross_d13;     ///< Off-diagonal 13th-order cross-coupling estimate
    double eigenvalue_min;///< Minimum diagonal (proxy for smallest eigenvalue)
    bool   positive_definite; ///< true when eigenvalue_min > 0
};

// ============================================================================
// SECTION 3b — SESSION 202 VARIANT-BRANCH RESULT STRUCTS
// Returned by the five new branch-derivation functions added in v1.3.0.
// Mirror the Python dataclasses in QCalcGeom.py (VDSBranchResult, etc.).
// ============================================================================

/**
 * @brief VDS variant-branch derivations: sensitivity, energy density, BH26-coupled.
 *
 * vds_prime      = ∂/∂z Li_{26}(z)|_{z=SSq} = Li_{25}(SSq)/SSq  (calibration sensitivity)
 * vds_density    = Li_{26}(SSq) × RHO_VAC_SCM_BSFG  [J/m³]  (energy-density form)
 * vds_k_weighted = Li_{25}(SSq) + 25·Li_{26}(SSq)            (VDS × BH26-ladder coupling)
 *
 * Tests: T61 (vds_prime ≈ 1.0), T62 (vds_density > 0)
 */
struct VDSBranchResult {
    double vds_li25;       ///< Li_{25}([SSq])  — polylogarithm one degree below VDS
    double vds_prime;      ///< d/dz Li_{26}(z)|_{z=SSq} = Li_{25}/SSq  ≈ 1.0
    double vds_density;    ///< Li_{26}(SSq) × 7.09e-37  [J/m³]
    double vds_k_weighted; ///< Li_{25}(SSq) + 25·Li_{26}(SSq)  — VDS×BH26 coupled amplitude
};

/**
 * @brief DVP variant-branch derivations: spectral sum, pair product, vorticity floor.
 *
 * zeta_sum       = Σ_{p>26, p≤p_max} a(p) where a(p) = SSq^{π(p)} / p^{26}
 * pair_product   = a(29) × a(31)  (double-vortex state coupling)
 * spectral_floor = a(p_max)       (Navier-Stokes vorticity lower bound)
 *
 * Tests: T63 (zeta_sum > 0), T64 (pair_product < a_29²)
 */
struct DVPBranchResult {
    double zeta_sum;        ///< Full DVP prime-vortex spectral sum
    int    n_primes_dvp;    ///< Count of DVP primes in (26, p_max]
    double pair_product;    ///< a(29) × a(31)  — double-vortex amplitude
    double spectral_floor;  ///< a(p_max)  — lowest DVP vortex amplitude
    double a_29;            ///< a(29) — first DVP term (dominant component)
};

/**
 * @brief BH26/DH26 variant-branch derivations: spectral ladder, Casimir, degeneracy, VDS bridge.
 *
 * spectral_sum   = Σ_{k=1}^{N} k(k+25)  [eigenvalue ladder sum; = 1760 for N=10]
 * casimir_energy = ℏ·f_{RR}/2 × Σ 1/λ_k  [J]  (quantum vacuum energy)
 * degeneracy_k1  = 26  (C(26,25): multiplicity of degree-1 on S^{25})
 * vds_coupling   = Σ_{k=1}^{N} λ_k^{-26}  (topology-to-VDS bridge: each eigenvalue raised to VDS power)
 *
 * Tests: T65 (spectral_sum==1760), T66 (casimir_energy>0), T67 (degeneracy_k1==26), T68 (vds_coupling>0)
 */
struct BH26BranchResult {
    double spectral_sum;    ///< Σ_{k=1}^{N} k(k+25)  — eigenvalue ladder sum
    double casimir_energy;  ///< ℏ·RERING_BB_HZ/2 × Σ 1/λ_k  [J]
    int    degeneracy_k1;   ///< 26 = C(26,25)  — multiplicity of degree-1 harmonic on S^{25}
    double vds_coupling;    ///< Σ_{k=1}^{N} λ_k^{-26}  — BH26→VDS topological bridge
    int    N;               ///< Number of eigenvalue levels summed
};

/**
 * @brief VDS×DVP coupled field coefficient and variant (calibration) branch.
 *
 * w_vds          = Li_{26}(SSq) / VDS_max  — normalised VDS weight in [0,1]
 * w_dvp          = zeta_sum / a(29)         — normalised DVP weight
 * joint_coeff    = sqrt(w_vds × w_dvp)      — geometric-mean field coupling
 * variant_branch = |w_vds − w_dvp|          — differential calibration magnitude
 *
 * Test: T69 (joint_coeff ≥ 0)
 */
struct VDSDVPCoupledResult {
    double w_vds;           ///< Normalised VDS weight
    double w_dvp;           ///< Normalised DVP weight
    double joint_coeff;     ///< sqrt(w_vds · w_dvp)  — coupled field amplitude
    double variant_branch;  ///< |w_vds − w_dvp|  — differential calibration gap
};

/**
 * @brief BH26×BSH resonance: BSH evaluated at a BH26 spectral frequency bin.
 *
 * freq_k         = RERING_BB_HZ / λ_k   [Hz]  — BH26 frequency bin
 * bsh_at_k       = BSH U_g2 at omega = 2π·freq_k  (phase-coherent buoyancy)
 * resonance      = bsh_at_k · cos(π·t_n)          (cross-resonance amplitude)
 * energy_density = resonance × RHO_VAC_SCM_BSFG    [J/m³]
 *
 * Test: T70 (energy_density > 0 at t_n=0)
 */
struct BH26BSHResonanceResult {
    double freq_k;          ///< f_k = RERING_BB_HZ / (k(k+25))  [Hz]
    double bsh_at_k;        ///< BSH U_g2 amplitude at BH26 frequency bin
    double resonance;       ///< bsh_at_k · cos(π·t_n)  — cross-resonance
    double energy_density;  ///< resonance × 7.09e-37  [J/m³]
};

// ============================================================================
// SECTION 4 — PUBLIC API DECLARATIONS
// Implemented in QCalcGeom.cpp (Phase B).
// Default parameters use canonical constants from Section 2.
// ============================================================================

/**
 * @brief Compute BSFG metric and curvature at (r, t_n).
 *
 * @param r     Radial coordinate [m] — must be > 0
 * @param t_n   Dimensionless phase parameter (integer or half-integer; cos(π t_n))
 * @return      BSFGMetricResult with all metric and Riemann tensor components
 *
 * Test coverage: T01, T02, T03, T04, T05, T06
 * Reference: CP4 #149 BSFGRiemannCurvatureAetherMetricCalculator
 */
BSFGMetricResult bsfg_metric(double r, double t_n);

/**
 * @brief Compute BSFG blinking horizon properties at phase t_n.
 *
 * @param t_n   Phase parameter (horizon present when |cos(π t_n)| > 0)
 * @return      BSFGHorizonResult with r_h, T_H, κ, and dimensionless ratio
 *
 * Test coverage: T07, T08
 * Reference: CP4 #156 BSFGBlackHoleSolutionHorizonCalculator
 */
BSFGHorizonResult bsfg_horizon(double t_n);

/**
 * @brief Compute BSFG Einstein field equation deviations at (r, t_n).
 *
 * @param r     Radial coordinate [m]
 * @param t_n   Phase parameter
 * @return      BSFGFieldEqResult with amp_factor, Λ_eff, Λ_ratio, ρ_vac_eff
 *
 * Test coverage: T28, T29, T32
 * Reference: CP4 #154 BSFGEinsteinTensorFieldEquationsCalculator
 */
BSFGFieldEqResult bsfg_field_equations(double r, double t_n);

/**
 * @brief Compute BSFG geodesic and orbital quantization at (r, t_n).
 *
 * @param r     Radial coordinate [m]
 * @param t_n   Phase parameter
 * @return      BSFGGeodesicResult with crossover radius, h_η, and δJ/J
 *
 * Test coverage: T09, T10, T34, T39
 * Reference: CP4 #157 BSFGBohrSommerfeldAetherQuantizationCalculator
 */
BSFGGeodesicResult bsfg_geodesic(double r, double t_n);

/**
 * @brief Compute BSFG holonomy and extra-dimension topology at (r, t_n).
 *
 * @param r             Radial coordinate [m]
 * @param t_n           Phase parameter
 * @param loop_area_m2  Area of the parallel-transport loop [m²]
 * @return              BSFGHolonomyResult with δφ, ω_{0r}, n_extra, exclusion bools
 *
 * Test coverage: T11, T12, T36
 * Reference: CP4 #155 BSFGHolonomyGroupTopologyCalculator
 */
BSFGHolonomyResult bsfg_holonomy(double r, double t_n, double loop_area_m2);

/**
 * @brief Evaluate the Vacuum Density Series VDS(SSq, n_terms).
 *
 * VDS(SSq, N) = Σ_{n=1}^{N}  SSq^n / n^{26}
 *
 * @param SSq     Convergence parameter (|SSq| < 1 required; default 0.57)
 * @param n_terms Number of terms to sum before checking convergence (default 200)
 * @return        VDSResult with partial sum, convergence flag, tail bound
 *
 * Test coverage: T13, T14, T15, T16
 * Reference: CP4 #83 (VDS component of ThreeNewNumberSystemsCalculator)
 */
VDSResult vds_series(double SSq = SSQ_DEFAULT, int n_terms = 200);

/**
 * @brief Compute Dipole Vortex Primes arithmetic constants.
 *
 * Computes 26! mod 113, the non-repeating property, and proplyd radius r_q.
 * Uses exact 128-bit integer arithmetic for 26! mod 113 to avoid float rounding.
 *
 * @return DVPResult with exact fac26_mod_113=12, non_repeating=true, r_q values
 *
 * Test coverage: T17, T18, T19, T20, T21
 * Reference: CP4 #83 (DVP component of ThreeNewNumberSystemsCalculator)
 */
DVPResult dvp_arithmetic();

/**
 * @brief Evaluate the Buoyancy Series Harmonics at the given parameters.
 *
 * U_g2 = f_Ub · (1 + SSq · Σ_{m=1}^{m_max} (1 − exp(−SSq · m)))
 *
 * @param f_Ub   Base frequency of Aether buoyancy [Hz] (default 3.3e7)
 * @param SSq    [SSq] convergence parameter (default 0.57)
 * @param omega  Angular frequency of system oscillation [rad/s] (default 2π·f_Ub)
 * @param t_n    Phase parameter (default 0.0)
 * @param m_max  Number of harmonic modes to sum (default 20)
 * @return       BSHResult with U_g2, H_m_max, saturated flag
 *
 * Test coverage: T22, T23, T24
 * Reference: CP4 #83 (BSH component of ThreeNewNumberSystemsCalculator)
 */
BSHResult bsh_harmonic(double f_Ub  = 3.3e7,
                        double SSq   = SSQ_DEFAULT,
                        double omega = 2.0 * M_PI * 3.3e7,
                        double t_n   = 0.0,
                        int    m_max = 20);

/**
 * @brief Return the BH26 Kaluza–Klein eigenvalue for mode k.
 *
 * λ_k = k(k + 25)  for the Laplacian on S^{25} (boundary of B^{26}).
 * Spectral bins map λ_k to frequencies via RERING_BB_HZ / λ_k.
 *
 * @param k   Mode index (k ≥ 1)
 * @return    BH26Result with lambda_k, freq_bin_hz, finite flag
 *
 * Test coverage: T25, T26, T27
 * Reference: CP4 #149 BH26 spectrum section
 */
BH26Result bh26_eigenvalue(int k);

/**
 * @brief Compute the BSFG buoyancy (Ubi) force coupling at (r, t_n).
 *
 * Implements SOURCE4 compute_Ubi_SOURCE4 formula:
 *   Ubi = −β_i · (G·M_⊙²/r²) · (Ω_g·M_bh/d_g) · wind_mod · U_UA · cos(π·t_n)
 *
 * where wind_mod = 1 + ε_sw · ρ_sw and all default parameters match the
 * SOURCE4 canonical constants defined in Section 2.
 *
 * Negative time (t_n < 0): cos(π·t_n) is even → Ubi(t_n) ≡ Ubi(−t_n).
 * Time-reversal symmetry is thus exact — source106 NegativeTimeModule property.
 *
 * @param r        Radial coordinate [m] (must be > 0)
 * @param t_n      Phase parameter (any real value, including t_n < 0)
 * @param beta_i   Buoyancy coupling β_i    (default BETA_I_BSFG = 0.6)
 * @param Omega_g  Galactic spin rate [rad/s] (default OMEGA_G_BSFG = 7.3e-16)
 * @param M_bh     BH mass [kg]              (default M_BH_BSFG = 8.15e36)
 * @param d_g      GC distance [m]           (default D_G_BSFG = 2.55e20)
 * @param epsilon_sw SW modulation factor    (default EPS_SW_BSFG = 0.001)
 * @param rho_sw   SW density [kg/m³]        (default RHO_SW_BSFG = 8e-21)
 * @param U_UA     Aether factor             (default U_UA_BSFG = 1.0)
 * @return         BSFGBuoyancyResult with Ubi, orbit_factor, cos_tn, sign flags
 *
 * Test coverage: T41–T50 (Session 151 Phase E)
 * Reference: SOURCE4 compute_Ubi_SOURCE4; source106 NegativeTimeModule
 */
BSFGBuoyancyResult bsfg_buoyancy(double r, double t_n,
    double beta_i     = BETA_I_BSFG,
    double Omega_g    = OMEGA_G_BSFG,
    double M_bh       = M_BH_BSFG,
    double d_g        = D_G_BSFG,
    double epsilon_sw = EPS_SW_BSFG,
    double rho_sw     = RHO_SW_BSFG,
    double U_UA       = U_UA_BSFG);

/**
 * @brief Compute the 26th-order polynomial derivative of (c / r^k).
 *
 * Implements the closed-form result:
 *   d^{26}/dr^{26} (c/r^k) = (k+25)!/(k-1)! · c / r^{k+26}
 *
 * The coefficient is computed as the Pochhammer product k*(k+1)*...*(k+25)
 * using double arithmetic (exact for k ≤ 40 given double's 53-bit mantissa
 * vs product magnitude).  Requires k ≥ 1 and r > 0.
 *
 * @param k  Power of inverse-r in c/r^k (k ≥ 1)
 * @param c  Scale coefficient
 * @param r  Radial coordinate [m] (must be > 0)
 * @return   Poly26Result with value, factorial_ratio, r_power, negligible flag
 *
 * Test coverage: T51–T56 (group POLY26)
 * Reference: Grok thread grok_share_79fdf5367d1.txt — 26th-order expansions;
 *            iterative differentiation: each step multiplies by -(k+m)/r
 */
Poly26Result poly26_derivative(int k, double c, double r);

/**
 * @brief Evaluate the UQFF compressed-field (UQFF_comp) matrix at (r, rho).
 *
 * Constructs the 3×3 UQFF_comp tensor diagonal and cross-coupling entries
 * using poly26_derivative.  Returns the minimum diagonal as eigenvalue_min.
 *
 * @param r    Radial coordinate [m] (must be > 0)
 * @param rho  Mass-energy density [kg/m³] (must be > 0)
 * @return     UQFFCompResult with m00/m11/m22 diagonals, cross_d13, eigenvalue_min
 *
 * Test coverage: T57–T60 (group UQFF-COMP)
 * Reference: Grok thread grok_share_79fdf5367d1.txt — UQFF_comp matrix
 */
UQFFCompResult uqff_comp_matrix(double r, double rho);

/**
 * @brief Compute VDS variant branches (sensitivity, energy density, BH26-coupled amplitude).
 *
 * @param SSq     Convergence parameter (default SSQ_DEFAULT = 0.57)
 * @param n_terms Polylogarithm truncation (default 200)
 * @return VDSBranchResult with vds_li25, vds_prime, vds_density, vds_k_weighted
 *
 * Tests: T61, T62.  Reference: CP4 #83 VDS + Session 202 derivations.
 */
VDSBranchResult vds_branches(double SSq = SSQ_DEFAULT, int n_terms = 200);

/**
 * @brief Compute DVP spectral sum, pair-product, and vorticity floor.
 *
 * Enumerates primes p > 26 up to p_max; for each prime p at prime-counting
 * position π(p): a(p) = [SSq]^{π(p)} / p^{26}.
 *
 * @param p_max  Upper limit of prime sieve (default 200)
 * @return DVPBranchResult with zeta_sum, n_primes_dvp, pair_product, spectral_floor, a_29
 *
 * Tests: T63, T64.  Reference: CP4 #83 DVP + Navier-Stokes vorticity bound.
 */
DVPBranchResult dvp_branches(int p_max = 200);

/**
 * @brief Compute BH26/DH26 spectral ladder: sum, Casimir energy, degeneracy, VDS bridge.
 *
 * λ_k = k(k+25) for the Laplacian on S^{25}.  Summed from k=1 to N.
 *
 * @param N  Number of eigenvalue levels to sum (default 10)
 * @return BH26BranchResult with spectral_sum, casimir_energy, degeneracy_k1, vds_coupling, N
 *
 * Tests: T65, T66, T67, T68.  Reference: CP4 #149 BH26 + Session 202.
 */
BH26BranchResult bh26_branches(int N = 10);

/**
 * @brief Compute VDS×DVP coupled field coefficient and calibration variant branch.
 *
 * Normalises both series against their respective maxima, then forms the
 * geometric-mean joint coefficient and the differential calibration gap.
 *
 * @param SSq     VDS convergence parameter (default SSQ_DEFAULT)
 * @param p_max   DVP sieve limit (default 200)
 * @param n_terms VDS truncation (default 200)
 * @return VDSDVPCoupledResult with w_vds, w_dvp, joint_coeff, variant_branch
 *
 * Test: T69.  Encodes "many ways to get from one place to the other".
 */
VDSDVPCoupledResult vds_dvp_coupled(double SSq    = SSQ_DEFAULT,
                                     int    p_max   = 200,
                                     int    n_terms = 200);

/**
 * @brief Evaluate BSH harmonics at a BH26 spectral frequency bin.
 *
 * For BH26 mode k: freq_k = RERING_BB_HZ / (k(k+25)) [Hz].
 * The BSH is evaluated with omega = 2π·freq_k at the given phase t_n.
 * Cross-resonance = bsh_at_k × cos(π·t_n).
 *
 * @param f_Ub  Base buoyancy frequency [Hz] (controls BSH H_m amplitude)
 * @param SSq   [SSq] convergence parameter (default SSQ_DEFAULT)
 * @param t_n   Phase parameter (default 0.0)
 * @param k     BH26 mode index k ≥ 1 (default 1)
 * @return BH26BSHResonanceResult with freq_k, bsh_at_k, resonance, energy_density
 *
 * Test: T70.  Reference: BH26 × BSH cross-resonance, Session 202.
 */
BH26BSHResonanceResult bh26_bsh_resonance(double f_Ub = 3.3e7,
                                            double SSq  = SSQ_DEFAULT,
                                            double t_n  = 0.0,
                                            int    k    = 1);

/**
 * @brief Run all 70 QCalcGeom requirements-boundary tests and print results.
 *
 * Phases A-E: 40 original BSFG/VDS/DVP/BSH/BH26/COSMO/CHALLENGE tests (T01-T40)
 * Phase E:    10 Session-151 negative-time + negative-buoyancy tests (T41-T50).
 * Phase G:    10 Session-151 26th-order expansion tests (T51-T60).
 * Phase H202: 10 Session-202 VDS/DVP/DH26 variant-branch + coupling tests (T61-T70).
 *
 * Expected result post Phase H202: 70/70 PASS.
 */
void runQCalcGeomTests();

// ============================================================================
// SECTION 5 — WOLFRAM WSTP SYMBOLIC INTERFACE (Phase C)
// Declared here; implemented in qcalcgeom_wolfram.h + wolfram_sources_bridge.cpp.
// Only compiled when USE_EMBEDDED_WOLFRAM is defined (Wolfram-only and
// Cosmic-Egg builds).
// ============================================================================

#ifdef USE_EMBEDDED_WOLFRAM
namespace geom_w {

/**
 * [W1] Symbolic Riemann tensor for the BSFG metric in variable r_sym.
 * Returns Wolfram Language result string from Simplify[RicciTensor[...]].
 * Maps to: W1 requirement in qcalcgeom_tests.cpp Wolfram block.
 */
std::string bsfg_riemann_symbolic(const std::string& r_sym);

/**
 * [W2] Exact VDS partial sum to 50 significant figures.
 * Returns Wolfram N[Sum[SSq^n/n^26, {n,1,n_terms}], 50].
 */
std::string vds_exact(int n_terms, const std::string& SSq_sym);

/**
 * [W3] Exact 26! mod p via Wolfram BigInteger arithmetic.
 * Returns Wolfram Mod[n!, p] as a string (verified = "12" for n=26, p=113).
 */
std::string dvp_factorial_mod(int n, int p);

/**
 * [W4] Solve BSFG Killing vector equations symbolically.
 * Returns Wolfram Solve[LieDerivative conditions, {K^μ}] result string.
 */
std::string bsfg_killing_vectors();

/**
 * [W5] Numerically integrate BSFG geodesic equations via Wolfram NDSolve.
 * @param r0    Initial radial coordinate [m]
 * @param t_max Integration endpoint in affine parameter
 * @return Wolfram NDSolve trajectory as string
 */
std::string bsfg_geodesic_wstp(double r0, double t_max);

/**
 * [W6] Compute BH26 Kaluza-Klein mass spectrum on R_compact[22].
 * Returns Wolfram T_{22} KK mass spectrum → BH26 λ_k bin string.
 */
std::string kkm_spectrum_wstp();

/**
 * [W7] Symbolic 26th-order derivative of c/r^k via Wolfram FullSimplify[D[...]].
 *
 * Sends to Wolfram:
 *   FullSimplify[D[c/r^k, {r,26}], Assumptions -> r > 0]
 *
 * Expected symbolic result: (k+25)!/(k-1)! * c / r^(k+26)
 * e.g. for k=1, c=1: 26! / r^27  →  "403291461126605635584000000/r^27"
 *
 * @param k  Power of inverse-r (k ≥ 1)
 * @param c  Scale coefficient (passed as WL literal string, e.g. "1", "1*^-22")
 * @return   Wolfram FullSimplify result string
 */
std::string poly26_symbolic(int k, const std::string& c_wl);

} // namespace geom_w
#endif // USE_EMBEDDED_WOLFRAM

} // namespace QCALCGEOM

#endif // QCALCGEOM_H
