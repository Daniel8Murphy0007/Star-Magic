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
 * Version  : 1.0.0
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
constexpr int    QCALCGEOM_VERSION_MINOR = 0;
constexpr int    QCALCGEOM_VERSION_PATCH = 0;
constexpr const char* QCALCGEOM_VERSION_STR = "1.0.0-S150";

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
 * @brief Run all 40 QCalcGeom requirements-boundary tests and print results.
 *
 * Calls the inline tests from qcalcgeom_tests.cpp when the test file is
 * #included as a translation unit alongside QCalcGeom.cpp.  In standalone
 * mode (QCALCGEOM_STANDALONE defined), qcalcgeom_tests.cpp's own main()
 * calls this instead.
 *
 * Expected result post Phase B: 40/40 PASS.
 * Current result (Phase A): 40/40 PASS (all tests self-contained, no QCalcGeom.cpp needed).
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

} // namespace geom_w
#endif // USE_EMBEDDED_WOLFRAM

} // namespace QCALCGEOM

#endif // QCALCGEOM_H
