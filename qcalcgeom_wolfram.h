/**
 * @file qcalcgeom_wolfram.h
 * @brief Phase C — WSTP symbolic bridge for QCALCGEOM::geom_w (W1-W6).
 *
 * Provides the function DEFINITIONS for the six declarations made in
 * QCalcGeom.h Section 5 under #ifdef USE_EMBEDDED_WOLFRAM:
 *
 *   W1  bsfg_riemann_symbolic(r_sym)       Simplify[R^r_{0r0}(r_sym)]
 *   W2  vds_exact(n_terms, SSq_sym)        N[Sum[SSq^n/n^26,{n,1,N}],50]
 *   W3  dvp_factorial_mod(n, p)            Mod[n!, p]   →  "12" for (26,113)
 *   W4  bsfg_killing_vectors()             ∂_r A00 Killing coefficient
 *   W5  bsfg_geodesic_wstp(r0, t_max)      NDSolve radial free-fall
 *   W6  kkm_spectrum_wstp()                BH26 λ_k = k(k+25) spectrum table
 *
 * This header is included exactly once by wolfram_sources_bridge.cpp
 * (appended at the end of the auto-generated file under
 * #ifdef USE_EMBEDDED_WOLFRAM).  The single-include-once pattern plus
 * the header guard below prevents ODR violations.
 *
 * All six functions delegate to WolframEvalToString() which is defined
 * in source174_wolfram_bridge_embedded.cpp and forward-declared here.
 * When the Wolfram kernel is unavailable at runtime each function returns
 * the string "KERNEL_NOT_AVAILABLE" transparently.
 *
 * Constants used match QCalcGeom.h Section 2 canonical values:
 *   ETA_BSFG     = 1e-22      (_S148_ETA)
 *   C_NUM_SOLAR  = 4.273e46   (_S148_C_FIELD)
 *   RERING_BB_HZ = 1.15e14    (_S146_RERING_BB)
 *
 * Author : Daniel T. Murphy
 * Session: 150 | March 27 2026
 * Phase  : C (follows Phase B commit 1a8c188)
 */

#ifndef QCALCGEOM_WOLFRAM_H
#define QCALCGEOM_WOLFRAM_H

#ifdef USE_EMBEDDED_WOLFRAM

#include "QCalcGeom.h"  // Pulls in geom_w declarations + QCALCGEOM namespace

#include <string>
#include <sstream>
#include <iomanip>

// ─────────────────────────────────────────────────────────────────────────────
// WSTP BRIDGE — forward declaration
//
// WolframEvalToString is defined in source174_wolfram_bridge_embedded.cpp.
// It wraps the EvaluatePacket / ToString / InputForm round-trip to the Wolfram
// kernel.  Returns "KERNEL_NOT_AVAILABLE" when the kernel cannot be reached.
// ─────────────────────────────────────────────────────────────────────────────
extern std::string WolframEvalToString(const std::string& code);

// ─────────────────────────────────────────────────────────────────────────────
// CANONICAL BSFG CONSTANTS (WL literal strings, match QCalcGeom.h Section 2)
// ─────────────────────────────────────────────────────────────────────────────
namespace {
    constexpr const char* WL_ETA_BSFG     = "1*^-22";       // ETA_BSFG
    constexpr const char* WL_C_NUM_SOLAR  = "4.273*^46";    // C_NUM_SOLAR
    constexpr const char* WL_RERING_BB    = "1.15*^14";     // RERING_BB_HZ
} // anonymous namespace

namespace QCALCGEOM {
namespace geom_w {

// ─────────────────────────────────────────────────────────────────────────────
// W1 — BSFG Riemann curvature component, symbolic form
//
// Computes R^r_{0r0} = ε''/2 − ε'^2/2 for the BSFG metric perturbation
//
//   ε(r) = η · C_n / r³     (A00(r) = 1 + ε, t_n = 0)
//   ε'(r) = −3η·C_n / r^4
//   ε''(r) = 12η·C_n / r^5
//
// Sends to Wolfram:
//   With[{eta=1*^-22, Cn=4.273*^46},
//     FullSimplify[12*eta*Cn/<r>^5/2 - (-3*eta*Cn/<r>^4)^2/2,
//                  Assumptions -> <r> > 0]]
//
// The result is the exact symbolic curvature coefficient; at r = R_SUN
// it evaluates numerically to ≈ 1.555e-19 m^-2 (reference R_R0R0_REF).
//
// @param r_sym  Wolfram variable name to use (default "r" when empty).
// @return       WL-simplified R^r_{0r0} expression string.
// ─────────────────────────────────────────────────────────────────────────────
std::string bsfg_riemann_symbolic(const std::string& r_sym)
{
    const std::string r = r_sym.empty() ? "r" : r_sym;

    std::ostringstream wl;
    wl << "With[{eta=" << WL_ETA_BSFG << ",Cn=" << WL_C_NUM_SOLAR << "},"
       << "FullSimplify["
       <<   "12*eta*Cn/" << r << "^5/2"
       <<   " - (-3*eta*Cn/" << r << "^4)^2/2,"
       <<   "Assumptions->" << r << ">0"
       << "]]";

    return WolframEvalToString(wl.str());
}

// ─────────────────────────────────────────────────────────────────────────────
// W2 — VDS partial sum to 50 significant figures
//
// Computes the Vector-Dimensional Sum to arbitrary precision:
//
//   VDS_N(SSq) = Sum[ SSq^n / n^26,  {n, 1, n_terms} ]
//
// Sends to Wolfram:
//   N[Sum[(SSq_sym)^n/n^26, {n,1,n_terms}], 50]
//
// For (n_terms=200, SSq_sym="0.57") the 50-digit result is the canonical
// high-precision VDS_200 value used as a convergence benchmark.
//
// @param n_terms   Number of terms in the partial sum.
// @param SSq_sym   Wolfram-literal SSq value (e.g. "0.57" or "57/100").
// @return          50-digit decimal string from N[..., 50].
// ─────────────────────────────────────────────────────────────────────────────
std::string vds_exact(int n_terms, const std::string& SSq_sym)
{
    std::ostringstream wl;
    wl << "N[Sum[(" << SSq_sym << ")^n/n^26,{n,1," << n_terms << "}],50]";
    return WolframEvalToString(wl.str());
}

// ─────────────────────────────────────────────────────────────────────────────
// W3 — Exact n! mod p via Wolfram BigInteger arithmetic
//
// Computes Mod[n!, p] using Wolfram's arbitrary-precision factorial.
// The iterative C++ dvp_arithmetic() uses sequential modular reduction
// which agrees with the exact BigInteger result.
//
// Verified identity: dvp_factorial_mod(26, 113) → "12"
//   ( 26! = 403291461126605635584000000;  403291461126605635584000000 mod 113 = 12 )
//
// Sends to Wolfram:
//   Mod[n!, p]
//
// @param n  Factorial argument (use 26 for DVP canonical check).
// @param p  Modulus prime (use 113 = DVP_PRIME for DVP canonical check).
// @return   String integer result (e.g. "12").
// ─────────────────────────────────────────────────────────────────────────────
std::string dvp_factorial_mod(int n, int p)
{
    std::ostringstream wl;
    wl << "Mod[" << n << "!," << p << "]";
    return WolframEvalToString(wl.str());
}

// ─────────────────────────────────────────────────────────────────────────────
// W4 — BSFG Killing vector metric coefficient (symbolic)
//
// For a static spherically-symmetric metric g = diag{A00, -A00, r², r²sin²θ}
// with A00(r) = 1 + ε(r), the Killing equation ∇_(μ)K_(ν) + ∇_(ν)K_(μ) = 0
// reduces to verifying ∂_r A00 ≠ 0 (non-flat) and that ∂_t A00 = 0 (static).
//
// This function computes ∂_r A00 symbolically — the fundamental Killing field
// coefficient that governs time-translation symmetry in the BSFG geometry:
//
//   ∂_r A00 = ε'(r) = −3·η·C_n / r^4
//
// Sends to Wolfram:
//   With[{eta=1*^-22, Cn=4.273*^46},
//     FullSimplify[D[1 + eta*Cn/r^3, r], Assumptions -> r > 0]]
//
// The result confirms the Killing structure: ∂_t K^t = 0 identically while
// ∂_r K^r ∝ r^{-4} encodes the BSFG radial gradient.
//
// @return  Wolfram-simplified ε'(r) string.
// ─────────────────────────────────────────────────────────────────────────────
std::string bsfg_killing_vectors()
{
    const std::string wl =
        "With[{eta=" + std::string(WL_ETA_BSFG) +
              ",Cn="  + std::string(WL_C_NUM_SOLAR) + "},"
        "FullSimplify["
          "D[1+eta*Cn/r^3,r],"
          "Assumptions->r>0"
        "]]";
    return WolframEvalToString(wl);
}

// ─────────────────────────────────────────────────────────────────────────────
// W5 — BSFG radial geodesic via Wolfram NDSolve
//
// Integrates the radial free-fall geodesic in the BSFG background:
//
//   d²r/dτ² = − ½ g^{rr} ∂_r g_{tt} (dt/dτ)²
//           ≈ − [3·η·C_n / (2·A00(r)·r^4)]    (dt/dτ = 1 normalisation)
//
// where A00(r) = 1 + η·C_n/r³.  Initial conditions r(0)=r0, r'(0)=0 (rest).
//
// Sends to Wolfram:
//   With[{eta=1*^-22, Cn=4.273*^46},
//     Module[{sol},
//       sol = NDSolve[{
//               r''[tau] == -(3*eta*Cn/(2*(1+eta*Cn/r[tau]^3)*r[tau]^4)),
//               r[0] == r0,  r'[0] == 0},
//             r, {tau,0,t_max},
//             MaxSteps->50000, WorkingPrecision->15];
//       ToString[sol[[1,1,2]], InputForm]]]
//
// @param r0     Initial radial coordinate [m]  (e.g. 6.96e8 for R_SUN).
// @param t_max  Integration endpoint in affine parameter τ.
// @return       NDSolve InterpolatingFunction expression as string.
// ─────────────────────────────────────────────────────────────────────────────
std::string bsfg_geodesic_wstp(double r0, double t_max)
{
    std::ostringstream wl;
    wl << std::scientific << std::setprecision(6);
    wl << "With[{eta=" << WL_ETA_BSFG << ",Cn=" << WL_C_NUM_SOLAR << "},"
       <<   "Module[{sol},"
       <<     "sol=NDSolve[{"
       <<       "r''[tau]==-(3*eta*Cn/(2*(1+eta*Cn/r[tau]^3)*r[tau]^4)),"
       <<       "r[0]==" << r0 << ","
       <<       "r'[0]==0"
       <<     "},"
       <<     "r,{tau,0," << t_max << "},"
       <<     "MaxSteps->50000,WorkingPrecision->15"
       <<   "];"
       <<   "ToString[sol[[1,1,2]],InputForm]"
       <<   "]]";
    return WolframEvalToString(wl.str());
}

// ─────────────────────────────────────────────────────────────────────────────
// W6 — BH26 Kaluza-Klein mass spectrum on T²² compact dimension
//
// The BH26 eigenvalue formula λ_k = k(k+25) is the spectrum of the scalar
// Laplacian restricted to the internal S^25 sphere in the 26-dimensional
// Kaluza-Klein tower (26+1D → 4D effective theory).  The corresponding
// physical frequency bin is:
//
//   f_k = RERING_BB_HZ / λ_k   [Hz]
//
// where RERING_BB_HZ = 1.15×10^14 Hz is the BSFG re-ring blackbody frequency
// from Session 146 (_S146_RERING_BB).
//
// Sends to Wolfram:
//   TableForm[
//     Table[{k, k*(k+25), N[1.15*^14/(k*(k+25)), 6]}, {k,1,10}],
//     TableHeadings -> {None, {"k", "lambda_k", "freq_Hz"}}]
//
// Returns the first 10 KK levels as a formatted Wolfram TableForm string.
//
// @return  String containing {k, λ_k, f_k} table for k = 1 … 10.
// ─────────────────────────────────────────────────────────────────────────────
std::string kkm_spectrum_wstp()
{
    // RERING_BB_HZ = 1.15e14 (WL_RERING_BB literal)
    const std::string wl =
        "TableForm["
          "Table[{"
            "k,"
            "k*(k+25),"
            "N[" + std::string(WL_RERING_BB) + "/(k*(k+25)),6]"
          "},{k,1,10}],"
          "TableHeadings->{None,{\"k\",\"lambda_k\",\"freq_Hz\"}}"
        "]";
    return WolframEvalToString(wl);
}

} // namespace geom_w
} // namespace QCALCGEOM

#endif // USE_EMBEDDED_WOLFRAM
#endif // QCALCGEOM_WOLFRAM_H
