// dpm_foundation.h - DPM Foundation Gravity Helpers (canonical implementation)
// DPM IS THE FOUNDATION. NEWTON (GM/r^2) IS THE EMERGENT OUTPUT. NEVER SWAP.
//
// === CANONICAL ONTOLOGY (immutable — Star-Magic.txt) ===
// 0_vacuum -> grad(UA) -> DPM_vortex -> mu_s -> Ug1[seed=DPM]
//   -> Ug_family[Ug1+Ug2+Ug3+Ug4+Ug4i]
//   -> [Ug_family + Um + FUBi + FUBii + UA_uv] -> F_U -> M -> GM/r^2 [LAST]
//
// DPM IS THE FOUNDATION. GM/r^2 IS THE LAST OBSERVABLE PROJECTION.
//
// Canonical Ug1 formula: Ug1 = mu_s * (M/R)   — NO G
//   where mu_s = B * R^3, grad(M_s/r) = M/R (mass gradient without Newton G)
//   G appears only at the FINAL downstream GM/r^2 projection step.

#ifndef DPM_FOUNDATION_H
#define DPM_FOUNDATION_H

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace DPM {

constexpr double G = 6.67430e-11;          // Gravitational constant [m^3/kg/s^2] — DOWNSTREAM ONLY
constexpr double rho_A = 7.09e-37;         // [SCm] vacuum density [J/m^3]
constexpr double rho_UA = 7.09e-36;        // [UA] vacuum density [J/m^3]

// Canonical Ug1 seed: mu_s * grad(M_s/r) — NO G in the seed
// B = magnetic field [T], R = body radius [m], M = mass [kg]
inline double seed_Ug1(double M, double R, double B = 1e-4) {
    if (R <= 0.0) return 0.0;
    double mu_s = B * R * R * R;           // Magnetic moment mu_s = B * R^3
    return mu_s * M / R;                   // grad(M_s/r) = M/R — no Newton G
}

// Backward-compat alias (new code must use seed_Ug1)
inline double emergent_Ug1(double M, double R, double B = 1e-4) {
    return seed_Ug1(M, R, B);
}

// Canonical Ug2 shell: charge-reactivity bubble using vacuum densities — NO G
inline double seed_Ug2(double M, double r, double R, double v_sw = 4e5) {
    if (R <= 0.0 || r <= 0.0) return 0.0;
    double V_body = (4.0 / 3.0) * M_PI * R * R * R;
    double Q_SCm = rho_A * V_body;
    double Q_UA = rho_UA * V_body;
    double E_react = rho_A * v_sw * v_sw / rho_UA;
    double R_b = R * 100.0;
    double S_rb = (r > R_b) ? 1.0 : 0.0;
    return (Q_SCm + Q_UA) * M / (r * r) * S_rb * E_react;
}

inline double emergent_Ug2(double M, double r, double R, double v_sw = 4e5) {
    return seed_Ug2(M, r, R, v_sw);
}

// DPM Ug4: vacuum concentration
inline double seed_Ug4(double rho_v = 7.09e-37, double C_conc = 1e30, double corr_B = 1.0) {
    return rho_v * C_conc * corr_B;
}

inline double emergent_Ug4(double rho_v = 7.09e-37, double C_conc = 1e30, double corr_B = 1.0) {
    return seed_Ug4(rho_v, C_conc, corr_B);
}

// DPM tidal gradient: 3 * Ug1 / r
inline double seed_tidal(double M, double R, double B = 1e-4) {
    return 3.0 * seed_Ug1(M, R, B) / R;
}

inline double emergent_tidal(double M, double R, double B = 1e-4) {
    return seed_tidal(M, R, B);
}

} // namespace DPM

#endif // DPM_FOUNDATION_H
