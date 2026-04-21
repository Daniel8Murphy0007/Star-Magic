// dpm_emergent.h - DPM-Emergent Gravity Helpers
//
// === CANONICAL ONTOLOGY LOCK (v1) — see also Star-Magic.txt and ARCHITECTURE_FLOW_DIAGRAM.md ===
// 1. Starting state: zero-mass vacuum — rho_UA = 0, rho_vac = |grad(UA)|, F_U(vacuum) = 0.
//    NO MASS exists at quantum cycle start.
// 2. Mass emergence precedes motion. DPM vortical dynamics -> Ug2 shell traps magnetics/
//    spawn material -> mass EMERGES -> only then does Ug1 look gravitational.
// 3. Fixed promotion order: Ug1 promotes the full family — Ug1 -> Ug2 + Ug3 + Ug4 (+ Ug4_i).
// 4. Gravity family is assembled simultaneously: Ug_family = Ug1 + Ug2 + Ug3 + Ug4 (+ Ug4_i).
// 5. Unified field follows: F_U = field(Ug_family, Ub, Um, A, Ui, E_react, t_n).
// 6. Operational modes (Compressed, Resonant, Superconductive, Buoyant) are downstream
//    simultaneous forms of F_U — not independent seed equations.
// 7. GM/r^2 is allowed only as a reduced observational projection AFTER mass emergence
//    and family assembly. It is NOT a seed or foundation term.
// ================================================================================================
//
// Ug1 = k1 * mu_s * grad(M_s/r) * exp(-alpha*t) * cos(pi*t_n) * (1+delta_def)
//   where mu_s = B * R^3 (magnetic dipole moment)
//         grad(M_s/r) = G * M_s / R^2 (mass gradient, NOT Newtonian force)

#ifndef DPM_EMERGENT_H
#define DPM_EMERGENT_H

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace DPM {

constexpr double G = 6.67430e-11;          // Gravitational constant [m^3/kg/s^2]
constexpr double rho_A = 7.09e-37;         // [SCm] vacuum density [J/m^3]
constexpr double rho_UA = 7.09e-36;        // [UA] vacuum density [J/m^3]

// DPM-emergent Ug1: magnetic moment x mass gradient
// B = magnetic field [T], R = body radius [m], M = mass [kg]
// Default B=1e-4 T (typical stellar surface field; override per system)
inline double emergent_Ug1(double M, double R, double B = 1e-4) {
    double mu_s = B * R * R * R;           // Magnetic moment mu_s = B * R^3
    double grad_M = G * M / (R * R);      // Mass gradient grad(M_s/r)
    return mu_s * grad_M;                 // DPM-emergent Ug1
}

// DPM-emergent Ug2: quantum shell trapping (dual charges x reactor energy)
inline double emergent_Ug2(double M, double r, double R, double v_sw = 4e5) {
    double V_body = (4.0 / 3.0) * M_PI * R * R * R;
    double Q_SCm = rho_A * V_body;
    double Q_UA = rho_UA * V_body;
    double E_react = rho_A * v_sw * v_sw / rho_UA;
    double R_b = R * 100.0;
    double S_rb = (r > R_b) ? 1.0 : 0.0;
    return (Q_SCm + Q_UA) * M / (r * r) * S_rb * E_react;
}

// DPM-emergent Ug4: vacuum concentration
inline double emergent_Ug4(double rho_v = 7.09e-37, double C_conc = 1e30, double corr_B = 1.0) {
    return rho_v * C_conc * corr_B;
}

// DPM tidal gradient: 3 * Ug1 / r
inline double emergent_tidal(double M, double R, double B = 1e-4) {
    return 3.0 * emergent_Ug1(M, R, B) / R;
}

} // namespace DPM

#endif // DPM_EMERGENT_H
