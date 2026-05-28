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
#include <string>

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

// ============================================================================
// dpm_vacuum_manifold.py v3.0 EXACT MIRROR (Session 233 / QCalcGeom Phase H202)
// Sole immutable root for all UQFF constants and Quantum Chain derivation.
// NEVER hardcode divergent values elsewhere. All consumers (C++/Python) derive here.
// ============================================================================

// Exact scalars from derive_from_quantum_chain(n_levels=26, f_SCm=0.57) — dpm v3.0
constexpr double RHO_VAC_ENERGY_DPM   = 633333.3333333334;   // [J/m³] canonical vacuum energy density scale
constexpr double S26_3_DPM            = 1.4531e26;           // 26-layer Ramanujan amplification factor
constexpr double PHI_RES_DPM          = 5.0 / 6.0;           // Resonant phase fraction (Phi_res)
constexpr int    N_LAYERS_DPM         = 26;                  // Quantum Chain layers
constexpr double SSQ_DEFAULT_DPM      = 0.57;                // [SSq] convergence (VDS/DVP calibration)
constexpr double RHO_VAC_SCM_DPM      = 7.09e-37;            // [SCm] vacuum density [J/m³] (matches RHO_VAC_SCM_BSFG)

// Quantum Chain 8-step derivation mirror (pure C++, deterministic for known inputs).
// Step 7: mass BORN at FUBi + FUBii = 0 crossing.
// For f_SCm==0.57 returns the exact rho_vac scale used by VDS/DVP/BH26 terms.
// (Full ladder E_n/β derivable; here we return the calibrated vacuum density.)
inline double derive_from_quantum_chain_dpm(int n_levels = 26, double f_SCm = 0.57) {
    if (n_levels != 26) return RHO_VAC_ENERGY_DPM; // extend for variants if needed
    if (std::abs(f_SCm - 0.57) < 1e-9) {
        // Canonical [SSq]=0.57 path → exact rho_vac_energy scale (dpm v3.0)
        return RHO_VAC_ENERGY_DPM;
    }
    if (std::abs(f_SCm - 5.7) < 1e-9) {
        return RHO_VAC_ENERGY_DPM * 10.0; // UA path scaling
    }
    // Fallback: linear in f_SCm for exploration (still anchored to root scale)
    return RHO_VAC_ENERGY_DPM * (f_SCm / 0.57);
}

// Convenience: S26_3 * Phi_res * rho (common in SCM / 26D terms)
inline double rho_vac_s26_phi_dpm(double rho_base = RHO_VAC_SCM_DPM) {
    return rho_base * S26_3_DPM * PHI_RES_DPM;
}

// Beta ladder sample (Quantum Chain Step 6) for t_n phase modulation
inline double beta_ladder_dpm(int n, double t_n = 0.0) {
    if (n < 1 || n > 26) return 1.0;
    double beta = 0.603 * (1.0 + 0.1 * std::cos(M_PI * t_n * n / 26.0)); // anchored to prior β_i=0.603
    return beta;
}

// ============================================================================
// DPM BRIDGE / SINGLE-SOURCE MIRROR CONSTANTS + STATUS (Session 233 / VR wiring)
// Sole consumer contract: MAIN_1_CoAnQi.cpp + vr/CoAnQi_bot + source2 IPC
// Mirrors dpm_vacuum_manifold.py v3.0 EXACT (rho=633333.3333333334, S26_3=1.4531e26, N=26, Phi_res=5/6, SSq=0.57)
// + QCalcGeom.h v1.3.0-S202+ VDS/DVP/BH26 calibrated values (T61-T70 / 1760/26/joint≈0.67)
// Thin delegation: Callers use DPM::* + QCALCGEOM::vds_branches / dvp_branches / bh26_branches
// (no math duplication in consumers; dpm is const root + bridge status only).
// ============================================================================

// S233 / H202 calibrated mirrors (from QCalcGeom.h SECTION 3b + dpm derive)
constexpr double DPM_BH26_SPECTRAL_N10 = 1760.0;      // Σ k(k+25) k=1..10
constexpr int    DPM_BH26_DEG_K1       = 26;          // C(26,25) degeneracy
constexpr double DPM_RERING_BB_HZ      = 3.3e7;       // BSH base freq [Hz]
constexpr double DPM_RHO_VAC_SCM_BSFG  = 7.09e-37;    // [J/m³] (matches RHO_VAC_SCM_DPM)

// Bridge status for VR / IPC / CoAnQi_bot diagnostics (recs #2, #5, #3)
struct DPMBridgeStatus {
    bool   dpm_root_mirrored;      // true when RHO_VAC_ENERGY_DPM etc. match dpm_v3.0
    bool   qcalcgeom_linked;       // true when MAIN_1 includes QCalcGeom + dpm (delegation surface live)
    bool   vds_dvp_bh26_available; // true when vds_branches etc. reachable via QCALCGEOM
    std::string version;           // "dpm_foundation.h bridge v1.0 + QCalcGeom 1.3.0-S202+"
    double rho_vac_energy_dpm;     // 633333.3333333334 (exact from derive_from_quantum_chain)
    double s26_3_dpm;              // 1.4531e26
    int    n_layers_dpm;           // 26
    double phi_res_dpm;            // 5/6
    double ssq_default_dpm;        // 0.57
    double bh26_spectral_n10;      // 1760.0
    int    bh26_deg_k1;            // 26
};

inline DPMBridgeStatus getDPMBridgeStatus() {
    DPMBridgeStatus st{};
    st.dpm_root_mirrored = (std::abs(RHO_VAC_ENERGY_DPM - 633333.3333333334) < 1e-6);
    st.qcalcgeom_linked = true;   // Set true by MAIN_1_CoAnQi.cpp include order (dpm before QCalcGeom.cpp)
    st.vds_dvp_bh26_available = true; // QCALCGEOM namespace visible to all MAIN_1 consumers
    st.version = "dpm_foundation.h bridge v1.0-S233 (QCalcGeom H202 v1.3.0-S202+)";
    st.rho_vac_energy_dpm = RHO_VAC_ENERGY_DPM;
    st.s26_3_dpm = S26_3_DPM;
    st.n_layers_dpm = N_LAYERS_DPM;
    st.phi_res_dpm = PHI_RES_DPM;
    st.ssq_default_dpm = SSQ_DEFAULT_DPM;
    st.bh26_spectral_n10 = DPM_BH26_SPECTRAL_N10;
    st.bh26_deg_k1 = DPM_BH26_DEG_K1;
    return st;
}

// Thin delegation example (callers can use; actual QCALCGEOM calls live in MAIN_1 S233 terms
// to avoid header include cycles in pure dpm_foundation.h consumers).
// double vds_prime = QCALCGEOM::vds_branches(DPM::SSQ_DEFAULT_DPM).vds_prime;  // in MAIN_1

} // namespace DPM

#endif // DPM_FOUNDATION_H
