/**
 * UQFFLearningAssessment_002.cpp
 * ─────────────────────────────────────────────────────────────────────────────
 * UQFF Learning Assessment — PhD Research Edition (v1.0)
 *
 * Systems [0]–[7]:
 *   [0] W2_TRIADIC       — Westerlund 2, DPM 4-Component Buoyancy Integral
 *   [1] SGRA_BUOY        — Sgr A* SMBH, Buoyancy-Dominant DPM
 *   [2] PILLARS_DPM      — Pillars of Creation, DPM Resonance ×10⁷
 *   [3] RESONANCE_26D    — Abstract 26-Layer Vacuum Resonance Framework
 *   [4] UTE2_SC          — UTe2 Topological Superconductor, δ_n Series
 *   [5] ANYONS_CERN      — Anyon Condensate (CERN 2025), Non-unitary TQFT
 *   [6] H_RES_FE56       — Fe-56 Nuclear Shell, H_res Periodic Table Model
 *   [7] VAC_HARMONICS    — Vacuum Density Series + Buoyancy Harmonics U_g2
 *
 * Audience  : PhD candidates in quantum physics, astrophysics, cosmology.
 * Prereq    : UQFFLearningAssessment_001.cpp (Research Student Edition).
 * Source    : grok_share_c020496d9e.txt + Daniel T. Murphy UQFF manuscript v4.80.
 * Author    : Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date      : 2026
 * Depends   : UQFF_LA002_Constants.h  (must be in same include path)
 *
 * Key features (PhD Edition beyond _001.cpp):
 *   1. DPM 4-component buoyancy integral: F_U_Bi_i ≈ ±10^{208..211} N
 *   2. Full 26-layer resonance sum R(t) = Σᵢ₌₁²⁶ [...]
 *   3. UTe2 topological SC δ_n series  (n = 1..9, f_topo = 0.20)
 *   4. Anyon F_UBii with Gaussian collapse exp(−δ_c²/2σ²)
 *   5. H_res nuclear shell model Z=1..118, Fe-56 reference
 *   6. Vacuum Density Series  Σ(1/n²⁶)·[SSq]^n  ≈  Li₂₆([SSq]) ≈ 0.570
 *   7. Buoyancy Harmonics  U_g2 = Σ H_m(1−e^{−[SSq]m})cos(ω t_n)
 *   8. Full U_m form with P_SCm, E_react, f_Heaviside, f_quasi, e^{−[SSq]}
 *   9. Dynamic [SSq](n,t) = log(ρ_SCm/ρ_UA)·n·e^{−(π−t_n)}
 *  10. Scale range: 1 fm (nuclear) → 1 Mpc (cosmic void) = 37 log-decades
 *
 * Integration:
 *   - Header-only (include guard pattern matching _001.cpp)
 *   - #include "UQFFLearningAssessment_002.cpp" from MAIN_1_CoAnQi.cpp
 *   - All UQFF 2.0 Self-Expanding Framework features present:
 *       setEnableLogging(), setDynamicParameter(), getDynamicParameter(),
 *       exportState(), cross_validate<OtherModuleT>()
 *
 * Reference files:
 *   UQFF_LA002_Constants.h      — PhD constants (F_rel, SSq, E_LEP, ...)
 *   UQFF_LA002_SystemParams.md  — Curriculum outline + learning objectives
 *   UQFF_LA002_PhysicsNotes.md  — Full equation reference sheet
 */

#ifndef UQFF_LEARNING_ASSESSMENT_002_H
#define UQFF_LEARNING_ASSESSMENT_002_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <limits>
#include <map>
#include <fstream>
#include <array>
#include <algorithm>
#include <numeric>

#include "UQFF_LA002_Constants.h"   // PhD-edition constants (namespace UQFF_LA002)

static constexpr double UQFF_LA002_PI     = 3.14159265358979323846;
static constexpr double UQFF_LA002_TWO_PI = 6.28318530717958647692;

// ── System index constants ────────────────────────────────────────────────────
static constexpr int SYS_W2_TRIADIC    = 0;  ///< Westerlund 2 — DPM 4-component
static constexpr int SYS_SGRA_BUOY     = 1;  ///< Sgr A* SMBH — buoyancy-dominant
static constexpr int SYS_PILLARS_DPM   = 2;  ///< Pillars of Creation — DPM resonance
static constexpr int SYS_RESONANCE_26D = 3;  ///< Abstract 26-layer resonance
static constexpr int SYS_UTE2_SC       = 4;  ///< UTe2 topological superconductor
static constexpr int SYS_ANYONS_CERN   = 5;  ///< Anyon condensate (CERN 2025)
static constexpr int SYS_H_RES_FE56    = 6;  ///< Fe-56 nuclear shell / H_res
static constexpr int SYS_VAC_HARMONICS = 7;  ///< Vacuum density series + harmonics
static constexpr int NUM_SYSTEMS_002   = 8;

// ─────────────────────────────────────────────────────────────────────────────
class UQFFLearningAssessment002 {
private:

    // ── Shared physical constants ─────────────────────────────────────────────
    const double G            = 6.6743e-11;   // N m² kg⁻²
    const double c_light      = 2.998e8;       // m s⁻¹
    const double hbar         = 1.0546e-34;    // J s
    const double Lambda       = 1.1e-52;       // m⁻²
    const double M_sun        = 1.989e30;      // kg
    const double m_e_kg       = 9.109e-31;     // kg  (electron)
    const double m_p_kg       = 1.6726e-27;    // kg  (proton)
    const double m_u_kg       = 1.66054e-27;   // kg  (atomic mass unit)
    const double mu_B_SI      = 9.274e-24;     // J T⁻¹  (Bohr magneton)
    const double h_planck     = 6.626e-34;     // J s
    const double t_Hubble     = 4.355e17;      // s  (13.8 Gyr)

    // ── UQFF standard buoyancy pipeline ─────────────────────────────────────
    const double beta_i       = 0.61;          // dimensionless
    const double omega_g      = 7.3e-16;       // rad s⁻¹
    const double U_UA         = 1.0e-11;       // dimensionless
    const double f_TRZ        = 0.1;           // TRZ fraction
    const double rho_vac_UA   = 7.09e-36;      // J m⁻³
    const double rho_vac_SCm  = 7.09e-37;      // J m⁻³
    const double B_crit       = 4.4e13;        // T

    // ── PhD edition constants ────────────────────────────────────────────────
    const double F_rel        = 4.31e33;       // N   (2024 LEP calibration)
    const double E_LEP        = 3.204e-8;      // J   (200 GeV)
    const double Q_wave       = 1.0e12;        // coupling amplifier
    const double SSq          = 0.57;          // [SSq] index
    const double gamma_d      = 5.787e-10;     // s⁻¹ (= 5e-5/86400)
    const double delta_k_eta  = 7.25e8;        // nuclear binding differential
    const double k_Ub_val     = 0.1;           // buoyancy calibration
    const double f_Ub_cal     = 2.20e7;        // calibrated f_Ub
    const double delta_c_any  = 1.686;         // anyon collapse threshold
    const double sigma_any    = 1.0;           // TQFT variance
    const double k_A_hres     = 1.0e-3;        // H_res amplitude k_A
    const double A_H_ref      = 1.0;           // hydrogen mass number reference
    const double k0_nuc       = 0.1;           // nuclear coupling baseline

    // ── Assessment metrics ───────────────────────────────────────────────────
    double diversity_score    = 1.0;   // 8/8 physical regimes
    double dynamic_score      = 1.0;   // 10 unique PhD processes
    double scalability_score  = 0.0;   // log-decades / 40 (computed at init)
    double coverage_score     = 1.0;   // 8/8 systems computed

    // ── UQFF 2.0 Self-Expanding Framework ───────────────────────────────────
    bool                        logging_enabled = false;
    std::map<std::string,double> dynamic_params;

    // ─────────────────────────────────────────────────────────────────────────
    // Per-system parameters
    // ─────────────────────────────────────────────────────────────────────────

    // [0] Westerlund 2 — DPM 4-component triadic
    double M_w2t, r_w2t, DPM_res_w2, F_LENR_w2, T_sf_w2;
    double M_ext_w2t, r_ext_w2t;

    // [1] Sgr A* — buoyancy-dominant DPM
    double M_sgra_b, r_sgra_b, DPM_res_sgra, F_LENR_sgra, T_sf_sgra;
    double M_ext_sgra_b, r_ext_sgra_b;

    // [2] Pillars of Creation — DPM resonance ×10⁷
    double M_pil_d, r_pil_d, DPM_res_pil, F_LENR_pil, T_sf_pil;
    double M_ext_pil_d, r_ext_pil_d;

    // [3] 26D resonance layer framework
    double M_res26, r_res26, N_layers_res26, T_sf_res26, SSq_res26;
    double M_ext_res26, r_ext_res26;

    // [4] UTe2 topological superconductor
    double M_ute2, r_ute2, B_ute2, f_topo_ute2;
    double M_ext_ute2, r_ext_ute2;

    // [5] Anyon condensate (CERN 2025)
    double M_any, r_any, E_anyons_loc, delta_c_loc, sigma_loc;
    double M_ext_any, r_ext_any;

    // [6] H_res nuclear shell model — Fe-56
    double M_fe56, r_fe56, Z_fe56, A_fe56, N_fe56;
    double M_ext_fe56, r_ext_fe56;

    // [7] Vacuum density series + buoyancy harmonics
    double M_vac, r_vac, N_harmonics_vac, gamma_vac_loc;
    double M_ext_vac, r_ext_vac;

    // ─────────────────────────────────────────────────────────────────────────
    // Private helpers
    // ─────────────────────────────────────────────────────────────────────────

    double g_base_for(int idx) const {
        switch (idx) {
            case SYS_W2_TRIADIC:    return G * M_w2t    / (r_w2t    * r_w2t);
            case SYS_SGRA_BUOY:     return G * M_sgra_b / (r_sgra_b * r_sgra_b);
            case SYS_PILLARS_DPM:   return G * M_pil_d  / (r_pil_d  * r_pil_d);
            case SYS_RESONANCE_26D: return G * M_res26  / (r_res26  * r_res26);
            case SYS_UTE2_SC:       return G * M_ute2   / (r_ute2   * r_ute2);
            case SYS_ANYONS_CERN:   return G * M_any    / (r_any    * r_any);
            case SYS_H_RES_FE56:    return G * M_fe56   / (r_fe56   * r_fe56);
            case SYS_VAC_HARMONICS: return G * M_vac    / (r_vac    * r_vac);
            default:                return 0.0;
        }
    }

    double get_M(int idx) const {
        switch (idx) {
            case SYS_W2_TRIADIC:    return M_w2t;
            case SYS_SGRA_BUOY:     return M_sgra_b;
            case SYS_PILLARS_DPM:   return M_pil_d;
            case SYS_RESONANCE_26D: return M_res26;
            case SYS_UTE2_SC:       return M_ute2;
            case SYS_ANYONS_CERN:   return M_any;
            case SYS_H_RES_FE56:    return M_fe56;
            case SYS_VAC_HARMONICS: return M_vac;
            default:                return 0.0;
        }
    }

    double get_r(int idx) const {
        switch (idx) {
            case SYS_W2_TRIADIC:    return r_w2t;
            case SYS_SGRA_BUOY:     return r_sgra_b;
            case SYS_PILLARS_DPM:   return r_pil_d;
            case SYS_RESONANCE_26D: return r_res26;
            case SYS_UTE2_SC:       return r_ute2;
            case SYS_ANYONS_CERN:   return r_any;
            case SYS_H_RES_FE56:    return r_fe56;
            case SYS_VAC_HARMONICS: return r_vac;
            default:                return 1.0;
        }
    }

    double get_DPM_res(int idx) const {
        switch (idx) {
            case SYS_W2_TRIADIC:   return DPM_res_w2;
            case SYS_SGRA_BUOY:    return DPM_res_sgra;
            case SYS_PILLARS_DPM:  return DPM_res_pil;
            default:               return 1.0;   // other systems: unscaled
        }
    }

    /// Effective magnetic field per system (for U_m computation)
    double get_B_eff(int idx) const {
        switch (idx) {
            case SYS_W2_TRIADIC:    return 1.0e-5;   // T  typical OB-wind ISM
            case SYS_SGRA_BUOY:     return 1.0e-3;   // T  Galactic Centre
            case SYS_PILLARS_DPM:   return 1.0e-5;   // T  ISM pillar
            case SYS_RESONANCE_26D: return 1.0e-5;   // T  reference
            case SYS_UTE2_SC:       return B_ute2;   // T  UTe2 lab field (16 T)
            case SYS_ANYONS_CERN:   return 1.0;       // T  experiment bench
            case SYS_H_RES_FE56:    return 30.0;      // T  hyperfine at Fe nucleus
            case SYS_VAC_HARMONICS: return 1.0e-9;    // T  IGM background
            default:                return 1.0e-5;
        }
    }

    const char* system_name(int idx) const {
        static const char* names[NUM_SYSTEMS_002] = {
            "Westerlund 2 [DPM Triadic]",
            "Sgr A* [Buoyancy Dominant]",
            "Pillars of Creation [DPM Resonance]",
            "26D Resonance Layer Framework",
            "UTe2 Topological Superconductor",
            "Anyon Condensate [CERN 2025]",
            "Fe-56 Nucleus [H_res Shell Model]",
            "Cosmic Void [Vacuum Harmonics]"
        };
        return (idx >= 0 && idx < NUM_SYSTEMS_002) ? names[idx] : "Unknown";
    }

    const char* regime_label(int idx) const {
        static const char* labels[NUM_SYSTEMS_002] = {
            "OB Star Cluster (6 pc, 30 kM\xe2\x98\x89)",
            "SMBH (1.27\xc3\x9710\xc2\xb9\xc2\xb0 m, 4\xc3\x9710\xe2\x81\xb6 M\xe2\x98\x89)",
            "Star-Forming Pillar (5 pc, ~10\xe2\x81\xb4 M\xe2\x98\x89)",
            "26 Orthogonal Vacuum States (Solar ref.)",
            "Condensed Matter — 1 nm crystal lattice",
            "Quasiparticle — fm scale (nuclear)",
            "Atomic Nucleus — 4.35 fm (Fe-56)",
            "Cosmic Void — 1 Mpc scale"
        };
        return (idx >= 0 && idx < NUM_SYSTEMS_002) ? labels[idx] : "Unknown";
    }

    const char* unique_physics_str(int idx) const {
        static const char* phys[NUM_SYSTEMS_002] = {
            "DPM integral F_U_Bi_i = int[LENR+DPM_mom+DPM_grav+...] dx ~2.11e208 N; "
            "DPM_res=1.67e3; F_LENR=1.56e36 N; x2=-1.35e172 m (quadratic); T_sf=2 Myr.",

            "F_U_Bi_i~-8.31e211 N (F_LENR dominant); DPM_res=1.67e7; NSC tidal; "
            "Kerr spin correction a*=0.3; QPO 20-min; GW back-reaction dOmega/dt.",

            "DPM_res=1.67e7 (UV-enhanced); F_U_Bi_i~2.11e212 N; triadic F_Ubi=9.79e-33 N; "
            "f_z,CGM~1.46e-73 ([SSq]-corrected); JWST NIRCam 2022 pillar density.",

            "R(t)=sum_{i=1}^{26} R_i*cos(omega_i*t); R_i=F_base*(1+SSq)*exp(-SSq*i/26); "
            "omega_i=2pi/T_sf*i*(1+SSq); [SSq]=0.57; 26 orthogonal vacuum quantum states.",

            "delta_{n,UTe2}=(2pi)^{n/6}*exp(-SSq*n/26)*(1+f_topo)*exp(-pi); "
            "B_threshold=16 T; f_topo=0.20; Andreev STM; non-unitary TQFT; n=1..9.",

            "F_UBii,anyons=-F_rel*(E_any/E_LEP)*Q_wave*g*exp(-delta_c^2/(2*sigma^2)); "
            "delta_c=1.686; sigma=1.0; Ising braiding E_any=1 MeV; 75% CERN 2025 alignment.",

            "H_res=A_res*sin(2pi*f_res*t)+U_dp*SC_m*k_nuc+S_shell; Z=26,A=56,N=30; "
            "S_shell=0.20 (near magic 28); Fe-56 binding 8.79 MeV/nucleon (global peak).",

            "U_g2=sum H_m*(1-exp(-SSq*m))*cos(omega*t_n); H_m=sum(1/k)*f_Ub; "
            "V_series=sum(1/n^26)*SSq^n~0.570=Li26(SSq); dynamic SSq(n,t)."
        };
        return (idx >= 0 && idx < NUM_SYSTEMS_002) ? phys[idx] : "";
    }

    void get_ext(int idx, double& M_ext, double& r_ext) const {
        switch (idx) {
            case SYS_W2_TRIADIC:    M_ext = M_ext_w2t;    r_ext = r_ext_w2t;    break;
            case SYS_SGRA_BUOY:     M_ext = M_ext_sgra_b; r_ext = r_ext_sgra_b; break;
            case SYS_PILLARS_DPM:   M_ext = M_ext_pil_d;  r_ext = r_ext_pil_d;  break;
            case SYS_RESONANCE_26D: M_ext = M_ext_res26;  r_ext = r_ext_res26;  break;
            case SYS_UTE2_SC:       M_ext = M_ext_ute2;   r_ext = r_ext_ute2;   break;
            case SYS_ANYONS_CERN:   M_ext = M_ext_any;    r_ext = r_ext_any;    break;
            case SYS_H_RES_FE56:    M_ext = M_ext_fe56;   r_ext = r_ext_fe56;   break;
            case SYS_VAC_HARMONICS: M_ext = M_ext_vac;    r_ext = r_ext_vac;    break;
            default: M_ext = M_sun; r_ext = 1.0e11; break;
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Initialise all system default parameters
    // ─────────────────────────────────────────────────────────────────────────
    void initializeDefaults() {
        // ── [0] Westerlund 2 — DPM triadic ───────────────────────────────────
        M_w2t      = 30000.0 * M_sun;  // 3.0e34 kg
        r_w2t      = 1.89e16;          // m  (6 pc cluster radius)
        DPM_res_w2 = 1.67e3;           // dimensionless
        F_LENR_w2  = 1.56e36;          // N
        T_sf_w2    = 6.307e13;         // s  (2 Myr)
        M_ext_w2t  = 4.0e6 * M_sun;   // SgrA* tidal (8 kpc distance)
        r_ext_w2t  = 2.469e20;         // m  (~8 kpc)

        // ── [1] Sgr A* — buoyancy dominant ───────────────────────────────────
        M_sgra_b     = 4.0e6 * M_sun;    // 7.956e36 kg
        r_sgra_b     = 1.27e10;           // m  (~0.09 AU, G2 cloud periastron)
        DPM_res_sgra = 1.67e7;            // NSC-enhanced
        F_LENR_sgra  = 6.16e39;           // N  (F_LENR >> F_rel at this scale)
        T_sf_sgra    = 3.156e11;          // s  (10 kyr accretion cycle)
        M_ext_sgra_b = 1.0e11 * M_sun;   // NSC mass
        r_ext_sgra_b = 9.461e16;          // m  (~3 pc, NSC half-light radius)

        // ── [2] Pillars of Creation — DPM resonance ───────────────────────────
        M_pil_d    = 10100.0 * M_sun;  // 2.009e34 kg
        r_pil_d    = 4.73e16;          // m  (5 pc pillar height)
        DPM_res_pil = 1.67e7;          // UV-ionisation enhanced
        F_LENR_pil  = 1.56e36;         // N
        T_sf_pil    = 3.156e13;        // s  (1 Myr)
        M_ext_pil_d = 40.0 * M_sun;   // Eta Carinae (η Car) + Trumpler 14
        r_ext_pil_d = 6.17e16;         // m  (~6.5 pc separation)

        // ── [3] 26D resonance layer ───────────────────────────────────────────
        M_res26        = 1.0 * M_sun;    // reference: 1 M☉
        r_res26        = 6.957e8;        // m  (R_sun — reference scale)
        N_layers_res26 = 26.0;           // 26 quantum vacuum states
        T_sf_res26     = 3.156e13;       // s  (1 Myr — resonance period)
        SSq_res26      = SSq;            // 0.57  (calibrated [SSq])
        M_ext_res26    = 4.0e6 * M_sun; // SgrA* (cosmological boundary)
        r_ext_res26    = 1.27e10;        // m  (ISCO scale)

        // ── [4] UTe2 topological superconductor ───────────────────────────────
        M_ute2      = m_e_kg;    // 9.109e-31 kg  (electron mass as reference)
        r_ute2      = 1.0e-9;   // m  (1 nm crystal unit cell)
        B_ute2      = 16.0;     // T  (superconducting threshold)
        f_topo_ute2 = 0.20;     // dimensionless (centre 0.1-0.3 range)
        M_ext_ute2  = m_p_kg;   // proton mass (nuclear context)
        r_ext_ute2  = 1.0e-10;  // m  (Bohr radius, atomic scale)

        // ── [5] Anyon condensate (CERN 2025) ─────────────────────────────────
        M_any        = m_e_kg;         // 9.109e-31 kg
        r_any        = 1.0e-15;        // m  (1 fm nuclear interaction scale)
        E_anyons_loc = 1.602e-13;      // J  (1 MeV Ising braiding energy)
        delta_c_loc  = 1.686;          // Press-Schechter collapse threshold
        sigma_loc    = 1.0;            // non-semisimple TQFT variance
        M_ext_any    = m_p_kg;         // proton (nuclear context)
        r_ext_any    = 1.0e-15;        // m  (fm scale)

        // ── [6] H_res nuclear shell — Fe-56 ──────────────────────────────────
        M_fe56     = 56.0 * m_u_kg;   // 9.299e-26 kg  (56 atomic mass units)
        r_fe56     = 4.35e-15;        // m  (4.35 fm nuclear radius, r = r₀ A^{1/3})
        Z_fe56     = 26.0;            // atomic number
        A_fe56     = 56.0;            // mass number
        N_fe56     = 30.0;            // neutron number (A-Z)
        M_ext_fe56 = m_e_kg;          // electron mass (atomic context)
        r_ext_fe56 = 5.292e-11;       // m  (Bohr radius a₀)

        // ── [7] Vacuum density series + buoyancy harmonics ───────────────────
        M_vac          = 1.0 * M_sun;      // cosmological reference (1 M☉)
        r_vac          = 3.086e22;          // m  (1 Mpc Hubble volume element)
        N_harmonics_vac = 26.0;            // harmonic series order
        gamma_vac_loc  = 5.787e-10;        // s⁻¹  (5×10⁻⁵ day⁻¹)
        M_ext_vac      = 1.0e14 * M_sun;  // galaxy cluster (Coma scale)
        r_ext_vac      = 3.086e23;          // m  (10 Mpc)

        // Compute scalability score
        double r_min   = r_any;         // 1e-15 m  (fm, anyons)
        double r_max1  = r_vac;         // 3.086e22 m (1 Mpc)
        scalability_score = std::min(std::log10(r_max1 / r_min) / 40.0, 1.0);
    }

public:
    // ─────────────────────────────────────────────────────────────────────────
    // Constructor / Destructor
    // ─────────────────────────────────────────────────────────────────────────
    UQFFLearningAssessment002() {
        initializeDefaults();
        if (logging_enabled)
            std::cout << "[LOG] UQFFLearningAssessment002 constructed ("
                      << NUM_SYSTEMS_002 << " systems).\n";
    }
    ~UQFFLearningAssessment002() = default;

    // ─────────────────────────────────────────────────────────────────────────
    // §A  Core UQFF Compute Methods (mirror _001.cpp interface)
    // ─────────────────────────────────────────────────────────────────────────

    /// PhD advancement score (0–100%).
    /// Weighted: diversity × dynamic × scalability × coverage, equal weights.
    double compute_advancement() const {
        double score = (diversity_score + dynamic_score
                        + scalability_score + coverage_score) / 4.0 * 100.0;
        if (logging_enabled)
            std::cout << "[ADVANC] PhD score = " << score << " %\n";
        return score;
    }

    /// g_base = G·M/r² for system idx.
    double compute_g_base(int idx) const {
        return g_base_for(idx);
    }

    /// Tier-1 buoyancy: U_bi = ½ × g_base.
    double compute_Ubi(int idx) const {
        return 0.5 * g_base_for(idx);
    }

    /// Tier-2 F_UBii: standard UQFF oscillatory coupling, scaled by DPM_res.
    /// For PhD systems using specialised methods (26D, UTe2, anyons, H_res,
    /// vacuum), this provides the standardised comparison baseline.
    double compute_F_UBii(int idx, double t) const {
        double g   = g_base_for(idx);
        double M   = get_M(idx);
        double r   = get_r(idx);
        double dpm = get_DPM_res(idx);   // 1.0 for systems [3]–[7]
        return -beta_i * g * omega_g * (M / r) * U_UA
               * std::cos(UQFF_LA002_PI * t) * dpm;
    }

    /// Tier-3: external-body tidal buoyancy F_Ub_i.
    double compute_Ub_i(int idx, double t) const {
        double g = g_base_for(idx);
        double M_ext, r_ext;
        get_ext(idx, M_ext, r_ext);
        return -beta_i * g * omega_g * (M_ext / r_ext) * U_UA
               * std::cos(UQFF_LA002_PI * t);
    }

    /// Log₁₀ of scale range across all 8 systems (fm to Mpc).
    double compute_scale_range() const {
        return std::log10(r_vac / r_any);   // log10(3.086e22 / 1e-15) ≈ 37.5
    }

    // ─────────────────────────────────────────────────────────────────────────
    // §B  PhD-Exclusive Compute Methods
    // ─────────────────────────────────────────────────────────────────────────

    /// §B.1  log₁₀|F_U_Bi_i| for the three DPM-integral systems.
    /// Returns the stored calibrated logarithmic scale of the DPM integral.
    double compute_DPM_log_scale(int idx) const {
        switch (idx) {
            case SYS_W2_TRIADIC:   return 208.32;   // log10(2.11e208)
            case SYS_SGRA_BUOY:    return 211.92;   // log10(8.31e211)
            case SYS_PILLARS_DPM:  return 212.32;   // log10(2.11e212)
            default:               return std::log10(std::abs(compute_F_UBii(idx, 0.0)));
        }
    }

    /// §B.2  Full 26-layer resonance sum R(t).
    /// R(t) = Σᵢ₌₁²⁶ R_{Ug1,i} · cos(ω_{Ug1,i} · t)
    /// R_{Ug1,i} = g_base · (1+SSq) · exp(−SSq·i/26)
    /// ω_{Ug1,i} = (2π/T_sf) · i · (1+SSq)
    double compute_26D_resonance_sum(int idx, double t) const {
        double F_base  = g_base_for(idx);
        double T_sf    = (idx == SYS_RESONANCE_26D) ? T_sf_res26
                       : (idx == SYS_W2_TRIADIC)    ? T_sf_w2
                       : (idx == SYS_PILLARS_DPM)   ? T_sf_pil
                       :                               T_sf_res26;
        double ssq_loc = (idx == SYS_RESONANCE_26D) ? SSq_res26 : SSq;
        double sum     = 0.0;
        for (int i = 1; i <= 26; ++i) {
            double R_i     = F_base * (1.0 + ssq_loc)
                           * std::exp(-ssq_loc * i / 26.0);
            double omega_i = UQFF_LA002_TWO_PI / T_sf * i * (1.0 + ssq_loc);
            sum += R_i * std::cos(omega_i * t);
        }
        if (logging_enabled)
            std::cout << "[26D] R(t=" << t << ") = " << sum << " N\n";
        return sum;
    }

    /// §B.3  UTe2 δ_n series coefficient at layer n (1-indexed).
    /// δ_{n,UTe2} = (2π)^{n/6} × exp(−[SSq]·n/26) × (1+f_topo) × exp(−π)
    double compute_delta_n_UTe2(int n) const {
        double phi_term   = std::pow(UQFF_LA002_TWO_PI, static_cast<double>(n) / 6.0);
        double decay_term = std::exp(-SSq * n / 26.0);
        double topo_boost = 1.0 + f_topo_ute2;
        double cosmic_exp = std::exp(-UQFF_LA002_PI);   // t_n → 0 baseline
        double result     = phi_term * decay_term * topo_boost * cosmic_exp;
        if (logging_enabled)
            std::cout << "[UTe2 delta_" << n << "] = " << result << "\n";
        return result;
    }

    /// §B.4  Anyon condensate F_UBii.
    /// F_UBii,anyons = −F_rel × (E_any/E_LEP) × Q_wave × g(r,t) × exp(−δ_c²/2σ²)
    double compute_F_UBii_anyons(double g_field = -1.0) const {
        if (g_field < 0.0) g_field = g_base_for(SYS_ANYONS_CERN);
        double gauss  = std::exp(-delta_c_loc * delta_c_loc / (2.0 * sigma_loc * sigma_loc));
        double result = -F_rel * (E_anyons_loc / E_LEP) * Q_wave * g_field * gauss;
        if (logging_enabled)
            std::cout << "[Anyons] F_UBii = " << result << " N (g=" << g_field << ")\n";
        return result;
    }

    /// §B.5  H_res nuclear shell amplitude (DC static approximation).
    /// Full H_res oscillates at nuclear transition frequency; here we return
    /// the time-averaged DC contribution: A_res + U_dp·SC_m·k_nuc + S_shell.
    /// Z = atomic number, A = mass number, N_n = neutron number.
    double compute_H_res(double Z, double A, double N_n) const {
        // Pairing delta: even-even → +0.1; odd-odd → -0.1; mixed → 0.0
        bool Z_even = (static_cast<int>(Z) % 2 == 0);
        bool N_even = (static_cast<int>(N_n) % 2 == 0);
        double delta_pair = (Z_even && N_even) ? 0.1 : ((!Z_even) && (!N_even)) ? -0.1 : 0.0;

        // Shell correction via proximity to magic numbers
        constexpr int magic[7] = {2, 8, 20, 28, 50, 82, 126};
        double S_shell = 0.0;
        for (int m : magic) {
            if (std::abs(static_cast<int>(Z) - m) <= 2)   S_shell += 0.1;
            if (std::abs(static_cast<int>(N_n) - m) <= 2) S_shell += 0.1;
        }

        double A_res  = k_A_hres * Z * (A / A_H_ref) * (1.0 + delta_pair);
        double k_nuc  = k0_nuc * (N_n / Z) * (1.0 + delta_pair);
        double U_dp   = k0_nuc * (A * A_H_ref) / std::max(A * A * 0.09, 1.0e-30);
        double SC_m   = 1.0;   // exact
        double H_result = A_res + U_dp * SC_m * k_nuc + S_shell;
        if (logging_enabled)
            std::cout << "[H_res Z=" << (int)Z << ",A=" << (int)A
                      << "] A_res=" << A_res << " S_shell=" << S_shell
                      << " -> H_res=" << H_result << "\n";
        return H_result;
    }

    /// §B.6  U_g2 buoyancy harmonics series (truncated at M_terms).
    /// U_g2 = Σ_{m=1}^{M} H_m × (1 − e^{−[SSq]m}) × cos(ω_{Ug2} × t_n)
    /// H_m  = Σ_{k=1}^m (1/k) × f_Ub
    /// t_n  = t_seconds / t_Hubble
    double compute_Ug2_harmonics(double t_seconds, int M_terms = 26) const {
        double t_n       = t_seconds / t_Hubble;
        double omega_ug2 = UQFF_LA002_TWO_PI / t_Hubble;   // 1.44e-17 rad/s
        double f_ub      = f_Ub_cal;   // 2.20e7 (calibrated)
        double result    = 0.0;
        double H_m       = 0.0;
        for (int m = 1; m <= M_terms; ++m) {
            H_m    += f_ub / static_cast<double>(m);
            double damp = 1.0 - std::exp(-SSq * m);
            result += H_m * damp * std::cos(omega_ug2 * t_n);
        }
        if (logging_enabled)
            std::cout << "[Ug2 harmonics, M=" << M_terms
                      << ", t_n=" << t_n << "] = " << result << "\n";
        return result;
    }

    /// §B.7  Vacuum Density Series (partial sum to N_terms).
    /// V_series = Σ_{n=1}^{N} (1/n²⁶) × [SSq]^n  →  Li₂₆([SSq]) ≈ 0.570
    double compute_Vacuum_Series(int N_terms = 26) const {
        double result = 0.0;
        double ssq_n  = SSq;
        for (int n = 1; n <= N_terms; ++n) {
            result += ssq_n / std::pow(static_cast<double>(n), 26.0);
            ssq_n  *= SSq;
        }
        if (logging_enabled)
            std::cout << "[VacSeries N=" << N_terms << "] = " << result << "\n";
        return result;
    }

    /// §B.8  Full U_m (single-j approximation with Heaviside correction).
    /// U_m = (μ_j/r_j) × (1 − e^{−γt}·cos(πt_n)) × P_SCm × E_react
    ///       × (1 + 10¹³·f_H) × (1+f_quasi) × e^{−[SSq]}
    double compute_Um_full(int idx, double t) const {
        double M    = get_M(idx);
        double r    = get_r(idx);
        double B_eff = get_B_eff(idx);
        double t_n  = t / t_Hubble;

        // Calibrated μ_j: Bohr magneton scale, field-enhanced
        double mu_j = mu_B_SI * B_eff / std::max(hbar * 1.0e-12, 1.0e-50);
        double gamma_eff = (idx == SYS_VAC_HARMONICS) ? gamma_vac_loc : gamma_d;

        double exp_term  = 1.0 - std::exp(-gamma_eff * t) * std::cos(UQFF_LA002_PI * t_n);
        double P_SCm     = 1.0;
        double E_react   = U_UA;    // reaction energy proxy
        double f_Heaviside = 0.01;
        double f_quasi   = 0.01;

        double Um_val = (mu_j / r) * exp_term * P_SCm * E_react
                        * (1.0 + 1.0e13 * f_Heaviside) * (1.0 + f_quasi)
                        * std::exp(-SSq);
        if (logging_enabled)
            std::cout << "[Um sys[" << idx << "] t=" << t << "] = " << Um_val << "\n";
        return Um_val;
    }

    /// §B.9  Dynamic [SSq](n, t_seconds): vacuum entanglement entropy at layer n.
    /// [SSq](n,t) = log(ρ_SCm/ρ_UA) × n × e^{−(π−t_n)}
    double compute_SSq_dynamic(int n, double t_seconds) const {
        double t_n      = t_seconds / t_Hubble;
        double log_ratio = std::log(rho_vac_SCm / rho_vac_UA);   // ≈ -2.303
        return log_ratio * n * std::exp(-(UQFF_LA002_PI - t_n));
    }

    // ─────────────────────────────────────────────────────────────────────────
    // §C  Print / Display Methods
    // ─────────────────────────────────────────────────────────────────────────

    /// Print the PhD-level student guide.
    void printPhDGuide() const {
        std::cout << "\n";
        std::cout << "================================================================================\n";
        std::cout << "  UQFF Learning Assessment 002 — PhD Research Edition (v1.0)\n";
        std::cout << "  Audience: PhD quantum physics / astrophysics / cosmology\n";
        std::cout << "================================================================================\n\n";

        std::cout << "FRAMEWORK OVERVIEW\n";
        std::cout << "------------------\n";
        std::cout << "  This assessment extends the Research Student Edition (_001)\n";
        std::cout << "  to 8 systems spanning 37+ log-decades of scale (1 fm -> 1 Mpc).\n";
        std::cout << "  10 unique PhD-level physics processes are evaluated:\n\n";
        std::cout << "    1. DPM 4-component buoyancy integral (F_U_Bi_i ~ 10^208..211 N)\n";
        std::cout << "    2. 26D resonance sum R(t) = sum_{i=1}^{26} R_i cos(omega_i t)\n";
        std::cout << "    3. UTe2 topological SC delta_n series (n=1..9, f_topo=0.20)\n";
        std::cout << "    4. Anyon F_UBii: Gaussian collapse exp(-delta_c^2/2sigma^2)\n";
        std::cout << "    5. H_res nuclear shell model (Fe-56 peak, magic numbers)\n";
        std::cout << "    6. Vacuum Density Series: sum(1/n^26)[SSq]^n = Li26([SSq])\n";
        std::cout << "    7. Buoyancy Harmonics: U_g2 = sum H_m (1-e^{-[SSq]m}) cos(...)\n";
        std::cout << "    8. Full Um: P_SCm x E_react x f_Heaviside x f_quasi x e^{-SSq}\n";
        std::cout << "    9. Dynamic [SSq](n,t) = log(rho_SCm/rho_UA) x n x e^{-(pi-t_n)}\n";
        std::cout << "   10. F_rel / E_LEP / Q_wave relativistic calibration (2024 LEP)\n\n";

        std::cout << "SYSTEM CATALOGUE\n";
        std::cout << "----------------\n";
        std::cout << std::left << std::setw(4)  << "Idx"
                  << std::setw(42) << "Name"
                  << std::setw(20) << "g_base (m/s^2)"
                  << std::setw(24) << "Primary PhD Method\n";
        std::cout << std::string(90, '-') << "\n";
        const char* methods[NUM_SYSTEMS_002] = {
            "DPM integral log10|F|", "DPM integral log10|F|", "DPM integral log10|F|",
            "26D resonance sum R(t)", "UTe2 delta_n series",   "Anyon F_UBii",
            "H_res(Z,A,N)", "Ug2 harmonics + VacSeries"
        };
        for (int i = 0; i < NUM_SYSTEMS_002; ++i) {
            std::cout << std::left
                      << std::setw(4)  << i
                      << std::setw(42) << system_name(i)
                      << std::setw(20) << std::scientific << std::setprecision(3)
                      << compute_g_base(i)
                      << std::setw(24) << methods[i] << "\n";
        }

        std::cout << "\nSCALE RANGE\n";
        std::cout << "-----------\n";
        std::cout << "  r_min = r_any  = " << r_any
                  << " m  [sys 5: anyon, nuclear scale]\n";
        std::cout << "  r_max = r_vac  = " << r_vac
                  << " m  [sys 7: 1 Mpc cosmic void]\n";
        std::cout << "  log10(r_max/r_min) = " << std::fixed << std::setprecision(2)
                  << compute_scale_range() << " decades  (37.49 expected)\n";

        std::cout << "\nKEY PhD CONSTANTS (UQFF PhD Edition)\n";
        std::cout << "--------------------------------------\n";
        std::cout << "  F_rel        = " << std::scientific << F_rel
                  << " N   (2024 LEP calibration)\n";
        std::cout << "  E_LEP        = " << E_LEP << " J   (200 GeV)\n";
        std::cout << "  Q_wave       = " << Q_wave << "    (coupling amplifier)\n";
        std::cout << "  [SSq]        = " << SSq    << "    (vacuum shell index)\n";
        std::cout << "  gamma_decay  = " << gamma_d << " s^-1  (5e-5/day)\n";
        std::cout << "  delta_c_any  = " << delta_c_any << "    (Press-Schechter)\n";
        std::cout << "  B_thr(UTe2)  = " << 16.0   << " T\n";
        std::cout << "  Li26([SSq])  = " << compute_Vacuum_Series(100) << "  (Vacuum Series partial)\n";

        std::cout << "\nUTe2 delta_n SERIES (n=1..9)\n";
        std::cout << "------------------------------\n";
        for (int n = 1; n <= 9; ++n)
            std::cout << "  delta_" << n << " = " << std::scientific
                      << std::setprecision(4) << compute_delta_n_UTe2(n) << "\n";

        std::cout << "\nH_res HIGHLIGHTS\n";
        std::cout << "-----------------\n";
        std::cout << "  Fe-56  (Z=26, A=56, N=30): H_res = " << std::fixed
                  << std::setprecision(4) << compute_H_res(Z_fe56, A_fe56, N_fe56)
                  << "  [peak binding 8.79 MeV/nucleon]\n";
        std::cout << "  O-16   (Z=8,  A=16, N=8):  H_res = "
                  << compute_H_res(8,16,8) << "  [doubly-magic]\n";
        std::cout << "  Pb-208 (Z=82, A=208,N=126): H_res = "
                  << compute_H_res(82,208,126) << "  [doubly-magic, terminus]\n";
        std::cout << "  H-1    (Z=1,  A=1,  N=0):  H_res = "
                  << compute_H_res(1,1,0) << "  [lightest, reference A_H=1]\n";

        std::cout << "\nDPM INTEGRAL LOG-SCALE RESULTS\n";
        std::cout << "--------------------------------\n";
        for (int i = 0; i <= 2; ++i)
            std::cout << "  sys[" << i << "] " << std::left << std::setw(32)
                      << system_name(i) << ": log10|F_U_Bi_i| = "
                      << std::fixed << std::setprecision(2)
                      << compute_DPM_log_scale(i) << "\n";

        std::cout << "\nADVANCEMENT METRIC (PhD)\n";
        std::cout << "--------------------------\n";
        std::cout << "  diversity    = " << std::fixed << std::setprecision(3)
                  << diversity_score * 100.0 << " %  (8/8 regimes)\n";
        std::cout << "  dynamic      = " << dynamic_score * 100.0
                  << " %  (10/10 PhD processes)\n";
        std::cout << "  scalability  = " << scalability_score * 100.0
                  << " %  (" << compute_scale_range() << " / 40 decades)\n";
        std::cout << "  coverage     = " << coverage_score * 100.0
                  << " %  (8/8 systems)\n";
        std::cout << "  TOTAL        = " << compute_advancement() << " %\n\n";
    }

    /// Print sorted comparative analysis table at time t (seconds).
    void printComparativeAnalysis(double t = 0.0) const {
        std::cout << "\n";
        std::cout << "================================================================================\n";
        std::cout << "  Comparative UQFF Buoyancy Analysis  (t = "
                  << std::scientific << t << " s)\n";
        std::cout << "================================================================================\n";
        std::cout << std::left
                  << std::setw(4)  << "Idx"
                  << std::setw(34) << "System"
                  << std::setw(16) << "g_base (m/s^2)"
                  << std::setw(16) << "U_bi (m/s^2)"
                  << std::setw(18) << "F_UBii (N,scaled)"
                  << std::setw(18) << "U_m (full)\n";
        std::cout << std::string(106, '-') << "\n";

        // Collect sorted order by g_base
        std::array<int, NUM_SYSTEMS_002> order;
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [this](int a, int b) {
            return std::abs(compute_g_base(a)) > std::abs(compute_g_base(b));
        });

        for (int i : order) {
            double g    = compute_g_base(i);
            double ubi  = compute_Ubi(i);
            double fii  = compute_F_UBii(i, t);
            double um   = compute_Um_full(i, t);
            std::cout << std::left
                      << std::setw(4)  << i
                      << std::setw(34) << system_name(i)
                      << std::setw(16) << std::scientific << std::setprecision(3) << g
                      << std::setw(16) << ubi
                      << std::setw(18) << fii
                      << std::setw(18) << um << "\n";
        }
        std::cout << "\n";
        std::cout << "  NOTE: F_UBii shown is Tier-2 standard (DPM-scaled for sys[0-2]).\n";
        std::cout << "  For calibrated DPM integral use compute_DPM_log_scale(idx).\n";
        std::cout << "  For anyon: compute_F_UBii_anyons(); 26D: compute_26D_resonance_sum().\n\n";
    }

    /// Print compact parameter dump for all 8 systems.
    void printParameters() const {
        std::cout << "\n--- UQFFLearningAssessment002 Parameters ---\n";
        std::cout << std::scientific << std::setprecision(4);
        auto row = [&](const char* tag, double v) {
            std::cout << "  " << std::left << std::setw(22) << tag << " = " << v << "\n";
        };
        std::cout << "[0] Westerlund 2 Triadic\n";
        row("M_w2t",       M_w2t);       row("r_w2t",     r_w2t);
        row("DPM_res_w2",  DPM_res_w2);  row("F_LENR_w2", F_LENR_w2);
        row("T_sf_w2",     T_sf_w2);
        std::cout << "[1] Sgr A* Buoyancy\n";
        row("M_sgra_b",     M_sgra_b);   row("r_sgra_b",   r_sgra_b);
        row("DPM_res_sgra", DPM_res_sgra); row("F_LENR_sgra", F_LENR_sgra);
        std::cout << "[2] Pillars DPM\n";
        row("M_pil_d",     M_pil_d);     row("r_pil_d",   r_pil_d);
        row("DPM_res_pil", DPM_res_pil); row("F_LENR_pil", F_LENR_pil);
        std::cout << "[3] 26D Resonance\n";
        row("M_res26",     M_res26);     row("r_res26",    r_res26);
        row("N_layers",    N_layers_res26); row("T_sf_res26", T_sf_res26);
        row("SSq_res26",   SSq_res26);
        std::cout << "[4] UTe2 SC\n";
        row("M_ute2",      M_ute2);      row("r_ute2",    r_ute2);
        row("B_ute2",      B_ute2);      row("f_topo_ute2", f_topo_ute2);
        std::cout << "[5] Anyons CERN\n";
        row("M_any",       M_any);       row("r_any",     r_any);
        row("E_anyons_loc", E_anyons_loc); row("delta_c_loc", delta_c_loc);
        std::cout << "[6] H_res Fe-56\n";
        row("M_fe56",      M_fe56);      row("r_fe56",    r_fe56);
        row("Z_fe56",      Z_fe56);      row("A_fe56",    A_fe56);
        row("N_fe56",      N_fe56);
        std::cout << "[7] Vacuum Harmonics\n";
        row("M_vac",       M_vac);       row("r_vac",     r_vac);
        row("N_harmonics", N_harmonics_vac); row("gamma_vac", gamma_vac_loc);
        std::cout << "---\n";
    }

    /// Run a representative PhD-level research state example.
    /// Evaluates all 8 systems at a 1-Myr reference time.
    void exampleAtResearchState() const {
        constexpr double t_1Myr = 3.156e13;   // s  (1 million years)
        printPhDGuide();
        printComparativeAnalysis(t_1Myr);

        std::cout << "\nPHD SPECIFIC RESULTS AT t = 1 Myr\n";
        std::cout << "-----------------------------------\n";

        // 26D resonance sum
        double R26 = compute_26D_resonance_sum(SYS_RESONANCE_26D, t_1Myr);
        std::cout << "  R(t=1Myr, 26D)         = " << std::scientific
                  << std::setprecision(4) << R26 << " N\n";

        // Anyon F_UBii
        double F_any = compute_F_UBii_anyons();
        std::cout << "  F_UBii,anyons (g_any)  = " << F_any << " N\n";

        // Vacuum series
        double V_ser = compute_Vacuum_Series(26);
        std::cout << "  Vac.Series(N=26)       = " << V_ser
                  << "  [Li26(0.57) approx]\n";

        // U_g2 harmonics at t=0 (cosmic reference)
        double Ug2  = compute_Ug2_harmonics(0.0, 26);
        std::cout << "  U_g2 harmonics(t=0)    = " << Ug2 << " [n.d.]\n";

        // Dynamic [SSq] at n=13, t=0
        double ssq_dyn = compute_SSq_dynamic(13, t_1Myr);
        std::cout << "  SSq_dynamic(n=13,1Myr) = " << ssq_dyn << "\n";

        // Full Um for Sgr A* at 1 Myr
        double Um_sgra = compute_Um_full(SYS_SGRA_BUOY, t_1Myr);
        std::cout << "  Um_full(SgrA*, 1Myr)   = " << Um_sgra << " (arb.)\n";

        std::cout << "\nDPM LOG-SCALE SUMMARY\n";
        std::cout << "---------------------\n";
        const char* sys_names[3] = {"W2 Triadic", "SgrA* Buoy", "Pillars DPM"};
        for (int i = 0; i < 3; ++i)
            std::cout << "  log10|F_U_Bi_i| (" << sys_names[i] << ") = "
                      << std::fixed << std::setprecision(1)
                      << compute_DPM_log_scale(i) << "\n";

        std::cout << "\n  PhD Advancement Score = "
                  << std::fixed << std::setprecision(2) << compute_advancement()
                  << " %\n\n";
    }

    // ─────────────────────────────────────────────────────────────────────────
    // §D  UQFF 2.0 Self-Expanding Framework — Variable Access
    // ─────────────────────────────────────────────────────────────────────────

    bool setVariable(const std::string& varName, double newValue) {
        // ── [0] W2 Triadic ────────────────────────────────────────────────────
        if      (varName == "M_w2t")        M_w2t        = newValue;
        else if (varName == "r_w2t")        r_w2t        = newValue;
        else if (varName == "DPM_res_w2")   DPM_res_w2   = newValue;
        else if (varName == "F_LENR_w2")    F_LENR_w2    = newValue;
        else if (varName == "T_sf_w2")      T_sf_w2      = newValue;
        else if (varName == "M_ext_w2t")    M_ext_w2t    = newValue;
        else if (varName == "r_ext_w2t")    r_ext_w2t    = newValue;
        // ── [1] SgrA* Buoyancy ────────────────────────────────────────────────
        else if (varName == "M_sgra_b")     M_sgra_b     = newValue;
        else if (varName == "r_sgra_b")     r_sgra_b     = newValue;
        else if (varName == "DPM_res_sgra") DPM_res_sgra = newValue;
        else if (varName == "F_LENR_sgra")  F_LENR_sgra  = newValue;
        else if (varName == "T_sf_sgra")    T_sf_sgra    = newValue;
        else if (varName == "M_ext_sgra_b") M_ext_sgra_b = newValue;
        else if (varName == "r_ext_sgra_b") r_ext_sgra_b = newValue;
        // ── [2] Pillars DPM ───────────────────────────────────────────────────
        else if (varName == "M_pil_d")      M_pil_d      = newValue;
        else if (varName == "r_pil_d")      r_pil_d      = newValue;
        else if (varName == "DPM_res_pil")  DPM_res_pil  = newValue;
        else if (varName == "F_LENR_pil")   F_LENR_pil   = newValue;
        else if (varName == "T_sf_pil")     T_sf_pil     = newValue;
        else if (varName == "M_ext_pil_d")  M_ext_pil_d  = newValue;
        else if (varName == "r_ext_pil_d")  r_ext_pil_d  = newValue;
        // ── [3] 26D Resonance ─────────────────────────────────────────────────
        else if (varName == "M_res26")         M_res26         = newValue;
        else if (varName == "r_res26")         r_res26         = newValue;
        else if (varName == "N_layers_res26")  N_layers_res26  = newValue;
        else if (varName == "T_sf_res26")      T_sf_res26      = newValue;
        else if (varName == "SSq_res26")       SSq_res26       = newValue;
        else if (varName == "M_ext_res26")     M_ext_res26     = newValue;
        else if (varName == "r_ext_res26")     r_ext_res26     = newValue;
        // ── [4] UTe2 SC ───────────────────────────────────────────────────────
        else if (varName == "M_ute2")       M_ute2       = newValue;
        else if (varName == "r_ute2")       r_ute2       = newValue;
        else if (varName == "B_ute2")       B_ute2       = newValue;
        else if (varName == "f_topo_ute2")  f_topo_ute2  = newValue;
        else if (varName == "M_ext_ute2")   M_ext_ute2   = newValue;
        else if (varName == "r_ext_ute2")   r_ext_ute2   = newValue;
        // ── [5] Anyons ────────────────────────────────────────────────────────
        else if (varName == "M_any")         M_any         = newValue;
        else if (varName == "r_any")         r_any         = newValue;
        else if (varName == "E_anyons_loc")  E_anyons_loc  = newValue;
        else if (varName == "delta_c_loc")   delta_c_loc   = newValue;
        else if (varName == "sigma_loc")     sigma_loc     = newValue;
        else if (varName == "M_ext_any")     M_ext_any     = newValue;
        else if (varName == "r_ext_any")     r_ext_any     = newValue;
        // ── [6] H_res Fe-56 ───────────────────────────────────────────────────
        else if (varName == "M_fe56")     M_fe56     = newValue;
        else if (varName == "r_fe56")     r_fe56     = newValue;
        else if (varName == "Z_fe56")     Z_fe56     = newValue;
        else if (varName == "A_fe56")     A_fe56     = newValue;
        else if (varName == "N_fe56")     N_fe56     = newValue;
        else if (varName == "M_ext_fe56") M_ext_fe56 = newValue;
        else if (varName == "r_ext_fe56") r_ext_fe56 = newValue;
        // ── [7] Vacuum Harmonics ──────────────────────────────────────────────
        else if (varName == "M_vac")           M_vac           = newValue;
        else if (varName == "r_vac")           r_vac           = newValue;
        else if (varName == "N_harmonics_vac") N_harmonics_vac = newValue;
        else if (varName == "gamma_vac_loc")   gamma_vac_loc   = newValue;
        else if (varName == "M_ext_vac")       M_ext_vac       = newValue;
        else if (varName == "r_ext_vac")       r_ext_vac       = newValue;
        // ── Assessment metrics ────────────────────────────────────────────────
        else if (varName == "diversity_score")   diversity_score   = newValue;
        else if (varName == "dynamic_score")     dynamic_score     = newValue;
        else if (varName == "scalability_score") scalability_score = newValue;
        else if (varName == "coverage_score")    coverage_score    = newValue;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'.\n";
            return false;
        }
        if (logging_enabled)
            std::cout << "[LOG] setVariable: " << varName << " = " << newValue << "\n";
        return true;
    }

    bool addToVariable(const std::string& varName, double delta) {
        return setVariable(varName, getVariable(varName) + delta);
    }

    bool subtractFromVariable(const std::string& varName, double delta) {
        return addToVariable(varName, -delta);
    }

    double getVariable(const std::string& varName) const {
        if      (varName == "M_w2t")        return M_w2t;
        else if (varName == "r_w2t")        return r_w2t;
        else if (varName == "DPM_res_w2")   return DPM_res_w2;
        else if (varName == "F_LENR_w2")    return F_LENR_w2;
        else if (varName == "T_sf_w2")      return T_sf_w2;
        else if (varName == "M_ext_w2t")    return M_ext_w2t;
        else if (varName == "r_ext_w2t")    return r_ext_w2t;
        else if (varName == "M_sgra_b")     return M_sgra_b;
        else if (varName == "r_sgra_b")     return r_sgra_b;
        else if (varName == "DPM_res_sgra") return DPM_res_sgra;
        else if (varName == "F_LENR_sgra")  return F_LENR_sgra;
        else if (varName == "T_sf_sgra")    return T_sf_sgra;
        else if (varName == "M_ext_sgra_b") return M_ext_sgra_b;
        else if (varName == "r_ext_sgra_b") return r_ext_sgra_b;
        else if (varName == "M_pil_d")      return M_pil_d;
        else if (varName == "r_pil_d")      return r_pil_d;
        else if (varName == "DPM_res_pil")  return DPM_res_pil;
        else if (varName == "F_LENR_pil")   return F_LENR_pil;
        else if (varName == "T_sf_pil")     return T_sf_pil;
        else if (varName == "M_ext_pil_d")  return M_ext_pil_d;
        else if (varName == "r_ext_pil_d")  return r_ext_pil_d;
        else if (varName == "M_res26")         return M_res26;
        else if (varName == "r_res26")         return r_res26;
        else if (varName == "N_layers_res26")  return N_layers_res26;
        else if (varName == "T_sf_res26")      return T_sf_res26;
        else if (varName == "SSq_res26")       return SSq_res26;
        else if (varName == "M_ext_res26")     return M_ext_res26;
        else if (varName == "r_ext_res26")     return r_ext_res26;
        else if (varName == "M_ute2")       return M_ute2;
        else if (varName == "r_ute2")       return r_ute2;
        else if (varName == "B_ute2")       return B_ute2;
        else if (varName == "f_topo_ute2")  return f_topo_ute2;
        else if (varName == "M_ext_ute2")   return M_ext_ute2;
        else if (varName == "r_ext_ute2")   return r_ext_ute2;
        else if (varName == "M_any")         return M_any;
        else if (varName == "r_any")         return r_any;
        else if (varName == "E_anyons_loc")  return E_anyons_loc;
        else if (varName == "delta_c_loc")   return delta_c_loc;
        else if (varName == "sigma_loc")     return sigma_loc;
        else if (varName == "M_ext_any")     return M_ext_any;
        else if (varName == "r_ext_any")     return r_ext_any;
        else if (varName == "M_fe56")     return M_fe56;
        else if (varName == "r_fe56")     return r_fe56;
        else if (varName == "Z_fe56")     return Z_fe56;
        else if (varName == "A_fe56")     return A_fe56;
        else if (varName == "N_fe56")     return N_fe56;
        else if (varName == "M_ext_fe56") return M_ext_fe56;
        else if (varName == "r_ext_fe56") return r_ext_fe56;
        else if (varName == "M_vac")           return M_vac;
        else if (varName == "r_vac")           return r_vac;
        else if (varName == "N_harmonics_vac") return N_harmonics_vac;
        else if (varName == "gamma_vac_loc")   return gamma_vac_loc;
        else if (varName == "M_ext_vac")       return M_ext_vac;
        else if (varName == "r_ext_vac")       return r_ext_vac;
        else if (varName == "diversity_score")   return diversity_score;
        else if (varName == "dynamic_score")     return dynamic_score;
        else if (varName == "scalability_score") return scalability_score;
        else if (varName == "coverage_score")    return coverage_score;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'.\n";
            return 0.0;
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // §E  UQFF 2.0 Self-Expanding Framework — Dynamic Params & State
    // ─────────────────────────────────────────────────────────────────────────

    void setEnableLogging(bool enabled) { logging_enabled = enabled; }
    bool getLoggingEnabled() const      { return logging_enabled; }

    void setDynamicParameter(const std::string& key, double value) {
        dynamic_params[key] = value;
        if (logging_enabled)
            std::cout << "[LOG] Dynamic param set: " << key << " = " << value << "\n";
    }

    double getDynamicParameter(const std::string& key) const {
        auto it = dynamic_params.find(key);
        if (it != dynamic_params.end()) return it->second;
        if (logging_enabled)
            std::cout << "[LOG] Dynamic param not found: " << key << "\n";
        return std::numeric_limits<double>::quiet_NaN();
    }

    void exportState(const std::string& filename = "UQFFLearningAssessment002_state.txt") const {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "Error: Cannot open '" << filename << "' for state export.\n";
            return;
        }
        ofs << "# UQFFLearningAssessment002 (PhD Research Edition) state export\n";
        ofs << std::scientific << std::setprecision(10);
        const char* names[] = {
            // [0] W2
            "M_w2t", "r_w2t", "DPM_res_w2", "F_LENR_w2", "T_sf_w2",
            "M_ext_w2t", "r_ext_w2t",
            // [1] SgrA*
            "M_sgra_b", "r_sgra_b", "DPM_res_sgra", "F_LENR_sgra", "T_sf_sgra",
            "M_ext_sgra_b", "r_ext_sgra_b",
            // [2] Pillars
            "M_pil_d", "r_pil_d", "DPM_res_pil", "F_LENR_pil", "T_sf_pil",
            "M_ext_pil_d", "r_ext_pil_d",
            // [3] 26D
            "M_res26", "r_res26", "N_layers_res26", "T_sf_res26", "SSq_res26",
            "M_ext_res26", "r_ext_res26",
            // [4] UTe2
            "M_ute2", "r_ute2", "B_ute2", "f_topo_ute2",
            "M_ext_ute2", "r_ext_ute2",
            // [5] Anyons
            "M_any", "r_any", "E_anyons_loc", "delta_c_loc", "sigma_loc",
            "M_ext_any", "r_ext_any",
            // [6] Fe-56
            "M_fe56", "r_fe56", "Z_fe56", "A_fe56", "N_fe56",
            "M_ext_fe56", "r_ext_fe56",
            // [7] Vac
            "M_vac", "r_vac", "N_harmonics_vac", "gamma_vac_loc",
            "M_ext_vac", "r_ext_vac",
            // metrics
            "diversity_score", "dynamic_score", "scalability_score", "coverage_score"
        };
        for (const char* n : names)
            ofs << n << " = " << getVariable(n) << "\n";
        for (const auto& kv : dynamic_params)
            ofs << "dynamic." << kv.first << " = " << kv.second << "\n";
        if (logging_enabled)
            std::cout << "[LOG] State exported to: " << filename << "\n";
    }

    /// Cross-validate g_base[sys_idx] against an external UQFF module.
    /// Store the other module's g beforehand: setDynamicParameter("cross_g", val).
    template<typename OtherModuleT>
    double cross_validate(const OtherModuleT& /*other*/, int sys_idx,
                          double /*t_seconds*/ = 0.0) const {
        double g_here  = compute_g_base(sys_idx);
        double g_other = 0.0;
        auto it = dynamic_params.find("cross_g");
        if (it != dynamic_params.end()) g_other = it->second;
        double frac = (g_here != 0.0)
                    ? std::abs(g_other - g_here) / std::abs(g_here)
                    : std::numeric_limits<double>::quiet_NaN();
        if (logging_enabled)
            std::cout << "[XVAL002] sys[" << sys_idx << "] g_here=" << g_here
                      << "  g_other=" << g_other
                      << "  frac_diff=" << frac << "\n";
        return frac;
    }
};

#endif // UQFF_LEARNING_ASSESSMENT_002_H
