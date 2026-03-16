/**
 * ================================================================================================
 * Header: UQFFLearningAssessment.h
 *
 * Description: C++ Module for UQFF Learning Assessment — Research Student Edition (v2.0)
 *              This is the ninth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations. Covers ALL SEVEN completed astrophysical systems
 *              in a unified research-student comparative framework:
 *                [0] SGR 1745-2900      — Magnetar (compact, GC proximity, D(t) burst, 18 terms)
 *                [1] SGR 0501+4516      — Magnetar (compact, Perseus arm, D(t) burst, 16 terms)
 *                [2] Sgr A*             — SMBH (galactic center, M(t) accretion, QPO, 15 terms)
 *                [3] Starbirth Tapestry — Young cluster (LMC NGC2014/2020, M(t) + wind, 12 terms)
 *                [4] Westerlund 2       — Super cluster (Carina, M(t) + OB winds, 12 terms)
 *                [5] Pillars of Creation— Star-forming pillars (Eagle Nebula M16, E(t) erosion, 12 terms)
 *                [6] Rings of Relativity— Einstein ring (GAL-CLUS-022058s, H(z), L_t lensing, 12 terms)
 *
 * Purpose: Research-student comparative framework for UQFF mastery assessment.
 *          - Aggregates canonical parameters from all 7 completed C++ modules (Sessions 63-70)
 *          - Computes g_base (Newtonian G*M/r^2), Ubi (buoyancy Tier-1), and advancement metrics
 *          - Metrics: diversity (7/7 regimes), dynamic processes (7/10 unique), scale span
 *            (~16.2 log-decades in r), coverage (7/7 planned modules complete)
 *          - Advancement score: weighted mean of 4 metrics -> 0-100%
 *          - Educational output: printStudentGuide() (unique physics per system, dynamic catalogue,
 *            scale analysis), printComparativeAnalysis() (sorted g-table across all 7 systems)
 *          - Full UQFF 2.0 Self-Expanding Framework (logging, dynamic_params, exportState,
 *            cross_validate<>)
 *
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: UQFFLearningAssessment assess;
 *              double adv = assess.compute_advancement();
 *              assess.printStudentGuide();
 *              assess.printComparativeAnalysis(1.0e6 * 3.156e7);  // at t = 1 Myr
 *
 * Key Features:
 *   - Canonical parameters sourced directly from all 7 completed UQFF C++ modules
 *   - 7 physical regimes: magnetar x2, SMBH, LMC cluster, MW super-cluster,
 *     MW star-forming pillars, cosmological Einstein ring
 *   - Scale range: r = 20 km (magnetar) to 10 kpc (Einstein ring) = 16.2 log-decades
 *   - Mass range: M = 1.4 M_sun (magnetar) to 1e14 M_sun (galaxy cluster) = ~20 dex
 *   - compute_advancement(): 4-metric weighted UQFF score (diversity, dynamics, scale, coverage)
 *   - compute_g_base(int): Newtonian g = G*M/r^2 for system 0 through 6
 *   - compute_Ubi(int): UQFF buoyancy Tier-1 = 0.5 * g_base
 *   - compute_F_UBii(int, double): UQFF buoyancy Tier-2 at time t
 *   - compute_Ub_i(int, double): UQFF buoyancy Tier-3 (external body) at time t
 *   - printStudentGuide(): educational walkthrough + unique physics catalogue
 *   - printComparativeAnalysis(double t): sorted g-table for all 7 systems
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: March 16, 2026 (Research Student Edition — replaced October 08, 2025 stub)
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef UQFF_LEARNING_ASSESSMENT_H
#define UQFF_LEARNING_ASSESSMENT_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <limits>
#include <map>
#include <fstream>
#include <array>
#include <algorithm>

static constexpr double UQFF_LA_PI = 3.14159265358979323846;

// System index constants
static constexpr int SYS_SGR1745    = 0;
static constexpr int SYS_SGR0501    = 1;
static constexpr int SYS_SGR_A_STAR = 2;
static constexpr int SYS_TAPESTRY   = 3;
static constexpr int SYS_WESTERLUND = 4;
static constexpr int SYS_PILLARS    = 5;
static constexpr int SYS_RINGS      = 6;
static constexpr int NUM_SYSTEMS    = 7;

class UQFFLearningAssessment {
private:
    // ── Shared physical constants ─────────────────────────────────────────────
    double G;               // Gravitational constant (m^3 kg^-1 s^-2)
    double c_light;         // Speed of light (m/s)
    double hbar;            // Reduced Planck constant (J·s)
    double Lambda;          // Cosmological constant (m^-2)
    double M_sun;           // Solar mass (kg)
    double f_TRZ;           // UQFF time-reversal factor (shared, = 0.1)
    double rho_vac_UA;      // UA vacuum density (J/m^3)
    double rho_vac_SCm;     // SCm vacuum density (J/m^3)
    double t_Hubble;        // Hubble time (s)

    // ── UQFF buoyancy pipeline (CP3/PAPER_198 canonical) ─────────────────────
    double beta_i;          // Buoyancy coupling constant (= 0.61)
    double omega_g;         // Galactic rotation rate (rad/s, = 7.3e-16)
    double U_UA;            // Unit charge aether parameter (C, = 1e-11)

    // ── Assessment metrics ────────────────────────────────────────────────────
    double diversity_score;   // Fraction of distinct physical regimes covered (0-1)
    double dynamic_score;     // Fraction of distinct dynamic-physics processes covered (0-1)
    double scalability_score; // Scale-span coverage (0-1, based on log10 r range / 20)
    double coverage_score;    // Module completion fraction (completed / planned, 0-1)

    // ── System 0: SGR 1745-2900 (Magnetar — GC proximity, 0.92 pc from Sgr A*) ──────────
    double M_sgr1745;           // 1.4 M_sun (kg)
    double r_sgr1745;           // 20 km neutron-star radius (m)
    double B_sgr1745;           // Magnetic field ~1e14 T (near quantum-critical)
    double D0_sgr1745;          // D(t) burst amplitude (m/s^2)
    double omega_D_sgr1745;     // Burst angular freq: 2*pi/11 s (rad/s)
    double tau_D_sgr1745;       // Burst decay timescale: 3.5 yr (s)
    double M_ext_sgr1745;       // External body: Sgr A* ~4e6 M_sun (kg)
    double r_ext_sgr1745;       // Sgr A* distance: 0.92 pc = 2.83e16 m

    // ── System 1: SGR 0501+4516 (Magnetar — Perseus arm, 10 kpc from GC) ────────────────
    double M_sgr0501;           // 1.4 M_sun (kg)
    double r_sgr0501;           // 20 km neutron-star radius (m)
    double B_sgr0501;           // Magnetic field ~1.9e14 T (SGR0501 canonical)
    double D0_sgr0501;          // D(t) burst amplitude (m/s^2)
    double omega_D_sgr0501;     // Burst angular freq: 2*pi/15 s (2008 outburst cadence, rad/s)
    double tau_D_sgr0501;       // Burst decay timescale: 3.5 yr (s)
    double M_ext_sgr0501;       // External body: Sgr A* ~4e6 M_sun (kg)
    double r_ext_sgr0501;       // Distance to GC: 10 kpc = 3.086e20 m

    // ── System 2: Sgr A* (SMBH — galactic center) ────────────────────────────────────────
    double M_smbh;              // Initial: 4e6 M_sun — M(t) grows by accretion (kg)
    double r_smbh;              // Schwarzschild-region reference radius: 1.27e10 m
    double M_dot_smbh;          // Accretion growth factor (dimensionless)
    double tau_acc_smbh;        // Accretion timescale: 1 Myr (s)
    double omega_D_smbh;        // QPO angular freq: 2*pi/1200 s (20-min Sgr A* QPO, rad/s)
    double tau_D_smbh;          // QPO envelope decay: 3 yr (s)
    double spin_smbh;           // Kerr spin factor (dimensionless, = 0.3)
    double M_ext_smbh;          // External body: NSC ~2.5e7 M_sun (kg)
    double r_ext_smbh;          // NSC distance: 3 pc = 9.258e16 m

    // ── System 3: Starbirth Tapestry (LMC — NGC 2014/2020 young cluster) ─────────────────
    double M_tapestry;          // Initial: 240 M_sun (kg)
    double r_tapestry;          // 10 ly = 9.461e16 m
    double M_dot_tapestry;      // SF mass-growth factor (~41.7, gas/stellar ratio)
    double tau_SF_tapestry;     // Star-formation timescale: 5 Myr (s)
    double rho_wind_tap;        // Wind density: 1e-21 kg/m^3
    double v_wind_tap;          // Wind velocity: 2e6 m/s
    double M_ext_tap;           // External body: R136 Tarantula core ~1e5 M_sun (kg)
    double r_ext_tap;           // R136 distance: 650 ly = 6.15e18 m

    // ── System 4: Westerlund 2 (super star cluster — Carina arm, ~4.16 kpc) ──────────────
    double M_wd2;               // Initial: 30,000 M_sun (kg)
    double r_wd2;               // 10 ly = 9.461e16 m
    double M_dot_wd2;           // SF mass-growth factor (~3.33, gas/stellar ratio)
    double tau_SF_wd2;          // Star-formation timescale: 2 Myr (s)
    double rho_wind_wd2;        // Dense OB-star wind: 1e-20 kg/m^3
    double v_wind_wd2;          // Wind velocity: 2e6 m/s
    double M_ext_wd2;           // External body: Sgr A* ~4e6 M_sun (kg)
    double r_ext_wd2;           // Wd2 to GC: 8 kpc = 2.47e20 m

    // ── System 5: Pillars of Creation (Eagle Nebula M16 — Sagittarius arm) ───────────────
    double M_pillars;           // Initial: 10,100 M_sun (kg)
    double r_pillars;           // 5 ly = 4.731e16 m
    double M_dot_pillars;       // SF mass-growth factor (~0.99)
    double tau_SF_pillars;      // Star-formation timescale: 1 Myr (s)
    double E_0_pillars;         // Erosion amplitude: 0.1
    double tau_E_pillars;       // Erosion timescale: 1 Myr (s)
    double rho_wind_pil;        // Wind density: 1e-21 kg/m^3
    double v_wind_pil;          // Wind velocity: 2e6 m/s
    double M_ext_pil;           // External body: Sgr A* ~4e6 M_sun (kg)
    double r_ext_pil;           // Pillars to GC: 6.85 kpc = 2.11e20 m

    // ── System 6: Rings of Relativity (GAL-CLUS-022058s — z=0.5 Einstein ring) ──────────
    double M_rings;             // Galaxy cluster: 1e14 M_sun; static (kg)
    double r_rings;             // Einstein radius: 10 kpc = 3.086e20 m
    double Hz_rings;            // H(z=0.5) = 2.42e-18 s^-1 (redshift-specific Hubble)
    double L_factor_rings;      // D_LS/D_S = 0.67 (lensing geometry)
    double z_rings;             // Lens redshift: 0.5
    double rho_wind_rings;      // Intracluster wind density: 1e-21 kg/m^3
    double v_wind_rings;        // ICM wind velocity: 2e6 m/s
    double M_ext_rings;         // External body: Virgo Cluster ~1.2e15 M_sun (kg)
    double r_ext_rings;         // Distance to Virgo: 16.5 Mpc = 5.09e23 m

    // ── UQFF 2.0 Self-Expanding Framework ────────────────────────────────────
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

    // ── Internal helpers (private) ────────────────────────────────────────────
    double g_base_for(int idx) const {
        switch (idx) {
            case SYS_SGR1745:    return (G * M_sgr1745)  / (r_sgr1745  * r_sgr1745);
            case SYS_SGR0501:    return (G * M_sgr0501)  / (r_sgr0501  * r_sgr0501);
            case SYS_SGR_A_STAR: return (G * M_smbh)     / (r_smbh     * r_smbh);
            case SYS_TAPESTRY:   return (G * M_tapestry) / (r_tapestry * r_tapestry);
            case SYS_WESTERLUND: return (G * M_wd2)      / (r_wd2      * r_wd2);
            case SYS_PILLARS:    return (G * M_pillars)  / (r_pillars  * r_pillars);
            case SYS_RINGS:      return (G * M_rings)    / (r_rings    * r_rings);
            default:             return 0.0;
        }
    }

    double get_M(int idx) const {
        switch (idx) {
            case SYS_SGR1745:    return M_sgr1745;
            case SYS_SGR0501:    return M_sgr0501;
            case SYS_SGR_A_STAR: return M_smbh;
            case SYS_TAPESTRY:   return M_tapestry;
            case SYS_WESTERLUND: return M_wd2;
            case SYS_PILLARS:    return M_pillars;
            case SYS_RINGS:      return M_rings;
            default:             return 0.0;
        }
    }

    double get_r(int idx) const {
        switch (idx) {
            case SYS_SGR1745:    return r_sgr1745;
            case SYS_SGR0501:    return r_sgr0501;
            case SYS_SGR_A_STAR: return r_smbh;
            case SYS_TAPESTRY:   return r_tapestry;
            case SYS_WESTERLUND: return r_wd2;
            case SYS_PILLARS:    return r_pillars;
            case SYS_RINGS:      return r_rings;
            default:             return 0.0;
        }
    }

    static const char* system_name(int idx) {
        static const char* names[NUM_SYSTEMS] = {
            "SGR 1745-2900 (Magnetar)",
            "SGR 0501+4516 (Magnetar)",
            "Sgr A* (SMBH)",
            "Starbirth Tapestry (LMC cluster)",
            "Westerlund 2 (super cluster)",
            "Pillars of Creation (Eagle Nebula M16)",
            "Rings of Relativity (Einstein ring)"
        };
        return (idx >= 0 && idx < NUM_SYSTEMS) ? names[idx] : "Unknown";
    }

    static const char* unique_physics_str(int idx) {
        static const char* phys[NUM_SYSTEMS] = {
            "D(t) burst [omega_D=2pi/11s, tau_D=3.5yr]; BH tidal [2G*M_BH*r/r_BH^3]; 18 terms; Sgr A* ext 0.92 pc",
            "D(t) burst [omega_D=2pi/15s, 2008 cadence]; GW back-reaction; galactic tidal; 16 terms; Sgr A* ext 10 kpc",
            "M(t) accretion growth; QPO burst [omega_D=2pi/1200s, 20-min]; NSC tidal; magnetic energy; 15 terms",
            "M(t) SF growth [tau=5Myr]; stellar wind [rho=1e-21]; R136 Tarantula ext 650ly; 12 terms; LMC/cosmological isolation",
            "M(t) SF growth [tau=2Myr]; OB-star wind [rho=1e-20, richest]; Sgr A* ext 8 kpc; 12 terms; Carina arm",
            "M(t) SF growth + E(t) erosion [corr_E=1-E0*exp(-t/tau)]; 4 simultaneous co-actions; Sgr A* ext 6.85 kpc; 12 terms",
            "Static M; H(z=0.5)=2.42e-18 s^-1 [H(z) — 1st cosmo module]; L_t=(GM/c^2*r)*Lfac -> corr_L=1+L_t; Virgo ext 16.5 Mpc; 12 terms"
        };
        return (idx >= 0 && idx < NUM_SYSTEMS) ? phys[idx] : "Unknown";
    }

    static const char* regime_label(int idx) {
        static const char* labels[NUM_SYSTEMS] = {
            "Compact / Nuclear-quantum",
            "Compact / Nuclear-quantum",
            "SMBH / Galactic center",
            "Young stellar cluster / LMC",
            "Super star cluster / MW Carina arm",
            "Star-forming pillars / MW Sagittarius arm",
            "Cosmological / Einstein ring"
        };
        return (idx >= 0 && idx < NUM_SYSTEMS) ? labels[idx] : "Unknown";
    }

    // External body: M_ext and r_ext for each system
    void get_ext(int idx, double& M_ext, double& r_ext) const {
        switch (idx) {
            case SYS_SGR1745:    M_ext = M_ext_sgr1745; r_ext = r_ext_sgr1745; break;
            case SYS_SGR0501:    M_ext = M_ext_sgr0501; r_ext = r_ext_sgr0501; break;
            case SYS_SGR_A_STAR: M_ext = M_ext_smbh;    r_ext = r_ext_smbh;    break;
            case SYS_TAPESTRY:   M_ext = M_ext_tap;     r_ext = r_ext_tap;     break;
            case SYS_WESTERLUND: M_ext = M_ext_wd2;     r_ext = r_ext_wd2;     break;
            case SYS_PILLARS:    M_ext = M_ext_pil;     r_ext = r_ext_pil;     break;
            case SYS_RINGS:      M_ext = M_ext_rings;   r_ext = r_ext_rings;   break;
            default:             M_ext = 0.0;           r_ext = 1.0;           break;
        }
    }

public:
    // ── Constructor ───────────────────────────────────────────────────────────
    UQFFLearningAssessment() { initializeDefaults(); }
    ~UQFFLearningAssessment() {}

    // ── Initialization ────────────────────────────────────────────────────────
    void initializeDefaults() {
        // Shared constants
        G           = 6.6743e-11;
        c_light     = 3.0e8;
        hbar        = 1.0546e-34;
        Lambda      = 1.1e-52;
        M_sun       = 1.989e30;
        f_TRZ       = 0.1;
        rho_vac_UA  = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        t_Hubble    = 13.8e9 * 3.15576e7;

        // UQFF buoyancy pipeline (CP3/PAPER_198 canonical)
        beta_i  = 0.61;
        omega_g = 7.3e-16;
        U_UA    = 1e-11;

        // Assessment metrics
        // diversity:   7 distinct physical regimes covered out of 7
        // dynamic:     7 unique time-dependent physics processes introduced
        // scalability: log10(r_max/r_min) = log10(3.086e20 / 2e4) ~= 16.2 / 20 = 0.81
        // coverage:    7/7 planned C++ modules complete
        diversity_score   = 7.0 / 7.0;
        dynamic_score     = 7.0 / 10.0;
        scalability_score = 16.2 / 20.0;
        coverage_score    = 7.0 / 7.0;

        // ── System 0: SGR 1745-2900 (Magnetar — 0.92 pc from Sgr A*) ──────────────────
        M_sgr1745       = 1.4 * M_sun;                     // 1.4 M_sun
        r_sgr1745       = 20.0e3;                           // 20 km
        B_sgr1745       = 1.0e14;                           // Near quantum-critical field (T)
        D0_sgr1745      = 1.0e10;                           // Burst amplitude (m/s^2)
        omega_D_sgr1745 = 2.0 * UQFF_LA_PI / 11.0;        // 11 s burst period (rad/s)
        tau_D_sgr1745   = 3.5 * 3.15576e7;                 // 3.5 yr burst decay (s)
        M_ext_sgr1745   = 4.0e6 * M_sun;                   // Sgr A* (kg)
        r_ext_sgr1745   = 2.83e16;                          // 0.92 pc (m)

        // ── System 1: SGR 0501+4516 (Magnetar — Perseus arm) ───────────────────────────
        M_sgr0501       = 1.4 * M_sun;                     // 1.4 M_sun
        r_sgr0501       = 20.0e3;                           // 20 km
        B_sgr0501       = 1.9e14;                           // SGR0501 canonical field (T)
        D0_sgr0501      = 1.0e10;                           // Burst amplitude (m/s^2)
        omega_D_sgr0501 = 2.0 * UQFF_LA_PI / 15.0;        // 15 s burst cadence (2008 outburst, rad/s)
        tau_D_sgr0501   = 3.5 * 3.15576e7;                 // 3.5 yr burst decay (s)
        M_ext_sgr0501   = 4.0e6 * M_sun;                   // Sgr A* (kg)
        r_ext_sgr0501   = 3.086e20;                         // 10 kpc (m)

        // ── System 2: Sgr A* (SMBH — galactic center) ──────────────────────────────────
        M_smbh          = 4.0e6 * M_sun;                   // 4e6 M_sun initial
        r_smbh          = 1.27e10;                          // Schwarzschild-region radius (m)
        M_dot_smbh      = 1.0e-4;                           // Accretion growth factor
        tau_acc_smbh    = 1.0e6 * 3.15576e7;               // 1 Myr accretion timescale (s)
        omega_D_smbh    = 2.0 * UQFF_LA_PI / 1200.0;      // 20-min QPO (rad/s)
        tau_D_smbh      = 3.0 * 3.15576e7;                 // 3-yr QPO envelope decay (s)
        spin_smbh       = 0.3;                              // Kerr spin parameter
        M_ext_smbh      = 2.5e7 * M_sun;                   // NSC ~2.5e7 M_sun (kg)
        r_ext_smbh      = 9.258e16;                         // NSC: 3 pc (m)

        // ── System 3: Starbirth Tapestry (LMC — NGC 2014/2020) ─────────────────────────
        M_tapestry      = 240.0 * M_sun;                   // 240 M_sun initial
        r_tapestry      = 9.461e16;                         // 10 ly (m)
        M_dot_tapestry  = 1.0e4 / 240.0;                   // SF growth factor ~41.7
        tau_SF_tapestry = 5.0e6 * 3.156e7;                 // 5 Myr (s)
        rho_wind_tap    = 1.0e-21;                          // Wind density (kg/m^3)
        v_wind_tap      = 2.0e6;                            // Wind velocity (m/s)
        M_ext_tap       = 1.0e5 * M_sun;                   // R136 Tarantula core (kg)
        r_ext_tap       = 6.15e18;                          // 650 ly (m)

        // ── System 4: Westerlund 2 (super star cluster — Carina arm) ───────────────────
        M_wd2           = 30000.0 * M_sun;                 // 30,000 M_sun initial
        r_wd2           = 9.461e16;                         // 10 ly (m)
        M_dot_wd2       = 1.0e5 / 30000.0;                 // SF growth factor ~3.33
        tau_SF_wd2      = 2.0e6 * 3.156e7;                 // 2 Myr (s)
        rho_wind_wd2    = 1.0e-20;                          // Dense OB-star wind (kg/m^3)
        v_wind_wd2      = 2.0e6;                            // Wind velocity (m/s)
        M_ext_wd2       = 4.0e6 * M_sun;                   // Sgr A* (kg)
        r_ext_wd2       = 2.47e20;                          // 8 kpc (m)

        // ── System 5: Pillars of Creation (Eagle Nebula M16 — Sagittarius arm) ─────────
        M_pillars       = 10100.0 * M_sun;                 // 10,100 M_sun initial
        r_pillars       = 4.731e16;                         // 5 ly (m)
        M_dot_pillars   = 1.0e4 / 10100.0;                 // SF growth factor ~0.99
        tau_SF_pillars  = 1.0e6 * 3.156e7;                 // 1 Myr (s)
        E_0_pillars     = 0.1;                              // Erosion amplitude
        tau_E_pillars   = 1.0e6 * 3.156e7;                 // 1 Myr erosion timescale (s)
        rho_wind_pil    = 1.0e-21;                          // Wind density (kg/m^3)
        v_wind_pil      = 2.0e6;                            // Wind velocity (m/s)
        M_ext_pil       = 4.0e6 * M_sun;                   // Sgr A* (kg)
        r_ext_pil       = 2.11e20;                          // 6.85 kpc (m)

        // ── System 6: Rings of Relativity (GAL-CLUS-022058s — z=0.5 Einstein ring) ─────
        M_rings         = 1.0e14 * M_sun;                  // Galaxy cluster: static (kg)
        r_rings         = 3.086e20;                         // Einstein radius: 10 kpc (m)
        Hz_rings        = 2.42e-18;                         // H(z=0.5) (s^-1)
        L_factor_rings  = 0.67;                             // D_LS/D_S lensing geometry
        z_rings         = 0.5;                              // Lens redshift
        rho_wind_rings  = 1.0e-21;                          // Intracluster wind (kg/m^3)
        v_wind_rings    = 2.0e6;                            // ICM wind velocity (m/s)
        M_ext_rings     = 2.387e45;                         // Virgo Cluster ~1.2e15 M_sun (kg)
        r_ext_rings     = 5.09e23;                          // 16.5 Mpc (m)

        logging_enabled = false;
    }

    // ── Core assessment computations ──────────────────────────────────────────

    // Advancement score (0-100%): weighted mean of 4 metrics
    double compute_advancement() const {
        double score = (diversity_score + dynamic_score + scalability_score + coverage_score)
                       / 4.0 * 100.0;
        if (logging_enabled)
            std::cout << "[LOG] compute_advancement()=" << score
                      << "  div=" << diversity_score << " dyn=" << dynamic_score
                      << " scale=" << scalability_score << " cov=" << coverage_score << std::endl;
        return score;
    }

    // Newtonian g = G*M/r^2 for system idx (0 = SGR1745 ... 6 = RINGS)
    double compute_g_base(int idx) const {
        if (idx < 0 || idx >= NUM_SYSTEMS) {
            std::cerr << "Error: system index " << idx
                      << " out of range [0," << (NUM_SYSTEMS - 1) << "].\n";
            return 0.0;
        }
        double g = g_base_for(idx);
        if (logging_enabled)
            std::cout << "[LOG] g_base[" << idx << " " << system_name(idx)
                      << "]=" << std::scientific << g << " m/s^2\n";
        return g;
    }

    // log10(r_max / r_min) across all 7 systems
    double compute_scale_range() const {
        double r_min = r_sgr1745;   // 20 km — magnetar surface
        double r_max = r_rings;     // 10 kpc — Einstein radius
        double span  = std::log10(r_max / r_min);
        if (logging_enabled)
            std::cout << "[LOG] scale_range: r_min=" << r_min
                      << " r_max=" << r_max << " log10=" << span << std::endl;
        return span;
    }

    // UQFF buoyancy Tier-1: term_Ubi = 0.5 * g_base (CP3/PAPER_198 static half-gravity)
    double compute_Ubi(int idx) const {
        return 0.5 * g_base_for(idx);
    }

    // UQFF buoyancy Tier-2: term_F_UBii = -beta_i * g * omega_g * (M/r) * U_UA * cos(pi*t)
    double compute_F_UBii(int idx, double t) const {
        double g  = g_base_for(idx);
        double Mt = get_M(idx);
        double ri = get_r(idx);
        return -beta_i * g * omega_g * (Mt / ri) * U_UA * std::cos(UQFF_LA_PI * t);
    }

    // UQFF buoyancy Tier-3: term_Ub_i = -beta_i * g * omega_g * (M_ext/r_ext) * U_UA * cos(pi*t)
    double compute_Ub_i(int idx, double t) const {
        double g = g_base_for(idx);
        double M_ext = 0.0, r_ext = 1.0;
        get_ext(idx, M_ext, r_ext);
        return -beta_i * g * omega_g * (M_ext / r_ext) * U_UA * std::cos(UQFF_LA_PI * t);
    }

    // ── Comparative analysis output ───────────────────────────────────────────

    // Print sorted g-value table for all 7 systems at time t (seconds)
    void printComparativeAnalysis(double t = 0.0, std::ostream& os = std::cout) const {
        os << "\n====================================================================\n";
        os << "  UQFF COMPARATIVE ANALYSIS — Research Student Module\n";
        os << "  t = " << std::scientific << std::setprecision(4) << t << " s\n";
        os << "====================================================================\n";

        // Header row
        os << std::left
           << std::setw(4)  << "Idx"
           << std::setw(44) << "System"
           << std::setw(15) << "g_base(m/s^2)"
           << std::setw(15) << "Ubi(m/s^2)"
           << std::setw(15) << "r(m)"
           << "M(kg)" << "\n";
        os << std::string(110, '-') << "\n";

        // Sort by g_base descending
        std::array<int, NUM_SYSTEMS> order;
        for (int i = 0; i < NUM_SYSTEMS; ++i) order[i] = i;
        std::sort(order.begin(), order.end(),
                  [&](int a, int b){ return g_base_for(a) > g_base_for(b); });

        for (int k = 0; k < NUM_SYSTEMS; ++k) {
            int i   = order[k];
            double gb  = g_base_for(i);
            double ubi = compute_Ubi(i);
            os << std::left
               << std::setw(4)  << i
               << std::setw(44) << system_name(i)
               << std::scientific << std::setprecision(4)
               << std::setw(15) << gb
               << std::setw(15) << ubi
               << std::setw(15) << get_r(i)
               << get_M(i) << "\n";
        }
        os << std::string(110, '-') << "\n";
        os << "  Scale range (log10 r_max/r_min): " << std::fixed << std::setprecision(2)
           << compute_scale_range() << " decades\n";
        os << "  Advancement score: " << std::fixed << std::setprecision(1)
           << compute_advancement() << "%\n";
        os << "====================================================================\n\n";
    }

    // Educational walkthrough: unique physics per system, dynamic catalogue, scale analysis
    void printStudentGuide(std::ostream& os = std::cout) const {
        os << "\n===============================================================\n";
        os << "  UQFF Research Student Guide — 7-System MUGE Series\n";
        os << "  Universal Quantum Field Framework | CP3/PAPER_198 Buoyancy\n";
        os << "===============================================================\n\n";

        os << "CORE ARCHITECTURE (applies to ALL systems)\n";
        os << "------------------------------------------\n";
        os << "  MUGE g(r,t) = sum of terms:\n";
        os << "    term1 = g_base * (1+Hz*t) * (1-B/B_crit) * [unique corr per system]\n";
        os << "    term2 = compute_Ug (Ug1+Ug2+Ug3+Ug4) * (1+f_TRZ)\n";
        os << "    term3 = Lambda*c^2/3  (cosmological constant)\n";
        os << "    term4 = EM correction with [UA] vacuum\n";
        os << "    term_q, term_fluid, term_osc, term_DM, term_wind (shared)\n";
        os << "  Buoyancy Pipeline (3 tiers, CP3/PAPER_198):\n";
        os << "    Tier-1  term_Ubi    = 0.5 * g_base           [static half-gravity]\n";
        os << "    Tier-2  term_F_UBii = -b_i * g * w_g * (M/r) * [UA] * cos(pi*t)\n";
        os << "    Tier-3  term_Ub_i   = -b_i * g * w_g * (M_ext/r_ext) * [UA] * cos(pi*t)\n";
        os << std::scientific << std::setprecision(3);
        os << "    b_i=" << beta_i << "  w_g=" << omega_g << " rad/s  [UA]=" << U_UA << " C\n\n";

        for (int i = 0; i < NUM_SYSTEMS; ++i) {
            os << "System [" << i << "] " << system_name(i) << "\n";
            os << "  Regime   : " << regime_label(i) << "\n";
            os << std::scientific << std::setprecision(4);
            os << "  M        : " << get_M(i) << " kg  ("
               << std::fixed << std::setprecision(0) << (get_M(i) / M_sun) << " M_sun)\n";
            os << std::scientific << std::setprecision(4);
            os << "  r        : " << get_r(i) << " m\n";
            os << "  g_base   : " << g_base_for(i) << " m/s^2\n";
            os << "  Ubi      : " << compute_Ubi(i) << " m/s^2\n";
            os << "  Unique   : " << unique_physics_str(i) << "\n\n";
        }

        os << "DYNAMIC PROCESSES CATALOGUE (7 unique time-dependent physics introduced)\n";
        os << "------------------------------------------------------------------------\n";
        os << "  1. D(t) burst   = D0*cos(omega_D*t)*exp(-t/tau_D)  [Systems 0,1 — magnetars]\n";
        os << "  2. QPO burst    = D0*cos(omega_D*t)*exp(-t/tau_D), omega_D=2pi/1200s [System 2 — Sgr A*]\n";
        os << "  3. M(t) growth  = M0*(1+Mdot*exp(-t/tau_SF))       [Systems 2,3,4,5 — accretion/SF]\n";
        os << "  4. E(t) erosion = E0*exp(-t/tau_E) -> corr_E=1-E(t)[System 5 — Pillars]\n";
        os << "  5. H(z) Hubble  = H0*sqrt(Om*(1+z)^3+OL)           [System 6 — Rings, z=0.5]\n";
        os << "  6. L_t lensing  = (G*M/c^2*r)*L_fac -> corr_L=1+L_t [System 6 — Rings]\n";
        os << "  7. Tidal force  = 2*G*M_ext*r/r_ext^3               [Systems 0,1,2 — compact]\n\n";

        os << "SCALE ANALYSIS\n";
        os << "--------------\n";
        os << std::scientific << std::setprecision(3);
        os << "  r_min = " << r_sgr1745 << " m  (magnetar, 20 km — nuclear quantum scale)\n";
        os << "  r_max = " << r_rings   << " m  (Einstein ring, 10 kpc — cosmological scale)\n";
        os << std::fixed << std::setprecision(2);
        os << "  log10(r_max/r_min) = " << compute_scale_range() << " decades\n\n";

        os << "ADVANCEMENT METRICS\n";
        os << "-------------------\n";
        os << std::fixed << std::setprecision(4);
        os << "  diversity_score   = " << diversity_score
           << "  (7/7 distinct physical regimes)\n";
        os << "  dynamic_score     = " << dynamic_score
           << "  (7/10 unique dynamic-physics processes)\n";
        os << "  scalability_score = " << scalability_score
           << "  (16.2/20 log-decades of r coverage)\n";
        os << "  coverage_score    = " << coverage_score
           << "  (7/7 planned C++ modules complete)\n";
        os << std::fixed << std::setprecision(1);
        os << "  OVERALL ADVANCEMENT : " << compute_advancement() << "%\n";
        os << "===============================================================\n\n";
    }

    // Compact parameter dump for all 7 systems
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "UQFFLearningAssessment (Research Student Edition)\n";
        os << "  G=" << G << "  c_light=" << c_light << "  hbar=" << hbar << "\n";
        os << "  beta_i=" << beta_i << "  omega_g=" << omega_g << "  U_UA=" << U_UA << "\n";
        for (int i = 0; i < NUM_SYSTEMS; ++i) {
            os << "  [" << i << "] " << system_name(i)
               << "  M=" << std::scientific << get_M(i)
               << "  r=" << get_r(i)
               << "  g_base=" << g_base_for(i) << "\n";
        }
        os << "  advancement=" << std::fixed << std::setprecision(1) << compute_advancement() << "%\n";
    }

    // Example: print guide + comparative analysis at t = 1 Myr
    double exampleAt1Myr() const {
        const double t = 1.0e6 * 3.156e7;
        printStudentGuide();
        printComparativeAnalysis(t);
        return compute_advancement();
    }

    // ── setVariable / getVariable / addToVariable / subtractFromVariable ──────
    bool setVariable(const std::string& varName, double newValue) {
        if      (varName == "G")               G               = newValue;
        else if (varName == "c_light")         c_light         = newValue;
        else if (varName == "hbar")            hbar            = newValue;
        else if (varName == "Lambda")          Lambda          = newValue;
        else if (varName == "f_TRZ")           f_TRZ           = newValue;
        else if (varName == "beta_i")          beta_i          = newValue;
        else if (varName == "omega_g")         omega_g         = newValue;
        else if (varName == "U_UA")            U_UA            = newValue;
        else if (varName == "diversity_score")   diversity_score   = newValue;
        else if (varName == "dynamic_score")     dynamic_score     = newValue;
        else if (varName == "scalability_score") scalability_score = newValue;
        else if (varName == "coverage_score")    coverage_score    = newValue;
        // SGR1745
        else if (varName == "M_sgr1745")       M_sgr1745       = newValue;
        else if (varName == "r_sgr1745")       r_sgr1745       = newValue;
        else if (varName == "B_sgr1745")       B_sgr1745       = newValue;
        else if (varName == "D0_sgr1745")      D0_sgr1745      = newValue;
        else if (varName == "M_ext_sgr1745")   M_ext_sgr1745   = newValue;
        else if (varName == "r_ext_sgr1745")   r_ext_sgr1745   = newValue;
        // SGR0501
        else if (varName == "M_sgr0501")       M_sgr0501       = newValue;
        else if (varName == "r_sgr0501")       r_sgr0501       = newValue;
        else if (varName == "B_sgr0501")       B_sgr0501       = newValue;
        else if (varName == "D0_sgr0501")      D0_sgr0501      = newValue;
        else if (varName == "M_ext_sgr0501")   M_ext_sgr0501   = newValue;
        else if (varName == "r_ext_sgr0501")   r_ext_sgr0501   = newValue;
        // SMBH Sgr A*
        else if (varName == "M_smbh")          M_smbh          = newValue;
        else if (varName == "r_smbh")          r_smbh          = newValue;
        else if (varName == "M_dot_smbh")      M_dot_smbh      = newValue;
        else if (varName == "spin_smbh")       spin_smbh       = newValue;
        else if (varName == "M_ext_smbh")      M_ext_smbh      = newValue;
        else if (varName == "r_ext_smbh")      r_ext_smbh      = newValue;
        // Tapestry
        else if (varName == "M_tapestry")      M_tapestry      = newValue;
        else if (varName == "r_tapestry")      r_tapestry      = newValue;
        else if (varName == "rho_wind_tap")    rho_wind_tap    = newValue;
        else if (varName == "v_wind_tap")      v_wind_tap      = newValue;
        else if (varName == "M_ext_tap")       M_ext_tap       = newValue;
        else if (varName == "r_ext_tap")       r_ext_tap       = newValue;
        // Westerlund 2
        else if (varName == "M_wd2")           M_wd2           = newValue;
        else if (varName == "r_wd2")           r_wd2           = newValue;
        else if (varName == "rho_wind_wd2")    rho_wind_wd2    = newValue;
        else if (varName == "v_wind_wd2")      v_wind_wd2      = newValue;
        else if (varName == "M_ext_wd2")       M_ext_wd2       = newValue;
        else if (varName == "r_ext_wd2")       r_ext_wd2       = newValue;
        // Pillars
        else if (varName == "M_pillars")       M_pillars       = newValue;
        else if (varName == "r_pillars")       r_pillars       = newValue;
        else if (varName == "E_0_pillars")     E_0_pillars     = newValue;
        else if (varName == "rho_wind_pil")    rho_wind_pil    = newValue;
        else if (varName == "v_wind_pil")      v_wind_pil      = newValue;
        else if (varName == "M_ext_pil")       M_ext_pil       = newValue;
        else if (varName == "r_ext_pil")       r_ext_pil       = newValue;
        // Rings
        else if (varName == "M_rings")         M_rings         = newValue;
        else if (varName == "r_rings")         r_rings         = newValue;
        else if (varName == "Hz_rings")        Hz_rings        = newValue;
        else if (varName == "L_factor_rings")  L_factor_rings  = newValue;
        else if (varName == "z_rings")         z_rings         = newValue;
        else if (varName == "M_ext_rings")     M_ext_rings     = newValue;
        else if (varName == "r_ext_rings")     r_ext_rings     = newValue;
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
        if      (varName == "G")               return G;
        else if (varName == "c_light")         return c_light;
        else if (varName == "hbar")            return hbar;
        else if (varName == "Lambda")          return Lambda;
        else if (varName == "f_TRZ")           return f_TRZ;
        else if (varName == "beta_i")          return beta_i;
        else if (varName == "omega_g")         return omega_g;
        else if (varName == "U_UA")            return U_UA;
        else if (varName == "diversity_score")   return diversity_score;
        else if (varName == "dynamic_score")     return dynamic_score;
        else if (varName == "scalability_score") return scalability_score;
        else if (varName == "coverage_score")    return coverage_score;
        else if (varName == "M_sgr1745")       return M_sgr1745;
        else if (varName == "r_sgr1745")       return r_sgr1745;
        else if (varName == "B_sgr1745")       return B_sgr1745;
        else if (varName == "D0_sgr1745")      return D0_sgr1745;
        else if (varName == "M_ext_sgr1745")   return M_ext_sgr1745;
        else if (varName == "r_ext_sgr1745")   return r_ext_sgr1745;
        else if (varName == "M_sgr0501")       return M_sgr0501;
        else if (varName == "r_sgr0501")       return r_sgr0501;
        else if (varName == "B_sgr0501")       return B_sgr0501;
        else if (varName == "D0_sgr0501")      return D0_sgr0501;
        else if (varName == "M_ext_sgr0501")   return M_ext_sgr0501;
        else if (varName == "r_ext_sgr0501")   return r_ext_sgr0501;
        else if (varName == "M_smbh")          return M_smbh;
        else if (varName == "r_smbh")          return r_smbh;
        else if (varName == "M_dot_smbh")      return M_dot_smbh;
        else if (varName == "spin_smbh")       return spin_smbh;
        else if (varName == "M_ext_smbh")      return M_ext_smbh;
        else if (varName == "r_ext_smbh")      return r_ext_smbh;
        else if (varName == "M_tapestry")      return M_tapestry;
        else if (varName == "r_tapestry")      return r_tapestry;
        else if (varName == "rho_wind_tap")    return rho_wind_tap;
        else if (varName == "v_wind_tap")      return v_wind_tap;
        else if (varName == "M_ext_tap")       return M_ext_tap;
        else if (varName == "r_ext_tap")       return r_ext_tap;
        else if (varName == "M_wd2")           return M_wd2;
        else if (varName == "r_wd2")           return r_wd2;
        else if (varName == "rho_wind_wd2")    return rho_wind_wd2;
        else if (varName == "v_wind_wd2")      return v_wind_wd2;
        else if (varName == "M_ext_wd2")       return M_ext_wd2;
        else if (varName == "r_ext_wd2")       return r_ext_wd2;
        else if (varName == "M_pillars")       return M_pillars;
        else if (varName == "r_pillars")       return r_pillars;
        else if (varName == "E_0_pillars")     return E_0_pillars;
        else if (varName == "rho_wind_pil")    return rho_wind_pil;
        else if (varName == "v_wind_pil")      return v_wind_pil;
        else if (varName == "M_ext_pil")       return M_ext_pil;
        else if (varName == "r_ext_pil")       return r_ext_pil;
        else if (varName == "M_rings")         return M_rings;
        else if (varName == "r_rings")         return r_rings;
        else if (varName == "Hz_rings")        return Hz_rings;
        else if (varName == "L_factor_rings")  return L_factor_rings;
        else if (varName == "z_rings")         return z_rings;
        else if (varName == "M_ext_rings")     return M_ext_rings;
        else if (varName == "r_ext_rings")     return r_ext_rings;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'.\n";
            return 0.0;
        }
    }

    // ── UQFF 2.0 Self-Expanding Framework ────────────────────────────────────

    void setEnableLogging(bool enabled) { logging_enabled = enabled; }
    bool getLoggingEnabled() const      { return logging_enabled; }

    void setDynamicParameter(const std::string& key, double value) {
        dynamic_params[key] = value;
        if (logging_enabled)
            std::cout << "[LOG] Dynamic param set: " << key << " = " << value << std::endl;
    }

    double getDynamicParameter(const std::string& key) const {
        auto it = dynamic_params.find(key);
        if (it != dynamic_params.end()) return it->second;
        if (logging_enabled)
            std::cout << "[LOG] Dynamic param not found: " << key << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }

    void exportState(const std::string& filename = "UQFFLearningAssessment_state.txt") const {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "Error: Cannot open '" << filename << "' for state export.\n";
            return;
        }
        ofs << "# UQFFLearningAssessment (Research Student Edition) state export\n";
        ofs << std::scientific << std::setprecision(10);
        const char* names[] = {
            "G","c_light","hbar","Lambda","f_TRZ","beta_i","omega_g","U_UA",
            "diversity_score","dynamic_score","scalability_score","coverage_score",
            "M_sgr1745","r_sgr1745","B_sgr1745","D0_sgr1745","M_ext_sgr1745","r_ext_sgr1745",
            "M_sgr0501","r_sgr0501","B_sgr0501","D0_sgr0501","M_ext_sgr0501","r_ext_sgr0501",
            "M_smbh","r_smbh","M_dot_smbh","spin_smbh","M_ext_smbh","r_ext_smbh",
            "M_tapestry","r_tapestry","rho_wind_tap","v_wind_tap","M_ext_tap","r_ext_tap",
            "M_wd2","r_wd2","rho_wind_wd2","v_wind_wd2","M_ext_wd2","r_ext_wd2",
            "M_pillars","r_pillars","E_0_pillars","rho_wind_pil","v_wind_pil","M_ext_pil","r_ext_pil",
            "M_rings","r_rings","Hz_rings","L_factor_rings","z_rings","M_ext_rings","r_ext_rings"
        };
        for (const char* n : names)
            ofs << n << " = " << getVariable(n) << "\n";
        for (const auto& kv : dynamic_params)
            ofs << "dynamic." << kv.first << " = " << kv.second << "\n";
        if (logging_enabled)
            std::cout << "[LOG] State exported to: " << filename << std::endl;
    }

    // Cross-validate g_base[sys_idx] against another UQFF module's g at the same time
    // Store the other module's g value via setDynamicParameter("cross_g", value) first.
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
            std::cout << "[XVAL] sys[" << sys_idx << "] g_here=" << g_here
                      << "  g_other=" << g_other
                      << "  frac_diff=" << frac << std::endl;
        return frac;
    }
};

#endif // UQFF_LEARNING_ASSESSMENT_H