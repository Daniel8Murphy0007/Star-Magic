/**
 * @file SystemCatalogue.cpp
 * @brief Implementation of UQFF master equation catalogue from Source10.cpp
 * @date 2025-11-14
 * @version 1.0
 * 
 * EXTRACTED FROM: Source10.cpp (lines 759-1009)
 * PURPOSE: Centralize all master equations and astrophysical system database
 * PART OF: Phase 1, Week 1 - Source10 extraction to Core/
 */

#include "SystemCatalogue.hpp"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>
#include <array> // MSVC requirement

namespace UQFFCatalogue {

// ============================================================================
// SYSTEM CATALOGUE DATABASE (26+ astrophysical systems)
// Extracted from Source10.cpp lines 759-820
// ============================================================================

// Helper macro to create SystemParams with field assignments
#define CREATE_SYSTEM(NAME, M_VAL, R_VAL, T_VAL, LX_VAL, B0_VAL, OMEGA_VAL, THETA_VAL, TIME_VAL, V_VAL, RHO_UA, RHO_SCM, DPM_S, DPM_M, DPM_G, K_LENR_V, K_ACT_V, K_DE_V, K_NEUT_V, SIGMA_V, K_REL_V, F_REL_V, K_VAC_V, K_THZ_V, OMEGA_THZ_V, N_FACT, C_SCALE, K_COND_V, WATER_V, H_ABUND_V, K_SPOOKY_V, STRING_V, Q_WAVE_V, DELTA_K, V_VOID_V, ALPHA_V, TERM1_V, TERM2_V, TERM3_V, TERM4_V, STD_V, DPM_L, RHO_A, RHO_L) \
    do { \
        SystemParams p; \
        p.name = (NAME); p.M = (M_VAL); p.r = (R_VAL); p.T = (T_VAL); p.L_X = (LX_VAL); p.B0 = (B0_VAL); p.omega0 = (OMEGA_VAL); \
        p.theta_deg = (THETA_VAL); p.t = (TIME_VAL); p.v = (V_VAL); \
        p.rho_vac_UA = (RHO_UA); p.rho_vac_SCm = (RHO_SCM); \
        p.DPM_stability = (DPM_S); p.DPM_momentum = (DPM_M); p.DPM_gravity = (DPM_G); \
        p.k_LENR = (K_LENR_V); p.k_act = (K_ACT_V); p.k_DE = (K_DE_V); p.k_neutron = (K_NEUT_V); p.sigma_n = (SIGMA_V); \
        p.k_rel = (K_REL_V); p.F_rel = (F_REL_V); \
        p.k_vac = (K_VAC_V); p.k_thz = (K_THZ_V); p.omega_thz = (OMEGA_THZ_V); \
        p.neutron_factor = (N_FACT); p.conduit_scale = (C_SCALE); \
        p.k_conduit = (K_COND_V); p.water_state = (WATER_V); p.H_abundance = (H_ABUND_V); \
        p.k_spooky = (K_SPOOKY_V); p.string_wave = (STRING_V); p.Q_wave = (Q_WAVE_V); \
        p.Delta_k_eta = (DELTA_K); p.V_void_fraction = (V_VOID_V); \
        p.alpha_i = (ALPHA_V); \
        p.term1 = (TERM1_V); p.term2 = (TERM2_V); p.term3 = (TERM3_V); p.term4 = (TERM4_V); \
        p.std_scale = (STD_V); p.DPM_life = (DPM_L); \
        p.rho_astro = (RHO_A); p.rho_LEP = (RHO_L); \
        catalogue[(NAME)] = p; \
    } while(0)

std::map<std::string, SystemParams> initializeSystemCatalogue() {
    std::map<std::string, SystemParams> catalogue;

    // NGC 1365 (Galaxy): Mass ~1e11 Msun, r ~1.54e21 m, T ~1e4 K, L_X ~1e33 W, B0 ~1e-9 T, omega0 ~1.95e-16 s^-1, v ~3e5 m/s
    CREATE_SYSTEM("NGC 1365", 1e11 * Msun, 1.54e21, 1e4, 1e33, 1e-9, 1.95e-16, 45.0, 1e15, 3e5, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    // Vela Pulsar: Mass 1.4 Msun, r 1e4 m, T 1e6 K, L_X 1e26 W, B0 3.4e8 T, omega0 70.6 s^-1, v 6.1e4 m/s
    CREATE_SYSTEM("Vela Pulsar", 1.4 * Msun, 1e4, 1e6, 1e26, 3.4e8, 70.6, 45.0, 1e10, 6.1e4, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-17, 1e-25);

    // ASASSN-14li (TDE): BH mass ~1e6 Msun ~1.989e36 kg, r ~3e9 m, T ~1e5 K, L_X ~1e37 W, B0 ~1e-3 T, omega0 ~0, v ~3e7 m/s
    CREATE_SYSTEM("ASASSN-14li", 1e6 * Msun, 3e9, 1e5, 1e37, 1e-3, 0.0, 45.0, 1e15, 3e7, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-17, 1e-25);

    // El Gordo (Cluster): Mass 2e15 Msun ~3.978e45 kg, r ~3.086e22 m, T 1.68e8 K, L_X 2.36e38 W, B0 ~1e-10 T, omega0 ~0, v ~1.3e6 m/s
    CREATE_SYSTEM("El Gordo", 2e15 * Msun, 3.086e22, 1.68e8, 2.36e38, 1e-10, 0.0, 45.0, 1e15, 1.3e6, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    // Magnetar SGR 1745-2900: Mass 1.4 Msun, r 1e4 m, T 1e6 K, L_X 1e28 W, B0 2e10 T, omega0 1.67 s^-1, v 1.3e5 m/s
    CREATE_SYSTEM("Magnetar SGR 1745-2900", 1.4 * Msun, 1e4, 1e6, 1e28, 2e10, 1.67, 45.0, 1e10, 1.3e5, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-17, 1e-25);

    // Tapestry of Blazing Starbirth NGC 2264 (Cluster): Mass ~500 Msun, r ~6.172e16 m, T ~1e4 K, L_X ~1e30 W, B0 ~1e-9 T, omega0 ~0, v ~1e4 m/s
    CREATE_SYSTEM("Tapestry of Blazing Starbirth NGC 2264", 500 * Msun, 6.172e16, 1e4, 1e30, 1e-9, 0.0, 45.0, 1e10, 1e4, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    // Westerlund 2 (Cluster): Mass ~1e4 Msun, r ~3.086e16 m, T ~1e4 K, L_X ~1e32 W, B0 ~1e-9 T, omega0 ~0, v ~5e3 m/s
    CREATE_SYSTEM("Westerlund 2", 1e4 * Msun, 3.086e16, 1e4, 1e32, 1e-9, 0.0, 45.0, 1e10, 5e3, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    // Pillars of Creation M16 (Nebula): Mass ~200 Msun, r ~3.086e16 m, T ~1e4 K, L_X ~1e30 W, B0 ~1e-8 T, omega0 ~0, v ~5e3 m/s
    CREATE_SYSTEM("Pillars of Creation M16", 200 * Msun, 3.086e16, 1e4, 1e30, 1e-8, 0.0, 45.0, 1e10, 5e3, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    // Rings of Relativity (Einstein rings): Typical lens mass ~1e12 Msun, r ~3.086e20 m, T ~1e4 K, L_X ~1e35 W, B0 ~1e-10 T, omega0 ~6.48e-16 s^-1, v ~2e5 m/s
    CREATE_SYSTEM("Rings of Relativity", 1e12 * Msun, 3.086e20, 1e4, 1e35, 1e-10, 6.48e-16, 45.0, 1e15, 2e5, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    // Chandra Archive Collection: Average values
    CREATE_SYSTEM("Chandra Archive Collection", 1e30, 1e16, 1e7, 1e30, 1e-9, 0.0, 45.0, 1e10, 1e6, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    // Cassiopeia (Cas A SNR): Mass ~4 Msun ejected, r ~1.54e17 m, T ~1e7 K, L_X ~1e30 W, B0 ~1e-9 T, omega0 ~0, v ~5e6 m/s
    CREATE_SYSTEM("Cassiopeia", 4 * Msun, 1.54e17, 1e7, 1e30, 1e-9, 0.0, 45.0, 1e10, 5e6, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    // New diversified systems from Chandra deepsearch
    CREATE_SYSTEM("3C273", 1e9 * Msun, 4.6e21, 1e7, 1e37, 1e-5, 1e-15, 45.0, 1e15, 2.7e8, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    CREATE_SYSTEM("Cen A AGN", 1e8 * Msun, 3e13, 1e7, 1e36, 1e-6, 1e-12, 45.0, 1e15, 3e7, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    CREATE_SYSTEM("UHZ1 AGN", 1e7 * Msun, 1e12, 1e8, 1e38, 1e-6, 1e-12, 45.0, 1e15, 3e7, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    CREATE_SYSTEM("Geminga", 1.4 * Msun, 1e4, 1e6, 1e26, 1.6e8, 26.5, 45.0, 1e10, 3.4e5, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-17, 1e-25);

    CREATE_SYSTEM("GW170817", 2.7 * Msun, 2e4, 1e10, 1e32, 1e11, 1e3, 45.0, 1e8, 6e7, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-17, 1e-25);

    // New integrated systems from Chandra deepsearch
    CREATE_SYSTEM("NGC 1068", 1e7 * Msun, 3e16, 1e7, 1e36, 1e-5, 1e-14, 45.0, 1e15, 1e6, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    CREATE_SYSTEM("PJ352-15", 1e9 * Msun, 4.6e21, 1e7, 1e37, 1e-5, 1e-15, 45.0, 1e15, 2.7e8, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    CREATE_SYSTEM("Quasar Survey (Typical)", 1e8 * Msun, 1e13, 1e7, 1e36, 1e-6, 1e-12, 45.0, 1e15, 3e8, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    CREATE_SYSTEM("GSN 069", 4e5 * Msun, 1e9, 1e5, 1e32, 1e8, 1e-13, 45.0, 1e15, 1e7, 7.09e-36, 7.09e-37, 0.01, 0.93, 1.0, 1e-10, 1e-6, 1e-30, 1e10, 1e-4, 1e-10, 4.30e33, 1e-30, 1e-10, 2 * PI * 1e12, 1.0, 10.0, 1e-22, 1.0, 10.0, 1e-30, 1e-10, 1.0, 7.25e8, 0.2, 0.01, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 1e-24, 1e-25);

    return catalogue;
}

#undef CREATE_SYSTEM

// ============================================================================
// HELPER FUNCTIONS
// Extracted from Source10.cpp lines 823-869
// ============================================================================

/**
 * @brief Compute E_cm scaling factor
 * Long-form: E_cm,astro,local,adj,eff,enhanced = E_LEP * sqrt(rho_astro / rho_LEP) * Q_wave
 */
double compute_E_cm(const SystemParams &p) {
    double sqrt_ratio = std::sqrt(p.rho_astro / p.rho_LEP); // Density ratio square root
    return E_LEP * sqrt_ratio * p.Q_wave;                   // Scaled energy
}

/**
 * @brief Compute DPM life span proportion
 * Proportion [(SCM):(UA'):(F_U_Bi_i):(Belly Button)] as a/b/c/d
 */
double dpm_life_proportion(const SystemParams &p) {
    double SCM = p.rho_vac_SCm;           // Superconductive magnetism density
    double UA_prime = p.rho_vac_UA;       // Adjusted UA
    double F_U_Bi = p.F_U_Bi_i;           // Buoyancy
    double belly_button = p.omega0 * p.r; // Placeholder for torque-related "belly button" term
    
    // Ratio (SCM / UA_prime) : (F_U_Bi / belly_button)
    double ratio1 = SCM / UA_prime;
    double ratio2 = F_U_Bi / belly_button;
    return ratio1 / ratio2; // Combined proportion
}

// ============================================================================
// MASTER EQUATIONS: F_U_Bi_i (Buoyancy Force)
// Extracted from Source10.cpp lines 836-903
// ============================================================================

/**
 * @brief Compute F_U_Bi_i (buoyancy force, integrating all UQFF terms)
 * 
 * PHYSICS TERMS INTEGRATED:
 * - Vacuum repulsion (UA - SCM density difference)
 * - THz shock (frequency ratio, neutron factor)
 * - Conduit (H abundance, water state)
 * - Spooky action (quantum wave normalization)
 * - LENR (1.25 THz resonance)
 * - Activation (300 Hz)
 * - Directed energy (L_X based)
 * - Resonance (omega0 oscillation)
 * - Neutron drop (sigma_n, neutron_factor)
 * - Relativistic coherence (F_rel)
 * - 26-layer amplification
 * - E_cm scaling
 * - Probabilistic Monte Carlo integration
 */
double F_U_Bi_i(const SystemParams &p) {
    std::srand(std::time(NULL)); // Seed for MC
    double randn = (std::rand() % 1000 / 1000.0 - 0.5) * 2 * std::sqrt(3) * p.std_scale; // Approx normal N(0,1)

    // Vacuum repulsion term
    double Delta_rho_vac = p.rho_vac_UA - p.rho_vac_SCm;
    double F_vac_rep = p.k_vac * Delta_rho_vac * p.M * p.v;

    // THz shock term
    double freq_ratio_sq = std::pow(p.omega_thz / p.omega0, 2);
    double F_thz_shock = p.k_thz * freq_ratio_sq * p.neutron_factor * p.conduit_scale;

    // Conduit term
    double material_interact = p.H_abundance * p.water_state;
    double F_conduit = p.k_conduit * material_interact * p.neutron_factor;

    // Spooky action term
    double wave_norm = p.string_wave / p.omega0;
    double F_spooky = p.k_spooky * wave_norm;

    // Additional core terms
    double LENR_term = p.k_LENR * (1.25e12);                          // Average 1.2-1.3 THz resonance
    double act_term = p.k_act * 300.0;                                // 300 Hz activation
    double DE_term = p.k_DE * (p.L_X / (4 * PI * p.r * p.r));         // Directed energy from luminosity
    double resonance_term = p.k_act * std::cos(p.omega0 * p.t);       // Simple resonance
    double neutron_term = p.k_neutron * p.sigma_n * p.neutron_factor; // Neutron drop
    double rel_term = p.k_rel * p.F_rel;                              // Relativistic coherence

    // Sum all terms
    double F_sum = F_vac_rep + F_thz_shock + F_conduit + F_spooky + 
                   LENR_term + act_term + DE_term + resonance_term + 
                   neutron_term + rel_term;

    // Apply layered scaling for 26 layers
    double layered_F = F_sum * layer_scale_factor;

    // Alternative buoyancy form from hydrogen calc (for small systems)
    double V_total = (4.0 / 3.0) * PI * std::pow(p.r, 3); // Assume spherical
    double V_void = p.V_void_fraction * V_total;
    double g_base = G * p.M / (p.r * p.r); // Base gravity
    double U_Bi_alt = 0.1 * p.Delta_k_eta * (p.rho_vac_UA / p.rho_vac_SCm) * 
                      (V_void / V_total) * g_base;
    
    // If negative, adjust sign
    if (layered_F > 0 && U_Bi_alt < 0)
        layered_F += U_Bi_alt * p.M; // To force, example integration

    // Diversify with GW ripple (currently zero in buoyancy)
    double gw_ripple = 0.0;
    layered_F += gw_ripple;

    // Probabilistic integration
    layered_F *= (1 + randn); // F = mean + std * randn

    // Multi-scale scalar refinement for E_cm
    double E_cm = compute_E_cm(p);

    return layered_F * E_cm; // Integrate E_cm scaling
}

// ============================================================================
// MASTER EQUATIONS: compressed_g (Gravity Field)
// Extracted from Source10.cpp lines 906-944
// ============================================================================

/**
 * @brief Compute compressed_g (sum over 26 quantum layers)
 * 
 * PHYSICS TERMS PER LAYER:
 * - Ug1_i: Dipole/spin term (E_DPM, rho_vac_UA, time-reversal zone)
 * - Ug2_i: Superconductor term (E_DPM, SCM, cosmological communication)
 * - Ug3_i: Resonance term (h_bar * omega0, quantum factor, cos oscillation)
 * - Ug4i_i: Modified Newtonian (G * M, alpha_i correction, SCM enhancement)
 * 
 * All terms summed over 26 layers with 1/i radius scaling
 */
double compressed_g(const SystemParams &p) {
    double g_total = 0.0;
    
    for (int i = 1; i <= num_layers; ++i) {
        // Long-form per layer
        double r_i = p.r / i;     // Scale radius
        double Q_i = i;           // Quantum factor
        double SCm_i = i * i;     // Superconductive magnetism
        double f_TRZ_i = 1.0 / i; // Time-reversal zone factor
        double f_Um_i = i;        // Cosmological communication

        // E_DPM,i
        double r_i_sq = r_i * r_i;
        double E_DPM_i = (h_bar * c / r_i_sq) * Q_i * SCm_i;

        // Ug1_i: Dipole/spin term
        double Ug1_i = (E_DPM_i / r_i_sq) * p.rho_vac_UA * f_TRZ_i;

        // Ug2_i: Superconductor term
        double Ug2_i = (E_DPM_i / r_i_sq) * SCm_i * f_Um_i;

        // Ug3_i: Resonance term
        double f_i = p.omega0 / (2 * PI); // Frequency
        double cos_term = std::cos(2 * PI * f_i * p.t);
        double Ug3_i = (h_bar * p.omega0 / 2.0) * Q_i * cos_term / r_i;

        // Ug4i_i: Modified Newtonian term
        double M_i = p.M / i; // Scaled mass
        double Ug4i_i = (G * M_i / r_i_sq) * (1.0 + p.alpha_i) * SCm_i;

        // Sum per layer
        double layer_g = Ug1_i + Ug2_i + Ug3_i + Ug4i_i;
        g_total += layer_g;
    }
    
    return g_total;
}

// ============================================================================
// RELATIVISTIC FUNCTIONS
// Extracted from Source10.cpp lines 947-970
// ============================================================================

/**
 * @brief Relativistic jet force
 */
double F_jet_rel(const SystemParams &p) {
    double gamma = 1.0 / std::sqrt(1 - std::pow(p.v / c, 2));
    return p.k_thz * std::pow(p.omega_thz / p.omega0, 2) * p.neutron_factor * 
           p.conduit_scale * (p.v / c) * gamma * gamma;
}

/**
 * @brief Relativistic acceleration energy
 */
double E_acc_rel(const SystemParams &p) {
    double beta = p.v / c;
    return (p.L_X / (4 * PI * p.r * p.r * c)) * (1 + beta); // Simplified from E_cm * term
}

/**
 * @brief Relativistic drag force
 */
double F_drag_rel(const SystemParams &p) {
    return p.k_vac * (p.rho_vac_UA - p.rho_vac_SCm) * p.M * p.v * 
           (std::pow(p.B0, 2) / (2 * 4 * PI * 1e-7)) / (p.rho_vac_UA * c);
}

/**
 * @brief Relativistic gravitational wave force
 */
double F_gw_rel(const SystemParams &p) {
    return 0.0; // No GW in buoyancy, set to zero
}

// ============================================================================
// VALIDATION PIPELINE (Simulation)
// Extracted from Source10.cpp lines 973-978
// ============================================================================

/**
 * @brief Validation Pipeline Simulation
 * Prints cross-ref suggestions (no real API)
 */
void validation_pipeline(const SystemParams &p) {
    std::cout << "Simulated Chandra/GW cross-ref for " << p.name << ":" << std::endl;
    std::cout << "Cross-ref L_X with GW strain: " << p.L_X * h_gw << " W (adjusted)" << std::endl;
    std::cout << "Suggest observation: JWST for buoyancy offset ~" << p.v / c * p.r << " m" << std::endl;
}

} // namespace UQFFCatalogue
