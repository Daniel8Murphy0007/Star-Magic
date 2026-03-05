/**
 * @file shared_constants.h
 * @brief Unified Physical Constants for UQFF Star-Magic Framework
 * 
 * This header provides ALL physical constants used across:
 * - MAIN_1_CoAnQi.cpp (C++ production calculator)
 * - source2.cpp (Qt6 GUI head program)
 * - CondensedPhysics.py (via shared_constants.py)
 * - index.js (JavaScript engine)
 * 
 * IMPORTANT: Any constant changes must be synchronized with:
 *   - shared_constants.py (Python)
 *   - index.js CONSTANTS object
 * 
 * Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
 * Framework: UQFF Star-Magic v3.0
 * Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
 */

#ifndef SHARED_CONSTANTS_H
#define SHARED_CONSTANTS_H

#include <cmath>

namespace UQFF {
namespace Constants {

// ═══════════════════════════════════════════════════════════════════════════════
// FUNDAMENTAL PHYSICAL CONSTANTS (CODATA 2018)
// ═══════════════════════════════════════════════════════════════════════════════

/// Speed of light in vacuum (m/s) - exact by definition
constexpr double c = 2.99792458e8;

/// Gravitational constant (m³ kg⁻¹ s⁻²)
constexpr double G = 6.67430e-11;

/// Reduced Planck constant ℏ (J·s)
constexpr double hbar = 1.054571817e-34;

/// Planck constant h (J·s)
constexpr double h_planck = 6.62607015e-34;

/// Elementary charge (C)
constexpr double e = 1.602176634e-19;

/// Boltzmann constant (J/K)
constexpr double k_B = 1.380649e-23;

/// Vacuum permittivity ε₀ (F/m)
constexpr double epsilon_0 = 8.8541878128e-12;

/// Vacuum permeability μ₀ (H/m) = 4π × 10⁻⁷
constexpr double mu_0 = 1.25663706212e-6;

/// Fine structure constant α ≈ 1/137
constexpr double alpha = 7.2973525693e-3;

/// Avogadro constant (mol⁻¹)
constexpr double N_A = 6.02214076e23;

/// Stefan-Boltzmann constant (W m⁻² K⁻⁴)
constexpr double sigma_SB = 5.670374419e-8;

/// Pi (mathematical constant)
constexpr double PI = 3.14159265358979323846;
constexpr double pi = PI;  // lowercase alias
constexpr double TWO_PI = 2.0 * PI;

// ═══════════════════════════════════════════════════════════════════════════════
// PARTICLE MASSES
// ═══════════════════════════════════════════════════════════════════════════════

/// Electron mass (kg)
constexpr double m_e = 9.1093837015e-31;

/// Proton mass (kg)
constexpr double m_p = 1.67262192369e-27;

/// Neutron mass (kg)
constexpr double m_n = 1.67492749804e-27;

/// Atomic mass unit (kg)
constexpr double u = 1.66053906660e-27;


// ═══════════════════════════════════════════════════════════════════════════════
// ASTROPHYSICAL CONSTANTS
// ═══════════════════════════════════════════════════════════════════════════════

/// Solar mass M☉ (kg)
constexpr double M_sun = 1.98892e30;

/// Solar radius R☉ (m)
constexpr double R_sun = 6.9634e8;

/// Solar luminosity L☉ (W)
constexpr double L_sun = 3.828e26;

/// Earth mass M⊕ (kg)
constexpr double M_earth = 5.9722e24;

/// Earth radius R⊕ (m)
constexpr double R_earth = 6.371e6;

/// Astronomical Unit (m)
constexpr double AU = 1.495978707e11;

/// Parsec (m)
constexpr double pc = 3.0856775814914e16;

/// Kiloparsec (m)
constexpr double kpc = 3.0856775814914e19;

/// Megaparsec (m)
constexpr double Mpc = 3.0856775814914e22;

/// Light-year (m)
constexpr double ly = 9.4607304725808e15;

/// Year (s)
constexpr double yr = 3.15576e7;

/// Day (s)
constexpr double day = 86400.0;


// ═══════════════════════════════════════════════════════════════════════════════
// COSMOLOGICAL CONSTANTS
// ═══════════════════════════════════════════════════════════════════════════════

/// Hubble constant H₀ (s⁻¹) - 70 km/s/Mpc
constexpr double H0 = 2.269e-18;

/// Cosmological constant Λ (m⁻²)
constexpr double Lambda = 1.1e-52;

/// Critical density ρ_crit (kg/m³) at z=0
constexpr double rho_crit = 9.47e-27;

/// CMB temperature T_CMB (K)
constexpr double T_CMB = 2.7255;


// ═══════════════════════════════════════════════════════════════════════════════
// UQFF-SPECIFIC CONSTANTS (Murphy Framework)
// ═══════════════════════════════════════════════════════════════════════════════

// ═══════════════════════════════════════════════════════════════════════════════
// VACUUM DENSITY GRADIENT SYSTEM
// ═══════════════════════════════════════════════════════════════════════════════
// 
// The UQFF framework uses TWO vacuum density scales that create a GRADIENT:
//
// 1. GRAVITATIONAL SCALE (rho_vac_UA): 7.09e-36 J/m³
//    - Used in: Ug1-4 equations, cosmological terms, UQFF buoyancy
//    - Represents: Cosmic vacuum energy density (dark energy scale)
//    - Derived from: Λc²/8πG ≈ 5.96e-27 J/m³ modified by SCm factor
//
// 2. FIELD SCALE (rho_vac_UA_field): 1e-27 J/m³
//    - Used in: Electric field terms, neutron production, magnetism
//    - Represents: Local field coupling vacuum density
//    - Derived from: Critical density ρ_c ≈ 9.47e-27 kg/m³
//
// GRADIENT RATIO: rho_vac_UA / rho_vac_UA_field = 7.09e-9
//    - This ~10^9 ratio creates coupling between gravitational and field sectors
//    - The gradient drives energy flow in DPM (Di-Pseudo-Monopole) interactions
//    - Mathematically: ∇ρ_vac = (ρ_UA - ρ_field) / L_transition
//
// DO NOT "FIX" THIS BY UNIFYING - THE GRADIENT IS INTENTIONAL PHYSICS
// ═══════════════════════════════════════════════════════════════════════════════

/// Universal Aether vacuum density ρ_UA (J/m³) - GRAVITATIONAL SCALE
/// Used in: UQFF gravity (Ug1-4), buoyancy (Ub_i), cosmological terms
constexpr double rho_vac_UA = 7.09e-36;

/// Universal Aether vacuum density ρ_UA_field (J/m³) - FIELD SCALE  
/// Used in: Electric field (E = U_m/ρ_field/r), neutron production (η)
constexpr double rho_vac_UA_field = 1e-27;

/// Vacuum density gradient ratio (dimensionless)
/// ∇ρ_ratio = rho_vac_UA / rho_vac_UA_field ≈ 7.09e-9
constexpr double rho_vac_gradient_ratio = rho_vac_UA / rho_vac_UA_field;

/// SCm vacuum density ρ_SCm (J/m³) - Superconductive medium
constexpr double rho_vac_SCm = 7.09e-37;

/// UQFF decay constant κ (day⁻¹)
constexpr double kappa = 0.0005;

/// UQFF solvability [SSq]
constexpr double SSq = 0.57;

/// SCm modulation factor H_SCm
constexpr double H_SCm = 0.99;

/// Universal Aether correction U_UA
constexpr double U_UA = 0.0001;

/// Eta coupling constant k_η
constexpr double k_eta = 1e-113;

/// Beta buoyancy factor β_i (dimensionless)
/// Physical interpretation: Uniform buoyancy opposition strength across all Ug ranges
/// Grok Thread 4e0ecf23: β_i = 0.6 uniformity reflects UQFF superposition principle
/// - Buoyancy opposes gravity uniformly regardless of Ug1-4 dominant range
/// - No range-specific tuning needed due to SCm-UA coupling universality
/// Ref: GrokThread_StarMagic_UnifiedFramework.py (Variable Documentation)
constexpr double beta_i = 0.603;

/// F₀ reference force (N) - UQFF normalization
constexpr double F0 = 1.83e71;

/// Number of magnetic strings (compact objects)
constexpr double num_strings = 1e9;


// ═══════════════════════════════════════════════════════════════════════════════
// UQFF COUPLING CONSTANTS (k_i coefficients)
// ═══════════════════════════════════════════════════════════════════════════════
//
// Physical interpretations from Grok Thread 4e0ecf23:
//
// k_1 = 1.5 (Ug1: Internal Dipole)
//    HIGHER value → emphasizes strong internal stellar irregularities
//    SCm modulation strengthens dipole distortion of gravitational field
//
// k_2 = 1.2 (Ug2: Heliosphere)
//    MODERATE value → balance between wind ram pressure and SCm envelope
//    Heliosphere acts as buffer, less dramatic than dipole or magnetic strings
//
// k_3 = 1.8 (Ug3: Magnetic Strings Disk)
//    HIGHEST value → magnetic strings have largest influence on rotation curves
//    Planetary/stellar cores trap SCm → strongest gravitational distortion
//    Explains why galaxy rotation curves deviate most from Newtonian at this scale
//
// k_4 = 1.0 (Ug4: Star-Black Hole)
//    BASELINE value → normalized reference for largest-scale interactions
//    No SCm penetration modulation (black holes have zero internal structure)
//
// Ref: GrokThread_StarMagic_UnifiedFramework.py (UQFF_VARIABLE_DOCUMENTATION)
// ═══════════════════════════════════════════════════════════════════════════════

/// LENR coupling k_LENR
constexpr double k_LENR = 1e-10;

/// Activation coupling k_act
constexpr double k_act = 1e-14;

/// Dark energy coupling k_DE
constexpr double k_DE = 1e-16;

/// Neutron coupling k_neutron
constexpr double k_neutron = 1e-20;

/// Relativistic coupling k_rel
constexpr double k_rel = 1e-12;

/// Vacuum coupling k_vac
constexpr double k_vac = 1e-10;

/// THz resonance coupling k_thz
constexpr double k_thz = 1e-15;

/// Conduit coupling k_conduit
constexpr double k_conduit = 1e-18;

/// Spooky (entanglement) coupling k_spooky
constexpr double k_spooky = 1e-20;


// ═══════════════════════════════════════════════════════════════════════════════
// UQFF FREQUENCY CONSTANTS
// ═══════════════════════════════════════════════════════════════════════════════

/// LENR resonance frequency ω_LENR (Hz)
constexpr double omega_LENR = 1.2e12;

/// Magnetar rotation frequency Ω_g (rad/s) - SGR 1745 reference
constexpr double Omega_g = 7.3e-16;

/// THz source frequency f_THz (Hz)
constexpr double f_THz = 1e12;


// ═══════════════════════════════════════════════════════════════════════════════
// MAGNETAR / EXTREME FIELD CONSTANTS
// ═══════════════════════════════════════════════════════════════════════════════

/// Critical (Schwinger) magnetic field B_crit (T)
constexpr double B_crit = 4.414e9;

/// Critical magnetar field (T) - SGR type reference
constexpr double B_crit_magnetar = 4.4e13;

/// Solar surface magnetic field (T)
constexpr double B_sun_surface = 1e-4;


// ═══════════════════════════════════════════════════════════════════════════════
// REFERENCE SYSTEM PARAMETERS (SGR 1745-2900)
// ═══════════════════════════════════════════════════════════════════════════════

/// SGR 1745 black hole mass M_bh (kg)
constexpr double M_bh_SGR1745 = 8.155e36;

/// SGR 1745 characteristic distance d_g (m)
constexpr double d_g_SGR1745 = 2.55e20;

/// SGR 1745 stellar wind density ρ_sw (kg/m³)
constexpr double rho_sw_SGR1745 = 8e-21;

/// SGR 1745 stellar wind efficiency ε_sw
constexpr double epsilon_sw_SGR1745 = 0.001;


// ═══════════════════════════════════════════════════════════════════════════════
// DERIVED CONSTANTS
// ═══════════════════════════════════════════════════════════════════════════════

/// Bohr magneton μ_B = eℏ/(2m_e) (J/T)
constexpr double mu_B = 9.2740100783e-24;

/// Electron Compton wavelength λ_C = h/(m_e c) (m)
constexpr double lambda_C = 2.42631023867e-12;

/// Classical electron radius r_e = e²/(4πε₀m_e c²) (m)
constexpr double r_e = 2.8179403262e-15;

/// Hartree energy E_H = m_e e⁴/(16π²ε₀²ℏ²) (J)
constexpr double E_H = 4.3597447222071e-18;

/// Bohr radius a₀ = 4πε₀ℏ²/(m_e e²) (m)
constexpr double a_0 = 5.29177210903e-11;

/// Schwarzschild radius of 1 M☉ = 2GM/c² (m)
constexpr double r_s_sun = 2953.0;


// ═══════════════════════════════════════════════════════════════════════════════
// QUANTUM GRAVITY / BSM CONSTANTS
// ═══════════════════════════════════════════════════════════════════════════════

/// Planck mass m_P (kg)
constexpr double m_P = 2.176434e-8;

/// Planck length l_P (m)
constexpr double l_P = 1.616255e-35;

/// Planck time t_P (s)
constexpr double t_P = 5.391247e-44;

/// Planck temperature T_P (K)
constexpr double T_P = 1.416784e32;

/// Hubble time t_H = 1/H₀ (s) ≈ 13.8 Gyr
constexpr double t_H = 4.35e17;


// ═══════════════════════════════════════════════════════════════════════════════
// INFORMATION PARADOX CONSTANTS (Batch 21)
// ═══════════════════════════════════════════════════════════════════════════════

/// Position uncertainty (26D compact dimensions)
constexpr double Delta_x_p = 1e-68;

/// Wavefunction integral normalization
constexpr double integral_psi = 2.176e-18;


// ═══════════════════════════════════════════════════════════════════════════════
// INFLATION/FORCE CHART EPOCH FRAMEWORK (Grok Thread 4e0ecf23 - March 4, 2026)
// ═══════════════════════════════════════════════════════════════════════════════
//
// 5-Epoch Cosmic Evolution Framework documenting WHEN Ug ranges become active
// in Universal History. This provides temporal context for UQFF calculations.
//
// Source: Grok Thread 4e0ecf23 "Star Magic: The Quest for Unity"
// URL: https://x.com/i/grok/share/4e0ecf23920b435cb3b2f410e93699b5
// Module: GrokThread_StarMagic_UnifiedFramework.py (857 lines)
// Analysis: GROK_THREAD_4E0ECF23_ANALYSIS.md
//
// Background: This is NOT new physics (all Ug1-4 already implemented in codebase).
// Instead, it provides VALIDATION CONTEXT showing which epochs activate which ranges.
//
// ═══════════════════════════════════════════════════════════════════════════════

namespace InflationForceChart {

/// Epoch 1: Fisile Nuclei/Nebular → Periodic Table Formation
/// Time: t = 1.0 - 1.9 (cosmic time scale)
/// SCm State: SCm (initial)
/// Ug Ranges: None active (pre-stellar)
/// Cosmic Structure: Periodic Table forming from fisile nuclei
constexpr int EPOCH_1_FISILE_NUCLEI = 1;
constexpr double EPOCH_1_TIME_START = 1.0;
constexpr double EPOCH_1_TIME_END = 1.9;

/// Epoch 2: Star/Planetary Atom → Ug1-Ug3 Activation
/// Time: t = 2.0 - 2.9
/// SCm State: SCm'' (second-order modulation)
/// Ug Ranges: Ug1, Ug2, Ug3 ACTIVE (star ignition)
/// Cosmic Structure: Stars and Planets (heliospheres form, planetary cores trap SCm)
constexpr int EPOCH_2_STAR_PLANETARY = 2;
constexpr double EPOCH_2_TIME_START = 2.0;
constexpr double EPOCH_2_TIME_END = 2.9;

/// Epoch 3: Galaxies/Quasar → Early Ug4
/// Time: t = 3.0 - 3.9
/// SCm State: SCm''' (third-order modulation)
/// Ug Ranges: Ug1, Ug2, Ug3, Early Ug4 (galaxy formation)
/// Cosmic Structure: Galaxies and Quasars
constexpr int EPOCH_3_GALAXIES_QUASAR = 3;
constexpr double EPOCH_3_TIME_START = 3.0;
constexpr double EPOCH_3_TIME_END = 3.9;

/// Epoch 4: Magnetar/SMBH → Ug4 DOMINANCE
/// Time: t = 4.0 - 4.9
/// SCm State: SCm'''' (fourth-order modulation)
/// Ug Ranges: ALL ACTIVE, Ug4 DOMINATES
/// Cosmic Structure: Magnetars and Supermassive Black Holes (Ug4 signature observable in Sagittarius A* orbits)
/// Validation: Gaia DR4 (2026) should show Ug4 signatures in stellar orbits around Sgr A*
constexpr int EPOCH_4_MAGNETAR_SMBH = 4;
constexpr double EPOCH_4_TIME_START = 4.0;
constexpr double EPOCH_4_TIME_END = 4.9;

/// Epoch 5: Globular Clusters → Stabilization
/// Time: t = 5.0 - 5.9
/// SCm State: SCm''''' (fifth-order modulation)
/// Ug Ranges: ALL ACTIVE, stabilized
/// Cosmic Structure: Globular Clusters (long-term equilibrium)
constexpr int EPOCH_5_GLOBULAR_CLUSTERS = 5;
constexpr double EPOCH_5_TIME_START = 5.0;
constexpr double EPOCH_5_TIME_END = 5.9;

/// Total number of epochs
constexpr int NUM_EPOCHS = 5;

/// F_U epoch calculation baseline (N)
/// F_U(t=0) = F_core + sum_{states=1 to 26} (Ui_state + F_p_state)
/// F_core = ℏ ω_LENR / (σ_n ρ_vac,[UA]) ~ 10^10 N
constexpr double F_U_EPOCH_CORE = 1e10;

} // namespace InflationForceChart


// ═══════════════════════════════════════════════════════════════════════════════
// DPM BIRTH SPHERE GEOMETRY (Grok Thread 4e0ecf23)
// ═══════════════════════════════════════════════════════════════════════════════
//
// Birth of Di-Pseudo-Monopole (DPM) Big Bang mechanism:
//   Sphere equation: (x - h)² + (y - k)² + (z - l)² = r²
//
// Interpretation: 26 quantum states = 26 centers in pre-Big Bang 26-shell EM field
// Each of the 26 quantum levels has a center (h_i, k_i, l_i) forming a geometric
// constellation that collapses into DPM during Big Bang.
//
// Pre-Big Bang: [SCm] and [UA] in vacuum → 26-shell EM field → DPM birth
//
// Ref: GrokThread_StarMagic_UnifiedFramework.py birth_of_dpm_sphere() function
// ═══════════════════════════════════════════════════════════════════════════════

namespace DPMGeometry {

/// Number of sphere centers (one per quantum level)
constexpr int NUM_DPM_CENTERS = 26;

/// DPM sphere characteristic radius (m) - Pre-Big Bang scale
/// Estimated from: sqrt(l_P * t_H) ~ 10^-18 m (Planck-Hubble geometric mean)
constexpr double DPM_SPHERE_RADIUS = 1e-18;

} // namespace DPMGeometry


// ═══════════════════════════════════════════════════════════════════════════════
// "BELLY BUTTON" COSMIC RESONANCE FACTOR (Grok Thread 4e0ecf23)
// ═══════════════════════════════════════════════════════════════════════════════
//
// Pre-Big Bang standing resonance factor:
//   - First foundational constant/source of electrostatic mechanism
//   - [SCm], [UA], electromagnetic, quantum envelope in 26-field ACP_massive
//   - -1/2 states as high energy superconductive barriers
//
// This is the "origin point" of the UQFF - the cosmic resonance that established
// the fundamental ratio a/b relating GM/r², e (elementary charge), and q (charge).
//
// Ref: GrokThread_StarMagic_UnifiedFramework.py BELLY_BUTTON_PARAMS
// ═══════════════════════════════════════════════════════════════════════════════

namespace BellyButtonResonance {

/// Pre-Big Bang resonance frequency (Hz) - estimated from Planck time
/// f_BB = 1/t_P ≈ 1.855 × 10^43 Hz
constexpr double PRE_BIG_BANG_RESONANCE_FREQ = 1.855e43;

/// 26-field envelope coupling (dimensionless)
/// Couples all 26 quantum levels to pre-Big Bang EM shell
constexpr double ACP_MASSIVE_COUPLING = 1.0;

/// High energy barrier for -1/2 states (J)
/// E_barrier = k_B * T_P (Planck temperature barrier)
constexpr double SUPERCONDUCTIVE_BARRIER_ENERGY = 1.9444e9; // k_B * T_P

} // namespace BellyButtonResonance

/// Thread 10220801 (March 5, 2026): Solar UQFF Calibration
/// 10 new CondensedPhysics2.py solar calculators — Ug1Solar through SolarAether
namespace SolarUQFFCalibration {

/// Sagittarius A* black hole mass (kg) — 2025 EHT GRAVITY+ measurement
/// M_bh = (4.297 ± 0.012) × 10^6 M_sun = 8.55 × 10^36 kg
constexpr double M_BH_SGRA_2025 = 8.55e36;

/// Solar cycle angular frequency (rad/s) — 11-year solar activity cycle
/// omega_c = 2π / (11 × 365.25 × 86400) ≈ 1.60 × 10⁻⁸ rad/s
constexpr double OMEGA_SOLAR_CYCLE = 1.6e-8;

/// Solar wind vacuum density at 1 AU (kg/m³)
/// rho_sw_1AU ≈ 8.4 × 10⁻²¹ kg/m³ (average quiet-sun conditions)
constexpr double RHO_SOLAR_WIND_1AU = 8.4e-21;

/// Solar coupling constants (dimensionless) — UQFF magnetic disk calibration
/// k1: poloidal coupling, k2: toroidal, k3: halo interface, k4: disk midplane
constexpr double K1_POLOIDAL_COUPLING = 1.5;
constexpr double K2_TOROIDAL_COUPLING = 1.2;
constexpr double K3_HALO_INTERFACE    = 1.8;
constexpr double K4_DISK_MIDPLANE     = 1.0;

/// Solar wind energy density correction (dimensionless)
/// epsilon_sw = 0.001 — small relativistic correction to Ug3 disk integral
constexpr double EPSILON_SOLAR_WIND = 0.001;

/// Solar wind vacuum polarization depth (dimensionless)
/// delta_sw = 0.01 — polarization of vacuum near heliosphere boundary
constexpr double DELTA_SOLAR_WIND = 0.01;

/// Distance from Solar System to Galactic Centre (m)
/// d_g = 26,000 ly = 2.55 × 10²⁰ m
constexpr double D_GALACTIC_CENTRE = 2.55e20;

/// AGN feedback fraction (dimensionless)
/// f_feedback = 0.20 — fraction of AGN luminosity coupled to ISM
constexpr double F_AGN_FEEDBACK = 0.2;

/// UQFF buoyancy calibration constant (dimensionless) — from Grok 4 Sept 2025
/// beta_i = 0.603 (calibrated against 121 astronomical systems)
constexpr double BETA_UQFF_CALIBRATION = 0.603;

/// Solar heliosphere outer bubble radius (m) — ~100 AU termination shock
/// R_bubble = 100 AU = 1.496 × 10¹³ m
constexpr double R_HELIOSPHERE_100AU = 1.496e13;

} // namespace SolarUQFFCalibration

/// StarMagicCanonical — canonical constants from Star Magic 14Apr2025 derivation
/// Thread 3a469fcc: https://x.com/i/grok/share/3a469fcc1af84841a645c923d15a1f8e
/// Source: Star Magic_14April2025.docx — Primary UQFF derivation by Daniel T. Murphy
namespace StarMagicCanonical {
    /// SCm reactivity decay rate kappa (day^-1) — governs E_react fall-off over stellar lifetime
    /// Calibrated from quasar jet observations; at t=0 E_react ~ 10^46 J/m^3
    constexpr double KAPPA_SCM_DECAY = 0.0005;

    /// B_SCm: superconductive interior field contribution to stellar dipole (T)
    /// Undetectable at surface (Qs = 0), drives Ug3 disk formation via CCW/CW rotation
    constexpr double B_SCM_SUPERCONDUCTIVE = 1.0e3;

    /// gamma_Um: near-lossless Um decay rate (day^-1) — refined from 1e-4 in thread 10220801
    /// Governs the magnetic string decay rate for the Um magnetism term
    constexpr double GAMMA_UM_LOSSLESS = 5.0e-5;

    /// delta_bh: Ug4 black hole field modulation factor (dimensionless)
    /// Scales the SCm volume contribution in Ug4 star-BH interaction
    constexpr double DELTA_BH_MODULATION = 0.1;

    /// P_core_planet: planetary core SCm/UA penetration factor (dimensionless)
    /// Fraction of stellar SCm/UA donated to planetary cores at creation epoch
    constexpr double P_CORE_PLANET = 1.0e-3;

    /// v_UA: Aether bulk velocity (m/s) — used in planetary core Hamiltonian H_UA term
    constexpr double V_UA_AETHER = 1.0e6;

    /// delta_omega: differential rotation amplitude (rad/s)
    /// Equatorial CCW vs coronal CW rotation produces Ug3 disk; delta_omega = 0.4e-6 rad/s
    constexpr double DELTA_OMEGA_ROTATION = 0.4e-6;

    /// omega_defect: Ug1 defect oscillation frequency (rad/s)
    /// Models surface irregularities in the internal dipole field
    constexpr double OMEGA_DEFECT = 0.001;

    /// delta_def_amp: Ug1 defect amplitude (dimensionless)
    constexpr double DELTA_DEFECT_AMP = 0.01;

    /// mu_0: vacuum permeability (H/m) — standard SI value
    constexpr double MU_0_VACUUM = 1.2566e-6;

    /// E_react_t0: reactor efficiency at t=0 (J/m^3 normalised)
    /// E_react_t0 = rho_SCm * v_SCm^2 / rho_A ~ 10^46 (solar calibration)
    constexpr double E_REACT_T0 = 1.0e46;

    /// mu_jet: quasar jet dynamic viscosity (Pa.s) — near-zero for SCm-Aether medium
    /// Used in Navier-Stokes quasar jet model (Millennium Problem connection)
    constexpr double MU_QUASAR_JET_VISCOSITY = 1.0e-35;
} // namespace StarMagicCanonical


} // namespace Constants
} // namespace UQFF

// For backward compatibility - expose commonly used constants at global scope
using UQFF::Constants::G;
using UQFF::Constants::c;
using UQFF::Constants::hbar;
using UQFF::Constants::M_sun;
using UQFF::Constants::epsilon_0;
using UQFF::Constants::mu_0;

#endif // SHARED_CONSTANTS_H
