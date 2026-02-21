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

/// Beta buoyancy factor β_i
constexpr double beta_i = 0.603;

/// F₀ reference force (N) - UQFF normalization
constexpr double F0 = 1.83e71;

/// Number of magnetic strings (compact objects)
constexpr double num_strings = 1e9;


// ═══════════════════════════════════════════════════════════════════════════════
// UQFF COUPLING CONSTANTS
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
