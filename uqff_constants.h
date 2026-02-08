// ============================================================================
// UQFF_CONSTANTS.H - Unified Physical Constants for Star-Magic Project
// ============================================================================
// Purpose: Shared physical constants between source2.cpp (GUI) and source4.cpp (compute)
// Created: January 20, 2026
// Watermark: ©2025-2026 Daniel T. Murphy, daniel.murphy00@gmail.com – All Rights Reserved
// ============================================================================

#pragma once

#ifndef UQFF_CONSTANTS_H
#define UQFF_CONSTANTS_H

namespace UQFF {

// ============================================================================
// FUNDAMENTAL PHYSICAL CONSTANTS (SI Units)
// ============================================================================

// Gravitational constant (m³/kg·s²) - CODATA 2018
constexpr double G = 6.67430e-11;

// Speed of light in vacuum (m/s) - exact by SI definition
constexpr double c = 2.99792458e8;

// Reduced Planck constant ℏ (J·s) - CODATA 2018
constexpr double hbar = 1.054571817e-34;

// Planck constant h (J·s)
constexpr double h_planck = 6.62607015e-34;

// Pi - mathematical constant
constexpr double PI = 3.14159265358979323846;

// Vacuum permeability μ₀ (H/m or N/A²)
constexpr double mu_0 = 1.25663706212e-6;

// Vacuum permittivity ε₀ (F/m)
constexpr double epsilon_0 = 8.8541878128e-12;

// Boltzmann constant k_B (J/K) - exact by SI definition
constexpr double k_B = 1.380649e-23;

// Elementary charge (C) - exact by SI definition
constexpr double e_charge = 1.602176634e-19;

// Electron mass (kg)
constexpr double m_electron = 9.1093837015e-31;

// Proton mass (kg)
constexpr double m_proton = 1.67262192369e-27;

// Stefan-Boltzmann constant σ (W/m²·K⁴)
constexpr double sigma_SB = 5.670374419e-8;

// Fine structure constant α (dimensionless)
constexpr double alpha_fine = 7.2973525693e-3;

// ============================================================================
// ASTRONOMICAL CONSTANTS
// ============================================================================

// Solar mass (kg) - IAU nominal value
constexpr double SUN_MASS_KG = 1.98892e30;
constexpr double M_sun = SUN_MASS_KG;  // Alias for UQFF equations

// Earth mass (kg)
constexpr double EARTH_MASS_KG = 5.9722e24;
constexpr double M_earth = EARTH_MASS_KG;  // Alias

// Jupiter mass (kg)
constexpr double JUPITER_MASS_KG = 1.89813e27;
constexpr double M_jupiter = JUPITER_MASS_KG;  // Alias

// Solar radius (m)
constexpr double SUN_RADIUS_M = 6.9634e8;
constexpr double R_sun = SUN_RADIUS_M;  // Alias

// Earth radius (m)
constexpr double EARTH_RADIUS_M = 6.371e6;
constexpr double R_earth = EARTH_RADIUS_M;  // Alias

// Astronomical Unit (m) - IAU 2012 exact
constexpr double AU_TO_METERS = 1.495978707e11;
constexpr double AU = AU_TO_METERS;  // Alias

// Parsec (m)
constexpr double PARSEC_TO_METERS = 3.0856775814913673e16;
constexpr double pc = PARSEC_TO_METERS;  // Alias

// Light year (m)
constexpr double LIGHT_YEAR_TO_METERS = 9.4607304725808e15;
constexpr double ly = LIGHT_YEAR_TO_METERS;  // Alias

// Kiloparsec (m)
constexpr double KPC_TO_METERS = 3.0856775814913673e19;
constexpr double kpc = KPC_TO_METERS;  // Alias

// Megaparsec (m)
constexpr double MPC_TO_METERS = 3.0856775814913673e22;
constexpr double Mpc = MPC_TO_METERS;  // Alias

// Julian year (days)
constexpr double JULIAN_YEAR_DAYS = 365.25;

// Sidereal day (seconds)
constexpr double SIDEREAL_DAY_SECONDS = 86164.0905;

// Age of universe (seconds) - Planck 2018
constexpr double AGE_OF_UNIVERSE_S = 4.35e17;

// Observable universe radius (m)
constexpr double OBSERVABLE_UNIVERSE_RADIUS_M = 8.8e26;

// Hubble constant H₀ (s⁻¹) - Planck 2018: 67.4 km/s/Mpc
constexpr double H0 = 2.184e-18;

// Critical density of universe (kg/m³)
constexpr double rho_critical = 9.47e-27;

// Cosmological constant Λ (m⁻²) - approximate
constexpr double Lambda_cosm = 1.1056e-52;

// ============================================================================
// UQFF-SPECIFIC CONSTANTS (Star-Magic Framework)
// ============================================================================

// SCm (String/Cosmic medium) velocity (m/s)
constexpr double v_SCm = 1.0e5;

// Aether density ρ_A (kg/m³) - speculative
constexpr double rho_A = 1e-21;

// Kappa coupling constant
constexpr double kappa = 1e-5;

// Critical magnetic field for magnetars (T)
constexpr double B_crit_magnetar = 4.4e13;

// Vacuum energy density (J/m³) - QFT estimate
constexpr double rho_vacuum = 1e-9;

// String tension parameter
constexpr double T_string = 1e-6;

// Quantum vacuum fluctuation scale
constexpr double delta_vacuum = 1e-15;

// ============================================================================
// MUGE (Multi-system Universal Gravity Equation) CONSTANTS
// ============================================================================

// Resonance parameters
constexpr double f_TRZ_default = 1e12;      // THz resonance frequency
constexpr double omega_quantum = 1e15;       // Quantum frequency (rad/s)
constexpr double omega_aether = 1e10;        // Aether resonance frequency (rad/s)
constexpr double omega_fluid = 1e6;          // Fluid frequency (rad/s)
constexpr double omega_exp = 1e3;            // Expansion frequency (rad/s)

// Dark matter halo parameters
constexpr double rho_DM_local = 0.3;         // Local DM density (GeV/cm³)
constexpr double r_s_NFW = 20e3 * pc;        // NFW scale radius (m)

// ============================================================================
// CONVERSION FACTORS
// ============================================================================

// Angle conversions
constexpr double DEG_TO_RAD = PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / PI;
constexpr double HOURS_TO_RAD = PI / 12.0;
constexpr double RAD_TO_HOURS = 12.0 / PI;

// Time conversions
constexpr double DAYS_TO_SECONDS = 86400.0;
constexpr double YEARS_TO_SECONDS = JULIAN_YEAR_DAYS * DAYS_TO_SECONDS;

// Energy conversions
constexpr double eV_TO_JOULES = 1.602176634e-19;
constexpr double JOULES_TO_eV = 1.0 / eV_TO_JOULES;

// ============================================================================
// VALIDATION BOUNDS
// ============================================================================

// Maximum physically reasonable values for validation
constexpr double MAX_MASS_KG = 1e54;              // Observable universe mass estimate
constexpr double MAX_DISTANCE_M = 8.8e26;         // Observable universe radius
constexpr double MAX_VELOCITY_MS = c;             // Speed of light limit
constexpr double MAX_TIME_S = 4.35e17;            // Age of universe
constexpr double MIN_MASS_KG = m_electron;        // Electron mass as minimum

} // namespace UQFF

#endif // UQFF_CONSTANTS_H
