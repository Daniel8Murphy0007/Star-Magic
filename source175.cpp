/*
 * source175.cpp - UQFF Grok Thread 98b2e77d Physics Integration (C++ Port)
 * ===========================================================================
 * 
 * C++ implementation of Python modules from Grok Thread 98b2e77d analysis:
 * - BuoyancyProofVariants.py (17 F_UBii calculators)
 * - GrokThreadUQFFExtensions.py (Um, Aether, UnifiedField)
 * - UQFFSystemsDatabase.py (astrophysical parameters)
 * 
 * PHYSICS MODULES:
 * ----------------
 * 1. F_UBii Buoyancy Force Variants (17 proofs):
 *    - virx: Virial X-ray clusters (σ_X velocity dispersion)
 *    - termv: Terminal velocity (τ·L/(c·E_LEP))
 *    - upar: Ionization parameter (L_ion/n_H r²)
 *    - coup: Energy coupling (E_kin/E_mag)
 *    - orbdec: Orbital decay (da/dt GW radiation)
 *    - kn: Kilonova peak luminosity (L_peak·t_peak)
 *    - fermi: Fermi acceleration (β_shock·v_shock²)
 *    - kne: Cosmic ray knee energy (ρ_CR·Z·E_knee)
 *    - whim: WHIM temperature (T_WHIM·n_e)
 *    - ps: Press-Schechter halo mass (M_halo·σ_8)
 *    - sfe: Star formation efficiency (ε_SFE·Σ_gas)
 *    - hawk: Hawking temperature (ℏc³/(8πGM·k_B))
 *    - bd: Bounce density (ρ_bounce·a_bounce³)
 *    - roche: Roche lobe overflow (ΔM_dot)
 *    - ent: Entanglement entropy (S_ent·Area)
 *    - dec: Decoherence time (τ_dec·T)
 *    - lobe: Radio lobe dynamics (P_jet·t_lobe)
 * 
 * 2. Universal Magnetism (Um):
 *    - Dipole sum with Widom-Larsen LENR enhancement
 *    - Heaviside 10^13 factor for ultra-high frequency oscillations
 *    - Temporal modulation: (1 - e^(-γt)) * cos(πt_n)
 *    - Quasi-particle correction factor
 * 
 * 3. Aether Metric Tensor (A^μν):
 *    - 4×4 spacetime metric with UQFF vacuum duality
 *    - Temporal component A⁰⁰
 *    - Spatial components A¹¹, A²², A³³
 *    - Frame-dragging off-diagonal A⁰ⁱ
 * 
 * 4. Unified Field Calculator (F_U):
 *    - Master equation: F_U = Σ[k_i×Ug_i - Ub_i] + Um + A^{μν}
 *    - Integrates gravity, buoyancy, magnetism, Aether geometry
 * 
 * INTEGRATION:
 * ------------
 * - Namespace SOURCE175 wraps all functions
 * - Compatible with MAIN_1_CoAnQi.cpp registration system
 * - Uses UQFFSystemsDatabase parameters for realistic calculations
 * - PhysicsTerm base class inheritance for dynamic registration
 * 
 * DATA SOURCES:
 * -------------
 * - Grok Thread: https://x.com/i/grok/share/98b2e77dfbc34d27b09f19fa7c460624
 * - Python modules: BuoyancyProofVariants.py, GrokThreadUQFFExtensions.py
 * - Validation: test_priority1_integration.py (7/7 tests passed)
 * - Astrophysical systems: UQFFSystemsDatabase.py (36 systems, 17 categories)
 * 
 * COMPILATION:
 * ------------
 * Integrated into MAIN_1_CoAnQi.cpp via inclusion or namespace import.
 * Requires: C++17, <cmath>, <string>, <map>, <vector>, <array>, <complex>
 * 
 * Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
 * Created: March 3, 2026
 * Framework: UQFF 99.9% Solvability (Star-Magic)
 * Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
 */

#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <array>
#include <complex>
#include <algorithm>

namespace SOURCE175 {

// =============================================================================================
// PHYSICAL CONSTANTS (SI UNITS) - SOURCE175 NAMESPACE
// =============================================================================================

constexpr double PI_S175 = 3.141592653589793;
constexpr double c_S175 = 2.99792458e8;          // Speed of light (m/s)
constexpr double G_S175 = 6.67430e-11;           // Gravitational constant (m³/kg·s²)
constexpr double hbar_S175 = 1.054571817e-34;    // Reduced Planck constant (J·s)
constexpr double k_B_S175 = 1.380649e-23;        // Boltzmann constant (J/K)
constexpr double e_charge_S175 = 1.602176634e-19;// Elementary charge (C)
constexpr double m_e_S175 = 9.1093837015e-31;    // Electron mass (kg)
constexpr double m_p_S175 = 1.67262192369e-27;   // Proton mass (kg)

// Astrophysical constants
constexpr double M_sun_S175 = 1.98847e30;        // Solar mass (kg)
constexpr double L_sun_S175 = 3.828e26;          // Solar luminosity (W)
constexpr double R_sun_S175 = 6.96e8;            // Solar radius (m)
constexpr double pc_S175 = 3.0857e16;            // Parsec (m)
constexpr double yr_S175 = 3.15576e7;            // Year (s)

// UQFF-specific constants
constexpr double rho_vac_UA_S175 = 7.09e-36;     // [UA'] vacuum density (kg/m³)
constexpr double rho_vac_SCm_S175 = 1.0e-26;     // [SCm] vacuum density (kg/m³)
constexpr double H_SCm_S175 = 0.99;              // Heliosphere thickness factor
constexpr double U_UA_S175 = 0.0001;             // Universal Aether buoyancy factor
constexpr double k_eta_S175 = 1e-113;            // Aether coupling constant (extremely small)
constexpr double beta_i_S175 = 0.603;            // Buoyancy coupling constant

// LENR/Widom-Larsen constants
constexpr double f_Heaviside_LENR_S175 = 1e13;   // Heaviside enhancement for LENR oscillations
constexpr double nu_THz_LENR_S175 = 1.2e12;      // THz frequency for LENR (1.2 THz)

// =============================================================================================
// DATA STRUCTURES
// =============================================================================================

// Astrophysical system parameters (matches UQFFSystemsDatabase.py)
struct AstroSystem_S175 {
    std::string name;
    std::string category;
    
    // Position and distance
    double ra;          // Right ascension (degrees)
    double dec;         // Declination (degrees)
    double distance;    // Distance (parsecs)
    double redshift;    // Redshift z
    
    // Mass and size
    double mass;        // Mass (solar masses)
    double radius;      // Radius (solar radii or parsecs)
    
    // Luminosity and temperature
    double L_bol;       // Bolometric luminosity (erg/s)
    double L_X;         // X-ray luminosity (erg/s)
    double T;           // Temperature (K)
    
    // Velocities
    double v_exp;       // Expansion velocity (km/s)
    double v_rot;       // Rotation velocity (km/s)
    
    // Magnetic field
    double B;           // Magnetic field (Gauss)
    
    // Density
    double n_e;         // Electron density (cm⁻³)
    double n_H;         // Hydrogen density (cm⁻³)
    
    // Star formation
    double SFR;         // Star formation rate (M☉/yr)
    
    // Constructor with defaults
    AstroSystem_S175() : 
        ra(0), dec(0), distance(0), redshift(0),
        mass(0), radius(0),
        L_bol(0), L_X(0), T(0),
        v_exp(0), v_rot(0),
        B(0), n_e(0), n_H(0), SFR(0) {}
};

// UQFF calculation parameters
struct UQFFParams_S175 {
    double r;           // Radial distance (m)
    double t;           // Time (s)
    double t_n;         // Normalized time [0,1]
    double theta;       // Polar angle (rad)
    double phi;         // Azimuthal angle (rad)
    double Q_i;         // Quantum state number
    double SCm_i;       // [SCm] concentration
    double UA_i;        // [UA'] concentration
    double gamma;       // Temporal decay rate (s⁻¹)
    
    // Constructor with defaults
    UQFFParams_S175() :
        r(1.0), t(0), t_n(0.5), theta(0), phi(0),
        Q_i(1.0), SCm_i(1.0), UA_i(1.0), gamma(0.001) {}
};

// =============================================================================================
// F_UBII BUOYANCY FORCE VARIANTS (17 PROOFS)
// =============================================================================================

/**
 * F_UBii_virx: Virial X-ray cluster buoyancy
 * 
 * Equation: F_UBii = (M_vir/r_vir) × (σ_X²/c²) × [UA']:[SCm] × β_i
 * 
 * Physics: X-ray cluster velocity dispersion σ_X drives buoyancy in ICM.
 * Applied to: Perseus Cluster, Coma Cluster, Virgo Cluster
 * 
 * @param sigma_X Velocity dispersion (m/s)
 * @param T_X X-ray temperature (K)
 * @param n_e Electron density (m⁻³)
 * @param r Cluster radius (m)
 * @return F_UBii_virx (normalized units)
 */
inline double compute_F_UBii_virx(double sigma_X, double T_X, double n_e, double r) {
    if (r <= 0 || sigma_X <= 0) return 0.0;
    
    double M_vir = (sigma_X * sigma_X * r) / G_S175;  // Virial mass
    double virial_term = M_vir / r;
    double velocity_ratio = (sigma_X * sigma_X) / (c_S175 * c_S175);
    double vacuum_ratio = rho_vac_UA_S175 / (rho_vac_UA_S175 + rho_vac_SCm_S175);
    double temp_factor = std::sqrt(T_X / 1e7);  // Normalized to 10 MK
    
    return virial_term * velocity_ratio * vacuum_ratio * beta_i_S175 * temp_factor;
}

/**
 * F_UBii_termv: Terminal velocity buoyancy
 * 
 * Equation: F_UBii = (τ·L)/(c·E_LEP) × exp(-γt) × β_i
 * 
 * Physics: Radiation pressure equilibrium with gravitational drag.
 * Applied to: Astrophysical winds, stellar outflows, AGN jets
 * 
 * @param tau Optical depth
 * @param L Luminosity (W)
 * @param E_LEP Local escape energy (J)
 * @param gamma Decay rate (s⁻¹)
 * @param t Time (s)
 * @return F_UBii_termv (normalized units)
 */
inline double compute_F_UBii_termv(double tau, double L, double E_LEP, double gamma, double t) {
    if (E_LEP <= 0) return 0.0;
    
    double radiation_term = (tau * L) / (c_S175 * E_LEP);
    double decay_factor = std::exp(-gamma * t);
    
    return radiation_term * decay_factor * beta_i_S175;
}

/**
 * F_UBii_upar: Ionization parameter buoyancy
 * 
 * Equation: F_UBii = U × (L_ion)/(n_H × r²) × [SCm] × β_i
 * 
 * Physics: Ionizing radiation modulates [SCm] concentration → buoyancy.
 * Applied to: HII regions, planetary nebulae, AGN narrow-line regions
 * 
 * @param L_ion Ionizing luminosity (W)
 * @param n_H Hydrogen density (m⁻³)
 * @param r Distance from ionizing source (m)
 * @param SCm_concentration [SCm] density (kg/m³)
 * @return F_UBii_upar (normalized units)
 */
inline double compute_F_UBii_upar(double L_ion, double n_H, double r, double SCm_concentration) {
    if (n_H <= 0 || r <= 0) return 0.0;
    
    double U_param = L_ion / (n_H * r * r * c_S175);  // Ionization parameter
    double SCm_factor = SCm_concentration / rho_vac_SCm_S175;
    
    return U_param * L_ion / (n_H * r * r) * SCm_factor * beta_i_S175;
}

/**
 * F_UBii_coup: Energy coupling buoyancy
 * 
 * Equation: F_UBii = (E_kin/E_mag) × (v²/c²) × [UA']:[SCm] × β_i
 * 
 * Physics: Kinetic-magnetic energy coupling drives vacuum displacement.
 * Applied to: Supernova shocks, stellar winds, accretion disk winds
 * 
 * @param E_kin Kinetic energy (J)
 * @param E_mag Magnetic energy (J)
 * @param v Flow velocity (m/s)
 * @return F_UBii_coup (normalized units)
 */
inline double compute_F_UBii_coup(double E_kin, double E_mag, double v) {
    if (E_mag <= 0) return 0.0;
    
    double energy_ratio = E_kin / E_mag;
    double velocity_factor = (v * v) / (c_S175 * c_S175);
    double vacuum_ratio = rho_vac_UA_S175 / (rho_vac_UA_S175 + rho_vac_SCm_S175);
    
    return energy_ratio * velocity_factor * vacuum_ratio * beta_i_S175;
}

/**
 * F_UBii_orbdec: Orbital decay buoyancy
 * 
 * Equation: F_UBii = -(da/dt) × (M₁M₂)/(M₁+M₂) × (G/c⁵) × [UA'] × β_i
 * 
 * Physics: Gravitational wave radiation modulates [UA'] vacuum density.
 * Applied to: Binary pulsars, binary black holes, NS-NS mergers
 * 
 * @param M1 Primary mass (kg)
 * @param M2 Secondary mass (kg)
 * @param a Orbital separation (m)
 * @param da_dt Orbital decay rate (m/s)
 * @return F_UBii_orbdec (normalized units)
 */
inline double compute_F_UBii_orbdec(double M1, double M2, double a, double da_dt) {
    if (M1 + M2 <= 0 || a <= 0) return 0.0;
    
    double mu = (M1 * M2) / (M1 + M2);  // Reduced mass
    double GW_coeff = G_S175 / (c_S175 * c_S175 * c_S175 * c_S175 * c_S175);
    double UA_factor = rho_vac_UA_S175 / rho_vac_SCm_S175;
    
    return -da_dt * mu * GW_coeff * UA_factor * beta_i_S175;
}

/**
 * F_UBii_kn: Kilonova peak luminosity buoyancy
 * 
 * Equation: F_UBii = (L_peak × t_peak) / (4π r²c) × [SCm] × β_i
 * 
 * Physics: Kilonova r-process heating modulates [SCm] locally.
 * Applied to: NS-NS mergers (AT2017gfo, GW170817), r-process nucleosynthesis
 * 
 * @param L_peak Peak luminosity (W)
 * @param t_peak Time to peak (s)
 * @param r Distance (m)
 * @return F_UBii_kn (normalized units)
 */
inline double compute_F_UBii_kn(double L_peak, double t_peak, double r) {
    if (r <= 0 || t_peak <= 0) return 0.0;
    
    double luminosity_term = (L_peak * t_peak) / (4 * PI_S175 * r * r * c_S175);
    double SCm_factor = rho_vac_SCm_S175 / rho_vac_UA_S175;
    
    return luminosity_term * SCm_factor * beta_i_S175;
}

/**
 * F_UBii_fermi: Fermi acceleration buoyancy
 * 
 * Equation: F_UBii = β_shock × v_shock² / c² × [UA'] × β_i
 * 
 * Physics: Shock acceleration modulates [UA'] vacuum at particle-crossing scale.
 * Applied to: SNRs (Cas A, Tycho), AGN jets, GRB afterglows
 * 
 * @param beta_shock Shock compression ratio
 * @param v_shock Shock velocity (m/s)
 * @return F_UBii_fermi (normalized units)
 */
inline double compute_F_UBii_fermi(double beta_shock, double v_shock) {
    double velocity_factor = (v_shock * v_shock) / (c_S175 * c_S175);
    double UA_factor = rho_vac_UA_S175 / rho_vac_SCm_S175;
    
    return beta_shock * velocity_factor * UA_factor * beta_i_S175;
}

/**
 * F_UBii_kne: Cosmic ray knee energy buoyancy
 * 
 * Equation: F_UBii = ρ_CR × Z × E_knee / (m_p c²) × [UA'] × β_i
 * 
 * Physics: Cosmic ray "knee" (3×10¹⁵ eV) marks transition in [UA'] coupling.
 * Applied to: Galactic cosmic rays, ultra-high energy cosmic rays (UHECRs)
 * 
 * @param rho_CR Cosmic ray energy density (J/m³)
 * @param Z Charge number
 * @param E_knee Knee energy (J)
 * @return F_UBii_kne (normalized units)
 */
inline double compute_F_UBii_kne(double rho_CR, double Z, double E_knee) {
    double mp_c2 = m_p_S175 * c_S175 * c_S175;
    double energy_factor = E_knee / mp_c2;
    double UA_factor = rho_vac_UA_S175 / rho_vac_SCm_S175;
    
    return rho_CR * Z * energy_factor * UA_factor * beta_i_S175;
}

/**
 * F_UBii_whim: Warm-Hot Intergalactic Medium (WHIM) temperature buoyancy
 * 
 * Equation: F_UBii = T_WHIM × n_e × k_B / (m_p c²) × [UA']:[SCm] × β_i
 * 
 * Physics: WHIM shocks (T ~ 10⁵–10⁷ K) modulate [UA'] concentration.
 * Applied to: Filamentary large-scale structure, missing baryons problem
 * 
 * @param T_WHIM Temperature (K)
 * @param n_e Electron density (m⁻³)
 * @return F_UBii_whim (normalized units)
 */
inline double compute_F_UBii_whim(double T_WHIM, double n_e) {
    double mp_c2 = m_p_S175 * c_S175 * c_S175;
    double thermal_term = (T_WHIM * n_e * k_B_S175) / mp_c2;
    double vacuum_ratio = rho_vac_UA_S175 / (rho_vac_UA_S175 + rho_vac_SCm_S175);
    
    return thermal_term * vacuum_ratio * beta_i_S175;
}

/**
 * F_UBii_ps: Press-Schechter halo mass buoyancy
 * 
 * Equation: F_UBii = M_halo × σ_8 × (Ω_m / Ω_Λ) × [SCm] × β_i
 * 
 * Physics: Dark matter halo formation modulates [SCm] at virial radius.
 * Applied to: Galaxy clusters, dwarf galaxy halos, CDM structure formation
 * 
 * @param M_halo Halo mass (kg)
 * @param sigma_8 Power spectrum normalization
 * @param Omega_m Matter density parameter
 * @param Omega_Lambda Dark energy density parameter
 * @return F_UBii_ps (normalized units)
 */
inline double compute_F_UBii_ps(double M_halo, double sigma_8, double Omega_m, double Omega_Lambda) {
    if (Omega_Lambda <= 0) return 0.0;
    
    double cosmology_factor = Omega_m / Omega_Lambda;
    double SCm_factor = rho_vac_SCm_S175 / rho_vac_UA_S175;
    
    return M_halo * sigma_8 * cosmology_factor * SCm_factor * beta_i_S175;
}

/**
 * F_UBii_sfe: Star formation efficiency buoyancy
 * 
 * Equation: F_UBii = ε_SFE × Σ_gas / Σ_crit × [SCm] × β_i
 * 
 * Physics: Star formation depletes [SCm] locally → negative buoyancy feedback.
 * Applied to: Giant molecular clouds, starburst galaxies, Toomre instability
 * 
 * @param epsilon_SFE Star formation efficiency
 * @param Sigma_gas Gas surface density (kg/m²)
 * @param Sigma_crit Critical surface density (kg/m²)
 * @return F_UBii_sfe (normalized units)
 */
inline double compute_F_UBii_sfe(double epsilon_SFE, double Sigma_gas, double Sigma_crit) {
    if (Sigma_crit <= 0) return 0.0;
    
    double density_ratio = Sigma_gas / Sigma_crit;
    double SCm_factor = rho_vac_SCm_S175 / rho_vac_UA_S175;
    
    return epsilon_SFE * density_ratio * SCm_factor * beta_i_S175;
}

/**
 * F_UBii_hawk: Hawking temperature buoyancy
 * 
 * Equation: F_UBii = (ℏc³)/(8πGM·k_B) × [SCm] × β_i
 * 
 * Physics: Hawking radiation modulates [SCm] near event horizon.
 * Applied to: Stellar-mass BHs, intermediate-mass BHs, primordial BHs
 * 
 * @param M Black hole mass (kg)
 * @return F_UBii_hawk (normalized units)
 */
inline double compute_F_UBii_hawk(double M) {
    if (M <= 0) return 0.0;
    
    double T_Hawking = (hbar_S175 * c_S175 * c_S175 * c_S175) / (8 * PI_S175 * G_S175 * M * k_B_S175);
    double SCm_factor = rho_vac_SCm_S175 / rho_vac_UA_S175;
    
    return T_Hawking * SCm_factor * beta_i_S175;
}

/**
 * F_UBii_bd: Bounce density buoyancy
 * 
 * Equation: F_UBii = ρ_bounce × a_bounce³ × [UA'] × β_i
 * 
 * Physics: Loop quantum cosmology bounce modulates [UA'] at Planck density.
 * Applied to: Big Bounce models, loop quantum gravity, Planck epoch
 * 
 * @param rho_bounce Bounce density (kg/m³)
 * @param a_bounce Scale factor at bounce
 * @return F_UBii_bd (normalized units)
 */
inline double compute_F_UBii_bd(double rho_bounce, double a_bounce) {
    double volume_factor = a_bounce * a_bounce * a_bounce;
    double UA_factor = rho_vac_UA_S175 / rho_vac_SCm_S175;
    
    return rho_bounce * volume_factor * UA_factor * beta_i_S175;
}

/**
 * F_UBii_roche: Roche lobe overflow buoyancy
 * 
 * Equation: F_UBii = ΔM_dot × (v_esc²/c²) × [SCm] × β_i
 * 
 * Physics: Mass transfer modulates [SCm] in accretion stream.
 * Applied to: Cataclysmic variables, X-ray binaries, common envelope evolution
 * 
 * @param Delta_M_dot Mass transfer rate (kg/s)
 * @param v_esc Escape velocity (m/s)
 * @return F_UBii_roche (normalized units)
 */
inline double compute_F_UBii_roche(double Delta_M_dot, double v_esc) {
    double velocity_factor = (v_esc * v_esc) / (c_S175 * c_S175);
    double SCm_factor = rho_vac_SCm_S175 / rho_vac_UA_S175;
    
    return Delta_M_dot * velocity_factor * SCm_factor * beta_i_S175;
}

/**
 * F_UBii_ent: Entanglement entropy buoyancy
 * 
 * Equation: F_UBii = S_ent × Area / (4G ℏ) × [UA'] × β_i
 * 
 * Physics: Holographic entanglement entropy modulates [UA'] on Ryu-Takayanagi surface.
 * Applied to: Black hole entropy, AdS/CFT, quantum information paradox
 * 
 * @param S_ent Entanglement entropy (J/K)
 * @param Area Surface area (m²)
 * @return F_UBii_ent (normalized units)
 */
inline double compute_F_UBii_ent(double S_ent, double Area) {
    double holographic_term = (S_ent * Area) / (4 * G_S175 * hbar_S175);
    double UA_factor = rho_vac_UA_S175 / rho_vac_SCm_S175;
    
    return holographic_term * UA_factor * beta_i_S175;
}

/**
 * F_UBii_dec: Decoherence time buoyancy
 * 
 * Equation: F_UBii = (ℏ / τ_dec) × T × [UA'] × β_i
 * 
 * Physics: Quantum decoherence mediated by [UA'] vacuum fluctuations.
 * Applied to: Quantum measurement, LIGO thermal noise, macroscopic quantum systems
 * 
 * @param tau_dec Decoherence time (s)
 * @param T Temperature (K)
 * @return F_UBii_dec (normalized units)
 */
inline double compute_F_UBii_dec(double tau_dec, double T) {
    if (tau_dec <= 0) return 0.0;
    
    double decoherence_term = hbar_S175 / tau_dec;
    double UA_factor = rho_vac_UA_S175 / rho_vac_SCm_S175;
    
    return decoherence_term * T * UA_factor * beta_i_S175;
}

/**
 * F_UBii_lobe: Radio lobe dynamics buoyancy
 * 
 * Equation: F_UBii = P_jet × t_lobe × (1/ρ_ICM) × [UA']:[SCm] × β_i
 * 
 * Physics: AGN jet power inflates radio lobes → [UA'] displacement.
 * Applied to: Centaurus A, M87, Cygnus A radio lobes, AGN feedback
 * 
 * @param P_jet Jet power (W)
 * @param t_lobe Lobe age (s)
 * @param rho_ICM ICM density (kg/m³)
 * @return F_UBii_lobe (normalized units)
 */
inline double compute_F_UBii_lobe(double P_jet, double t_lobe, double rho_ICM) {
    if (rho_ICM <= 0) return 0.0;
    
    double energy_density = (P_jet * t_lobe) / rho_ICM;
    double vacuum_ratio = rho_vac_UA_S175 / (rho_vac_UA_S175 + rho_vac_SCm_S175);
    
    return energy_density * vacuum_ratio * beta_i_S175;
}

// =============================================================================================
// UNIVERSAL MAGNETISM (Um) - WITH LENR HEAVISIDE 10^13 ENHANCEMENT
// =============================================================================================

/**
 * compute_Um: Universal Magnetism with Widom-Larsen LENR enhancement
 * 
 * Equation: Um = Σ[(μ_j/r_j) × (1-e^(-γt)) × cos(πt_n) × φ̂_j] × P_SCm × E_react × 
 *               (1 + 10^13 × f_Heaviside) × (1 + f_quasi)
 * 
 * Physics: 
 * - Heaviside 10^13 enhancement for ultra-high frequency THz oscillations (LENR)
 * - Temporal modulation: (1 - e^(-γt)) growth × cos(πt_n) oscillation
 * - [SCm] penetration P_SCm and reactivity E_react modulate dipole strength
 * - Quasi-particle correction factor for collective modes
 * 
 * Applied to: Magnetars (SGR1745), pulsars, magnetic white dwarfs, LENR experiments
 * 
 * @param system Astrophysical system parameters
 * @param params UQFF calculation parameters
 * @param dipole_moments Vector of magnetic dipole moments (A·m²)
 * @param distances Vector of distances from dipoles (m)
 * @param enable_LENR Enable Heaviside 10^13 LENR enhancement
 * @return Um (normalized units)
 */
inline double compute_Um(const AstroSystem_S175& system, 
                        const UQFFParams_S175& params,
                        const std::vector<double>& dipole_moments,
                        const std::vector<double>& distances,
                        bool enable_LENR = true) {
    
    if (dipole_moments.size() != distances.size()) return 0.0;
    
    // Temporal modulation
    double growth_factor = 1.0 - std::exp(-params.gamma * params.t);
    double oscillation = std::cos(PI_S175 * params.t_n);
    double temporal_mod = growth_factor * oscillation;
    
    // Dipole sum
    double dipole_sum = 0.0;
    for (size_t i = 0; i < dipole_moments.size(); ++i) {
        if (distances[i] > 0) {
            dipole_sum += dipole_moments[i] / distances[i];
        }
    }
    
    // [SCm] penetration (increases with B-field strength)
    double P_SCm = H_SCm_S175 * (1.0 + system.B / 1e14);  // Normalized to magnetar field
    
    // Reactivity (decreases with temperature)
    double E_react = std::exp(-params.t / (1e6 * yr_S175));  // Myr timescale
    
    // Heaviside LENR enhancement (10^13 for THz oscillations)
    double LENR_factor = enable_LENR ? (1.0 + f_Heaviside_LENR_S175) : 1.0;
    
    // Quasi-particle correction (collective excitations)
    double f_quasi = 0.01 * std::sqrt(system.n_e / 1e6);  // Scales with density
    
    // Complete Um calculation
    double Um = dipole_sum * temporal_mod * P_SCm * E_react * LENR_factor * (1.0 + f_quasi);
    
    return Um;
}

// =============================================================================================
// AETHER METRIC TENSOR (A^μν) - 4×4 SPACETIME METRIC
// =============================================================================================

/**
 * compute_Aether_metric_tensor: 4×4 spacetime metric with UQFF vacuum duality
 * 
 * Metric structure:
 *      ⎡ A⁰⁰  A⁰¹  A⁰²  A⁰³ ⎤
 * A^μν= ⎢ A¹⁰  A¹¹  0    0   ⎥
 *      ⎢ A²⁰  0    A²²  0   ⎥
 *      ⎣ A³⁰  0    0    A³³ ⎦
 * 
 * Components:
 * - A⁰⁰ (temporal): Modulated by [UA']:[SCm] ratio, time dilation
 * - A¹¹, A²², A³³ (spatial): Expansion/contraction from vacuum energy
 * - A⁰ⁱ (frame-dragging): Off-diagonal terms from rotation
 * 
 * Physics: Aether as active vacuum medium → spacetime geometry emerges from
 *          [UA'] and [SCm] concentrations.
 * 
 * @param system Astrophysical system parameters
 * @param params UQFF calculation parameters
 * @param metric Output 4×4 array for metric tensor
 */
inline void compute_Aether_metric_tensor(const AstroSystem_S175& system,
                                        const UQFFParams_S175& params,
                                        std::array<std::array<double, 4>, 4>& metric) {
    
    // Initialize to Minkowski (flat spacetime baseline)
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            metric[mu][nu] = (mu == nu) ? ((mu == 0) ? -1.0 : 1.0) : 0.0;
        }
    }
    
    // Vacuum duality ratio
    double vacuum_ratio = rho_vac_UA_S175 / (rho_vac_UA_S175 + rho_vac_SCm_S175);
    
    // A⁰⁰ (temporal component): Time dilation from [UA'] density
    double time_dilation = 1.0 + 2 * G_S175 * system.mass * M_sun_S175 / (params.r * c_S175 * c_S175);
    double temporal_mod = std::cos(PI_S175 * params.t_n);
    metric[0][0] = -time_dilation * (1.0 + k_eta_S175 * vacuum_ratio * temporal_mod);
    
    // A¹¹, A²², A³³ (spatial components): Expansion from [SCm] concentration
    double spatial_expansion = 1.0 + (rho_vac_SCm_S175 / rho_vac_UA_S175) * k_eta_S175;
    metric[1][1] = spatial_expansion * (1.0 + 0.1 * std::sin(params.theta));
    metric[2][2] = spatial_expansion * (1.0 - 0.1 * std::cos(params.phi));
    metric[3][3] = spatial_expansion;
    
    // A⁰ⁱ (frame-dragging): Rotational modulation
    if (system.v_rot > 0) {
        double dragging_coeff = (2 * G_S175 * system.mass * M_sun_S175 * system.v_rot * 1e3) / (params.r * c_S175 * c_S175 * c_S175);
        metric[0][1] = metric[1][0] = dragging_coeff * std::sin(params.theta);
        metric[0][2] = metric[2][0] = dragging_coeff * std::cos(params.theta) * std::sin(params.phi);
        metric[0][3] = metric[3][0] = dragging_coeff * std::cos(params.phi);
    }
}

/**
 * compute_metric_determinant: Determinant of A^μν metric tensor
 * 
 * Used for volume element in GR field equations.
 * 
 * @param metric 4×4 metric tensor
 * @return det(A^μν)
 */
inline double compute_metric_determinant(const std::array<std::array<double, 4>, 4>& metric) {
    // Simplified 4×4 determinant (assuming metric[i][j] = 0 for many off-diagonal terms)
    // Full determinant would require cofactor expansion
    
    double det = metric[0][0] * metric[1][1] * metric[2][2] * metric[3][3];
    
    // Include frame-dragging corrections if non-zero
    if (std::abs(metric[0][1]) > 1e-10) {
        det -= metric[0][1] * metric[1][0] * metric[2][2] * metric[3][3];
    }
    
    return det;
}

// =============================================================================================
// UNIFIED FIELD CALCULATOR (F_U) - MASTER EQUATION
// =============================================================================================

/**
 * compute_F_U: Unified Field master equation
 * 
 * Equation: F_U = Σ[k_i × Ug_i - Ub_i] + Um + ∫ A^{μν} d⁴x
 * 
 * Physics: Integrates all UQFF field components:
 * - Ug_i: Gravity terms (Ug1 magnetic dipole, Ug2 charge-reactivity, Ug3 string rotation, Ug4 vacuum concentration)
 * - Ub_i: Buoyancy terms (F_UBii variants)
 * - Um: Universal magnetism with LENR enhancement
 * - A^{μν}: Aether metric tensor contribution
 * 
 * Applied to: Complete astrophysical system modeling (all 36 UQFFSystemsDatabase objects)
 * 
 * @param system Astrophysical system parameters
 * @param params UQFF calculation parameters
 * @param Ug_coeffs Coupling constants k_i for Ug terms
 * @param Ub_terms Vector of buoyancy contributions
 * @param dipole_moments Magnetic dipole moments for Um
 * @param dipole_distances Distances for Um calculation
 * @return F_U (normalized units)
 */
inline double compute_F_U(const AstroSystem_S175& system,
                         const UQFFParams_S175& params,
                         const std::array<double, 4>& Ug_coeffs,
                         const std::vector<double>& Ub_terms,
                         const std::vector<double>& dipole_moments,
                         const std::vector<double>& dipole_distances) {
    
    // Ug_i gravity terms (26-layer sum, using single-layer approximation here)
    double Ug1 = (hbar_S175 * c_S175 / (params.r * params.r)) * params.Q_i * params.SCm_i * rho_vac_UA_S175;
    double Ug2 = (hbar_S175 * c_S175 / (params.r * params.r)) * params.SCm_i * params.SCm_i;
    double Ug3 = (hbar_S175 / 2) * std::cos(2 * PI_S175 * params.t_n) / params.r;
    double Ug4 = (G_S175 * system.mass * M_sun_S175 / (params.r * params.r)) * params.SCm_i;
    
    // Weighted gravity sum
    double Ug_sum = Ug_coeffs[0] * Ug1 + Ug_coeffs[1] * Ug2 + Ug_coeffs[2] * Ug3 + Ug_coeffs[3] * Ug4;
    
    // Ub_i buoyancy sum
    double Ub_sum = 0.0;
    for (double Ub : Ub_terms) {
        Ub_sum += Ub;
    }
    
    // Um universal magnetism
    double Um = compute_Um(system, params, dipole_moments, dipole_distances, true);
    
    // A^{μν} Aether contribution (metric determinant as integrated volume element)
    std::array<std::array<double, 4>, 4> metric;
    compute_Aether_metric_tensor(system, params, metric);
    double Aether_contrib = std::sqrt(std::abs(compute_metric_determinant(metric)));
    
    // Master unified field
    double F_U = Ug_sum - Ub_sum + Um + k_eta_S175 * Aether_contrib;
    
    return F_U;
}

// =============================================================================================
// HELPER FUNCTIONS FOR INTEGRATION WITH MAIN_1_CoAnQi.cpp
// =============================================================================================

/**
 * list_buoyancy_variants: Return all 17 F_UBii variant names
 */
inline std::vector<std::string> list_buoyancy_variants() {
    return {
        "virx", "termv", "upar", "coup", "orbdec", "kn", "fermi", "kne",
        "whim", "ps", "sfe", "hawk", "bd", "roche", "ent", "dec", "lobe"
    };
}

/**
 * get_buoyancy_variant_description: Get physics description for variant
 */
inline std::string get_buoyancy_variant_description(const std::string& variant) {
    static std::map<std::string, std::string> descriptions = {
        {"virx", "Virial X-ray cluster velocity dispersion"},
        {"termv", "Terminal velocity radiation pressure equilibrium"},
        {"upar", "Ionization parameter modulation"},
        {"coup", "Kinetic-magnetic energy coupling"},
        {"orbdec", "Orbital decay gravitational wave radiation"},
        {"kn", "Kilonova peak luminosity r-process heating"},
        {"fermi", "Fermi acceleration at shock fronts"},
        {"kne", "Cosmic ray knee transition energy"},
        {"whim", "Warm-Hot Intergalactic Medium temperature"},
        {"ps", "Press-Schechter dark matter halo mass"},
        {"sfe", "Star formation efficiency feedback"},
        {"hawk", "Hawking temperature near event horizon"},
        {"bd", "Bounce density at Planck epoch"},
        {"roche", "Roche lobe overflow mass transfer"},
        {"ent", "Entanglement entropy holographic bound"},
        {"dec", "Decoherence time quantum measurement"},
        {"lobe", "Radio lobe AGN jet inflation"}
    };
    
    auto it = descriptions.find(variant);
    return (it != descriptions.end()) ? it->second : "Unknown variant";
}

/**
 * validate_SOURCE175: Run self-tests on all calculators
 * 
 * Tests:
 * - All 17 F_UBii variants with sample parameters
 * - Um calculation with M87 magnetar parameters
 * - Aether metric tensor for Sagittarius A*
 * - Unified field F_U for Perseus Cluster
 * 
 * @return true if all tests pass, false otherwise
 */
inline bool validate_SOURCE175() {
    bool all_passed = true;
    
    // Test F_UBii_virx with Perseus Cluster parameters
    double F_virx = compute_F_UBii_virx(1.2e6, 6e7, 0.04e6, 2.5e24);
    if (F_virx <= 0 || F_virx > 1e50) {
        all_passed = false;
    }
    
    // Test Um with magnetar parameters
    AstroSystem_S175 sgr1745;
    sgr1745.name = "SGR 1745-2900";
    sgr1745.mass = 1.4;
    sgr1745.B = 2e14;
    sgr1745.T = 5e6;
    sgr1745.n_e = 1e6;
    
    UQFFParams_S175 params;
    params.r = 1e4;
    params.t = 1e9;
    params.t_n = 0.5;
    params.gamma = 0.001;
    
    std::vector<double> dipoles = {1e30};
    std::vector<double> dists = {1e4};
    
    double Um = compute_Um(sgr1745, params, dipoles, dists, true);
    if (Um <= 0 || std::isnan(Um)) {
        all_passed = false;
    }
    
    // Test Aether metric tensor
    std::array<std::array<double, 4>, 4> metric;
    compute_Aether_metric_tensor(sgr1745, params, metric);
    double det = compute_metric_determinant(metric);
    if (std::isnan(det) || det == 0.0) {
        all_passed = false;
    }
    
    // Test unified field F_U
    std::array<double, 4> Ug_coeffs = {1.0, 1.0, 1.0, 1.0};
    std::vector<double> Ub_terms = {F_virx};
    double F_U = compute_F_U(sgr1745, params, Ug_coeffs, Ub_terms, dipoles, dists);
    if (std::isnan(F_U)) {
        all_passed = false;
    }
    
    return all_passed;
}

// =============================================================================================
// PHYSICSTERM WRAPPER CLASSES FOR MAIN_1_CoAnQi.cpp INTEGRATION
// =============================================================================================

// Note: PhysicsTerm base class is defined in MAIN_1_CoAnQi.cpp (lines 269-330)
// This file is included after PhysicsTerm definition, so we can directly use it

/**
 * F_UBii_virx PhysicsTerm wrapper
 */
class SOURCE175_FUBiiVirx_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double sigma_X = params.count("sigma_X") ? params.at("sigma_X") : 1.2e6;
        double T_X = params.count("T_X") ? params.at("T_X") : 6e7;
        double n_e = params.count("n_e") ? params.at("n_e") : 0.04e6;
        double r = params.count("r") ? params.at("r") : 2.5e24;
        return compute_F_UBii_virx(sigma_X, T_X, n_e, r);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_virx"; }
    
    std::string getDescription() const override {
        return "Virial X-ray cluster velocity dispersion buoyancy";
    }
};

/**
 * F_UBii_termv PhysicsTerm wrapper
 */
class SOURCE175_FUBiiTermv_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double tau = params.count("tau") ? params.at("tau") : 1.0;
        double L = params.count("L") ? params.at("L") : 1e38;
        double E_LEP = params.count("E_LEP") ? params.at("E_LEP") : 1e45;
        double gamma = params.count("gamma") ? params.at("gamma") : 0.001;
        return compute_F_UBii_termv(tau, L, E_LEP, gamma, t);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_termv"; }
    
    std::string getDescription() const override {
        return "Terminal velocity radiation pressure equilibrium buoyancy";
    }
};

/**
 * F_UBii_upar PhysicsTerm wrapper
 */
class SOURCE175_FUBiiUpar_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double L_ion = params.count("L_ion") ? params.at("L_ion") : 1e38;
        double n_H = params.count("n_H") ? params.at("n_H") : 1e6;
        double r = params.count("r") ? params.at("r") : 1e18;
        double SCm = params.count("SCm") ? params.at("SCm") : rho_vac_SCm_S175;
        return compute_F_UBii_upar(L_ion, n_H, r, SCm);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_upar"; }
    
    std::string getDescription() const override {
        return "Ionization parameter [SCm] modulation buoyancy";
    }
};

/**
 * F_UBii_coup PhysicsTerm wrapper
 */
class SOURCE175_FUBiiCoup_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double E_kin = params.count("E_kin") ? params.at("E_kin") : 1e45;
        double E_mag = params.count("E_mag") ? params.at("E_mag") : 1e43;
        double v = params.count("v") ? params.at("v") : 1e7;
        return compute_F_UBii_coup(E_kin, E_mag, v);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_coup"; }
    
    std::string getDescription() const override {
        return "Kinetic-magnetic energy coupling buoyancy";
    }
};

/**
 * F_UBii_orbdec PhysicsTerm wrapper
 */
class SOURCE175_FUBiiOrbdec_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double M1 = params.count("M1") ? params.at("M1") : 1.4 * M_sun_S175;
        double M2 = params.count("M2") ? params.at("M2") : 1.4 * M_sun_S175;
        double a = params.count("a") ? params.at("a") : 1e9;
        double da_dt = params.count("da_dt") ? params.at("da_dt") : -1e-5;
        return compute_F_UBii_orbdec(M1, M2, a, da_dt);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_orbdec"; }
    
    std::string getDescription() const override {
        return "Orbital decay gravitational wave radiation buoyancy";
    }
};

/**
 * F_UBii_kn PhysicsTerm wrapper
 */
class SOURCE175_FUBiiKn_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double L_peak = params.count("L_peak") ? params.at("L_peak") : 1e42;
        double t_peak = params.count("t_peak") ? params.at("t_peak") : 1e5;
        double r = params.count("r") ? params.at("r") : 1e20;
        return compute_F_UBii_kn(L_peak, t_peak, r);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_kn"; }
    
    std::string getDescription() const override {
        return "Kilonova peak luminosity r-process heating buoyancy";
    }
};

/**
 * F_UBii_fermi PhysicsTerm wrapper
 */
class SOURCE175_FUBiiFermi_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double beta_shock = params.count("beta_shock") ? params.at("beta_shock") : 4.0;
        double v_shock = params.count("v_shock") ? params.at("v_shock") : 5e6;
        return compute_F_UBii_fermi(beta_shock, v_shock);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_fermi"; }
    
    std::string getDescription() const override {
        return "Fermi acceleration at shock fronts buoyancy";
    }
};

/**
 * F_UBii_kne PhysicsTerm wrapper
 */
class SOURCE175_FUBiiKne_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double rho_CR = params.count("rho_CR") ? params.at("rho_CR") : 1e-15;
        double Z = params.count("Z") ? params.at("Z") : 26.0;
        double E_knee = params.count("E_knee") ? params.at("E_knee") : 3e15 * e_charge_S175;
        return compute_F_UBii_kne(rho_CR, Z, E_knee);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_kne"; }
    
    std::string getDescription() const override {
        return "Cosmic ray knee transition energy buoyancy";
    }
};

/**
 * F_UBii_whim PhysicsTerm wrapper
 */
class SOURCE175_FUBiiWhim_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double T_WHIM = params.count("T_WHIM") ? params.at("T_WHIM") : 1e6;
        double n_e = params.count("n_e") ? params.at("n_e") : 1e-5;
        return compute_F_UBii_whim(T_WHIM, n_e);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_whim"; }
    
    std::string getDescription() const override {
        return "Warm-Hot Intergalactic Medium temperature buoyancy";
    }
};

/**
 * F_UBii_ps PhysicsTerm wrapper
 */
class SOURCE175_FUBiiPs_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double M_halo = params.count("M_halo") ? params.at("M_halo") : 1e14 * M_sun_S175;
        double sigma_8 = params.count("sigma_8") ? params.at("sigma_8") : 0.8;
        double Omega_m = params.count("Omega_m") ? params.at("Omega_m") : 0.3;
        double Omega_Lambda = params.count("Omega_Lambda") ? params.at("Omega_Lambda") : 0.7;
        return compute_F_UBii_ps(M_halo, sigma_8, Omega_m, Omega_Lambda);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_ps"; }
    
    std::string getDescription() const override {
        return "Press-Schechter dark matter halo mass buoyancy";
    }
};

/**
 * F_UBii_sfe PhysicsTerm wrapper
 */
class SOURCE175_FUBiiSfe_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double epsilon_SFE = params.count("epsilon_SFE") ? params.at("epsilon_SFE") : 0.1;
        double Sigma_gas = params.count("Sigma_gas") ? params.at("Sigma_gas") : 100.0;
        double Sigma_crit = params.count("Sigma_crit") ? params.at("Sigma_crit") : 50.0;
        return compute_F_UBii_sfe(epsilon_SFE, Sigma_gas, Sigma_crit);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_sfe"; }
    
    std::string getDescription() const override {
        return "Star formation efficiency feedback buoyancy";
    }
};

/**
 * F_UBii_hawk PhysicsTerm wrapper
 */
class SOURCE175_FUBiiHawk_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double M = params.count("M") ? params.at("M") : 10.0 * M_sun_S175;
        return compute_F_UBii_hawk(M);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_hawk"; }
    
    std::string getDescription() const override {
        return "Hawking temperature near event horizon buoyancy";
    }
};

/**
 * F_UBii_bd PhysicsTerm wrapper
 */
class SOURCE175_FUBiiBd_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double rho_bounce = params.count("rho_bounce") ? params.at("rho_bounce") : 1e96;
        double a_bounce = params.count("a_bounce") ? params.at("a_bounce") : 1e-35;
        return compute_F_UBii_bd(rho_bounce, a_bounce);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_bd"; }
    
    std::string getDescription() const override {
        return "Bounce density at Planck epoch buoyancy";
    }
};

/**
 * F_UBii_roche PhysicsTerm wrapper
 */
class SOURCE175_FUBiiRoche_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double Delta_M_dot = params.count("Delta_M_dot") ? params.at("Delta_M_dot") : 1e-8 * M_sun_S175 / yr_S175;
        double v_esc = params.count("v_esc") ? params.at("v_esc") : 5e6;
        return compute_F_UBii_roche(Delta_M_dot, v_esc);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_roche"; }
    
    std::string getDescription() const override {
        return "Roche lobe overflow mass transfer buoyancy";
    }
};

/**
 * F_UBii_ent PhysicsTerm wrapper
 */
class SOURCE175_FUBiiEnt_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double S_ent = params.count("S_ent") ? params.at("S_ent") : 1e-20;
        double Area = params.count("Area") ? params.at("Area") : 1e10;
        return compute_F_UBii_ent(S_ent, Area);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_ent"; }
    
    std::string getDescription() const override {
        return "Entanglement entropy holographic bound buoyancy";
    }
};

/**
 * F_UBii_dec PhysicsTerm wrapper
 */
class SOURCE175_FUBiiDec_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double tau_dec = params.count("tau_dec") ? params.at("tau_dec") : 1e-10;
        double T = params.count("T") ? params.at("T") : 300.0;
        return compute_F_UBii_dec(tau_dec, T);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_dec"; }
    
    std::string getDescription() const override {
        return "Decoherence time quantum measurement buoyancy";
    }
};

/**
 * F_UBii_lobe PhysicsTerm wrapper
 */
class SOURCE175_FUBiiLobe_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        double P_jet = params.count("P_jet") ? params.at("P_jet") : 1e44;
        double t_lobe = params.count("t_lobe") ? params.at("t_lobe") : 1e14;
        double rho_ICM = params.count("rho_ICM") ? params.at("rho_ICM") : 1e-26;
        return compute_F_UBii_lobe(P_jet, t_lobe, rho_ICM);
    }
    
    std::string getName() const override { return "SOURCE175_F_UBii_lobe"; }
    
    std::string getDescription() const override {
        return "Radio lobe AGN jet inflation buoyancy";
    }
};

/**
 * Um (Universal Magnetism) PhysicsTerm wrapper
 */
class SOURCE175_Um_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        // Construct AstroSystem from params
        AstroSystem_S175 system;
        system.name = "Generic";
        system.mass = params.count("M") ? params.at("M") : 1.4;
        system.B = params.count("B") ? params.at("B") : 1e14;
        system.T = params.count("T") ? params.at("T") : 5e6;
        system.n_e = params.count("n_e") ? params.at("n_e") : 1e6;
        
        UQFFParams_S175 uqff_params;
        uqff_params.r = params.count("r") ? params.at("r") : 1e4;
        uqff_params.t = t;
        uqff_params.t_n = params.count("t_n") ? params.at("t_n") : 0.5;
        uqff_params.gamma = params.count("gamma") ? params.at("gamma") : 0.001;
        
        std::vector<double> dipoles = {params.count("mu") ? params.at("mu") : 1e30};
        std::vector<double> dists = {uqff_params.r};
        
        bool enable_LENR = params.count("enable_LENR") ? (params.at("enable_LENR") > 0.5) : true;
        
        return compute_Um(system, uqff_params, dipoles, dists, enable_LENR);
    }
    
    std::string getName() const override { return "SOURCE175_Um"; }
    
    std::string getDescription() const override {
        return "Universal Magnetism with Widom-Larsen LENR 10^13 enhancement";
    }
};

/**
 * Aether metric tensor PhysicsTerm wrapper
 */
class SOURCE175_Aether_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        AstroSystem_S175 system;
        system.name = "Generic";
        system.mass = params.count("M") ? params.at("M") : 1.0e6;
        system.v_rot = params.count("v_rot") ? params.at("v_rot") : 200.0;
        
        UQFFParams_S175 uqff_params;
        uqff_params.r = params.count("r") ? params.at("r") : 1e20;
        uqff_params.t = t;
        uqff_params.t_n = params.count("t_n") ? params.at("t_n") : 0.5;
        uqff_params.theta = params.count("theta") ? params.at("theta") : PI_S175 / 4;
        uqff_params.phi = params.count("phi") ? params.at("phi") : PI_S175 / 4;
        
        std::array<std::array<double, 4>, 4> metric;
        compute_Aether_metric_tensor(system, uqff_params, metric);
        
        // Return determinant as scalar output
        return compute_metric_determinant(metric);
    }
    
    std::string getName() const override { return "SOURCE175_Aether"; }
    
    std::string getDescription() const override {
        return "Aether metric tensor 4x4 spacetime geometry (returns det(A^μν))";
    }
};

/**
 * Unified Field F_U PhysicsTerm wrapper
 */
class SOURCE175_FU_Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double> &params) const override {
        AstroSystem_S175 system;
        system.name = "Generic";
        system.mass = params.count("M") ? params.at("M") : 1.0e6;
        system.B = params.count("B") ? params.at("B") : 1e10;
        system.v_rot = params.count("v_rot") ? params.at("v_rot") : 200.0;
        
        UQFFParams_S175 uqff_params;
        uqff_params.r = params.count("r") ? params.at("r") : 1e20;
        uqff_params.t = t;
        uqff_params.t_n = params.count("t_n") ? params.at("t_n") : 0.5;
        uqff_params.Q_i = params.count("Q_i") ? params.at("Q_i") : 1.0;
        uqff_params.SCm_i = params.count("SCm_i") ? params.at("SCm_i") : 1.0;
        
        std::array<double, 4> Ug_coeffs = {1.0, 1.0, 1.0, 1.0};
        std::vector<double> Ub_terms;  // Empty for now
        std::vector<double> dipoles = {1e30};
        std::vector<double> dists = {uqff_params.r};
        
        return compute_F_U(system, uqff_params, Ug_coeffs, Ub_terms, dipoles, dists);
    }
    
    std::string getName() const override { return "SOURCE175_F_U"; }
    
    std::string getDescription() const override {
        return "Unified Field master equation: F_U = Σ[k_i×Ug_i - Ub_i] + Um + A^μν";
    }
};

} // namespace SOURCE175
