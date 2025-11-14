/**
 * ================================================================================================
 * SystemCatalogue.hpp - UQFF Master Equation Catalogue
 * ================================================================================================
 * 
 * Extracted from: Source10.cpp (lines 695-990)
 * Purpose: Central repository for all UQFF system parameters and physics equations
 * Phase: 1, Week 1 - Foundation Extraction
 * 
 * Contains:
 * - SystemParams struct: 63+ physics terms for UQFF calculations
 * - System database: 26+ astrophysical systems (SN 1006, Eta Carinae, etc.)
 * - Core calculations: F_U_Bi_i (buoyancy), compressed_g (gravity field)
 * 
 * Dependencies:
 * - C++17 standard library
 * - <cmath> for physics calculations
 * - <map>, <string> for system catalogue
 * 
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef SYSTEMCATALOGUE_HPP
#define SYSTEMCATALOGUE_HPP

#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <ctime>
#include <ctime>

// Physical constants (from Source10.cpp)
const double PI = 3.141592653589793;
const double c = 299792458.0;            // Speed of light (m/s)
const double G = 6.67430e-11;            // Gravitational constant (m^3 kg^-1 s^-2)
const double h_bar = 1.054571817e-34;    // Reduced Planck constant (J·s)
const double Msun = 1.989e30;            // Solar mass (kg)
const double E_LEP = 208e9 * 1.602e-19;  // LEP energy (J)
const double h_gw = 1e-21;               // Gravitational wave strain

// UQFF layer constants
const int num_layers = 26;               // 26 quantum layers
const double layer_scale_factor = 1e12;  // Trillion-scale amplification

namespace UQFFCatalogue {

/**
 * SystemParams: Complete parameter set for UQFF calculations
 * 
 * Structure contains all 63+ physics terms required for:
 * - Buoyancy force calculations (F_U_Bi_i)
 * - Compressed gravity field (compressed_g)
 * - Multi-scale quantum-relativistic analysis
 * 
 * Fields organized by category:
 * - Basic: name, M, r, T, L_X, B0, omega0, theta_deg, t, v
 * - Vacuum: rho_vac_UA, rho_vac_SCm
 * - DPM: DPM_stability, DPM_momentum, DPM_gravity, DPM_life
 * - UQFF terms: k_LENR, k_act, k_DE, k_neutron, sigma_n, k_rel, F_rel
 * - Vacuum/THz: k_vac, k_thz, omega_thz, neutron_factor, conduit_scale
 * - Conduit: k_conduit, water_state, H_abundance
 * - Quantum: k_spooky, string_wave, Q_wave
 * - Hydrogen: Delta_k_eta, V_void_fraction
 * - Adjustment: alpha_i
 * - Results: F_U_Bi_i (computed), term1-4, std_scale
 * - Density: rho_astro, rho_LEP
 */
struct SystemParams {
    // Basic system identification and properties
    std::string name;                // System name (e.g., "SN 1006", "Eta Carinae")
    double M;                        // Mass (kg)
    double r;                        // Characteristic radius (m)
    double T;                        // Temperature (K)
    double L_X;                      // X-ray luminosity (W)
    double B0;                       // Magnetic field strength (T)
    double omega0;                   // Characteristic angular frequency (rad/s)
    double theta_deg;                // Angle parameter (degrees, default 45.0)
    double t;                        // Time parameter (s)
    double v;                        // Velocity (m/s)
    
    // Vacuum energy densities
    double rho_vac_UA;               // Universal Aether vacuum density (J/m³, default 7.09e-36)
    double rho_vac_SCm;              // Superconductive magnetism vacuum density (J/m³, default 7.09e-37)
    
    // DPM (Dipole Momentum) parameters
    double DPM_stability;            // Stability factor (default 0.01)
    double DPM_momentum;             // Momentum factor (default 0.93)
    double DPM_gravity;              // Gravity factor (default 1.0)
    
    // UQFF coupling constants
    double k_LENR;                   // Low-Energy Nuclear Reactions coupling (default 1e-10)
    double k_act;                    // Activation frequency coupling (default 1e-6)
    double k_DE;                     // Directed energy coupling (default 1e-30)
    double k_neutron;                // Neutron drop coupling (default 1e10)
    double sigma_n;                  // Neutron cross-section (default 1e-4)
    double k_rel;                    // Relativistic coupling (default 1e-10)
    double F_rel;                    // Relativistic force from LEP (N, default 4.30e33)
    
    // Vacuum repulsion and THz shock
    double k_vac;                    // Vacuum repulsion coupling (default 1e-30)
    double k_thz;                    // THz shock coupling (default 1e-10)
    double omega_thz;                // THz frequency (rad/s, default 2π×1e12)
    double neutron_factor;           // Neutron stability factor (1.0 stable, 0.0 unstable)
    double conduit_scale;            // Material abundance scale (default 10.0)
    
    // Conduit parameters
    double k_conduit;                // Conduit coupling (default 1e-22)
    double water_state;              // Water phase state (1.0 stable incompressible, variable in plasma)
    double H_abundance;              // Hydrogen abundance factor (default 10.0)
    
    // Quantum wave parameters
    double k_spooky;                 // Spooky action coupling (default 1e-30)
    double string_wave;              // String wave parameter (default 1e-10)
    double Q_wave;                   // Wave quality factor (default 1.0)
    
    // Hydrogen calculation parameters
    double Delta_k_eta;              // Eta parameter from hydrogen calc (default 7.25e8)
    double V_void_fraction;          // Void fraction (default 0.2)
    
    // Adjustment factor
    double alpha_i;                  // Layer adjustment factor (default 0.01)
    
    // Computed results (stored after calculation)
    double F_U_Bi_i;                 // Buoyancy force (N) - COMPUTED
    double term1;                    // Intermediate term 1
    double term2;                    // Intermediate term 2
    double term3;                    // Intermediate term 3
    double term4;                    // Intermediate term 4
    double std_scale;                // Probabilistic Monte Carlo standard deviation scale (default 0.1)
    double DPM_life;                 // DPM lifespan parameter
    
    // Density scaling
    double rho_astro;                // Astrophysical density (g/cm³, default 1e-17)
    double rho_LEP;                  // LEP reference density (g/cm³, default 1e-25)
    
    // Default constructor with standard UQFF values
    SystemParams() : 
        name("Unknown"),
        M(1e30), r(1e16), T(1e7), L_X(1e30), B0(1e-9), omega0(0.0),
        theta_deg(45.0), t(1e10), v(1e6),
        rho_vac_UA(7.09e-36), rho_vac_SCm(7.09e-37),
        DPM_stability(0.01), DPM_momentum(0.93), DPM_gravity(1.0),
        k_LENR(1e-10), k_act(1e-6), k_DE(1e-30), k_neutron(1e10), sigma_n(1e-4),
        k_rel(1e-10), F_rel(4.30e33),
        k_vac(1e-30), k_thz(1e-10), omega_thz(2.0 * PI * 1e12),
        neutron_factor(1.0), conduit_scale(10.0),
        k_conduit(1e-22), water_state(1.0), H_abundance(10.0),
        k_spooky(1e-30), string_wave(1e-10), Q_wave(1.0),
        Delta_k_eta(7.25e8), V_void_fraction(0.2),
        alpha_i(0.01),
        F_U_Bi_i(0.0), term1(0.0), term2(0.0), term3(0.0), term4(0.0),
        std_scale(0.1), DPM_life(0.0),
        rho_astro(1e-17), rho_LEP(1e-25)
    {}
};

// System Catalogue: Database of astrophysical systems with parameters
// Initialized in SystemCatalogue.cpp
extern std::map<std::string, SystemParams> systems;

// Core UQFF calculation functions (implemented in SystemCatalogue.cpp)

/**
 * compute_E_cm: Compute center-of-mass energy scaling
 * 
 * Long-form: E_cm,astro,local,adj,eff,enhanced = E_LEP * sqrt(ρ_astro / ρ_LEP) * Q_wave
 * 
 * @param p System parameters
 * @return Scaled energy (J)
 */
double compute_E_cm(const SystemParams& p);

/**
 * dpm_life_proportion: Compute DPM lifespan proportion
 * 
 * Calculates: [(SCM):(UA'):(F_U_Bi_i):(Belly Button)] proportion
 * 
 * @param p System parameters
 * @return Lifespan proportion (dimensionless)
 */
double dpm_life_proportion(const SystemParams& p);

/**
 * F_U_Bi_i: Compute indexed buoyancy force
 * 
 * Integrates all UQFF terms:
 * - Vacuum repulsion (k_vac * Δρ_vac * M * v)
 * - THz shock (k_thz * (ω_thz/ω_0)² * neutron_factor * conduit_scale)
 * - Conduit (k_conduit * H_abundance * water_state * neutron_factor)
 * - Spooky action (k_spooky * string_wave/ω_0)
 * - LENR (k_LENR * 1.25e12 Hz)
 * - Activation (k_act * 300 Hz)
 * - Directed energy (k_DE * L_X/(4πr²))
 * - Resonance (k_act * cos(ω_0 * t))
 * - Neutron drop (k_neutron * σ_n * neutron_factor)
 * - Relativistic (k_rel * F_rel)
 * 
 * Applies 26-layer amplification and E_cm scaling
 * Includes probabilistic Monte Carlo integration
 * 
 * @param p System parameters
 * @return Buoyancy force (N)
 */
double F_U_Bi_i(const SystemParams& p);

/**
 * compressed_g: Compute compressed gravity field (26-layer sum)
 * 
 * Sums over 26 quantum layers:
 * g(r,t) = Σ(i=1 to 26) [Ug1_i + Ug2_i + Ug3_i + Ug4i_i]
 * 
 * Where:
 * - Ug1_i: Dipole/spin term (ℏc/r_i² * Q_i * SCm_i * ρ_vac * f_TRZ)
 * - Ug2_i: Superconductor quality (E_DPM/r_i² * SCm_i * f_Um)
 * - Ug3_i: Resonance/magnetic disk (ℏω_0/2 * Q_i * cos(2πf_i*t) / r_i)
 * - Ug4i_i: Adjusted Newtonian (G*M_i/r_i² * (1+α_i) * SCm_i)
 * 
 * @param p System parameters
 * @return Effective gravity field (m/s²)
 */
double compressed_g(const SystemParams& p);

/**
 * F_jet_rel: Relativistic jet force
 * @param p System parameters
 * @return Jet force (N)
 */
double F_jet_rel(const SystemParams& p);

/**
 * E_acc_rel: Relativistic acceleration energy
 * @param p System parameters
 * @return Acceleration energy (J)
 */
double E_acc_rel(const SystemParams& p);

/**
 * F_drag_rel: Relativistic drag force
 * @param p System parameters
 * @return Drag force (N)
 */
double F_drag_rel(const SystemParams& p);

/**
 * F_gw_rel: Gravitational wave force (placeholder)
 * @param p System parameters
 * @return GW force (N, currently 0.0)
 */
double F_gw_rel(const SystemParams& p);

/**
 * validation_pipeline: Simulate cross-reference validation
 * Prints Chandra/GW/JWST observation suggestions
 * @param p System parameters
 */
void validation_pipeline(const SystemParams& p);

/**
 * initializeSystemCatalogue: Load system catalogue (returns map)
 * Populates and returns a map with 20 astrophysical systems
 * @return Map of system name to SystemParams
 */
std::map<std::string, SystemParams> initializeSystemCatalogue();

/**
 * initialize_systems: Load system catalogue
 * Populates the global systems map with 26+ astrophysical systems
 * Called once during engine initialization
 */
void initialize_systems();

/**
 * get_system: Retrieve system parameters by name
 * @param name System name (e.g., "SN 1006")
 * @return SystemParams reference (throws if not found)
 */
SystemParams& get_system(const std::string& name);

/**
 * list_systems: Get all system names
 * @return Vector of system names
 */
std::vector<std::string> list_systems();

} // namespace UQFFCatalogue

#endif // SYSTEMCATALOGUE_HPP
