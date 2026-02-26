// uqff_wormhole_transverse_time.h
// UQFF Wormhole Transverse Time Module
// Computes traversal time τ_UQFF through stable wormholes

#ifndef UQFF_WORMHOLE_TRANSVERSE_TIME_H
#define UQFF_WORMHOLE_TRANSVERSE_TIME_H

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <iostream>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <chrono>
#include <random>
#include <sstream>
#include <iomanip>

/**
 * WormholeTransverseTime: UQFF wormhole traversal time calculator
 * 
 * Derivation: Modulates GR τ ≈ l/c ~ r_throat/c = 2GM/c³ with UQFF elements.
 * 
 * Step 1: τ_base = 2GM/c³ from r_throat = 2GM/c²
 * Step 2: v_eff = c(1 - ρ_SCm/ρ_UA), ρ ratio ≈ 0.1 reduces to 0.9c
 * Step 3: τ' = τ_base / (1 - ρ_SCm/ρ_UA)
 * Step 4: τ'' = τ'(1 + f_TRZ), f_TRZ ≈ 0.1 adds drag
 * Step 5: τ_UQFF = τ'' × exp(U_m / (k_B × T_eff))
 * 
 * Full equation:
 *   τ_UQFF = [2GM/c³] / (1 - ρ_SCm/ρ_UA) × (1 + f_TRZ) × exp(U_m / (k_B × T_eff))
 * 
 * Numerical example: Sgr A* (4×10⁶ M_sun)
 *   τ_base ≈ 1×10⁻⁴ s
 *   ρ ratio = 0.1 → τ' ≈ 1.11×10⁻⁴ s
 *   f_TRZ = 0.1 → τ'' ≈ 1.22×10⁻⁴ s
 *   U_m/(k_B T_eff) = 1 → τ_UQFF ≈ 3.3×10⁻⁴ s (longer traversal)
 * 
 * Physics: Traversable stable wormholes with UQFF-induced delays
 */
class WormholeTransverseTime {
private:
    // Physical constants
    double G;                       // Gravitational constant: 6.6743e-11 m³ kg⁻¹ s⁻²
    double c;                       // Speed of light: 2.998e8 m/s
    double hbar;                    // Reduced Planck: 1.0545718e-34 J·s
    double k_B;                     // Boltzmann constant: 1.380649e-23 J/K

    // UQFF vacuum densities
    double rho_vac_SCm;             // Superconductive vacuum: 7.09e-37 J/m³
    double rho_vac_UA;              // Aether vacuum: 7.09e-36 J/m³

    // UQFF parameters
    double f_TRZ;                   // Time-reversal zone factor: 0.1
    double U_m;                     // Magnetic barrier energy (J)
    double mu_j;                    // Magnetic permeability factor
    double gamma;                   // Temporal decay rate
    double t_n;                     // Normalized time parameter

    // UQFF scaling factors (from uqff_wormhole_formation)
    double kappa_UQFF;              // Energy reduction factor: 1e-60
    double lambda_UQFF;             // Magnetic scaling factor: 1e-9
    double T_eff_floor;             // Effective temperature floor: 1e16 K

    // Stochastic generator
    std::mt19937 rng;
    std::normal_distribution<double> noise_dist;

    // Add-on modifications for τ_UQFF(M)
    std::vector<std::function<double(double)>> additional_mods;

    // Derivation explanations
    std::vector<std::string> explanations;

public:
    // Constructor with optional random seed
    WormholeTransverseTime(unsigned int seed = 0);

    // Core computation methods
    double compute_tau_base(double M);
    double compute_v_eff();
    double compute_tau_prime(double tau_base);
    double compute_tau_double_prime(double tau_prime);
    double compute_T_H(double M);
    double compute_T_eff(double T_H);
    double compute_U_m(double M);
    double compute_tau_UQFF(double tau_double_prime, double U_m_val, double T_eff);
    double compute_full_tau(double M, double noise_level = 0.0);

    // Self-expansion: Add custom modification
    void add_mod(std::function<double(double)> mod);

    // Self-update: Load parameters from file
    void update_from_file(const std::string& config_file);

    // Self-simulate: Compute τ_UQFF over mass range
    void simulate_over_mass(double M_start, double M_end, double dM, const std::string& output_file = "");

    // Display derivation explanations
    void display_explanations();

    // Generate full derivation string
    std::string generate_derivation(double M);

    // Run validation tests
    int run_tests();

    // Getters for current parameters
    double get_f_TRZ() const { return f_TRZ; }
    double get_U_m() const { return U_m; }
    double get_rho_ratio() const { return rho_vac_SCm / rho_vac_UA; }
};

#endif // UQFF_WORMHOLE_TRANSVERSE_TIME_H
