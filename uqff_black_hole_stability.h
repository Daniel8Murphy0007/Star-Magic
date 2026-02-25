/**
 * @file uqff_black_hole_stability.h
 * @brief UQFF Black Hole Stability Module
 * 
 * Models black hole lifetime and stability in UQFF framework.
 * 
 * Standard Hawking Instability:
 *   τ_standard = 5120πG²M³/(ℏc⁴)
 *   Solar mass: τ ≈ 10⁶⁷ years >> universe age (practically stable)
 *   Primordial BH (10¹² kg): τ ~ Hubble time (evaporating now)
 * 
 * UQFF Stability Enhancement:
 *   τ_UQFF = τ_standard × 1/(1-f_TRZ) × (ρ_UA/ρ_SCm) × exp(U_m/(k_BT_H))
 * 
 *   Step 1: (1-f_TRZ)⁻¹ ≈ 1.11 - Time-reversal reverses 10% of emissions
 *   Step 2: ρ_UA/ρ_SCm ≈ 10 - Aether density suppresses pair separation
 *   Step 3: exp(U_m/k_BT_H) - Magnetic strings block radiation exponentially
 * 
 * Result: Massive BHs become "infinitely" stable in UQFF framework.
 * Sgr A* (4×10⁶ M☉): τ_UQFF → ∞ (stable aether-superstructure)
 * 
 * Reference: Star Magic UQFF Framework (Feb 2026)
 * Author: Daniel Murphy / Star Magic Team
 */

#ifndef UQFF_BLACK_HOLE_STABILITY_H
#define UQFF_BLACK_HOLE_STABILITY_H

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <chrono>
#include <random>
#include <map>
#include <sstream>
#include <iomanip>

// M_PI fallback for Windows MSVC
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @class UQFFBlackHoleStability
 * @brief UQFF-enhanced black hole lifetime calculator
 * 
 * Features:
 * - Standard Hawking evaporation lifetime τ ∝ M³
 * - UQFF stability enhancement via f_TRZ, ρ_UA/ρ_SCm, U_m
 * - Self-expanding: add custom scaling factors
 * - Self-updating: load parameters from config
 * - Self-simulating: τ over mass range
 */
class UQFFBlackHoleStability {
private:
    // Physical constants and UQFF parameters
    std::map<std::string, double> params;

    // Stochastic generator for perturbations
    std::mt19937 rng;
    std::normal_distribution<double> noise_dist;

    // Add-on scalings for τ (function of M, T_H)
    std::vector<std::function<double(double, double)>> additional_scalings;

public:
    /**
     * @brief Constructor with optional random seed
     * @param seed Random seed for stochastic perturbations
     */
    UQFFBlackHoleStability(unsigned int seed = 
        static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count()));

    // ═══════════════════════════════════════════════════════════════════════════
    // PARAMETER ACCESS
    // ═══════════════════════════════════════════════════════════════════════════
    
    /**
     * @brief Get parameter value
     * @param key Parameter name
     * @return Parameter value (0 if not found)
     */
    double get_param(const std::string& key) const;
    
    /**
     * @brief Set parameter value
     * @param key Parameter name
     * @param value New value
     */
    void set_param(const std::string& key, double value);

    // ═══════════════════════════════════════════════════════════════════════════
    // LIFETIME COMPUTATION (Step-by-Step)
    // ═══════════════════════════════════════════════════════════════════════════
    
    /**
     * @brief Compute Hawking temperature
     * @param M Black hole mass [kg]
     * @return T_H = ℏc³/(8πGMk_B) [K]
     */
    double compute_T_H(double M);
    
    /**
     * @brief Standard Hawking evaporation lifetime
     * @param M Black hole mass [kg]
     * @return τ_standard = 5120πG²M³/(ℏc⁴) [s]
     * 
     * Derivation: From dM/dt = -L_H/c² and L_H ∝ M⁻²
     */
    double compute_tau_standard(double M);
    
    /**
     * @brief Step 1: Time-reversal enhancement
     * @param tau_std Standard lifetime from compute_tau_standard
     * @return τ' = τ/(1 - f_TRZ) [s]
     * 
     * Negentropic reversal increases lifetime by reversing ~10%
     */
    double compute_tau_prime(double tau_std);
    
    /**
     * @brief Step 2: Aether density enhancement
     * @param tau_prime τ' from Step 1
     * @return τ'' = τ' × (ρ_UA/ρ_SCm) [s]
     * 
     * Dense [UA] aether suppresses pair fluctuations at horizon
     */
    double compute_tau_double_prime(double tau_prime);
    
    /**
     * @brief Step 3: Magnetic string barrier
     * @param tau_double_prime τ'' from Step 2
     * @param T_H Hawking temperature [K]
     * @return τ_UQFF = τ'' × exp(U_m/(k_B T_H)) [s]
     * 
     * Exponential suppression of evaporation for cold horizons
     */
    double compute_tau_UQFF(double tau_double_prime, double T_H);
    
    /**
     * @brief Full UQFF lifetime with all enhancements + noise + custom scalings
     * @param M Black hole mass [kg]
     * @param noise_level Stochastic noise amplitude (default: 0.001)
     * @return τ_UQFF_full [s]
     */
    double compute_full_tau(double M, double noise_level = 0.001);

    /**
     * @brief Compute stability enhancement factor relative to Hawking
     * @param M Black hole mass [kg]
     * @return τ_UQFF / τ_standard (> 1 means more stable)
     */
    double compute_stability_factor(double M);

    // ═══════════════════════════════════════════════════════════════════════════
    // SELF-EXPANDING / SELF-UPDATING / SELF-SIMULATING
    // ═══════════════════════════════════════════════════════════════════════════
    
    /**
     * @brief Self-expand: Add custom scaling factor
     * @param scaling Function(M, T_H) returning multiplicative factor
     */
    void add_scaling(std::function<double(double, double)> scaling);
    
    /**
     * @brief Self-update: Load parameters from config file
     * @param config_file Path to config file (key=value format)
     */
    void update_from_file(const std::string& config_file);
    
    /**
     * @brief Self-simulate: Compute τ_UQFF over mass range
     * @param M_start Starting mass [kg]
     * @param M_end Ending mass [kg]
     * @param dM Mass step [kg]
     * @param output_file Output CSV file (empty = console output)
     */
    void simulate_over_mass(double M_start, double M_end, double dM, 
                            const std::string& output_file = "");
    
    /**
     * @brief Generate long-form equation string with substituted values
     * @param M Black hole mass [kg]
     * @return Formatted equation string
     */
    std::string long_form_equation(double M);
};

#endif // UQFF_BLACK_HOLE_STABILITY_H
