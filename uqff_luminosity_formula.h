/**
 * @file uqff_luminosity_formula.h
 * @brief UQFF Luminosity Formula for Black Hole Radiation
 * 
 * Extends Hawking L_H = ℏc⁶/(15360πG²M²) with UQFF elements:
 * 
 * Derivation Steps:
 * Step 1: Horizon A=4πr_s², r_s=2GM/c², L_H=σAT_H⁴, σ=π²k_B⁴/(60ℏ³c²), T_H=ℏc³/(8πGMk_B)
 * Step 2: L' = L_H × (1 - f_TRZ), f_TRZ≈0.1 negentropic subtraction
 * Step 3: L'' = L' × (1 - ρ_SCm/ρ_UA), aether ratio damps ~90%
 * Step 4: L_UQFF = L'' × exp(-U_m/(k_B T_H)), magnetic string barrier
 * Step 5: Full L_UQFF = ℏc⁶/(15360πG²M²) × (1-f_TRZ) × (1-ρ_SCm/ρ_UA) × exp(-U_m/(k_BT_H))
 * 
 * Numerical Example: M=8.552e36 kg (Sgr A*), L_H≈1e-5 W, L_UQFF≈3e-6 W (reduced)
 * 
 * Advances UQFF: Quantifies aether/superconductive damping; stable BHs;
 * links THz/q-scope validation.
 * 
 * Reference: Star Magic UQFF Framework (Feb 2026)
 * Author: Daniel Murphy / Star Magic Team
 */

#ifndef UQFF_LUMINOSITY_FORMULA_H
#define UQFF_LUMINOSITY_FORMULA_H

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <chrono>
#include <random>
#include <map>

// M_PI fallback for Windows MSVC
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @class UQFFLuminosityFormula
 * @brief UQFF-extended Hawking radiation luminosity calculator
 * 
 * Features:
 * - Standard Hawking luminosity L_H
 * - Time-reversal correction (f_TRZ)
 * - Aether-superconductive damping (ρ_SCm/ρ_UA)
 * - Magnetic string barrier (U_m)
 * - Self-expanding: add custom dampings
 * - Self-updating: load parameters from config
 * - Self-simulating: L over mass range
 */
class UQFFLuminosityFormula {
private:
    // Physical constants and UQFF parameters
    std::map<std::string, double> params;

    // Stochastic generator for perturbations
    std::mt19937 rng;
    std::normal_distribution<double> noise_dist;

    // Add-on dampings for L (function of M, T_H)
    std::vector<std::function<double(double, double)>> additional_dampings;

public:
    /**
     * @brief Constructor with optional random seed
     * @param seed Random seed for stochastic perturbations
     */
    UQFFLuminosityFormula(unsigned int seed = 
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
    // LUMINOSITY COMPUTATION (Step-by-Step)
    // ═══════════════════════════════════════════════════════════════════════════
    
    /**
     * @brief Compute Hawking temperature
     * @param M Black hole mass [kg]
     * @return T_H = ℏc³/(8πGMk_B) [K]
     */
    double compute_T_H(double M);
    
    /**
     * @brief Step 1: Standard Hawking luminosity
     * @param M Black hole mass [kg]
     * @return L_H = ℏc⁶/(15360πG²M²) [W]
     * 
     * From Stefan-Boltzmann law: L = σAT⁴
     * With A = 4πr_s² = 16π(GM/c²)², T = T_H
     */
    double compute_L_H(double M);
    
    /**
     * @brief Step 2: Time-reversal correction
     * @param L_H Hawking luminosity from Step 1
     * @return L' = L_H × (1 - f_TRZ) [W]
     * 
     * Negentropic reversal reduces emission by f_TRZ ≈ 10%
     */
    double compute_L_prime(double L_H);
    
    /**
     * @brief Step 3: Aether-superconductive damping
     * @param L_prime L' from Step 2
     * @return L'' = L' × (1 - ρ_SCm/ρ_UA) [W]
     * 
     * Dense [UA] vacuum suppresses radiation by ~90%
     */
    double compute_L_double_prime(double L_prime);
    
    /**
     * @brief Step 4: Magnetic string barrier
     * @param L_double_prime L'' from Step 3
     * @param T_H Hawking temperature [K]
     * @return L_UQFF = L'' × exp(-U_m/(k_B T_H)) [W]
     * 
     * Exponential suppression by magnetic string energy barrier
     */
    double compute_L_UQFF(double L_double_prime, double T_H);
    
    /**
     * @brief Full UQFF luminosity with all corrections + noise + custom dampings
     * @param M Black hole mass [kg]
     * @param noise_level Stochastic noise amplitude (default: 0.001)
     * @return L_UQFF_full [W]
     */
    double compute_full_L(double M, double noise_level = 0.001);

    // ═══════════════════════════════════════════════════════════════════════════
    // SELF-EXPANDING / SELF-UPDATING / SELF-SIMULATING
    // ═══════════════════════════════════════════════════════════════════════════
    
    /**
     * @brief Self-expand: Add custom damping factor
     * @param damping Function(M, T_H) returning multiplicative factor
     * 
     * Example: Add quantum coherence modulation
     *   formula.add_damping([](double M, double T_H) { return 1.0 + 0.01/(M*T_H); });
     */
    void add_damping(std::function<double(double, double)> damping);
    
    /**
     * @brief Self-update: Load parameters from config file
     * @param config_file Path to config file (key=value format)
     */
    void update_from_file(const std::string& config_file);
    
    /**
     * @brief Self-simulate: Compute L_UQFF over mass range
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

#endif // UQFF_LUMINOSITY_FORMULA_H
