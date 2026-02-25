/**
 * @file uqff_framework.h
 * @brief UQFF Framework for Astronomical Modeling - Header
 * 
 * Class implementing the UQFF framework for astronomical modeling.
 * Captures core principles, math structure, equations in comments.
 * Allows self-update (param file load), self-expand (add terms to MUGE),
 * self-simulate (evolve g over t).
 * 
 * Core Principles:
 * - Gravity emergent from quantum-superconductive in aether
 * - [UA] superfluid ρ≈7.09e-36 J/m³ modulates energy/light
 * - [SCm] Type-II B_crit≈10^11 T for extremes
 * - f_TRZ=0.1 negentropic (echo contractions/flares)
 * - U_m km-scale magnetic/THz strings
 * - F_env tidal/rad/SN/BH/eta
 * - Quantum ħ/√(Δx Δp) + ψ coherence
 * - Scalable reactor-galaxy via U_g1-4
 * 
 * Math (MUGE):
 *   g(r,t) = [G M(t)/r(t)²] (1+H(t,z)) (1-B(t)/B_crit) (1+F_env(t))
 *          + ∑U_g + U_i + Λc²/3
 *          + [ħ/√(Δx·Δp)] ∫ψ* H ψ dV (2π/t_Hub)
 *          + ρ_fluid V g_local
 *          + (M_vis + M_DM) (δ_ρ/ρ + 3GM(t)/r(t)³)
 * 
 * Reference: Star Magic UQFF Framework (Feb 2026)
 * Cross-platform: Windows (MSVC) + Linux (g++)
 */

#ifndef UQFF_FRAMEWORK_H
#define UQFF_FRAMEWORK_H

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <map>
#include <chrono>
#include <random>

// M_PI fallback for Windows MSVC
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @class UQFFFramework
 * @brief Core UQFF implementation with MUGE gravity computation
 * 
 * Features:
 * - Full MUGE computation with all 10+ terms
 * - Quantum coherence modeling near horizons
 * - Self-expanding: add custom physics terms
 * - Self-updating: load parameters from config file
 * - Self-simulating: time evolution of gravity
 */
class UQFFFramework {
private:
    // Stochastic generator for perturbations
    std::mt19937 rng;
    std::normal_distribution<double> noise_dist;

    // Add-on terms (functions of r, t)
    std::vector<std::function<double(double, double)>> additional_terms;

    // Framework explanations
    std::vector<std::string> explanations;

    // Initialize explanations from core principles
    void init_explanations();

public:
    // Parameters map (public for direct access and calibration)
    std::map<std::string, double> params;

    /**
     * @brief Constructor with optional random seed
     * @param seed Random seed for stochastic terms (default: system time)
     */
    UQFFFramework(unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count());

    /**
     * @brief Compute quantum coherence contribution
     * @param r Radial distance [m]
     * @param t Time [s]
     * @return Quantum coherence term contribution
     * @brief Time-reversal correction factor
     * @param base_value Input value to correct
     * @return Corrected value with f_TRZ applied
     * 
     * Applies negentropic correction for time-reversal zone effects.
     * Enhancement: value × (1 + f_TRZ), Suppression: value × (1 - f_TRZ)
     */
    double time_reversal_correction(double base_value);
    
    /**
     * @brief Compute effective sigma with aether damping
     * @return σ_eff = σ × (1 - ρ_SCm/ρ_UA) [m]
     * 
     * Aether vacuum densities modify the coherence length scale.
     */
    double compute_sigma_effective();
    
    /**
     * @brief Compute wavefunction normalization amplitude
     * @return A = (√(2π) σ_eff)^(-1/2)
     * 
     * Ensures ∫|ψ|² dr = 1 for Gaussian wavepacket.
     */
    double compute_normalization_amplitude();
    
    /**
     * @brief Compute full UQFF quantum coherence
     * @param r Radial distance from center [m]
     * @param t Time [s]
     * @return C_UQFF coherence measure
     * 
     * Full UQFF derivation:
     *   ψ(r,t) = A exp(-(r-r_h)²/2σ_eff²) exp(-i 2πft(1+f_TRZ))
     *   C_UQFF = (ℏ²/2m σ_eff²) × |cos(2πft(1+f_TRZ))| × exp(-U_m/k_BT) × gaussian
     * 
     * Set use_full_uqff_coherence=0 for simple model (original).
     */
    double quantum_coherence(double r, double t);

    /**
     * @brief Compute full MUGE gravity
     * @param r Radial distance [m]
     * @param t Time [s]
     * @param noise_level Optional stochastic noise amplitude (default: 0.001)
     * @return Total gravity g(r,t) [m/s²]
     * 
     * Full equation:
     *   g(r,t) = base × (1+H) × (1-B/B_crit) × (1+F_env)
     *          + ∑U_g + U_i + Λc²/3 + quantum + fluid + DM
     */
    double compute_MUGE(double r, double t, double noise_level = 0.001);

    /**
     * @brief Add custom physics term (self-expand)
     * @param term Function of (r, t) returning gravity contribution
     */
    void add_term(std::function<double(double, double)> term);

    /**
     * @brief Load parameters from config file (self-update)
     * @param config_file Path to config file (key=value format)
     */
    void update_from_file(const std::string& config_file);

    /**
     * @brief Simulate gravity evolution over time (self-simulate)
     * @param r Fixed radial distance [m]
     * @param t_start Start time [s]
     * @param t_end End time [s]
     * @param dt Time step [s]
     * @param output_file Optional CSV output file
     */
    void simulate_evolution(double r, double t_start, double t_end, double dt, 
                           const std::string& output_file = "");

    /**
     * @brief Display framework explanations
     */
    void display_explanations();

    /**
     * @brief Set a single parameter
     * @param key Parameter name
     * @param value Parameter value
     */
    void set_param(const std::string& key, double value);

    /**
     * @brief Get a parameter value
     * @param key Parameter name
     * @return Parameter value (0.0 if not found)
     */
    double get_param(const std::string& key) const;

    /**
     * @brief Check if parameter exists
     * @param key Parameter name
     * @return true if parameter exists
     */
    bool has_param(const std::string& key) const;

    /**
     * @brief Get number of additional terms
     * @return Count of custom terms added
     */
    size_t term_count() const;

    /**
     * @brief Clear all additional terms
     */
    void clear_terms();

    /**
     * @brief Export parameters to file
     * @param output_file Path to output file
     */
    void export_params(const std::string& output_file) const;
};

#endif // UQFF_FRAMEWORK_H
