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
#include <tuple>

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
    
    // ═══════════════════════════════════════════════════════════════════════════
    // AETHER SUPERFLUID DYNAMICS (GPE-derived)
    // ═══════════════════════════════════════════════════════════════════════════
    
    /**
     * @brief Compute healing length (coherence length)
     * @return ξ = √(ℏ²/(2mgρ)) [m]
     * 
     * Determines vortex core size. Cosmic scale: ~10^18 m.
     */
    double compute_healing_length();
    
    /**
     * @brief Compute quantized vortex circulation
     * @param n Vortex quantum number (integer)
     * @return Γ = 2πℏn/m [m²/s]
     */
    double compute_vortex_circulation(int n = 1);
    
    /**
     * @brief Compute vortex energy per unit length
     * @return U_m = (πρℏ²/m) × ln(L/ξ) [J/m]
     */
    double compute_vortex_energy();
    
    /**
     * @brief Compute superfluid velocity from phase gradient
     * @param grad_theta Phase gradient ∇θ [1/m]
     * @return v_s = (ℏ/m) × ∇θ [m/s]
     */
    double compute_superfluid_velocity(double grad_theta);
    
    /**
     * @brief Compute time-reversal modified interaction strength
     * @return g_TRZ = g × (1 - f_TRZ)
     */
    double compute_g_TRZ();
    
    /**
     * @brief Compute quantum pressure term
     * @param rho Density ρ
     * @param grad2_sqrt_rho Laplacian of √ρ
     * @return P_Q = -(ℏ²/2m) × (∇²√ρ/√ρ) [J/m³]
     */
    double compute_quantum_pressure(double rho, double grad2_sqrt_rho);
    
    /**
     * @brief Compute GPE effective potential
     * @param rho Local density |ψ|²
     * @return V_eff = V_ext + g_TRZ × ρ - μ [J]
     */
    double compute_GPE_potential(double rho);
    
    /**
     * @brief Compute Meissner-like expulsion factor at [SCm]
     * @return Factor (1 - ρ_SCm/ρ_UA) for supercurrent
     */
    double compute_meissner_factor();
    
    /**
     * @brief Compute superfluid density time evolution
     * @param rho_0 Initial density
     * @param t Time [s]
     * @return ρ(t) with f_TRZ stabilization
     */
    double compute_superfluid_density_time(double rho_0, double t);
    
    // ═══════════════════════════════════════════════════════════════════════════
    // SUPERFLUID VORTEX DYNAMICS (Kelvin-Helmholtz derived)
    // ═══════════════════════════════════════════════════════════════════════════
    
    /**
     * @brief Compute vortex core size from chemical potential
     * @return ξ = ℏ/√(2mμ) [m]
     */
    double compute_vortex_core_size_mu();
    
    /**
     * @brief Compute azimuthal vortex velocity field
     * @param r Radial distance from vortex center [m]
     * @return v(r) = (ℏn/mr) [m/s] (1/r decay)
     */
    double compute_vortex_velocity_field(double r);
    
    /**
     * @brief Compute self-induced velocity for curved vortex
     * @return v_s,self ≈ (κ/4πR)ln(R/ξ) [m/s]
     */
    double compute_self_induced_velocity();
    
    /**
     * @brief Compute 2D point-vortex Hamiltonian
     * @param positions Vector of (x,y) vortex positions [m]
     * @return H = -(ρκ²/4π)Σln|r_i-r_j| [J]
     */
    double compute_point_vortex_hamiltonian(const std::vector<std::pair<double, double>>& positions);
    
    /**
     * @brief Compute UQFF circulation quantum
     * @return κ_UQFF = h/m_eff [m²/s]
     */
    double compute_kappa_UQFF();
    
    /**
     * @brief Compute UQFF vortex energy with f_TRZ damping
     * @return E_v,UQFF = E_v × (1 - f_TRZ) [J/m]
     */
    double compute_vortex_energy_UQFF();
    
    /**
     * @brief Compute [SCm] damped vortex velocity
     * @param v_s Input velocity [m/s]
     * @return v_s × (1 - B/B_crit) [m/s]
     */
    double compute_vortex_velocity_SCm(double v_s);
    
    /**
     * @brief Compute vortex energy with U_m line tension
     * @return E_v + U_m×ln(b/ξ) [J/m]
     */
    double compute_vortex_energy_with_tension();
    
    /**
     * @brief Compute full UQFF vortex velocity (Step 7 formula)
     * @param r Current radial position [m]
     * @param vortex_positions Vector of other vortex positions (relative) [m]
     * @return Full vortex velocity with self, mutual, tension contributions [m/s]
     */
    double compute_full_vortex_velocity_UQFF(double r, const std::vector<std::pair<double, double>>& vortex_positions);
    
    // ═══════════════════════════════════════════════════════════════════════════
    // VORTEX QUANTIZATION (From wavefunction single-valuedness)
    // ═══════════════════════════════════════════════════════════════════════════
    
    /**
     * @brief Compute effective vortex quantum number with f_TRZ modulation
     * @return n_eff = n × (1 - f_TRZ)
     * 
     * Damps multi-quanta vortices negentropically.
     */
    double compute_n_effective();
    
    /**
     * @brief Compute full UQFF circulation quantum (Step 10 formula)
     * @return κ_UQFF = (h/m_eff) × n × (1-f_TRZ) × (1-B/B_crit) [m²/s]
     * 
     * Includes all modulations: quantization, TRZ damping, SCm pinning.
     * Cosmic scale: ~10^{34} m²/s for m_eff ~ 10^{-68} kg.
     */
    double compute_kappa_UQFF_full();
    
    /**
     * @brief Compute UQFF vortex core size
     * @return ξ_UQFF = ξ × √(ρ_UA/ρ_SCm) [m]
     * 
     * Aether density ratio enlarges cores by ~3.16×.
     */
    double compute_xi_UQFF();
    
    // ═══════════════════════════════════════════════════════════════════════════
    // GPE VORTEX WAVEFUNCTION HELPERS
    // ═══════════════════════════════════════════════════════════════════════════
    
    /**
     * @brief Compute vortex wavefunction amplitude
     * @param r Radial distance from vortex center [m]
     * @return |ψ| = √(μ/g) × tanh(r/ξ)
     * 
     * Core depletion: tanh profile from 0 at center to 1 far away.
     */
    double compute_vortex_wavefunction_amplitude(double r);
    
    /**
     * @brief Compute vortex density profile
     * @param r Radial distance [m]
     * @return ρ(r) = |ψ(r)|²
     */
    double compute_vortex_density(double r);
    
    /**
     * @brief Compute GPE kinetic energy density in k-space
     * @param k_squared k² = kx² + ky² [1/m²]
     * @return K = ℏ²k²/(2m) [J]
     * 
     * Used in split-step Fourier method.
     */
    double compute_GPE_kinetic_energy_density(double k_squared);
    
    /**
     * @brief Compute GPE interaction potential
     * @param density Local density ρ = |ψ|²
     * @return V = g_TRZ × ρ - μ + V_ext [J]
     */
    double compute_GPE_interaction_potential(double density);
    
    // ═══════════════════════════════════════════════════════════════════════════
    // MULTI-VORTEX DYNAMICS
    // ═══════════════════════════════════════════════════════════════════════════
    
    /**
     * @brief Compute multi-vortex wavefunction amplitude
     * @param x X coordinate [m or dimensionless]
     * @param y Y coordinate [m or dimensionless]
     * @param vortices Vector of (x_v, y_v, n) tuples (n>0: vortex, n<0: antivortex)
     * @return |ψ| = √(μ/g) × Π_j tanh(|r - r_j|/ξ)
     */
    double compute_multi_vortex_amplitude(double x, double y,
        const std::vector<std::tuple<double, double, int>>& vortices);
    
    /**
     * @brief Compute multi-vortex density profile
     * @param x X coordinate
     * @param y Y coordinate
     * @param vortices Vector of (x_v, y_v, n) tuples
     * @return ρ(x,y) = |ψ(x,y)|² showing depletions at each core
     */
    double compute_multi_vortex_density(double x, double y,
        const std::vector<std::tuple<double, double, int>>& vortices);
    
    /**
     * @brief Compute pairwise vortex interaction energy
     * @param x1, y1 Position of vortex 1
     * @param n1 Winding number of vortex 1
     * @param x2, y2 Position of vortex 2
     * @param n2 Winding number of vortex 2
     * @return E_pair = -(ρ_s κ₁ κ₂ / 4π) × ln(r₁₂/ξ) [J/m]
     * 
     * Same-sign: repulsion (orbiting), opposite-sign: attraction (annihilation)
     */
    double compute_vortex_pair_interaction_energy(
        double x1, double y1, int n1,
        double x2, double y2, int n2);
    
    /**
     * @brief Compute total interaction energy for multi-vortex system
     * @param vortices Vector of (x_v, y_v, n) tuples
     * @return E_total = Σ_{i<j} E_pair(i,j) [J/m]
     */
    double compute_total_vortex_interaction_energy(
        const std::vector<std::tuple<double, double, int>>& vortices);
    
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
     * @brief Simulate mass evolution with accretion/decay (self-simulate)
     * @param M_start Initial mass [kg]
     * @param t_start Start time [s]
     * @param t_end End time [s]
     * @param dt Time step [s]
     * @param output_file Optional CSV output file
     * 
     * Models: dM/dt = M_dot - M_decay_rate × M × (1 - f_TRZ)
     * UQFF enhancement: f_TRZ damps evaporation negentropically.
     */
    void simulate_mass_evolution(double M_start, double t_start, double t_end, double dt,
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
