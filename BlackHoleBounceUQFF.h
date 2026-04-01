#ifndef BLACK_HOLE_BOUNCE_UQFF_H
#define BLACK_HOLE_BOUNCE_UQFF_H

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @file BlackHoleBounceUQFF.h
 * @class BlackHoleBounceUQFF
 * @brief UQFF extension of Loop Quantum Gravity (LQG) black hole bounces and Loop Quantum Cosmology (LQC).
 *
 * SOURCE: grok_share_fc21e30c24b4.txt — "60. Loop Quantum Gravity Bounce_cpp_12Jan2026"
 * PAPER: 658 | Session 172 | April 1–2, 2026
 *
 * ─── THEORY ───────────────────────────────────────────────────────────────────
 * Loop Quantum Gravity (LQG) quantises space-time into discrete "loops" — quanta
 * of area and volume — that prevent the infinite curvature blow-ups (singularities)
 * of classical GR at the Big Bang and inside black holes.
 *
 * Loop Quantum Cosmology (LQC) applies LQG to the whole universe. The classical
 * Friedmann equation
 *     (ȧ/a)² = (8πG/3)ρ  −  kc²/a²
 * is replaced by the LQC-modified version
 *     (ȧ/a)² = (8πG/3)ρ (1 − ρ/ρ_c)
 * where ρ_c = 0.41 ρ_Pl ≈ 5.16 × 10^{96} kg/m³ is the critical density.
 * As ρ → ρ_c the expansion rate ȧ/a → 0 ⟹ classical singularity is replaced
 * by a smooth "Big Bounce": contraction → bounce → expansion.
 *
 * Minimum scale factor at bounce:
 *     a_min ~ (ħG/c³)^{1/2}  ≈  Planck length ≈ 1.6 × 10^{-35} m
 *
 * Near-bounce approximation (flat, matter-dominated):
 *     a(t) ≈ a_min · cosh(t / t_Pl),   t_Pl = (ħG/c⁵)^{1/2} ≈ 5.39 × 10^{-44} s
 *
 * For black holes, "black bounces" in LQG-inspired models replace the singularity
 * with a bounce to a white hole or another universe (Ashtekar, Olmedo & Singh 2018;
 * Modesto 2025).
 *
 * ─── UQFF EXTENSION ──────────────────────────────────────────────────────────
 * The UQFF extends LQC by raising the effective critical density via the
 * Vacuum Density Series (VDS, PAPER_646) ratio:
 *     ρ_c,UQFF = ρ_c · (1 + ρ_vac_UA / ρ_vac_SCm)  ≈  ρ_c × 11
 * This pushes the bounce to smaller a_min and higher energies, consistent with
 * the five quantum variables of PAPER_657 (Knowledge Base_7).
 *
 * Buoyancy Harmonics Series (PAPER_648) contribution:
 *     The aether density gradient [UA]/[SCm] creates a buoyant restoring force
 *     that prevents singularity formation even before ρ → ρ_c.
 *
 * ─── SELF-EXPANDING FRAMEWORK ────────────────────────────────────────────────
 * Dynamic extra terms can be added at runtime via add_extra_term().
 * Parameters are adjustable at runtime or loaded from a key=value file.
 * Time-evolution simulation is available via simulate().
 *
 * Validated constants:
 *   ρ_vac_UA  = 7.09 × 10^{-36} J/m³
 *   ρ_vac_SCm = 7.09 × 10^{-37} J/m³
 *   f_TRZ     = 0.1
 *   κ         = 0.0005  day^{-1}
 *   [SSq]     = 0.57
 */

class BlackHoleBounceUQFF {
private:
    // ── Physical constants ────────────────────────────────────────────────────
    double G    = 6.6743e-11;          // Gravitational constant [m³ kg^{-1} s^{-2}]
    double c    = 2.998e8;             // Speed of light [m/s]
    double hbar = 1.0546e-34;          // Reduced Planck constant [J s]
    double k_B  = 1.380649e-23;        // Boltzmann constant [J/K]

    // ── Derived Planck quantities ─────────────────────────────────────────────
    double rho_Planck;   // ρ_Pl = c^5 / (ħ G²)  ≈ 5.16 × 10^{96} kg/m³
    double t_Planck;     // t_Pl = sqrt(ħ G / c^5) ≈ 5.39 × 10^{-44} s
    double a_min;        // a_min = sqrt(ħ G / c³) ≈ 1.616 × 10^{-35} m

    // ── LQC critical density ──────────────────────────────────────────────────
    double rho_c_factor = 0.41;        // ρ_c = 0.41 ρ_Pl (Ashtekar & Singh 2011)
    double rho_c;                      // LQC critical density [kg/m³]

    // ── UQFF parameters (Vacuum Density Series PAPER_646) ────────────────────
    double rho_vac_UA  = 7.09e-36;    // [UA] vacuum energy density [J/m³]
    double rho_vac_SCm = 7.09e-37;    // [SCm] vacuum energy density [J/m³]
    double f_TRZ       = 0.1;          // Time-reversal correction factor
    double kappa       = 0.0005;       // κ calibration [day^{-1}]
    double SSq         = 0.57;         // [SSq] calibration constant

    // ── Curvature flag ────────────────────────────────────────────────────────
    int k_curv = 0;                    // 0=flat, +1=closed, -1=open

    // ── Dynamic extra terms (self-expanding framework) ────────────────────────
    std::vector<std::function<double(double /*a*/, double /*rho*/)>> extra_terms;

    // ── Stochastic perturbation ───────────────────────────────────────────────
    std::mt19937 rng;
    std::normal_distribution<double> noise_dist{0.0, 1.0};

public:
    // ── Constructor ───────────────────────────────────────────────────────────
    explicit BlackHoleBounceUQFF(
        unsigned int seed = static_cast<unsigned int>(
            std::chrono::system_clock::now().time_since_epoch().count()));

    // ── Core equations ────────────────────────────────────────────────────────

    /**
     * @brief LQC-modified Friedmann: (ȧ/a)² = (8πG/3)ρ(1 − ρ/ρ_c) − kc²/a²
     * @param rho  Energy density [kg/m³]
     * @param a    Scale factor [m]
     * @return (ȧ/a)² [s^{-2}]
     */
    double compute_LQC_friedmann(double rho, double a) const;

    /**
     * @brief Classical Friedmann (k=0): (ȧ/a)² = (8πG/3)ρ
     * @param rho  Energy density [kg/m³]
     * @return (ȧ/a)² [s^{-2}]
     */
    double compute_classical_friedmann(double rho) const;

    /**
     * @brief UQFF-extended critical density:
     *        ρ_c,UQFF = ρ_c · (1 + ρ_vac_UA / ρ_vac_SCm)
     * @return ρ_c,UQFF [kg/m³]
     */
    double compute_rho_c_UQFF() const;

    /**
     * @brief Near-bounce approximation: a(t) ≈ a_min · cosh(t / t_Pl)
     * @param t  Time offset from bounce [s]
     * @return Scale factor [m]
     */
    double compute_scale_factor_near_bounce(double t) const;

    /**
     * @brief UQFF-corrected scale factor near bounce:
     *        a(t,UQFF) = a_min · cosh(t / t_Pl) · (1 + f_TRZ · ρ_vac_UA / ρ_vac_SCm)^{1/3}
     * @param t  Time offset from bounce [s]
     * @return UQFF scale factor [m]
     */
    double compute_scale_factor_UQFF(double t) const;

    /**
     * @brief Effective equation-of-state UQFF correction at bounce:
     *        w_eff = −1 + (1+f_TRZ)(ρ_vac_UA/ρ_vac_SCm)·κ·[SSq]
     * @return Dimensionless effective EOS
     */
    double compute_effective_eos() const;

    /**
     * @brief Rate of density change near bounce (numerical derivative):
     *        drho/dt ≈ −3H(1+w)ρ where H = ȧ/a
     * @param rho  Current density [kg/m³]
     * @param a    Current scale factor [m]
     * @return drho/dt [kg m^{-3} s^{-1}]
     */
    double compute_rho_rate(double rho, double a) const;

    // ── Self-expanding interface ───────────────────────────────────────────────

    /**
     * @brief Register a new additive term f(a, rho) to the modified Friedmann.
     */
    void add_extra_term(std::function<double(double, double)> term);

    // ── Self-update interface ─────────────────────────────────────────────────

    /**
     * @brief Load parameters from key=value file.
     *        Keys: G, c, hbar, k_B, rho_c_factor, rho_vac_UA, rho_vac_SCm,
     *              f_TRZ, kappa, SSq, k_curv
     */
    void update_from_file(const std::string& filename);

    // ── Self-simulate interface ───────────────────────────────────────────────

    /**
     * @brief Evolve scale factor a(t) and density ρ(t) over a time range.
     *        Uses simple Euler integration of the LQC Friedmann equation.
     * @param a0          Initial scale factor [m]
     * @param rho0        Initial density [kg/m³]
     * @param t_start     Start time [s]
     * @param t_end       End time [s]
     * @param dt          Time step [s]
     * @param output_file CSV path (empty → console)
     */
    void simulate(double a0, double rho0, double t_start, double t_end,
                  double dt, const std::string& output_file = "") const;

    // ── Accessors ─────────────────────────────────────────────────────────────
    double get_rho_c()       const { return rho_c; }
    double get_rho_c_UQFF()  const { return compute_rho_c_UQFF(); }
    double get_t_Planck()    const { return t_Planck; }
    double get_a_min()       const { return a_min; }
    double get_rho_Planck()  const { return rho_Planck; }

    void print_explanations() const;
};

#endif // BLACK_HOLE_BOUNCE_UQFF_H
