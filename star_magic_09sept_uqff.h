// star_magic_09sept_uqff.h
// UQFF 2.0 Header — Star Magic_09Sept2025 Module
// Session 100 — PAPER_368 / PAPER_369 / PAPER_370
// ©2025 Daniel T. Murphy — All Rights Reserved

#pragma once

// Forward declarations and public API for STAR_MAGIC_09SEPT_UQFF_MODULE.cpp
// Include this header to access the StarMagic09Sept namespace from other modules.

namespace StarMagic09Sept {

    // ---- CelestialBody (defined in .cpp) ----
    struct CelestialBody;

    // ---- Factory functions ----
    CelestialBody make_Sun();
    CelestialBody make_Earth();
    CelestialBody make_Jupiter();
    CelestialBody make_Neptune();

    // ---- Core equation functions ----
    double compute_Ereact(const CelestialBody& b, double t);
    double compute_Ug1(const CelestialBody& b, double r, double t, double tn);
    double compute_Ug2(const CelestialBody& b, double r, double t, double tn);
    double compute_Ug3(const CelestialBody& b, double r, double t, double tn);

    /**
     * PAPER_368 — Ug4 Vacuum Energy ΛCDM Galactic BH Coupling
     *
     * Ug4 = k4 × ρ_v × C_conc × Mbh/dg × exp(−α×t) × cos(π×tn) × (1+f_feedback)
     *
     * Parameters (canonical value):
     *   k4        = 2.0              Vacuum coupling constant
     *   rho_v     = 6×10⁻²⁷ kg/m³  ΛCDM dark-energy mass density
     *   C_conc    = 1.0              Vacuum concentration factor
     *   Mbh       = 8.15×10³⁶ kg    Galactic centre BH mass (SgrA*, EHT 2022)
     *   dg        = 2.55×10²⁰ m     Distance to galactic centre
     *   alpha     = 0.001 day⁻¹     Non-linear time decay
     *   f_feedback= 0.1             AGN feedback coupling
     *
     * Numerical: Ug4(t=0, tn=0) ≈ 4.22×10⁻¹⁰ m/s²
     *
     * DISTINCT from Ug4VacuumMediatedCalculator (Thread f3c55f52):
     *   - f3c55f52: k4=1×10⁻⁴⁰, rho in J/m³, [SCm] multiplier
     *   - This form: k4=2.0, rho in kg/m³ (ΛCDM ρ_DE), C_concentration multiplier
     *
     * @param t  time in days
     * @param tn negative-time factor (days)
     * @return Ug4 in m/s²
     */
    double compute_Ug4(double t, double tn);

    double compute_Ubi(double Ugi, double tn);
    double compute_Um(const CelestialBody& b, double t, double tn, double rj);
    double compute_A_mu_nu_trace(double tn);
    double compute_FU(const CelestialBody& b, double r, double t, double tn);

    // ---- PAPER_369: Navier-Stokes Quasar Jet ----
    /**
     * PAPER_369 — Navier-Stokes Stable Fluids UQFF Quasar Jet
     *
     * Simulates 2D incompressible flow on N×N grid using Jos Stam (1999)
     * "Stable Fluids" method. SCm velocity is used as jet force input.
     *
     * Step sequence:  diffuse → project → advect → project
     * Jet forcing:    v[i, N/2] += v_SCm/1e7  for i ∈ [N/4, 3N/4]
     * Grid:           N=32, dt=0.1, visc=0.0001
     *
     * @param initial_velocity  v_SCm in m/s (default 1×10⁸)
     * @param steps             number of time steps (default 10)
     * @param verbose           print velocity field to os
     * @param os                output stream
     * @return mean velocity magnitude (normalised grid units)
     */
    double simulate_quasar_jet(double initial_velocity, int steps,
                                bool verbose, std::ostream& os);

    // ---- PAPER_370: Multi-body Pcore Planetary Scaling Law ----
    /**
     * PAPER_370 — Multi-body Solar Pcore Planetary Scaling Law
     *
     * Pcore = 1.0    for stellar bodies (Sun); full SCm core penetration
     * Pcore = 1×10⁻³ for planetary bodies (Earth, Jupiter, Neptune)
     *        → partial SCm blocking by dense solid/liquid core
     *
     * omega_c = 2π / T_orbital   (planets) — FIRST UQFF orbital-cycle bridge
     * omega_c = 2π / T_solar_cycle (Sun)
     *
     * Neptune (72K, 164.8yr) = FIRST UQFF ice giant / frozen planet module.
     *
     * Surface gravity validation:
     *   Sun=274, Earth=9.8, Jupiter=24.8, Neptune=11.2 m/s²  (all ±0.5%)
     *
     * NOTE: β_i discrepancy — thread source uses 0.6; UQFF canonical is 0.61.
     *       All pipeline classes use β_i=0.61.
     */

    // ---- UQFF 2.0 Module Class ----
    class StarMagic09SeptModule;

    // ---- Standalone self-test ----
    void run_selftest(std::ostream& os);

} // namespace StarMagic09Sept
