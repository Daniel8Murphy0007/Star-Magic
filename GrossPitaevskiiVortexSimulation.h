#ifndef GROSS_PITAEVSKII_VORTEX_SIMULATION_H
#define GROSS_PITAEVSKII_VORTEX_SIMULATION_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class GrossPitaevskiiVortexSimulation
 * @brief UQFF Gross-Pitaevskii Vortex Simulation — PAPER_681 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~21332
 *
 * Full 1D Gross-Pitaevskii equation for the [UA] Aether wavefunction
 * in radial (spherical) coordinates around a UQFF black hole:
 *
 *   i hbar dpsi/dt = [ -hbar^2/(2 m_UA) (1/r^2 d/dr r^2 d/dr)
 *                    + V_grav(r)
 *                    + g_UA |psi|^2
 *                    + U_m(r,t) ] psi
 *
 * where V_grav(r) = -G M m_UA / r  (gravitational well)
 *       U_m(r,t) from UQFF magnetic string term
 *
 * Imaginary-time propagation to find ground state.
 * Real-time propagation to simulate vortex dynamics.
 * Normalisation: integral |psi|^2 r^2 dr = N_UA * V_sphere
 *
 * Output: density profile |psi(r)|^2, phase phi(r), current j(r).
 */
class GrossPitaevskiiVortexSimulation {
public:
    static constexpr double G       = 6.6743e-11;
    static constexpr double C       = 2.998e8;
    static constexpr double HBAR    = 1.0546e-34;
    static constexpr double K_B     = 1.380649e-23;
    static constexpr double PI      = 3.14159265358979323846;
    static constexpr double M_SUN   = 1.989e30;
    static constexpr double RHO_UA  = 7.09e-36;
    static constexpr double RHO_SCM = 7.09e-37;
    static constexpr double F_TRZ   = 0.1;
    static constexpr double KAPPA   = 0.0005;
    static constexpr double SSQ     = 0.57;
    static constexpr double MU_J    = 3.38e23;
    static constexpr double GAMMA   = 5.0e-5 / 86400.0;
    static constexpr double T_N_REF = 1.0e8;

    static constexpr double H0_SI = 2.184e-18;
    static constexpr double M_UA  = HBAR * H0_SI / (C * C);
    static constexpr double G_UA  = 1.0e-10;
    static constexpr double N_UA  = RHO_UA / (HBAR * H0_SI / (C * C));

    GrossPitaevskiiVortexSimulation();
    // Compute single GP time step (imaginary time)
    std::vector<double> imaginary_time_step(
        const std::vector<double>& psi,
        const std::vector<double>& r_grid,
        double M_bh, double dt, int n_steps) const;
    double compute_chemical_potential(const std::vector<double>& psi,
                                      const std::vector<double>& r_grid,
                                      double M_bh) const;
    void simulate(double r_min, double r_max, int nr, double M_bh,
                  int n_steps, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // GROSS_PITAEVSKII_VORTEX_SIMULATION_H
