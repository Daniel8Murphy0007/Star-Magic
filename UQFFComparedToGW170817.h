#ifndef UQFF_COMPARED_TO_GW170817_H
#define UQFF_COMPARED_TO_GW170817_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFComparedToGW170817
 * @brief UQFF vs LIGO/Virgo GW170817 (NS-NS merger) — PAPER_675 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~19205
 *
 * GW170817: Binary neutron star coalescence, Aug 17, 2017.
 *   m1 = 1.36 Msun, m2 = 1.17 Msun
 *   Luminosity distance: ~40 Mpc = 1.234e24 m
 *   Peak frequency: ~1500 Hz
 *   Multi-messenger: GRB170817A detected 1.7 s after merger
 *
 * UQFF modifications:
 *   T_H,NS = HBAR c^3 / (8 pi G M_NS k_B)   (characteristic NS Hawking T)
 *   h_UQFF = h_GR * (1-f_TRZ) * S_SCm * S_Um
 *   GRB delay: Delta_t_UQFF = 1.7 s * (1 + f_TRZ * rho_UA/rho_SCm)
 *              -> probes aether propagation speed
 *   tidal deformability Lambda modified by rho_UA suppression
 *   kappa_tidal_UQFF = kappa_tidal_GR * (1 - f_TRZ * rho_SCm/rho_UA)
 */
class UQFFComparedToGW170817 {
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

    static constexpr double M1_GW170817 = 1.36 * 1.989e30;
    static constexpr double M2_GW170817 = 1.17 * 1.989e30;
    static constexpr double D_GW170817  = 1.234e24;  // 40 Mpc in m
    static constexpr double F_PEAK      = 1500.0;    // Hz
    static constexpr double DT_GRB      = 1.7;       // s

    UQFFComparedToGW170817();
    double compute_chirp_mass(double m1, double m2) const;
    double compute_h_GR(double Mc, double d, double f) const;
    double compute_h_UQFF(double m1, double m2, double d, double f,
                          double t=1e8) const;
    double compute_grb_delay_UQFF() const;
    double compute_tidal_UQFF(double kappa_GR) const;
    void simulate_inspiral(double f_start, double f_peak, double df,
                           const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // UQFF_COMPARED_TO_GW170817_H
