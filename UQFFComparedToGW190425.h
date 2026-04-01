#ifndef UQFF_COMPARED_TO_GW190425_H
#define UQFF_COMPARED_TO_GW190425_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFComparedToGW190425
 * @brief UQFF vs LIGO/Virgo GW190425 (heavy NS-NS) — PAPER_676 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~19417
 *
 * GW190425: Binary NS merger April 25, 2019. Heaviest NS binary known.
 *   Total mass: 3.4 Msun (m1~1.9, m2~1.5 Msun)
 *   Distance: ~159 Mpc = 4.9e24 m
 *   Only LIGO-Livingston + Virgo (LIGO-Hanford offline)
 *   No EM counterpart detected (high distance, low inclination sensitivity)
 *
 * UQFF analysis:
 *   Mc_GW190425 = (1.9*1.5*Msun^2)^0.6 / (3.4*Msun)^0.2
 *   h_UQFF(f) = h_GR(f) * (1-f_TRZ) * S_SCm * S_Um
 *   Mass increase over GW170817 (~2.5x higher total mass)
 *   -> UQFF suppression factors scale with T_H and r_s
 *   post-merger frequency f_pm ~ 2.5 kHz: UQFF phase shift probes aether
 *   UQFF upper-limit on ejecta mass through U_m suppression
 */
class UQFFComparedToGW190425 {
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

    static constexpr double M1_DEFAULT = 1.9 * 1.989e30;
    static constexpr double M2_DEFAULT = 1.5 * 1.989e30;
    static constexpr double D_DEFAULT  = 4.9e24;
    static constexpr double F_PM       = 2500.0;

    UQFFComparedToGW190425();
    double compute_chirp_mass(double m1, double m2) const;
    double compute_h_GR(double Mc, double d, double f) const;
    double compute_h_UQFF(double m1, double m2, double d, double f,
                          double t=1e8) const;
    double compute_post_merger_phase(double m1, double m2) const;
    double compute_ejecta_limit_UQFF(double m1, double m2) const;
    void simulate_inspiral(double f_start, double f_end, double df,
                           const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // UQFF_COMPARED_TO_GW190425_H
