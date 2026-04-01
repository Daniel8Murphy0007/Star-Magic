#ifndef UQFF_HAWKING_TEMPERATURE_MODULATION_H
#define UQFF_HAWKING_TEMPERATURE_MODULATION_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFHawkingTemperatureModulation
 * @brief UQFF Hawking Temperature Modulation — PAPER_683 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~21788
 *
 * Standard Hawking temperature: T_H = hbar c^3 / (8 pi G M k_B)
 *
 * UQFF modulates T_H through three channels:
 *   1. f_TRZ time-reversal boost:     T_UQFF = T_H * (1 + f_TRZ)
 *   2. Density ratio suppression:     T_UQFF = T_H * (1 - rho_SCm/rho_UA)
 *   3. U_m magnetic string modulation:
 *      T_UQFF = T_H * (1 + f_TRZ) * (1 - rho_SCm/rho_UA) * (1 + U_m/k_B T_H)
 *
 *   Combined (first-order expansion):
 *   T_UQFF = T_H * [1 + f_TRZ + U_m/k_B T_H - rho_SCm/rho_UA * (1 + f_TRZ)]
 *
 * Spectral modulation:
 *   n(omega)_UQFF = 1 / (exp(hbar omega / k_B T_UQFF) - 1)
 *   -> shifted Planck spectrum peak wavelength
 *   lambda_max_UQFF = hbar c / (2.82 k_B T_UQFF)  (Wien law)
 *
 * Time-dependence: T_UQFF(t) includes time-varying U_m(r,t) term.
 */
class UQFFHawkingTemperatureModulation {
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

    UQFFHawkingTemperatureModulation();
    double compute_T_H(double M) const;
    double compute_T_UQFF(double M, double r=0.0, double t=1e8) const;
    double compute_planck_spectrum(double omega, double T) const;
    double compute_wien_peak(double T) const;
    double compute_T_modulation_factor(double M, double r, double t) const;
    void simulate_T_vs_M(double M_start, double M_end, double dM,
                         const std::string& out="") const;
    void simulate_spectrum(double M, double omega_min, double omega_max,
                           double d_omega, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // UQFF_HAWKING_TEMPERATURE_MODULATION_H
