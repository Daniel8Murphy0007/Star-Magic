#ifndef UQFF_PRIMORDIAL_BH_EVAPORATION_H
#define UQFF_PRIMORDIAL_BH_EVAPORATION_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFPrimordialBHEvaporation
 * @brief UQFF Primordial Black Hole Evaporation — PAPER_684 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~22011
 *
 * Primordial BHs (PBHs) form in the early universe from density
 * fluctuations. Standard Hawking evaporation:
 *   dM/dt = -hbar c^4 / (15360 pi G^2 M^2)  [Schwarzschild PBH]
 *
 * UQFF modifies evaporation amplitude by suppression chain:
 *   dM/dt_UQFF = dM/dt * (1-f_TRZ) * (rho_SCm/rho_UA) * exp(-U_m/k_B T_H)
 *   -> factor ~0.033 slower than standard
 *
 * PBH mass at formation: M_f = rho_rad(t_f) * (4/3) pi (c t_f / 2)^3
 *   For t_f ~ 1e-23 s: M_f ~ 1e10 kg  (minimum viable DM mass in UQFF)
 *
 * Analytic lifetime: tau_UQFF = 5120 pi G^2 M_f^3 / (hbar c^4)
 *                              * 1/(1-f_TRZ) * (rho_UA/rho_SCm) * exp(U_m/k_B T_H)
 *
 * Euler / RK4 numerical integration of M(t) from M_f to 0.
 * Burst parameters: dN_gamma/dE at end of life — UQFF delays burst epoch.
 */
class UQFFPrimordialBHEvaporation {
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

    UQFFPrimordialBHEvaporation();
    double compute_dM_dt_std(double M) const;
    double compute_dM_dt_UQFF(double M, double t) const;
    double compute_tau_std(double M) const;
    double compute_tau_UQFF(double M) const;
    double compute_M_formation(double t_form) const;
    std::vector<std::pair<double,double>> evolve_rk4(double M0, double t_end,
                                                     double dt) const;
    void simulate_evaporation(double M0, double dt, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // UQFF_PRIMORDIAL_BH_EVAPORATION_H
