#ifndef UQFF_STABILITY_NUMERICALLY_FOR_SgrA_H
#define UQFF_STABILITY_NUMERICALLY_FOR_SgrA_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFStabilityNumericallyForSgrA
 * @brief UQFF Numerical Stability Analysis for Sgr A* — PAPER_682 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~21554
 *
 * Sgr A* parameters:
 *   M_SgrA = 4.297e6 Msun = 8.548e36 kg
 *   Distance: 8.178 kpc = 2.523e20 m
 *
 * Numerical stability analysis in UQFF:
 *   1. Perturbation expansion: delta_M(t) / M_0 = epsilon * exp(i omega t)
 *      omega = omega_R + i omega_I
 *      stable: omega_I < 0  (damped perturbation)
 *
 *   2. UQFF stability equation:
 *      omega_UQFF = omega_GR * (1 + f_TRZ * rho_UA/rho_SCm * U_m/k_B T_H)
 *      -> enhanced damping rate
 *
 *   3. Lyapunov exponent:
 *      lambda_UQFF = lambda_GR * (RHO_SCM/RHO_UA) * exp(-U_m/k_B T_H)
 *      lambda < 0: stable fixed point
 *
 *   4. RK4 integration of perturbed BH mass evolution:
 *      dM/dt = -(1-f_TRZ) * L_H(M) / c^2 * (rho_SCm/rho_UA) * exp(-U_m/k_B T_H)
 */
class UQFFStabilityNumericallyForSgrA {
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

    static constexpr double M_SGRA   = 4.297e6 * 1.989e30;
    static constexpr double D_SGRA   = 2.523e20;

    UQFFStabilityNumericallyForSgrA();
    double compute_omega_I_UQFF(double M) const;
    double compute_lyapunov_UQFF(double M) const;
    double compute_dM_dt_UQFF(double M, double t) const;
    std::vector<std::pair<double,double>> rk4_evolution(double M0, double t_end,
                                                        double dt) const;
    bool is_stable(double M) const;
    void simulate_stability(double t_end, double dt, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // UQFF_STABILITY_NUMERICALLY_FOR_SgrA_H
