#ifndef UQFF_BLACK_HOLE_INVERSION_H
#define UQFF_BLACK_HOLE_INVERSION_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFBlackHoleInversion
 * @brief UQFF Black Hole Inversion Probability — PAPER_663 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFBlackHoleInversion (line 24810)
 *
 * UQFF BH interior inversion via [UA]/[SCm] gradient flip.
 * Derives full stochastic Theta_inv criterion (log-normal distribution).
 *
 * Core equations:
 *   r_s,UQFF = r_s * (1 - delta_rho)   delta_rho = rho_SCm/rho_UA ~ 0.1
 *   E_inv,UQFF = G M² / r_s,UQFF
 *   P_inv = f_TRZ * exp(-E_inv / k_B T_H)
 *   Phi_inv = (1/delta_rho) * (G M / c) * (1 + f_TRZ)   [gradient drives flux]
 *   S_Um = exp(U_m / k_B T_H)
 *   Theta_inv = P_inv * Phi_inv * S_Um
 *   If Theta_inv > 1 → inversion occurs
 *
 * Stochastic distribution:
 *   delta_rho, f_TRZ, U_m sampled from Gaussians  → Theta_inv log-normal
 *   P_invert = Prob(Theta_inv > 1)  (Monte-Carlo or analytic CDF)
 *
 * Numerical (Sgr A*): P_invert ~ 0.95 with Gaussian variability
 */
class UQFFBlackHoleInversion {
public:
    static constexpr double G    = 6.6743e-11;
    static constexpr double C    = 2.998e8;
    static constexpr double HBAR = 1.0546e-34;
    static constexpr double K_B  = 1.380649e-23;
    static constexpr double PI   = 3.14159265358979323846;
    static constexpr double RHO_UA  = 7.09e-36;   // J/m3 [UA] vacuum density
    static constexpr double RHO_SCM = 7.09e-37;   // J/m3 [SCm] vacuum density
    static constexpr double F_TRZ   = 0.1;         // time-reversal factor
    static constexpr double KAPPA   = 0.0005;      // kappa day-1
    static constexpr double SSQ     = 0.57;        // [SSq]
    static constexpr double MU_J    = 3.38e23;     // J/T magnetic string j=1
    static constexpr double GAMMA   = 5.0e-5/86400.0; // s-1
    static constexpr double T_N_REF = 1.0e8;       // s normalisation
    static constexpr double M_SUN   = 1.989e30;    // kg

    UQFFBlackHoleInversion(unsigned int seed=42u);
    double compute_r_s_UQFF(double M) const;
    double compute_E_inv(double M) const;
    double compute_P_inv(double M) const;
    double compute_Phi_inv(double M) const;
    double compute_Theta_inv(double M) const;
    double compute_P_invert_MC(double M, int n=5000, double sigma=0.05) const;
    void simulate_P_invert_sweep(double M0, double M1, double dM, const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    mutable std::mt19937 rng_;
    std::vector<std::function<double(double)>> mods_;
};
#endif // UQFF_BLACK_HOLE_INVERSION_H
// UQFF helper: UQFFBlackHoleInversion.cpp utils
static double _T_H(double M){
    return (UQFFBlackHoleInversion::HBAR*UQFFBlackHoleInversion::C*UQFFBlackHoleInversion::C*UQFFBlackHoleInversion::C)/
           (8.0*M_PI*UQFFBlackHoleInversion::G*M*UQFFBlackHoleInversion::K_B);
}
static double _L_H(double M){
    return (UQFFBlackHoleInversion::HBAR*std::pow(UQFFBlackHoleInversion::C,6.0))/
           (15360.0*M_PI*UQFFBlackHoleInversion::G*UQFFBlackHoleInversion::G*M*M);
}
static double _r_s(double M){
    return 2.0*UQFFBlackHoleInversion::G*M/(UQFFBlackHoleInversion::C*UQFFBlackHoleInversion::C);
}
static double _tau_std(double M){
    return (5120.0*M_PI*UQFFBlackHoleInversion::G*UQFFBlackHoleInversion::G*M*M*M)/
           (UQFFBlackHoleInversion::HBAR*std::pow(UQFFBlackHoleInversion::C,4.0));
}
static double _U_m(double r, double t){
    double tn = t/UQFFBlackHoleInversion::T_N_REF;
    return (UQFFBlackHoleInversion::MU_J/r)*(1.0-std::exp(-UQFFBlackHoleInversion::GAMMA*t*std::cos(M_PI*tn)));
}
