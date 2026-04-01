#ifndef UQFF_EVAPORATION_TIMESCALE_H
#define UQFF_EVAPORATION_TIMESCALE_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFEvaporationTimescale
 * @brief UQFF BH Evaporation Timescale — PAPER_672 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFEvaporationTimescale (line 26358)
 *
 * Standard evaporation timescale:
 *   tau_Hawking = 5120 pi G² M³ / (hbar c⁴)
 *   For M = M☉: tau ~ 2.1e74 years (much > universe age)
 *   For M = 1e11 kg:  tau ~ 3.4e-7 s  (instant)
 *   Boundary mass M_evap (tau = t_Hubble): M_evap ~ 1e11 kg
 *
 * UQFF timescale:
 *   tau_UQFF = tau_Hawking * (1/(1-f_TRZ)) * (rho_UA/rho_SCm) * exp(U_m/k_B T_H)
 *   Factor ~ 30 (same as stability proofs)
 *   UQFF boundary mass: M_evap,UQFF = M_evap * (1/factor)^(1/3) ~ M_evap / 3.1
 *
 * Universe-age crossings:
 *   Standard: M_cross,std  = (hbar c⁴ t_H / 5120 pi G²)^(1/3)
 *   UQFF:     M_cross,UQFF = M_cross,std / (factor)^(1/3)
 *
 * Sensitivity to U_m: tau_UQFF / tau_std vs U_m/k_B T_H
 */
class UQFFEvaporationTimescale {
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

    static constexpr double T_HUBBLE = 4.35e17; // s
    UQFFEvaporationTimescale();
    double compute_tau_standard(double M) const;
    double compute_tau_UQFF(double M) const;
    double compute_factor(double M) const;
    double compute_M_cross_standard() const;
    double compute_M_cross_UQFF() const;
    void sensitivity_Um_sweep(double Um_min, double Um_max, double dUm,
                              double M, const std::string& out="") const;
    void simulate_timescale_sweep(double M0, double M1, int n_pts, const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
};
#endif // UQFF_EVAPORATION_TIMESCALE_H
// UQFF helper: UQFFEvaporationTimescale.cpp utils
static double _T_H(double M){
    return (UQFFEvaporationTimescale::HBAR*UQFFEvaporationTimescale::C*UQFFEvaporationTimescale::C*UQFFEvaporationTimescale::C)/
           (8.0*M_PI*UQFFEvaporationTimescale::G*M*UQFFEvaporationTimescale::K_B);
}
static double _L_H(double M){
    return (UQFFEvaporationTimescale::HBAR*std::pow(UQFFEvaporationTimescale::C,6.0))/
           (15360.0*M_PI*UQFFEvaporationTimescale::G*UQFFEvaporationTimescale::G*M*M);
}
static double _r_s(double M){
    return 2.0*UQFFEvaporationTimescale::G*M/(UQFFEvaporationTimescale::C*UQFFEvaporationTimescale::C);
}
static double _tau_std(double M){
    return (5120.0*M_PI*UQFFEvaporationTimescale::G*UQFFEvaporationTimescale::G*M*M*M)/
           (UQFFEvaporationTimescale::HBAR*std::pow(UQFFEvaporationTimescale::C,4.0));
}
static double _U_m(double r, double t){
    double tn = t/UQFFEvaporationTimescale::T_N_REF;
    return (UQFFEvaporationTimescale::MU_J/r)*(1.0-std::exp(-UQFFEvaporationTimescale::GAMMA*t*std::cos(M_PI*tn)));
}
