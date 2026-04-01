#ifndef UQFF_DM_DT_DERIVATION_H
#define UQFF_DM_DT_DERIVATION_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFDMDtDerivation
 * @brief UQFF dM/dt Full Derivation — PAPER_671 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFDMDtFromUQFFSteps (line 26116)
 *
 * Full step-by-step derivation of dM/dt from UQFF mechanics:
 *
 * Standard Hawking:
 *   dM/dt = -L_H / c²  = -hbar c⁴ / (15360 pi G² M² c²)
 *
 * UQFF Step 1 — f_TRZ negentropic reversal:
 *   (dM/dt)' = dM/dt * (1 - f_TRZ)         [suppresses evaporation ~10%]
 *
 * UQFF Step 2 — Aether density quench:
 *   (dM/dt)'' = (dM/dt)' * (rho_SCm / rho_UA)  [x0.1 further suppression]
 *
 * UQFF Step 3 — U_m string anchor:
 *   (dM/dt)_UQFF = (dM/dt)'' * exp(-U_m / k_B T_H)  [exponential suppression]
 *
 * Combined:
 *   (dM/dt)_UQFF = dM/dt * (1-f_TRZ) * (rho_SCm/rho_UA) * exp(-U_m/k_B T_H)
 *   Suppression factor: 0.9 * 0.1 * e^(-1) ~ 0.033  (30x slower evaporation)
 *
 * Time-integrated:
 *   M(t) via Euler solver; also analytic approximation for near-constant suppressor:
 *   M(t) ≈ (M0³ - 3 * |dM/dt|_UQFF / c² * t * c²)^(1/3)
 */
class UQFFDMDtDerivation {
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

    UQFFDMDtDerivation();
    double compute_dM_dt_standard(double M) const;
    double compute_dM_dt_step1(double M) const;
    double compute_dM_dt_step2(double M) const;
    double compute_dM_dt_UQFF(double M) const;
    double compute_suppression_factor(double M) const;
    double compute_M_at_t(double M0, double t) const;  // analytic approx
    void simulate_evaporation(double M0, double t0, double t1, double dt, const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
};
#endif // UQFF_DM_DT_DERIVATION_H
// UQFF helper: UQFFDMDtDerivation.cpp utils
static double _T_H(double M){
    return (UQFFDMDtDerivation::HBAR*UQFFDMDtDerivation::C*UQFFDMDtDerivation::C*UQFFDMDtDerivation::C)/
           (8.0*M_PI*UQFFDMDtDerivation::G*M*UQFFDMDtDerivation::K_B);
}
static double _L_H(double M){
    return (UQFFDMDtDerivation::HBAR*std::pow(UQFFDMDtDerivation::C,6.0))/
           (15360.0*M_PI*UQFFDMDtDerivation::G*UQFFDMDtDerivation::G*M*M);
}
static double _r_s(double M){
    return 2.0*UQFFDMDtDerivation::G*M/(UQFFDMDtDerivation::C*UQFFDMDtDerivation::C);
}
static double _tau_std(double M){
    return (5120.0*M_PI*UQFFDMDtDerivation::G*UQFFDMDtDerivation::G*M*M*M)/
           (UQFFDMDtDerivation::HBAR*std::pow(UQFFDMDtDerivation::C,4.0));
}
static double _U_m(double r, double t){
    double tn = t/UQFFDMDtDerivation::T_N_REF;
    return (UQFFDMDtDerivation::MU_J/r)*(1.0-std::exp(-UQFFDMDtDerivation::GAMMA*t*std::cos(M_PI*tn)));
}
