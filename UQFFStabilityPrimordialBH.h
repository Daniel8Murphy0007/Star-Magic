#ifndef UQFF_STABILITY_PRIMORDIAL_BH_H
#define UQFF_STABILITY_PRIMORDIAL_BH_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFStabilityPrimordialBH
 * @brief UQFF Primordial BH Stability Analysis — PAPER_668 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFStabilityForPrimordialBHs (line 25798)
 *
 * PBH stability analysis in UQFF context:
 *  Standard: τ_std = 5120 π G² M³ / (ħ c⁴)
 *  M < 1e12 kg → τ < universe age → evaporate
 *  M ~ 1e10-1e15 g constrained by gamma-ray observations
 *
 * UQFF step-by-step:
 *   1. τ_std (Hawking)
 *   2. τ' = τ / (1 - f_TRZ)      x1.11  negentropic
 *   3. τ'' = τ' * (rho_UA/rho_SCm)  x10
 *   4. τ_UQFF = τ'' * exp(U_m/k_B T_H) x2.7
 *   Total: ~30x
 *
 * Numerical check M=1e12 kg:
 *   τ_std ~ 1e10 yr; τ_UQFF ~ 3e11 yr > universe age → UQFF promotes DM candidate
 *
 * Mass categories: stable (τ_UQFF > t_H), marginal (0.1 t_H < τ < t_H), evaporating
 * t_H = 4.35e17 s (Hubble time)
 */
class UQFFStabilityPrimordialBH {
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
    UQFFStabilityPrimordialBH();
    double compute_tau_std(double M) const;
    double compute_tau_UQFF(double M) const;
    std::string classify(double M) const;    // stable/marginal/evaporating
    double pbh_dm_window_min_mass_uqff() const;  // minimum DM-viable mass
    void simulate_mass_stability_map(double M0, double M1, int n_pts, const std::string& out="") const;
    void add_modifier(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
};
#endif // UQFF_STABILITY_PRIMORDIAL_BH_H
// UQFF helper: UQFFStabilityPrimordialBH.cpp utils
static double _T_H(double M){
    return (UQFFStabilityPrimordialBH::HBAR*UQFFStabilityPrimordialBH::C*UQFFStabilityPrimordialBH::C*UQFFStabilityPrimordialBH::C)/
           (8.0*M_PI*UQFFStabilityPrimordialBH::G*M*UQFFStabilityPrimordialBH::K_B);
}
static double _L_H(double M){
    return (UQFFStabilityPrimordialBH::HBAR*std::pow(UQFFStabilityPrimordialBH::C,6.0))/
           (15360.0*M_PI*UQFFStabilityPrimordialBH::G*UQFFStabilityPrimordialBH::G*M*M);
}
static double _r_s(double M){
    return 2.0*UQFFStabilityPrimordialBH::G*M/(UQFFStabilityPrimordialBH::C*UQFFStabilityPrimordialBH::C);
}
static double _tau_std(double M){
    return (5120.0*M_PI*UQFFStabilityPrimordialBH::G*UQFFStabilityPrimordialBH::G*M*M*M)/
           (UQFFStabilityPrimordialBH::HBAR*std::pow(UQFFStabilityPrimordialBH::C,4.0));
}
static double _U_m(double r, double t){
    double tn = t/UQFFStabilityPrimordialBH::T_N_REF;
    return (UQFFStabilityPrimordialBH::MU_J/r)*(1.0-std::exp(-UQFFStabilityPrimordialBH::GAMMA*t*std::cos(M_PI*tn)));
}
