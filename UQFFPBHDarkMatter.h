#ifndef UQFF_PBH_DARK_MATTER_H
#define UQFF_PBH_DARK_MATTER_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFPBHDarkMatter
 * @brief UQFF Primordial Black Hole Dark Matter — PAPER_661 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFPBHDarkMatterImplications (line 22222)
 * Audit: SESSION_172_AUDIT_HELPER.md — PBH Dark Matter section
 *
 * Standard Hawking evaporation: tau_std = 5120 pi G² M³ / (hbar c⁴)
 * PBH with M < 1e12 kg: tau < age of universe → evaporate → gamma rays → constrain DM
 * M ~ 1e10-1e15 g: constrained window
 *
 * UQFF elevates lifetime by ~11x (same factor as LQC critical density):
 *   Step 1: tau_std (Hawking)
 *   Step 2: tau' = tau_std / (1 - f_TRZ)          x1.11  (negentropic)
 *   Step 3: tau'' = tau' * (rho_UA/rho_SCm)        x10    (aether suppression)
 *   Step 4: tau_UQFF = tau'' * exp(U_m / k_B T_H)  x2.7   (string anchor)
 *   Total factor: ~30x (PBHs in "evaporation window" become viable DM)
 *
 *   f_PBH_UQFF = Omega_PBH / Omega_DM  (enhanced by tau factor)
 *   UQFF reopens DM window: M ~ 1e10-1e15 g now stable
 */
class UQFFPBHDarkMatter {
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

    UQFFPBHDarkMatter();
    double compute_tau_standard(double M) const;
    double compute_tau_prime(double tau_std) const;
    double compute_tau_double_prime(double tau_prime) const;
    double compute_tau_UQFF(double M) const;
    double compute_f_PBH_enhancement(double M) const;
    double compute_T_H(double M) const;
    bool is_DM_candidate(double M) const;
    void simulate_lifetime_sweep(double M_start, double M_end, double dM, const std::string& out="") const;
    void add_modifier(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
};
#endif // UQFF_PBH_DARK_MATTER_H
// UQFF helper: UQFFPBHDarkMatter.cpp utils
static double _T_H(double M){
    return (UQFFPBHDarkMatter::HBAR*UQFFPBHDarkMatter::C*UQFFPBHDarkMatter::C*UQFFPBHDarkMatter::C)/
           (8.0*M_PI*UQFFPBHDarkMatter::G*M*UQFFPBHDarkMatter::K_B);
}
static double _L_H(double M){
    return (UQFFPBHDarkMatter::HBAR*std::pow(UQFFPBHDarkMatter::C,6.0))/
           (15360.0*M_PI*UQFFPBHDarkMatter::G*UQFFPBHDarkMatter::G*M*M);
}
static double _r_s(double M){
    return 2.0*UQFFPBHDarkMatter::G*M/(UQFFPBHDarkMatter::C*UQFFPBHDarkMatter::C);
}
static double _tau_std(double M){
    return (5120.0*M_PI*UQFFPBHDarkMatter::G*UQFFPBHDarkMatter::G*M*M*M)/
           (UQFFPBHDarkMatter::HBAR*std::pow(UQFFPBHDarkMatter::C,4.0));
}
static double _U_m(double r, double t){
    double tn = t/UQFFPBHDarkMatter::T_N_REF;
    return (UQFFPBHDarkMatter::MU_J/r)*(1.0-std::exp(-UQFFPBHDarkMatter::GAMMA*t*std::cos(M_PI*tn)));
}
