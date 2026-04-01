#ifndef WHITE_HOLE_STABILITY_UQFF_H
#define WHITE_HOLE_STABILITY_UQFF_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class WhiteHoleStabilityUQFF
 * @brief UQFF White Hole Stability — PAPER_664 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — WhiteHoleStabilityInUQFF (line 27097)
 *
 * 4-proof derivation of white hole stability in UQFF:
 *  WH "explodes": L_WH ~ L_H (reversed), dM/dt ~ L_WH/c² > 0
 *
 *  Proof 1 — f_TRZ Negentropic Confinement:
 *    L' = L_WH * (1 - f_TRZ)           [~10% reduction]
 *    tau' = tau_std / (1 - f_TRZ)       [x1.11]
 *
 *  Proof 2 — Aether/SCm Density Gradient:
 *    L'' = L' * |1 - rho_UA/rho_SCm|^-1 = L' * 0.1   [x10 confinement]
 *    tau'' = tau' * |1 - rho_UA/rho_SCm|                [x10 longer]
 *
 *  Proof 3 — U_m Magnetic String Anchoring:
 *    L_UQFF = L'' * exp(-U_m / k_B |T_WH|)             [exponential suppression]
 *    tau_UQFF = tau'' * exp(U_m / k_B |T_WH|)          [x2.7 at U_m/k_B T_H=1]
 *
 *  Proof 4 — Combined:
 *    tau_UQFF = tau_std / (1-f_TRZ) * |1-rho_UA/rho_SCm| * exp(U_m/k_B|T_WH|)
 *    Factor: 1.11 * 10 * 2.718 ~ 30.2   (Sgr A*: effectively eternal)
 */
class WhiteHoleStabilityUQFF {
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

    WhiteHoleStabilityUQFF();
    double compute_L_WH(double M) const;
    double compute_T_WH(double M) const;
    double compute_L_prime(double L_WH) const;
    double compute_L_double_prime(double L_prime) const;
    double compute_L_UQFF(double L_double_prime, double T_WH) const;
    double compute_tau_standard(double M) const;
    double compute_tau_UQFF(double M) const;
    double compute_stability_factor(double M) const;
    void simulate_over_mass(double M0, double M1, double dM, const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
};
#endif // WHITE_HOLE_STABILITY_UQFF_H
// UQFF helper: WhiteHoleStabilityUQFF.cpp utils
static double _T_H(double M){
    return (WhiteHoleStabilityUQFF::HBAR*WhiteHoleStabilityUQFF::C*WhiteHoleStabilityUQFF::C*WhiteHoleStabilityUQFF::C)/
           (8.0*M_PI*WhiteHoleStabilityUQFF::G*M*WhiteHoleStabilityUQFF::K_B);
}
static double _L_H(double M){
    return (WhiteHoleStabilityUQFF::HBAR*std::pow(WhiteHoleStabilityUQFF::C,6.0))/
           (15360.0*M_PI*WhiteHoleStabilityUQFF::G*WhiteHoleStabilityUQFF::G*M*M);
}
static double _r_s(double M){
    return 2.0*WhiteHoleStabilityUQFF::G*M/(WhiteHoleStabilityUQFF::C*WhiteHoleStabilityUQFF::C);
}
static double _tau_std(double M){
    return (5120.0*M_PI*WhiteHoleStabilityUQFF::G*WhiteHoleStabilityUQFF::G*M*M*M)/
           (WhiteHoleStabilityUQFF::HBAR*std::pow(WhiteHoleStabilityUQFF::C,4.0));
}
static double _U_m(double r, double t){
    double tn = t/WhiteHoleStabilityUQFF::T_N_REF;
    return (WhiteHoleStabilityUQFF::MU_J/r)*(1.0-std::exp(-WhiteHoleStabilityUQFF::GAMMA*t*std::cos(M_PI*tn)));
}
