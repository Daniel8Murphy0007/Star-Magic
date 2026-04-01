#ifndef WHITE_HOLE_RADIATION_UQFF_H
#define WHITE_HOLE_RADIATION_UQFF_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class WhiteHoleRadiationUQFF
 * @brief UQFF White Hole Radiation — PAPER_660 | Session 172 | April 2, 2026
 *
 * Source: grok_share_fc21e30c24b4.txt — WhiteHoleRadiation class (line 26585)
 *
 * 6-step derivation (time-reversed Hawking in UQFF):
 *  Step 1: L_WH ~ -L_H  (time-reversed expulsion, magnitude = L_H)
 *  Step 2: UQFF inversion: r_s,UQFF = r_s(1-rho_SCm/rho_UA)
 *  Step 3: f_TRZ boost: L_WH' = L_H * (1+f_TRZ)
 *  Step 4: Aether amplification: L_WH'' = L_WH' * (rho_UA/rho_SCm)  ~x10
 *  Step 5: U_m channeling: L_WH,UQFF = L_WH'' * exp(U_m / k_B T_H)
 *  Step 6: Full formula:
 *    L_WH,UQFF = (hbar c^6)/(15360 pi G^2 M^2)*(1+f_TRZ)*(rho_UA/rho_SCm)*exp(U_m/k_B T_H)
 * Numerical (Sgr A*): L_WH,UQFF ~ 3e-3 W  (vs L_H ~1e-29 W)
 */
class WhiteHoleRadiationUQFF {
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

    WhiteHoleRadiationUQFF();
    double compute_L_H(double M) const;
    double compute_T_H(double M) const;
    double compute_r_s_UQFF(double M) const;
    double compute_L_WH_prime(double L_H) const;
    double compute_L_WH_double_prime(double L_WH_prime) const;
    double compute_L_WH_UQFF(double M, double r=0.0, double t=0.0) const;
    void simulate_over_M(double M_start, double M_end, double dM, const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
    mutable std::mt19937 rng_{42u};
};
#endif // WHITE_HOLE_RADIATION_UQFF_H
// UQFF helper: WhiteHoleRadiationUQFF.cpp utils
static double _T_H(double M){
    return (WhiteHoleRadiationUQFF::HBAR*WhiteHoleRadiationUQFF::C*WhiteHoleRadiationUQFF::C*WhiteHoleRadiationUQFF::C)/
           (8.0*M_PI*WhiteHoleRadiationUQFF::G*M*WhiteHoleRadiationUQFF::K_B);
}
static double _L_H(double M){
    return (WhiteHoleRadiationUQFF::HBAR*std::pow(WhiteHoleRadiationUQFF::C,6.0))/
           (15360.0*M_PI*WhiteHoleRadiationUQFF::G*WhiteHoleRadiationUQFF::G*M*M);
}
static double _r_s(double M){
    return 2.0*WhiteHoleRadiationUQFF::G*M/(WhiteHoleRadiationUQFF::C*WhiteHoleRadiationUQFF::C);
}
static double _tau_std(double M){
    return (5120.0*M_PI*WhiteHoleRadiationUQFF::G*WhiteHoleRadiationUQFF::G*M*M*M)/
           (WhiteHoleRadiationUQFF::HBAR*std::pow(WhiteHoleRadiationUQFF::C,4.0));
}
static double _U_m(double r, double t){
    double tn = t/WhiteHoleRadiationUQFF::T_N_REF;
    return (WhiteHoleRadiationUQFF::MU_J/r)*(1.0-std::exp(-WhiteHoleRadiationUQFF::GAMMA*t*std::cos(M_PI*tn)));
}
