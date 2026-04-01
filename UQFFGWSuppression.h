#ifndef UQFF_GW_SUPPRESSION_H
#define UQFF_GW_SUPPRESSION_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFGWSuppression
 * @brief UQFF Gravitational Wave Power Suppression — PAPER_666 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFGWSuppression (line 23812)
 *
 * Standard GW power (Peters formula for circular orbit):
 *   P_GW = 32/5 * G^4/c^5 * m1^2 * m2^2 * (m1+m2) / r^5
 *
 * UQFF suppression mechanisms:
 *   S_UA  = 1 - (rho_UA / rho_critical)    [aether absorption]
 *   S_SCm = exp(-rho_SCm * r_s / (k_B T_H)) [superconductive damping]
 *   S_TRZ = (1 - f_TRZ)                    [negentropic reversal]
 *   S_Um  = exp(-U_m / (omega_GW * c))      [string impedance; omega_GW = c/r_s]
 *
 *   P_GW_UQFF = P_GW * S_UA * S_SCm * S_TRZ * S_Um
 *
 * Strain suppression: h_UQFF = h_standard * sqrt(P_GW_UQFF / P_GW)
 * Validation: GW150914 — compare UQFF strain vs LIGO observed
 */
class UQFFGWSuppression {
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

    static constexpr double RHO_CRIT = 9.47e-27; // kg/m3 cosmological critical density
    UQFFGWSuppression();
    double compute_P_GW_standard(double m1, double m2, double r) const;
    double compute_S_UA() const;
    double compute_S_SCm(double M) const;
    double compute_S_TRZ() const;
    double compute_S_Um(double M) const;
    double compute_P_GW_UQFF(double m1, double m2, double r) const;
    double compute_h_suppression(double m1, double m2, double r) const;
    void simulate_P_sweep(double r_min, double r_max, double dr,
                          double m1, double m2, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // UQFF_GW_SUPPRESSION_H
// UQFF helper: UQFFGWSuppression.cpp utils
static double _T_H(double M){
    return (UQFFGWSuppression::HBAR*UQFFGWSuppression::C*UQFFGWSuppression::C*UQFFGWSuppression::C)/
           (8.0*M_PI*UQFFGWSuppression::G*M*UQFFGWSuppression::K_B);
}
static double _L_H(double M){
    return (UQFFGWSuppression::HBAR*std::pow(UQFFGWSuppression::C,6.0))/
           (15360.0*M_PI*UQFFGWSuppression::G*UQFFGWSuppression::G*M*M);
}
static double _r_s(double M){
    return 2.0*UQFFGWSuppression::G*M/(UQFFGWSuppression::C*UQFFGWSuppression::C);
}
static double _tau_std(double M){
    return (5120.0*M_PI*UQFFGWSuppression::G*UQFFGWSuppression::G*M*M*M)/
           (UQFFGWSuppression::HBAR*std::pow(UQFFGWSuppression::C,4.0));
}
static double _U_m(double r, double t){
    double tn = t/UQFFGWSuppression::T_N_REF;
    return (UQFFGWSuppression::MU_J/r)*(1.0-std::exp(-UQFFGWSuppression::GAMMA*t*std::cos(M_PI*tn)));
}
