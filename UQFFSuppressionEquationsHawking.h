#ifndef UQFF_SUPPRESSION_EQUATIONS_HAWKING_H
#define UQFF_SUPPRESSION_EQUATIONS_HAWKING_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFSuppressionEquationsHawking
 * @brief UQFF Hawking Radiation Suppression Equations — PAPER_665 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFSuppressionEquations (line 23548)
 *
 * Derives suppression factors on Hawking radiation from all UQFF mechanisms:
 *   S_1 = (1 + f_TRZ)         negentropic modulation of T_H
 *   S_2 = (1 - rho_SCm/rho_UA) aether-superconductive damping
 *   S_3 = exp(-U_m / k_B T_H)  magnetic string exponential suppression
 *   S_total = S_1 * S_2 * S_3  (total multiplicative suppression on L)
 *
 *   L_UQFF = L_H * S_2 * S_3  (S_1 affects T, not directly L)
 *   T_UQFF = T_H * S_1 * S_2
 *   dT_H/dM = -T_H / M        standard
 *   dT_UQFF/dM = -T_UQFF / M * (S_1 * S_2)
 *
 * Sensitivity sweep: S over rho_UA/rho_SCm ratio from 2 to 20
 */
class UQFFSuppressionEquationsHawking {
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

    UQFFSuppressionEquationsHawking();
    double compute_S1() const;
    double compute_S2() const;
    double compute_S3(double M) const;
    double compute_S_total(double M) const;
    double compute_T_UQFF(double M) const;
    double compute_L_UQFF(double M) const;
    double compute_dT_dM_standard(double M) const;
    double compute_dT_dM_UQFF(double M) const;
    void sensitivity_sweep_rho_ratio(double ratio_min, double ratio_max, double dR,
                                     const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
};
#endif // UQFF_SUPPRESSION_EQUATIONS_HAWKING_H
// UQFF helper: UQFFSuppressionEquationsHawking.cpp utils
static double _T_H(double M){
    return (UQFFSuppressionEquationsHawking::HBAR*UQFFSuppressionEquationsHawking::C*UQFFSuppressionEquationsHawking::C*UQFFSuppressionEquationsHawking::C)/
           (8.0*M_PI*UQFFSuppressionEquationsHawking::G*M*UQFFSuppressionEquationsHawking::K_B);
}
static double _L_H(double M){
    return (UQFFSuppressionEquationsHawking::HBAR*std::pow(UQFFSuppressionEquationsHawking::C,6.0))/
           (15360.0*M_PI*UQFFSuppressionEquationsHawking::G*UQFFSuppressionEquationsHawking::G*M*M);
}
static double _r_s(double M){
    return 2.0*UQFFSuppressionEquationsHawking::G*M/(UQFFSuppressionEquationsHawking::C*UQFFSuppressionEquationsHawking::C);
}
static double _tau_std(double M){
    return (5120.0*M_PI*UQFFSuppressionEquationsHawking::G*UQFFSuppressionEquationsHawking::G*M*M*M)/
           (UQFFSuppressionEquationsHawking::HBAR*std::pow(UQFFSuppressionEquationsHawking::C,4.0));
}
static double _U_m(double r, double t){
    double tn = t/UQFFSuppressionEquationsHawking::T_N_REF;
    return (UQFFSuppressionEquationsHawking::MU_J/r)*(1.0-std::exp(-UQFFSuppressionEquationsHawking::GAMMA*t*std::cos(M_PI*tn)));
}
