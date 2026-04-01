#ifndef UQFF_COMPARED_TO_GW150914_H
#define UQFF_COMPARED_TO_GW150914_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFComparedToGW150914
 * @brief UQFF Waveform vs LIGO GW150914 — PAPER_669 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFComparedToGW150914 (line 25030)
 *
 * GW150914: First direct GW detection (Sep 14, 2015, LIGO)
 *   m1 = 36 M☉, m2 = 29 M☉, M_final = 62 M☉, radiated = 3 M☉ c²
 *   Peak frequency: f_peak ~ 150 Hz, strain: h ~ 1e-21
 *   Distance: d ~ 410 Mpc
 *
 * UQFF modifications to strain:
 *   h_UQFF(t) = h_GR(t) * (1 - f_TRZ) * S_SCm * exp(-U_m * omega / c²)
 *   omega = 2 pi f_GW  (instantaneous angular frequency)
 *   S_SCm = exp(-rho_SCm * lambda_GW)   lambda_GW = c/f_GW
 *   Phase shift: dphi_UQFF = dphi_GR + kappa * f_TRZ * t
 *
 * Chirp mass: M_chirp = (m1*m2)^(3/5) / (m1+m2)^(1/5)
 * Inspiral frequency evolution: df/dt = 96/5 * pi^(8/3) * (G M_chirp)^(5/3) / c^5 * f^(11/3)
 */
class UQFFComparedToGW150914 {
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

    // GW150914 canonical parameters
    static constexpr double M1_GW150914 = 36.0*1.989e30;
    static constexpr double M2_GW150914 = 29.0*1.989e30;
    static constexpr double D_GW150914  = 410.0*3.086e22; // 410 Mpc in m
    static constexpr double F_PEAK      = 150.0;          // Hz
    UQFFComparedToGW150914();
    double compute_chirp_mass(double m1, double m2) const;
    double compute_h_GR(double Mc, double d, double f) const;
    double compute_h_UQFF(double Mc, double d, double f) const;
    double compute_df_dt(double Mc, double f) const;
    double compute_dphi_UQFF(double f, double t) const;
    double compute_S_SCm(double f) const;
    void simulate_inspiral(double f_start, double f_end, double df,
                           double m1, double m2, double d, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // UQFF_COMPARED_TO_GW150914_H
// UQFF helper: UQFFComparedToGW150914.cpp utils
static double _T_H(double M){
    return (UQFFComparedToGW150914::HBAR*UQFFComparedToGW150914::C*UQFFComparedToGW150914::C*UQFFComparedToGW150914::C)/
           (8.0*M_PI*UQFFComparedToGW150914::G*M*UQFFComparedToGW150914::K_B);
}
static double _L_H(double M){
    return (UQFFComparedToGW150914::HBAR*std::pow(UQFFComparedToGW150914::C,6.0))/
           (15360.0*M_PI*UQFFComparedToGW150914::G*UQFFComparedToGW150914::G*M*M);
}
static double _r_s(double M){
    return 2.0*UQFFComparedToGW150914::G*M/(UQFFComparedToGW150914::C*UQFFComparedToGW150914::C);
}
static double _tau_std(double M){
    return (5120.0*M_PI*UQFFComparedToGW150914::G*UQFFComparedToGW150914::G*M*M*M)/
           (UQFFComparedToGW150914::HBAR*std::pow(UQFFComparedToGW150914::C,4.0));
}
static double _U_m(double r, double t){
    double tn = t/UQFFComparedToGW150914::T_N_REF;
    return (UQFFComparedToGW150914::MU_J/r)*(1.0-std::exp(-UQFFComparedToGW150914::GAMMA*t*std::cos(M_PI*tn)));
}
