#ifndef LISA_VS_LIGO_COMPARISONS_H
#define LISA_VS_LIGO_COMPARISONS_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class LISAVsLIGOComparisons
 * @brief UQFF LISA vs LIGO Cross-Comparison — PAPER_678 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~20057
 *
 * Quantitative comparison of UQFF corrections in the LISA band (0.1 mHz – 1 Hz)
 * versus the LIGO band (10 – 2000 Hz):
 *
 *   LIGO band: f ~ 10-2000 Hz, U_m suppression dominant (S_Um small at high f)
 *   LISA band: f ~ 1e-4 – 1 Hz, S_UA_LISA suppression dominant (long arm)
 *
 *   Suppression ratio: R_supp(f) = h_UQFF(f) / h_GR(f)
 *   LIGO regime:  R ~ (1-f_TRZ) * S_SCm  (U_m term ~ 1 at low f)
 *   LISA regime:  R ~ (1-f_TRZ) * S_UA * S_SCm
 *
 *   Crossover frequency: f_cross where S_Um = S_UA
 *   f_cross = c^2 * rho_UA * L_LISA / (2 pi U_m * k_B * T_H)^{-1}
 *
 *   Sensitivity noise floor UQFF correction:
 *   S_n,UQFF(f) = S_n,GR(f) * (1 + F_TRZ * rho_UA/rho_SCm * f^(-2/3))
 */
class LISAVsLIGOComparisons {
public:
    static constexpr double G       = 6.6743e-11;
    static constexpr double C       = 2.998e8;
    static constexpr double HBAR    = 1.0546e-34;
    static constexpr double K_B     = 1.380649e-23;
    static constexpr double PI      = 3.14159265358979323846;
    static constexpr double M_SUN   = 1.989e30;
    static constexpr double RHO_UA  = 7.09e-36;
    static constexpr double RHO_SCM = 7.09e-37;
    static constexpr double F_TRZ   = 0.1;
    static constexpr double KAPPA   = 0.0005;
    static constexpr double SSQ     = 0.57;
    static constexpr double MU_J    = 3.38e23;
    static constexpr double GAMMA   = 5.0e-5 / 86400.0;
    static constexpr double T_N_REF = 1.0e8;

    static constexpr double L_LISA = 2.5e9;
    static constexpr double L_LIGO = 4.0e3;

    LISAVsLIGOComparisons();
    double compute_R_supp_LIGO(double M, double f, double t=1e8) const;
    double compute_R_supp_LISA(double M, double f, double T_eff=2.73) const;
    double compute_crossover_freq(double M, double T_eff=2.73) const;
    double compute_Sn_UQFF(double Sn_GR, double f) const;
    void compare_sweep(double M, double d, double f_ligo_lo, double f_ligo_hi,
                       double f_lisa_lo, double f_lisa_hi,
                       const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // LISA_VS_LIGO_COMPARISONS_H
