#ifndef UQFF_BLACK_HOLE_ACCRETION_MODEL_H
#define UQFF_BLACK_HOLE_ACCRETION_MODEL_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFBlackHoleAccretionModel
 * @brief UQFF BH Accretion Model — PAPER_670 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFBlackHoleAccretionModel (line 25268)
 *
 * Bondi accretion in UQFF:
 *   Mdot_Bondi = 4 pi lambda_B (G M)² rho_inf / cs³
 *   lambda_B = 1/4 (polytropic index gamma_ad=5/3)
 *   cs = sound speed = sqrt(gamma_ad k_B T_inf / m_p)  [proton mass m_p]
 *
 * UQFF modifications:
 *   rho_eff = rho_inf + rho_UA - rho_SCm  (aether mass contribution)
 *   S_TRZ   = (1 + f_TRZ)                 (negentropic inflow boost)
 *   S_Um    = 1 - exp(-U_m / k_B T_inf)   (magnetic string impedance)
 *   Mdot_UQFF = Mdot_Bondi * (rho_eff/rho_inf) * S_TRZ * S_Um
 *
 * Eddington accretion:
 *   L_Edd = 4 pi G M m_p c / sigma_T   [sigma_T Thompson cross-section]
 *   Mdot_Edd = L_Edd / (eta c²)         [eta=0.1 efficiency]
 *   Mdot_UQFF vs. Mdot_Edd ratio
 *
 * Simulation: M(t) evolution from dM/dt = Mdot_UQFF - Mdot_Hawking
 */
class UQFFBlackHoleAccretionModel {
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

    static constexpr double M_P = 1.6726e-27;   // kg proton mass
    static constexpr double SIGMA_T = 6.652e-29; // m2 Thompson
    static constexpr double ETA = 0.1;           // radiative efficiency
    UQFFBlackHoleAccretionModel();
    double compute_cs(double T_inf) const;
    double compute_Mdot_Bondi(double M, double rho_inf, double T_inf) const;
    double compute_Mdot_UQFF(double M, double rho_inf, double T_inf) const;
    double compute_L_Edd(double M) const;
    double compute_Mdot_Edd(double M) const;
    double compute_Eddington_ratio(double M, double rho_inf, double T_inf) const;
    void simulate_M_evolution(double M0, double t0, double t1, double dt,
                              double rho_inf, double T_inf, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // UQFF_BLACK_HOLE_ACCRETION_MODEL_H
// UQFF helper: UQFFBlackHoleAccretionModel.cpp utils
static double _T_H(double M){
    return (UQFFBlackHoleAccretionModel::HBAR*UQFFBlackHoleAccretionModel::C*UQFFBlackHoleAccretionModel::C*UQFFBlackHoleAccretionModel::C)/
           (8.0*M_PI*UQFFBlackHoleAccretionModel::G*M*UQFFBlackHoleAccretionModel::K_B);
}
static double _L_H(double M){
    return (UQFFBlackHoleAccretionModel::HBAR*std::pow(UQFFBlackHoleAccretionModel::C,6.0))/
           (15360.0*M_PI*UQFFBlackHoleAccretionModel::G*UQFFBlackHoleAccretionModel::G*M*M);
}
static double _r_s(double M){
    return 2.0*UQFFBlackHoleAccretionModel::G*M/(UQFFBlackHoleAccretionModel::C*UQFFBlackHoleAccretionModel::C);
}
static double _tau_std(double M){
    return (5120.0*M_PI*UQFFBlackHoleAccretionModel::G*UQFFBlackHoleAccretionModel::G*M*M*M)/
           (UQFFBlackHoleAccretionModel::HBAR*std::pow(UQFFBlackHoleAccretionModel::C,4.0));
}
static double _U_m(double r, double t){
    double tn = t/UQFFBlackHoleAccretionModel::T_N_REF;
    return (UQFFBlackHoleAccretionModel::MU_J/r)*(1.0-std::exp(-UQFFBlackHoleAccretionModel::GAMMA*t*std::cos(M_PI*tn)));
}
