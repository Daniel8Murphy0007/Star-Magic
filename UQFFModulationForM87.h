#ifndef UQFF_MODULATION_FOR_M87_H
#define UQFF_MODULATION_FOR_M87_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFModulationForM87
 * @brief UQFF Modulation for M87 SMBH — PAPER_686 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~22436
 *
 * M87 (Messier 87 / NGC 4486):
 *   SMBH mass: M87* = 6.5e9 Msun = 1.293e40 kg (EHT April 2019)
 *   Distance: 16.4 Mpc = 5.06e23 m
 *   Shadow radius: r_shadow = 3*sqrt(3) * G*M/c^2 ~ 2e13 m
 *   Jet power: P_jet ~ 10^44 erg/s = 10^37 W
 *
 * UQFF modifications for M87*:
 *   1. Hawking T: T_UQFF = T_H * (1+f_TRZ) * (1-rho_SCm/rho_UA)
 *   2. Shadow size: r_sh,UQFF = r_sh * (1 + f_TRZ * rho_UA/rho_SCm)^0.5
 *   3. Jet power: P_jet,UQFF = P_jet,BZ * (1+f_TRZ) * (rho_UA/rho_SCm)^0.5
 *      Blandford-Znajek: P_BZ ~ kappa * Phi_BH^2 * Omega_H^2 / (4 pi c)
 *   4. Ring brightness: B_UQFF = B_GR * (rho_UA/rho_SCm)^(f_TRZ/4)
 *   5. U_m coupling to accretion disk: T_acc,UQFF = T_acc * (1 + U_m/k_B T_H)
 */
class UQFFModulationForM87 {
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

    static constexpr double M_M87     = 6.5e9 * 1.989e30;
    static constexpr double D_M87     = 5.06e23;
    static constexpr double P_JET_GR  = 1.0e37;   // W (Blandford-Znajek estimate)
    static constexpr double B_FIELD   = 1.0e3;     // T near horizon

    UQFFModulationForM87();
    double compute_T_H(double M=M_M87) const;
    double compute_T_UQFF(double M=M_M87) const;
    double compute_shadow_radius(double M=M_M87) const;
    double compute_shadow_UQFF(double M=M_M87) const;
    double compute_jet_power_UQFF(double P_BZ=P_JET_GR) const;
    double compute_ring_brightness_UQFF(double B_GR) const;
    double compute_T_accretion_UQFF(double T_acc, double M=M_M87) const;
    void simulate_T_evolution(double M_start, double M_end, double dM,
                              const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // UQFF_MODULATION_FOR_M87_H
