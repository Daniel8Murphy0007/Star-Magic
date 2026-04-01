#ifndef UQFF_COMPARED_TO_LIGO_DATA_H
#define UQFF_COMPARED_TO_LIGO_DATA_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFComparedToLIGOData
 * @brief UQFF vs General LIGO Dataset — PAPER_674 | Session 173 | April 1, 2026
 * Source: grok_share_fc21e30c24b4.txt line ~18514
 *
 * Compares UQFF-modified GW strain and phase to the generalised LIGO O1/O2/O3
 * dataset. UQFF corrections enter through three suppression factors:
 *
 *   h_UQFF(f) = h_GR(f) * (1 - f_TRZ) * exp(-rho_SCm * r_s / k_B T_H)
 *                        * exp(-U_m * 2*pi*f / c^2)
 *
 *   Delta_phi_UQFF = Delta_phi_GR + kappa * f_TRZ * t_coal
 *
 * LIGO sensitivity band: 10–2000 Hz.
 * Chirp mass: Mc = (m1*m2)^(3/5) / (m1+m2)^(1/5)
 * GW luminosity (Peters): P_GW = -32/5 * G^4/c^5 * m1^2*m2^2*(m1+m2)/a^5
 * UQFF modifies via S_UA * S_SCm * S_TRZ * S_Um factors.
 * Self-simulate: sweep frequency 10–2000 Hz, output h_UQFF vs h_GR.
 */
class UQFFComparedToLIGOData {
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

    UQFFComparedToLIGOData();
    double compute_chirp_mass(double m1, double m2) const;
    double compute_h_GR(double Mc, double d, double f) const;
    double compute_S_SCm(double M, double T_H) const;
    double compute_S_Um(double M, double r, double t, double f) const;
    double compute_h_UQFF(double m1, double m2, double d, double f,
                          double r=0.0, double t=1e8) const;
    double compute_delta_phi(double t_coal) const;
    void simulate_frequency_sweep(double f_start, double f_end, double df,
                                  double m1, double m2, double d,
                                  const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // UQFF_COMPARED_TO_LIGO_DATA_H
