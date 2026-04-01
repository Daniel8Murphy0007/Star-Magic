#ifndef UQFF_PBH_DARK_MATTER_IMPLICATIONS_H
#define UQFF_PBH_DARK_MATTER_IMPLICATIONS_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class UQFFPBHDarkMatterImplications
 * @brief UQFF Primordial BH Dark Matter Implications — PAPER_685 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~22222
 *
 * Cosmological dark matter fraction from stable PBHs in UQFF:
 *
 *   Standard: PBHs with M < M_crit_std evaporate before today
 *             M_crit_std = (hbar c^4 * t_age / (5120 pi G^2))^(1/3)
 *             M_crit_std ~ 5e11 kg  (only heavier PBHs survive)
 *
 *   UQFF: M_crit_UQFF lowered by factor ~0.033 suppression -> smaller PBHs survive
 *         M_crit_UQFF = M_crit_std * (tau_ratio)^(-1/3)
 *         where tau_ratio = tau_UQFF/tau_std ~ 30
 *         -> M_crit_UQFF ~ M_crit_std / 30^(1/3) ~ 0.32 * M_crit_std
 *
 *   DM fraction:
 *     f_PBH = rho_PBH / rho_DM
 *     rho_PBH ~ M_f * n_PBH  (n_PBH from inflationary spectrum)
 *     UQFF enhances f_PBH by keeping lower-mass PBHs alive:
 *     f_PBH_UQFF = f_PBH_std * 30^(2/3)  (more PBHs survive)
 *
 *   Observational constraints: microlensing, CMB distortions, GW background
 *   UQFF shifts the viable mass window: 1e10 – 1e17 kg fully open in UQFF.
 */
class UQFFPBHDarkMatterImplications {
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

    static constexpr double T_AGE = 4.34e17;  // 13.8 Gyr in s

    UQFFPBHDarkMatterImplications();
    double compute_M_crit_std() const;
    double compute_M_crit_UQFF() const;
    double compute_f_PBH_boost() const;
    double compute_rho_PBH_UQFF(double M_f, double n_PBH) const;
    bool is_DM_viable_std(double M) const;
    bool is_DM_viable_UQFF(double M) const;
    void scan_mass_window(double M_min, double M_max, int N,
                          const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // UQFF_PBH_DARK_MATTER_IMPLICATIONS_H
