#ifndef M87_MASS_EVOLUTION_SIMULATION_H
#define M87_MASS_EVOLUTION_SIMULATION_H
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

/**
 * @class M87MassEvolutionSimulation
 * @brief M87* SMBH Mass Evolution in UQFF — PAPER_687 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~22623
 *
 * M87*: 6.5e9 Msun at D=16.4 Mpc. Combines UQFF accretion with UQFF-suppressed
 * Hawking evaporation for a full mass evolution simulation.
 *
 * Coupled equations:
 *   dM/dt_acc  = Mdot_Bondi_UQFF  (accretion from ISM)
 *   dM/dt_evap = dM/dt_Hawking_UQFF  (suppressed Hawking)
 *   dM/dt_jet  = -P_jet_UQFF / c^2  (jet mass loss)
 *   dM/dt_total = dM_acc - |dM_evap| - dM_jet
 *
 * Bondi-UQFF:
 *   Mdot_Bondi = 4 pi G^2 M^2 rho_inf / (c_s^3)
 *   Mdot_UQFF  = Mdot_Bondi * (rho_inf + rho_UA - rho_SCm) / rho_inf * (1+f_TRZ)
 *
 * Jet (Blandford-Znajek UQFF):
 *   P_jet = kappa_BZ * Phi_BH^2 * Omega_H^2 / (4 pi c) * (1+f_TRZ) * sqrt(rho_UA/rho_SCm)
 *   Omega_H ~ c/(4 G M) (maximal spin approximation)
 *
 * Simulate: M(t) over 14 Gyr (age of universe). Compare GR vs UQFF evolution.
 */
class M87MassEvolutionSimulation {
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

    static constexpr double M_M87    = 6.5e9 * 1.989e30;
    static constexpr double RHO_ISM  = 1.67e-25;   // kg/m^3 M87 environment
    static constexpr double T_ISM    = 1.0e7;       // K hot plasma
    static constexpr double T_HUBBLE = 4.34e17;     // s (13.8 Gyr)

    M87MassEvolutionSimulation();
    double compute_sound_speed_ISM(double T_ISM) const;
    double compute_Mdot_Bondi_UQFF(double M, double rho_inf, double T_inf) const;
    double compute_dM_dt_evap_UQFF(double M, double t) const;
    double compute_jet_power_UQFF(double M) const;
    double compute_dM_dt_total(double M, double t,
                               double rho_inf=RHO_ISM,
                               double T_inf=T_ISM) const;
    std::vector<std::pair<double,double>> evolve(double M0, double t_end, double dt,
                                                 double rho_inf=RHO_ISM,
                                                 double T_inf=T_ISM) const;
    void simulate_over_hubble(double dt, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
};
#endif // M87_MASS_EVOLUTION_SIMULATION_H
