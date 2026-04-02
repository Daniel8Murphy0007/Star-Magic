#ifndef WESTERLUND2_STAR_CLUSTER_UQFF_H
#define WESTERLUND2_STAR_CLUSTER_UQFF_H
// STANDALONE_WESTERLUND2STARCLUSTERUQFF
// PAPER_709: Westerlund 2 Young Super Star Cluster -- Full UQFF MUGE
// Milky Way Carina, 15000 ly, ~3000 stars, 2 Myr age, M_stars~30000 M_sun
// Source: grok_share_ba508f76c8e.txt entry #85 | Session 175
#include <vector>
#include <string>
#include <cmath>

class Westerlund2StarClusterUQFF {
public:

    // UQFF Universal Constants
    static constexpr double G        = 6.6743e-11;  // m^3 kg^-1 s^-2
    static constexpr double c        = 3.0e8;       // m/s
    static constexpr double hbar     = 1.0546e-34;  // J*s
    static constexpr double mu_0     = 1.2566e-6;   // H/m
    static constexpr double k_B      = 1.3806e-23;  // J/K
    static constexpr double M_sun    = 1.989e30;    // kg
    static constexpr double kpc      = 3.086e19;    // m (1 kpc)
    static constexpr double Mpc      = 3.086e22;    // m (1 Mpc)
    static constexpr double rho_UA   = 7.09e-36;    // J/m^3 Universal Aether
    static constexpr double rho_SCm  = 7.09e-37;    // J/m^3 Schwarzschild condensate
    static constexpr double f_TRZ    = 0.1;         // time-reversal zone factor
    static constexpr double kappa    = 5.0e-4;      // UQFF calibration day^-1
    static constexpr double SSq      = 0.57;        // superstring quenching
    static constexpr double mu_J     = 3.38e23;     // J*m magnetic string coupling
    static constexpr double Lambda   = 1.1e-52;     // m^-2 cosmological constant
    static constexpr double H_0      = 2.269e-18;   // s^-1 Hubble constant
    static constexpr double t_H      = 4.35e17;     // s Hubble time

    // Westerlund 2 parameters
    double M_init     = 3.0e4  * 1.989e30; // kg initial stellar mass
    double r_wd2      = 9.461e16;           // m cluster radius (10 ly)
    double M_0_gas    = 3.333;              // fractional gas-to-star mass ratio
    double tau_SF     = 6.312e13;           // s SF timescale (2 Myr)
    double B_wd2      = 1.0e-5;             // T cluster magnetic field
    double v_gas      = 1.0e5;              // m/s gas velocity
    double q_p        = 1.602e-19;          // C proton charge
    double H0_wd2     = 2.184e-18;          // s^-1 H_0 for Westerlund calc

    // Self-simulation state
    double time_step  = 3.156e13;           // 1 Myr
    mutable double curr_t = 3.156e13;       // start at 1 Myr

    Westerlund2StarClusterUQFF() = default;

    // M_dot: rapid star formation from initial gas reservoir
    double M_dot(double t) const;
    double M_t(double t) const;

    // EM aether term
    double a_EM() const;

    // Full Westerlund 2 MUGE
    // g_W2 = G*M(t)/r^2*(1+H0*t)*(1-B/B_crit)
    //       + (Ug1+Ug4)*(1+f_TRZ) + Lambda*c^2/3 + a_EM
    double g_Westerlund2(double t) const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // WESTERLUND2_STAR_CLUSTER_UQFF_H
