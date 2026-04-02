#ifndef NGC3603_STAR_CLUSTER_2_UQFF_H
#define NGC3603_STAR_CLUSTER_2_UQFF_H
// STANDALONE_NGC3603STARCLUSTER2UQFF
// PAPER_705: NGC 3603 Star Cluster Variant 2 -- Refined MUGE Analysis
// Extreme star-forming region, Carina arm, 20,000 ly, 400,000 M_sun cluster
// Source: grok_share_ba508f76c8e.txt entry #90 | Session 175
#include <vector>
#include <string>
#include <cmath>

class NGC3603StarCluster2UQFF {
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

    // NGC 3603 parameters (variant 2 emphasis on Bok globule secondary SF)
    double M_initial  = 4.0e5 * 1.989e30; // kg initial cluster mass (400k M_sun)
    double r_clust    = 8.998e15;          // m cluster radius (19/2 ly)
    double M_0        = 0.10;              // fractional secondary SF mass fraction
    double tau_SF     = 3.156e13;          // s SF timescale (1 Myr)
    double P_0        = 0.10;              // fractional stellar-wind pressure reduction
    double tau_exp    = 3.156e13;          // s cavity expansion timescale
    double B_clust    = 1.0e-5;            // T cluster magnetic field
    double v_gas      = 1.0e5;             // m/s gas velocity
    double q_p        = 1.602e-19;         // C proton charge
    double rho_gas    = 1.0e-20;           // kg/m^3 cluster gas density
    double v_wind_OB  = 2.0e6;             // m/s OB star wind velocity

    // Self-simulation state
    double time_step  = 1.578e13;          // 0.5 Myr
    mutable double curr_t = 1.578e13;      // start at 0.5 Myr

    NGC3603StarCluster2UQFF() = default;

    // M_dot: secondary star formation growth exponent
    double M_dot(double t) const;

    // M_t: total mass at time t
    double M_t(double t) const;

    // Stellar wind pressure reduction factor (1 - P(t))
    double one_minus_P(double t) const;

    // EM aether term
    double a_EM() const;

    // Full NGC 3603 v2 MUGE
    // g = G*M(t)/r^2 * (1+H0*t) * (1-P(t)) * (1+f_TRZ) + a_EM
    double g_NGC3603v2(double t) const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // NGC3603_STAR_CLUSTER_2_UQFF_H
