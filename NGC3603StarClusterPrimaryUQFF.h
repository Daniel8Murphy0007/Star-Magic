#ifndef NGC3603_STAR_CLUSTER_PRIMARY_UQFF_H
#define NGC3603_STAR_CLUSTER_PRIMARY_UQFF_H
// STANDALONE_NGC3603STARCLUSTERPRIMARYUQFF
// PAPER_706: NGC 3603 Star Cluster (Primary Analysis) -- UQFF MUGE
// Milky Way Carina arm, most massive young cluster, ~1 Myr age
// Source: grok_share_ba508f76c8e.txt entry #89 | Session 175
#include <vector>
#include <string>
#include <cmath>

class NGC3603StarClusterPrimaryUQFF {
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

    // NGC 3603 primary parameters (document-base variant)
    double M_base     = 4.0e5 * 1.989e30; // kg initial cluster mass
    double r_clust    = 9.0e15;            // m cluster radius (slightly larger)
    double M_gas_frac = 0.30;              // fraction of mass as unformed gas
    double tau_SF     = 3.156e13;          // s SF timescale  
    double alpha_wind = 0.05;              // wind fractional suppression
    double B_clust    = 1.2e-5;            // T cluster magnetic field
    double v_gas      = 1.2e5;             // m/s cluster gas velocity
    double q_p        = 1.602e-19;         // C proton charge
    double z_ngc      = 3.0e-4;           // redshift (20,000 ly)
    double rho_gas    = 1.2e-20;           // kg/m^3 gas density

    // Self-simulation state
    double time_step  = 3.156e13;          // 1 Myr
    mutable double curr_t = 5.0e12;        // 0.16 Myr initial age

    NGC3603StarClusterPrimaryUQFF() = default;

    // Effective mass including gas conversion
    // M_eff(t) = M_base * (1 + M_gas_frac * exp(-t/tau_SF))
    double M_eff(double t) const;

    // Gravitational acceleration
    double g_grav(double t) const;

    // EM aether term
    double a_EM() const;

    // Hubble factor at z~3e-4
    double H_z() const;

    // Full primary MUGE
    // g = g_grav*(1+H_z*t)*(1-alpha_wind)*(1+f_TRZ) + a_EM
    double g_NGC3603primary(double t) const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // NGC3603_STAR_CLUSTER_PRIMARY_UQFF_H
