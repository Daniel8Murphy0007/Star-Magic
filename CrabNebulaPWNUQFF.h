#ifndef CRAB_NEBULA_PWN_UQFF_H
#define CRAB_NEBULA_PWN_UQFF_H
// STANDALONE_CRABNEBULAPWNUQFF
// PAPER_694: Crab Nebula (M1 / NGC 1952) — SNR + Pulsar Wind Nebula UQFF
// 6500 ly, 11 ly diameter, Crab Pulsar 33 ms spin period
// Source: grok_share_ba508f76c8e.txt (#100) | Session 174
#include <cmath>
#include <string>

class CrabNebulaPWNUQFF {
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

    // Crab Nebula / Pulsar parameters
    double M_ejecta   = 4.6  * M_sun;    // supernova ejecta mass
    double E_SN       = 1.0e44;          // SN explosion energy J
    double R_SNR      = 5.5  * 3.086e16; // SNR radius (5.5 ly in m)
    double P_pulsar   = 0.033;           // Crab pulsar spin period s
    double B_pulsar   = 3.8e8;           // pulsar surface B field T
    double Omega_psr  = 2*M_PI/0.033;   // angular velocity rad/s
    double I_psr      = 1.0e38;          // moment of inertia kg*m^2
    double age_SNR    = 975.0 * 3.156e7; // age since SN 1054 (s)
    double rho_PWN    = 1.67e-25;        // PWN density kg/m^3

    double time_step = 3.156e10; // ~1000 yr
    mutable double curr_t = 0.0;

    CrabNebulaPWNUQFF() = default;

    // Pulsar spin-down luminosity
    // L_sd = I * Omega * |dOmega/dt|
    // dOmega/dt = -B^2 * R^6 * Omega^3 / (6*I*c^3) for magnetic dipole
    double spin_down_luminosity() const;

    // Synchrotron cooling time for electrons in PWN
    // t_sync = 9*m_e^3*c^5 / (4*e^4*B^2*gamma_e)
    double synchrotron_cooling(double gamma_e) const;

    // UQFF-modified SNR expansion velocity
    // v_SNR = sqrt(2*E_SN*(1-f_TRZ) / M_ejecta) * (rho_SCm/rho_UA correction)
    double v_SNR_UQFF() const;

    // Magnetic energy density in PWN
    double B_energy_PWN() const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // CRAB_NEBULA_PWN_UQFF_H
