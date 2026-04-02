#ifndef NGC2525_WITH_SUPERNOVAE_SN2018GV_H
#define NGC2525_WITH_SUPERNOVAE_SN2018GV_H
// STANDALONE_NGC2525WITHSUPERNOVAESN2018GV
// PAPER_697: NGC 2525 — Barred Spiral with Type Ia Supernova SN 2018gv
// Barred spiral, 60 Mly, Hubble Time-lapse SN 2018gv observations
// Source: grok_share_ba508f76c8e.txt (#87/88) | Session 174
#include <cmath>
#include <string>

class NGC2525WithSupernovaeSN2018gv {
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

    // NGC 2525 parameters
    double M_galaxy    = 3.0e10 * M_sun;  // total galaxy mass
    double R_bar       = 8.0e3  * kpc;    // bar length
    double distance    = 60.0   * Mpc;    // to Earth
    double SFR_gal     = 1.5;             // star formation rate M_sun/yr

    // SN 2018gv Type Ia parameters
    double L_SN_peak   = 1.5e43;          // peak luminosity W
    double t_rise      = 14.0 * 86400.0;  // rise time s
    double t_decline   = 60.0 * 86400.0;  // 60-day decline
    double E_SN        = 1.0e44;          // kinetic energy J
    double M_Ni56      = 0.6 * M_sun;     // Ni-56 mass (powers light curve)

    double time_step = 86400.0;  // 1 day steps for SN
    mutable double curr_t = 0.0;

    NGC2525WithSupernovaeSN2018gv() = default;

    // Phillips relation: Delta_m_15 for Type Ia distance
    // M_B ~ -19.3 + 0.74*(Delta_m_15 - 1.1) [mag]
    double absolute_magnitude(double delta_m15) const;

    // Arnett's rule: L_peak ~ E_Ni / t_rise * radioactive_decay
    // L(t) = L_peak * exp(-t/t_rise) * (1 - exp(-t_decline/t)) for simplified model
    double L_light_curve(double t) const;

    // UQFF-modified SN luminosity
    // L_UQFF = L_SN * (1 + rho_SCm/rho_UA) * (1 - f_TRZ)
    double L_SN_UQFF(double t) const;

    // Host galaxy rotation curve (bar potential)
    double v_circ_barred(double R) const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // NGC2525_WITH_SUPERNOVAE_SN2018GV_H
