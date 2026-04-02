#ifndef M51_WHIRLPOOL_TIDAL_INTERACTION_H
#define M51_WHIRLPOOL_TIDAL_INTERACTION_H
// STANDALONE_M51WHIRLPOOLTIDALINTERACTION
// PAPER_692: M51 Whirlpool Galaxy — Tidal Arm Formation and UQFF
// Interacting spiral pair: M51a (NGC 5194) + M51b (NGC 5195)
// Source: grok_share_ba508f76c8e.txt | Session 174
#include <cmath>
#include <string>

class M51WhirlpoolTidalInteraction {
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

    // M51 system parameters
    double M_main     = 1.6e11 * M_sun;  // M51a primary mass
    double M_comp     = 2.0e10 * M_sun;  // M51b companion mass
    double d_sep      = 7.5e3 * kpc;     // current separation
    double R_spiral   = 15.0e3 * kpc;    // spiral arm radius
    double SFR        = 3.0;             // star formation rate M_sun/yr
    double z_M51      = 0.002;           // redshift
    double rho_ISM    = 1.67e-21;        // ISM density kg/m^3

    double time_step  = 3.156e13;
    mutable double curr_t = 0.0;

    M51WhirlpoolTidalInteraction() = default;

    // Tidal force on main disk from companion
    // F_tidal = 2*G*M_comp*R_d / d_sep^3 (differential tidal force)
    double tidal_force(double R_d) const;

    // Spiral arm pitch angle induced by tidal interaction
    double pitch_angle() const;

    // UQFF-modified MUGE for M51
    // g_M51 = G*(M_main+M_comp)/r^2 * (1+H(z)) * (1-B/B_crit) + rho_ISM*V*g_loc
    double g_M51_UQFF(double r) const;

    // Star formation efficiency (Kennicutt-Schmidt + UQFF)
    double SFE_UQFF() const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // M51_WHIRLPOOL_TIDAL_INTERACTION_H
