#ifndef NGC7635_BUBBLE_NEBULA_H
#define NGC7635_BUBBLE_NEBULA_H
// STANDALONE_NGC7635BUBBLENEBULA
// PAPER_695: NGC 7635 — Bubble Nebula, Stellar Wind Driven Cavity
// O-star (BD+60°2522) wind bubble, 7100 ly, ~10 ly diameter
// Source: grok_share_ba508f76c8e.txt (#91) | Session 174
#include <cmath>
#include <string>

class NGC7635BubbleNebula {
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

    // NGC 7635 / BD+60 2522 parameters
    double M_star      = 44.0  * M_sun;   // driving O-star mass
    double L_star      = 4.0e5 * 3.83e26; // stellar luminosity (4e5 L_sun)
    double M_dot_wind  = 1.0e-5 * M_sun / 3.156e7; // mass loss rate kg/s
    double v_wind      = 2.0e6;           // terminal wind velocity m/s
    double R_bubble    = 5.0  * 3.086e16; // bubble radius (5 ly in m)
    double rho_ISM     = 3.34e-24;        // ISM density kg/m^3
    double t_age       = 1.0e5 * 3.156e7; // bubble age ~100 kyr in s

    double time_step = 3.156e12;  // ~100 kyr steps
    mutable double curr_t = 0.0;

    NGC7635BubbleNebula() = default;

    // Wind mechanical luminosity: L_w = 0.5 * M_dot * v_wind^2
    double L_wind() const;

    // Bubble expansion velocity (Weaver et al. 1977)
    // R_b(t) = (L_w * t^3 / rho_ISM)^(1/5) * 0.88
    double R_bubble_expansion(double t) const;

    // UQFF stellar wind modification
    // v_wind_UQFF = v_wind * (1 + f_TRZ) * sqrt(rho_UA/rho_SCm)
    double v_wind_UQFF() const;

    // Shock compression ratio (strong shock: rho_2/rho_1 = (gamma+1)/(gamma-1))
    double shock_compression_ratio() const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // NGC7635_BUBBLE_NEBULA_H
