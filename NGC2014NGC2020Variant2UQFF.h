#ifndef NGC2014_NGC2020_VARIANT2_UQFF_H
#define NGC2014_NGC2020_VARIANT2_UQFF_H
// STANDALONE_NGC2014NGC2020VARIANT2UQFF
// PAPER_711: NGC 2014 + NGC 2020 Variant 2 -- WR-dominated UQFF Analysis
// Second analysis with Wolf-Rayet cone structure and oxygen gas at 11000 C
// Source: grok_share_ba508f76c8e.txt entry #103 | Session 175
#include <vector>
#include <string>
#include <cmath>

class NGC2014NGC2020Variant2UQFF {
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

    // NGC 2014/2020 Variant 2 -- Wolf-Rayet emphasis
    double M_WR       = 200e3 * 1.989e30;  // kg WR-dominated region mass
    double r_cone     = 1.2e16;             // m  WR cone radius
    double L_WR       = 2.0e5 * 3.826e26;  // W  WR star luminosity (200k L_sun)
    double v_eject    = 3.0e6;              // m/s WR ejecta velocity
    double T_O_gas    = 11000.0;            // K  oxygen gas temperature
    double rho_WR     = 5.0e-22;           // kg/m^3 WR wind density
    double B_WR       = 1.5e-6;             // T  WR magnetic field
    double v_gas      = 2.0e5;              // m/s gas velocity
    double q_p        = 1.602e-19;          // C  proton charge
    double tau_WR     = 3.0e15;             // s  WR lifetime (~100 Myr)
    double M_0_WR     = 50.0;              // fractional SF from WR winds
    double tau_SF_WR  = 2.0e14;            // s  WR-driven SF timescale

    // Self-simulation state
    double time_step  = 7.89e13;
    mutable double curr_t = 7.89e13;

    NGC2014NGC2020Variant2UQFF() = default;

    double M_WR_t(double t) const;
    double a_WR_radiation() const;
    double a_EM_WR() const;

    // Full WR cone MUGE variant 2
    // g_v2 = G*M_WR(t)/r^2*(1+H0*t)*(1+f_TRZ) + a_WR_rad + a_EM_WR
    double g_WRcone(double t) const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // NGC2014_NGC2020_VARIANT2_UQFF_H
