#ifndef NGC2014_NGC2020_STARFORMING_UQFF_H
#define NGC2014_NGC2020_STARFORMING_UQFF_H
// STANDALONE_NGC2014NGC2020STARFORMINGUQFF
// PAPER_710: NGC 2014 + NGC 2020 (Tapestry of Blazing Starbirth) UQFF
// LMC star-forming nebula pair, 163,000 ly, Wolf-Rayet + OB star dynamics
// Source: grok_share_ba508f76c8e.txt entry #84 | Session 175
#include <vector>
#include <string>
#include <cmath>

class NGC2014NGC2020StarformingUQFF {
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

    // NGC 2014/2020 parameters (LMC Tapestry of Blazing Starbirth)
    double M_initial  = 240.0 * 1.989e30;  // kg initial OB/WR stellar mass
    double r_region   = 9.461e16;           // m region span (10 ly)
    double M_0_SF     = 41.67;              // fractional secondary SF (10^4/240)
    double tau_SF     = 1.578e14;           // s SF timescale (5 Myr)
    double B_LMC      = 1.0e-6;             // T LMC magnetic field
    double v_gas      = 1.0e5;              // m/s gas velocity
    double q_p        = 1.602e-19;          // C proton charge
    double H0_LMC     = 2.184e-18;          // s^-1 Hubble constant

    // Self-simulation state
    double time_step  = 7.89e13;            // 2.5 Myr
    mutable double curr_t = 7.89e13;        // start at 2.5 Myr

    NGC2014NGC2020StarformingUQFF() = default;

    double M_dot(double t) const;
    double M_t(double t) const;
    double a_EM() const;

    // Full LMC starbirth MUGE
    // g_TB = G*M(t)/r^2*(1+H0*t)*(1-B/Bc)
    //       + (Ug1+Ug4)*(1+f_TRZ) + Lambda*c^2/3 + a_EM
    double g_Starbirth(double t) const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // NGC2014_NGC2020_STARFORMING_UQFF_H
