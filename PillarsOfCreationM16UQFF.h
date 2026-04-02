#ifndef PILLARS_OF_CREATION_M16_UQFF_H
#define PILLARS_OF_CREATION_M16_UQFF_H
// STANDALONE_PILLARSOFCREATIONM16UQFF
// PAPER_708: Pillars of Creation (Eagle Nebula M16) -- Full UQFF MUGE
// Star-forming pillars 6500-7000 ly, 4-5 ly length, cool molecular H2
// Source: grok_share_ba508f76c8e.txt entry #86 | Session 175
#include <vector>
#include <string>
#include <cmath>

class PillarsOfCreationM16UQFF {
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

    // M16 Pillars of Creation parameters
    double M_initial  = 1.01e4 * 1.989e30; // kg initial pillar mass (10100 M_sun)
    double r_pillar   = 4.731e16;           // m pillar height (5 ly half-span)
    double M_0_SF     = 0.9901;             // fractional secondary SF mass
    double tau_SF     = 3.156e13;           // s SF timescale (1 Myr)
    double E_0        = 0.10;               // fractional erosion rate
    double tau_erode  = 3.156e13;           // s erosion timescale
    double B_pillar   = 1.0e-6;             // T pillar magnetic field
    double v_gas      = 1.0e5;              // m/s gas velocity
    double q_p        = 1.602e-19;          // C proton charge
    double z_M16      = 1.2e-3;             // redshift (6500 ly)

    // Self-simulation state
    double time_step  = 1.578e13;           // 0.5 Myr
    mutable double curr_t = 1.578e13;       // start at 0.5 Myr

    PillarsOfCreationM16UQFF() = default;

    double M_dot(double t) const;
    double M_t(double t) const;

    // Erosion: (1 - E(t)) = 1 - E_0*exp(-t/tau_erode)
    double one_minus_E(double t) const;

    double H_z() const;
    double a_EM() const;

    // Full Pillars MUGE
    // g_P = G*M(t)/r^2*(1+H_z*t)*(1-B/B_crit)*(1-E(t)) + Ug_sum*(1+f_TRZ) + Lambda*c^2/3 + a_EM
    double g_Pillars(double t) const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // PILLARS_OF_CREATION_M16_UQFF_H
