#ifndef HORSEHEAD_NEBULA_BARNARD33_UQFF_H
#define HORSEHEAD_NEBULA_BARNARD33_UQFF_H
// STANDALONE_HORSEHEADNEBULABARNARD33UQFF
// PAPER_704: Horsehead Nebula (Barnard 33) -- Infrared UQFF Evolution
// Dark nebula in Orion Molecular Cloud, 1500 ly distance
// Source: grok_share_ba508f76c8e.txt entry #93 | Session 175
#include <vector>
#include <string>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class HorseheadNebulaBarnard33UQFF {
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

    // Horsehead Nebula parameters
    double M_neb      = 120.0 * 1.989e30; // kg  mass (100 gas + 20 young stars M_sun)
    double r_neb      = 1.182e16;         // m   half-span of 2.5 ly nebula
    double L_star     = 1.0e5 * 3.826e26; // W   Sigma Orionis luminosity
    double rho_neb    = 1.0e-21;          // kg/m^3 nebula density
    double B_neb      = 1.0e-5;           // T   nebula magnetic field
    double v_gas      = 1.0e5;            // m/s gas velocity
    double tau_erode  = 1.578e14;         // s   erosion timescale (5 Myr)
    double E_0        = 0.2;              // fractional erosion rate
    double z_neb      = 3.0e-4;          // redshift (1500 ly)
    double q_p        = 1.602e-19;        // C proton charge

    // Self-simulation state
    double time_step  = 3.156e13;         // 1 Myr
    mutable double curr_t = 3.156e13;     // start at 1 Myr

    HorseheadNebulaBarnard33UQFF() = default;

    // Base gravitational acceleration
    // g_grav = G*M_neb/r^2
    double g_grav() const;

    // Erosion factor (1 - E(t))
    // E(t) = E_0*(1 - exp(-t/tau_erode))
    double one_minus_E(double t) const;

    // Hubble correction at z~3e-4
    double H_z() const;

    // Radiation pressure from Sigma Orionis
    // P_rad = (L_star/(4*pi*r^2*c)) * (rho_neb/m_H)
    double a_rad() const;

    // EM aether term
    // a_EM = q_p*v_gas*B_neb*(1+rho_UA/rho_SCm)*1e-12
    double a_EM() const;

    // Full MUGE for Horsehead Nebula
    // g_HH = g_grav*(1+H_z*t)*(1-E(t))*(1+f_TRZ) + a_rad + a_EM
    double g_Horsehead(double t) const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // HORSEHEAD_NEBULA_BARNARD33_UQFF_H
