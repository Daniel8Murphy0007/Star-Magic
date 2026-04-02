#ifndef NGC2525_BARRED_SPIRAL_2_UQFF_H
#define NGC2525_BARRED_SPIRAL_2_UQFF_H
// STANDALONE_NGC2525BARREDSPIRAL2UQFF
// PAPER_707: NGC 2525 Barred Spiral Galaxy (Variant 2) -- UQFF MUGE
// Host galaxy of Type Ia supernova SN2018gv, second analysis variant
// Source: grok_share_ba508f76c8e.txt entry #88 | Session 175
#include <vector>
#include <string>
#include <cmath>

class NGC2525BarredSpiral2UQFF {
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

    // NGC 2525 variant 2 parameters
    double M_stellar  = 3.0e10 * 1.989e30; // kg stellar mass
    double M_DM       = 9.0e10 * 1.989e30; // kg dark matter halo
    double r_gal      = 30.0e3 * 3.086e19; // m galaxy radius (30 kpc)
    double z_ngc      = 0.0051;             // redshift (721 Mly)
    double L_SN       = 1.0e43;             // W SN Ia peak luminosity
    double tau_SN     = 3.156e8;            // s SN decay timescale (10 yr)
    double v_bar      = 1.5e5;              // m/s bar streaming velocity
    double B_bar      = 3.0e-10;            // T bar magnetic field
    double q_p        = 1.602e-19;          // C proton charge

    // Self-simulation state
    double time_step  = 3.156e13;           // 1 Myr
    mutable double curr_t = 0.0;

    NGC2525BarredSpiral2UQFF() = default;

    // Total mass
    double M_tot() const;

    // SN Ia feedback impulse at t
    // a_SN(t) = L_SN/(M_tot*r_gal*c) * exp(-t/tau_SN)
    double a_SN(double t) const;

    // Hubble parameter at z=0.0051
    double H_z() const;

    // EM aether bar streaming term
    double a_EM_bar() const;

    // Full variant 2 MUGE
    // g_v2(t) = G*M_tot/r^2*(1+H_z*t)*(1+f_TRZ) + a_SN(t) + a_EM_bar
    double g_NGC2525v2(double t) const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // NGC2525_BARRED_SPIRAL_2_UQFF_H
