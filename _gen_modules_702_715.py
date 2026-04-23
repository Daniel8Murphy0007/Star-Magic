from dpm_helpers import dpm_ug1_seed, dpm_ug2_shell
"""
Session 175 -- Generate C++ .h and .cpp standalone modules for PAPER_702-715.
Source: grok_share_ba508f76c8e.txt
14 new modules: Saturn, NGC1275, HorseheadB33, NGC3603x2, NGC2525v2,
                PillarsM16, Westerlund2, NGC2014-2020, variants, KB17-19
"""
import os

ROOT = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
os.chdir(ROOT)

UQFF_CONSTS = """
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
"""

SELF_METHODS = """\

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);
"""

# ---------------------------------------------------------------------------
# PAPER_702: Saturn Ring System UQFF
# ---------------------------------------------------------------------------
def gen_702():
    cls   = "SaturnRingSystemUQFF"
    guard = "SATURN_RING_SYSTEM_UQFF_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_702: Saturn Ring System -- Master UQFF Gravity Equation
// Source: grok_share_ba508f76c8e.txt entry #98 | Session 175
#include <vector>
#include <string>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class {cls} {{
public:
{UQFF_CONSTS}
    // Saturn system parameters
    double M_Saturn   = 5.683e26;       // kg  Saturn mass
    double r_eq       = 6.0268e7;       // m   equatorial radius
    double M_ring     = 1.5e19;         // kg  ring system mass
    double r_ring     = 7.0e7;          // m   main ring radius
    double M_Sun      = 1.989e30;       // kg  solar mass (for orbit gravity)
    double r_orbit    = 1.427e12;       // m   Saturn orbital radius
    double v_wind     = 500.0;          // m/s atmospheric wind speed
    double rho_atm    = 2.0e-4;         // kg/m^3 atmosphere density
    double B_Sat      = 1.0e-7;         // T   Saturn magnetic field
    double q_p        = 1.602e-19;      // C   proton charge (EM term)
    double v_charged  = 1.0e4;          // m/s charged particle velocity
    double z_Sat      = 0.0;            // Saturn redshift (solar system)

    // Self-simulation state
    double time_step  = 3.156e7;        // 1 yr
    mutable double curr_t = 0.0;

    {cls}() = default;

    // Solar gravitational acceleration at Saturn orbit
    // g_solar = G*M_Sun/r_orbit^2 * (1 + H_0*t) * (1 + f_TRZ)
    double g_solar(double t) const;

    // Saturn self-gravity at equatorial radius
    // g_Saturn = G*M_Saturn/r_eq^2
    double g_self() const;

    // Ring tidal contribution
    // T_ring = G*M_ring/r_ring^2
    double T_ring() const;

    // Atmospheric wind dynamic pressure
    // a_wind = rho_atm * v_wind^2 / M_Saturn * V_atmosphere
    double a_wind() const;

    // Electromagnetic aether term
    // a_EM = q_p * v_charged * B_Sat * (1 + rho_UA/rho_SCm) * 1e-12
    double a_EM() const;

    // Full Saturn MUGE
    // g_Saturn = g_solar + g_self + T_ring + a_wind + a_EM
    double g_Saturn(double t) const;

    std::string primary_equation() const;
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>
#include <fstream>

// Solar gravitational pull at Saturn orbit, with Hubble + TRZ correction
// g_solar = G*M_Sun/r_orbit^2 * (1 + H_0*t) * (1 + f_TRZ)
double {cls}::g_solar(double t) const {{
    double g_base = dpm_ug1_seed(M_Sun, r_orbit);
    return g_base * (1.0 + H_0 * t) * (1.0 + f_TRZ);
}}

// Saturn self-gravity at equatorial surface
// g_self = G*M_Saturn/r_eq^2  => ~10.44 m/s^2
double {cls}::g_self() const {{
    return dpm_ug1_seed(M_Saturn, r_eq);
}}

// Ring tidal acceleration from main ring system
// T_ring = G*M_ring/r_ring^2 ~ 2.043e-7 m/s^2
double {cls}::T_ring() const {{
    return dpm_ug1_seed(M_ring, r_ring);
}}

// Wind dynamic pressure converted to acceleration
// a_wind ~ rho_atm * v_wind^2 scaled ~ 2.5e-7 m/s^2
double {cls}::a_wind() const {{
    double P_wind = rho_atm * v_wind * v_wind;
    // Scale by 1/M_Saturn * V_eff (effective column volume ~ 1e18 m^3)
    double V_eff = 1.0e18;
    return P_wind * V_eff / M_Saturn;
}}

// Electromagnetic aether term with [UA]/[SCm] correction
// a_EM = q_p*(v x B)*(1 + rho_UA/rho_SCm)*1e-12
//      = q_p*v_charged*B_Sat / m_proton * 11 * 1e-12
double {cls}::a_EM() const {{
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_charged * B_Sat / m_p;
    double correction = 1.0 + rho_UA / rho_SCm;  // = 11
    return a_bare * correction * 1.0e-12;
}}

// Full Saturn MUGE
// g_Saturn(t) = g_solar(t) + g_self + T_ring + a_wind + a_EM
// Numerical example (t=0): ~ 10.44 m/s^2 (surface gravity)
double {cls}::g_Saturn(double t) const {{
    return g_solar(t) + g_self() + T_ring() + a_wind() + a_EM();
}}

std::string {cls}::primary_equation() const {{
    return "g_Saturn(t) = G*M_Sun/r_orbit^2*(1+H0*t)*(1+f_TRZ) + G*M_Sat/r_eq^2 + G*M_ring/r_ring^2 + rho*v_wind^2/M_eff + q*(vxB)*(1+rho_UA/rho_SCm)*1e-12";
}}

void {cls}::self_update() {{
    curr_t += time_step;
    r_orbit *= 1.0 + 1.0e-10 * time_step; // very slow orbital drift
}}

void {cls}::self_expand() {{
    r_ring *= 1.001;
    M_ring *= 0.9999; // ring mass slowly dissipates
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        double g = g_Saturn(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e7 << " yr"
                  << "  g_Saturn=" << g << " m/s^2\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} sat;
    std::cout << "Saturn Ring System UQFF Simulation\\n";
    std::cout << sat.primary_equation() << "\\n\\n";
    std::cout << "g_Saturn(t=0) = " << sat.g_Saturn(0.0) << " m/s^2\\n";
    sat.simulate(5);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_703: NGC 1275 Magnetic Monster UQFF
# ---------------------------------------------------------------------------
def gen_703():
    cls   = "NGC1275MagneticMonsterUQFF"
    guard = "NGC1275_MAGNETIC_MONSTER_UQFF_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_703: NGC 1275 (Magnetic Monster / Perseus A) -- UQFF Master Gravity
// Type 1.5 Seyfert galaxy in Perseus Cluster, 237 Mly, SMBH 800e6 M_sun
// Source: grok_share_ba508f76c8e.txt entry #94 | Session 175
#include <vector>
#include <string>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class {cls} {{
public:
{UQFF_CONSTS}
    // NGC 1275 system parameters
    double M_SMBH     = 8.0e8 * 1.989e30;  // SMBH mass (800e6 M_sun)
    double M_stellar  = 1.0e12 * 1.989e30; // stellar mass
    double r_gal      = 9.46e20;           // m galaxy radius (100 kly)
    double z_ngc      = 0.0176;            // NGC 1275 redshift
    double B_fil      = 1.0e-8;            // T  filament magnetic field
    double v_HVS      = 3.0e6;             // m/s merger velocity
    double tau_BH     = 3.156e15;          // s  BH feedback timescale (100 Myr)
    double F_0        = 0.1;               // fractional feedback reduction
    double rho_ICM    = 1.0e-24;           // kg/m^3 intracluster medium density
    double V_fil      = 1.42e50;           // m^3 filament volume
    double M_fil      = 1.0e6 * 1.989e30;  // kg filament mass (1e6 M_sun)
    double q_p        = 1.602e-19;         // C proton charge

    // Self-simulation state
    double time_step  = 1.578e15;          // ~50 Myr
    mutable double curr_t = 1.578e15;      // start at 50 Myr

    {cls}() = default;

    // Base gravitational acceleration
    // g_grav = G*(M_SMBH + M_stellar)/r^2
    double g_grav(double r) const;

    // Black hole feedback factor
    // F_BH(t) = F_0*(1 - exp(-t/tau_BH))
    // Returns (1 - F_BH) to reduce g_grav
    double one_minus_FBH(double t) const;

    // Hubble expansion correction for z=0.0176
    // H(z) = H_0 * sqrt(0.3*(1+z)^3 + 0.7)
    double H_z() const;

    // Magnetic filament support acceleration
    // a_fil = B_fil^2/(2*mu_0) * V_fil / M_fil
    double a_fil() const;

    // Electromagnetic aether term
    // a_EM = q_p*v_HVS*B_fil*(1 + rho_UA/rho_SCm)*1e-12
    double a_EM() const;

    // Full NGC 1275 MUGE
    // g_NGC1275 = g_grav*(1+H_z*t)*(1-F_BH)*(1+f_TRZ) + a_fil + a_EM
    double g_NGC1275(double r, double t) const;

    std::string primary_equation() const;
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>

// Base gravitational acceleration
double {cls}::g_grav(double r) const {{
    double M_tot = M_SMBH + M_stellar;
    return dpm_ug1_seed(M_tot, r);
}}

// Hubble parameter at z=0.0176
// H(z) = H_0*sqrt(0.3*(1+z)^3 + 0.7)
double {cls}::H_z() const {{
    double z = z_ngc;
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + z, 3) + 0.7);
}}

// BH feedback reduction: (1 - F_BH(t))
double {cls}::one_minus_FBH(double t) const {{
    double F = F_0 * (1.0 - std::exp(-t / tau_BH));
    return 1.0 - F;
}}

// Filament magnetic support acceleration (upward, stabilising)
// Energy density: B^2/(2*mu_0) * V_fil = 5.649e39 N  => / M_fil => 2.840e3 m/s^2
// Scaled by 1e-12 for macroscopic system
double {cls}::a_fil() const {{
    double energy = (B_fil * B_fil) / (2.0 * mu_0) * V_fil;
    return energy / M_fil * 1.0e-12;
}}

// EM aether term
// a_EM = q_p * v_HVS * B_fil / m_p * (1 + rho_UA/rho_SCm) * 1e-12
// ~ 3.16e-5 m/s^2
double {cls}::a_EM() const {{
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_HVS * B_fil / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}}

// Full MUGE for NGC 1275
// g_NGC1275(r,t) = g_grav*(1+H*t)*(1-F_BH)*(1+f_TRZ) + a_fil + a_EM
// At t=50 Myr, r=r_gal: ~ 3.16e-5 m/s^2
double {cls}::g_NGC1275(double r, double t) const {{
    double Ht = H_z() * t;
    double g_base = g_grav(r) * (1.0 + Ht) * one_minus_FBH(t) * (1.0 + f_TRZ);
    return g_base + a_fil() + a_EM();
}}

std::string {cls}::primary_equation() const {{
    return "g_NGC1275(r,t)=G*(M_SMBH+M_star)/r^2*(1+H_z*t)*(1-F_BH(t))*(1+f_TRZ)+a_fil+q*(vxB)*(1+rho_UA/rho_SCm)*1e-12";
}}

void {cls}::self_update() {{
    curr_t += time_step;
    M_stellar *= (1.0 - 1.0e-5 * time_step / 3.156e7); // slow mass loss
}}

void {cls}::self_expand() {{
    r_gal *= 1.001;
    B_fil *= 0.999;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        double g = g_NGC1275(r_gal, t);
        std::cout << "Step " << step << "  t=" << t / 3.156e7
                  << " yr  g=" << g << " m/s^2\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} ngc;
    std::cout << "NGC 1275 Magnetic Monster UQFF\\n";
    std::cout << ngc.primary_equation() << "\\n\\n";
    double g = ngc.g_NGC1275(ngc.r_gal, 1.578e15);
    std::cout << "g(50 Myr) = " << g << " m/s^2\\n";
    ngc.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_704: Horsehead Nebula Barnard 33 UQFF
# ---------------------------------------------------------------------------
def gen_704():
    cls   = "HorseheadNebulaBarnard33UQFF"
    guard = "HORSEHEAD_NEBULA_BARNARD33_UQFF_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_704: Horsehead Nebula (Barnard 33) -- Infrared UQFF Evolution
// Dark nebula in Orion Molecular Cloud, 1500 ly distance
// Source: grok_share_ba508f76c8e.txt entry #93 | Session 175
#include <vector>
#include <string>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class {cls} {{
public:
{UQFF_CONSTS}
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

    {cls}() = default;

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
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>

double {cls}::g_grav() const {{
    return dpm_ug1_seed(M_neb, r_neb);
}}

// Hubble expansion at z~3e-4
double {cls}::H_z() const {{
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + z_neb, 3) + 0.7);
}}

// Erosion factor: as stars erode the nebula
// E(t) = E_0*(1 - exp(-t/tau)) => 1-E = 1 - 0.03626 at t=1 Myr
double {cls}::one_minus_E(double t) const {{
    double E = E_0 * (1.0 - std::exp(-t / tau_erode));
    return 1.0 - E;
}}

// Radiation pressure acceleration
// P_rad = (L/(4*pi*r^2*c)) * (rho/m_H) ~ 4.347e-5 m/s^2
double {cls}::a_rad() const {{
    const double m_H = 1.67e-27;
    double flux = L_star / (4.0 * M_PI * r_neb * r_neb * c);
    return flux * (rho_neb / m_H);
}}

// EM aether term ~ 1.053e-3 m/s^2
double {cls}::a_EM() const {{
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_neb / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}}

// Full MUGE
// g_HH = g_grav*(1+H_z*t)*(1-E(t))*(1+f_TRZ) + a_rad + a_EM ~ 1.097e-3 m/s^2
double {cls}::g_Horsehead(double t) const {{
    double g_base = g_grav() * (1.0 + H_z() * t) * one_minus_E(t) * (1.0 + f_TRZ);
    return g_base + a_rad() + a_EM();
}}

std::string {cls}::primary_equation() const {{
    return "g_HH(t)=G*M/r^2*(1+H_z*t)*(1-E(t))*(1+f_TRZ)+P_rad+(q*v_gas*B)*(1+rho_UA/rho_SCm)*1e-12";
}}

void {cls}::self_update() {{
    curr_t += time_step;
    M_neb *= (1.0 - E_0 * 0.01); // gradual mass erosion
}}

void {cls}::self_expand() {{
    r_neb *= 1.001;
    rho_neb *= 0.999;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        double g = g_Horsehead(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} hh;
    std::cout << "Horsehead Nebula Barnard 33 UQFF\\n";
    std::cout << hh.primary_equation() << "\\n\\n";
    std::cout << "g(1 Myr) = " << hh.g_Horsehead(3.156e13) << " m/s^2\\n";
    hh.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_705: NGC 3603 Star Cluster Variant 2 UQFF
# ---------------------------------------------------------------------------
def gen_705():
    cls   = "NGC3603StarCluster2UQFF"
    guard = "NGC3603_STAR_CLUSTER_2_UQFF_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_705: NGC 3603 Star Cluster Variant 2 -- Refined MUGE Analysis
// Extreme star-forming region, Carina arm, 20,000 ly, 400,000 M_sun cluster
// Source: grok_share_ba508f76c8e.txt entry #90 | Session 175
#include <vector>
#include <string>
#include <cmath>

class {cls} {{
public:
{UQFF_CONSTS}
    // NGC 3603 parameters (variant 2 emphasis on Bok globule secondary SF)
    double M_initial  = 4.0e5 * 1.989e30; // kg initial cluster mass (400k M_sun)
    double r_clust    = 8.998e15;          // m cluster radius (19/2 ly)
    double M_0        = 0.10;              // fractional secondary SF mass fraction
    double tau_SF     = 3.156e13;          // s SF timescale (1 Myr)
    double P_0        = 0.10;              // fractional stellar-wind pressure reduction
    double tau_exp    = 3.156e13;          // s cavity expansion timescale
    double B_clust    = 1.0e-5;            // T cluster magnetic field
    double v_gas      = 1.0e5;             // m/s gas velocity
    double q_p        = 1.602e-19;         // C proton charge
    double rho_gas    = 1.0e-20;           // kg/m^3 cluster gas density
    double v_wind_OB  = 2.0e6;             // m/s OB star wind velocity

    // Self-simulation state
    double time_step  = 1.578e13;          // 0.5 Myr
    mutable double curr_t = 1.578e13;      // start at 0.5 Myr

    {cls}() = default;

    // M_dot: secondary star formation growth exponent
    double M_dot(double t) const;

    // M_t: total mass at time t
    double M_t(double t) const;

    // Stellar wind pressure reduction factor (1 - P(t))
    double one_minus_P(double t) const;

    // EM aether term
    double a_EM() const;

    // Full NGC 3603 v2 MUGE
    // g = G*M(t)/r^2 * (1+H0*t) * (1-P(t)) * (1+f_TRZ) + a_EM
    double g_NGC3603v2(double t) const;

    std::string primary_equation() const;
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>

// M_dot: Bok globule secondary star formation mass fraction
// M_dot(t) = M_0 * exp(-t/tau_SF)
double {cls}::M_dot(double t) const {{
    return M_0 * std::exp(-t / tau_SF);
}}

// M_t: M_initial * (1 + M_dot(t))
double {cls}::M_t(double t) const {{
    return M_initial * (1.0 + M_dot(t));
}}

// Stellar wind cavity pressure reduction
// P(t) = P_0 * exp(-t/tau_exp) => (1 - P) at t=0.5 Myr ~ 0.9394
double {cls}::one_minus_P(double t) const {{
    return 1.0 - P_0 * std::exp(-t / tau_exp);
}}

// EM aether term ~ 1.053e-3 m/s^2
double {cls}::a_EM() const {{
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_clust / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}}

// Full MUGE: g = G*M(t)/r^2 * (1+H0*t) * (1-P(t)) * (1+f_TRZ) + a_EM
// At t=0.5 Myr: ~ 1.053e-3 m/s^2
double {cls}::g_NGC3603v2(double t) const {{
    double M   = M_t(t);
    double g_n = dpm_ug1_seed(M, r_clust);
    return g_n * (1.0 + H_0 * t) * one_minus_P(t) * (1.0 + f_TRZ) + a_EM();
}}

std::string {cls}::primary_equation() const {{
    return "g_NGC3603v2(t)=G*M(t)/r^2*(1+H0*t)*(1-P(t))*(1+f_TRZ)+q_p*v_gas*B/m_p*(1+rho_UA/rho_SCm)*1e-12";
}}

void {cls}::self_update() {{
    curr_t += time_step;
}}

void {cls}::self_expand() {{
    r_clust  *= 1.001;
    M_initial *= 1.0001;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        double g = g_NGC3603v2(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} n3;
    std::cout << "NGC 3603 Star Cluster v2 UQFF\\n";
    std::cout << n3.primary_equation() << "\\n\\n";
    std::cout << "g(0.5 Myr) = " << n3.g_NGC3603v2(1.578e13) << " m/s^2\\n";
    n3.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_706: NGC 3603 Star Cluster Primary UQFF
# ---------------------------------------------------------------------------
def gen_706():
    cls   = "NGC3603StarClusterPrimaryUQFF"
    guard = "NGC3603_STAR_CLUSTER_PRIMARY_UQFF_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_706: NGC 3603 Star Cluster (Primary Analysis) -- UQFF MUGE
// Milky Way Carina arm, most massive young cluster, ~1 Myr age
// Source: grok_share_ba508f76c8e.txt entry #89 | Session 175
#include <vector>
#include <string>
#include <cmath>

class {cls} {{
public:
{UQFF_CONSTS}
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

    {cls}() = default;

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
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>

double {cls}::M_eff(double t) const {{
    return M_base * (1.0 + M_gas_frac * std::exp(-t / tau_SF));
}}

double {cls}::g_grav(double t) const {{
    return G * M_eff(t) / (r_clust * r_clust);
}}

double {cls}::H_z() const {{
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + z_ngc, 3) + 0.7);
}}

double {cls}::a_EM() const {{
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_clust / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}}

// g_NGC3603primary = g_grav*(1+H*t)*(1-alpha_wind)*(1+f_TRZ) + a_EM
double {cls}::g_NGC3603primary(double t) const {{
    double g_base = g_grav(t) * (1.0 + H_z() * t)
                              * (1.0 - alpha_wind)
                              * (1.0 + f_TRZ);
    return g_base + a_EM();
}}

std::string {cls}::primary_equation() const {{
    return "g_NGC3603p(t)=G*M_eff(t)/r^2*(1+H_z*t)*(1-alpha_wind)*(1+f_TRZ)+q*v*B/m_p*(1+rho_UA/rho_SCm)*1e-12";
}}

void {cls}::self_update() {{
    curr_t += time_step;
}}

void {cls}::self_expand() {{
    r_clust   *= 1.001;
    M_gas_frac *= 0.999; // gas fraction consumed over time
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        double g = g_NGC3603primary(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} n3p;
    std::cout << "NGC 3603 Primary UQFF\\n";
    std::cout << n3p.primary_equation() << "\\n\\n";
    std::cout << "g(1 Myr) = " << n3p.g_NGC3603primary(3.156e13) << " m/s^2\\n";
    n3p.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_707: NGC 2525 Barred Spiral Variant 2 UQFF
# ---------------------------------------------------------------------------
def gen_707():
    cls   = "NGC2525BarredSpiral2UQFF"
    guard = "NGC2525_BARRED_SPIRAL_2_UQFF_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_707: NGC 2525 Barred Spiral Galaxy (Variant 2) -- UQFF MUGE
// Host galaxy of Type Ia supernova SN2018gv, second analysis variant
// Source: grok_share_ba508f76c8e.txt entry #88 | Session 175
#include <vector>
#include <string>
#include <cmath>

class {cls} {{
public:
{UQFF_CONSTS}
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

    {cls}() = default;

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
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>

double {cls}::M_tot() const {{
    return M_stellar + M_DM;
}}

double {cls}::H_z() const {{
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + z_ngc, 3) + 0.7);
}}

// SN Ia feedback: energy impulse decaying over ~10 yr
// a_SN = L_SN / (c * M_tot * r_gal) * exp(-t/tau_SN)
double {cls}::a_SN(double t) const {{
    double a_peak = L_SN / (c * M_tot() * r_gal);
    return a_peak * std::exp(-t / tau_SN);
}}

// EM bar streaming with aether correction
double {cls}::a_EM_bar() const {{
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_bar * B_bar / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}}

// Full MUGE v2
double {cls}::g_NGC2525v2(double t) const {{
    double g_grav = G * M_tot() / (r_gal * r_gal);
    double g_core = g_grav * (1.0 + H_z() * t) * (1.0 + f_TRZ);
    return g_core + a_SN(t) + a_EM_bar();
}}

std::string {cls}::primary_equation() const {{
    return "g_NGC2525v2(t)=G*(M_s+M_DM)/r^2*(1+H_z*t)*(1+f_TRZ)+L_SN/(c*M*r)*exp(-t/tau_SN)+q*v_bar*B/(m_p)*(1+rho_UA/rho_SCm)*1e-12";
}}

void {cls}::self_update() {{
    curr_t += time_step;
}}

void {cls}::self_expand() {{
    r_gal *= 1.0001;
    M_DM  *= 1.00005;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        double g = g_NGC2525v2(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} n;
    std::cout << "NGC 2525 Barred Spiral v2 UQFF\\n";
    std::cout << n.primary_equation() << "\\n\\n";
    std::cout << "g(0) = " << n.g_NGC2525v2(0.0) << " m/s^2\\n";
    n.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_708: Pillars of Creation M16 UQFF
# ---------------------------------------------------------------------------
def gen_708():
    cls   = "PillarsOfCreationM16UQFF"
    guard = "PILLARS_OF_CREATION_M16_UQFF_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_708: Pillars of Creation (Eagle Nebula M16) -- Full UQFF MUGE
// Star-forming pillars 6500-7000 ly, 4-5 ly length, cool molecular H2
// Source: grok_share_ba508f76c8e.txt entry #86 | Session 175
#include <vector>
#include <string>
#include <cmath>

class {cls} {{
public:
{UQFF_CONSTS}
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

    {cls}() = default;

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
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>

double {cls}::M_dot(double t) const {{
    return M_0_SF * std::exp(-t / tau_SF);
}}

double {cls}::M_t(double t) const {{
    return M_initial * (1.0 + M_dot(t));
}}

// E(t) = E_0 * exp(-t/tau_erode) => 1-E at t=0.5 Myr ~ 0.9394
double {cls}::one_minus_E(double t) const {{
    return 1.0 - E_0 * std::exp(-t / tau_erode);
}}

double {cls}::H_z() const {{
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + z_M16, 3) + 0.7);
}}

double {cls}::a_EM() const {{
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_pillar / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}}

// g_Pillars = G*M(t)/r^2*(1+H_z*t)*(1-B/B_crit)*(1-E(t))
//           + (Ug1+Ug4)*(1+f_TRZ) + Lambda*c^2/3 + a_EM
// ~ 1.053e-4 m/s^2 at t=0.5 Myr
double {cls}::g_Pillars(double t) const {{
    double M   = M_t(t);
    double g_n = dpm_ug1_seed(M, r_pillar);
    const double B_crit = 1.0e11; // T
    double g_core = g_n * (1.0 + H_z() * t) * (1.0 - B_pillar / B_crit)
                        * one_minus_E(t);
    // Ug1 + Ug4 = 2 * Ug1 (approximately)
    double Ug1 = g_n;
    double Ug4 = Ug1 * (1.0 - B_pillar / B_crit);
    double Ug_sum = (Ug1 + Ug4) * (1.0 + f_TRZ);
    double g_lambda = Lambda * c * c / 3.0;
    return g_core + Ug_sum + g_lambda + a_EM();
}}

std::string {cls}::primary_equation() const {{
    return "g_Pillars(t)=G*M(t)/r^2*(1+H_z*t)*(1-B/Bc)*(1-E(t))+(Ug1+Ug4)*(1+f_TRZ)+Lambda*c^2/3+q*v*B/m*(1+rho_UA/rho_SCm)*1e-12";
}}

void {cls}::self_update() {{
    curr_t += time_step;
    M_initial *= one_minus_E(curr_t);
}}

void {cls}::self_expand() {{
    r_pillar *= 0.999; // pillar shortens as it erodes
    B_pillar *= 1.001;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        double g = g_Pillars(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} p;
    std::cout << "Pillars of Creation M16 UQFF\\n";
    std::cout << p.primary_equation() << "\\n\\n";
    std::cout << "g(0.5 Myr) = " << p.g_Pillars(1.578e13) << " m/s^2\\n";
    p.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_709: Westerlund 2 Star Cluster UQFF
# ---------------------------------------------------------------------------
def gen_709():
    cls   = "Westerlund2StarClusterUQFF"
    guard = "WESTERLUND2_STAR_CLUSTER_UQFF_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_709: Westerlund 2 Young Super Star Cluster -- Full UQFF MUGE
// Milky Way Carina, 15000 ly, ~3000 stars, 2 Myr age, M_stars~30000 M_sun
// Source: grok_share_ba508f76c8e.txt entry #85 | Session 175
#include <vector>
#include <string>
#include <cmath>

class {cls} {{
public:
{UQFF_CONSTS}
    // Westerlund 2 parameters
    double M_init     = 3.0e4  * 1.989e30; // kg initial stellar mass
    double r_wd2      = 9.461e16;           // m cluster radius (10 ly)
    double M_0_gas    = 3.333;              // fractional gas-to-star mass ratio
    double tau_SF     = 6.312e13;           // s SF timescale (2 Myr)
    double B_wd2      = 1.0e-5;             // T cluster magnetic field
    double v_gas      = 1.0e5;              // m/s gas velocity
    double q_p        = 1.602e-19;          // C proton charge
    double H0_wd2     = 2.184e-18;          // s^-1 H_0 for Westerlund calc

    // Self-simulation state
    double time_step  = 3.156e13;           // 1 Myr
    mutable double curr_t = 3.156e13;       // start at 1 Myr

    {cls}() = default;

    // M_dot: rapid star formation from initial gas reservoir
    double M_dot(double t) const;
    double M_t(double t) const;

    // EM aether term
    double a_EM() const;

    // Full Westerlund 2 MUGE
    // g_W2 = G*M(t)/r^2*(1+H0*t)*(1-B/B_crit)
    //       + (Ug1+Ug4)*(1+f_TRZ) + Lambda*c^2/3 + a_EM
    double g_Westerlund2(double t) const;

    std::string primary_equation() const;
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>

// M_dot: gas reservoir rapidly converts to stars
// M_dot(t) = M_0_gas * exp(-t/tau_SF)
// At t=1 Myr: M_dot = 3.333*0.6065 = 2.021 => M(t) = M_init * 3.021
double {cls}::M_dot(double t) const {{
    return M_0_gas * std::exp(-t / tau_SF);
}}

double {cls}::M_t(double t) const {{
    return M_init * (1.0 + M_dot(t));
}}

double {cls}::a_EM() const {{
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_wd2 / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}}

// Full Westerlund 2 MUGE
// g_W2 ~ 1.053e-3 m/s^2 at t=1 Myr
double {cls}::g_Westerlund2(double t) const {{
    double M   = M_t(t);
    double g_n = dpm_ug1_seed(M, r_wd2);
    const double B_crit = 1.0e11;
    double g_core = g_n * (1.0 + H0_wd2 * t) * (1.0 - B_wd2 / B_crit);
    double Ug1   = g_n;
    double Ug4   = Ug1 * (1.0 - B_wd2 / B_crit);
    double Ug_sum = (Ug1 + Ug4) * (1.0 + f_TRZ);
    double g_lambda = Lambda * c * c / 3.0;
    return g_core + Ug_sum + g_lambda + a_EM();
}}

std::string {cls}::primary_equation() const {{
    return "g_W2(t)=G*M(t)/r^2*(1+H0*t)*(1-B/Bc)+(Ug1+Ug4)*(1+f_TRZ)+Lambda*c^2/3+q*v*B/m*(1+rho_UA/rho_SCm)*1e-12";
}}

void {cls}::self_update() {{
    curr_t += time_step;
    M_init *= 1.0001; // slight mass gain from fallback
}}

void {cls}::self_expand() {{
    r_wd2 *= 1.001;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        double g = g_Westerlund2(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} wd2;
    std::cout << "Westerlund 2 Star Cluster UQFF\\n";
    std::cout << wd2.primary_equation() << "\\n\\n";
    std::cout << "g(1 Myr) = " << wd2.g_Westerlund2(3.156e13) << " m/s^2\\n";
    wd2.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_710: NGC 2014 - NGC 2020 Tapestry of Blazing Starbirth UQFF
# ---------------------------------------------------------------------------
def gen_710():
    cls   = "NGC2014NGC2020StarformingUQFF"
    guard = "NGC2014_NGC2020_STARFORMING_UQFF_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_710: NGC 2014 + NGC 2020 (Tapestry of Blazing Starbirth) UQFF
// LMC star-forming nebula pair, 163,000 ly, Wolf-Rayet + OB star dynamics
// Source: grok_share_ba508f76c8e.txt entry #84 | Session 175
#include <vector>
#include <string>
#include <cmath>

class {cls} {{
public:
{UQFF_CONSTS}
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

    {cls}() = default;

    double M_dot(double t) const;
    double M_t(double t) const;
    double a_EM() const;

    // Full LMC starbirth MUGE
    // g_TB = G*M(t)/r^2*(1+H0*t)*(1-B/Bc)
    //       + (Ug1+Ug4)*(1+f_TRZ) + Lambda*c^2/3 + a_EM
    double g_Starbirth(double t) const;

    std::string primary_equation() const;
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>

// Rapid secondary star formation from WR + OB stellar winds
// M_dot(t) = 41.67 * exp(-t/tau_SF); M(t=2.5 Myr) ~ 1.254e34 kg
double {cls}::M_dot(double t) const {{
    return M_0_SF * std::exp(-t / tau_SF);
}}

double {cls}::M_t(double t) const {{
    return M_initial * (1.0 + M_dot(t));
}}

double {cls}::a_EM() const {{
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_LMC / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}}

// g_Starbirth ~ 1.053e-4 m/s^2 at t=2.5 Myr
double {cls}::g_Starbirth(double t) const {{
    double M   = M_t(t);
    double g_n = dpm_ug1_seed(M, r_region);
    const double B_crit = 1.0e11;
    double g_core = g_n * (1.0 + H0_LMC * t) * (1.0 - B_LMC / B_crit);
    double Ug1 = g_n;
    double Ug4 = Ug1 * (1.0 - B_LMC / B_crit);
    double Ug_sum = (Ug1 + Ug4) * (1.0 + f_TRZ);
    double g_lambda = Lambda * c * c / 3.0;
    return g_core + Ug_sum + g_lambda + a_EM();
}}

std::string {cls}::primary_equation() const {{
    return "g_Starbirth(t)=G*M(t)/r^2*(1+H0*t)*(1-B/Bc)+(Ug1+Ug4)*(1+f_TRZ)+Lambda*c^2/3+q*v*B/m*(1+rho_UA/rho_SCm)*1e-12";
}}

void {cls}::self_update() {{
    curr_t += time_step;
}}

void {cls}::self_expand() {{
    r_region *= 1.001;
    M_initial *= 1.001;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        double g = g_Starbirth(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} tb;
    std::cout << "NGC 2014/2020 Tapestry of Blazing Starbirth UQFF\\n";
    std::cout << tb.primary_equation() << "\\n\\n";
    std::cout << "g(2.5 Myr) = " << tb.g_Starbirth(7.89e13) << " m/s^2\\n";
    tb.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_711: NGC 2014 - NGC 2020 Variant 2 UQFF
# ---------------------------------------------------------------------------
def gen_711():
    cls   = "NGC2014NGC2020Variant2UQFF"
    guard = "NGC2014_NGC2020_VARIANT2_UQFF_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_711: NGC 2014 + NGC 2020 Variant 2 -- WR-dominated UQFF Analysis
// Second analysis with Wolf-Rayet cone structure and oxygen gas at 11000 C
// Source: grok_share_ba508f76c8e.txt entry #103 | Session 175
#include <vector>
#include <string>
#include <cmath>

class {cls} {{
public:
{UQFF_CONSTS}
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

    {cls}() = default;

    double M_WR_t(double t) const;
    double a_WR_radiation() const;
    double a_EM_WR() const;

    // Full WR cone MUGE variant 2
    // g_v2 = G*M_WR(t)/r^2*(1+H0*t)*(1+f_TRZ) + a_WR_rad + a_EM_WR
    double g_WRcone(double t) const;

    std::string primary_equation() const;
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>

// WR mass decreases as wind ejects material
double {cls}::M_WR_t(double t) const {{
    return M_WR * (1.0 + M_0_WR * std::exp(-t / tau_SF_WR));
}}

// WR radiation pressure on cone gas
// a_rad = L_WR / (4*pi*r_cone^2 * c) * (rho_WR/m_H)
double {cls}::a_WR_radiation() const {{
    const double m_H = 1.67e-27;
    double flux = L_WR / (4.0 * M_PI * r_cone * r_cone * c);
    return flux * (rho_WR / m_H);
}}

double {cls}::a_EM_WR() const {{
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_WR / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}}

// WR cone MUGE variant 2
double {cls}::g_WRcone(double t) const {{
    double M  = M_WR_t(t);
    double g_n = dpm_ug1_seed(M, r_cone);
    double g_core = g_n * (1.0 + H_0 * t) * (1.0 + f_TRZ);
    return g_core + a_WR_radiation() + a_EM_WR();
}}

std::string {cls}::primary_equation() const {{
    return "g_WRcone(t)=G*M_WR(t)/r^2*(1+H0*t)*(1+f_TRZ)+L_WR/(4pi*r^2*c)*(rho_WR/m_H)+q*v*B/m*(1+rho_UA/rho_SCm)*1e-12";
}}

void {cls}::self_update() {{
    curr_t += time_step;
}}

void {cls}::self_expand() {{
    r_cone  *= 1.002;
    rho_WR  *= 0.998;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        double g = g_WRcone(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} wrc;
    std::cout << "NGC 2014/2020 WR-variant 2 UQFF\\n";
    std::cout << wrc.primary_equation() << "\\n\\n";
    std::cout << "g(2.5 Myr) = " << wrc.g_WRcone(7.89e13) << " m/s^2\\n";
    wrc.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_712: Pillars of Creation M16 Variant 2 UQFF
# ---------------------------------------------------------------------------
def gen_712():
    cls   = "PillarsOfCreationM16v2UQFF"
    guard = "PILLARS_OF_CREATION_M16_V2_UQFF_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_712: Pillars of Creation M16 Variant 2 -- Post-Supernova Shock Analysis
// With 450,000 mph protostar jets and shockwave disruption scenario
// Source: grok_share_ba508f76c8e.txt entry #99 | Session 175
#include <vector>
#include <string>
#include <cmath>

class {cls} {{
public:
{UQFF_CONSTS}
    // M16 Pillars v2 parameters -- supernova shockwave scenario
    double M_initial  = 1.01e4 * 1.989e30; // kg initial pillar mass
    double r_pillar   = 4.731e16;           // m  pillar height (5 ly)
    double M_0_SF     = 0.9901;             // secondary SF fraction
    double tau_SF     = 3.156e13;           // s  SF timescale
    double E_0_shock  = 0.15;               // fractional erosion (higher: SN shockwave)
    double tau_shock  = 1.893e14;           // s  shock dissipation (6000 yr)
    double v_jet      = 2.0e5;              // m/s protostar jet velocity (450k mph)
    double L_jet      = 1.0e28;             // W  jet luminosity
    double B_jet      = 2.0e-6;             // T  jet magnetic field  
    double q_p        = 1.602e-19;          // C  proton charge
    double z_M16      = 1.2e-3;             // redshift

    // Self-simulation state
    double time_step  = 1.578e13;
    mutable double curr_t = 1.578e13;

    {cls}() = default;

    double M_dot(double t) const;
    double M_t(double t) const;

    // SN shockwave erosion (1 - E_shock(t))
    double one_minus_E_shock(double t) const;

    double H_z() const;
    double a_jet() const;
    double a_EM_jet() const;

    // Full v2 MUGE with shock+jet terms
    double g_Pillars_v2(double t) const;

    std::string primary_equation() const;
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>

double {cls}::M_dot(double t) const {{
    return M_0_SF * std::exp(-t / tau_SF);
}}

double {cls}::M_t(double t) const {{
    return M_initial * (1.0 + M_dot(t));
}}

// Supernova shockwave disruption: stronger initial erosion
// E_shock(t) = E_0_shock * exp(-t/tau_shock)
double {cls}::one_minus_E_shock(double t) const {{
    return 1.0 - E_0_shock * std::exp(-t / tau_shock);
}}

double {cls}::H_z() const {{
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + z_M16, 3) + 0.7);
}}

// Protostar jet radiation pressure
// a_jet = L_jet / (c * M_t) ~ jet forcing per unit mass
double {cls}::a_jet() const {{
    return L_jet / (c * M_initial);
}}

double {cls}::a_EM_jet() const {{
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_jet * B_jet / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}}

// Full v2 MUGE: g_Pv2 = g_core*(1-E_shock) + (Ug1+Ug4)*(1+f_TRZ) + a_jet + a_EM
double {cls}::g_Pillars_v2(double t) const {{
    double M   = M_t(t);
    double g_n = dpm_ug1_seed(M, r_pillar);
    const double B_crit = 1.0e11;
    double g_core = g_n * (1.0 + H_z() * t) * (1.0 - B_jet / B_crit)
                        * one_minus_E_shock(t);
    double Ug1 = g_n;
    double Ug4 = Ug1;
    double Ug_sum = (Ug1 + Ug4) * (1.0 + f_TRZ);
    double g_lambda = Lambda * c * c / 3.0;
    return g_core + Ug_sum + g_lambda + a_jet() + a_EM_jet();
}}

std::string {cls}::primary_equation() const {{
    return "g_Pv2(t)=G*M(t)/r^2*(1+H*t)*(1-B/Bc)*(1-E_shock)+(Ug1+Ug4)*(1+f_TRZ)+Lambda*c^2/3+L_jet/(c*M)+q*v_jet*B/m*(1+rho_UA/rho_SCm)*1e-12";
}}

void {cls}::self_update() {{
    curr_t += time_step;
}}

void {cls}::self_expand() {{
    r_pillar *= 0.9995; // SN shock collapses pillars faster
    L_jet    *= 0.999;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        double g = g_Pillars_v2(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} pv2;
    std::cout << "Pillars of Creation M16 v2 (Post-SN) UQFF\\n";
    std::cout << pv2.primary_equation() << "\\n\\n";
    std::cout << "g(0.5 Myr) = " << pv2.g_Pillars_v2(1.578e13) << " m/s^2\\n";
    pv2.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_713: UQFF Knowledge Base KB19 -- THz Signals 1-50
# ---------------------------------------------------------------------------
def gen_713():
    cls   = "UQFFKnowledgeBaseKB19"
    guard = "UQFF_KNOWLEDGE_BASE_KB19_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_713: UQFF Knowledge Base 19 -- THz Q-Scope Signal Analysis (Sets 10-50)
// 50 THz signals captured from Earth's core q-scope, 1.2-1.3 THz resonance range
// Source: grok_share_ba508f76c8e.txt entry #83 | Session 175
#include <vector>
#include <string>
#include <cmath>

class {cls} {{
public:
{UQFF_CONSTS}
    // THz q-scope signal parameters (Sets 10-50, Signals 1-50)
    double f_THz_meas = 1.246e12;   // Hz   measured THz frequency (1.246 THz)
    double omega_THz  = 7.853e12;   // rad/s angular frequency (2*pi*f_THz)
    double dA_range   = 6.205;      // A    differential amperage range (+/-3.102 A)
    double V_pp_max   = 1.00;       // V    peak-to-peak voltage maximum (Ch1)
    double V_eff_max  = 0.35;       // V    effective (RMS) voltage maximum
    double Z_line     = 50.0;       // Ohm  line impedance
    double dT_step    = 13.0;       // s    time step between images
    double tau_flow   = 52.0;       // s    ACE/DCE flow reversal period
    int    n_signals  = 50;         // total signals in set

    // UQFF-integrated signal parameters
    double mu_j_THz   = 1.0e-30;   // J*m magnetic string coupling at THz scale
    double rho_vac_UA = 7.09e-36;  // J/m^3 [UA] vacuum density
    double rho_vac_SCm = 7.09e-37; // J/m^3 [SCm] vacuum density

    // Self-simulation state
    double time_step   = 13.0;      // 13 s per signal
    mutable double curr_t = 0.0;
    mutable int    curr_sig = 1;

    {cls}() = default;

    // Signal power from voltage and impedance
    // P = V_eff^2 / Z
    double signal_power(double V_eff) const;

    // U_m magnetic string dynamics at f_THz
    // U_m = mu_j_THz * omega_THz * (1 - exp(-gamma*t*cos(pi*t_n))) * P_SCm
    double U_m_THz(double t) const;

    // Energy density contribution from THz bundle
    // rho_E = P / (c * A_eff),  A_eff ~ 1e-4 m^2
    double rho_E_THz(double V_eff) const;

    // Phase inversion function (ACE/DCE flow reversal)
    // phi_inv(t) = cos(2*pi*t/tau_flow)
    double phi_inversion(double t) const;

    // Full UQFF THz bundle integral
    // I_THz = sum_j [U_m_j * phi_inversion(j*dT_step)]
    double I_THz_bundle(int N) const;

    std::string primary_equation() const;
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>
#include <vector>

// Signal power P = V_eff^2 / Z ~ 0.00245 W
double {cls}::signal_power(double V_eff) const {{
    return V_eff * V_eff / Z_line;
}}

// U_m at THz frequency
// gamma ~ kappa for UQFF calibration
double {cls}::U_m_THz(double t) const {{
    double gamma = kappa;
    double t_n   = 0.0;
    double P_SCm = rho_vac_SCm * 1.0e-3; // approximate phase space volume
    return mu_j_THz * omega_THz * (1.0 - std::exp(-gamma * t * std::cos(M_PI * t_n))) * P_SCm;
}}

// Energy density from THz bundle
double {cls}::rho_E_THz(double V_eff) const {{
    const double A_eff = 1.0e-4; // m^2 effective area
    return signal_power(V_eff) / (c * A_eff);
}}

// ACE/DCE flow reversal phase
double {cls}::phi_inversion(double t) const {{
    return std::cos(2.0 * M_PI * t / tau_flow);
}}

// Bundle integral over N signals
double {cls}::I_THz_bundle(int N) const {{
    double sum = 0.0;
    for (int j = 0; j < N; j++) {{
        double t_j = j * dT_step;
        sum += U_m_THz(t_j) * phi_inversion(t_j);
    }}
    return sum;
}}

std::string {cls}::primary_equation() const {{
    return "I_THz = sum_j[mu_j*omega_THz*(1-exp(-kappa*t*cos(pi*t_n)))*P_SCm]*cos(2pi*t/tau_flow); U_m=[U_m:SM_m / Ug1=UQFF_g+SM_g]^SCm";
}}

void {cls}::self_update() {{
    curr_t += time_step;
    curr_sig++;
    if (curr_sig > n_signals) curr_sig = 1;
}}

void {cls}::self_expand() {{
    // Extend analysis window
    n_signals += 10;
}}

void {cls}::simulate(int num_steps) {{
    std::cout << "THz Bundle Integral (KB19, " << n_signals << " signals):\\n";
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        double I = I_THz_bundle(n_signals);
        double phi = phi_inversion(t);
        std::cout << "Step " << step
                  << "  t=" << t << " s"
                  << "  I_THz=" << I
                  << "  phi=" << phi << "\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} kb19;
    std::cout << "UQFF KB19 -- THz Q-Scope Signals 1-50\\n";
    std::cout << kb19.primary_equation() << "\\n\\n";
    std::cout << "P(V_eff=0.35) = " << kb19.signal_power(0.35) << " W\\n";
    std::cout << "I_bundle(50) = " << kb19.I_THz_bundle(50) << "\\n";
    kb19.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_714: UQFF Knowledge Base KB18 -- THz Signals 41-50
# ---------------------------------------------------------------------------
def gen_714():
    cls   = "UQFFKnowledgeBaseKB18"
    guard = "UQFF_KNOWLEDGE_BASE_KB18_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_714: UQFF Knowledge Base 18 -- THz Q-Scope Signal Analysis (Set 50: Signals 41-50)
// 10 signals: 16:48:23-16:50:20 (117 s interval), 200ns/div, 500mV/div
// Source: grok_share_ba508f76c8e.txt entry #82 | Session 175
#include <vector>
#include <string>
#include <cmath>

class {cls} {{
public:
{UQFF_CONSTS}
    // Set 50 (Signals 41-50) parameters
    double f_THz        = 1.246e12;   // Hz   THz signal frequency
    double omega_THz    = 7.853e12;   // rad/s
    double time_div     = 200.0e-9;   // s    time per division
    double V_div        = 0.5;        // V    voltage per division
    double dA_range     = 6.205;      // A    differential amperage
    double dT_interv    = 13.0;       // s    inter-signal interval
    double total_time   = 117.0;      // s    set duration
    int    n_signals_set = 10;        // signals in Set 50

    // Observed amplitude data [Ch1_mV, Ch2_mV] for signals 41-50
    double V_ch1[10] = {0.60, 0.65, 0.60, 0.55, 0.50, 0.60, 0.55, 0.50, 0.50, 0.50};
    double V_ch2[10] = {0.35, 0.40, 0.35, 0.30, 0.35, 0.40, 0.35, 0.30, 0.35, 0.35};

    double rho_vac_UA  = 7.09e-36;
    double rho_vac_SCm = 7.09e-37;

    // Self-simulation state
    double time_step    = 13.0;
    mutable double curr_t = 0.0;
    mutable int curr_idx  = 0;

    {cls}() = default;

    // Signal energy from voltage and impedance (Z=50 Ohm)
    double signal_energy(int idx) const;

    // U_m contribution for signal i
    double U_m_signal(int idx, double t) const;

    // Phase analysis: detect flow reversals
    // Returns +1 (normal), -1 (inverted), 0 (chaotic)
    int flow_state(int idx) const;

    // Cumulative bundle sum for set
    double bundle_sum() const;

    std::string primary_equation() const;
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>
#include <algorithm>

// Signal energy E = V_eff^2 / Z * dT_interval
double {cls}::signal_energy(int idx) const {{
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    double V_eff = V_ch1[idx] / std::sqrt(2.0);
    const double Z = 50.0;
    return V_eff * V_eff / Z * dT_interv;
}}

// U_m for signal at index idx, time offset t
double {cls}::U_m_signal(int idx, double t) const {{
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    double A  = V_ch1[idx];
    double t_j = idx * dT_interv;
    double gamma = kappa;
    return mu_J * omega_THz * A * (1.0 - std::exp(-gamma * (t + t_j)));
}}

// Flow state detection from Ch2 inversion
// Signals 44,47-48: chaotic (0); 45,46,49-50: inverted (-1); else normal (+1)
int {cls}::flow_state(int idx) const {{
    // Chaotic: indices 3 (sig44), 6 (sig47), 7 (sig48)
    if (idx == 3 || idx == 6 || idx == 7) return 0;
    // Inverted: 4 (sig45), 5 (sig46), 8 (sig49), 9 (sig50)
    if (idx == 4 || idx == 5 || idx == 8 || idx == 9) return -1;
    return 1;
}}

// Bundle sum: weighted by energy and flow states
double {cls}::bundle_sum() const {{
    double sum = 0.0;
    for (int i = 0; i < n_signals_set; i++) {{
        int fs = flow_state(i);
        if (fs != 0) {{
            sum += signal_energy(i) * fs;
        }}
    }}
    return sum;
}}

std::string {cls}::primary_equation() const {{
    return "KB18_bundle = sum_i[V_eff_i^2/Z*dT*flow_state(i)]; U_m_i = mu_J*omega_THz*V_i*(1-exp(-kappa*(t+t_i)))";
}}

void {cls}::self_update() {{
    curr_t += time_step;
    curr_idx = (curr_idx + 1) % n_signals_set;
}}

void {cls}::self_expand() {{
    // No geometric expansion for KB; extend time window
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        int    idx = step % n_signals_set;
        double E   = signal_energy(idx);
        double Um  = U_m_signal(idx, t);
        int    fs  = flow_state(idx);
        std::cout << "Sig " << (idx + 41)
                  << "  t=" << (curr_t + idx * dT_interv) << " s"
                  << "  E=" << E << " J  U_m=" << Um
                  << "  flow=" << (fs > 0 ? "normal" : fs < 0 ? "inverted" : "chaotic") << "\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} kb18;
    std::cout << "UQFF KB18 -- THz Signals 41-50 Analysis\\n";
    std::cout << kb18.primary_equation() << "\\n\\n";
    std::cout << "Bundle sum = " << kb18.bundle_sum() << " J\\n";
    kb18.simulate(5);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_715: UQFF Knowledge Base KB17 -- THz Signals 31-40
# ---------------------------------------------------------------------------
def gen_715():
    cls   = "UQFFKnowledgeBaseKB17"
    guard = "UQFF_KNOWLEDGE_BASE_KB17_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_715: UQFF Knowledge Base 17 -- THz Q-Scope Signal Analysis (Set 40: Signals 31-40)
// 10 signals: 16:46:13-16:48:10 (117 s interval), 200ns/div, 500mV/div
// Source: grok_share_ba508f76c8e.txt entry #81 | Session 175
#include <vector>
#include <string>
#include <cmath>

class {cls} {{
public:
{UQFF_CONSTS}
    // Set 40 (Signals 31-40) parameters
    double f_THz         = 1.246e12;  // Hz   signal frequency
    double omega_THz     = 7.853e12;  // rad/s
    double time_div      = 200.0e-9;  // s    time per division
    double V_div         = 0.5;       // V    voltage per division
    double dA_range      = 6.205;     // A    differential amperage  
    double dT_interv     = 13.0;      // s    inter-signal interval
    double total_time    = 117.0;     // s    set duration
    int    n_signals_set = 10;        // signals in Set 40

    // Observed amplitudes [Ch1_V, Ch2_V] for signals 31-40
    double V_ch1[10] = {0.60, 0.65, 0.60, 0.55, 0.50, 0.60, 0.55, 0.50, 0.50, 0.50};
    double V_ch2[10] = {0.35, 0.40, 0.35, 0.30, 0.35, 0.40, 0.35, 0.30, 0.35, 0.35};

    double rho_vac_UA  = 7.09e-36;
    double rho_vac_SCm = 7.09e-37;

    // UQFF resonance contribution
    double omega_i     = 7.85e12;   // rad/s primary THz mode
    double kappa_THz   = 5.0e-4;   // calibration constant

    // Self-simulation state
    double time_step    = 13.0;
    mutable double curr_t = 0.0;
    mutable int curr_idx  = 0;

    {cls}() = default;

    // Signal power (W)
    double signal_power(int idx) const;

    // U_bi adjustment for q-scope signal i
    // U_bi_i = rho_vac_UA * omega_i * V_ch1[i] * phi_inv(t)
    double U_bi_signal(int idx, double t) const;

    // Ug1 thread strength from frequency clustering
    // Ug1_i = mu_J * omega_i * V_ch1[i]^2
    double Ug1_signal(int idx) const;

    // Phase inversion state
    // -1 = inverted (signals 35,36,39,40), 0 = chaotic (34,37,38), +1 = normal
    int flow_state(int idx) const;

    // Total Ug1 thread strength across set
    double total_Ug1_thread() const;

    std::string primary_equation() const;
{SELF_METHODS}
}};

#endif // {guard}
"""
    cpp = f"""\
// STANDALONE_{cls.upper()}
#include "{cls}.h"
#include <iostream>

double {cls}::signal_power(int idx) const {{
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    double V_eff = V_ch1[idx] / std::sqrt(2.0);
    return V_eff * V_eff / 50.0;
}}

// U_bi: buoyancy adjustment for q-scope observing Earth's core THz hole
// U_bi ~ rho_UA * omega_i * V_ch1 * cos(2pi*t/tau) 
double {cls}::U_bi_signal(int idx, double t) const {{
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    const double tau_flow = 52.0; // s
    double phi = std::cos(2.0 * M_PI * t / tau_flow);
    return rho_vac_UA * omega_i * V_ch1[idx] * phi;
}}

// Ug1 thread strength: mu_J * omega * A^2
double {cls}::Ug1_signal(int idx) const {{
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    double A = V_ch1[idx];
    return mu_J * omega_i * A * A;
}}

// Classify flow state for signals 31-40 (0-indexed)
// Chaotic: 3 (sig34), 6 (sig37), 7 (sig38)
// Inverted: 4 (sig35), 5 (sig36), 8 (sig39), 9 (sig40)
int {cls}::flow_state(int idx) const {{
    if (idx == 3 || idx == 6 || idx == 7) return 0;
    if (idx == 4 || idx == 5 || idx == 8 || idx == 9) return -1;
    return 1;
}}

// Total Ug1 thread strength across all 10 signals in set
double {cls}::total_Ug1_thread() const {{
    double sum = 0.0;
    for (int i = 0; i < n_signals_set; i++) {{
        sum += Ug1_signal(i);
    }}
    return sum;
}}

std::string {cls}::primary_equation() const {{
    return "Ug1_thread = sum_i[mu_J*omega_THz*V_i^2]; U_bi_i = rho_UA*omega_i*V_i*cos(2pi*t/tau_flow); [U_m:SM_m/Ug1=UQFF_g+SM_g]^SCm";
}}

void {cls}::self_update() {{
    curr_t += time_step;
    curr_idx = (curr_idx + 1) % n_signals_set;
}}

void {cls}::self_expand() {{
    // Extend UQFF resonance model
    kappa_THz *= 1.001;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; step++) {{
        double t = curr_t + step * time_step;
        int    idx = step % n_signals_set;
        double P   = signal_power(idx);
        double Ug1 = Ug1_signal(idx);
        double Ubi = U_bi_signal(idx, t);
        std::cout << "Sig " << (idx + 31)
                  << "  P=" << P << " W  Ug1=" << Ug1
                  << "  U_bi=" << Ubi
                  << "  flow=" << (flow_state(idx) > 0 ? "normal" : flow_state(idx) < 0 ? "inv" : "chaotic") << "\\n";
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} kb17;
    std::cout << "UQFF KB17 -- THz Signals 31-40 Analysis\\n";
    std::cout << kb17.primary_equation() << "\\n\\n";
    std::cout << "Total Ug1 thread = " << kb17.total_Ug1_thread() << " J*m\\n";
    kb17.simulate(5);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# Main: generate all 14 module pairs
# ---------------------------------------------------------------------------
GENERATORS = [
    gen_702, gen_703, gen_704, gen_705, gen_706, gen_707,
    gen_708, gen_709, gen_710, gen_711, gen_712, gen_713,
    gen_714, gen_715,
]

if __name__ == "__main__":
    count = 0
    for gen_fn in GENERATORS:
        cls, h_content, cpp_content = gen_fn()
        h_path   = os.path.join(ROOT, f"{cls}.h")
        cpp_path = os.path.join(ROOT, f"{cls}.cpp")
        with open(h_path,   "w", encoding="utf-8") as f:
            f.write(h_content)
        with open(cpp_path, "w", encoding="utf-8") as f:
            f.write(cpp_content)
        print(f"  Written: {cls}.h + {cls}.cpp")
        count += 2
    print(f"\nDone. {count} files created.")
