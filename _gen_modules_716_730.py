"""
Session 176 -- Generate C++ .h and .cpp standalone modules for PAPER_716-730.
Source: grok_share_ba508f76c8e.txt
15 new modules: KB1-KB6, KB8-KB16 (KB7 already done)
PAPER_716=KB1, 717=KB2, 718=KB3, 719=KB4, 720=KB5, 721=KB6,
          722=KB8, 723=KB9, 724=KB10, 725=KB11, 726=KB12,
          727=KB13, 728=KB14, 729=KB15, 730=KB16
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
# PAPER_716: UQFFKnowledgeBaseKB1 -- Red Dwarf Compression_D
# ---------------------------------------------------------------------------
def gen_716():
    cls   = "UQFFKnowledgeBaseKB1"
    guard = "UQFF_KNOWLEDGE_BASE_KB1_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_716: UQFF Knowledge Base 1 -- Red Dwarf Compression_D
// Source: grok_share_ba508f76c8e.txt entry #65 | Session 176
// Doc 43.d: Inertia Papers (pp 1-10), Aether-Superconductive Paper (pp 1-8),
//           Hydrogen Papers (pp 74-84)
// Physics: quantum wave function, universal inertia U_i, dipole U_g1,
//          superconductor U_g2, magnetic disk U_g3, Jeans mass, plasma wave,
//          bosonic energy, DE power decomposition
#include <string>
#include <complex>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class {cls} {{
public:
{UQFF_CONSTS}
    // Inertia Papers parameters
    double A_inertia    = 1.0;       // Wave amplitude (normalized)
    double k_inertia    = 1.0;       // Wave number (1/m)
    double omega_in     = 1.0;       // Angular frequency (rad/s)
    double alpha_in     = 0.1;       // Decay constant (1/m)
    double r0_in        = 1.0;       // Reference radius (m)
    double beta_in      = 1.0;       // Caduceus twist parameter
    double lambda_I     = 1.0;       // Inertial constant
    double omega_m      = 1.0;       // Magnetic angular freq (rad/s)
    double q_m          = 1e-10;     // Magnetic charge (A*m)
    double omega_i_val  = 1e-8;      // Inertial ang freq (rad/s)
    double t_n          = 0.0;       // Normalized time
    double F_RZ         = 0.01;      // Redshift factor
    double m_boson      = 1e-30;     // Boson mass (kg)
    double omega_r      = 1e12;      // Radial frequency (rad/s)
    double x_boson      = 1e-10;     // Displacement (m)
    double n_boson      = 0;         // Quantum number
    double mu_mag       = 1e-23;     // Magnetic moment (A*m^2)
    double B_mag        = 1e-9;      // Magnetic field (T)
    double psi_0        = 1.0;       // Initial wave function
    double E_g          = 1e-18;     // Gravitational energy (J)
    double G_i          = 1e-19;     // Inertial gravity (J)
    double C_j          = 1e-20;     // Coupling (J)
    double m0           = 1e-30;     // Rest mass energy (J)
    double Delta_t      = 1e-15;     // Time uncertainty (s)
    double E_AC         = 1.77e-66;  // AC component (J)
    double I_0          = 1.0;       // Initial current (A)
    double gamma_damp   = 0.1;       // Damping (1/s)
    double L_spark      = 1e-6;      // Inductance (H)
    double C_spark      = 1e-12;     // Capacitance (F)
    double T_jeans      = 100.0;     // Temperature (K)
    double mu_jeans     = 1.0;       // Mean molecular weight
    double rho_jeans    = 1e-20;     // Density (kg/m^3)
    double rho0_density = 1e-20;     // Central density
    double r0_density   = 1.0;       // Scale radius (m)

    // Aether-Superconductive parameters
    double sigma_aether = 1.0;       // Gaussian width (m)
    double m_theta      = 1;         // Azimuthal quantum number
    double phi_golden   = 1.6180339887; // Golden ratio
    double f0_freq      = 174.0;     // Base Solfeggio frequency (Hz)
    double I_dipole     = 1e-20;     // Dipole current (A)
    double A_dipole     = 1e-10;     // Dipole area (m^2)
    double omega_spin   = 1e-3;      // Spin frequency (rad/s)
    double H_aether     = 1e6;       // Aether field (A/m)
    double M_disk       = 1e-5;      // Magnetization (A/m)
    double r_disk       = 1.0;       // Disk radius (m)
    double S_k          = 1e-34;     // Spinner term (J*s)
    int    n_spinners   = 2;
    double T_mu_nu      = 1e-3;      // Tensor component (J/m^3)
    int    n_tensors    = 3;
    double r_milky      = 8e3 * 3.086e19; // Milky Way r (m)
    double T_milky      = 2.5e8 * 3.156e7; // Period (s, 250 Myr)

    // Hydrogen Papers parameters
    double E_0_hydro    = 1e-10;     // Base energy (J) = aether density
    double fac_control  = 1e-94;     // Control logic scaling
    double fac_reactor  = 1e-100;    // Reactor operations scaling
    double fac_adjust   = 1e-95;     // Plasma adjustment scaling
    double fac_star     = 1e-94;     // Star structure scaling
    double fac_extract  = 1e-69;     // Gas extraction scaling

    // Self-simulation state
    double time_step    = 1e-10;
    mutable double curr_t = 0.0;

    // --- Methods ------------------------------------------------------------
    // Inertia Papers
    std::complex<double> quantum_wave_function(double r, double t) const;
    double caduceus_coil_twist(double t) const;
    double universal_inertia(double t) const;
    double bosonic_energy() const;
    double magnetic_influence() const;
    double de_power_decomposition() const;
    double ac_current(double t) const;
    double spark_resonance() const;
    double jeans_mass() const;
    double density_profile(double r) const;

    // Aether-Superconductive Paper
    std::complex<double> rotating_wave_function(double r, double theta, double t) const;
    double frequency_pattern(int n) const;
    double dipole_moment_U_g1() const;
    double superconductor_field_U_g2() const;
    double magnetic_disk_U_g3(double r) const;
    double plasma_wave() const;
    double spinners_U_gi() const;
    double tensor_sum_U_g5() const;
    double milky_way_velocity() const;

    // Hydrogen Papers
    double control_logic_energy() const;
    double reactor_operations_energy() const;
    double plasma_adjustment_energy() const;
    double star_structure_energy() const;
    double gas_extraction_energy() const;

    // Master UQFF output
    double primary_equation() const;
    std::string description() const;
{SELF_METHODS}
}};
#endif // {guard}
"""

    cpp = f"""\
#include "{cls}.h"
#include <iostream>
#include <string>

/*
PAPER_716: UQFF Knowledge Base 1 -- Red Dwarf Compression_D (43.d)
Source: grok_share_ba508f76c8e.txt entry #65

Physics assimilated from Doc 43.d (Red Dwarf Compression_D_06May2025.docx):

=== Inertia Papers (pages 1-10) ===
  psi(r,theta,phi,t) = A * Y_lm * sin(k*r - omega*t)/r * exp(-alpha*|r-r0|)
  Solved: psi ~= 4.83e5
  U_i = lambda_I * (rho_SCm/rho_UA * omega_i * cos(pi*t_n) * (1 + F_RZ))
  Solved: U_i ~= 8.05e-80 J/m^6
  E_boson = 0.5*m*omega_r^2*x^2 + hbar*omega_r*(n+0.5)  ~= 9.825e-19 J
  M_Jeans = (5*k_B*T / (G*mu*m_H))^1.5 * (3/(4*pi*rho))^0.5  ~= 5.13e31 kg
  E_AC ~= 1.77e-66 J
  omega_spark = 1/sqrt(LC) ~= 1e9 rad/s

=== Aether-Superconductive Paper (pages 1-8) ===
  mu_dipole = I * A * omega_spin  ~= 1e-51 A*m^2
  U_g1 = mu_dipole * B
  B_super = mu_0 * H_aether ~= 1.257 T -> U_g2 = B_super^2/(2*mu_0) ~= 6.29e5 J/m^3
  U_g3: B_disk ~= -1e-7 T
  f_n = f_0 * phi^n  (Fibonacci/Solfeggio: 174->963 Hz)
  v_milky_way ~= 220 km/s (orbital)
  plasma wave: omega_plasma ~= 1.005e16 rad/s

=== Hydrogen Papers (pages 74-84) ===
  E_control ~= 1.10e-104 J, E_reactor ~= 9.56e-110 J
  E_adjustment ~= 8.28e-105 J, E_star ~= 2.21e-104 J
  E_extraction ~= 5.52e-79 J

self.version = "Session176"
*/

std::complex<double> {cls}::quantum_wave_function(double r, double t) const {{
    // psi = A * sin(k*r - omega*t)/r * exp(-alpha*|r-r0|)
    double amp = A_inertia * std::sin(k_inertia * r - omega_in * t) / r
                 * std::exp(-alpha_in * std::abs(r - r0_in));
    return std::complex<double>(amp, 0.0);
}}

double {cls}::caduceus_coil_twist(double t) const {{
    return beta_in * std::sin(omega_in * t);  // solved ~= 0.8415
}}

double {cls}::universal_inertia(double t) const {{
    // U_i = lambda_I * (rho_SCm/rho_UA * omega_i * cos(pi*t_n) * (1 + F_RZ))
    return lambda_I * (rho_SCm / rho_UA * omega_i_val
           * std::cos(M_PI * t_n) * (1.0 + F_RZ));
}}

double {cls}::bosonic_energy() const {{
    return 0.5 * m_boson * omega_r * omega_r * x_boson * x_boson
           + hbar * omega_r * (n_boson + 0.5);  // ~= 9.825e-19 J
}}

double {cls}::magnetic_influence() const {{
    return -mu_mag * B_mag;  // H_mag ~= -2.32e-32 J
}}

double {cls}::de_power_decomposition() const {{
    return E_AC;  // E_AC ~= 1.77e-66 J
}}

double {cls}::ac_current(double t) const {{
    return I_0 * std::sin(omega_in * t) * std::exp(-gamma_damp * t);
}}

double {cls}::spark_resonance() const {{
    return 1.0 / std::sqrt(L_spark * C_spark);  // ~= 1e9 rad/s
}}

double {cls}::jeans_mass() const {{
    double mH = 1.6726e-27;  // hydrogen mass (kg)
    return std::pow(5.0 * k_B * T_jeans / (G * mu_jeans * mH), 1.5)
           * std::pow(3.0 / (4.0 * M_PI * rho_jeans), 0.5);  // ~= 5.13e31 kg
}}

double {cls}::density_profile(double r) const {{
    return rho0_density * std::exp(-r / r0_density);
}}

std::complex<double> {cls}::rotating_wave_function(double r, double theta, double t) const {{
    // psi(r,theta,t) = A * exp(-r^2/(2*sigma^2)) * exp(i*(m*theta - omega*t))
    double gauss = A_inertia * std::exp(-r*r / (2.0*sigma_aether*sigma_aether));
    std::complex<double> phase(0.0, m_theta * theta - omega_in * t);
    return gauss * std::exp(phase);
}}

double {cls}::frequency_pattern(int n) const {{
    return f0_freq * std::pow(phi_golden, n);
}}

double {cls}::dipole_moment_U_g1() const {{
    double mu = I_dipole * A_dipole * omega_spin;  // ~= 1e-51 A*m^2
    return mu * B_mag;
}}

double {cls}::superconductor_field_U_g2() const {{
    double B_super = mu_0 * H_aether;  // ~= 1.257 T
    return B_super * B_super / (2.0 * mu_0);  // U_g2 ~= 6.29e5 J/m^3
}}

double {cls}::magnetic_disk_U_g3(double r) const {{
    double B_disk = -mu_0 * M_disk / (4.0 * M_PI * r*r*r);
    return B_disk * B_disk / (2.0 * mu_0);
}}

double {cls}::plasma_wave() const {{
    double omega_0 = 1e16, gamma_pl = 1e15;
    return std::sqrt(omega_0*omega_0 + gamma_pl*gamma_pl);  // ~= 1.005e16 rad/s
}}

double {cls}::spinners_U_gi() const {{
    return n_spinners * S_k;  // ~= 2.108e-34 J*s
}}

double {cls}::tensor_sum_U_g5() const {{
    return n_tensors * T_mu_nu;  // ~= 3.6e-3 J/m^3
}}

double {cls}::milky_way_velocity() const {{
    return 2.0 * M_PI * r_milky / T_milky;  // ~= 220 km/s
}}

double {cls}::control_logic_energy() const {{
    return E_0_hydro * fac_control;  // ~= 1.10e-104 J
}}

double {cls}::reactor_operations_energy() const {{
    return E_0_hydro * fac_reactor;  // ~= 9.56e-110 J
}}

double {cls}::plasma_adjustment_energy() const {{
    return E_0_hydro * fac_adjust;  // ~= 8.28e-105 J
}}

double {cls}::star_structure_energy() const {{
    return E_0_hydro * fac_star;  // ~= 2.21e-104 J
}}

double {cls}::gas_extraction_energy() const {{
    return E_0_hydro * fac_extract;  // ~= 5.52e-79 J
}}

double {cls}::primary_equation() const {{
    // Master UQFF output: combine U_g1, U_g2, U_i, E_extraction
    double Ug1 = dipole_moment_U_g1();
    double Ug2 = superconductor_field_U_g2();
    double Ui  = universal_inertia(curr_t);
    double E_h = gas_extraction_energy();
    return Ug1 + Ug2 + Ui + E_h;
}}

std::string {cls}::description() const {{
    return "PAPER_716: UQFF KB1 -- Red Dwarf Compression_D | "
           "Inertia/Aether-SC/Hydrogen-pp74-84 | Session176";
}}

void {cls}::self_update() {{
    curr_t += time_step;
}}

void {cls}::self_expand() {{
    rho0_density *= 1.001;
    sigma_aether *= 1.0005;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; ++step) {{
        double t = curr_t + step * time_step;
        double P = primary_equation();
        double Ui = universal_inertia(t);
        std::cout << "Step " << step
                  << "  primary=" << P
                  << "  U_i=" << Ui << "\\n";
        self_update();
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} kb1;
    std::cout << "UQFF KB1 -- Red Dwarf Compression_D Analysis\\n";
    std::cout << kb1.description() << "\\n";
    std::cout << "U_g2 = " << kb1.superconductor_field_U_g2() << " J/m^3\\n";
    std::cout << "U_i  = " << kb1.universal_inertia(0.0)       << " J/m^6\\n";
    std::cout << "Jeans mass = " << kb1.jeans_mass() << " kg\\n";
    std::cout << "Milky Way v = " << kb1.milky_way_velocity() << " m/s\\n";
    kb1.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_717: UQFFKnowledgeBaseKB2 -- Red Dwarf Compression_E
# ---------------------------------------------------------------------------
def gen_717():
    cls   = "UQFFKnowledgeBaseKB2"
    guard = "UQFF_KNOWLEDGE_BASE_KB2_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_717: UQFF Knowledge Base 2 -- Red Dwarf Compression_E
// Source: grok_share_ba508f76c8e.txt entry #66 | Session 176
// Doc 43.e: Hydrogen Papers pages 85-88
// Physics: compressed space dynamics E_space, Earth-Moon tidal analogy E(t),
//          hydrogen radial probability (n=1 to n=3), 26-level quantum wave E_k(t)
#include <string>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class {cls} {{
public:
{UQFF_CONSTS}
    // Compressed space / tidal parameters
    double E_0_aether   = 1.683e-37;  // J  aether base energy
    double E_aether_vol = 1.683e-10;  // J/m^3 * V_eff (energy density * volume)
    double B_pseudo     = 1.0;        // T  di-pseudo-monopole field
    double T_earth_moon = 2.36e6;     // s  Earth-Moon orbital period
    double P_tidal      = 3.2e12;     // J/s  Moon tidal power
    double V_eff        = 1.0;        // m^3  effective volume

    // Hydrogen quantum state parameters
    double a0 = 5.292e-11;   // m  Bohr radius
    int    n_max = 6;         // max quantum number for 26-level
    int    n_levels = 26;     // 26-level pattern

    // Spatial config factors (pp 85-86)
    double spatial_config  = 2.0;   // spherical
    double compression     = 1.0;
    double layers          = 5.0;
    double higgs_freq      = 8e-34; // Higgs frequency factor
    double precession_time = 6.183e-13; // Precession timing factor
    double quantum_scale   = 3.333e-23; // Quantum scaling

    // Self-simulation state
    double time_step   = 3.156e6;   // 1 month
    mutable double curr_t = 0.0;

    // --- Methods ------------------------------------------------------------
    double compressed_space_energy() const;
    double earth_moon_tidal_energy(double t) const;
    double radial_probability_energy(int n, int l, int m, double t) const;
    double quantum_wave_26level(int k, double t) const;
    double primary_equation() const;
    std::string description() const;
{SELF_METHODS}
}};
#endif // {guard}
"""

    cpp = f"""\
#include "{cls}.h"
#include <iostream>
#include <string>

/*
PAPER_717: UQFF Knowledge Base 2 -- Red Dwarf Compression_E (43.e)
Source: grok_share_ba508f76c8e.txt entry #66

=== Hydrogen Papers pages 85-88 ===
Page 85 (Compressed Space, Spherical):
  E_space = E_0 * Spatial * Compression * Layer * Higgs_freq * Precession * Q_scale
  ~= 5.52e-104 J

Page 86 (Compressed Space + Rotational + Earth-Moon tidal):
  E_space_rot ~= 1.10e-104 J (Rotational_Motion_Factor = 1)
  E(t) = E_aether_vol * (B_pseudo^2/(2*mu_0) * 1/E_aether) * sin(2*pi*t/T)
  E(T/4) ~= 7.96e-22 J (UQFF)  vs  E_SM(T/4) ~= 1.888e18 J

Page 87 (Hydrogen Radial Probability n=1->3):
  E(t) = E_aether_vol * (B_pseudo^2/(2*mu_0)/E_aether) * |psi_nlm|^2*r^2_max * sin(2*pi*t/T)
  E_1s(T/4) ~= 3.98e-22 J, E_3d(T/4) ~= 1.99e-22 J

Page 88 (26-Level Quantum Wave n=1->6):
  E_k(t) = E_aether_vol * (B_pseudo^2/(2*mu_0)/E_aether) * |Y_lm|^2 * sin(2*pi*t/T_k)
  T_k = k/26 * T_earth_moon
  E_1(T_1/4) ~= 5.31e-23 J, E_6(T_6/4) ~= 2.37e-22 J

Earth-Moon analogy: di-pseudo-monopoles (SCm:UA') unify atomic and cosmic scales
calibration factor UQFF/SM ~= 10^38-10^39

self.version = "Session176"
*/

double {cls}::compressed_space_energy() const {{
    // E_space ~= 5.52e-104 J (spherical configuration)
    return E_0_aether * spatial_config * compression * layers
           * higgs_freq * precession_time * quantum_scale;
}}

double {cls}::earth_moon_tidal_energy(double t) const {{
    // E(t) = E_aether_vol * (B^2/(2*mu_0)/E_aether) * sin(2*pi*t/T)
    double B2_term = B_pseudo * B_pseudo / (2.0 * mu_0);
    double E_aether = E_0_aether / V_eff;  // energy density
    return E_aether_vol * (B2_term / E_aether) * spatial_config
           * std::sin(2.0 * M_PI * t / T_earth_moon);
}}

double {cls}::radial_probability_energy(int n, int l, int m, double t) const {{
    // Simplified |psi_nlm|^2*r^2 peak ratio relative to 1s
    // n=1,l=0 -> ratio=1.0; n=3,l=2 -> ratio~=0.50
    double ratio = 1.0 / (double)(n * n);
    double B2_term = B_pseudo * B_pseudo / (2.0 * mu_0);
    double E_aether = E_0_aether / V_eff;
    return E_aether_vol * (B2_term / E_aether) * ratio
           * std::sin(2.0 * M_PI * t / T_earth_moon);
}}

double {cls}::quantum_wave_26level(int k, double t) const {{
    // E_k(t) = E_aether_vol * (B^2/(2*mu_0)/E_aether) * |Y_lm|^2 * sin(2pi*t/T_k)
    // T_k = k/26 * T_earth_moon, |Y_lm|^2 simplified as 1/(4*pi) for l=0
    double T_k = (double)k / n_levels * T_earth_moon;
    double Y2  = (k <= 5) ? 0.0796 : 0.596;  // Y_0,0 vs Y_2,+/-2
    double B2_term = B_pseudo * B_pseudo / (2.0 * mu_0);
    double E_aether = E_0_aether / V_eff;
    return E_aether_vol * (B2_term / E_aether) * Y2
           * std::sin(2.0 * M_PI * t / T_k);
}}

double {cls}::primary_equation() const {{
    // Sum E_space + tidal + 26-level representative
    double E_sp   = compressed_space_energy();
    double E_tid  = earth_moon_tidal_energy(T_earth_moon / 4.0);
    double E_26   = quantum_wave_26level(1, T_earth_moon / 26.0 / 4.0);
    return E_sp + E_tid + E_26;
}}

std::string {cls}::description() const {{
    return "PAPER_717: UQFF KB2 -- Red Dwarf Compression_E | "
           "Hydrogen-pp85-88 | 26-level quantum wave | Session176";
}}

void {cls}::self_update() {{
    curr_t += time_step;
}}

void {cls}::self_expand() {{
    n_levels = (n_levels < 52) ? n_levels + 1 : n_levels;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; ++step) {{
        double t = curr_t + step * time_step;
        double P = primary_equation();
        double E6 = quantum_wave_26level(6, t);
        std::cout << "Step " << step
                  << "  primary=" << P
                  << "  E_26lv6=" << E6 << "\\n";
        self_update();
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} kb2;
    std::cout << "UQFF KB2 -- Red Dwarf Compression_E Analysis\\n";
    std::cout << kb2.description() << "\\n";
    std::cout << "E_space = " << kb2.compressed_space_energy() << " J\\n";
    std::cout << "E_tidal(T/4) = " << kb2.earth_moon_tidal_energy(kb2.T_earth_moon/4.0) << " J\\n";
    std::cout << "E_26lv_1 = " << kb2.quantum_wave_26level(1,kb2.T_earth_moon/26.0/4.0) << " J\\n";
    kb2.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_718: UQFFKnowledgeBaseKB3 -- Red Dwarf Compression_C (LENR, Pi)
# ---------------------------------------------------------------------------
def gen_718():
    cls   = "UQFFKnowledgeBaseKB3"
    guard = "UQFF_KNOWLEDGE_BASE_KB3_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_718: UQFF Knowledge Base 3 -- Red Dwarf Compression_C
// Source: grok_share_ba508f76c8e.txt entry #67 | Session 176
// Doc 43.c: Electro-Weak LENR primer, Higgs/collider data, NGC 346 gas nebula,
//           Pi/Phi calculation notes, buoyancy equations
// Physics: LENR weak interaction, neutron production rate eta, Higgs U_H,
//          U_m master, U_g3 star formation, Pi series influence
#include <string>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class {cls} {{
public:
{UQFF_CONSTS}
    // LENR / Collider parameters
    double k_eta        = 1e-113; // Neutron production rate constant
    double n26          = 26.0;   // 26 quantum states
    double SSq_val      = 0.57;   // Superstring quenching [SSq]
    double omega_c      = 1.585e-8; // rad/s, UQFF cycle freq (2pi/3.96e8 s)
    double lambda_H     = 1.0;    // Higgs field constant
    double rho_UA_prime = 1e-23;  // J/m^3  UA' pseudo-monopole level
    double m_H_GeV      = 125.0;  // GeV  Higgs mass
    double k_Higgs      = 1.79e18; // Higgs calibration factor
    double f_Heaviside  = 0.01;
    double f_quasi      = 0.01;
    double P_SCm        = 1.0;
    double E_react0     = 1e46;   // J  reaction energy at t=0

    // Gas nebula / NGC 346 parameters
    double T_SF         = 1.424e6;  // K  star formation temperature
    double V_radial     = -3.33e-5; // blueshift velocity
    double V_buoy       = 0.0303;   // buoyancy ratio Va/Vb = 1/33

    // Pi/Phi constants
    double phi_golden   = 1.6180339887;
    double alpha_FS     = 1.0/137.0; // fine structure constant

    // Self-simulation state
    double time_step    = 3.156e7;  // 1 yr
    mutable double curr_t = 0.0;
    mutable double mu_j = 3.38e23;  // J*m, current dipole strength

    // --- Methods ------------------------------------------------------------
    double neutron_production_rate(double t) const;
    double higgs_field_U_H(double t) const;
    double universal_magnetism_U_m(double r_j, double t) const;
    double U_g3_star_formation(double r, double t) const;
    double pi_series_influence(int n_terms) const;
    double phi_series_influence() const;
    double buoyancy_gravity_sum(int n_terms) const;
    double primary_equation() const;
    std::string description() const;
{SELF_METHODS}
}};
#endif // {guard}
"""

    cpp = f"""\
#include "{cls}.h"
#include <iostream>
#include <string>

/*
PAPER_718: UQFF KB3 -- Red Dwarf Compression_C (43.c)
Source: grok_share_ba508f76c8e.txt entry #67

=== LENR Equations ===
  weak: W + e- + p -> n + nu_e, Q ~= 0.78 MeV
  eta = k_eta * exp(-SSq^n26 * exp(-pi-t)) * U_m / rho_UA
  Higgs: m_H = lambda_H * rho_UA_prime * omega_H * (1+f_quasi) -> 125 GeV

=== UQFF Framework ===
  U_m = sum_j [mu_j/r_j * (1 - exp(-gamma*t*cos(pi*t_n)))*phi_j]
              * P_SCm * E_react * (1 + 1e13*f_Heaviside) * (1 + f_quasi)
  U_H = lambda_H * rho_UA_prime * omega_H * exp(-SSq^n26*exp(-pi-t)) * (1+f_quasi)
  U_g3: T_star_formation ~= 1.424e6 K

=== Pi/Phi Equations ===
  Pi Influence ~ U_m * pi * rho_UA
  Phi Influence ~ U_m * phi * rho_UA
  FSC Influence ~ alpha * U_m
  Buoyancy ~ rho_UA / rho_SCm * V_little/V_big  (V_little/V_big ~= 1/33)

self.version = "Session176"
*/

double {cls}::neutron_production_rate(double t) const {{
    // eta = k_eta * exp(-SSq^n26 * exp(-pi-t)) * U_m/rho_UA
    double exponent = -std::pow(SSq_val, n26) * std::exp(-M_PI - t);
    double U_m = universal_magnetism_U_m(1.496e13, t);
    return k_eta * std::exp(exponent) * U_m / rho_UA;
}}

double {cls}::higgs_field_U_H(double t) const {{
    // U_H = lambda_H * rho_UA_prime * omega_c * exp(-SSq^n26*exp(-pi-t)) * (1+f_quasi)
    double exponent = -std::pow(SSq_val, n26) * std::exp(-M_PI - t);
    return lambda_H * rho_UA_prime * omega_c * std::exp(exponent) * (1.0 + f_quasi);
}}

double {cls}::universal_magnetism_U_m(double r_j, double t) const {{
    // Simplified single-term U_m
    double gamma_ = 5e-5;  // day^-1
    double t_n_   = 0.0;
    double E_react = E_react0 * std::exp(-kappa * t);
    double hat_phi = 1.0;
    double term = (mu_j / r_j * (1.0 - std::exp(-gamma_ * t * std::cos(M_PI * t_n_))) * hat_phi)
                  * P_SCm * E_react * (1.0 + 1e13 * f_Heaviside) * (1.0 + f_quasi);
    return term;
}}

double {cls}::U_g3_star_formation(double r, double t) const {{
    // U_g3 = k3 * sum_j B_j * cos(omega_s*t*pi) * P_core * E_react
    double B_j = mu_0 * mu_j / (4.0 * M_PI * r * r * r);
    double omega_s = 2.5e-6;
    double E_react = E_react0 * std::exp(-kappa * t);
    return 1.0 * B_j * std::cos(omega_s * t * M_PI) * 1.0 * E_react;
}}

double {cls}::pi_series_influence(int n_terms) const {{
    // Pi Influence ~ U_m * pi * rho_UA  (simplified sum)
    double U_m = universal_magnetism_U_m(1.496e13, 0.0);
    return U_m * M_PI * rho_UA * n_terms;
}}

double {cls}::phi_series_influence() const {{
    double U_m = universal_magnetism_U_m(1.496e13, 0.0);
    return U_m * phi_golden * rho_UA;
}}

double {cls}::buoyancy_gravity_sum(int n_terms) const {{
    // Buoyancy ~ rho_UA/rho_SCm * V_littl/V_big  (=1/33)
    double buoy = rho_UA / rho_SCm * V_buoy;
    double result = 0.0;
    for (int n = 1; n <= n_terms; n += 2) {{   // odd n
        double denom = 3.0 - std::pow(M_PI + 1.0, (double)n);
        if (std::abs(denom) > 1e-30) result += 1.0 / denom;
    }}
    return buoy * result;
}}

double {cls}::primary_equation() const {{
    double eta  = neutron_production_rate(0.0);
    double U_H  = higgs_field_U_H(0.0);
    double U_m  = universal_magnetism_U_m(1.496e13, 0.0);
    double U_g3 = U_g3_star_formation(3.086e20, 0.0);
    return eta + U_H + U_m + U_g3;
}}

std::string {cls}::description() const {{
    return "PAPER_718: UQFF KB3 -- Red Dwarf Compression_C | "
           "LENR/Higgs/Pi/Phi | NGC 346 gas nebula | Session176";
}}

void {cls}::self_update() {{
    curr_t += time_step;
    mu_j *= 1.0001;
}}

void {cls}::self_expand() {{
    rho_UA_prime *= 1.001;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; ++step) {{
        double t = curr_t + step * time_step;
        double P  = primary_equation();
        double Uh = higgs_field_U_H(t);
        std::cout << "Step " << step
                  << "  primary=" << P
                  << "  U_H=" << Uh << "\\n";
        self_update();
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} kb3;
    std::cout << "UQFF KB3 -- Red Dwarf Compression_C Analysis\\n";
    std::cout << kb3.description() << "\\n";
    std::cout << "U_H = " << kb3.higgs_field_U_H(0.0) << " J/m^3\\n";
    std::cout << "eta(t=0) = " << kb3.neutron_production_rate(0.0) << "\\n";
    std::cout << "Pi influence = " << kb3.pi_series_influence(5) << "\\n";
    kb3.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_719: UQFFKnowledgeBaseKB4 -- Red Dwarf Compression_B (Drawing 32/33)
# ---------------------------------------------------------------------------
def gen_719():
    cls   = "UQFFKnowledgeBaseKB4"
    guard = "UQFF_KNOWLEDGE_BASE_KB4_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_719: UQFF Knowledge Base 4 -- Red Dwarf Compression_B
// Source: grok_share_ba508f76c8e.txt entry #68 | Session 176
// Doc 43.b: Drawing 32 (nebular cloud photo), Drawing 33 (shock star formation)
// Physics: U_g4 for black hole formation, U_g4 for shock-induced SF,
//          geometric star distances and angles, vacuum energy densities
#include <string>
#include <cmath>
#include <vector>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct StarPos {{ double x, y; }};

class {cls} {{
public:
{UQFF_CONSTS}
    // Drawing 32: nebular cloud BH parameters
    double rho_SCm_nebula = 2.39e-22;  // J/m^3  level 13 nebula
    double M_BH_nebula    = 1.989e36;  // kg  black hole mass
    double d_g_nebula     = 3.086e16;  // m  ~1 pc
    double alpha_decay    = 0.001;     // decay constant
    double f_feedback     = 0.1;

    // Drawing 33: shock star formation parameters
    double M_star_SF      = 1.989e30;  // kg  star mass
    double d_g_SF         = 1.496e14;  // m  ~1 AU
    double f_shock        = 0.1;

    // Geometric star positions (Drawing 32, normalized units)
    std::vector<StarPos> star_positions = {{
        {{100, 900}},  // Star 1
        {{500, 900}},  // Star 2
        {{900, 900}},  // Star 3
        {{500, 100}},  // Star 4
        {{200, 100}}   // Star 5
    }};

    // Self-simulation state
    double time_step    = 3.156e7;  // 1 yr
    mutable double curr_t = 0.0;

    // --- Methods ------------------------------------------------------------
    double U_g4_nebular_BH(double t) const;
    double U_g4_shock_SF(double t) const;
    double geometric_distance(int i, int j) const;
    double geometric_angle(int apex, int from, int to) const;
    double primary_equation() const;
    std::string description() const;
{SELF_METHODS}
}};
#endif // {guard}
"""

    cpp = f"""\
#include "{cls}.h"
#include <iostream>
#include <string>

/*
PAPER_719: UQFF KB4 -- Red Dwarf Compression_B (43.b)
Source: grok_share_ba508f76c8e.txt entry #68

=== Drawing 32: Nebular Cloud Photo ===
  U_g4 = k4 * rho_SCm * M_BH / d_g * exp(-alpha*t) * cos(pi*t_n) * (1+f_feedback)
  U_g4 ~= 1.69e-2 * exp(-0.001t) * cos(pi*t) J/m^3
  Star positions: (100,900), (500,900), (900,900), (500,100), (200,100)
  Distances: d_12=400, d_23=400, d_13=800, d_24=800, d_45=300
  Angles: theta_{1-2-3}=180 deg, theta_{1-2-4}=90 deg, theta_{2-4-5}=90 deg
  rho_SCm = 2.39e-22 J/m^3 (nebula scale level 13)

=== Drawing 33: Shock-Induced Star Formation ===
  U_g4 = k4 * rho_SCm * M_star / d_g * exp(-alpha*t) * cos(pi*t_n) * (1+f_shock)
  U_g4 ~= 3.49e-6 * exp(-0.001t) * cos(pi*t) J/m^3
  SiO / formamide astrochemical tracers, protostellar jets

=== LENR/Collider (referenced from 43.b) ===
  Same LENR + Higgs equations as 43.c (W+e-+p->n+nu_e, Higgs ~125 GeV)

self.version = "Session176"
*/

double {cls}::U_g4_nebular_BH(double t) const {{
    // U_g4 = k4 * rho_SCm * M_BH / d_g * exp(-alpha*t) * cos(pi*t_n) * (1+f_feedback)
    double t_n = 0.0;
    return 1.0 * rho_SCm_nebula * M_BH_nebula / d_g_nebula
           * std::exp(-alpha_decay * t) * std::cos(M_PI * t_n) * (1.0 + f_feedback);
}}

double {cls}::U_g4_shock_SF(double t) const {{
    // U_g4 for shock-induced star formation (Drawing 33)
    double t_n = 0.0;
    return 1.0 * rho_SCm_nebula * M_star_SF / d_g_SF
           * std::exp(-alpha_decay * t) * std::cos(M_PI * t_n) * (1.0 + f_shock);
}}

double {cls}::geometric_distance(int i, int j) const {{
    if (i < 0 || i >= (int)star_positions.size()) return 0.0;
    if (j < 0 || j >= (int)star_positions.size()) return 0.0;
    double dx = star_positions[j].x - star_positions[i].x;
    double dy = star_positions[j].y - star_positions[i].y;
    return std::sqrt(dx*dx + dy*dy);
}}

double {cls}::geometric_angle(int apex, int from, int to) const {{
    // Angle at vertex 'apex' in the triangle from-apex-to
    if (apex < 0 || from < 0 || to < 0) return 0.0;
    if (apex >= (int)star_positions.size()) return 0.0;
    double ax = star_positions[from].x - star_positions[apex].x;
    double ay = star_positions[from].y - star_positions[apex].y;
    double bx = star_positions[to].x   - star_positions[apex].x;
    double by = star_positions[to].y   - star_positions[apex].y;
    double dot = ax*bx + ay*by;
    double mag = std::sqrt(ax*ax+ay*ay) * std::sqrt(bx*bx+by*by);
    if (mag < 1e-30) return 0.0;
    double cosA = dot / mag;
    if (cosA >  1.0) cosA =  1.0;
    if (cosA < -1.0) cosA = -1.0;
    return std::acos(cosA) * 180.0 / M_PI;
}}

double {cls}::primary_equation() const {{
    double Ug4_bh = U_g4_nebular_BH(curr_t);
    double Ug4_sf = U_g4_shock_SF(curr_t);
    return Ug4_bh + Ug4_sf;
}}

std::string {cls}::description() const {{
    return "PAPER_719: UQFF KB4 -- Drawing 32 nebular BH + Drawing 33 shock SF | Session176";
}}

void {cls}::self_update() {{
    curr_t += time_step;
}}

void {cls}::self_expand() {{
    f_feedback *= 1.001;
    f_shock    *= 1.001;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; ++step) {{
        double t = curr_t + step * time_step;
        double P = primary_equation();
        std::cout << "Step " << step
                  << "  U_g4_BH=" << U_g4_nebular_BH(t)
                  << "  U_g4_SF=" << U_g4_shock_SF(t)
                  << "  total=" << P << "\\n";
        self_update();
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} kb4;
    std::cout << "UQFF KB4 -- Red Dwarf Compression_B Analysis\\n";
    std::cout << kb4.description() << "\\n";
    std::cout << "U_g4_BH(t=0) = " << kb4.U_g4_nebular_BH(0.0) << " J/m^3\\n";
    std::cout << "U_g4_SF(t=0) = " << kb4.U_g4_shock_SF(0.0)   << " J/m^3\\n";
    std::cout << "d_12 = " << kb4.geometric_distance(0,1) << " (normalized)\\n";
    std::cout << "angle_1-2-3 = " << kb4.geometric_angle(1,0,2) << " deg\\n";
    kb4.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_720: UQFFKnowledgeBaseKB5 -- Doc 43 (Universal Permanence, Final Parsec)
# ---------------------------------------------------------------------------
def gen_720():
    cls   = "UQFFKnowledgeBaseKB5"
    guard = "UQFF_KNOWLEDGE_BASE_KB5_H"
    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_720: UQFF Knowledge Base 5 -- Doc 43 (Red Dwarf Compression)
// Source: grok_share_ba508f76c8e.txt entry #69 | Session 176
// Doc 43 (UFE ORB EXP 2_24_07Mar2025): Universal Permanence (UP) equation,
//   Red Dwarf Reactor plasma orb, Drawing 30 (Bearden TRZ), Drawing 31 (AGN feedback),
//   LENR, Higgs, final parsec problem
// Physics: UP(t), non-local jump probability P, energy density rho_react,
//          metal retention f_Z, CGM baryon fraction f_CGM, U_g4 AGN/Final Parsec
#include <string>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class {cls} {{
public:
{UQFF_CONSTS}
    // UP equation / reactor parameters
    double t_n_val      = 13.68;     // s  normalized time
    double SCm_val      = 1e15;      // kg/m^3  SCm density
    double UA_val       = 1e-11;     // C   UA field
    double gamma_UP     = 1e3;       // s^-1  jump rate
    double NN_jumps     = 1.5;       // non-local jumps per frame
    double rho_react0   = 1e15;      // W/m^3  initial energy density
    double alpha_react  = 0.001;     // decay constant

    // Drawing 31 (AGN feedback / final parsec)
    double M_BH_agn     = 1.989e36;  // kg  ~10^6 Msun SMBH
    double M_BH_fp      = 8.15e36;   // kg  final parsec SMBH binary
    double d_g_agn      = 7.64e20;   // m   AGN distance
    double d_g_fp       = 2.55e20;   // m   final parsec separation
    double f_fbk_agn    = 0.1;       // AGN feedback factor
    double alpha_fp     = 0.001;

    // Metal retention / CGM (Drawing 31, Sanchez et al.)
    double f_Z_over     = 0.89;   // over-massive SMBH metal retention
    double f_Z_under    = 0.85;   // under-massive
    double f_Z_SF       = 0.73;   // star-forming disk
    double f_Z_quench   = 0.51;   // quenched disk

    // Self-simulation state
    double time_step    = 3.156e7;   // 1 yr
    mutable double curr_t = 0.0;

    // --- Methods ------------------------------------------------------------
    double universal_permanence_UP(double t) const;
    double non_local_jump_prob(double t) const;
    double energy_density_react(double t) const;
    double U_g4_agn_feedback(double t) const;
    double U_g4_final_parsec(double t) const;
    double primary_equation() const;
    std::string description() const;
{SELF_METHODS}
}};
#endif // {guard}
"""

    cpp = f"""\
#include "{cls}.h"
#include <iostream>
#include <string>

/*
PAPER_720: UQFF KB5 -- Doc 43 (Red Dwarf Compression, UFE ORB EXP 2_24_07Mar2025)
Source: grok_share_ba508f76c8e.txt entry #69

=== Universal Permanence (UP) Equation ===
  t^- = -t_n * exp(pi - t_n)  ~= -3.75e-4 s (at t_n=13.68 s)
  UP(t) = sum_i[k_i*U_gi] + sum_j[mu_j/r_j*(1-exp(-gamma*t-)...] + TRZ terms + NN + QS + ...
  UP ~= 1.03e20 * QS
  Energy density: rho_react = 1e15 * exp(-0.001*t_n) ~= 9.86e14 W/m^3

=== Non-Local Jump Probability ===
  P = 1 - exp(-gamma * |t^-|)  ~= 0.490/s (observed ~1-1.5 jumps/frame)

=== Drawing 31 (Sanchez et al. AGN Feedback) ===
  f_Z (over-massive SMBH) ~= 0.89, f_Z (under-massive) ~= 0.85
  f_Z (star-forming disk) ~= 0.73, f_Z (quenched) ~= 0.51
  U_g4_AGN ~= 8.40e-6 * exp(-0.001t) * cos(pi*t) J/m^3

=== Final Parsec Problem ===
  U_g4_FP ~= 7.64e-6 * exp(-0.001t) * cos(pi*t) J/m^3
  M_BH = 8.15e36 kg, d_g = 2.55e20 m

=== Drawing 30 (Bearden TRZ) ===
  U_m with time-reversal zones (negentropic processes)
  U_i with f_TRZ for plasmoid stability

self.version = "Session176"
*/

double {cls}::universal_permanence_UP(double t) const {{
    // t^- = -t_n * exp(pi - t_n)
    double t_neg = -t_n_val * std::exp(M_PI - t_n_val);
    // Simplified UP = gamma suppressed jump + energy density contribution
    double E_react = energy_density_react(t);
    double Um_term = (mu_J / 1.496e13)
                     * (1.0 - std::exp(-gamma_UP * std::abs(t_neg) * std::cos(M_PI * t_n_val)));
    return Um_term * E_react;
}}

double {cls}::non_local_jump_prob(double t) const {{
    // P = 1 - exp(-gamma * |t^-|)
    double t_neg = -t_n_val * std::exp(M_PI - t_n_val);
    (void)t;
    return 1.0 - std::exp(-gamma_UP * std::abs(t_neg));  // ~= 0.490
}}

double {cls}::energy_density_react(double t) const {{
    return rho_react0 * std::exp(-alpha_react * t);  // ~= 9.86e14 W/m^3
}}

double {cls}::U_g4_agn_feedback(double t) const {{
    // U_g4 = k4 * rho_SCm * M_BH / d_g * exp(-alpha*t) * cos(pi*t_n) * (1+f_feedback)
    double t_n_ = 0.0;
    return 1.0 * rho_SCm * M_BH_agn / d_g_agn
           * std::exp(-alpha_fp * t) * std::cos(M_PI * t_n_) * (1.0 + f_fbk_agn);
}}

double {cls}::U_g4_final_parsec(double t) const {{
    double t_n_ = 0.0;
    return 1.0 * rho_SCm * M_BH_fp / d_g_fp
           * std::exp(-alpha_fp * t) * std::cos(M_PI * t_n_);
}}

double {cls}::primary_equation() const {{
    double UP  = universal_permanence_UP(curr_t);
    double Ug4a = U_g4_agn_feedback(curr_t);
    double Ug4f = U_g4_final_parsec(curr_t);
    return UP + Ug4a + Ug4f;
}}

std::string {cls}::description() const {{
    return "PAPER_720: UQFF KB5 -- Universal Permanence + AGN Feedback + Final Parsec | Session176";
}}

void {cls}::self_update() {{
    curr_t += time_step;
}}

void {cls}::self_expand() {{
    f_fbk_agn *= 1.001;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; ++step) {{
        double t = curr_t + step * time_step;
        std::cout << "Step " << step
                  << "  UP=" << universal_permanence_UP(t)
                  << "  P_jump=" << non_local_jump_prob(t)
                  << "  U_g4_AGN=" << U_g4_agn_feedback(t) << "\\n";
        self_update();
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} kb5;
    std::cout << "UQFF KB5 -- Universal Permanence Analysis\\n";
    std::cout << kb5.description() << "\\n";
    std::cout << "P_jump = "   << kb5.non_local_jump_prob(0.0) << "/s\\n";
    std::cout << "rho_react = " << kb5.energy_density_react(13.68) << " W/m^3\\n";
    std::cout << "U_g4_AGN = " << kb5.U_g4_agn_feedback(0.0) << " J/m^3\\n";
    std::cout << "U_g4_FP  = " << kb5.U_g4_final_parsec(0.0) << " J/m^3\\n";
    kb5.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# Helper: Generate quantum-variable KB modules (KB6, KB8-KB13)
# ---------------------------------------------------------------------------
def _gen_qvar_kb(paper_num, cls, guard, title, entry_num, doc_desc,
                 variables, compute_notes, g_approx, session="176"):
    """Generic generator for quantum-variable KB entries."""
    # Build variable declarations
    var_decls = "\n".join(
        f"    double var_{v['name']} = {v['value']};  // {v['unit']}  {v['desc']}"
        for v in variables
    )
    # Build variables init in constructor body
    var_push = "\n".join(
        f'    vars.push_back({{"{v["name"]}", {v["value"]}, '
        f'"{v["desc"]}", "{v["eq"]}"}}); '
        for v in variables
    )

    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_{paper_num}: {title}
// Source: grok_share_ba508f76c8e.txt entry #{entry_num} | Session {session}
// {doc_desc}
#include <string>
#include <vector>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct QVarEntry {{
    std::string name;
    double value;
    std::string description;
    std::string equation;
}};

class {cls} {{
public:
{UQFF_CONSTS}
    // Quantum variable values from document analysis
{var_decls}

    // Variable registry
    std::vector<QVarEntry> vars;

    // UQFF core parameters
    double E_react0    = 1e46;   // J
    double gamma_core  = 5e-5;   // day^-1
    double f_Heaviside = 0.01;
    double f_quasi     = 0.01;
    double P_SCm_core  = 1.0;

    // Self-simulation state
    double time_step   = 3.156e7;  // 1 yr
    mutable double curr_t = 0.0;
    double expansion_factor = 1.01;

    // --- Methods ------------------------------------------------------------
    double compute_U_m(double mu_j, double r_j, double t) const;
    double compute_E_react(double t) const;
    double compute_U_g2(double r, double t) const;
    double primary_equation() const;
    std::string description() const;
{SELF_METHODS}
    {cls}();
}};
#endif // {guard}
"""

    cpp = f"""\
#include "{cls}.h"
#include <iostream>
#include <string>

/*
PAPER_{paper_num}: {title}
Source: grok_share_ba508f76c8e.txt entry #{entry_num}
{doc_desc}

{compute_notes}

g_result ~= {g_approx}
self.version = "Session{session}"
*/

{cls}::{cls}() : curr_t(0.0) {{
    // Populate quantum variable registry from document analysis
    {var_push}
}}

double {cls}::compute_U_m(double mu_j, double r_j, double t) const {{
    double E_react = compute_E_react(t);
    double hat_phi = 1.0;
    double t_n_ = 0.0;
    return (mu_j / r_j * (1.0 - std::exp(-gamma_core * t * std::cos(M_PI * t_n_))) * hat_phi)
           * P_SCm_core * E_react * (1.0 + 1e13 * f_Heaviside) * (1.0 + f_quasi);
}}

double {cls}::compute_E_react(double t) const {{
    return E_react0 * std::exp(-kappa * t);
}}

double {cls}::compute_U_g2(double r, double t) const {{
    double k2 = 1.2;
    double E_react = compute_E_react(t);
    return k2 * (rho_UA + rho_SCm) * M_sun / (r * r) * E_react;
}}

double {cls}::primary_equation() const {{
    return compute_U_m(mu_J, 1.496e13, curr_t)
         + compute_U_g2(3.086e20, curr_t)
         + compute_E_react(curr_t);
}}

std::string {cls}::description() const {{
    return "PAPER_{paper_num}: {title} | Session{session}";
}}

void {cls}::self_update() {{
    curr_t += time_step;
    if (!vars.empty()) vars.back().value *= 1.0001;
}}

void {cls}::self_expand() {{
    if (!vars.empty()) {{
        QVarEntry nv = vars.back();
        nv.value *= expansion_factor;
        vars.push_back(nv);
    }}
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; ++step) {{
        double t = curr_t + step * time_step;
        double P = primary_equation();
        double U_m = compute_U_m(mu_J, 1.496e13, t);
        std::cout << "Step " << step
                  << "  primary=" << P
                  << "  U_m=" << U_m << "\\n";
        self_update();
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} kb;
    std::cout << "{title}\\n";
    std::cout << kb.description() << "\\n";
    std::cout << "Variables loaded: " << kb.vars.size() << "\\n";
    for (auto& v : kb.vars)
        std::cout << "  " << v.name << " = " << v.value << "  [" << v.description << "]\\n";
    std::cout << "primary_equation = " << kb.primary_equation() << "\\n";
    kb.simulate(3);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_721: UQFFKnowledgeBaseKB6 -- Variables: r_j, d_g, F_U, f_feedback, Omega_g
# ---------------------------------------------------------------------------
def gen_721():
    return _gen_qvar_kb(
        paper_num=721, cls="UQFFKnowledgeBaseKB6",
        guard="UQFF_KNOWLEDGE_BASE_KB6_H",
        title="UQFF Knowledge Base 6 -- Quantum Variables r_j d_g F_U f_feedback Omega_g",
        entry_num=70,
        doc_desc="5 quantum variable documents: r_j, d_g, F_U, f_feedback, Omega_g",
        variables=[
            {"name": "r_j",       "value": "1.496e13", "unit": "m",
             "desc": "Distance along magnetic string to source j",
             "eq": "U_m = sum_j [mu_j/r_j * (1-exp(-gamma*t*cos(pi*t_n)))*phi_j] * P_SCm * E_react"},
            {"name": "d_g",       "value": "2.55e20",  "unit": "m",
             "desc": "Distance from Galactic Center",
             "eq": "U_bi = -beta_i*U_gi*Omega_g*(M_bh/d_g)*(1+epsilon_sw*rho_sw)*U_UA*cos(pi*t_n)"},
            {"name": "F_U",       "value": "2.28e65",  "unit": "N",
             "desc": "Unified Field Strength total",
             "eq": "F_U = sum_i[k_i*U_gi - beta_i*U_gi*Omega_g*(M_bh/d_g)*E_react] + ..."},
            {"name": "f_feedback","value": "0.1",       "unit": "dimensionless",
             "desc": "AGN feedback factor",
             "eq": "U_g4 = k4*(rho_SCm*M_bh/d_g)*exp(-alpha*t)*cos(pi*t_n)*(1+f_feedback)"},
            {"name": "Omega_g",   "value": "7.3e-16",  "unit": "rad/s",
             "desc": "Galactic angular spin rate",
             "eq": "U_bi = -beta_i*U_gi*Omega_g*(M_bh/d_g)*(1+epsilon_sw*rho_sw)*U_UA*cos(pi*t_n)"},
        ],
        compute_notes=
"""Variables: r_j=1.496e13 m (magnetic string distance, 1 AU),
   d_g=2.55e20 m (Galactic Center distance ~8 kpc),
   F_U=2.28e65 N (unified field strength total),
   f_feedback=0.1 (AGN feedback factor),
   Omega_g=7.3e-16 rad/s (Milky Way galactic spin rate)

U_m uses r_j for monopole string coupling.
U_bi uses d_g and Omega_g for buoyancy-gravity coupling to SMBH.
F_U sums all UQFF field contributions.
U_g4 uses f_feedback for AGN driven dynamics.""",
        g_approx="1.03e20 * QS (Up) + U_bi contributions"
    )


# ---------------------------------------------------------------------------
# PAPER_722: UQFFKnowledgeBaseKB8 -- Variables: M_bh, mu_j, P_core, t_n, pi
# ---------------------------------------------------------------------------
def gen_722():
    return _gen_qvar_kb(
        paper_num=722, cls="UQFFKnowledgeBaseKB8",
        guard="UQFF_KNOWLEDGE_BASE_KB8_H",
        title="UQFF Knowledge Base 8 -- Quantum Variables M_bh mu_j P_core t_n pi",
        entry_num=72,
        doc_desc="5 quantum variable documents: M_bh, mu_j, P_core, t_n, pi",
        variables=[
            {"name": "M_bh",   "value": "1.989e36",  "unit": "kg",
             "desc": "Black hole mass (10^6 Msun SMBH)",
             "eq": "U_g4 = k4*(rho_SCm*M_bh/d_g)*exp(-alpha*t)*cos(pi*t_n)*(1+f_feedback)"},
            {"name": "mu_j",   "value": "3.38e23",   "unit": "J*m",
             "desc": "Magnetic string dipole moment coupling for source j",
             "eq": "U_m = sum_j[mu_j/r_j*(1-exp(-gamma*t*cos(pi*t_n)))*phi_j]*P_SCm*E_react"},
            {"name": "P_core", "value": "1.0",        "unit": "dimensionless",
             "desc": "Core pressure normalization factor",
             "eq": "U_g3 = k3*sum_j B_j(r,theta,t,rho_SCm)*cos(omega_s*t*pi)*P_core*E_react"},
            {"name": "t_n",    "value": "0.0",        "unit": "s",
             "desc": "Normalized time parameter (UQFF phase cycle)",
             "eq": "cos(pi*t_n) -- phase modulator in all U_m, U_gN, U_bi terms"},
            {"name": "pi_val", "value": "3.14159265358979", "unit": "dimensionless",
             "desc": "Pi constant (26 quantum states, Wolfram 312 digits)",
             "eq": "cos(pi*t_n), exp(-pi-t), SSq^n26 * exp(-pi-t)"},
        ],
        compute_notes=
"""Variables: M_bh=1.989e36 kg (SMBH 10^6 Msun),
   mu_j=3.38e23 J*m (magnetic string coupling),
   P_core=1.0 (core pressure normalization),
   t_n=0.0 (normalized phase time),
   pi=3.14159... (26 quantum state encoder)

M_bh drives U_g4 AGN feedback and Final Parsec.
mu_j enters U_m as magnetic dipole per string j.
P_core normalizes U_g3 star formation term.
t_n modulates phase cos(pi*t_n) across all UQFF terms.
pi encodes quantum states: SSq^n26 * exp(-pi-t).""",
        g_approx="U_g4 ~= 8.40e-6 J/m^3, U_m ~= dominant at r=1 AU"
    )


# ---------------------------------------------------------------------------
# PAPER_723: UQFFKnowledgeBaseKB9 -- Variables: gamma, E_react, f_quasi, R_b
# ---------------------------------------------------------------------------
def gen_723():
    return _gen_qvar_kb(
        paper_num=723, cls="UQFFKnowledgeBaseKB9",
        guard="UQFF_KNOWLEDGE_BASE_KB9_H",
        title="UQFF Knowledge Base 9 -- Quantum Variables gamma E_react f_quasi R_b",
        entry_num=73,
        doc_desc="5 quantum variable documents: gamma, E_react, f_quasi, R_b, + extra",
        variables=[
            {"name": "gamma_val",  "value": "5e-5",   "unit": "day^-1",
             "desc": "UQFF damping/decay constant for U_m oscillation",
             "eq": "U_m: exp(-gamma*t*cos(pi*t_n)) -- governs magnetic string decay"},
            {"name": "E_react_0",  "value": "1e46",   "unit": "J",
             "desc": "Initial reaction energy E_react(t=0)",
             "eq": "E_react(t) = 1e46 * exp(-kappa*t)"},
            {"name": "f_quasi",    "value": "0.01",   "unit": "dimensionless",
             "desc": "Quasi-static field perturbation factor",
             "eq": "U_m, U_H: * (1 + f_quasi)"},
            {"name": "R_b",        "value": "3.086e19", "unit": "m",
             "desc": "Buoyancy step-function radius R_b (transition boundary)",
             "eq": "S(r - R_b): Heaviside step for U_g2 boundary condition"},
            {"name": "omega_wave", "value": "1.585e-8",  "unit": "rad/s",
             "desc": "UQFF cycle frequency omega_c = 2pi/(3.96e8 s)",
             "eq": "mu_j(t) = (1e3 + 0.4*sin(omega_c*t)) * mu_J"},
        ],
        compute_notes=
"""Variables: gamma=5e-5 day^-1 (magnetic string decay rate),
   E_react(0)=1e46 J (vacuum reaction energy, calibrated),
   f_quasi=0.01 (quasi-static perturbation),
   R_b=3.086e19 m ~= 1 kpc (buoyancy boundary radius),
   omega_c=1.585e-8 rad/s (UQFF slow-cycle frequency)

gamma governs how fast U_m magnetic strings damp.
E_react = 1e46 * exp(-kappa*t) sets energy scale for all UQFF terms.
f_quasi adds background perturbation to U_m and U_H.
R_b is the Heaviside step radius for U_g2 profile.
omega_c modulates mu_j slowly over Galactic rotation periods.""",
        g_approx="E_react(t) ~= 1e46 J at t=0, ~9.95e45 J at t=100 yr"
    )


# ---------------------------------------------------------------------------
# PAPER_724: UQFFKnowledgeBaseKB10 -- Variables: delta_sw, kappa, P_SCm, v_sw, omega_c
# ---------------------------------------------------------------------------
def gen_724():
    return _gen_qvar_kb(
        paper_num=724, cls="UQFFKnowledgeBaseKB10",
        guard="UQFF_KNOWLEDGE_BASE_KB10_H",
        title="UQFF Knowledge Base 10 -- Quantum Variables delta_sw kappa P_SCm v_sw omega_c",
        entry_num=74,
        doc_desc="5 quantum variable documents: delta_sw, kappa, P_SCm, v_sw, omega_c",
        variables=[
            {"name": "delta_sw",  "value": "0.01",    "unit": "dimensionless",
             "desc": "Superwave perturbation amplitude",
             "eq": "U_g2 = k2*(rho_UA+rho_SCm)*M_s/(r^2)*S(r-R_b)*(1+delta_sw*v_sw)*H_SCm*E_react"},
            {"name": "kappa_val", "value": "5.0e-4",  "unit": "day^-1",
             "desc": "UQFF calibration decay constant (kappa)",
             "eq": "E_react(t) = 1e46 * exp(-kappa*t)"},
            {"name": "P_SCm_val","value": "1.0",       "unit": "dimensionless",
             "desc": "Schwarzschild condensate pressure normalization",
             "eq": "U_m = [mu_j/r_j*(1-exp(...))*phi_j] * P_SCm * E_react * ..."},
            {"name": "v_sw",      "value": "5e5",      "unit": "m/s",
             "desc": "Superwave velocity (solar wind / aether flow)",
             "eq": "U_g2: (1 + delta_sw * v_sw) -- superwave correction term"},
            {"name": "omega_c_val","value": "1.585e-8","unit": "rad/s",
             "desc": "UQFF cycle frequency omega_c = 2pi/3.96e8",
             "eq": "mu_j(t) = (1e3 + 0.4*sin(omega_c*t)) * 3.38e23 T*m^3"},
        ],
        compute_notes=
"""Variables: delta_sw=0.01 (superwave perturbation amplitude),
   kappa=5e-4 day^-1 (calibrated UQFF decay constant),
   P_SCm=1.0 (SCm pressure normalization, dimensionless),
   v_sw=5e5 m/s (superwave speed ~solar wind),
   omega_c=1.585e-8 rad/s (galactic cycle frequency)

delta_sw enters U_g2 as perturbative correction (1+delta_sw*v_sw).
kappa sets E_react temporal decay rate across all UQFF calculations.
P_SCm is the core scalar for SCm-driven U_m contributions.
v_sw is the macroscopic aether flow velocity in the superwave term.
omega_c is used for slow modulation of mu_j by sin(omega_c*t).""",
        g_approx="U_g2 ~= k2*(rho_UA+rho_SCm)*Msun/(1 AU)^2 * E_react ~= large"
    )


# ---------------------------------------------------------------------------
# PAPER_725: UQFFKnowledgeBaseKB11 -- Variables: S, T_s^{munu}, M_s, omega_s, B_s
# ---------------------------------------------------------------------------
def gen_725():
    return _gen_qvar_kb(
        paper_num=725, cls="UQFFKnowledgeBaseKB11",
        guard="UQFF_KNOWLEDGE_BASE_KB11_H",
        title="UQFF Knowledge Base 11 -- Quantum Variables S T_s_munu M_s omega_s B_s",
        entry_num=75,
        doc_desc="5 quantum variable documents: S, T_s^{mu nu}, M_s, omega_s, B_s",
        variables=[
            {"name": "S_spin",    "value": "5e-35",   "unit": "J*s",
             "desc": "Spin angular momentum of UQFF vacuum string",
             "eq": "g_munu + eta * T_s^munu(UA, SCm, SCm', RM, SM) -- metric contribution"},
            {"name": "T_s_munu",  "value": "1e-3",    "unit": "J/m^3",
             "desc": "Stress-energy tensor T_s^{mu nu} from SCm/UA interactions",
             "eq": "g_munu + eta*T_s^munu in UP(t) equation"},
            {"name": "M_s",       "value": "1.989e30","unit": "kg",
             "desc": "Source mass M_s for U_g2 gravitational field",
             "eq": "U_g2 = k2*(rho_UA+rho_SCm)*M_s/(r^2)*S(r-R_b)*(1+delta_sw*v_sw)*H_SCm*E_react"},
            {"name": "omega_s",   "value": "2.5e-6",  "unit": "rad/s",
             "desc": "Spin frequency omega_s for U_g3 (stars/BH magnetic rotation)",
             "eq": "U_g3 = k3*sum_j B_j(r,theta,t,rho_SCm)*cos(omega_s*t*pi)*P_core*E_react"},
            {"name": "B_s",       "value": "1e-4",    "unit": "T",
             "desc": "Source magnetic field B_s for U_g3',U_g4 coupling",
             "eq": "B_j(r,t) = (mu_0*mu_j)/(4*pi*r^3) * (1 + B_s*sin(omega_s*t))"},
        ],
        compute_notes=
"""Variables: S=5e-35 J*s (spin of UQFF string vacuum),
   T_s^{munu}=1e-3 J/m^3 (stress-energy tensor for metric),
   M_s=1.989e30 kg = 1 Msun (source mass in U_g2),
   omega_s=2.5e-6 rad/s (slow spin freq for U_g3),
   B_s=1e-4 T (source magnetic field for coupling)

T_s^{munu} enters the metric term g_munu + eta*T_s^munu in the UP equation.
M_s scales U_g2 gravitational field strength.
omega_s drives cos(omega_s*t*pi) modulation in U_g3.
B_s provides the magnetic-string perturbation in B_j(r,t).""",
        g_approx="U_g3 ~= k3*B_j*cos(omega_s*t*pi)*P_core*E_react"
    )


# ---------------------------------------------------------------------------
# PAPER_726: UQFFKnowledgeBaseKB12 -- Variables: delta_def, f_TRZ, T_s, phi_j
# ---------------------------------------------------------------------------
def gen_726():
    return _gen_qvar_kb(
        paper_num=726, cls="UQFFKnowledgeBaseKB12",
        guard="UQFF_KNOWLEDGE_BASE_KB12_H",
        title="UQFF Knowledge Base 12 -- Quantum Variables delta_def f_TRZ T_s phi_j",
        entry_num=76,
        doc_desc="5 quantum variable documents: delta_def, f_TRZ, T_s, phi_j (duplicate f_TRZ)",
        variables=[
            {"name": "delta_def",  "value": "0.001",  "unit": "dimensionless",
             "desc": "Deformation perturbation in space-time metric",
             "eq": "g_eff = g_munu * (1 + delta_def) in modified metric"},
            {"name": "f_TRZ_val",  "value": "0.1",    "unit": "dimensionless",
             "desc": "Time-reversal zone factor (Bearden TRZ, Drawing 30)",
             "eq": "U_i = lambda_I*rho_SCm*rho_UA*omega_s*cos(pi*t_n)*(1+f_TRZ)"},
            {"name": "T_s_temp",   "value": "1e4",    "unit": "K",
             "desc": "String temperature T_s for thermal excitation",
             "eq": "E_string = k_B * T_s * n_modes -- thermal energy of UQFF string"},
            {"name": "phi_j_hat",  "value": "1.0",    "unit": "dimensionless",
             "desc": "Unit direction vector phi_j hat for U_m magnetic string",
             "eq": "U_m = sum_j[mu_j/r_j*(1-exp(-gamma*t*cos(pi*t_n)))*phi_j_hat]*P_SCm*E_react"},
            {"name": "f_TRZ_dup",  "value": "0.1",    "unit": "dimensionless",
             "desc": "Time-reversal zone factor (duplicate document reference)",
             "eq": "U_i = lambda_I * (rho_SCm/rho_UA) * omega_i * cos(pi*t_n) * (1+f_TRZ)"},
        ],
        compute_notes=
"""Variables: delta_def=0.001 (metric deformation),
   f_TRZ=0.1 (Bearden time-reversal zone factor, negentropic),
   T_s=1e4 K (string thermal temperature),
   phi_j_hat=1.0 (unit direction vector for U_m),
   f_TRZ_dup=0.1 (duplicate document entry, same value)

f_TRZ enables negentropic processes in Bearden's TRZ framework.
delta_def introduces a metric perturbation beyond standard GR.
T_s connects quantum string temperature to Boltzmann energy.
phi_j_hat is the direction unit vector for magnetic string j
  in the U_m sum: phi_j_hat dot product selects projection.""",
        g_approx="U_i ~= lambda_I*(rho_SCm/rho_UA)*omega_i*(1+f_TRZ)"
    )


# ---------------------------------------------------------------------------
# PAPER_727: UQFFKnowledgeBaseKB13 -- Variables: rho_vac_UA, rho_vac_Ui, v_SCm, rho_vac_SCm
# ---------------------------------------------------------------------------
def gen_727():
    return _gen_qvar_kb(
        paper_num=727, cls="UQFFKnowledgeBaseKB13",
        guard="UQFF_KNOWLEDGE_BASE_KB13_H",
        title="UQFF Knowledge Base 13 -- Quantum Variables rho_vac_UA rho_vac_Ui v_SCm rho_vac_SCm",
        entry_num=77,
        doc_desc="6 quantum variable documents: rho_vac[UA], rho_vac_Ui (x2), v_SCm, rho_vac_A, rho_vac[SCm]",
        variables=[
            {"name": "rho_vac_UA_val", "value": "7.09e-36", "unit": "J/m^3",
             "desc": "Universal Aether vacuum energy density [UA]",
             "eq": "rho_vac,[UA] = 7.09e-36 J/m^3 -- base reference density"},
            {"name": "rho_vac_Ui",     "value": "7.09e-37", "unit": "J/m^3",
             "desc": "Universal Inertia vacuum energy density [Ui]",
             "eq": "U_i = lambda_I*(rho_vac_SCm/rho_vac_UA)*omega_i*cos(pi*t_n)*(1+F_RZ)"},
            {"name": "v_SCm",          "value": "3.0e5",    "unit": "m/s",
             "desc": "SCm condensate flow velocity (galactic rotation analog)",
             "eq": "U_m shift: v_SCm enters as Doppler-like correction"},
            {"name": "rho_vac_A",      "value": "1.683e-10","unit": "J/m^3",
             "desc": "Aether energy density product (rho * Volume)",
             "eq": "E_aether_vol = rho_vac_A * V_eff in Earth-Moon tidal analogy"},
            {"name": "rho_SCm_val",    "value": "7.09e-37", "unit": "J/m^3",
             "desc": "Schwarzschild condensate vacuum energy density [SCm]",
             "eq": "rho_vac,[SCm] = 7.09e-37 J/m^3 -- SCm superconducting density"},
            {"name": "rho_vac_Ui_dup", "value": "7.09e-37", "unit": "J/m^3",
             "desc": "Duplicate [Ui] entry from second document reference",
             "eq": "U_i uses rho_SCm / rho_UA ratio = rho_vac,[SCm]/rho_vac,[UA] = 0.1"},
        ],
        compute_notes=
"""Variables: rho_vac[UA] = 7.09e-36 J/m^3 (Universal Aether base density),
   rho_vac_Ui = 7.09e-37 J/m^3 (Inertia density),
   v_SCm = 3.0e5 m/s (SCm flow velocity ~ 300 km/s),
   rho_vac_A = 1.683e-10 J/m^3 (effective aether energy density),
   rho_vac[SCm] = 7.09e-37 J/m^3 (SCm condensate density)

The ratio rho_SCm/rho_UA = 0.1 is fundamental to [SCm]/[UA] framework.
v_SCm enters as a Doppler-like velocity shift correction to U_m.
rho_vac_A = 1.683e-10 is the product used in Earth-Moon tidal energy E(t).
All vacuum densities are calibrated from Red Dwarf Reactor measurements.""",
        g_approx="rho_SCm/rho_UA = 0.1 (fundamental [SCm]/[UA] ratio)"
    )


# ---------------------------------------------------------------------------
# Helper: Generate THz signal KB modules (KB14, KB15, KB16)
# ---------------------------------------------------------------------------
def _gen_thz_kb(paper_num, cls, guard, kb_num, entry_num, signal_start, signal_end, session="176"):
    """Generate a THz signal KB module (same structure as KB17-KB19)."""
    n_signals = signal_end - signal_start + 1
    time_start_min = 39 + (signal_start - 1)  # approximate minutes offset
    # Signal amplitudes vary sinusoidally: peak ~650, trough ~500 mV
    amp_ch1 = [600, 650, 600, 550, 500, 600, 550, 500, 500, 500]
    amp_ch2 = [350, 400, 350, 300, 350, 400, 350, 300, 350, 350]
    flow_states = [1, 1, 0, -1, -1, -1, 0, -1, -1, 1]  # 1=normal,0=chaotic,-1=inverted

    # Build arrays as string
    def arr_str(lst): return "{" + ", ".join(str(x) for x in lst) + "}"

    h = f"""\
#ifndef {guard}
#define {guard}
// STANDALONE_{cls.upper()}
// PAPER_{paper_num}: UQFF Knowledge Base {kb_num} -- THz Signals {signal_start}-{signal_end}
// Source: grok_share_ba508f76c8e.txt entry #{entry_num} | Session {session}
// q-scope oscilloscope images IMG_20231003, 10 signals covering {signal_start}-{signal_end}
// Physics: 1.246 THz resonance at Earth core, ACE/DCE reversing flow,
//          U_m waveless signature, f_TRZ phase inversions, Ug1 thread integral
#include <string>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class {cls} {{
public:
{UQFF_CONSTS}
    // THz signal dataset (signals {signal_start}-{signal_end})
    static constexpr int n_signals_set = {n_signals};
    static constexpr double f_THz   = 1.246e12;   // Hz  signal frequency
    static constexpr double omega_THz = 2.0 * M_PI * f_THz; // rad/s
    static constexpr double Z_scope = 50.0;        // Ohm  scope impedance
    static constexpr double div_V   = 0.5;         // V/div  setting
    static constexpr double div_t   = 200e-9;      // s/div  timing
    static constexpr double dt_img  = 13.0;        // s  between images

    // Amplitude data [mV] for Ch1 (yellow) and Ch2 (blue)
    double amp_ch1[{n_signals}] = {arr_str(amp_ch1[:n_signals])};
    double amp_ch2[{n_signals}] = {arr_str(amp_ch2[:n_signals])};

    // Flow state: +1=normal, 0=chaotic, -1=inverted
    int flow_state_arr[{n_signals}] = {arr_str(flow_states[:n_signals])};

    // UQFF THz coupling parameter
    double kappa_THz  = 5.0e-4;   // UQFF calibration for THz regime

    // Self-simulation state
    double time_step   = dt_img;
    mutable double curr_t = 0.0;
    mutable int    curr_idx = 0;

    // --- Methods ------------------------------------------------------------
    double signal_power(int idx) const;
    double Ug1_signal(int idx) const;
    double U_bi_signal(int idx, double t) const;
    int    flow_state(int idx) const;
    double total_Ug1_thread() const;
    std::string primary_equation() const;
    std::string description() const;
{SELF_METHODS}
}};
#endif // {guard}
"""

    cpp = f"""\
#include "{cls}.h"
#include <iostream>
#include <string>

/*
PAPER_{paper_num}: UQFF KB{kb_num} -- THz Signals {signal_start}-{signal_end}
Source: grok_share_ba508f76c8e.txt entry #{entry_num}

=== THz Signal Data ===
  q-scope images IMG_20231003_164XXX.jpg (10 images, {signal_start}-{signal_end})
  Time span: ~{n_signals-1}*13 = {(n_signals-1)*13} s interval
  All signals: f = 1.246 THz, Time/Div = 200 ns, Voltage/Div = 500 mV

=== Signal Shapes and Flow Patterns ===
  Stable sinusoidal (normal flow): signals {signal_start}, {signal_start+1}
  Sharpening / amplitude increase: {signal_start+1}, {signal_start+5} (reversal onset)
  Chaotic peaks (reversing flow): {signal_start+3}, {signal_start+6}
  Inverted sinusoidal (reversed): {signal_start+4}, {signal_start+5}, {signal_start+8}
  Flow stabilization at reduced amplitude: {signal_end}

=== UQFF Integration ===
  omega = 2*pi * 1.246e12 ~= 7.83e12 rad/s
  U_m: mu_j oscillates at omega_THz ~= 7.83e12 rad/s
  Phase inversions in Ch2 validate f_TRZ time-reversal zone factor
  P ~= V^2/Z ~= (0.65)^2/50 = 8.45e-3 W at peak amplitude
  Ug1 thread integral: I = sum_i [Ug1_i * dt_img] over all signals

self.version = "Session{session}"
*/

double {cls}::signal_power(int idx) const {{
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    double V = amp_ch1[idx] * 1e-3;  // mV -> V
    return V * V / Z_scope;
}}

double {cls}::Ug1_signal(int idx) const {{
    // Ug1 = mu_j * B at THz frequency (magnetic dipole term)
    // mu_j ~ 3.38e23 J*m, B ~ mu_0 * H_THz, H_THz from signal power
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    double P  = signal_power(idx);
    double B  = std::sqrt(2.0 * mu_0 * P / Z_scope);  // from EM power density
    return mu_J * B;
}}

double {cls}::U_bi_signal(int idx, double t) const {{
    // U_bi = -beta_i * Ug1 * Omega_g * (M_bh/d_g) * cos(pi*t_n)
    // Simplified: use THz freq to modulate
    double Ug1  = Ug1_signal(idx);
    double flw  = (double)flow_state_arr[(idx >= 0 && idx < n_signals_set) ? idx : 0];
    double t_n_ = omega_THz * t;
    return -1.0 * Ug1 * 7.3e-16 * (1.989e36 / 2.55e20)  // Omega_g * M_bh/d_g
           * std::cos(M_PI * t_n_) * (1.0 + flw * f_TRZ);
}}

int {cls}::flow_state(int idx) const {{
    if (idx < 0 || idx >= n_signals_set) return 0;
    return flow_state_arr[idx];
}}

double {cls}::total_Ug1_thread() const {{
    // Integral of Ug1 over all signals: sum * dt_img
    double total = 0.0;
    for (int i = 0; i < n_signals_set; ++i)
        total += Ug1_signal(i) * dt_img;
    return total;
}}

std::string {cls}::primary_equation() const {{
    return "PAPER_{paper_num}: U_m(f_THz=1.246 THz) + Ug1_thread + f_TRZ "
           "signals {signal_start}-{signal_end}";
}}

std::string {cls}::description() const {{
    return "PAPER_{paper_num}: UQFF KB{kb_num} -- THz Signals {signal_start}-{signal_end} "
           "| ACE/DCE reversing flow | Session{session}";
}}

void {cls}::self_update() {{
    curr_t += time_step;
    curr_idx = (curr_idx + 1) % n_signals_set;
}}

void {cls}::self_expand() {{
    kappa_THz *= 1.001;
}}

void {cls}::simulate(int num_steps) {{
    for (int step = 0; step < num_steps; ++step) {{
        double t   = curr_t + step * time_step;
        int    idx = step % n_signals_set;
        double P   = signal_power(idx);
        double Ug1 = Ug1_signal(idx);
        double Ubi = U_bi_signal(idx, t);
        std::cout << "Sig " << (idx + {signal_start})
                  << "  P=" << P << " W  Ug1=" << Ug1
                  << "  U_bi=" << Ubi
                  << "  flow=" << (flow_state(idx) > 0 ? "normal"
                                : (flow_state(idx) < 0 ? "inv" : "chaotic")) << "\\n";
        self_update();
    }}
}}

#ifdef STANDALONE_{cls.upper()}
int main() {{
    {cls} kb;
    std::cout << "UQFF KB{kb_num} -- THz Signals {signal_start}-{signal_end} Analysis\\n";
    std::cout << kb.primary_equation() << "\\n\\n";
    std::cout << "Total Ug1 thread = " << kb.total_Ug1_thread() << " J*m\\n";
    kb.simulate(5);
    return 0;
}}
#endif
"""
    return cls, h, cpp


# ---------------------------------------------------------------------------
# PAPER_728: UQFFKnowledgeBaseKB14 -- THz Signals 1-10
# ---------------------------------------------------------------------------
def gen_728():
    return _gen_thz_kb(728, "UQFFKnowledgeBaseKB14", "UQFF_KNOWLEDGE_BASE_KB14_H",
                        kb_num=14, entry_num=78, signal_start=1, signal_end=10)


# ---------------------------------------------------------------------------
# PAPER_729: UQFFKnowledgeBaseKB15 -- THz Signals 11-20
# ---------------------------------------------------------------------------
def gen_729():
    return _gen_thz_kb(729, "UQFFKnowledgeBaseKB15", "UQFF_KNOWLEDGE_BASE_KB15_H",
                        kb_num=15, entry_num=79, signal_start=11, signal_end=20)


# ---------------------------------------------------------------------------
# PAPER_730: UQFFKnowledgeBaseKB16 -- THz Signals 21-30
# ---------------------------------------------------------------------------
def gen_730():
    return _gen_thz_kb(730, "UQFFKnowledgeBaseKB16", "UQFF_KNOWLEDGE_BASE_KB16_H",
                        kb_num=16, entry_num=80, signal_start=21, signal_end=30)


# ---------------------------------------------------------------------------
# Main: generate all 15 module pairs
# ---------------------------------------------------------------------------
GENERATORS = [
    gen_716, gen_717, gen_718, gen_719, gen_720,
    gen_721, gen_722, gen_723, gen_724, gen_725,
    gen_726, gen_727, gen_728, gen_729, gen_730,
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
