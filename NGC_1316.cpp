#include "NGC_1316.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <complex>
#include <cmath>

/*
PAPER_731: NGC 1316 Merger Evolution MUGE
Source: grok_share_ba508f76c8e.txt entry #64 | Session 177

=== NGC 1316 (Fornax A) / "Hubble Spies Cosmic Dust Bunnies" ===
ACS Observations (March 2003): F435W/F555W/F814W colour composite
  - Complex dust lanes, patches, red star clusters detected
  - Fewer low-mass globular clusters (<1e5 M_sun) in inner regions
    consistent with tidal disruption during spiral-galaxy mergers
  - Tidal features (loops, plumes) from mergers 1-3 Gyr ago
  - SMBH ~1e8 M_sun powering radio lobes (Fornax A)

=== Master Universal Gravity Equation (MUGE) for NGC 1316 ===
g_NGC1316(r,t) =
  (G*M(t)) / (r(t)^2) * (1 + H(t,z)) * (1 - B_AGN/B_crit) * (1 + F_env(t))
  + (U_g1 + U_g2 + U_g3' + U_g4) + U_i
  + (Lambda*c^2/3)
  + (hbar/sqrt(dx*dp)) * |psi_dust(r,t)|^2 * (2*pi/t_H)
  + rho_dust * V * g_local
  + (M_visible + M_DM) * (delta_rho/rho + 3*G*M/r^3)

=== Parameter Table ===
  M(t)       = M_visible + M_DM + M_spiral*exp(-t/tau_merge)
             = 3.5e11*M_sun + 1.5e11*M_sun + 1e10*M_sun * exp(-t/3.156e16)
             solved at t=0: M = 5.01e11 M_sun

  r(t)       = r_0 + v_r*t  = 46 kpc + 1e3*t
             solved at t=0: r = 46 kpc = 1.420e24 m

  H(z=0.005) = H_0 * sqrt(0.3*(1.005)^3 + 0.7)  ~= 2.274e-18 s^-1

  F_env(t)   = G*M_spiral/d^2 + k_cluster*M_cluster
             = 6.6743e-11 * 1.989e40 / (1.543e24)^2 + 1e-12 * 1.989e36
             ~= 5.56e-19 + 1.99e24  (dominated by cluster term)

  U_g1       = I*A*omega_spin * B_AGN
             with I=1e22 A, A=1e16 m^2, omega_spin=1e-3 rad/s, B_AGN=1e-4 T
             = 1e22 * 1e16 * 1e-3 * 1e-4 = 1e31  (A*m^2*T = J*rad)

  U_g2       = (mu_0*H_aether)^2 / (2*mu_0)
             = mu_0 * H_aether^2 / 2
             = 1.2566e-6 * (1e-5)^2 / 2 ~= 6.28e-17 J/m^3

  U_g3'      = G*M_spiral / d_spiral^2
             = 6.6743e-11 * 1.989e40 / (1.543e24)^2  ~= 5.56e-19

  U_g4       = k_4 * E_react(t) = 1.0 * 1e46 * exp(-0.0005*t)
             solved at t=0: U_g4 = 1.0e46 J/m^3

  U_i        = lambda_I * (rho_SCm/rho_UA) * omega_i * cos(pi*t_n) * (1+F_RZ)
             = 1.0 * (7.09e-37/7.09e-36) * 1e-8 * cos(0) * 1.01
             = 0.1 * 1e-8 * 1.01 ~= 1.01e-9

  psi_dust   = A*exp(-r^2/(2*sigma^2)) * exp(-i*omega_psi*t)
             |psi_dust|^2 = A^2 * exp(-r^2/sigma^2)
             at r=46 kpc, sigma=2 kpc: ~= A^2 * exp(-529) ~= negligible

  Lambda*c^2/3 = 1.1e-52 * 9e16 / 3 ~= 3.3e-36

  rho_dust*V*g_local:
    g_local = G*M/r^2 ~= 6.6743e-11*5e11*1.989e30 / (1.42e24)^2
            ~= 3.30e-11 m/s^2
    rho_dust*V*g_local = 1e-21 * 1e51 * 3.30e-11 ~= 3.30e19

  (M_vis+M_DM)*(delta_rho/rho + 3*G*M/r^3):
    = 5e11*1.989e30 * (1e-5 + 3*6.6743e-11*5e11*1.989e30/(1.42e24)^3)
    ~= dominant perturbation term

self.version = "Session177"
*/

// ---------------------------------------------------------------------------
double NGC_1316::M_t(double t) const {
    // M(t) = M_visible + M_DM + M_spiral * exp(-t/tau_merge)
    return M_visible + M_DM + M_spiral * std::exp(-t / tau_merge);
}

double NGC_1316::r_t(double t) const {
    // r(t) = r_0 + v_r * t
    return r_0 + v_r * t;
}

double NGC_1316::H_tz(double z) const {
    // H(z) = H_0 * sqrt(Omega_m*(1+z)^3 + Omega_Lambda)
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + z, 3.0) + 0.7);
}

double NGC_1316::F_env(double t) const {
    // Tidal force from merger progenitor
    double F_tidal   = G * M_spiral / (d_spiral * d_spiral);
    // Star cluster disruption
    double F_cluster = k_cluster * M_cluster;
    return F_tidal + F_cluster;
}

double NGC_1316::U_g1(double t) const {
    // U_g1 = mu_dipole * B_AGN
    // mu_dipole = I * A * omega_spin  (AGN jet dipole moment)
    double I_agn     = 1.0e22;   // A  current
    double A_area    = 1.0e16;   // m^2 loop area
    double mu_dipole = I_agn * A_area * omega_spin;
    return mu_dipole * B_AGN;
}

double NGC_1316::U_g2(double t) const {
    // U_g2 = B_super^2 / (2*mu_0),  B_super = mu_0*H_aether
    double B_super = mu_0 * H_aether;
    return (B_super * B_super) / (2.0 * mu_0);
}

double NGC_1316::U_g3_prime(double t) const {
    // U_g3' = G*M_spiral / d_spiral^2  (merger remnant gravitational influence)
    return G * M_spiral / (d_spiral * d_spiral);
}

double NGC_1316::U_g4(double t) const {
    // U_g4 = k_4 * E_react(t),  E_react(t) = 1e46 * exp(-kappa*t)
    double E_react = 1.0e46 * std::exp(-kappa * t);
    return k_4 * E_react;
}

double NGC_1316::U_i(double t) const {
    // U_i = lambda_I * (rho_SCm/rho_UA) * omega_i * cos(pi*t_n) * (1+F_RZ)
    return lambda_I * (rho_SCm / rho_UA) * omega_i
           * std::cos(M_PI * t_n) * (1.0 + F_RZ);
}

double NGC_1316::psi_integral(double r, double t) const {
    // |psi_dust(r,t)|^2 = A^2 * exp(-r^2/sigma^2)
    // quantum wavefunction for dust lanes (Gaussian envelope)
    double A_psi = 1.0e-10;  // normalised amplitude
    double psi2  = A_psi * A_psi * std::exp(-r * r / (sigma_dust * sigma_dust));
    // quantum energy term weight
    double dx = 1.0e-10;     // position uncertainty (m)
    double dp = 1.0e-20;     // momentum uncertainty (kg*m/s)
    double prefac = hbar / std::sqrt(dx * dp);
    return prefac * psi2 * (2.0 * M_PI / t_H);
}

double NGC_1316::g_NGC1316(double r, double t) const {
    // Full MUGE for NGC 1316
    double M      = M_t(t);
    double r_eff  = r_t(t);
    double H      = H_tz(z_redshift);
    double F_e    = F_env(t);
    double g_loc  = G * M / (r_eff * r_eff);  // local Newtonian reference

    // Core Newtonian + UQFF modulation
    double core   = (G * M / (r * r))
                    * (1.0 + H)
                    * (1.0 - B_AGN / B_crit)
                    * (1.0 + F_e);

    // UQFF potential terms
    double uqff   = U_g1(t) + U_g2(t) + U_g3_prime(t) + U_g4(t) + U_i(t);

    // Cosmological constant contribution
    double cosm   = (Lambda * c * c) / 3.0;

    // Quantum dust lane integral
    double quant  = psi_integral(r, t);

    // Fluid / AGN dust drag
    double fluid  = rho_dust * Vol_galaxy * g_loc;

    // Dark matter + visible density perturbation
    double dm     = (M_visible + M_DM)
                    * (delta_rho_rho + (3.0 * G * M) / (r * r * r));

    return core + uqff + cosm + quant + fluid + dm;
}

std::string NGC_1316::primary_equation_str(double r, double t) const {
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(6);
    double val = g_NGC1316(r, t);
    ss << "g_NGC1316(r=" << r/kpc << " kpc, t=" << t << " s) = " << val << " m/s^2\n";
    ss << "  M(t)   = " << M_t(t)      << " kg\n";
    ss << "  r(t)   = " << r_t(t)/kpc  << " kpc\n";
    ss << "  H(z)   = " << H_tz(z_redshift) << " s^-1\n";
    ss << "  F_env  = " << F_env(t)    << "\n";
    ss << "  U_g1   = " << U_g1(t)     << " J\n";
    ss << "  U_g2   = " << U_g2(t)     << " J/m^3\n";
    ss << "  U_g3'  = " << U_g3_prime(t) << "\n";
    ss << "  U_g4   = " << U_g4(t)     << " J/m^3\n";
    ss << "  U_i    = " << U_i(t)      << " J/m^6\n";
    ss << "  Lambda*c^2/3 = " << (Lambda*c*c/3.0) << " m^-2\n";
    ss << "  psi quantum  = " << psi_integral(r, t) << "\n";
    return ss.str();
}

std::string NGC_1316::description() const {
    return "PAPER_731: NGC 1316 Merger Evolution MUGE | "
           "Fornax A elliptical galaxy | merger dynamics / AGN / dust lanes / DM | Session177";
}

void NGC_1316::self_update() {
    curr_t += time_step;
}

void NGC_1316::self_expand() {
    rho_dust   *= 1.0001;   // dust lane density grows with time
    sigma_dust *= 1.0005;   // dust distribution widens
}

void NGC_1316::simulate(int num_steps) {
    std::cout << "=== NGC 1316 Self-Simulation ===\n";
    double r_sim = r_0;   // start at galaxy edge
    for (int step = 0; step < num_steps; ++step) {
        double t = curr_t + step * time_step;
        double g = g_NGC1316(r_sim, t);
        std::cout << "Step " << std::setw(3) << step
                  << "  t=" << std::scientific << t
                  << "  g=" << g
                  << "  M=" << M_t(t) << " kg\n";
        self_update();
    }
}

// ---------------------------------------------------------------------------
#ifdef STANDALONE_NGC1316
int main() {
    NGC_1316 ngc;
    std::cout << "======================================================\n";
    std::cout << ngc.description() << "\n";
    std::cout << "======================================================\n";

    // Point evaluation at galaxy edge, t=0
    double r_test = ngc.r_0;
    double t_test = 0.0;
    std::cout << ngc.primary_equation_str(r_test, t_test) << "\n";

    // Radial profile at t = 100 Myr
    double t_100Myr = 3.156e15;  // s
    std::cout << "--- Radial profile at t = 100 Myr ---\n";
    for (int i = 1; i <= 10; ++i) {
        double r = i * 10.0e3 * NGC_1316::kpc;  // 10-100 kpc steps
        double g = ngc.g_NGC1316(r, t_100Myr);
        std::cout << "  r=" << std::setw(6) << (i*10) << " kpc  g=" 
                  << std::scientific << g << " m/s^2\n";
    }

    // Self-simulation: 5 steps
    std::cout << "\n";
    ngc.simulate(5);

    return 0;
}
#endif
