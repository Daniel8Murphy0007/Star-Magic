// STANDALONE_NGC1316MUGECALCULATION
#include "NGC1316MUGECalculation.h"
#include <iostream>
#include <fstream>

// M(t): Total mass with merger decay
// M(t) = (M_visible + M_DM) + M_spiral * exp(-t/tau_merge)
double NGC1316MUGECalculation::M_t(double t) const {
    return (M_visible + M_DM) + M_spiral * std::exp(-t / tau_merge);
}

// r(t): Dynamical radius (expanding with 1 km/s velocity)
double NGC1316MUGECalculation::r_t(double t) const {
    return r_0 + 1.0e3 * t;
}

// H(t,z): Hubble parameter at redshift z
// H(z) = H_0 * sqrt(0.3*(1+z)^3 + 0.7)  [flat LCDM]
double NGC1316MUGECalculation::H_tz(double zz) const {
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + zz, 3) + 0.7);
}

// F_env(t): Environmental force (tidal + cluster disruption)
// F_tidal = G * M_spiral / d^2  [merger remnant tidal forcing]
// F_cluster = k_cluster * M_cluster  [globular cluster tidal stripping]
double NGC1316MUGECalculation::F_env(double t) const {
    double F_tidal   = G * M_spiral / (d_spiral * d_spiral);
    double F_cluster = k_cluster * M_cluster;
    return F_tidal + F_cluster;
}

// U_g1: AGN magnetic dipole (Blandford-Znajek mechanism)
// mu_dipole = I * A * omega_spin  => approximated by mu_J * omega_spin
// U_g1 = mu_dipole * B_AGN
double NGC1316MUGECalculation::U_g1(double t) const {
    double mu_dipole = mu_J * omega_spin;
    return mu_dipole * B_AGN;
}

// U_g2: AGN superconductor (aether jet field)
// B_super = mu_0 * H_aether  (H_aether ~ 1e-5 A/m)
// U_g2 = B_super^2 / (2 * mu_0)
double NGC1316MUGECalculation::U_g2() const {
    double H_aether = 1.0e-5;
    double B_super  = mu_0 * H_aether;
    return (B_super * B_super) / (2.0 * mu_0);
}

// U_g3': Merger remnant gravitational influence
// U_g3' = G * M_spiral / d_spiral^2
double NGC1316MUGECalculation::U_g3_prime() const {
    return G * M_spiral / (d_spiral * d_spiral);
}

// U_g4: Reactive vacuum energy (exponentially decaying)
// U_g4 = k_4 * E_react(t),  E_react = 1e46 * exp(-0.0005 * t)
double NGC1316MUGECalculation::U_g4(double t) const {
    double E_react = 1.0e46 * std::exp(-0.0005 * t);
    return 1.0 * E_react;
}

// U_i: UQFF vacuum oscillation integral
// U_i = lambda_I * (rho_SCm/rho_UA) * omega_i * cos(pi*t_n) * (1 + f_TRZ)
double NGC1316MUGECalculation::U_i(double t) const {
    double t_n = 0.0;
    return 1.0 * (rho_SCm / rho_UA) * omega_i * std::cos(M_PI * t_n) * (1.0 + f_TRZ);
}

// psi_integral: Dust lane wavefunction contribution
// psi_dust = A * exp(-r^2/(2*sigma^2)) * exp(i(m*theta - omega*t))
// Real part: A * exp(-r^2/(2*sigma^2)) * cos(omega*t)
// Integral approximated: magnitude * volume element
double NGC1316MUGECalculation::psi_integral(double r, double t) const {
    double A = 1.0e-10;
    double psi_mag = A * std::exp(-r * r / (2.0 * sigma_dust * sigma_dust));
    double phase   = std::cos(omega_i * t);
    return psi_mag * phase * V_gal;
}

// g_NGC1316: Full Master Universal Gravity Equation for NGC 1316
// g(r,t) = [G*M(t)/r(t)^2] * [1+H(z)] * [1-B/B_crit] * [1+F_env]
//         + (U_g1 + U_g2 + U_g3' + U_g4) + U_i + Lambda*c^2/3
//         + hbar/sqrt(dx*dp) * psi_integral * (2*pi/t_H)
//         + rho_dust*V*g0 + (M_vis+M_DM)*(delta_rho/rho + 3*G*M/r^3)
double NGC1316MUGECalculation::g_NGC1316(double r, double t) const {
    double M   = M_t(t);
    double rt  = r_t(t);
    double H   = H_tz(z);
    double F   = F_env(t);
    double B_crit = 1.0e11;

    // Core Newtonian + UQFF modifiers
    double g_core = (G * M / (rt * rt)) * (1.0 + H) * (1.0 - B_AGN / B_crit) * (1.0 + F);

    // UQFF potential terms
    double U_sum = U_g1(t) + U_g2() + U_g3_prime() + U_g4(t);

    // Oscillation integral
    double U_osc = U_i(t);

    // Cosmological term
    double g_lambda = Lambda * c * c / 3.0;

    // Wavefunction Hamiltonian term (approximate)
    double dx = 1.0e-10, dp = hbar / (2.0 * dx);
    double g_psi = (hbar / std::sqrt(dx * dp)) * psi_integral(r, t) * (2.0 * M_PI / t_H);

    // Dust fluid dynamics
    double g0        = G * M / (r * r);
    double g_dust    = rho_dust * V_gal * g0;

    // Dark matter density perturbation
    double delta_rho = 1.0e-5;
    double g_dm      = (M_visible + M_DM) * (delta_rho + 3.0 * G * M / (r * r * r));

    return g_core + U_sum + U_osc + g_lambda + g_psi + g_dust + g_dm;
}

std::string NGC1316MUGECalculation::primary_equation() const {
    return "g_NGC1316(r,t) = [G*M(t)/r^2]*[1+H(z)]*[1-B/B_crit]*[1+F_env] + U_g1+U_g2+U_g3p+U_g4 + U_i + Lambda*c^2/3 + psi_H_term + rho_dust*V*g + (M_vis+M_DM)*(delta_rho/rho+3GM/r^3)";
}

void NGC1316MUGECalculation::self_update() {
    curr_t += time_step;
    // NGC 1316 evolves: radius expands dynamically
    r_0 += 1.0e3 * time_step;
}

void NGC1316MUGECalculation::self_expand() {
    // Increase simulation resolution
    sigma_dust *= 1.01;
    M_cluster  *= 1.001;
}

void NGC1316MUGECalculation::simulate(int num_steps) {
    std::vector<double> r_pts;
    for (int i = 1; i <= 100; i++) r_pts.push_back(i * 1.0e3 * kpc);
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        double g_max = 0.0, r_max = r_pts[0];
        for (double r : r_pts) {
            double g = g_NGC1316(r, t);
            if (g > g_max) { g_max = g; r_max = r; }
        }
        std::cout << "Step " << step << "  t=" << t/3.156e7 << " yr"
                  << "  g_peak=" << g_max << " m/s^2 at r=" << r_max/kpc << " kpc\n";
    }
}

#ifdef STANDALONE_NGC1316MUGECALCULATION
int main() {
    NGC1316MUGECalculation ngc;
    std::cout << "NGC 1316 MUGE Simulation\n";
    std::cout << ngc.primary_equation() << "\n\n";
    // Test at 10 kpc, t=100 Myr
    double r_test = 10.0e3 * kpc;
    double t_test = 3.156e15;
    double g = ngc.g_NGC1316(r_test, t_test);
    std::cout << "g(10 kpc, 100 Myr) = " << g << " m/s^2\n";
    ngc.simulate(5);
    return 0;
}
#endif
