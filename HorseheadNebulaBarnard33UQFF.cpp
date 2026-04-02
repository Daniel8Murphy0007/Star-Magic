// STANDALONE_HORSEHEADNEBULABARNARD33UQFF
#include "HorseheadNebulaBarnard33UQFF.h"
#include <iostream>

double HorseheadNebulaBarnard33UQFF::g_grav() const {
    return G * M_neb / (r_neb * r_neb);
}

// Hubble expansion at z~3e-4
double HorseheadNebulaBarnard33UQFF::H_z() const {
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + z_neb, 3) + 0.7);
}

// Erosion factor: as stars erode the nebula
// E(t) = E_0*(1 - exp(-t/tau)) => 1-E = 1 - 0.03626 at t=1 Myr
double HorseheadNebulaBarnard33UQFF::one_minus_E(double t) const {
    double E = E_0 * (1.0 - std::exp(-t / tau_erode));
    return 1.0 - E;
}

// Radiation pressure acceleration
// P_rad = (L/(4*pi*r^2*c)) * (rho/m_H) ~ 4.347e-5 m/s^2
double HorseheadNebulaBarnard33UQFF::a_rad() const {
    const double m_H = 1.67e-27;
    double flux = L_star / (4.0 * M_PI * r_neb * r_neb * c);
    return flux * (rho_neb / m_H);
}

// EM aether term ~ 1.053e-3 m/s^2
double HorseheadNebulaBarnard33UQFF::a_EM() const {
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_neb / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}

// Full MUGE
// g_HH = g_grav*(1+H_z*t)*(1-E(t))*(1+f_TRZ) + a_rad + a_EM ~ 1.097e-3 m/s^2
double HorseheadNebulaBarnard33UQFF::g_Horsehead(double t) const {
    double g_base = g_grav() * (1.0 + H_z() * t) * one_minus_E(t) * (1.0 + f_TRZ);
    return g_base + a_rad() + a_EM();
}

std::string HorseheadNebulaBarnard33UQFF::primary_equation() const {
    return "g_HH(t)=G*M/r^2*(1+H_z*t)*(1-E(t))*(1+f_TRZ)+P_rad+(q*v_gas*B)*(1+rho_UA/rho_SCm)*1e-12";
}

void HorseheadNebulaBarnard33UQFF::self_update() {
    curr_t += time_step;
    M_neb *= (1.0 - E_0 * 0.01); // gradual mass erosion
}

void HorseheadNebulaBarnard33UQFF::self_expand() {
    r_neb *= 1.001;
    rho_neb *= 0.999;
}

void HorseheadNebulaBarnard33UQFF::simulate(int num_steps) {
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        double g = g_Horsehead(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\n";
    }
}

#ifdef STANDALONE_HORSEHEADNEBULABARNARD33UQFF
int main() {
    HorseheadNebulaBarnard33UQFF hh;
    std::cout << "Horsehead Nebula Barnard 33 UQFF\n";
    std::cout << hh.primary_equation() << "\n\n";
    std::cout << "g(1 Myr) = " << hh.g_Horsehead(3.156e13) << " m/s^2\n";
    hh.simulate(3);
    return 0;
}
#endif
