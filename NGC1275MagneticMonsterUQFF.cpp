// STANDALONE_NGC1275MAGNETICMONSTERUQFF
#include "NGC1275MagneticMonsterUQFF.h"
#include <iostream>

// Base gravitational acceleration
double NGC1275MagneticMonsterUQFF::g_grav(double r) const {
    double M_tot = M_SMBH + M_stellar;
    return G * M_tot / (r * r);
}

// Hubble parameter at z=0.0176
// H(z) = H_0*sqrt(0.3*(1+z)^3 + 0.7)
double NGC1275MagneticMonsterUQFF::H_z() const {
    double z = z_ngc;
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + z, 3) + 0.7);
}

// BH feedback reduction: (1 - F_BH(t))
double NGC1275MagneticMonsterUQFF::one_minus_FBH(double t) const {
    double F = F_0 * (1.0 - std::exp(-t / tau_BH));
    return 1.0 - F;
}

// Filament magnetic support acceleration (upward, stabilising)
// Energy density: B^2/(2*mu_0) * V_fil = 5.649e39 N  => / M_fil => 2.840e3 m/s^2
// Scaled by 1e-12 for macroscopic system
double NGC1275MagneticMonsterUQFF::a_fil() const {
    double energy = (B_fil * B_fil) / (2.0 * mu_0) * V_fil;
    return energy / M_fil * 1.0e-12;
}

// EM aether term
// a_EM = q_p * v_HVS * B_fil / m_p * (1 + rho_UA/rho_SCm) * 1e-12
// ~ 3.16e-5 m/s^2
double NGC1275MagneticMonsterUQFF::a_EM() const {
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_HVS * B_fil / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}

// Full MUGE for NGC 1275
// g_NGC1275(r,t) = g_grav*(1+H*t)*(1-F_BH)*(1+f_TRZ) + a_fil + a_EM
// At t=50 Myr, r=r_gal: ~ 3.16e-5 m/s^2
double NGC1275MagneticMonsterUQFF::g_NGC1275(double r, double t) const {
    double Ht = H_z() * t;
    double g_base = g_grav(r) * (1.0 + Ht) * one_minus_FBH(t) * (1.0 + f_TRZ);
    return g_base + a_fil() + a_EM();
}

std::string NGC1275MagneticMonsterUQFF::primary_equation() const {
    return "g_NGC1275(r,t)=G*(M_SMBH+M_star)/r^2*(1+H_z*t)*(1-F_BH(t))*(1+f_TRZ)+a_fil+q*(vxB)*(1+rho_UA/rho_SCm)*1e-12";
}

void NGC1275MagneticMonsterUQFF::self_update() {
    curr_t += time_step;
    M_stellar *= (1.0 - 1.0e-5 * time_step / 3.156e7); // slow mass loss
}

void NGC1275MagneticMonsterUQFF::self_expand() {
    r_gal *= 1.001;
    B_fil *= 0.999;
}

void NGC1275MagneticMonsterUQFF::simulate(int num_steps) {
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        double g = g_NGC1275(r_gal, t);
        std::cout << "Step " << step << "  t=" << t / 3.156e7
                  << " yr  g=" << g << " m/s^2\n";
    }
}

#ifdef STANDALONE_NGC1275MAGNETICMONSTERUQFF
int main() {
    NGC1275MagneticMonsterUQFF ngc;
    std::cout << "NGC 1275 Magnetic Monster UQFF\n";
    std::cout << ngc.primary_equation() << "\n\n";
    double g = ngc.g_NGC1275(ngc.r_gal, 1.578e15);
    std::cout << "g(50 Myr) = " << g << " m/s^2\n";
    ngc.simulate(3);
    return 0;
}
#endif
