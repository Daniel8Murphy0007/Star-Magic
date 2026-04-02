// STANDALONE_NGC3603STARCLUSTERPRIMARYUQFF
#include "NGC3603StarClusterPrimaryUQFF.h"
#include <iostream>

double NGC3603StarClusterPrimaryUQFF::M_eff(double t) const {
    return M_base * (1.0 + M_gas_frac * std::exp(-t / tau_SF));
}

double NGC3603StarClusterPrimaryUQFF::g_grav(double t) const {
    return G * M_eff(t) / (r_clust * r_clust);
}

double NGC3603StarClusterPrimaryUQFF::H_z() const {
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + z_ngc, 3) + 0.7);
}

double NGC3603StarClusterPrimaryUQFF::a_EM() const {
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_clust / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}

// g_NGC3603primary = g_grav*(1+H*t)*(1-alpha_wind)*(1+f_TRZ) + a_EM
double NGC3603StarClusterPrimaryUQFF::g_NGC3603primary(double t) const {
    double g_base = g_grav(t) * (1.0 + H_z() * t)
                              * (1.0 - alpha_wind)
                              * (1.0 + f_TRZ);
    return g_base + a_EM();
}

std::string NGC3603StarClusterPrimaryUQFF::primary_equation() const {
    return "g_NGC3603p(t)=G*M_eff(t)/r^2*(1+H_z*t)*(1-alpha_wind)*(1+f_TRZ)+q*v*B/m_p*(1+rho_UA/rho_SCm)*1e-12";
}

void NGC3603StarClusterPrimaryUQFF::self_update() {
    curr_t += time_step;
}

void NGC3603StarClusterPrimaryUQFF::self_expand() {
    r_clust   *= 1.001;
    M_gas_frac *= 0.999; // gas fraction consumed over time
}

void NGC3603StarClusterPrimaryUQFF::simulate(int num_steps) {
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        double g = g_NGC3603primary(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\n";
    }
}

#ifdef STANDALONE_NGC3603STARCLUSTERPRIMARYUQFF
int main() {
    NGC3603StarClusterPrimaryUQFF n3p;
    std::cout << "NGC 3603 Primary UQFF\n";
    std::cout << n3p.primary_equation() << "\n\n";
    std::cout << "g(1 Myr) = " << n3p.g_NGC3603primary(3.156e13) << " m/s^2\n";
    n3p.simulate(3);
    return 0;
}
#endif
