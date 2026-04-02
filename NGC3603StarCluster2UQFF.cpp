// STANDALONE_NGC3603STARCLUSTER2UQFF
#include "NGC3603StarCluster2UQFF.h"
#include <iostream>

// M_dot: Bok globule secondary star formation mass fraction
// M_dot(t) = M_0 * exp(-t/tau_SF)
double NGC3603StarCluster2UQFF::M_dot(double t) const {
    return M_0 * std::exp(-t / tau_SF);
}

// M_t: M_initial * (1 + M_dot(t))
double NGC3603StarCluster2UQFF::M_t(double t) const {
    return M_initial * (1.0 + M_dot(t));
}

// Stellar wind cavity pressure reduction
// P(t) = P_0 * exp(-t/tau_exp) => (1 - P) at t=0.5 Myr ~ 0.9394
double NGC3603StarCluster2UQFF::one_minus_P(double t) const {
    return 1.0 - P_0 * std::exp(-t / tau_exp);
}

// EM aether term ~ 1.053e-3 m/s^2
double NGC3603StarCluster2UQFF::a_EM() const {
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_clust / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}

// Full MUGE: g = G*M(t)/r^2 * (1+H0*t) * (1-P(t)) * (1+f_TRZ) + a_EM
// At t=0.5 Myr: ~ 1.053e-3 m/s^2
double NGC3603StarCluster2UQFF::g_NGC3603v2(double t) const {
    double M   = M_t(t);
    double g_n = G * M / (r_clust * r_clust);
    return g_n * (1.0 + H_0 * t) * one_minus_P(t) * (1.0 + f_TRZ) + a_EM();
}

std::string NGC3603StarCluster2UQFF::primary_equation() const {
    return "g_NGC3603v2(t)=G*M(t)/r^2*(1+H0*t)*(1-P(t))*(1+f_TRZ)+q_p*v_gas*B/m_p*(1+rho_UA/rho_SCm)*1e-12";
}

void NGC3603StarCluster2UQFF::self_update() {
    curr_t += time_step;
}

void NGC3603StarCluster2UQFF::self_expand() {
    r_clust  *= 1.001;
    M_initial *= 1.0001;
}

void NGC3603StarCluster2UQFF::simulate(int num_steps) {
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        double g = g_NGC3603v2(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\n";
    }
}

#ifdef STANDALONE_NGC3603STARCLUSTER2UQFF
int main() {
    NGC3603StarCluster2UQFF n3;
    std::cout << "NGC 3603 Star Cluster v2 UQFF\n";
    std::cout << n3.primary_equation() << "\n\n";
    std::cout << "g(0.5 Myr) = " << n3.g_NGC3603v2(1.578e13) << " m/s^2\n";
    n3.simulate(3);
    return 0;
}
#endif
