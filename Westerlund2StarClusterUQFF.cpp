// STANDALONE_WESTERLUND2STARCLUSTERUQFF
#include "Westerlund2StarClusterUQFF.h"
#include <iostream>

// M_dot: gas reservoir rapidly converts to stars
// M_dot(t) = M_0_gas * exp(-t/tau_SF)
// At t=1 Myr: M_dot = 3.333*0.6065 = 2.021 => M(t) = M_init * 3.021
double Westerlund2StarClusterUQFF::M_dot(double t) const {
    return M_0_gas * std::exp(-t / tau_SF);
}

double Westerlund2StarClusterUQFF::M_t(double t) const {
    return M_init * (1.0 + M_dot(t));
}

double Westerlund2StarClusterUQFF::a_EM() const {
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_wd2 / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}

// Full Westerlund 2 MUGE
// g_W2 ~ 1.053e-3 m/s^2 at t=1 Myr
double Westerlund2StarClusterUQFF::g_Westerlund2(double t) const {
    double M   = M_t(t);
    double g_n = G * M / (r_wd2 * r_wd2);
    const double B_crit = 1.0e11;
    double g_core = g_n * (1.0 + H0_wd2 * t) * (1.0 - B_wd2 / B_crit);
    double Ug1   = g_n;
    double Ug4   = Ug1 * (1.0 - B_wd2 / B_crit);
    double Ug_sum = (Ug1 + Ug4) * (1.0 + f_TRZ);
    double g_lambda = Lambda * c * c / 3.0;
    return g_core + Ug_sum + g_lambda + a_EM();
}

std::string Westerlund2StarClusterUQFF::primary_equation() const {
    return "g_W2(t)=G*M(t)/r^2*(1+H0*t)*(1-B/Bc)+(Ug1+Ug4)*(1+f_TRZ)+Lambda*c^2/3+q*v*B/m*(1+rho_UA/rho_SCm)*1e-12";
}

void Westerlund2StarClusterUQFF::self_update() {
    curr_t += time_step;
    M_init *= 1.0001; // slight mass gain from fallback
}

void Westerlund2StarClusterUQFF::self_expand() {
    r_wd2 *= 1.001;
}

void Westerlund2StarClusterUQFF::simulate(int num_steps) {
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        double g = g_Westerlund2(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\n";
    }
}

#ifdef STANDALONE_WESTERLUND2STARCLUSTERUQFF
int main() {
    Westerlund2StarClusterUQFF wd2;
    std::cout << "Westerlund 2 Star Cluster UQFF\n";
    std::cout << wd2.primary_equation() << "\n\n";
    std::cout << "g(1 Myr) = " << wd2.g_Westerlund2(3.156e13) << " m/s^2\n";
    wd2.simulate(3);
    return 0;
}
#endif
