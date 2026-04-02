// STANDALONE_NGC2014NGC2020STARFORMINGUQFF
#include "NGC2014NGC2020StarformingUQFF.h"
#include <iostream>

// Rapid secondary star formation from WR + OB stellar winds
// M_dot(t) = 41.67 * exp(-t/tau_SF); M(t=2.5 Myr) ~ 1.254e34 kg
double NGC2014NGC2020StarformingUQFF::M_dot(double t) const {
    return M_0_SF * std::exp(-t / tau_SF);
}

double NGC2014NGC2020StarformingUQFF::M_t(double t) const {
    return M_initial * (1.0 + M_dot(t));
}

double NGC2014NGC2020StarformingUQFF::a_EM() const {
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_LMC / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}

// g_Starbirth ~ 1.053e-4 m/s^2 at t=2.5 Myr
double NGC2014NGC2020StarformingUQFF::g_Starbirth(double t) const {
    double M   = M_t(t);
    double g_n = G * M / (r_region * r_region);
    const double B_crit = 1.0e11;
    double g_core = g_n * (1.0 + H0_LMC * t) * (1.0 - B_LMC / B_crit);
    double Ug1 = g_n;
    double Ug4 = Ug1 * (1.0 - B_LMC / B_crit);
    double Ug_sum = (Ug1 + Ug4) * (1.0 + f_TRZ);
    double g_lambda = Lambda * c * c / 3.0;
    return g_core + Ug_sum + g_lambda + a_EM();
}

std::string NGC2014NGC2020StarformingUQFF::primary_equation() const {
    return "g_Starbirth(t)=G*M(t)/r^2*(1+H0*t)*(1-B/Bc)+(Ug1+Ug4)*(1+f_TRZ)+Lambda*c^2/3+q*v*B/m*(1+rho_UA/rho_SCm)*1e-12";
}

void NGC2014NGC2020StarformingUQFF::self_update() {
    curr_t += time_step;
}

void NGC2014NGC2020StarformingUQFF::self_expand() {
    r_region *= 1.001;
    M_initial *= 1.001;
}

void NGC2014NGC2020StarformingUQFF::simulate(int num_steps) {
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        double g = g_Starbirth(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\n";
    }
}

#ifdef STANDALONE_NGC2014NGC2020STARFORMINGUQFF
int main() {
    NGC2014NGC2020StarformingUQFF tb;
    std::cout << "NGC 2014/2020 Tapestry of Blazing Starbirth UQFF\n";
    std::cout << tb.primary_equation() << "\n\n";
    std::cout << "g(2.5 Myr) = " << tb.g_Starbirth(7.89e13) << " m/s^2\n";
    tb.simulate(3);
    return 0;
}
#endif
