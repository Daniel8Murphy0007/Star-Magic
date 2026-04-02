// STANDALONE_PILLARSOFCREATIONM16UQFF
#include "PillarsOfCreationM16UQFF.h"
#include <iostream>

double PillarsOfCreationM16UQFF::M_dot(double t) const {
    return M_0_SF * std::exp(-t / tau_SF);
}

double PillarsOfCreationM16UQFF::M_t(double t) const {
    return M_initial * (1.0 + M_dot(t));
}

// E(t) = E_0 * exp(-t/tau_erode) => 1-E at t=0.5 Myr ~ 0.9394
double PillarsOfCreationM16UQFF::one_minus_E(double t) const {
    return 1.0 - E_0 * std::exp(-t / tau_erode);
}

double PillarsOfCreationM16UQFF::H_z() const {
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + z_M16, 3) + 0.7);
}

double PillarsOfCreationM16UQFF::a_EM() const {
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_pillar / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}

// g_Pillars = G*M(t)/r^2*(1+H_z*t)*(1-B/B_crit)*(1-E(t))
//           + (Ug1+Ug4)*(1+f_TRZ) + Lambda*c^2/3 + a_EM
// ~ 1.053e-4 m/s^2 at t=0.5 Myr
double PillarsOfCreationM16UQFF::g_Pillars(double t) const {
    double M   = M_t(t);
    double g_n = G * M / (r_pillar * r_pillar);
    const double B_crit = 1.0e11; // T
    double g_core = g_n * (1.0 + H_z() * t) * (1.0 - B_pillar / B_crit)
                        * one_minus_E(t);
    // Ug1 + Ug4 = 2 * Ug1 (approximately)
    double Ug1 = g_n;
    double Ug4 = Ug1 * (1.0 - B_pillar / B_crit);
    double Ug_sum = (Ug1 + Ug4) * (1.0 + f_TRZ);
    double g_lambda = Lambda * c * c / 3.0;
    return g_core + Ug_sum + g_lambda + a_EM();
}

std::string PillarsOfCreationM16UQFF::primary_equation() const {
    return "g_Pillars(t)=G*M(t)/r^2*(1+H_z*t)*(1-B/Bc)*(1-E(t))+(Ug1+Ug4)*(1+f_TRZ)+Lambda*c^2/3+q*v*B/m*(1+rho_UA/rho_SCm)*1e-12";
}

void PillarsOfCreationM16UQFF::self_update() {
    curr_t += time_step;
    M_initial *= one_minus_E(curr_t);
}

void PillarsOfCreationM16UQFF::self_expand() {
    r_pillar *= 0.999; // pillar shortens as it erodes
    B_pillar *= 1.001;
}

void PillarsOfCreationM16UQFF::simulate(int num_steps) {
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        double g = g_Pillars(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\n";
    }
}

#ifdef STANDALONE_PILLARSOFCREATIONM16UQFF
int main() {
    PillarsOfCreationM16UQFF p;
    std::cout << "Pillars of Creation M16 UQFF\n";
    std::cout << p.primary_equation() << "\n\n";
    std::cout << "g(0.5 Myr) = " << p.g_Pillars(1.578e13) << " m/s^2\n";
    p.simulate(3);
    return 0;
}
#endif
