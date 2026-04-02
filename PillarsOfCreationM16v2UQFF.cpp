// STANDALONE_PILLARSOFCREATIONM16V2UQFF
#include "PillarsOfCreationM16v2UQFF.h"
#include <iostream>

double PillarsOfCreationM16v2UQFF::M_dot(double t) const {
    return M_0_SF * std::exp(-t / tau_SF);
}

double PillarsOfCreationM16v2UQFF::M_t(double t) const {
    return M_initial * (1.0 + M_dot(t));
}

// Supernova shockwave disruption: stronger initial erosion
// E_shock(t) = E_0_shock * exp(-t/tau_shock)
double PillarsOfCreationM16v2UQFF::one_minus_E_shock(double t) const {
    return 1.0 - E_0_shock * std::exp(-t / tau_shock);
}

double PillarsOfCreationM16v2UQFF::H_z() const {
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + z_M16, 3) + 0.7);
}

// Protostar jet radiation pressure
// a_jet = L_jet / (c * M_t) ~ jet forcing per unit mass
double PillarsOfCreationM16v2UQFF::a_jet() const {
    return L_jet / (c * M_initial);
}

double PillarsOfCreationM16v2UQFF::a_EM_jet() const {
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_jet * B_jet / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}

// Full v2 MUGE: g_Pv2 = g_core*(1-E_shock) + (Ug1+Ug4)*(1+f_TRZ) + a_jet + a_EM
double PillarsOfCreationM16v2UQFF::g_Pillars_v2(double t) const {
    double M   = M_t(t);
    double g_n = G * M / (r_pillar * r_pillar);
    const double B_crit = 1.0e11;
    double g_core = g_n * (1.0 + H_z() * t) * (1.0 - B_jet / B_crit)
                        * one_minus_E_shock(t);
    double Ug1 = g_n;
    double Ug4 = Ug1;
    double Ug_sum = (Ug1 + Ug4) * (1.0 + f_TRZ);
    double g_lambda = Lambda * c * c / 3.0;
    return g_core + Ug_sum + g_lambda + a_jet() + a_EM_jet();
}

std::string PillarsOfCreationM16v2UQFF::primary_equation() const {
    return "g_Pv2(t)=G*M(t)/r^2*(1+H*t)*(1-B/Bc)*(1-E_shock)+(Ug1+Ug4)*(1+f_TRZ)+Lambda*c^2/3+L_jet/(c*M)+q*v_jet*B/m*(1+rho_UA/rho_SCm)*1e-12";
}

void PillarsOfCreationM16v2UQFF::self_update() {
    curr_t += time_step;
}

void PillarsOfCreationM16v2UQFF::self_expand() {
    r_pillar *= 0.9995; // SN shock collapses pillars faster
    L_jet    *= 0.999;
}

void PillarsOfCreationM16v2UQFF::simulate(int num_steps) {
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        double g = g_Pillars_v2(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\n";
    }
}

#ifdef STANDALONE_PILLARSOFCREATIONM16V2UQFF
int main() {
    PillarsOfCreationM16v2UQFF pv2;
    std::cout << "Pillars of Creation M16 v2 (Post-SN) UQFF\n";
    std::cout << pv2.primary_equation() << "\n\n";
    std::cout << "g(0.5 Myr) = " << pv2.g_Pillars_v2(1.578e13) << " m/s^2\n";
    pv2.simulate(3);
    return 0;
}
#endif
