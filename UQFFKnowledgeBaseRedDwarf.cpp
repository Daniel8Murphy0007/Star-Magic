// STANDALONE_UQFFKNOWLEDGEBASEREDDWARF
#include "UQFFKnowledgeBaseRedDwarf.h"
#include <iostream>
#include <numeric>

// |psi|^2 ~ A^2 * exp(-2*alpha*(r-r0)) * sin^2(k*r - omega*t) / r^2
double UQFFKnowledgeBaseRedDwarf::psi_magnitude(double r, double t) const {
    double r0   = 1.0e-10;
    double A    = 1.0;
    double sinK = std::sin(k_wave * r - omega_wave * t);
    double exp_d = std::exp(-2.0 * alpha_decay * std::abs(r - r0));
    return A * A * exp_d * sinK * sinK / (r * r + 1e-60);
}

// O_hat * psi approximation: lambda_I * omega * |psi|
double UQFFKnowledgeBaseRedDwarf::inertial_operator(double r, double t) const {
    return lambda_I * omega_wave * psi_magnitude(r, t);
}

// E_solfeggio = sum hbar * 2*pi * f_n (quantum energy of each Solfeggio harmonic)
double UQFFKnowledgeBaseRedDwarf::E_solfeggio_sum() const {
    double total = 0.0;
    for (double f : solfeggio) total += hbar * 2.0 * M_PI * f;
    return total;
}

// B_pseudo = mu_0 * q_m / (4*pi*r^2)
double UQFFKnowledgeBaseRedDwarf::B_pseudo(double r) const {
    return mu_0 * q_m / (4.0 * M_PI * r * r);
}

// P_DE_cosm = rho_SCm * c^2 * V / t_H
double UQFFKnowledgeBaseRedDwarf::P_DE_cosmological(double V_horizon) const {
    return rho_SCm * c * c * V_horizon / t_H;
}

// phi_twist = beta * sin(omega * t)  (KBI caduceus coil)
double UQFFKnowledgeBaseRedDwarf::phi_twist(double t) const {
    double beta = 0.9999;  // near-unity twist parameter
    return beta * std::sin(omega_wave * t);
}

// g_KB_composite: Combines inertial operator, solfeggio, B_pseudo, and UQFF
double UQFFKnowledgeBaseRedDwarf::g_KB_composite(double r, double t) const {
    double g_inertial = inertial_operator(r, t);
    double g_solfeggio = E_solfeggio_sum() * (rho_SCm / rho_UA);
    double g_B_pseudo  = B_pseudo(r) * k_wave;
    double g_DE        = P_DE * (1.0 + f_TRZ) / (k_B * 2.73);  // normalized by CMB
    return g_inertial + g_solfeggio + g_B_pseudo + g_DE;
}

std::string UQFFKnowledgeBaseRedDwarf::primary_equation() const {
    return "g_KB = O_psi(r,t) + E_solfeggio*(rho_SCm/rho_UA) + B_pseudo*k + P_DE*(1+f_TRZ)/(kB*T_CMB)";
}

void UQFFKnowledgeBaseRedDwarf::self_update() { curr_t += time_step; nu_plasma *= 1.001; }
void UQFFKnowledgeBaseRedDwarf::self_expand() { lambda_I *= 1.01; }
void UQFFKnowledgeBaseRedDwarf::simulate(int num_steps) {
    double r_test = 1.0e-12;  // quantum scale
    for (int i = 0; i < num_steps; i++) {
        self_update();
        std::cout << "t=" << curr_t << " s  g_KB=" << g_KB_composite(r_test, curr_t) << "\n";
    }
}

#ifdef STANDALONE_UQFFKNOWLEDGEBASEREDDWARF
int main() {
    UQFFKnowledgeBaseRedDwarf kb;
    std::cout << "UQFF Knowledge Base (Red Dwarf Compression Assimilation)\n";
    std::cout << kb.primary_equation() << "\n";
    std::cout << "E_solfeggio = " << kb.E_solfeggio_sum() << " J\n";
    std::cout << "P_DE_cosm = " << kb.P_DE_cosmological(4.0e26*4.0e26*4.0e26) << " W\n";
    kb.simulate(3);
    return 0;
}
#endif
