// STANDALONE_UQFFKNOWLEDGEBASEKB19
#include "UQFFKnowledgeBaseKB19.h"
#include <iostream>
#include <vector>

// Signal power P = V_eff^2 / Z ~ 0.00245 W
double UQFFKnowledgeBaseKB19::signal_power(double V_eff) const {
    return V_eff * V_eff / Z_line;
}

// U_m at THz frequency
// gamma ~ kappa for UQFF calibration
double UQFFKnowledgeBaseKB19::U_m_THz(double t) const {
    double gamma = kappa;
    double t_n   = 0.0;
    double P_SCm = rho_vac_SCm * 1.0e-3; // approximate phase space volume
    return mu_j_THz * omega_THz * (1.0 - std::exp(-gamma * t * std::cos(M_PI * t_n))) * P_SCm;
}

// Energy density from THz bundle
double UQFFKnowledgeBaseKB19::rho_E_THz(double V_eff) const {
    const double A_eff = 1.0e-4; // m^2 effective area
    return signal_power(V_eff) / (c * A_eff);
}

// ACE/DCE flow reversal phase
double UQFFKnowledgeBaseKB19::phi_inversion(double t) const {
    return std::cos(2.0 * M_PI * t / tau_flow);
}

// Bundle integral over N signals
double UQFFKnowledgeBaseKB19::I_THz_bundle(int N) const {
    double sum = 0.0;
    for (int j = 0; j < N; j++) {
        double t_j = j * dT_step;
        sum += U_m_THz(t_j) * phi_inversion(t_j);
    }
    return sum;
}

std::string UQFFKnowledgeBaseKB19::primary_equation() const {
    return "I_THz = sum_j[mu_j*omega_THz*(1-exp(-kappa*t*cos(pi*t_n)))*P_SCm]*cos(2pi*t/tau_flow); U_m=[U_m:SM_m / Ug1=UQFF_g+SM_g]^SCm";
}

void UQFFKnowledgeBaseKB19::self_update() {
    curr_t += time_step;
    curr_sig++;
    if (curr_sig > n_signals) curr_sig = 1;
}

void UQFFKnowledgeBaseKB19::self_expand() {
    // Extend analysis window
    n_signals += 10;
}

void UQFFKnowledgeBaseKB19::simulate(int num_steps) {
    std::cout << "THz Bundle Integral (KB19, " << n_signals << " signals):\n";
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        double I = I_THz_bundle(n_signals);
        double phi = phi_inversion(t);
        std::cout << "Step " << step
                  << "  t=" << t << " s"
                  << "  I_THz=" << I
                  << "  phi=" << phi << "\n";
    }
}

#ifdef STANDALONE_UQFFKNOWLEDGEBASEKB19
int main() {
    UQFFKnowledgeBaseKB19 kb19;
    std::cout << "UQFF KB19 -- THz Q-Scope Signals 1-50\n";
    std::cout << kb19.primary_equation() << "\n\n";
    std::cout << "P(V_eff=0.35) = " << kb19.signal_power(0.35) << " W\n";
    std::cout << "I_bundle(50) = " << kb19.I_THz_bundle(50) << "\n";
    kb19.simulate(3);
    return 0;
}
#endif
