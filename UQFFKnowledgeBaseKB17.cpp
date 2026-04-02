// STANDALONE_UQFFKNOWLEDGEBASEKB17
#include "UQFFKnowledgeBaseKB17.h"
#include <iostream>

double UQFFKnowledgeBaseKB17::signal_power(int idx) const {
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    double V_eff = V_ch1[idx] / std::sqrt(2.0);
    return V_eff * V_eff / 50.0;
}

// U_bi: buoyancy adjustment for q-scope observing Earth's core THz hole
// U_bi ~ rho_UA * omega_i * V_ch1 * cos(2pi*t/tau) 
double UQFFKnowledgeBaseKB17::U_bi_signal(int idx, double t) const {
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    const double tau_flow = 52.0; // s
    double phi = std::cos(2.0 * M_PI * t / tau_flow);
    return rho_vac_UA * omega_i * V_ch1[idx] * phi;
}

// Ug1 thread strength: mu_J * omega * A^2
double UQFFKnowledgeBaseKB17::Ug1_signal(int idx) const {
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    double A = V_ch1[idx];
    return mu_J * omega_i * A * A;
}

// Classify flow state for signals 31-40 (0-indexed)
// Chaotic: 3 (sig34), 6 (sig37), 7 (sig38)
// Inverted: 4 (sig35), 5 (sig36), 8 (sig39), 9 (sig40)
int UQFFKnowledgeBaseKB17::flow_state(int idx) const {
    if (idx == 3 || idx == 6 || idx == 7) return 0;
    if (idx == 4 || idx == 5 || idx == 8 || idx == 9) return -1;
    return 1;
}

// Total Ug1 thread strength across all 10 signals in set
double UQFFKnowledgeBaseKB17::total_Ug1_thread() const {
    double sum = 0.0;
    for (int i = 0; i < n_signals_set; i++) {
        sum += Ug1_signal(i);
    }
    return sum;
}

std::string UQFFKnowledgeBaseKB17::primary_equation() const {
    return "Ug1_thread = sum_i[mu_J*omega_THz*V_i^2]; U_bi_i = rho_UA*omega_i*V_i*cos(2pi*t/tau_flow); [U_m:SM_m/Ug1=UQFF_g+SM_g]^SCm";
}

void UQFFKnowledgeBaseKB17::self_update() {
    curr_t += time_step;
    curr_idx = (curr_idx + 1) % n_signals_set;
}

void UQFFKnowledgeBaseKB17::self_expand() {
    // Extend UQFF resonance model
    kappa_THz *= 1.001;
}

void UQFFKnowledgeBaseKB17::simulate(int num_steps) {
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        int    idx = step % n_signals_set;
        double P   = signal_power(idx);
        double Ug1 = Ug1_signal(idx);
        double Ubi = U_bi_signal(idx, t);
        std::cout << "Sig " << (idx + 31)
                  << "  P=" << P << " W  Ug1=" << Ug1
                  << "  U_bi=" << Ubi
                  << "  flow=" << (flow_state(idx) > 0 ? "normal" : flow_state(idx) < 0 ? "inv" : "chaotic") << "\n";
    }
}

#ifdef STANDALONE_UQFFKNOWLEDGEBASEKB17
int main() {
    UQFFKnowledgeBaseKB17 kb17;
    std::cout << "UQFF KB17 -- THz Signals 31-40 Analysis\n";
    std::cout << kb17.primary_equation() << "\n\n";
    std::cout << "Total Ug1 thread = " << kb17.total_Ug1_thread() << " J*m\n";
    kb17.simulate(5);
    return 0;
}
#endif
