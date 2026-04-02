// STANDALONE_UQFFKNOWLEDGEBASEKB18
#include "UQFFKnowledgeBaseKB18.h"
#include <iostream>
#include <algorithm>

// Signal energy E = V_eff^2 / Z * dT_interval
double UQFFKnowledgeBaseKB18::signal_energy(int idx) const {
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    double V_eff = V_ch1[idx] / std::sqrt(2.0);
    const double Z = 50.0;
    return V_eff * V_eff / Z * dT_interv;
}

// U_m for signal at index idx, time offset t
double UQFFKnowledgeBaseKB18::U_m_signal(int idx, double t) const {
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    double A  = V_ch1[idx];
    double t_j = idx * dT_interv;
    double gamma = kappa;
    return mu_J * omega_THz * A * (1.0 - std::exp(-gamma * (t + t_j)));
}

// Flow state detection from Ch2 inversion
// Signals 44,47-48: chaotic (0); 45,46,49-50: inverted (-1); else normal (+1)
int UQFFKnowledgeBaseKB18::flow_state(int idx) const {
    // Chaotic: indices 3 (sig44), 6 (sig47), 7 (sig48)
    if (idx == 3 || idx == 6 || idx == 7) return 0;
    // Inverted: 4 (sig45), 5 (sig46), 8 (sig49), 9 (sig50)
    if (idx == 4 || idx == 5 || idx == 8 || idx == 9) return -1;
    return 1;
}

// Bundle sum: weighted by energy and flow states
double UQFFKnowledgeBaseKB18::bundle_sum() const {
    double sum = 0.0;
    for (int i = 0; i < n_signals_set; i++) {
        int fs = flow_state(i);
        if (fs != 0) {
            sum += signal_energy(i) * fs;
        }
    }
    return sum;
}

std::string UQFFKnowledgeBaseKB18::primary_equation() const {
    return "KB18_bundle = sum_i[V_eff_i^2/Z*dT*flow_state(i)]; U_m_i = mu_J*omega_THz*V_i*(1-exp(-kappa*(t+t_i)))";
}

void UQFFKnowledgeBaseKB18::self_update() {
    curr_t += time_step;
    curr_idx = (curr_idx + 1) % n_signals_set;
}

void UQFFKnowledgeBaseKB18::self_expand() {
    // No geometric expansion for KB; extend time window
}

void UQFFKnowledgeBaseKB18::simulate(int num_steps) {
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        int    idx = step % n_signals_set;
        double E   = signal_energy(idx);
        double Um  = U_m_signal(idx, t);
        int    fs  = flow_state(idx);
        std::cout << "Sig " << (idx + 41)
                  << "  t=" << (curr_t + idx * dT_interv) << " s"
                  << "  E=" << E << " J  U_m=" << Um
                  << "  flow=" << (fs > 0 ? "normal" : fs < 0 ? "inverted" : "chaotic") << "\n";
    }
}

#ifdef STANDALONE_UQFFKNOWLEDGEBASEKB18
int main() {
    UQFFKnowledgeBaseKB18 kb18;
    std::cout << "UQFF KB18 -- THz Signals 41-50 Analysis\n";
    std::cout << kb18.primary_equation() << "\n\n";
    std::cout << "Bundle sum = " << kb18.bundle_sum() << " J\n";
    kb18.simulate(5);
    return 0;
}
#endif
