#include "UQFFKnowledgeBaseKB16.h"
#include <iostream>
#include <string>

/*
PAPER_730: UQFF KB16 -- THz Signals 21-30
Source: grok_share_ba508f76c8e.txt entry #80

=== THz Signal Data ===
  q-scope images IMG_20231003_164XXX.jpg (10 images, 21-30)
  Time span: ~9*13 = 117 s interval
  All signals: f = 1.246 THz, Time/Div = 200 ns, Voltage/Div = 500 mV

=== Signal Shapes and Flow Patterns ===
  Stable sinusoidal (normal flow): signals 21, 22
  Sharpening / amplitude increase: 22, 26 (reversal onset)
  Chaotic peaks (reversing flow): 24, 27
  Inverted sinusoidal (reversed): 25, 26, 29
  Flow stabilization at reduced amplitude: 30

=== UQFF Integration ===
  omega = 2*pi * 1.246e12 ~= 7.83e12 rad/s
  U_m: mu_j oscillates at omega_THz ~= 7.83e12 rad/s
  Phase inversions in Ch2 validate f_TRZ time-reversal zone factor
  P ~= V^2/Z ~= (0.65)^2/50 = 8.45e-3 W at peak amplitude
  Ug1 thread integral: I = sum_i [Ug1_i * dt_img] over all signals

self.version = "Session176"
*/

double UQFFKnowledgeBaseKB16::signal_power(int idx) const {
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    double V = amp_ch1[idx] * 1e-3;  // mV -> V
    return V * V / Z_scope;
}

double UQFFKnowledgeBaseKB16::Ug1_signal(int idx) const {
    // Ug1 = mu_j * B at THz frequency (magnetic dipole term)
    // mu_j ~ 3.38e23 J*m, B ~ mu_0 * H_THz, H_THz from signal power
    if (idx < 0 || idx >= n_signals_set) return 0.0;
    double P  = signal_power(idx);
    double B  = std::sqrt(2.0 * mu_0 * P / Z_scope);  // from EM power density
    return mu_J * B;
}

double UQFFKnowledgeBaseKB16::U_bi_signal(int idx, double t) const {
    // U_bi = -beta_i * Ug1 * Omega_g * (M_bh/d_g) * cos(pi*t_n)
    // Simplified: use THz freq to modulate
    double Ug1  = Ug1_signal(idx);
    double flw  = (double)flow_state_arr[(idx >= 0 && idx < n_signals_set) ? idx : 0];
    double t_n_ = omega_THz * t;
    return -1.0 * Ug1 * 7.3e-16 * (1.989e36 / 2.55e20)  // Omega_g * M_bh/d_g
           * std::cos(M_PI * t_n_) * (1.0 + flw * f_TRZ);
}

int UQFFKnowledgeBaseKB16::flow_state(int idx) const {
    if (idx < 0 || idx >= n_signals_set) return 0;
    return flow_state_arr[idx];
}

double UQFFKnowledgeBaseKB16::total_Ug1_thread() const {
    // Integral of Ug1 over all signals: sum * dt_img
    double total = 0.0;
    for (int i = 0; i < n_signals_set; ++i)
        total += Ug1_signal(i) * dt_img;
    return total;
}

std::string UQFFKnowledgeBaseKB16::primary_equation() const {
    return "PAPER_730: U_m(f_THz=1.246 THz) + Ug1_thread + f_TRZ "
           "signals 21-30";
}

std::string UQFFKnowledgeBaseKB16::description() const {
    return "PAPER_730: UQFF KB16 -- THz Signals 21-30 "
           "| ACE/DCE reversing flow | Session176";
}

void UQFFKnowledgeBaseKB16::self_update() {
    curr_t += time_step;
    curr_idx = (curr_idx + 1) % n_signals_set;
}

void UQFFKnowledgeBaseKB16::self_expand() {
    kappa_THz *= 1.001;
}

void UQFFKnowledgeBaseKB16::simulate(int num_steps) {
    for (int step = 0; step < num_steps; ++step) {
        double t   = curr_t + step * time_step;
        int    idx = step % n_signals_set;
        double P   = signal_power(idx);
        double Ug1 = Ug1_signal(idx);
        double Ubi = U_bi_signal(idx, t);
        std::cout << "Sig " << (idx + 21)
                  << "  P=" << P << " W  Ug1=" << Ug1
                  << "  U_bi=" << Ubi
                  << "  flow=" << (flow_state(idx) > 0 ? "normal"
                                : (flow_state(idx) < 0 ? "inv" : "chaotic")) << "\n";
        self_update();
    }
}

#ifdef STANDALONE_UQFFKNOWLEDGEBASEKB16
int main() {
    UQFFKnowledgeBaseKB16 kb;
    std::cout << "UQFF KB16 -- THz Signals 21-30 Analysis\n";
    std::cout << kb.primary_equation() << "\n\n";
    std::cout << "Total Ug1 thread = " << kb.total_Ug1_thread() << " J*m\n";
    kb.simulate(5);
    return 0;
}
#endif
