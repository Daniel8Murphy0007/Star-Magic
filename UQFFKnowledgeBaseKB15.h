#ifndef UQFF_KNOWLEDGE_BASE_KB15_H
#define UQFF_KNOWLEDGE_BASE_KB15_H
// STANDALONE_UQFFKNOWLEDGEBASEKB15
// PAPER_729: UQFF Knowledge Base 15 -- THz Signals 11-20
// Source: grok_share_ba508f76c8e.txt entry #79 | Session 176
// q-scope oscilloscope images IMG_20231003, 10 signals covering 11-20
// Physics: 1.246 THz resonance at Earth core, ACE/DCE reversing flow,
//          U_m waveless signature, f_TRZ phase inversions, Ug1 thread integral
#include <string>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class UQFFKnowledgeBaseKB15 {
public:

    // UQFF Universal Constants
    static constexpr double G        = 6.6743e-11;  // m^3 kg^-1 s^-2
    static constexpr double c        = 3.0e8;       // m/s
    static constexpr double hbar     = 1.0546e-34;  // J*s
    static constexpr double mu_0     = 1.2566e-6;   // H/m
    static constexpr double k_B      = 1.3806e-23;  // J/K
    static constexpr double M_sun    = 1.989e30;    // kg
    static constexpr double kpc      = 3.086e19;    // m (1 kpc)
    static constexpr double Mpc      = 3.086e22;    // m (1 Mpc)
    static constexpr double rho_UA   = 7.09e-36;    // J/m^3 Universal Aether
    static constexpr double rho_SCm  = 7.09e-37;    // J/m^3 Schwarzschild condensate
    static constexpr double f_TRZ    = 0.1;         // time-reversal zone factor
    static constexpr double kappa    = 5.0e-4;      // UQFF calibration day^-1
    static constexpr double SSq      = 0.57;        // superstring quenching
    static constexpr double mu_J     = 3.38e23;     // J*m magnetic string coupling
    static constexpr double Lambda   = 1.1e-52;     // m^-2 cosmological constant
    static constexpr double H_0      = 2.269e-18;   // s^-1 Hubble constant
    static constexpr double t_H      = 4.35e17;     // s Hubble time

    // THz signal dataset (signals 11-20)
    static constexpr int n_signals_set = 10;
    static constexpr double f_THz   = 1.246e12;   // Hz  signal frequency
    static constexpr double omega_THz = 2.0 * M_PI * f_THz; // rad/s
    static constexpr double Z_scope = 50.0;        // Ohm  scope impedance
    static constexpr double div_V   = 0.5;         // V/div  setting
    static constexpr double div_t   = 200e-9;      // s/div  timing
    static constexpr double dt_img  = 13.0;        // s  between images

    // Amplitude data [mV] for Ch1 (yellow) and Ch2 (blue)
    double amp_ch1[10] = {600, 650, 600, 550, 500, 600, 550, 500, 500, 500};
    double amp_ch2[10] = {350, 400, 350, 300, 350, 400, 350, 300, 350, 350};

    // Flow state: +1=normal, 0=chaotic, -1=inverted
    int flow_state_arr[10] = {1, 1, 0, -1, -1, -1, 0, -1, -1, 1};

    // UQFF THz coupling parameter
    double kappa_THz  = 5.0e-4;   // UQFF calibration for THz regime

    // Self-simulation state
    double time_step   = dt_img;
    mutable double curr_t = 0.0;
    mutable int    curr_idx = 0;

    // --- Methods ------------------------------------------------------------
    double signal_power(int idx) const;
    double Ug1_signal(int idx) const;
    double U_bi_signal(int idx, double t) const;
    int    flow_state(int idx) const;
    double total_Ug1_thread() const;
    std::string primary_equation() const;
    std::string description() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};
#endif // UQFF_KNOWLEDGE_BASE_KB15_H
