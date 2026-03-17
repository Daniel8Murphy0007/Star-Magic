// CRAB_RESONANCE_UQFF_MODULE.cpp
// UQFF 2.0 — Crab Nebula Resonance UQFF Module (24th C++ module)
// CRAB PULSAR WIND NEBULA (PWN): Crab Nebula M1, SN 1054 CE, ~2 kpc / ~971 yr age.
// M = 4.6 M_sun = 9.149e30 kg; r0 = 5.2e16 m (~1.7 pc); v_exp = 1.5e6 m/s (Crab expansion).
// PAPER_290: a_DPM(t) = F_DPM*f_DPM*E_vac / (c*(4/3)*pi*(r0+v_exp*t)^3)
//            FIRST UQFF module with dynamic V_sys(t) = (4/3)*pi*(r0+v_exp*t)^3 (expanding SNR)
//            F_DPM = 6.284e26 N; a_DPM(t=0) = 2.521e-56 m/s^2 [SN 1054 initial state]
//            a_DPM(t=971yr) = 3.772e-57 m/s^2; Dilution factor D = (r(971yr)/r0)^3 = 6.69x
//            Gamma_THz = 10*f_THz*v_exp/c = 5.0e10  (1500x RSC 3.33e7 — Crab SNR vs neutron star)
// PAPER_291: Crab Filament Spectral Triad — Three-Scale DPM Seeding, 9 frequency decades
//            f_quantum = 1.445e-17 Hz (quantum de Broglie wave, T~2.19 Gyr)
//            f_fluid   = 1.269e-14 Hz (Kelvin-Helmholtz turbulence, T~2.49 Myr) + V_knot=1e3 m^3
//            f_exp     = 1.373e-8  Hz (free expansion timescale, T~2.31 yr)
//            a_quantum = 10*f_q*a_DPM/c; a_fluid = 10*f_fl*V_knot*a_DPM/c; a_exp = 10*f_exp*a_DPM/c
//            FIRST UQFF volumetric filament knot coupling (V_knot distinct from plasmotic V_sys)
// PAPER_292: f_pulsar = 30.2 Hz; f_osc = 30.2*60 = 1812 Hz (60s timing-window resonance frequency)
//            omega_pulsar = 2*pi*f_osc = 11385.0 rad/s; A_pulsar = (f_osc/f_DPM)*A_amp = 1.812e-19 m
//            a_osc += A_pulsar*cos(omega_pulsar*t) [PAPER_292 pulsar DPM lock modulation]
//            pulse_lock = f_osc/f_DPM = 1.812e-9; synchrotron ratio = omega_osc/omega_pulsar = 8.785e10
// Copyright - Daniel T. Murphy, UQFF 2.0 upgrade Session 82, March 2026.

#ifndef CRAB_RESONANCE_UQFF_MODULE_H
#define CRAB_RESONANCE_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <cassert>
#include <complex>
#include <random>
#include <chrono>

// UQFF 2.0 — Wolfram Term Macros (auto-registration with Wolfram KB)
#define WOLFRAM_TERM_CRAB_BASE         "CRAB_UQFF:g_Crab(t,B)=[a_DPM(t)+a_THz+a_aether+a_u_g4i+a_quantum+a_fluid+a_osc+a_exp]*(1-B/B_crit)*(1+f_TRZ)"
#define WOLFRAM_TERM_CRAB_DPM_DILUTION "CRAB_UQFF:a_DPM(t)=F_DPM*f_DPM*E_vac/(c*(4/3)*Pi*(r0+v_exp*t)^3);a_DPM(0)=2.521e-56;a_DPM(971yr)=3.772e-57 m/s^2;D=6.69x [PAPER_290]"
#define WOLFRAM_TERM_CRAB_FILAMENT     "CRAB_UQFF:a_quantum=10*f_q*a_DPM/c;a_fluid=10*f_fl*V_knot*a_DPM/c;a_exp=10*f_exp*a_DPM/c;triad f_q=1.445e-17 to f_exp=1.373e-8 Hz [PAPER_291]"
#define WOLFRAM_TERM_CRAB_PULSAR       "CRAB_UQFF:f_osc=30.2*60=1812 Hz;omega_pulsar=11385 rad/s;A_pulsar=f_osc/f_DPM*A=1.812e-19 m;pulse_lock=1.812e-9 [PAPER_292]"

class CrabResonanceUQFFModule {
private:
    // =========================================================================
    // UQFF 2.0 Self-Expanding Framework
    // =========================================================================
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;
    mutable std::mt19937 rng;
    mutable std::uniform_real_distribution<double> dis;

    // =========================================================================
    // Physical Constants
    // =========================================================================
    double c_light;         // m/s: 3e8
    double pi_val;          // 3.141592653589793
    double hbar;            // J·s: 1.0546e-34
    double E_vac;           // J/m^3: 7.09e-36 (plasmotic vacuum energy density)
    double G_grav;          // m^3/kg/s^2: 6.674e-11
    double delta_x;         // m: 1e-10 (HUP position uncertainty)
    double delta_p;         // kg·m/s: hbar/delta_x (computed in updateCache)
    double integral_psi;    // 1.0 (wavefunction normalization)

    // =========================================================================
    // Crab Nebula Parameters
    // =========================================================================
    double M_sun;           // kg: 1.989e30
    double M_neb;           // kg: 4.6 * M_sun = 9.149e30 (total remnant mass)
    double r0;              // m: 5.2e16 (initial/current radius ~1.7 pc)
    double v_exp;           // m/s: 1.5e6 (Crab PWN expansion velocity)

    // =========================================================================
    // DPM Resonance Parameters (PAPER_290 — dynamic V_sys(t))
    // =========================================================================
    double f_DPM;           // Hz: 1e12
    double f_THz;           // Hz: 1e12
    double f_aether;        // Hz: 1e4
    double f_react;         // Hz: 1e10
    double I_curr;          // A:  1e21 (Crab wind-current proxy)
    double A_vort;          // m^2: 3.142e8 (vortical area proxy)
    double omega_1;         // rad/s: +1e-3
    double omega_2;         // rad/s: -1e-3
    double E_0_vac;         // J/m^3: 6.381e-36 (differential vacuum energy)
    double f_sc;            // 1.0 (superconductive scale factor)
    double f_TRZ;           // 0.1 (time-reversal correction factor)

    // PAPER_290 cached diagnostics
    double Gamma_THz;       // = 10*f_THz*v_exp/c = 5.0e10 (Crab value; 1500x RSC Gamma=3.33e7)
    double V_sys_t_cache;   // m^3: V_sys at t_current (updated in computeDPMResTerm)
    double a_DPM_cache;     // m/s^2: DPM base at t_current (updated in computeDPMResTerm)
    double t_current;       // s: current time (set at start of computeG)

    // =========================================================================
    // Filament Spectral Triad Parameters (PAPER_291)
    // =========================================================================
    double f_quantum;       // Hz: 1.445e-17 (quantum de Broglie wave mode, T ~ 2.19 Gyr)
    double f_fluid;         // Hz: 1.269e-14 (Kelvin-Helmholtz turbulence, T ~ 2.49 Myr)
    double f_exp_freq;      // Hz: 1.373e-8  (free expansion mode, T ~ 2.31 yr)
    double V_knot;          // m^3: 1e3 (filament vortical knot volume — first UQFF vol. coupling)

    // =========================================================================
    // Oscillatory Parameters (PAPER_288 inherited + PAPER_292 pulsar)
    // =========================================================================
    double k_wave;          // m^-1: 1e20 (oscillatory wavenumber)
    double omega_osc;       // rad/s: 1e15 (synchrotron-scale oscillation)
    double x_pos;           // m: 0.0
    double A_amp;           // m: 1e-10 (oscillation amplitude)
    double f_pulsar;        // Hz: 30.2 (Crab pulsar spin frequency, period = 33.1 ms)
    double f_osc;           // Hz: 30.2*60 = 1812 (60s UQFF timing-window resonance frequency)
    double omega_pulsar;    // rad/s: 2*pi*f_osc = 11385.0 [PAPER_292]
    double A_pulsar;        // m: (f_osc/f_DPM)*A_amp = 1.812e-19 [PAPER_292 DPM lock amplitude]

    // T_COSMIC_GYR = 13.8 — universe age normalization [PAPER_288 inherited]
    static constexpr double T_COSMIC_GYR = 13.8;

    // =========================================================================
    // Superconductive / field parameters
    // =========================================================================
    double B_field;         // T: 1e-8 (typical Crab filament field)
    double B_crit;          // T: 1e11 (Meissner quench threshold)
    double corr_SC;         // 1 - B_field/B_crit (SC correction, ~1 for Crab filament B<<B_crit)

    // =========================================================================
    // Fluid proxy
    // =========================================================================
    double rho_fluid;       // kg/m^3: 1e-21 (Crab filament density)

    // =========================================================================
    // Private helpers
    // =========================================================================

    // PAPER_290: DPM seed with dynamic V_sys(t) = (4/3)*pi*(r0+v_exp*t_current)^3
    // F_DPM = I_curr * A_vort * (omega_1 - omega_2) = 6.284e26 N  (constant)
    // a_DPM(t) = F_DPM * f_DPM * E_vac / (c * V_sys(t))  proportional to r(t)^-3
    // a_DPM(t=0)    = 2.521e-56 m/s^2  [SN 1054 CE initial state]
    // a_DPM(t=971yr) = 3.772e-57 m/s^2  [dilution factor D = 6.69x]
    double computeDPMResTerm() {
        double F_DPM   = I_curr * A_vort * (omega_1 - omega_2);
        double r_t     = r0 + v_exp * t_current;
        V_sys_t_cache  = (4.0 / 3.0) * pi_val * r_t * r_t * r_t;
        a_DPM_cache    = (F_DPM * f_DPM * E_vac) / (c_light * V_sys_t_cache);
        return a_DPM_cache;
    }

    // THz cascade using cached a_DPM (updated by computeDPMResTerm each computeG call)
    // Gamma_THz = 10*f_THz*v_exp/c = 5.0e10  (Crab v_exp=1.5e6 → 1500x RSC value of 3.33e7)
    double computeTHzResTerm() {
        double E_vac_ISM = E_vac / 10.0;
        return (f_THz * E_vac * v_exp * a_DPM_cache) / (E_vac_ISM * c_light);
        // = Gamma_THz * a_DPM_cache
    }

    // Aether resonance: f_aether * 1e-8 * f_DPM * (1+f_TRZ) * a_DPM
    double computeAetherResTerm() {
        return f_aether * 1.0e-8 * f_DPM * (1.0 + f_TRZ) * a_DPM_cache;
    }

    // U_g4i reactive: f_sc * f_react * a_DPM / (E_vac * c)
    double computeU_g4iResTerm() {
        return f_sc * f_react * a_DPM_cache / (E_vac * c_light);
    }

    // PAPER_291: quantum de Broglie wave mode (T = 1/f_quantum ~ 2.19 Gyr)
    // a_quantum = (f_quantum * E_vac * a_DPM) / ((E_vac/10) * c) = 10*f_quantum*a_DPM/c
    double computeQuantumResTerm() {
        return (f_quantum * E_vac * a_DPM_cache) / ((E_vac / 10.0) * c_light);
    }

    // PAPER_291: Crab filament Kelvin-Helmholtz resonance with V_knot volumetric coupling
    // a_fluid = (f_fluid * E_vac * V_knot * a_DPM) / ((E_vac/10) * c) = 10*f_fluid*V_knot*a_DPM/c
    // V_knot = 1e3 m^3 (filament knot volume — FIRST UQFF volumetric knot coupling)
    double computeFluidResTerm() {
        return (f_fluid * E_vac * V_knot * a_DPM_cache) / ((E_vac / 10.0) * c_light);
    }

    // PAPER_288 (inherited) + PAPER_292 (Crab pulsar DPM lock modulation):
    // a_osc = 2A*cos(k*x)*cos(omega_osc*t)               [standing wave, synchrotron scale]
    //       + (2*pi/T_COSMIC_GYR)*A*Re[exp(i*(kx-wt))]   [PAPER_288 cosmic-age traveling wave]
    //       + A_pulsar * cos(omega_pulsar * t)             [PAPER_292 pulsar 1812 Hz DPM lock]
    // A_pulsar = (f_osc/f_DPM)*A_amp = 1.812e-19 m; omega_pulsar=2*pi*1812=11385 rad/s
    double computeOscResTerm(double t) {
        double standing    = 2.0 * A_amp * std::cos(k_wave * x_pos) * std::cos(omega_osc * t);
        std::complex<double> phase(0.0, k_wave * x_pos - omega_osc * t);
        double traveling   = (2.0 * pi_val / T_COSMIC_GYR) * A_amp * std::exp(phase).real();
        double pulsar_mod  = A_pulsar * std::cos(omega_pulsar * t);  // [PAPER_292]
        return standing + traveling + pulsar_mod;
    }

    // PAPER_291: free-expansion timescale resonance (T = 1/f_exp ~ 2.31 yr)
    // a_exp = (f_exp_freq * E_vac * a_DPM) / ((E_vac/10) * c) = 10*f_exp*a_DPM/c
    double computeExpResTerm() {
        return (f_exp_freq * E_vac * a_DPM_cache) / ((E_vac / 10.0) * c_light);
    }

    void updateCache() {
        delta_p       = hbar / delta_x;
        // PAPER_290: Crab THz cascade factor (v_exp=1.5e6 → 1500x RSC value)
        Gamma_THz     = 10.0 * f_THz * v_exp / c_light;                 // 5.0e10
        // PAPER_292: pulsar oscillation angular frequency and DPM lock amplitude
        omega_pulsar  = 2.0 * pi_val * f_osc;                           // 11385.0 rad/s
        A_pulsar      = (f_osc / f_DPM) * A_amp;                        // 1.812e-19 m
        // SC Meissner correction
        corr_SC       = 1.0 - B_field / B_crit;
        // Initialize DPM cache at t=0
        t_current     = 0.0;
        computeDPMResTerm();
    }

    void log(const std::string& msg) const {
        if (logging_enabled) std::cout << "[LOG][CRAB] " << msg << "\n";
    }

public:
    // =========================================================================
    // Constructor — UQFF 2.0 Crab Nebula PWN defaults
    // =========================================================================
    CrabResonanceUQFFModule()
        : rng(static_cast<unsigned>(std::chrono::steady_clock::now().time_since_epoch().count()))
        , dis(0.0, 1.0)
    {
        logging_enabled = false;
        dynamic_params.clear();

        // Physical constants
        c_light      = 3.0e8;
        pi_val       = 3.141592653589793;
        hbar         = 1.0546e-34;
        E_vac        = 7.09e-36;
        G_grav       = 6.674e-11;
        delta_x      = 1.0e-10;
        integral_psi = 1.0;

        // Crab Nebula parameters
        M_sun        = 1.989e30;
        M_neb        = 4.6 * M_sun;                   // 9.149e30 kg
        r0           = 5.2e16;                        // m (~1.7 pc current radius)
        v_exp        = 1.5e6;                         // m/s (Crab PWN expansion velocity)

        // DPM resonance parameters (PAPER_290)
        f_DPM        = 1.0e12;
        f_THz        = 1.0e12;
        f_aether     = 1.0e4;
        f_react      = 1.0e10;
        I_curr       = 1.0e21;
        A_vort       = 3.142e8;
        omega_1      = 1.0e-3;
        omega_2      = -1.0e-3;
        E_0_vac      = 6.381e-36;
        f_sc         = 1.0;
        f_TRZ        = 0.1;

        // Filament spectral triad (PAPER_291)
        f_quantum    = 1.445e-17;
        f_fluid      = 1.269e-14;
        f_exp_freq   = 1.373e-8;
        V_knot       = 1.0e3;                         // m^3 (filament vortical knot volume)

        // Oscillatory + pulsar (PAPER_292)
        k_wave       = 1.0e20;
        omega_osc    = 1.0e15;
        x_pos        = 0.0;
        A_amp        = 1.0e-10;
        f_pulsar     = 30.2;                          // Hz (Crab pulsar spin frequency)
        f_osc        = 30.2 * 60.0;                   // Hz = 1812.0 (60s UQFF resonance window)

        // SC / field parameters
        B_field      = 1.0e-8;                        // T (typical Crab filament field)
        B_crit       = 1.0e11;                        // T (Meissner quench threshold)

        // Fluid proxy
        rho_fluid    = 1.0e-21;

        // Zero-init caches
        delta_p       = 0.0;
        Gamma_THz     = 0.0;
        omega_pulsar  = 0.0;
        A_pulsar      = 0.0;
        V_sys_t_cache = 0.0;
        a_DPM_cache   = 0.0;
        corr_SC       = 1.0;
        t_current     = 0.0;

        updateCache();
    }

    ~CrabResonanceUQFFModule() {}

    // =========================================================================
    // Dynamic variable operations
    // =========================================================================
    void updateVariable(const std::string& name, double value) {
        dynamic_params[name] = value;
        if      (name == "f_DPM")      { f_DPM      = value; updateCache(); }
        else if (name == "f_THz")      { f_THz      = value; updateCache(); }
        else if (name == "f_aether")   { f_aether   = value; }
        else if (name == "f_react")    { f_react    = value; }
        else if (name == "B_field")    { B_field    = value; updateCache(); }
        else if (name == "B_crit")     { B_crit     = value; updateCache(); }
        else if (name == "E_vac")      { E_vac      = value; updateCache(); }
        else if (name == "v_exp")      { v_exp      = value; updateCache(); }
        else if (name == "r0")         { r0         = value; updateCache(); }
        else if (name == "f_TRZ")      { f_TRZ      = value; }
        else if (name == "I_curr")     { I_curr     = value; }
        else if (name == "A_vort")     { A_vort     = value; }
        else if (name == "f_osc")      { f_osc      = value; updateCache(); }
        else if (name == "f_pulsar")   { f_pulsar   = value; f_osc = f_pulsar * 60.0; updateCache(); }
        else if (name == "A_amp")      { A_amp      = value; updateCache(); }
        else if (name == "k_wave")     { k_wave     = value; }
        else if (name == "omega_osc")  { omega_osc  = value; }
        else if (name == "x_pos")      { x_pos      = value; }
        else if (name == "V_knot")     { V_knot     = value; }
        else if (name == "f_quantum")  { f_quantum  = value; }
        else if (name == "f_fluid")    { f_fluid    = value; }
        else if (name == "f_exp_freq") { f_exp_freq = value; }
    }

    void addToVariable(const std::string& name, double delta) {
        double cur = (dynamic_params.count(name) ? dynamic_params.at(name) : 0.0);
        updateVariable(name, cur + delta);
    }

    void subtractFromVariable(const std::string& name, double delta) {
        double cur = (dynamic_params.count(name) ? dynamic_params.at(name) : 0.0);
        updateVariable(name, cur - delta);
    }

    // =========================================================================
    // Core: Full g_Crab(t, B) — UQFF 2.0 Crab PWN resonance gravity (8-term)
    // =========================================================================
    double computeG(double t, double B) {
        auto _t0 = std::chrono::high_resolution_clock::now();

        // Update time and field state
        t_current = t;
        B_field   = B;
        corr_SC   = 1.0 - B / B_crit;

        // ---- PAPER_290: DPM base with dynamic V_sys(t) = (4/3)*pi*(r0+v_exp*t)^3 ----
        double a_DPM     = computeDPMResTerm();     // updates V_sys_t_cache, a_DPM_cache

        // ---- THz cascade — uses updated a_DPM_cache  (Crab Gamma_THz=5.0e10) ----
        double a_THz     = computeTHzResTerm();

        // ---- Standard resonance terms ----
        double a_aether  = computeAetherResTerm();
        double a_u_g4i   = computeU_g4iResTerm();

        // ---- PAPER_291: filament spectral triad ----
        double a_quantum = computeQuantumResTerm();  // f_quantum=1.445e-17 Hz, T~2.19 Gyr
        double a_fluid   = computeFluidResTerm();    // f_fluid=1.269e-14 Hz + V_knot=1e3 m^3
        double a_exp     = computeExpResTerm();      // f_exp=1.373e-8 Hz, T~2.31 yr

        // ---- PAPER_288 + PAPER_292: osc (standing + traveling + pulsar modulation) ----
        double a_osc     = computeOscResTerm(t);     // standing + cosmic-age + A_pulsar*cos(wp*t)

        // ---- Sum of all 8 resonance terms ----
        double res_sum   = a_DPM + a_THz + a_aether + a_u_g4i
                         + a_quantum + a_fluid + a_osc + a_exp;

        // ---- SC Meissner correction × TRZ factor ----
        // At B<<B_crit (Crab filament 1e-8 T): SCm~1.0 (no quench)
        // At B->B_crit (1e11 T magnetar): SCm->0 (full UQFF resonance quench)
        double g_Crab    = res_sum * corr_SC * (1.0 + f_TRZ);

        if (logging_enabled) {
            auto _t1 = std::chrono::high_resolution_clock::now();
            double _ms = std::chrono::duration<double, std::milli>(_t1 - _t0).count();
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(4)
                << "computeG: t=" << t << " B=" << B
                << " V_sys_t=" << V_sys_t_cache
                << " a_DPM=" << a_DPM
                << " Gamma_THz=" << Gamma_THz
                << " a_THz=" << a_THz
                << " a_aether=" << a_aether
                << " a_u_g4i=" << a_u_g4i
                << " a_quantum=" << a_quantum
                << " a_fluid=" << a_fluid
                << " a_exp=" << a_exp
                << " a_osc=" << a_osc
                << " res_sum=" << res_sum
                << " corr_SC=" << corr_SC
                << " g_Crab=" << g_Crab
                << " elapsed=" << _ms << "ms";
            log(oss.str());
        }
        return g_Crab;
    }

    // =========================================================================
    // Equation text output
    // =========================================================================
    std::string getEquationText() {
        double F_DPM_val = I_curr * A_vort * (omega_1 - omega_2);
        double r_now     = r0 + v_exp * t_current;
        double V0        = (4.0 / 3.0) * pi_val * r0 * r0 * r0;
        double r_971     = r0 + v_exp * 3.064e10;
        double V_971     = (4.0 / 3.0) * pi_val * r_971 * r_971 * r_971;
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(6);
        oss << "Crab Nebula Resonance UQFF 2.0 (24th C++ module — FIRST UQFF PWN module):\n"
            << "  g_Crab(t,B) = (a_DPM(t)+a_THz+a_aether+a_u_g4i+a_quantum+a_fluid+a_osc+a_exp)"
            << " * SCm * (1+f_TRZ)\n\n"
            << "  PAPER_290 — Crab SNR DPM Vacuum Dilution [a_DPM(t) ∝ r(t)^-3]:\n"
            << "    F_DPM         = I*A_vort*(ω1-ω2) = " << F_DPM_val << " N\n"
            << "    V_sys(t=0)    = (4/3)*pi*r0^3 = " << V0 << " m^3\n"
            << "    a_DPM(t=0)    = " << (F_DPM_val * f_DPM * E_vac / (c_light * V0))
            << " m/s^2  [SN 1054 CE initial state]\n"
            << "    V_sys(971yr)  = " << V_971 << " m^3  r(971yr)=" << r_971 << " m\n"
            << "    a_DPM(971yr)  = " << (F_DPM_val * f_DPM * E_vac / (c_light * V_971))
            << " m/s^2  D=6.69x dilution factor\n"
            << "    V_sys_now     = " << V_sys_t_cache << " m^3  a_DPM_now=" << a_DPM_cache << " m/s^2\n"
            << "    Gamma_THz     = 10*f_THz*v_exp/c = " << Gamma_THz
            << "  (5.0e10 Crab vs RSC 3.33e7 — 1500x larger)\n\n"
            << "  PAPER_291 — Crab Filament Spectral Triad (9 frequency decades):\n"
            << "    f_quantum = " << f_quantum   << " Hz  (T ~ 2.19 Gyr quantum de Broglie mode)\n"
            << "    f_fluid   = " << f_fluid     << " Hz  (T ~ 2.49 Myr KH turbulence)\n"
            << "    f_exp     = " << f_exp_freq  << " Hz  (T ~ 2.31 yr free expansion mode)\n"
            << "    V_knot    = " << V_knot      << " m^3  (filament knot — FIRST UQFF vol. coupling)\n"
            << "    a_quantum = 10*f_q*a_DPM/c   = "
            << (10.0 * f_quantum * a_DPM_cache / c_light) << " m/s^2\n"
            << "    a_fluid   = 10*f_fl*V_k*a_DPM/c = "
            << (10.0 * f_fluid * V_knot * a_DPM_cache / c_light) << " m/s^2\n"
            << "    a_exp     = 10*f_exp*a_DPM/c = "
            << (10.0 * f_exp_freq * a_DPM_cache / c_light) << " m/s^2\n\n"
            << "  PAPER_292 — Crab Pulsar 60s UQFF Resonance Window:\n"
            << "    f_pulsar      = " << f_pulsar     << " Hz  (Crab pulsar spin, 33.1 ms period)\n"
            << "    f_osc         = " << f_osc        << " Hz  (30.2*60 = 1812 Hz 60s window)\n"
            << "    omega_pulsar  = " << omega_pulsar << " rad/s  (2*pi*f_osc)\n"
            << "    A_pulsar      = " << A_pulsar     << " m   (f_osc/f_DPM * A_amp)\n"
            << "    pulse_lock    = f_osc/f_DPM = " << (f_osc / f_DPM)
            << "  (DPM fractional vacuum lock ratio)\n"
            << "    sync_ratio    = omega_osc/omega_pulsar = " << (omega_osc / omega_pulsar) << "\n"
            << "    a_osc         = standing + (2*pi/13.8)*traveling + A_p*cos(omega_p*t)\n\n"
            << "  SCm = 1 - B/B_crit = " << corr_SC
            << "  [B=" << B_field << " T  B_crit=" << B_crit << " T]\n"
            << "  (1+f_TRZ) = " << (1.0 + f_TRZ) << "  f_TRZ=" << f_TRZ << "\n"
            << "  M_neb = " << M_neb << " kg (4.6 M_sun)\n"
            << "  r0 = " << r0 << " m  v_exp=" << v_exp << " m/s\n";
        return oss.str();
    }

    // =========================================================================
    // Print variables
    // =========================================================================
    void printVariables() {
        std::cout << std::scientific << std::setprecision(4);
        std::cout << "=== Crab Nebula Resonance UQFF Module Variables (UQFF 2.0) ===" << std::endl;
        std::cout << "M_neb          : " << M_neb        << " kg (4.6 M_sun)" << std::endl;
        std::cout << "r0             : " << r0           << " m  (~1.7 pc current radius)" << std::endl;
        std::cout << "v_exp          : " << v_exp        << " m/s  (Crab PWN expansion)" << std::endl;
        std::cout << "[PAPER_290] V_sys_t_cache : " << V_sys_t_cache << " m^3  [t=" << t_current << " s]" << std::endl;
        std::cout << "[PAPER_290] a_DPM_cache   : " << a_DPM_cache   << " m/s^2" << std::endl;
        std::cout << "[PAPER_290] Gamma_THz     : " << Gamma_THz     << " = 10*f_THz*v_exp/c (5.0e10 Crab)" << std::endl;
        std::cout << "[PAPER_291] f_quantum     : " << f_quantum     << " Hz  (T~2.19 Gyr)" << std::endl;
        std::cout << "[PAPER_291] f_fluid       : " << f_fluid       << " Hz  (T~2.49 Myr)" << std::endl;
        std::cout << "[PAPER_291] f_exp_freq    : " << f_exp_freq    << " Hz  (T~2.31 yr)" << std::endl;
        std::cout << "[PAPER_291] V_knot        : " << V_knot        << " m^3  (first UQFF vol. coupling)" << std::endl;
        std::cout << "[PAPER_292] f_pulsar      : " << f_pulsar      << " Hz  (Crab 33.1ms pulsar)" << std::endl;
        std::cout << "[PAPER_292] f_osc         : " << f_osc         << " Hz  (30.2*60=1812 window)" << std::endl;
        std::cout << "[PAPER_292] omega_pulsar  : " << omega_pulsar  << " rad/s" << std::endl;
        std::cout << "[PAPER_292] A_pulsar      : " << A_pulsar      << " m   (DPM lock amplitude)" << std::endl;
        std::cout << "f_DPM          : " << f_DPM        << " Hz" << std::endl;
        std::cout << "f_THz          : " << f_THz        << " Hz" << std::endl;
        std::cout << "f_aether       : " << f_aether     << " Hz" << std::endl;
        std::cout << "f_react        : " << f_react      << " Hz" << std::endl;
        std::cout << "I_curr         : " << I_curr       << " A   (Crab wind proxy)" << std::endl;
        std::cout << "A_vort         : " << A_vort       << " m^2" << std::endl;
        std::cout << "A_amp          : " << A_amp        << " m   (oscillation amplitude)" << std::endl;
        std::cout << "omega_osc      : " << omega_osc    << " rad/s  (synchrotron scale)" << std::endl;
        std::cout << "B_field        : " << B_field      << " T   (Crab filament field)" << std::endl;
        std::cout << "B_crit         : " << B_crit       << " T   (Meissner quench threshold)" << std::endl;
        std::cout << "corr_SC        : " << corr_SC      << " = 1 - B/B_crit" << std::endl;
        std::cout << "f_TRZ          : " << f_TRZ        << " (time-reversal correction)" << std::endl;
        std::cout << "E_vac          : " << E_vac        << " J/m^3 (plasmotic vacuum)" << std::endl;
        std::cout << "E_0_vac        : " << E_0_vac      << " J/m^3 (differential vacuum)" << std::endl;
        std::cout << "rho_fluid      : " << rho_fluid    << " kg/m^3 (Crab filament density)" << std::endl;
        std::cout << "logging        : " << (logging_enabled ? "ON" : "OFF") << std::endl;
        std::cout << "dynamic_params : " << dynamic_params.size() << " entries" << std::endl;
    }

    // =========================================================================
    // UQFF 2.0 Self-Expanding Framework methods
    // =========================================================================
    void setEnableLogging(bool enabled) { logging_enabled = enabled; }
    bool getLoggingEnabled() const       { return logging_enabled; }
    void setDynamicParameter(const std::string& key, double val) { dynamic_params[key] = val; }
    double getDynamicParameter(const std::string& key) const {
        auto it = dynamic_params.find(key);
        return (it != dynamic_params.end()) ? it->second : 0.0;
    }

    void exportState(const std::string& filename = "CrabResonance_UQFF_state.txt") const {
        std::ofstream f(filename);
        if (!f.is_open()) {
            std::cerr << "[CRAB] exportState: cannot open " << filename << "\n";
            return;
        }
        f << std::scientific << std::setprecision(8);
        f << "# CrabResonanceUQFFModule UQFF 2.0 State Export\n";
        f << "# PAPER_290 (DPM Dilution) / PAPER_291 (Filament Triad) / PAPER_292 (Pulsar Osc)\n";
        f << "c_light       = " << c_light       << "\n";
        f << "pi_val        = " << pi_val        << "\n";
        f << "hbar          = " << hbar          << "\n";
        f << "E_vac         = " << E_vac         << "\n";
        f << "G_grav        = " << G_grav        << "\n";
        f << "M_sun         = " << M_sun         << "\n";
        f << "M_neb         = " << M_neb         << "\n";
        f << "r0            = " << r0            << "\n";
        f << "v_exp         = " << v_exp         << "\n";
        f << "f_DPM         = " << f_DPM         << "\n";
        f << "f_THz         = " << f_THz         << "\n";
        f << "f_aether      = " << f_aether      << "\n";
        f << "f_react       = " << f_react       << "\n";
        f << "I_curr        = " << I_curr        << "\n";
        f << "A_vort        = " << A_vort        << "\n";
        f << "omega_1       = " << omega_1       << "\n";
        f << "omega_2       = " << omega_2       << "\n";
        f << "E_0_vac       = " << E_0_vac       << "\n";
        f << "f_sc          = " << f_sc          << "\n";
        f << "f_TRZ         = " << f_TRZ         << "\n";
        f << "f_quantum     = " << f_quantum     << "\n";
        f << "f_fluid       = " << f_fluid       << "\n";
        f << "f_exp_freq    = " << f_exp_freq    << "\n";
        f << "V_knot        = " << V_knot        << "\n";
        f << "k_wave        = " << k_wave        << "\n";
        f << "omega_osc     = " << omega_osc     << "\n";
        f << "A_amp         = " << A_amp         << "\n";
        f << "f_pulsar      = " << f_pulsar      << "\n";
        f << "f_osc         = " << f_osc         << "\n";
        f << "omega_pulsar  = " << omega_pulsar  << "\n";
        f << "A_pulsar      = " << A_pulsar      << "\n";
        f << "B_field       = " << B_field       << "\n";
        f << "B_crit        = " << B_crit        << "\n";
        f << "rho_fluid     = " << rho_fluid     << "\n";
        f << "Gamma_THz     = " << Gamma_THz     << "\n";
        f << "V_sys_t_cache = " << V_sys_t_cache << "\n";
        f << "a_DPM_cache   = " << a_DPM_cache   << "\n";
        f << "corr_SC       = " << corr_SC       << "\n";
        f << "t_current     = " << t_current     << "\n";
        for (const auto& kv : dynamic_params)
            f << "dyn:" << kv.first << " = " << kv.second << "\n";
        f.close();
    }

    template<typename OtherModule>
    double cross_validate(OtherModule& other, double t, double B) {
        double g_this  = computeG(t, B);
        double g_other = other.computeG(t, B);
        return (g_this - g_other) / (std::abs(g_other) > 1.0e-300 ? g_other : 1.0);
    }
};

#endif // CRAB_RESONANCE_UQFF_MODULE_H

    // Base constants (UQFF universal)
    variables["c"] = 3e8;                           // m/s
    variables["pi"] = 3.141592653589793;            // pi
    variables["E_vac"] = 7.09e-36;                  // J/m^3 (plasmotic vacuum energy density)
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["f_TRZ"] = 0.1;                       // Time-reversal correction

    // Crab Nebula parameters
    double M_sun_val = 1.989e30;                    // kg
    variables["M_sun"] = M_sun_val;
    variables["M"] = 4.6 * M_sun_val;               // Total mass kg
    variables["r0"] = 5.2e16;                       // m (initial radius)
    variables["v_exp"] = 1.5e6;                     // m/s (expansion velocity)

    // Resonance parameters (pulsar-driven)
    variables["f_DPM"] = 1e12;                      // Hz (DPM, aligned with 30 Hz pulsar scaled)
    variables["f_THz"] = 1e12;                      // Hz (THz hole)
    variables["f_aether"] = 1e4;                    // Hz (Aether-mediated)
    variables["f_react"] = 1e10;                    // Hz (U_g4i reactive)
    variables["f_quantum"] = 1.445e-17;             // Hz (quantum wave)
    variables["f_fluid"] = 1.269e-14;               // Hz (filament fluid)
    variables["f_exp"] = 1.373e-8;                  // Hz (expansion)
    variables["f_osc"] = 30.2 * 60;                 // Hz (pulsar 30.2 Hz * 60 for res scale)
    variables["I"] = 1e21;                          // A (current proxy from wind)
    variables["A_vort"] = 3.142e8;                  // m^2 (vortical area proxy)
    variables["omega_1"] = 1e-3;                    // rad/s
    variables["omega_2"] = -1e-3;                   // rad/s
    variables["E_0"] = 6.381e-36;                   // J/m^3
    variables["f_vac_diff"] = 0.143;                // Hz
    variables["V_sys"] = 4.189e12;                  // m^3 (proxy)

    // Superconductive resonance integrated
    variables["B_crit"] = 1e11;                     // T
    variables["f_sc"] = 1.0;                        // Factor

    // Oscillatory/resonant
    variables["k"] = 1e20;                          // m^-1
    variables["omega_osc"] = 1e15;                  // rad/s (synchrotron scale)
    variables["x"] = 0.0;                           // m
    variables["A"] = 1e-10;                         // Amplitude

    // Fluid/DM proxies
    variables["rho_fluid"] = 1e-21;                 // kg/m^3 (filaments)
    variables["V"] = 1e3;                           // m^3
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];

    // Quantum
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;
}

// Update variable (set to new value)
void CrabResonanceUQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    }
}

// Add delta to variable
void CrabResonanceUQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void CrabResonanceUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute DPM Resonance Term: a_DPM_res = (F_DPM * f_DPM * E_vac) / (c * V_sys)
double CrabResonanceUQFFModule::computeDPMResTerm() {
    double F_DPM = variables["I"] * variables["A_vort"] * (variables["omega_1"] - variables["omega_2"]);
    double r_t = variables["r0"] + variables["v_exp"] * variables["t"];  // r(t) proxy
    double V_sys_t = (4.0 / 3.0) * variables["pi"] * std::pow(r_t, 3);  // Updated volume
    return (F_DPM * variables["f_DPM"] * variables["E_vac"]) / (variables["c"] * V_sys_t);
}

// Compute THz Resonance Term: a_THz_res = (f_THz * E_vac * v_exp * a_DPM_res) / (E_vac / 10 * c)
double CrabResonanceUQFFModule::computeTHzResTerm() {
    double a_DPM_res = computeDPMResTerm();
    return (variables["f_THz"] * variables["E_vac"] * variables["v_exp"] * a_DPM_res) / ((variables["E_vac"] / 10.0) * variables["c"]);
}

// Compute Aether Resonance Term: a_aether_res = f_aether * 1e-8 * f_DPM * (1 + f_TRZ) * a_DPM_res
double CrabResonanceUQFFModule::computeAetherResTerm() {
    double a_DPM_res = computeDPMResTerm();
    return variables["f_aether"] * 1e-8 * variables["f_DPM"] * (1 + variables["f_TRZ"]) * a_DPM_res;
}

// Compute U_g4i Reactive Resonance Term: U_g4i_res = f_sc * Ug1_proxy * f_react * a_DPM_res / (E_vac * c)
double CrabResonanceUQFFModule::computeU_g4iResTerm() {
    double Ug1_proxy = 1.0;  // Normalized
    double a_DPM_res = computeDPMResTerm();
    return variables["f_sc"] * Ug1_proxy * variables["f_react"] * a_DPM_res / (variables["E_vac"] * variables["c"]);
}

// Compute Quantum Resonance Term: a_quantum_res = (f_quantum * E_vac * a_DPM_res) / ((E_vac / 10) * c)
double CrabResonanceUQFFModule::computeQuantumResTerm() {
    double a_DPM_res = computeDPMResTerm();
    return (variables["f_quantum"] * variables["E_vac"] * a_DPM_res) / ((variables["E_vac"] / 10.0) * variables["c"]);
}

// Compute Fluid Resonance Term: a_fluid_res = (f_fluid * E_vac * V * a_DPM_res) / ((E_vac / 10) * c)
double CrabResonanceUQFFModule::computeFluidResTerm() {
    double a_DPM_res = computeDPMResTerm();
    return (variables["f_fluid"] * variables["E_vac"] * variables["V"] * a_DPM_res) / ((variables["E_vac"] / 10.0) * variables["c"]);
}

// Compute Oscillatory Resonance Term: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
double CrabResonanceUQFFModule::computeOscResTerm(double t) {
    double cos_term = 2 * variables["A"] * std::cos(variables["k"] * variables["x"]) * std::cos(variables["omega_osc"] * t);
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega_osc"] * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);
    return cos_term + exp_factor * real_exp;
}

// Compute Expansion Resonance Term: a_exp_res = (f_exp * E_vac * a_DPM_res) / ((E_vac / 10) * c)
double CrabResonanceUQFFModule::computeExpResTerm() {
    double a_DPM_res = computeDPMResTerm();
    return (variables["f_exp"] * variables["E_vac"] * a_DPM_res) / ((variables["E_vac"] / 10.0) * variables["c"]);
}

// Compute SC Resonance Integrated: (1 - B / B_crit) * f_sc
double CrabResonanceUQFFModule::computeSCResIntegrated(double B) {
    return (1.0 - (B / variables["B_crit"])) * variables["f_sc"];
}

// Full g_UQFF Resonance: Sum resonance terms * SC * (1 + f_TRZ)
double CrabResonanceUQFFModule::computeG(double t, double B) {
    variables["t"] = t;
    double a_DPM_res = computeDPMResTerm();
    double a_THz_res = computeTHzResTerm();
    double a_aether_res = computeAetherResTerm();
    double a_u_g4i_res = computeU_g4iResTerm();
    double a_quantum_res = computeQuantumResTerm();
    double a_fluid_res = computeFluidResTerm();
    double a_osc_res = computeOscResTerm(t);
    double a_exp_res = computeExpResTerm();
    double sc_int = computeSCResIntegrated(B);
    double tr_factor = 1.0 + variables["f_TRZ"];
    double res_sum = a_DPM_res + a_THz_res + a_aether_res + a_u_g4i_res + a_quantum_res + a_fluid_res + a_osc_res + a_exp_res;
    return res_sum * sc_int * tr_factor;
}

// Get equation text (descriptive)
std::string CrabResonanceUQFFModule::getEquationText() {
    return "g_Crab_Res(t, B) = [a_DPM_res + a_THz_res + a_aether_res + U_g4i_res + a_quantum_res + a_fluid_res + a_osc_res + a_exp_res] * SC_int * (1 + f_TRZ)\n"
           "Where:\n"
           "- a_DPM_res = (F_DPM * f_DPM * E_vac) / (c * V_sys(t)); F_DPM = I * A * (ω1 - ω2); V_sys(t) = 4/3 π r(t)^3, r(t)=r0 + v_exp t\n"
           "- a_THz_res = (f_THz * E_vac * v_exp * a_DPM_res) / (E_vac/10 * c)\n"
           "- a_aether_res = f_aether * 1e-8 * f_DPM * (1 + f_TRZ) * a_DPM_res\n"
           "- U_g4i_res = f_sc * Ug1_proxy * f_react * a_DPM_res / (E_vac * c)\n"
           "- a_quantum_res = (f_quantum * E_vac * a_DPM_res) / (E_vac/10 * c)\n"
           "- a_fluid_res = (f_fluid * E_vac * V * a_DPM_res) / (E_vac/10 * c)\n"
           "- a_osc_res = 2 A cos(k x) cos(ω t) + (2π / 13.8) A Re[exp(i (k x - ω t))]\n"
           "- a_exp_res = (f_exp * E_vac * a_DPM_res) / (E_vac/10 * c)\n"
           "- SC_int = (1 - B / B_crit) * f_sc\n"
           "Special Terms: UQFF resonance via plasmotic vacuum; Aether replaces dark energy; no SM terms; pulsar-driven f_osc.\n"
           "Solutions: At t=971 yr, B=1e-8 T, g ≈ 1e-40 m/s² (resonance micro-scale, wind proxy).\n"
           "Adaptations: Resonance focus for Crab wisps/shocks per Hubble/Chandra.";
}

// Print variables
void CrabResonanceUQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}