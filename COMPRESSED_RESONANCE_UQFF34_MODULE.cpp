// COMPRESSED_RESONANCE_UQFF34_MODULE.cpp
// UQFF 2.0 Multi-System Dual-Channel Compressed+Resonance Module — Session 92 (Systems 26-28,30-32,34)
// 34th C++ UQFF module — SECOND UQFF Dual-Channel (Compressed+Resonance) architecture
// Systems: 26=Universe Diameter, 27=Hydrogen Atom, 28=Hydrogen PToE Resonance,
//          30=Lagoon Nebula M8, 31=Spirals+SN Ia, 32=NGC 6302 Bug Nebula, 34=Orion Nebula M42
// g_CR34 = (a_DPM+a_THz+a_vac_diff+a_super + a_aether+a_u_g4i+a_osc+a_quantum+a_fluid+a_exp)
//          * (1-B/B_crit) * (1+f_TRZ)   [10-term dual-channel co-sum]
// PAPER_320: xi_span=1e35 — 7-system DPM force density spectral atlas (atomic to cosmic)
// PAPER_321: V_f_crossover=5.43e28 m^3/Hz — FIRST UQFF cross-channel dominance reversal threshold
// PAPER_322: a_THz Orion/Lagoon ratio=8.59 — FIRST UQFF intra-HII THz geometric amplification differential
// Stub bugs fixed: system_id not member->current_system_id; E_0 uninitialized->6.381e-36;
//   V_sys Orion stub=6.132e51->canonical 6.887e51; computeFullUQFF34 unimplemented->implemented
// Watermark: Copyright - Daniel T. Murphy, upgraded UQFF 2.0 Session 92 (March 18, 2026)

#ifndef COMPRESSED_RESONANCE_UQFF34_MODULE_H
#define COMPRESSED_RESONANCE_UQFF34_MODULE_H

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>

// ---------------------------------------------------------------------------
// WOLFRAM_TERM macros x4 — Session 92, PAPER_320-322
// ---------------------------------------------------------------------------
#define WOLFRAM_TERM_CR34_BASE \
    "g_CR34(t,B,sys)=(a_DPM+a_THz+a_vac_diff+a_super+a_aether+a_u_g4i+a_osc+a_quantum+a_fluid+a_exp)*(1-B/B_crit)*(1+f_TRZ);10-term dual-channel 7-system co-sum;34th C++ module [Session 92]"

#define WOLFRAM_TERM_CR34_DPM_SPECTRAL_ATLAS \
    "f_density=I*A_vort*d_omega/V_sys[N/m^3];Universe=1.500e-10;H_Atom=1.500e25;xi_span=1e35;Orion=9.12 N/m^3 HII balance;FIRST UQFF 35-order DPM force density spectral atlas 7 systems [PAPER_320]"

#define WOLFRAM_TERM_CR34_CHANNEL_CROSSOVER \
    "V_f_crossover=hbar/(E_0*f_vac_diff*E_vac*c)=5.43e28 m^3/Hz;V_sys/f_react>crossover:compressed dominant;H-Atom resonance-dominant 69 orders below;Universe compressed 44 orders above;FIRST UQFF cross-channel dominance reversal [PAPER_321]"

#define WOLFRAM_TERM_CR34_HII_THz_DIFFERENTIAL \
    "a_THz_34/a_THz_30=8.59;Orion/Lagoon same f_DPM=1e11/f_THz=1e11/v_exp=1e4;ratio=(A_vort_34*V_30)/(A_vort_30*V_34)=8.59;FIRST UQFF intra-HII THz geometric amplification differential [PAPER_322]"

// ---------------------------------------------------------------------------
// Per-system parameter struct
// ---------------------------------------------------------------------------
struct CR34SystemParams {
    double f_DPM, I_curr, A_vort, omega_diff, v_exp, V_sys;
    double f_THz, f_vac_diff, f_super;
    double f_aether, f_react, f_quantum, f_fluid, f_exp_freq;
    double omega_osc, k_wave, A_amp, V_fluid;
    const char* name;
};

// ===========================================================================
// Class declaration
// ===========================================================================
class CompressedResonanceUQFF34Module {
private:
    static constexpr double C_LIGHT      = 3.0e8;
    static constexpr double HBAR         = 1.0546e-34;
    static constexpr double E_VAC        = 7.09e-36;
    static constexpr double E_0          = 6.381e-36;   // reduced vacuum [PAPER_294 pattern]
    static constexpr double PI           = 3.141592653589793;
    static constexpr double T_COSMIC_GYR = 13.8;

    // Active system parameters
    double f_DPM, f_THz, f_vac_diff, f_super;
    double I_curr, A_vort, omega_diff, v_exp, V_sys;
    double f_aether, f_react, f_quantum, f_fluid, f_exp_freq;
    double omega_osc, k_wave, A_amp, V_fluid;
    double f_TRZ, B_crit, f_sc;
    int    current_system_id;

    // Cached quantities
    double F_DPM_cache, a_DPM_cache, Gamma_THz_cache;
    double A_sc_cache, a_vac_diff_cache, f_density_cache;

    // UQFF 2.0 runtime
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

    static const CR34SystemParams SYSTEMS[7];

    void updateCache();
    double computeDPMTerm();
    double computeTHzTerm();
    double computeVacDiffTerm();
    double computeSuperTerm();
    double computeAetherTerm();
    double computeU_g4iTerm();
    double computeOscTerm(double t);
    double computeQuantumTerm();
    double computeFluidTerm();
    double computeExpTerm();

public:
    explicit CompressedResonanceUQFF34Module(int system_id = 34);

    void setSystem(int system_id);

    // Main API (legacy + UQFF 2.0)
    double computeCompressedResTerm(double t, double B = 1e-5);
    double computeCompressed(int system_id);
    double computeResonance(int system_id, double t);
    double computeFullUQFF34(int system_id, double t, double B = 1e-5);
    double computeG(double t, double B = 1e-5) { return computeCompressedResTerm(t, B); }

    // UQFF 2.0
    void   setEnableLogging(bool e)                    { logging_enabled = e; }
    bool   getLoggingEnabled()               const     { return logging_enabled; }
    void   setDynamicParameter(const std::string& n, double v) { dynamic_params[n] = v; }
    double getDynamicParameter(const std::string& n)   const   {
        auto it = dynamic_params.find(n);
        return (it != dynamic_params.end()) ? it->second : std::numeric_limits<double>::quiet_NaN();
    }

    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    void exportState(const std::string& filename =
                     "CompressedResonanceUQFF34_state.txt") const;

    template <typename OtherModule>
    double cross_validate(OtherModule& other, double t, double B) {
        double g_self  = computeCompressedResTerm(t, B);
        double g_other = other.computeG(t, B);
        double ratio   = (g_other != 0.0) ? g_self / g_other : 0.0;
        if (logging_enabled)
            std::cout << "[CR34 XVAL] g_self=" << std::scientific << g_self
                      << " g_other=" << g_other << " ratio=" << ratio << "\n";
        return ratio;
    }

    std::string getEquationText(int system_id = -1) const;
    void printVariables() const;
};

// ===========================================================================
// System parameter table
// ===========================================================================
const CR34SystemParams CompressedResonanceUQFF34Module::SYSTEMS[7] = {
    // sys26: Universe Diameter
    {1e9,  1e24, 3.142e52, 2e-6,  1e8,    4.189e80,
     1e9,  0.143, 1.411e13, 1e3, 1e7,  1.445e-17, 1.269e-14, 1.373e-8,
     1e14, 1e17, 1e-9, 1e3,   "Universe Diameter"},
    // sys27: Hydrogen Atom
    {1e15, 1e18, 3.142e-21, 2e-3, 2.2e6, 4.189e-31,
     1e15, 0.143, 1.411e16, 1e4, 1e10, 1.445e-17, 1.269e-14, 1.373e-8,
     2.47e15, 1e11, 1e-10, 4.189e-31, "Hydrogen Atom"},
    // sys28: Hydrogen PToE Resonance
    {1e15, 1e18, 3.142e-21, 2e-3, 2.2e6, 4.189e-31,
     1e15, 0.143, 1.411e16, 1e4, 1e10, 1.445e-17, 1.269e-14, 1.373e-8,
     2.47e15, 1e11, 1e-10, 4.189e-31, "Hydrogen PToE Resonance"},
    // sys30: Lagoon Nebula M8
    {1e11, 1e20, 3.142e35, 2e-2, 1e4,   5.913e53,
     1e11, 0.143, 1.411e15, 1e2, 1e9,  1.445e-17, 1.269e-14, 1.373e-8,
     1e14, 1e15, 1e-9, 1e9,   "Lagoon Nebula M8"},
    // sys31: Spirals+SN Ia
    {1e10, 1e22, 3.142e41, 2e-1, 2e5,   1.543e64,
     1e10, 0.143, 1.411e14, 1e1, 1e8,  1.445e-17, 1.269e-14, 1.373e-8,
     1e13, 1e16, 1e-8, 1e12,  "Spirals+SN Ia"},
    // sys32: NGC 6302 Bug Nebula
    {1e12, 1e20, 3.142e32, 2e-3, 2.68e5, 1.458e48,
     1e12, 0.143, 1.411e16, 1e4, 1e10, 1.445e-17, 1.269e-14, 1.373e-8,
     1e15, 1e20, 1e-10, 1e3,  "NGC 6302 Bug Nebula"},
    // sys34: Orion Nebula M42 (V_sys canonical=(4/3)pi*(1.18e17)^3=6.887e51; stub 6.132e51 FIXED)
    {1e11, 1e20, 3.142e34, 2e-2, 1e4,   6.887e51,
     1e11, 0.143, 1.411e15, 1e2, 1e9,  1.445e-17, 1.269e-14, 1.373e-8,
     1e14, 1e15, 1e-9, 1e9,   "Orion Nebula M42/NGC 1976"},
};

static int cr34_sys_idx(int sid) {
    static const int IDS[7] = {26,27,28,30,31,32,34};
    for (int i = 0; i < 7; ++i) if (IDS[i] == sid) return i;
    return 6;
}

// ===========================================================================
// Implementation
// ===========================================================================

void CompressedResonanceUQFF34Module::updateCache() {
    F_DPM_cache      = I_curr * A_vort * omega_diff;
    a_DPM_cache      = (F_DPM_cache * f_DPM * E_VAC) / (C_LIGHT * V_sys);
    Gamma_THz_cache  = 10.0 * f_THz * v_exp / C_LIGHT;
    A_sc_cache       = (HBAR * f_super * f_DPM) / (E_VAC * C_LIGHT);   // [PAPER_295 pattern]
    a_vac_diff_cache = (E_0 * f_vac_diff * V_sys * a_DPM_cache) / HBAR; // [PAPER_294 pattern]
    f_density_cache  = (V_sys > 0.0) ? F_DPM_cache / V_sys : 0.0;       // N/m^3 [PAPER_320]
}

void CompressedResonanceUQFF34Module::setSystem(int sid) {
    int idx = cr34_sys_idx(sid);
    const CR34SystemParams& p = SYSTEMS[idx];
    current_system_id = sid;
    f_DPM      = p.f_DPM;      f_THz      = p.f_THz;    f_vac_diff = p.f_vac_diff;
    f_super    = p.f_super;    I_curr     = p.I_curr;   A_vort     = p.A_vort;
    omega_diff = p.omega_diff; v_exp      = p.v_exp;    V_sys      = p.V_sys;
    f_aether   = p.f_aether;   f_react    = p.f_react;  f_quantum  = p.f_quantum;
    f_fluid    = p.f_fluid;    f_exp_freq = p.f_exp_freq;
    omega_osc  = p.omega_osc;  k_wave     = p.k_wave;   A_amp      = p.A_amp;
    V_fluid    = p.V_fluid;
    updateCache();
}

CompressedResonanceUQFF34Module::CompressedResonanceUQFF34Module(int system_id)
    : f_TRZ(0.1), B_crit(1.0e11), f_sc(1.0),
      F_DPM_cache(0), a_DPM_cache(0), Gamma_THz_cache(0),
      A_sc_cache(0), a_vac_diff_cache(0), f_density_cache(0),
      current_system_id(system_id), logging_enabled(false)
{
    setSystem(system_id);
}

// Compressed channel
double CompressedResonanceUQFF34Module::computeDPMTerm()     { return a_DPM_cache; }
double CompressedResonanceUQFF34Module::computeTHzTerm()     { return Gamma_THz_cache * a_DPM_cache; }
double CompressedResonanceUQFF34Module::computeVacDiffTerm() { return a_vac_diff_cache; }
double CompressedResonanceUQFF34Module::computeSuperTerm()   { return A_sc_cache * a_DPM_cache; }

// Resonance channel
double CompressedResonanceUQFF34Module::computeAetherTerm() {
    return f_aether * 1.0e-8 * f_DPM * (1.0 + f_TRZ) * a_DPM_cache;
}
double CompressedResonanceUQFF34Module::computeU_g4iTerm() {
    return f_sc * f_react * a_DPM_cache / (E_VAC * C_LIGHT);
}
double CompressedResonanceUQFF34Module::computeOscTerm(double t) {
    double a_stand = 2.0 * A_amp * std::cos(k_wave * 0.0) * std::cos(omega_osc * t);
    std::complex<double> phase(0.0, k_wave * 0.0 - omega_osc * t);
    double a_trav  = (2.0 * PI / T_COSMIC_GYR) * A_amp * std::exp(phase).real();
    return a_stand + a_trav;
}
double CompressedResonanceUQFF34Module::computeQuantumTerm() {
    return (f_quantum * E_VAC * a_DPM_cache) / ((E_VAC / 10.0) * C_LIGHT);
}
double CompressedResonanceUQFF34Module::computeFluidTerm() {
    return (f_fluid * E_VAC * V_fluid * a_DPM_cache) / ((E_VAC / 10.0) * C_LIGHT);
}
double CompressedResonanceUQFF34Module::computeExpTerm() {
    return (f_exp_freq * E_VAC * a_DPM_cache) / ((E_VAC / 10.0) * C_LIGHT);
}

// Main computation
double CompressedResonanceUQFF34Module::computeCompressedResTerm(double t, double B) {
    const double a_DPM      = computeDPMTerm();
    const double a_THz      = computeTHzTerm();
    const double a_vac_diff = computeVacDiffTerm();
    const double a_super    = computeSuperTerm();
    const double a_comp     = a_DPM + a_THz + a_vac_diff + a_super;

    const double a_aether  = computeAetherTerm();
    const double a_u_g4i   = computeU_g4iTerm();
    const double a_osc     = computeOscTerm(t);
    const double a_quantum = computeQuantumTerm();
    const double a_fluid   = computeFluidTerm();
    const double a_exp     = computeExpTerm();
    const double a_res     = a_aether + a_u_g4i + a_osc + a_quantum + a_fluid + a_exp;

    const double SCm  = 1.0 - B / B_crit;
    const double g_CR = (a_comp + a_res) * SCm * (1.0 + f_TRZ);

    // [PAPER_321] channel dominance ratio — crossover at V_sys/f_react = 5.43e28
    const double R_CR = (a_res != 0.0) ? a_comp / a_res : 0.0;

    if (logging_enabled) {
        std::cout << std::scientific << std::setprecision(6)
                  << "[CR34][sys=" << current_system_id << "] "
                  << SYSTEMS[cr34_sys_idx(current_system_id)].name << "\n"
                  << "  [COMP]  a_DPM="      << a_DPM
                  << "  a_THz="      << a_THz
                  << "  a_vac_diff=" << a_vac_diff << "  a_super=" << a_super
                  << "  => a_comp=" << a_comp << "\n"
                  << "  [RES]   a_aether="   << a_aether
                  << "  a_u_g4i=" << a_u_g4i << "  a_osc=" << a_osc
                  << "  a_quantum=" << a_quantum << "  a_fluid=" << a_fluid
                  << "  a_exp=" << a_exp   << "  => a_res=" << a_res << "\n"
                  << "  [P320]  f_density="  << f_density_cache << " N/m^3\n"
                  << "  [P321]  R_CR(comp/res)=" << R_CR
                  << "  V_f=" << V_sys/f_react << " (crossover=5.43e28)\n"
                  << "  [RESULT] g_CR34=" << g_CR << " m/s^2\n";
    }
    return g_CR;
}

double CompressedResonanceUQFF34Module::computeCompressed(int sid) {
    setSystem(sid);
    return computeDPMTerm() + computeTHzTerm() + computeVacDiffTerm() + computeSuperTerm();
}
double CompressedResonanceUQFF34Module::computeResonance(int sid, double t) {
    setSystem(sid);
    return computeAetherTerm() + computeU_g4iTerm() + computeOscTerm(t)
         + computeQuantumTerm() + computeFluidTerm() + computeExpTerm();
}
double CompressedResonanceUQFF34Module::computeFullUQFF34(int sid, double t, double B) {
    setSystem(sid);
    return computeCompressedResTerm(t, B);
}

// exportState
void CompressedResonanceUQFF34Module::exportState(const std::string& filename) const {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) { std::cerr << "[CR34] Cannot open: " << filename << "\n"; return; }
    ofs << std::scientific << std::setprecision(10);
    ofs << "# CompressedResonanceUQFF34Module — UQFF 2.0 State (Session 92)\n";
    ofs << "# Papers: PAPER_320 PAPER_321 PAPER_322\n";
    ofs << "# 34th C++ UQFF module — 2nd dual-channel, 7-system atomic-to-cosmic\n";
    ofs << "current_system_id = " << current_system_id
        << "  # " << SYSTEMS[cr34_sys_idx(current_system_id)].name << "\n";
    ofs << "C_LIGHT     = " << C_LIGHT       << "\n";
    ofs << "HBAR        = " << HBAR          << "\n";
    ofs << "E_VAC       = " << E_VAC         << "\n";
    ofs << "E_0         = " << E_0           << "  # 6.381e-36 J/m^3 reduced vacuum\n";
    ofs << "f_DPM       = " << f_DPM         << "\n";
    ofs << "f_THz       = " << f_THz         << "\n";
    ofs << "f_vac_diff  = " << f_vac_diff    << "  # 0.143 Hz T_vac=6.993s\n";
    ofs << "f_super     = " << f_super       << "\n";
    ofs << "I_curr      = " << I_curr        << "\n";
    ofs << "A_vort      = " << A_vort        << "\n";
    ofs << "omega_diff  = " << omega_diff    << "\n";
    ofs << "v_exp       = " << v_exp         << "\n";
    ofs << "V_sys       = " << V_sys         << "\n";
    ofs << "f_aether    = " << f_aether      << "\n";
    ofs << "f_react     = " << f_react       << "\n";
    ofs << "f_quantum   = " << f_quantum     << "\n";
    ofs << "f_fluid     = " << f_fluid       << "\n";
    ofs << "f_exp_freq  = " << f_exp_freq    << "\n";
    ofs << "omega_osc   = " << omega_osc     << "\n";
    ofs << "k_wave      = " << k_wave        << "\n";
    ofs << "A_amp       = " << A_amp         << "\n";
    ofs << "V_fluid     = " << V_fluid       << "\n";
    ofs << "f_TRZ       = " << f_TRZ         << "\n";
    ofs << "B_crit      = " << B_crit        << "\n";
    ofs << "f_sc        = " << f_sc          << "\n";
    ofs << "# --- Cached ---\n";
    ofs << "F_DPM       = " << F_DPM_cache         << "\n";
    ofs << "a_DPM       = " << a_DPM_cache         << "\n";
    ofs << "Gamma_THz   = " << Gamma_THz_cache     << "\n";
    ofs << "A_sc        = " << A_sc_cache          << "  # [PAPER_295 pattern]\n";
    ofs << "a_vac_diff  = " << a_vac_diff_cache    << "  # [PAPER_294 pattern]\n";
    ofs << "f_density   = " << f_density_cache     << "  # N/m^3 [PAPER_320]\n";
    ofs << "# PAPER_320: xi_span=1e35 H_Atom/Universe DPM force density spectral atlas\n";
    ofs << "# PAPER_321: V_f_crossover=5.43e28 m^3/Hz compressed<->resonance dominance\n";
    ofs << "# PAPER_322: a_THz_Orion/a_THz_Lagoon=8.59 intra-HII THz geometric differential\n";
    for (const auto& dp : dynamic_params) ofs << "dyn:" << dp.first << " = " << dp.second << "\n";
    ofs.close();
}

// Legacy compat
void CompressedResonanceUQFF34Module::updateVariable(const std::string& name, double value) {
    dynamic_params[name] = value;
}
void CompressedResonanceUQFF34Module::addToVariable(const std::string& name, double delta) {
    dynamic_params[name] = getDynamicParameter(name) + delta;
}
void CompressedResonanceUQFF34Module::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

std::string CompressedResonanceUQFF34Module::getEquationText(int sid) const {
    if (sid < 0) sid = current_system_id;
    std::string sys_name(SYSTEMS[cr34_sys_idx(sid)].name);
    return
        "=== WOLFRAM_TERM_CR34_BASE ===\n" WOLFRAM_TERM_CR34_BASE "\n\n"
        "=== WOLFRAM_TERM_CR34_DPM_SPECTRAL_ATLAS ===\n" WOLFRAM_TERM_CR34_DPM_SPECTRAL_ATLAS "\n\n"
        "=== WOLFRAM_TERM_CR34_CHANNEL_CROSSOVER ===\n" WOLFRAM_TERM_CR34_CHANNEL_CROSSOVER "\n\n"
        "=== WOLFRAM_TERM_CR34_HII_THz_DIFFERENTIAL ===\n" WOLFRAM_TERM_CR34_HII_THz_DIFFERENTIAL "\n\n"
        "Active system: " + sys_name + " (id=" + std::to_string(sid) + ")\n"
        "  Compressed: a_DPM + a_THz + a_vac_diff + a_super\n"
        "  Resonance:  a_aether + a_u_g4i + a_osc + a_quantum + a_fluid + a_exp\n"
        "  g_CR34 = (a_comp + a_res) * (1-B/B_crit) * (1+f_TRZ)\n";
}

void CompressedResonanceUQFF34Module::printVariables() const {
    std::cout << std::scientific << std::setprecision(6)
              << "CompressedResonanceUQFF34Module — UQFF 2.0 (Session 92, 34th C++ module)\n"
              << "  Active: sys" << current_system_id
              << " (" << SYSTEMS[cr34_sys_idx(current_system_id)].name << ")\n"
              << "  [Comp]   f_DPM=" << f_DPM << " f_THz=" << f_THz
              << " f_vac_diff=" << f_vac_diff << " f_super=" << f_super << "\n"
              << "  [Vortex] I=" << I_curr << " A_vort=" << A_vort
              << " omega_diff=" << omega_diff << " V_sys=" << V_sys << "\n"
              << "  [Cache]  a_DPM=" << a_DPM_cache
              << " Gamma_THz=" << Gamma_THz_cache
              << " f_density=" << f_density_cache << " N/m^3 [P320]\n"
              << "  [Cache]  A_sc=" << A_sc_cache
              << " a_vac_diff=" << a_vac_diff_cache << "\n"
              << "  [SC/TRZ] B_crit=" << B_crit << " f_sc=" << f_sc
              << " f_TRZ=" << f_TRZ << "\n"
              << "  [P321]   V_f=" << V_sys/f_react
              << " crossover=5.43e28 dominant="
              << (V_sys/f_react > 5.43e28 ? "compressed" : "resonance") << "\n";
    if (!dynamic_params.empty())
        for (const auto& dp : dynamic_params)
            std::cout << "  [dyn] " << dp.first << "=" << dp.second << "\n";
}

#endif // COMPRESSED_RESONANCE_UQFF34_MODULE_H

// Set system-specific variables (system_id: 26=Universe, 27=Hydrogen, 28=PToE H, 30=Lagoon, 31=Spirals SN, 32=NGC6302, 34=Orion)
void CompressedResonanceUQFF34Module::setSystemVariables(int system_id) {
    switch (system_id) {
        case 26:  // Universe Diameter
            variables["f_DPM"] = 1e9; variables["I"] = 1e24; variables["A_vort"] = 3.142e52; variables["omega_1"] = 1e-6; variables["omega_2"] = -1e-6;
            variables["v_exp"] = 1e8; variables["V_sys"] = 4.189e80; variables["f_THz"] = 1e9; variables["f_vac_diff"] = 0.143; variables["f_super"] = 1.411e13;
            variables["f_aether"] = 1e3; variables["f_react"] = 1e7; variables["f_quantum"] = 1.445e-17; variables["f_fluid"] = 1.269e-14; variables["f_exp"] = 1.373e-8;
            variables["f_osc"] = 4.57e11; variables["k"] = 1e17; variables["omega_osc"] = 1e14; variables["x"] = 0.0; variables["A"] = 1e-9;
            variables["rho_fluid"] = 8.6e-27; variables["V"] = 1e3; variables["delta_rho"] = 0.1 * variables["rho_fluid"]; variables["rho"] = variables["rho_fluid"];
            variables["Delta_x"] = 1e-10; variables["Delta_p"] = variables["hbar"] / variables["Delta_x"]; variables["integral_psi"] = 1.0;
            break;
        case 27:  // Hydrogen Atom
            variables["f_DPM"] = 1e15; variables["I"] = 1e18; variables["A_vort"] = 3.142e-21; variables["omega_1"] = 1e-3; variables["omega_2"] = -1e-3;
            variables["v_exp"] = 2.2e6; variables["V_sys"] = 4.189e-31; variables["f_THz"] = 1e15; variables["f_vac_diff"] = 0.143; variables["f_super"] = 1.411e16;
            variables["f_aether"] = 1e4; variables["f_react"] = 1e10; variables["f_quantum"] = 1.445e-17; variables["f_fluid"] = 1.269e-14; variables["f_exp"] = 1.373e-8;
            variables["f_osc"] = 2.47e15; variables["k"] = 1e11; variables["omega_osc"] = 2.47e15; variables["x"] = 0.0; variables["A"] = 1e-10;
            variables["rho_fluid"] = 1e-25; variables["V"] = 4.189e-31; variables["delta_rho"] = 0.1 * variables["rho_fluid"]; variables["rho"] = variables["rho_fluid"];
            variables["Delta_x"] = 5.29e-11; variables["Delta_p"] = variables["hbar"] / variables["Delta_x"]; variables["integral_psi"] = 1.0;
            break;
        case 28:  // Hydrogen PToE Resonance
            variables["f_DPM"] = 1e15; variables["I"] = 1e18; variables["A_vort"] = 3.142e-21; variables["omega_1"] = 1e-3; variables["omega_2"] = -1e-3;
            variables["v_exp"] = 2.2e6; variables["V_sys"] = 4.189e-31; variables["f_THz"] = 1e15; variables["f_vac_diff"] = 0.143; variables["f_super"] = 1.411e16;
            variables["f_aether"] = 1e4; variables["f_react"] = 1e10; variables["f_quantum"] = 1.445e-17; variables["f_fluid"] = 1.269e-14; variables["f_exp"] = 1.373e-8;
            variables["f_osc"] = 2.47e15; variables["k"] = 1e11; variables["omega_osc"] = 2.47e15; variables["x"] = 0.0; variables["A"] = 1e-10;
            variables["rho_fluid"] = 1e-25; variables["V"] = 4.189e-31; variables["delta_rho"] = 0.1 * variables["rho_fluid"]; variables["rho"] = variables["rho_fluid"];
            variables["Delta_x"] = 5.29e-11; variables["Delta_p"] = variables["hbar"] / variables["Delta_x"]; variables["integral_psi"] = 1.0;
            break;
        case 30:  // Lagoon Nebula
            variables["f_DPM"] = 1e11; variables["I"] = 1e20; variables["A_vort"] = 3.142e35; variables["omega_1"] = 1e-2; variables["omega_2"] = -1e-2;
            variables["v_exp"] = 1e4; variables["V_sys"] = 5.913e53; variables["f_THz"] = 1e11; variables["f_vac_diff"] = 0.143; variables["f_super"] = 1.411e15;
            variables["f_aether"] = 1e2; variables["f_react"] = 1e9; variables["f_quantum"] = 1.445e-17; variables["f_fluid"] = 1.269e-14; variables["f_exp"] = 1.373e-8;
            variables["f_osc"] = 4.57e13; variables["k"] = 1e15; variables["omega_osc"] = 1e14; variables["x"] = 0.0; variables["A"] = 1e-9;
            variables["rho_fluid"] = 1e-20; variables["V"] = 1e9; variables["delta_rho"] = 0.1 * variables["rho_fluid"]; variables["rho"] = variables["rho_fluid"];
            variables["Delta_x"] = 1e-10; variables["Delta_p"] = variables["hbar"] / variables["Delta_x"]; variables["integral_psi"] = 1.0;
            break;
        case 31:  // Spirals and Supernovae
            variables["f_DPM"] = 1e10; variables["I"] = 1e22; variables["A_vort"] = 3.142e41; variables["omega_1"] = 1e-1; variables["omega_2"] = -1e-1;
            variables["v_exp"] = 2e5; variables["V_sys"] = 1.543e64; variables["f_THz"] = 1e10; variables["f_vac_diff"] = 0.143; variables["f_super"] = 1.411e14;
            variables["f_aether"] = 1e1; variables["f_react"] = 1e8; variables["f_quantum"] = 1.445e-17; variables["f_fluid"] = 1.269e-14; variables["f_exp"] = 1.373e-8;
            variables["f_osc"] = 4.57e12; variables["k"] = 1e16; variables["omega_osc"] = 1e13; variables["x"] = 0.0; variables["A"] = 1e-8;
            variables["rho_fluid"] = 1e-21; variables["V"] = 1e12; variables["delta_rho"] = 0.1 * variables["rho_fluid"]; variables["rho"] = variables["rho_fluid"];
            variables["Delta_x"] = 1e-10; variables["Delta_p"] = variables["hbar"] / variables["Delta_x"]; variables["integral_psi"] = 1.0;
            break;
        case 32:  // NGC 6302
            variables["f_DPM"] = 1e12; variables["I"] = 1e20; variables["A_vort"] = 3.142e32; variables["omega_1"] = 1e-3; variables["omega_2"] = -1e-3;
            variables["v_exp"] = 2.68e5; variables["V_sys"] = 1.458e48; variables["f_THz"] = 1e12; variables["f_vac_diff"] = 0.143; variables["f_super"] = 1.411e16;
            variables["f_aether"] = 1e4; variables["f_react"] = 1e10; variables["f_quantum"] = 1.445e-17; variables["f_fluid"] = 1.269e-14; variables["f_exp"] = 1.373e-8;
            variables["f_osc"] = 4.57e14; variables["k"] = 1e20; variables["omega_osc"] = 1e15; variables["x"] = 0.0; variables["A"] = 1e-10;
            variables["rho_fluid"] = 1e-21; variables["V"] = 1e3; variables["delta_rho"] = 0.1 * variables["rho_fluid"]; variables["rho"] = variables["rho_fluid"];
            variables["Delta_x"] = 1e-10; variables["Delta_p"] = variables["hbar"] / variables["Delta_x"]; variables["integral_psi"] = 1.0;
            break;
        case 34:  // Orion Nebula
            variables["f_DPM"] = 1e11; variables["I"] = 1e20; variables["A_vort"] = 3.142e34; variables["omega_1"] = 1e-2; variables["omega_2"] = -1e-2;
            variables["v_exp"] = 1e4; variables["V_sys"] = 6.132e51; variables["f_THz"] = 1e11; variables["f_vac_diff"] = 0.143; variables["f_super"] = 1.411e15;
            variables["f_aether"] = 1e2; variables["f_react"] = 1e9; variables["f_quantum"] = 1.445e-17; variables["f_fluid"] = 1.269e-14; variables["f_exp"] = 1.373e-8;
            variables["f_osc"] = 4.57e13; variables["k"] = 1e15; variables["omega_osc"] = 1e14; variables["x"] = 0.0; variables["A"] = 1e-9;
            variables["rho_fluid"] = 1e-20; variables["V"] = 1e9; variables["delta_rho"] = 0.1 * variables["rho_fluid"]; variables["rho"] = variables["rho_fluid"];
            variables["Delta_x"] = 1e-10; variables["Delta_p"] = variables["hbar"] / variables["Delta_x"]; variables["integral_psi"] = 1.0;
            break;
        default:
            std::cerr << "Unknown system_id: " << system_id << std::endl;
            break;
    }
}

// Update variable (set to new value)
void CompressedResonanceUQFF34Module::updateVariable(const std::string& name, double value) {
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
void CompressedResonanceUQFF34Module::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void CompressedResonanceUQFF34Module::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute Compressed Term: Sum streamlined DPM + THz + vac_diff + super
double CompressedResonanceUQFF34Module::computeCompressedTerm() {
    double F_DPM = variables["I"] * variables["A_vort"] * (variables["omega_1"] - variables["omega_2"]);
    double a_DPM = (F_DPM * variables["f_DPM"] * variables["E_vac"]) / (variables["c"] * variables["V_sys"]);
    double a_THz = (variables["f_THz"] * variables["E_vac"] * variables["v_exp"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
    double a_vac_diff = (variables["E_0"] * variables["f_vac_diff"] * variables["V_sys"] * a_DPM) / variables["hbar"];
    double a_super = (variables["hbar"] * variables["f_super"] * variables["f_DPM"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
    return a_DPM + a_THz + a_vac_diff + a_super;
}

// Compute Resonance Term: Sum aether + U_g4i + osc + quantum + fluid + exp
double CompressedResonanceUQFF34Module::computeResonanceTerm(double t) {
    double a_DPM = (variables["I"] * variables["A_vort"] * (variables["omega_1"] - variables["omega_2"]) * variables["f_DPM"] * variables["E_vac"]) / (variables["c"] * variables["V_sys"]);
    double a_aether = variables["f_aether"] * 1e-8 * variables["f_DPM"] * (1 + variables["f_TRZ"]) * a_DPM;
    double Ug1_proxy = 1.0;
    double a_u_g4i = variables["f_sc"] * Ug1_proxy * variables["f_react"] * a_DPM / (variables["E_vac"] * variables["c"]);
    double cos_term = 2 * variables["A"] * std::cos(variables["k"] * variables["x"]) * std::cos(variables["omega_osc"] * t);
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega_osc"] * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);
    double a_osc = cos_term + exp_factor * real_exp;
    double a_quantum = (variables["f_quantum"] * variables["E_vac"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
    double a_fluid = (variables["f_fluid"] * variables["E_vac"] * variables["V"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
    double a_exp = (variables["f_exp"] * variables["E_vac"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
    return a_aether + a_u_g4i + a_osc + a_quantum + a_fluid + a_exp;
}

// Compute SC Integrated: (1 - B / B_crit) * f_sc
double CompressedResonanceUQFF34Module::computeSCIntegrated(double B) {
    return (1.0 - (B / variables["B_crit"])) * variables["f_sc"];
}

// Full Compressed + Resonance with SC: (compressed + resonance) * SC * (1 + f_TRZ)
double CompressedResonanceUQFF34Module::computeCompressedResTerm(double t, double B) {
    setSystemVariables(system_id);  // Wait, this is in the method, but system_id not passed. Wait, for this code, assume it's set externally or add parameter.
    // Note: In usage, set system first.
    double comp = computeCompressedTerm();
    double res = computeResonanceTerm(t);
    double sc_int = computeSCIntegrated(B);
    double tr_factor = 1.0 + variables["f_TRZ"];
    return (comp + res) * sc_int * tr_factor;
}

// Get equation text (descriptive)
std::string CompressedResonanceUQFF34Module::getEquationText(int system_id) {
    std::string sys_name;
    switch (system_id) {
        case 26: sys_name = "Universe Diameter"; break;
        case 27: sys_name = "Hydrogen Atom"; break;
        case 28: sys_name = "Hydrogen PToE Resonance"; break;
        case 30: sys_name = "Lagoon Nebula"; break;
        case 31: sys_name = "Spirals and Supernovae"; break;
        case 32: sys_name = "NGC 6302"; break;
        case 34: sys_name = "Orion Nebula"; break;
        default: sys_name = "Unknown"; break;
    }
    return "Compressed Terms: a_comp = a_DPM + a_THz + a_vac_diff + a_super (scaled for " + sys_name + ")\n"
           "Resonance Terms: a_res = a_aether + U_g4i + a_osc + a_quantum + a_fluid + a_exp\n"
           "Full: g_comp_res = (a_comp + a_res) * SC_int * (1 + f_TRZ)\n"
           "Where SC_int = (1 - B / B_crit) * f_sc\n"
           "Special Terms: UQFF compressed/resonance via plasmotic vacuum; no SM; for system " + std::to_string(system_id) + " (" + sys_name + ").\n"
           "Solutions: See doc for system-specific g ~1e-33 to 1e35 m/s² (micro to macro scale).\n"
           "Adaptations: Frequencies scaled per system (e.g., f_DPM=1e9 for Universe, 1e15 for Hydrogen).";
}

// Print variables
void CompressedResonanceUQFF34Module::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}