// COMPRESSED_RESONANCE_UQFF34b_MODULE.cpp
// UQFF 2.0 Multi-System Dual-Channel Compressed+Resonance Module — Session 93 (Systems 18-20,22-24)
// 35th C++ UQFF module — CR34b VARIANT: 11-term dual-channel, 6-system galaxy-to-planetary
// Systems: 18=Sombrero Galaxy M104, 19=Andromeda Galaxy M31, 20=Universe Diameter,
//          22=Saturn, 23=M16 Eagle Nebula, 24=Crab Nebula M1
// g_CR34b = (a_DPM+a_THz+a_vac_diff+a_super + a_aether_res+a_u_g4i+a_osc+a_quantum_freq+
//             a_aether_freq+a_fluid_rho+a_exp)
//           * (1-B/B_crit) * (1+f_TRZ)   [11-term dual-channel co-sum; +1 vs CR34]
// Variant of CR34: adds a_aether_freq (11th term) and rho-weighted fluid term
// PAPER_323: F_AETHER=1.576e-35 Hz — 11th UQFF term: vacuum aether frequency mode, 5.253e-43*a_DPM
// PAPER_324: Saturn FIRST planetary dual-channel — g_vac_diff=1.29e-2 m/s^2 dominant
// PAPER_325: rho-ISM fluid coupling — f_fluid*rho_ISM=1.269e-35 kg/m^3/Hz density-enhanced fluid term
// Stub bugs fixed: V_sys Orion 6.132e51->6.887e51 (N/A here); omega_diff=2e-3 for all systems;
//   computeG stub unimplemented->11-term UQFF 2.0; variables[] map->UQFF 2.0 named typed members
// Watermark: Copyright - Daniel T. Murphy, upgraded UQFF 2.0 Session 93 (March 18, 2026)

#ifndef COMPRESSED_RESONANCE_UQFF34b_MODULE_H
#define COMPRESSED_RESONANCE_UQFF34b_MODULE_H

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>

// ---------------------------------------------------------------------------
// WOLFRAM_TERM macros x4 — Session 93, PAPER_323-325
// ---------------------------------------------------------------------------
#define WOLFRAM_TERM_CR34b_BASE \
    "g_CR34b(t,B,sys)=(a_DPM+a_THz+a_vac_diff+a_super+a_aether_res+a_u_g4i+a_osc+a_quantum_freq+a_aether_freq+a_fluid_rho+a_exp)*(1-B/B_crit)*(1+f_TRZ);11-term dual-channel 6-system galaxy-to-planetary;35th C++ module [Session 93]"

#define WOLFRAM_TERM_CR34b_AETHER_FREQ \
    "a_aether_freq=F_AETHER*E_VAC_neb*a_DPM/(E_VAC_ISM*c)=5.253e-43*a_DPM;F_AETHER=1.576e-35 Hz;smallest UQFF coupling coefficient;cosmological aether frequency mode distinct from a_aether_res;FIRST UQFF vacuum aether freq term [PAPER_323]"

#define WOLFRAM_TERM_CR34b_SATURN \
    "Saturn sys22 V_sys=9.184e23 m^3;FIRST planetary dual-channel;g_vac_diff=1.29e-2 m/s^2 dominant;f_DPM=1e12 Hz microwave regime;fills planetary gap in xi_span between atomic(4.189e-31) and nebular(5.913e50) [PAPER_324]"

#define WOLFRAM_TERM_CR34b_FLUID_RHO \
    "a_fluid_rho=f_fluid*E_VAC_neb*V_fluid*rho_ISM*a_DPM/(E_VAC_ISM*c);f_fluid*rho_ISM=1.269e-35 kg/m^3/Hz ISM density coupling;CR34b enhancement over CR34 which omits rho;reduces to CR34 fluid term when rho_ISM=1 [PAPER_325]"

// ---------------------------------------------------------------------------
// Per-system parameter struct (19 doubles + name — adds rho_fluid over CR34)
// ---------------------------------------------------------------------------
struct CR34bSystemParams {
    double f_DPM, I_curr, A_vort, omega_diff, v_exp, V_sys;
    double f_THz, f_vac_diff, f_super;
    double f_aether, f_react, f_quantum, f_fluid, f_exp_freq;
    double omega_osc, k_wave, A_amp, V_fluid, rho_fluid;
    const char* name;
};

// ===========================================================================
// Class declaration
// ===========================================================================
class CompressedResonanceUQFF34bModule {
private:
    static constexpr double C_LIGHT      = 3.0e8;
    static constexpr double HBAR         = 1.0546e-34;
    static constexpr double E_VAC_NEB    = 7.09e-36;    // nebular vacuum energy density
    static constexpr double E_VAC_ISM    = 7.09e-37;    // ISM vacuum (E_VAC_NEB/10) [CR34b]
    static constexpr double E_0          = 6.381e-36;   // reduced vacuum [PAPER_294 pattern]
    static constexpr double F_AETHER     = 1.576e-35;   // vacuum aether frequency [PAPER_323]
    static constexpr double LAMBDA       = 1.1e-52;     // cosmological constant
    static constexpr double PI           = 3.141592653589793;
    static constexpr double T_COSMIC_GYR = 13.8;

    // Active system parameters
    double f_DPM, f_THz, f_vac_diff, f_super;
    double I_curr, A_vort, omega_diff, v_exp, V_sys;
    double f_aether, f_react, f_quantum, f_fluid, f_exp_freq;
    double omega_osc, k_wave, A_amp, V_fluid, rho_fluid;
    double f_TRZ, B_crit, f_sc;
    int    current_system_id;

    // Cached quantities
    double F_DPM_cache, a_DPM_cache, Gamma_THz_cache;
    double A_sc_cache, a_vac_diff_cache, a_aether_freq_cache, f_density_cache;

    // UQFF 2.0 runtime
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

    static const CR34bSystemParams SYSTEMS[6];

    void updateCache();
    double computeDPMTerm();
    double computeTHzTerm();
    double computeVacDiffTerm();
    double computeSuperTerm();
    double computeAetherResTerm();
    double computeU_g4iTerm();
    double computeOscTerm(double t);
    double computeQuantumFreqTerm();
    double computeAetherFreqTerm();
    double computeFluidRhoTerm();
    double computeExpTerm();

public:
    explicit CompressedResonanceUQFF34bModule(int system_id = 18);

    void setSystem(int system_id);

    // Main API (legacy + UQFF 2.0)
    double computeG(double t, double B = 1e-5);
    double computeCompressedResTerm(double t, double B = 1e-5) { return computeG(t, B); }
    double computeCompressed(int system_id);
    double computeResonance(int system_id, double t);
    double computeFullUQFF34b(int system_id, double t, double B = 1e-5);

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
                     "CompressedResonanceUQFF34b_state.txt") const;

    template <typename OtherModule>
    double cross_validate(OtherModule& other, double t, double B) {
        double g_self  = computeG(t, B);
        double g_other = other.computeG(t, B);
        double ratio   = (g_other != 0.0) ? g_self / g_other : 0.0;
        if (logging_enabled)
            std::cout << "[CR34b XVAL] g_self=" << std::scientific << g_self
                      << " g_other=" << g_other << " ratio=" << ratio << "\n";
        return ratio;
    }

    std::string getEquationText(int system_id = -1) const;
    void printVariables() const;
};

// ===========================================================================
// System parameter table — 6 systems, all omega_diff=2e-3 (stub common default)
// ===========================================================================
const CR34bSystemParams CompressedResonanceUQFF34bModule::SYSTEMS[6] = {
    // sys18: Sombrero Galaxy M104 (V_sys=4.434e62; A_vort=7.487e42 precise)
    {1e10, 1e22, 7.487e42, 2e-3, 2e5,    4.434e62,
     1e10, 0.143, 1.411e14, 1e1, 1e8,  1.445e-17, 1.269e-14, 1.373e-8,
     1e13, 1e16, 1e-8, 1e12, 1e-21,   "Sombrero Galaxy M104"},
    // sys19: Andromeda Galaxy M31 (V_sys=1.543e64; A_vort=3.142e42 pi-placeholder noted)
    {1e10, 1e23, 3.142e42, 2e-3, 3e5,    1.543e64,
     1e10, 0.143, 1.411e14, 1e1, 1e8,  1.445e-17, 1.269e-14, 1.373e-8,
     1e13, 1e16, 1e-8, 1e12, 1e-21,   "Andromeda Galaxy M31"},
    // sys20: Universe Diameter (V_sys=3.568e80; note CR34 sys26 uses 4.189e80; see PAPER_325 note)
    {1e9,  1e24, 3.142e52, 2e-3, 1e8,    3.568e80,
     1e9,  0.143, 1.411e13, 1e3, 1e7,  1.445e-17, 1.269e-14, 1.373e-8,
     1e14, 1e17, 1e-9, 1e3,  8.6e-27, "Universe Diameter"},
    // sys22: Saturn (V_sys=9.184e23; FIRST planetary dual-channel [PAPER_324]; rho_fluid=ISM proxy)
    {1e12, 1e19, 3.142e15, 2e-3, 5e3,    9.184e23,
     1e12, 0.143, 1.411e16, 1e4, 1e10, 1.445e-17, 1.269e-14, 1.373e-8,
     1e15, 1e20, 1e-10, 1e3, 1e-21,   "Saturn"},
    // sys23: M16 Eagle Nebula (V_sys=5.913e53; A_vort=3.142e35 pi-placeholder; same scale as Lagoon)
    {1e11, 1e20, 3.142e35, 2e-3, 1e4,    5.913e53,
     1e11, 0.143, 1.411e15, 1e2, 1e9,  1.445e-17, 1.269e-14, 1.373e-8,
     1e14, 1e15, 1e-9, 1e9,  1e-20,   "M16 Eagle Nebula"},
    // sys24: Crab Nebula M1 (V_sys=5.913e50; A_vort=3.142e32 pi-placeholder; v_exp=1.34e6 observed)
    {1e12, 1e21, 3.142e32, 2e-3, 1.34e6, 5.913e50,
     1e12, 0.143, 1.411e16, 1e4, 1e10, 1.445e-17, 1.269e-14, 1.373e-8,
     1e15, 1e20, 1e-10, 1e3, 1e-21,   "Crab Nebula M1"},
};

static int cr34b_sys_idx(int sid) {
    static const int IDS[6] = {18,19,20,22,23,24};
    for (int i = 0; i < 6; ++i) if (IDS[i] == sid) return i;
    return 0;
}

// ===========================================================================
// Implementation
// ===========================================================================

void CompressedResonanceUQFF34bModule::updateCache() {
    F_DPM_cache         = I_curr * A_vort * omega_diff;
    a_DPM_cache         = (F_DPM_cache * f_DPM * E_VAC_NEB) / (C_LIGHT * V_sys);
    Gamma_THz_cache     = 10.0 * f_THz * v_exp / C_LIGHT;
    A_sc_cache          = (HBAR * f_super * f_DPM) / (E_VAC_NEB * C_LIGHT); // [PAPER_295 pattern]
    a_vac_diff_cache    = (E_0 * f_vac_diff * V_sys * a_DPM_cache) / HBAR;  // [PAPER_294 pattern]
    a_aether_freq_cache = (F_AETHER * E_VAC_NEB * a_DPM_cache) / (E_VAC_ISM * C_LIGHT); // [PAPER_323]
    f_density_cache     = (V_sys > 0.0) ? F_DPM_cache / V_sys : 0.0;        // N/m^3 [PAPER_320 ext.]
}

void CompressedResonanceUQFF34bModule::setSystem(int sid) {
    int idx = cr34b_sys_idx(sid);
    const CR34bSystemParams& p = SYSTEMS[idx];
    current_system_id = sid;
    f_DPM      = p.f_DPM;      f_THz      = p.f_THz;    f_vac_diff = p.f_vac_diff;
    f_super    = p.f_super;    I_curr     = p.I_curr;   A_vort     = p.A_vort;
    omega_diff = p.omega_diff; v_exp      = p.v_exp;    V_sys      = p.V_sys;
    f_aether   = p.f_aether;   f_react    = p.f_react;  f_quantum  = p.f_quantum;
    f_fluid    = p.f_fluid;    f_exp_freq = p.f_exp_freq;
    omega_osc  = p.omega_osc;  k_wave     = p.k_wave;   A_amp      = p.A_amp;
    V_fluid    = p.V_fluid;    rho_fluid  = p.rho_fluid;
    updateCache();
}

CompressedResonanceUQFF34bModule::CompressedResonanceUQFF34bModule(int system_id)
    : f_TRZ(0.1), B_crit(1.0e11), f_sc(1.0),
      F_DPM_cache(0), a_DPM_cache(0), Gamma_THz_cache(0),
      A_sc_cache(0), a_vac_diff_cache(0), a_aether_freq_cache(0), f_density_cache(0),
      current_system_id(system_id), logging_enabled(false)
{
    setSystem(system_id);
}

// Compressed channel
double CompressedResonanceUQFF34bModule::computeDPMTerm()     { return a_DPM_cache; }
double CompressedResonanceUQFF34bModule::computeTHzTerm()     { return Gamma_THz_cache * a_DPM_cache; }
double CompressedResonanceUQFF34bModule::computeVacDiffTerm() { return a_vac_diff_cache; }
double CompressedResonanceUQFF34bModule::computeSuperTerm()   { return A_sc_cache * a_DPM_cache; }

// Resonance channel
double CompressedResonanceUQFF34bModule::computeAetherResTerm() {
    return f_aether * 1.0e-8 * f_DPM * (1.0 + f_TRZ) * a_DPM_cache;
}
double CompressedResonanceUQFF34bModule::computeU_g4iTerm() {
    return f_sc * f_react * a_DPM_cache / (E_VAC_NEB * C_LIGHT);
}
double CompressedResonanceUQFF34bModule::computeOscTerm(double t) {
    double a_stand = 2.0 * A_amp * std::cos(k_wave * 0.0) * std::cos(omega_osc * t);
    std::complex<double> phase(0.0, k_wave * 0.0 - omega_osc * t);
    double a_trav  = (2.0 * PI / T_COSMIC_GYR) * A_amp * std::exp(phase).real();
    return a_stand + a_trav;
}
double CompressedResonanceUQFF34bModule::computeQuantumFreqTerm() {
    return (f_quantum * E_VAC_NEB * a_DPM_cache) / (E_VAC_ISM * C_LIGHT);
}
double CompressedResonanceUQFF34bModule::computeAetherFreqTerm() {
    // 11th term — vacuum aether frequency mode [PAPER_323]
    return a_aether_freq_cache;
}
double CompressedResonanceUQFF34bModule::computeFluidRhoTerm() {
    // rho-weighted fluid term; CR34 uses volumetric only [PAPER_325]
    return (f_fluid * E_VAC_NEB * V_fluid * rho_fluid * a_DPM_cache) / (E_VAC_ISM * C_LIGHT);
}
double CompressedResonanceUQFF34bModule::computeExpTerm() {
    return (f_exp_freq * E_VAC_NEB * a_DPM_cache) / (E_VAC_ISM * C_LIGHT);
}

// Main computation — 11-term dual-channel
double CompressedResonanceUQFF34bModule::computeG(double t, double B) {
    const double a_DPM          = computeDPMTerm();
    const double a_THz          = computeTHzTerm();
    const double a_vac_diff     = computeVacDiffTerm();
    const double a_super        = computeSuperTerm();
    const double a_comp         = a_DPM + a_THz + a_vac_diff + a_super;

    const double a_aether_res   = computeAetherResTerm();
    const double a_u_g4i        = computeU_g4iTerm();
    const double a_osc          = computeOscTerm(t);
    const double a_quantum_freq = computeQuantumFreqTerm();
    const double a_aether_freq  = computeAetherFreqTerm();   // 11th term [PAPER_323]
    const double a_fluid_rho    = computeFluidRhoTerm();     // rho-weighted [PAPER_325]
    const double a_exp          = computeExpTerm();
    const double a_res          = a_aether_res + a_u_g4i + a_osc + a_quantum_freq
                                + a_aether_freq + a_fluid_rho + a_exp;

    const double SCm     = 1.0 - B / B_crit;
    const double g_CR34b = (a_comp + a_res) * SCm * (1.0 + f_TRZ);

    // [PAPER_324] Saturn planetary gap; [PAPER_325] rho-ISM coupling
    const double R_CR = (a_res != 0.0) ? a_comp / a_res : 0.0;

    if (logging_enabled) {
        std::cout << std::scientific << std::setprecision(6)
                  << "[CR34b][sys=" << current_system_id << "] "
                  << SYSTEMS[cr34b_sys_idx(current_system_id)].name << "\n"
                  << "  [COMP]  a_DPM="      << a_DPM
                  << "  a_THz="      << a_THz
                  << "  a_vac_diff=" << a_vac_diff  << "  a_super=" << a_super
                  << "  => a_comp=" << a_comp << "\n"
                  << "  [RES]   a_aether_res=" << a_aether_res
                  << "  a_u_g4i=" << a_u_g4i  << "  a_osc=" << a_osc
                  << "  a_q_freq=" << a_quantum_freq
                  << "  a_aeth_freq=" << a_aether_freq << " [P323]\n"
                  << "          a_fluid_rho=" << a_fluid_rho << " [P325]"
                  << "  a_exp=" << a_exp << "  => a_res=" << a_res << "\n"
                  << "  [P323]  a_aether_freq=" << a_aether_freq
                  << " coeff=" << (a_DPM_cache != 0.0 ? a_aether_freq/a_DPM_cache : 0.0) << "\n"
                  << "  [P324]  f_density=" << f_density_cache << " N/m^3"
                  << "  V_sys=" << V_sys << "\n"
                  << "  [P325]  rho_fluid=" << rho_fluid
                  << "  f_fluid*rho=" << f_fluid*rho_fluid << " kg/m^3/Hz\n"
                  << "  [CH]    R_CR(comp/res)=" << R_CR << "\n"
                  << "  [RESULT] g_CR34b=" << g_CR34b << " m/s^2\n";
    }
    return g_CR34b;
}

double CompressedResonanceUQFF34bModule::computeCompressed(int sid) {
    setSystem(sid);
    return computeDPMTerm() + computeTHzTerm() + computeVacDiffTerm() + computeSuperTerm();
}
double CompressedResonanceUQFF34bModule::computeResonance(int sid, double t) {
    setSystem(sid);
    return computeAetherResTerm() + computeU_g4iTerm() + computeOscTerm(t)
         + computeQuantumFreqTerm() + computeAetherFreqTerm()
         + computeFluidRhoTerm() + computeExpTerm();
}
double CompressedResonanceUQFF34bModule::computeFullUQFF34b(int sid, double t, double B) {
    setSystem(sid);
    return computeG(t, B);
}

// exportState
void CompressedResonanceUQFF34bModule::exportState(const std::string& filename) const {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) { std::cerr << "[CR34b] Cannot open: " << filename << "\n"; return; }
    ofs << std::scientific << std::setprecision(10);
    ofs << "# CompressedResonanceUQFF34bModule — UQFF 2.0 State (Session 93)\n";
    ofs << "# Papers: PAPER_323 PAPER_324 PAPER_325\n";
    ofs << "# 35th C++ UQFF module — CR34b variant: 11-term, 6-system galaxy-to-planetary\n";
    ofs << "# CR34b adds: a_aether_freq (11th term), rho-weighted fluid, E_VAC_NEB/ISM split\n";
    ofs << "current_system_id = " << current_system_id
        << "  # " << SYSTEMS[cr34b_sys_idx(current_system_id)].name << "\n";
    ofs << "C_LIGHT      = " << C_LIGHT       << "\n";
    ofs << "HBAR         = " << HBAR          << "\n";
    ofs << "E_VAC_NEB    = " << E_VAC_NEB     << "  # nebular vacuum energy density\n";
    ofs << "E_VAC_ISM    = " << E_VAC_ISM     << "  # E_VAC_NEB/10 ISM proxy\n";
    ofs << "E_0          = " << E_0           << "  # 6.381e-36 reduced vacuum\n";
    ofs << "F_AETHER     = " << F_AETHER      << "  # 1.576e-35 Hz vacuum aether freq [PAPER_323]\n";
    ofs << "LAMBDA       = " << LAMBDA        << "  # 1.1e-52 cosmological constant\n";
    ofs << "f_DPM        = " << f_DPM         << "\n";
    ofs << "f_THz        = " << f_THz         << "\n";
    ofs << "f_vac_diff   = " << f_vac_diff    << "  # 0.143 Hz T_vac=6.993s\n";
    ofs << "f_super      = " << f_super       << "\n";
    ofs << "I_curr       = " << I_curr        << "\n";
    ofs << "A_vort       = " << A_vort        << "\n";
    ofs << "omega_diff   = " << omega_diff    << "\n";
    ofs << "v_exp        = " << v_exp         << "\n";
    ofs << "V_sys        = " << V_sys         << "\n";
    ofs << "f_aether     = " << f_aether      << "\n";
    ofs << "f_react      = " << f_react       << "\n";
    ofs << "f_quantum    = " << f_quantum     << "\n";
    ofs << "f_fluid      = " << f_fluid       << "\n";
    ofs << "f_exp_freq   = " << f_exp_freq    << "\n";
    ofs << "omega_osc    = " << omega_osc     << "\n";
    ofs << "k_wave       = " << k_wave        << "\n";
    ofs << "A_amp        = " << A_amp         << "\n";
    ofs << "V_fluid      = " << V_fluid       << "\n";
    ofs << "rho_fluid    = " << rho_fluid     << "  # ISM density proxy kg/m^3 [PAPER_325]\n";
    ofs << "f_TRZ        = " << f_TRZ         << "\n";
    ofs << "B_crit       = " << B_crit        << "\n";
    ofs << "f_sc         = " << f_sc          << "\n";
    ofs << "# --- Cached ---\n";
    ofs << "F_DPM        = " << F_DPM_cache          << "\n";
    ofs << "a_DPM        = " << a_DPM_cache          << "\n";
    ofs << "Gamma_THz    = " << Gamma_THz_cache      << "\n";
    ofs << "A_sc         = " << A_sc_cache           << "  # [PAPER_295 pattern]\n";
    ofs << "a_vac_diff   = " << a_vac_diff_cache     << "  # [PAPER_294 pattern]\n";
    ofs << "a_aether_freq= " << a_aether_freq_cache  << "  # [PAPER_323] 11th term\n";
    ofs << "f_density    = " << f_density_cache      << "  # N/m^3 [PAPER_320 ext.]\n";
    ofs << "# PAPER_323: F_AETHER=1.576e-35 Hz; a_aether_freq=5.253e-43*a_DPM; 11th UQFF term\n";
    ofs << "# PAPER_324: Saturn V_sys=9.184e23 m^3; FIRST planetary dual-channel; g_vac_diff=1.29e-2\n";
    ofs << "# PAPER_325: a_fluid_rho=f_fluid*E_VAC_neb*V*rho*a_DPM/(E_VAC_ISM*c); f_fluid*rho=1.269e-35\n";
    for (const auto& dp : dynamic_params) ofs << "dyn:" << dp.first << " = " << dp.second << "\n";
    ofs.close();
}

// Legacy compat
void CompressedResonanceUQFF34bModule::updateVariable(const std::string& name, double value) {
    dynamic_params[name] = value;
}
void CompressedResonanceUQFF34bModule::addToVariable(const std::string& name, double delta) {
    dynamic_params[name] = getDynamicParameter(name) + delta;
}
void CompressedResonanceUQFF34bModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

std::string CompressedResonanceUQFF34bModule::getEquationText(int sid) const {
    if (sid < 0) sid = current_system_id;
    std::string sys_name(SYSTEMS[cr34b_sys_idx(sid)].name);
    return
        "=== WOLFRAM_TERM_CR34b_BASE ===\n" WOLFRAM_TERM_CR34b_BASE "\n\n"
        "=== WOLFRAM_TERM_CR34b_AETHER_FREQ ===\n" WOLFRAM_TERM_CR34b_AETHER_FREQ "\n\n"
        "=== WOLFRAM_TERM_CR34b_SATURN ===\n" WOLFRAM_TERM_CR34b_SATURN "\n\n"
        "=== WOLFRAM_TERM_CR34b_FLUID_RHO ===\n" WOLFRAM_TERM_CR34b_FLUID_RHO "\n\n"
        "Active system: " + sys_name + " (id=" + std::to_string(sid) + ")\n"
        "  Compressed: a_DPM + a_THz + a_vac_diff + a_super\n"
        "  Resonance:  a_aether_res + a_u_g4i + a_osc + a_quantum_freq"
              " + a_aether_freq + a_fluid_rho + a_exp\n"
        "  g_CR34b = (a_comp + a_res) * (1-B/B_crit) * (1+f_TRZ)\n";
}

void CompressedResonanceUQFF34bModule::printVariables() const {
    std::cout << std::scientific << std::setprecision(6)
              << "CompressedResonanceUQFF34bModule — UQFF 2.0 (Session 93, 35th C++ module)\n"
              << "  Active: sys" << current_system_id
              << " (" << SYSTEMS[cr34b_sys_idx(current_system_id)].name << ")\n"
              << "  [Comp]   f_DPM=" << f_DPM << " f_THz=" << f_THz
              << " f_vac_diff=" << f_vac_diff << " f_super=" << f_super << "\n"
              << "  [Vortex] I=" << I_curr << " A_vort=" << A_vort
              << " omega_diff=" << omega_diff << " V_sys=" << V_sys << "\n"
              << "  [Cache]  a_DPM=" << a_DPM_cache
              << " Gamma_THz=" << Gamma_THz_cache
              << " f_density=" << f_density_cache << " N/m^3\n"
              << "  [Cache]  A_sc=" << A_sc_cache
              << " a_vac_diff=" << a_vac_diff_cache << "\n"
              << "  [P323]   a_aether_freq=" << a_aether_freq_cache
              << " F_AETHER=" << F_AETHER << "\n"
              << "  [P324]   Saturn? " << (current_system_id==22 ? "YES — first planetary" : "no")
              << " V_sys=" << V_sys << "\n"
              << "  [P325]   rho_fluid=" << rho_fluid
              << " f_fluid*rho=" << f_fluid*rho_fluid << " kg/m^3/Hz\n"
              << "  [SC/TRZ] B_crit=" << B_crit << " f_sc=" << f_sc
              << " f_TRZ=" << f_TRZ << "\n";
    if (!dynamic_params.empty())
        for (const auto& dp : dynamic_params)
            std::cout << "  [dyn] " << dp.first << "=" << dp.second << "\n";
}

#endif // COMPRESSED_RESONANCE_UQFF34b_MODULE_H
