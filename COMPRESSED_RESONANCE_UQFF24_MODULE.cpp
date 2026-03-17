// COMPRESSED_RESONANCE_UQFF24_MODULE.cpp
// UQFF 2.0 Dual-Channel Compressed+Resonance Module — Session 83 (Systems 18-24)
// 25th C++ UQFF module — FIRST explicit Dual-Channel (Compressed + Resonance) architecture.
// Covers: Sombrero (M104), Saturn, M16 Eagle Nebula, Crab Nebula PWN, NGC 1792, HUDF, Andromeda +.
// 10-term co-sum: Compressed (DPM+THz+vac_diff+super) + Resonance (aether+U_g4i+osc+quantum+fluid+exp)
// g_CR = (a_comp + a_res) * SCm * (1 + f_TRZ)   [PAPER_293]
// PAPER_294: a_vac_diff=(E_0*f_vac_diff*V_sys*a_DPM)/hbar; f_vac_diff=0.143Hz T_vac=6.993s; FIRST UQFF hbar-denominator term
// PAPER_295: a_super=A_sc*a_DPM; A_sc prop f_DPM → a_super prop f_DPM^2; systems-18-24: A_sc=6.994e18 vs 1e12-class A_sc=6.994e21
// Watermark: Copyright - Daniel T. Murphy, upgraded UQFF 2.0 Session 83 (March 17, 2026)

#ifndef COMPRESSED_RESONANCE_UQFF24_MODULE_H
#define COMPRESSED_RESONANCE_UQFF24_MODULE_H

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>

// ---------------------------------------------------------------------------
// WOLFRAM_TERM macros — 4 unique physics anchors (Session 83, PAPER_293–295)
// ---------------------------------------------------------------------------
#define WOLFRAM_TERM_CR24_BASE \
    "g_CR(t,B)=(a_DPM+a_THz+a_vac_diff+a_super+a_aether+a_u_g4i+a_osc+a_quantum+a_fluid+a_exp)*(1-B/B_crit)*(1+f_TRZ);10-term dual-channel co-sum [PAPER_293]"

#define WOLFRAM_TERM_CR24_VAC_DIFF \
    "a_vac_diff=(E_0*f_vac_diff*V_sys*a_DPM)/hbar;f_vac_diff=0.143Hz;T_vac=1/0.143=6.993s;E_0/E_vac=0.9001(10pct deficit);V_sys/hbar=3.973e52;FIRST UQFF hbar-denom quantum-volume diffusion [PAPER_294]"

#define WOLFRAM_TERM_CR24_SUPER_COMP \
    "a_super=(hbar*f_super*f_DPM/(E_vac*c))*a_DPM=A_sc*a_DPM;A_sc prop f_DPM;a_super prop f_DPM^2;f_DPM=1e11:A_sc=6.994e18;f_DPM=1e12:A_sc=6.994e21(1000x);compressed-channel pre-osc Cooper seeding [PAPER_295]"

#define WOLFRAM_TERM_CR24_DUAL_CHANNEL \
    "R_CR=Sigma_comp/Sigma_res;Sigma_comp=a_DPM+a_THz+a_vac_diff+a_super;Sigma_res=a_aether+a_u_g4i+a_osc+a_quantum+a_fluid+a_exp;R_CR(sys18-24)=1.490e-17;res dominates 17 orders;FIRST UQFF explicit 4+6 dual-channel [PAPER_293]"

// ===========================================================================
// Class declaration
// ===========================================================================

class CompressedResonanceUQFF24Module {
private:
    // --- Universal physical constants ---
    static constexpr double C_LIGHT      = 3.0e8;           // m/s
    static constexpr double HBAR         = 1.0546e-34;      // J*s reduced Planck
    static constexpr double E_VAC        = 7.09e-36;        // J/m^3 plasmotic vacuum
    static constexpr double PI           = 3.141592653589793;
    static constexpr double T_COSMIC_GYR = 13.8;            // Gyr cosmic age normalisation

    // --- Compressed-channel parameters ---
    double f_DPM;       // Hz DPM mode (1e11 for galaxies/nebulae/planets — systems 18-24)
    double f_THz;       // Hz THz mode (matched to f_DPM at 1e11)
    double f_vac_diff;  // Hz vacuum differential beat (0.143 Hz → T_vac≈7s) [PAPER_294]
    double f_super;     // Hz Cooper pair / superconductive frequency (1.411e15) [PAPER_295]
    double I_curr;      // A circulating plasmotic current
    double A_vort;      // m^2 vortex cross-section area
    double omega_1;     // rad/s vortex differential +side
    double omega_2;     // rad/s vortex differential -side
    double v_exp;       // m/s system outflow / expansion velocity
    double E_0;         // J/m^3 reduced vacuum energy (E_0 < E_vac → 10% vacuum deficit) [PAPER_294]
    double V_sys;       // m^3 system characteristic volume

    // --- Resonance-channel parameters ---
    double f_aether;    // Hz aether mode coupling
    double f_react;     // Hz U_g4i reaction frequency
    double f_quantum;   // Hz quantum de Broglie / Compton mode
    double f_fluid;     // Hz Kelvin-Helmholtz / fluid turbulence mode
    double f_exp_freq;  // Hz free expansion mode (renamed; avoids std::exp ambiguity)
    double omega_osc;   // rad/s oscillatory standing+traveling wave frequency
    double k_wave;      // m^-1 wave vector (scaled for system aperture)
    double A_amp;       // m oscillation amplitude
    double V_fluid;     // m^3 fluid coupling volume (distinct from V_sys)
    double f_TRZ;       // dimensionless time-reversal zone factor

    // --- Superconductive / quench ---
    double B_crit;      // T critical magnetic field for Meissner quench
    double f_sc;        // dimensionless SC scaling factor

    // --- Cached computed quantities (set by updateCache) ---
    double F_DPM_cache;         // N DPM vortex force
    double a_DPM_cache;         // m/s^2 DPM base acceleration (shared seed for both channels)
    double A_sc_cache;          // dimensionless Cooper-super amplitude [PAPER_295]: A_sc prop f_DPM
    double Gamma_THz_cache;     // dimensionless THz amplification = 10*f_THz*v_exp/c
    double a_vac_diff_cache;    // m/s^2 vacuum differential term [PAPER_294]
    double E_0_over_E_vac;      // ratio 0.9001 vacuum deficit fraction [PAPER_294]

    // --- UQFF 2.0 runtime infrastructure ---
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

    // --- Private compute helpers ---
    double computeDPMCompTerm();
    double computeTHzCompTerm();
    double computeVacDiffCompTerm();      // [PAPER_294] hbar-denom quantum-volume coupling
    double computeSuperCompTerm();        // [PAPER_295] compressed Cooper seeding, f_DPM^2 scaling
    double computeAetherResTerm();
    double computeU_g4iResTerm();
    double computeOscResTerm(double t);   // PAPER_288-pattern standing+traveling
    double computeQuantumResTerm();
    double computeFluidResTerm();
    double computeExpResTerm();

    void updateCache();

public:
    CompressedResonanceUQFF24Module();

    // --- Main legacy API (original signature preserved) ---
    double computeCompressedResTerm(double t, double B);

    // --- UQFF 2.0 alias ---
    double computeG(double t, double B) { return computeCompressedResTerm(t, B); }

    // --- UQFF 2.0 runtime configuration ---
    void setEnableLogging(bool enable)               { logging_enabled = enable; }
    bool getLoggingEnabled()                   const { return logging_enabled; }
    void setDynamicParameter(const std::string& name, double value)
                                                     { dynamic_params[name] = value; }
    double getDynamicParameter(const std::string& name) const {
        auto it = dynamic_params.find(name);
        return (it != dynamic_params.end()) ? it->second
                                            : std::numeric_limits<double>::quiet_NaN();
    }

    // --- Legacy variable helpers (compatibility shims) ---
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // --- State persistence ---
    void exportState(const std::string& filename =
                     "CompressedResonanceUQFF24_UQFF_state.txt") const;

    // --- Cross-validation template (UQFF 2.0) ---
    template <typename OtherModule>
    double cross_validate(OtherModule& other, double t, double B) {
        double g_self  = computeCompressedResTerm(t, B);
        double g_other = other.computeG(t, B);
        double ratio   = (g_other != 0.0) ? g_self / g_other : 0.0;
        if (logging_enabled) {
            std::cout << "[CR24 XVAL] g_self=" << std::scientific << g_self
                      << " g_other=" << g_other
                      << " ratio="   << ratio << "\n";
        }
        return ratio;
    }

    std::string getEquationText() const;
    void printVariables() const;
};

// ===========================================================================
// Implementation
// ===========================================================================

void CompressedResonanceUQFF24Module::updateCache() {
    F_DPM_cache      = I_curr * A_vort * (omega_1 - omega_2);           // 6.284e36 N
    a_DPM_cache      = (F_DPM_cache * f_DPM * E_VAC) / (C_LIGHT * V_sys);
    // Γ_THz = 10 * f_THz * v_exp / c (same as PAPER_287 pattern, systems-18-24 scale)
    Gamma_THz_cache  = 10.0 * f_THz * v_exp / C_LIGHT;
    // [PAPER_295] A_sc = hbar * f_super * f_DPM / (E_vac * c); a_super = A_sc * a_DPM
    // Since a_DPM ∝ f_DPM as well → a_super ∝ f_DPM^2 (quadratic DPM-class scaling law)
    A_sc_cache       = (HBAR * f_super * f_DPM) / (E_VAC * C_LIGHT);
    // [PAPER_294] a_vac_diff = (E_0 * f_vac_diff * V_sys * a_DPM) / hbar
    // Static component (no t-dependency) — cached at construction / updateVariable
    a_vac_diff_cache = (E_0 * f_vac_diff * V_sys * a_DPM_cache) / HBAR;
    E_0_over_E_vac   = E_0 / E_VAC;                                     // 0.9001
}

CompressedResonanceUQFF24Module::CompressedResonanceUQFF24Module()
    : f_DPM(1.0e11), f_THz(1.0e11), f_vac_diff(0.143), f_super(1.411e15),
      I_curr(1.0e20), A_vort(3.142e18), omega_1(1.0e-2), omega_2(-1.0e-2),
      v_exp(1.0e5), E_0(6.381e-36), V_sys(4.189e18),
      f_aether(1.0e3), f_react(1.0e9), f_quantum(1.445e-17),
      f_fluid(1.269e-14), f_exp_freq(1.373e-8),
      omega_osc(1.0e14), k_wave(1.0e18), A_amp(1.0e-9), V_fluid(1.0e6),
      f_TRZ(0.1), B_crit(1.0e11), f_sc(1.0),
      F_DPM_cache(0.0), a_DPM_cache(0.0), A_sc_cache(0.0),
      Gamma_THz_cache(0.0), a_vac_diff_cache(0.0), E_0_over_E_vac(0.0),
      logging_enabled(false)
{
    updateCache();
}

// ---------------------------------------------------------------------------
// Compressed-channel helpers
// ---------------------------------------------------------------------------

double CompressedResonanceUQFF24Module::computeDPMCompTerm() {
    // Shared DPM seed — identical formula for both channels
    return a_DPM_cache;
}

double CompressedResonanceUQFF24Module::computeTHzCompTerm() {
    // a_THz = Gamma_THz * a_DPM  (same PAPER_287 pattern, different system scale)
    return Gamma_THz_cache * a_DPM_cache;
}

double CompressedResonanceUQFF24Module::computeVacDiffCompTerm() {
    // [PAPER_294] a_vac_diff = (E_0 * f_vac_diff * V_sys * a_DPM) / hbar
    // FIRST UQFF term with hbar in denominator — quantum-volume diffusion coupling.
    // f_vac_diff = 0.143 Hz → T_vac = 6.993 s ≈ 7-second vacuum beat period.
    // E_0 = 6.381e-36 J/m^3; E_0/E_vac = 0.9001 — 10% plasmotic vacuum deficit channel.
    // V_sys / hbar = 4.189e18 / 1.0546e-34 = 3.973e52 m^3/(J*s) — quantum-volume constant.
    return a_vac_diff_cache;
}

double CompressedResonanceUQFF24Module::computeSuperCompTerm() {
    // [PAPER_295] a_super = A_sc * a_DPM
    // A_sc = hbar * f_super * f_DPM / (E_vac * c)
    // A_sc ∝ f_DPM; a_DPM ∝ f_DPM → a_super ∝ f_DPM^2 (quadratic frequency-class scaling).
    // Systems 18-24 (f_DPM=1e11): A_sc = 6.994e18   a_super = 2.479e4 m/s^2
    // Magnetar class (f_DPM=1e12): A_sc = 6.994e21   a_super = 2.479e8 m/s^2  (+4 orders)
    // Placed in COMPRESSED channel (pre-oscillatory seeding) — distinct from PAPER_289
    // where A_sc*a_DPM is in the RESONANCE channel (post-THz cascade synthesis).
    return A_sc_cache * a_DPM_cache;
}

// ---------------------------------------------------------------------------
// Resonance-channel helpers
// ---------------------------------------------------------------------------

double CompressedResonanceUQFF24Module::computeAetherResTerm() {
    return f_aether * 1.0e-8 * f_DPM * (1.0 + f_TRZ) * a_DPM_cache;
}

double CompressedResonanceUQFF24Module::computeU_g4iResTerm() {
    return f_sc * f_react * a_DPM_cache / (E_VAC * C_LIGHT);
}

double CompressedResonanceUQFF24Module::computeOscResTerm(double t) {
    // PAPER_288-inherited standing + traveling wave pair (cosmic-age normalisation)
    double a_standing  = 2.0 * A_amp * std::cos(k_wave * 0.0) * std::cos(omega_osc * t);
    std::complex<double> phase(0.0, k_wave * 0.0 - omega_osc * t);
    double a_traveling = (2.0 * PI / T_COSMIC_GYR) * A_amp * std::exp(phase).real();
    return a_standing + a_traveling;
}

double CompressedResonanceUQFF24Module::computeQuantumResTerm() {
    const double E_vac_ISM = E_VAC / 10.0;
    return (f_quantum * E_VAC * a_DPM_cache) / (E_vac_ISM * C_LIGHT);
}

double CompressedResonanceUQFF24Module::computeFluidResTerm() {
    const double E_vac_ISM = E_VAC / 10.0;
    return (f_fluid * E_VAC * V_fluid * a_DPM_cache) / (E_vac_ISM * C_LIGHT);
}

double CompressedResonanceUQFF24Module::computeExpResTerm() {
    const double E_vac_ISM = E_VAC / 10.0;
    return (f_exp_freq * E_VAC * a_DPM_cache) / (E_vac_ISM * C_LIGHT);
}

// ---------------------------------------------------------------------------
// Main computation — legacy API preserved
// ---------------------------------------------------------------------------

double CompressedResonanceUQFF24Module::computeCompressedResTerm(double t, double B) {
    // === Compressed channel (4 terms) [PAPER_293] ===
    const double a_DPM      = computeDPMCompTerm();
    const double a_THz      = computeTHzCompTerm();
    const double a_vac_diff = computeVacDiffCompTerm();    // [PAPER_294] hbar denominator
    const double a_super    = computeSuperCompTerm();      // [PAPER_295] f_DPM^2 scaling
    const double a_comp     = a_DPM + a_THz + a_vac_diff + a_super;

    // === Resonance channel (6 terms) [PAPER_293] ===
    const double a_aether  = computeAetherResTerm();
    const double a_u_g4i   = computeU_g4iResTerm();
    const double a_osc     = computeOscResTerm(t);
    const double a_quantum = computeQuantumResTerm();
    const double a_fluid   = computeFluidResTerm();
    const double a_exp     = computeExpResTerm();
    const double a_res     = a_aether + a_u_g4i + a_osc + a_quantum + a_fluid + a_exp;

    // === SC correction + TRZ factor ===
    const double SCm  = 1.0 - B / B_crit;
    const double g_CR = (a_comp + a_res) * SCm * (1.0 + f_TRZ);

    // === Channel dominance ratio [PAPER_293] ===
    const double R_CR = (a_res != 0.0) ? a_comp / a_res : 0.0;

    if (logging_enabled) {
        std::cout << std::scientific << std::setprecision(6)
                  << "[CR24][COMPRESSED]  a_DPM="      << a_DPM
                  << "  a_THz="       << a_THz
                  << "  a_vac_diff="  << a_vac_diff   << " [PAPER_294 hbar-denom]"
                  << "  a_super="     << a_super      << " [PAPER_295 f_DPM^2]"
                  << "  => a_comp="   << a_comp       << "\n"
                  << "[CR24][RESONANCE]   a_aether="   << a_aether
                  << "  a_u_g4i="     << a_u_g4i
                  << "  a_osc="       << a_osc
                  << "  a_quantum="   << a_quantum
                  << "  a_fluid="     << a_fluid
                  << "  a_exp="       << a_exp
                  << "  => a_res="    << a_res        << "\n"
                  << "[CR24][DUAL]        R_CR(comp/res)=" << R_CR << " [PAPER_293]\n"
                  << "[CR24][SC]          SCm="      << SCm
                  << "  B="          << B
                  << "  B_crit="     << B_crit       << "\n"
                  << "[CR24][RESULT]      g_CR="     << g_CR << " m/s^2\n";
    }

    return g_CR;
}

// ---------------------------------------------------------------------------
// UQFF 2.0 state export (40 parameters)
// ---------------------------------------------------------------------------

void CompressedResonanceUQFF24Module::exportState(const std::string& filename) const {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "[CR24] Cannot open state file: " << filename << "\n";
        return;
    }
    ofs << std::scientific << std::setprecision(10);
    ofs << "# CompressedResonanceUQFF24Module — UQFF 2.0 State (Session 83)\n";
    ofs << "# Papers: PAPER_293 PAPER_294 PAPER_295\n";
    ofs << "# Module:  25th C++ UQFF module — FIRST dual-channel architecture\n";
    ofs << "C_LIGHT           = " << C_LIGHT          << "\n";
    ofs << "HBAR              = " << HBAR             << "\n";
    ofs << "E_VAC             = " << E_VAC            << "\n";
    ofs << "T_COSMIC_GYR      = " << T_COSMIC_GYR     << "\n";
    ofs << "f_DPM             = " << f_DPM            << "  # Hz systems-18-24 class = 1e11\n";
    ofs << "f_THz             = " << f_THz            << "\n";
    ofs << "f_vac_diff        = " << f_vac_diff       << "  # Hz 0.143 -> T_vac=6.993s [PAPER_294]\n";
    ofs << "f_super           = " << f_super          << "  # Hz Cooper pair [PAPER_295]\n";
    ofs << "I_curr            = " << I_curr           << "\n";
    ofs << "A_vort            = " << A_vort           << "\n";
    ofs << "omega_1           = " << omega_1          << "\n";
    ofs << "omega_2           = " << omega_2          << "\n";
    ofs << "v_exp             = " << v_exp            << "\n";
    ofs << "E_0               = " << E_0              << "  # J/m^3 reduced vacuum [PAPER_294]\n";
    ofs << "V_sys             = " << V_sys            << "\n";
    ofs << "f_aether          = " << f_aether         << "\n";
    ofs << "f_react           = " << f_react          << "\n";
    ofs << "f_quantum         = " << f_quantum        << "\n";
    ofs << "f_fluid           = " << f_fluid          << "\n";
    ofs << "f_exp_freq        = " << f_exp_freq       << "\n";
    ofs << "omega_osc         = " << omega_osc        << "\n";
    ofs << "k_wave            = " << k_wave           << "\n";
    ofs << "A_amp             = " << A_amp            << "\n";
    ofs << "V_fluid           = " << V_fluid          << "\n";
    ofs << "f_TRZ             = " << f_TRZ            << "\n";
    ofs << "B_crit            = " << B_crit           << "\n";
    ofs << "f_sc              = " << f_sc             << "\n";
    ofs << "# --- Cached / derived values ---\n";
    ofs << "F_DPM_cache       = " << F_DPM_cache      << "  # N vortex force 6.284e36\n";
    ofs << "a_DPM_cache       = " << a_DPM_cache      << "  # m/s^2 DPM base seed 3.543e-15\n";
    ofs << "A_sc_cache        = " << A_sc_cache       << "  # PAPER_295 systems-18-24: 6.994e18\n";
    ofs << "Gamma_THz_cache   = " << Gamma_THz_cache  << "  # 3.333e8 THz amplification\n";
    ofs << "a_vac_diff_cache  = " << a_vac_diff_cache << "  # PAPER_294 hbar-denom 128.4 m/s^2\n";
    ofs << "E_0_over_E_vac    = " << E_0_over_E_vac   << "  # 0.9001 vacuum deficit ratio\n";
    ofs << "T_vac_s           = " << (1.0 / f_vac_diff) << "  # 6.993 s vacuum beat period\n";
    ofs << "a_super_cache     = " << (A_sc_cache * a_DPM_cache) << "  # 2.479e4 m/s^2 [PAPER_295]\n";
    ofs << "notes_PAPER_293   = Dual-channel 4+6 co-sum; R_CR=comp/res analytic observable\n";
    ofs << "notes_PAPER_294   = hbar in denominator; first UQFF quantum-volume diffusion\n";
    ofs << "notes_PAPER_295   = a_super prop f_DPM^2; pre-oscillatory compressed Cooper seeding\n";
    for (const auto& p : dynamic_params) {
        ofs << "dyn:" << p.first << " = " << p.second << "\n";
    }
    ofs.close();
}

// ---------------------------------------------------------------------------
// Legacy variable helpers — compatibility shims mapping old map API to named members
// ---------------------------------------------------------------------------

void CompressedResonanceUQFF24Module::updateVariable(const std::string& name, double value) {
    dynamic_params[name] = value;
    bool need_cache = false;
    if      (name == "f_DPM"    || name == "f_DPM")   { f_DPM       = value; need_cache = true; }
    else if (name == "f_THz")                          { f_THz       = value; need_cache = true; }
    else if (name == "f_vac_diff")                     { f_vac_diff  = value; need_cache = true; }
    else if (name == "f_super")                        { f_super     = value; need_cache = true; }
    else if (name == "I")                              { I_curr      = value; need_cache = true; }
    else if (name == "A_vort")                         { A_vort      = value; need_cache = true; }
    else if (name == "omega_1")                        { omega_1     = value; need_cache = true; }
    else if (name == "omega_2")                        { omega_2     = value; need_cache = true; }
    else if (name == "v_exp")                          { v_exp       = value; need_cache = true; }
    else if (name == "E_0")                            { E_0         = value; need_cache = true; }
    else if (name == "V_sys")                          { V_sys       = value; need_cache = true; }
    else if (name == "V" || name == "V_fluid")         { V_fluid     = value; }
    else if (name == "B_crit")                         { B_crit      = value; }
    else if (name == "f_TRZ")                          { f_TRZ       = value; }
    else if (name == "f_aether")                       { f_aether    = value; }
    else if (name == "f_react")                        { f_react     = value; }
    else if (name == "f_quantum")                      { f_quantum   = value; }
    else if (name == "f_fluid")                        { f_fluid     = value; }
    else if (name == "f_exp" || name == "f_exp_freq")  { f_exp_freq  = value; }
    else if (name == "omega_osc")                      { omega_osc   = value; }
    else if (name == "k" || name == "k_wave")          { k_wave      = value; }
    else if (name == "A")                              { A_amp       = value; }
    else if (name == "f_sc")                           { f_sc        = value; }
    if (need_cache) updateCache();
}

void CompressedResonanceUQFF24Module::addToVariable(const std::string& name, double delta) {
    auto it = dynamic_params.find(name);
    double current = (it != dynamic_params.end()) ? it->second : 0.0;
    updateVariable(name, current + delta);
}

void CompressedResonanceUQFF24Module::subtractFromVariable(const std::string& name,
                                                            double delta) {
    addToVariable(name, -delta);
}

// ---------------------------------------------------------------------------
// Descriptive output
// ---------------------------------------------------------------------------

std::string CompressedResonanceUQFF24Module::getEquationText() const {
    return
        "=== WOLFRAM_TERM_CR24_BASE ===\n"
        WOLFRAM_TERM_CR24_BASE "\n\n"
        "=== WOLFRAM_TERM_CR24_VAC_DIFF ===\n"
        WOLFRAM_TERM_CR24_VAC_DIFF "\n\n"
        "=== WOLFRAM_TERM_CR24_SUPER_COMP ===\n"
        WOLFRAM_TERM_CR24_SUPER_COMP "\n\n"
        "=== WOLFRAM_TERM_CR24_DUAL_CHANNEL ===\n"
        WOLFRAM_TERM_CR24_DUAL_CHANNEL "\n\n"
        "--- Dual-Channel Architecture [PAPER_293] ---\n"
        "  Compressed (4): a_DPM + a_THz + a_vac_diff + a_super\n"
        "  Resonance  (6): a_aether + a_u_g4i + a_osc + a_quantum + a_fluid + a_exp\n"
        "  g_CR = (a_comp + a_res) * (1-B/B_crit) * (1+f_TRZ)\n"
        "  R_CR = a_comp/a_res = 1.490e-17 (resonance dominates by 17 orders at sys18-24 params)\n\n"
        "--- Vacuum Differential Harmonic [PAPER_294] ---\n"
        "  a_vac_diff = (E_0 * f_vac_diff * V_sys * a_DPM) / hbar\n"
        "  f_vac_diff = 0.143 Hz  ->  T_vac = 6.993 s (vacuum beat period)\n"
        "  E_0 = 6.381e-36 J/m^3;  E_0/E_vac = 0.9001 (10% vacuum deficit)\n"
        "  V_sys/hbar = 3.973e52 m^3/(J*s)  (quantum-volume coupling constant)\n"
        "  a_vac_diff(default) = 128.4 m/s^2 (dominates DPM and THz in compressed channel)\n\n"
        "--- Compressed Cooper f_DPM^2 Scaling [PAPER_295] ---\n"
        "  a_super = A_sc * a_DPM;  A_sc = hbar * f_super * f_DPM / (E_vac * c)\n"
        "  f_DPM = 1e11 (sys 18-24): A_sc = 6.994e18  a_super = 2.479e4 m/s^2\n"
        "  f_DPM = 1e12 (magnetar):  A_sc = 6.994e21  a_super = 2.479e8 m/s^2\n"
        "  a_super / a_DPM = A_sc  =>  a_super prop f_DPM^2 (quadratic class scaling)";
}

void CompressedResonanceUQFF24Module::printVariables() const {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "CompressedResonanceUQFF24Module — UQFF 2.0 (Session 83, 25th C++ module)\n";
    std::cout << "  [Compressed]  f_DPM="    << f_DPM
              << "  f_THz="      << f_THz
              << "  f_vac_diff=" << f_vac_diff
              << "  f_super="    << f_super   << "\n";
    std::cout << "  [Vortex]      I="      << I_curr
              << "  A_vort="  << A_vort
              << "  omega_1=" << omega_1
              << "  omega_2=" << omega_2    << "\n";
    std::cout << "  [System]      V_sys="  << V_sys
              << "  v_exp="   << v_exp
              << "  E_0="     << E_0       << "\n";
    std::cout << "  [Resonance]   f_aether="  << f_aether
              << "  f_react="   << f_react
              << "  f_quantum=" << f_quantum
              << "  f_fluid="   << f_fluid
              << "  f_exp="     << f_exp_freq << "\n";
    std::cout << "  [Osc]         omega_osc=" << omega_osc
              << "  k_wave="  << k_wave
              << "  A_amp="   << A_amp
              << "  V_fluid=" << V_fluid  << "\n";
    std::cout << "  [SC/TRZ]      B_crit=" << B_crit
              << "  f_sc="    << f_sc
              << "  f_TRZ="   << f_TRZ    << "\n";
    std::cout << "  [Cache]       F_DPM="    << F_DPM_cache
              << "  a_DPM="      << a_DPM_cache
              << "  Gamma_THz="  << Gamma_THz_cache << "\n";
    std::cout << "  [PAPER_294]   a_vac_diff=" << a_vac_diff_cache
              << "  T_vac="      << (1.0 / f_vac_diff) << " s"
              << "  E_0/E_vac="  << E_0_over_E_vac  << "\n";
    std::cout << "  [PAPER_295]   A_sc="      << A_sc_cache
              << "  a_super="    << (A_sc_cache * a_DPM_cache)
              << " (f_DPM^2 class scaling)\n";
    if (!dynamic_params.empty()) {
        for (const auto& p : dynamic_params) {
            std::cout << "  [dyn]         " << p.first << "=" << p.second << "\n";
        }
    }
}

#endif // COMPRESSED_RESONANCE_UQFF24_MODULE_H