// ============================================================
// NGC6302_RESONANCE_UQFF_MODULE.cpp
// Session 90 — Star-Magic UQFF 2.0 Self-Expanding Framework
// 32nd C++ UQFF module — FIRST Resonance-Channel companion to
//   a UQFF Bipolar Planetary Nebula module (NGC 6302 "Bug Nebula")
// 11-term resonance co-sum (DPM/THz/VacDiff/SuperFreq/AetherRes/
//   U_g4i/QuantumFreq/AetherFreq/FluidFreq/Osc/ExpFreq)
// NGC 6302 resonance params:
//   r=1.42e16 m (~1.5 ly lobe half-span); rho=1e-21 kg/m³
//   v_exp=2.68e5 m/s (268 km/s, HST bipolar lobe expansion)
//   I=1e20 A (wind current proxy); f_DPM=1e12 Hz (wind-aligned)
//   E_vac_neb=7.09e-36 J/m³; E_vac_ISM=7.09e-37 J/m³
// Unique Physics computed:
//   PAPER_314: F_DPM=1.267e50 N (PN lobe macro-antenna DPM force;
//              ratio to compact F_DPM=6.284e36 N: 2.017e13; 13 orders)
//   PAPER_315: Gamma_THz=8.939e9; a_THz=2.232e-21 m/s²;
//              r_cross=3.280 km (VacDiff-THz crossover radius;
//              a_vac_diff/a_THz at PN scale=8.118e37, 38 orders)
//   PAPER_316: A_sc=6.994e21; a_super=1.747e-9 m/s²;
//              confirms PAPER_295 f_DPM=1e12 class prediction
// Stub fixes applied:
//   - A_area/A_osc collision resolved (A_area=pi*r² for DPM;
//     A_osc=1e-10 m canonical quantum amplitude for Osc term)
//   - B/B_crit explicit members (vs hardcoded 1e-5/1e11)
//   - Ug1_proxy uses M_proxy member (vs hardcoded constants)
//   - VacDiff canonical: (E_0*V_sys/hbar)*a_DPM [f_vac_diff cancels]
//   - std::map replaced by typed private members + dynamic_params
// WOLFRAM_TERM x4: NGC6302_RES_BASE / NGC6302_RES_DPM_LOBE /
//                  NGC6302_RES_THz_EXPANSION / NGC6302_RES_COOPER_SC
// Copyright: Daniel T. Murphy. Session 90, March 17, 2026.
// ============================================================

#ifndef NGC6302_RESONANCE_UQFF_MODULE_H
#define NGC6302_RESONANCE_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <fstream>
#include <limits>

// ── WOLFRAM_TERM Integration Tags ─────────────────────────────
#define WOLFRAM_TERM_NGC6302_RES_BASE          "NGC6302_RES_BASE"
#define WOLFRAM_TERM_NGC6302_RES_DPM_LOBE      "NGC6302_RES_DPM_LOBE"
#define WOLFRAM_TERM_NGC6302_RES_THz_EXPANSION "NGC6302_RES_THz_EXPANSION"
#define WOLFRAM_TERM_NGC6302_RES_COOPER_SC     "NGC6302_RES_COOPER_SC"

// ──────────────────────────────────────────────────────────────
class NGC6302ResonanceUQFFModule {
public:
    // ── Static Physical Constants ─────────────────────────────
    static constexpr double G_NEWTON    = 6.6743e-11;   // m³/(kg·s²)
    static constexpr double C_LIGHT     = 2.998e8;      // m/s
    static constexpr double HBAR        = 1.0546e-34;   // J·s
    static constexpr double PI          = 3.14159265358979;
    static constexpr double E_VAC_NEB   = 7.09e-36;     // J/m³ (plasmotic, nebula)
    static constexpr double E_VAC_ISM   = 7.09e-37;     // J/m³ (ISM)
    static constexpr double VAC_RATIO   = 10.0;         // E_VAC_NEB / E_VAC_ISM

private:
    // ── NGC 6302 Lobe Geometry ────────────────────────────────
    double r;           // m  — lobe half-span radius
    double V_sys;       // m³ — lobe volume (4/3 π r³ = 1.199e49)
    double A_area;      // m² — lobe cross-section (π r² = 6.333e32)
    double rho_fluid;   // kg/m³ — lobe density
    double v_exp;       // m/s — bipolar expansion velocity (HST)

    // ── DPM Parameters ────────────────────────────────────────
    double I_wind;      // A   — current proxy from winds
    double omega_1;     // rad/s — DPM angular freq +
    double omega_2;     // rad/s — DPM angular freq -
    double delta_omega; // rad/s — omega_1 - omega_2 = 2e-3
    double f_DPM;       // Hz  — DPM intrinsic frequency (wind-aligned)

    // ── Resonance Frequencies ─────────────────────────────────
    double f_THz;          // Hz — THz hole resonance
    double f_vac_diff;     // Hz — vacuum differential (cancels in VacDiff formula)
    double f_super;        // Hz — Cooper-DPM superconductive
    double f_aether;       // Hz — Aether-mediated resonance
    double f_react;        // Hz — U_g4i reactive
    double f_quantum;      // Hz — quantum wave
    double f_Aether_eff;   // Hz — Aether effect (effective, very low f)
    double f_fluid;        // Hz — fluid resonance
    double f_exp;          // Hz — cosmic expansion resonance

    // ── Vacuum / Field Parameters ─────────────────────────────
    double E_0;            // J/m³ — vacuum differential energy density
    double Lambda_res;     // m⁻²  — Aether proxy
    double B;              // T    — lobe magnetic field
    double B_crit;         // T    — critical field (Meissner)
    double f_TRZ;          // dimensionless — time-reversal factor
    double f_sc;           // dimensionless — superconductive factor
    double V_fluid;        // m³  — fluid volume element for fluid resonance term

    // ── Oscillation Parameters (PAPER_288 canonical) ──────────
    // A_osc is SEPARATE from A_area to avoid the stub naming collision
    double A_osc;          // m   — oscillation amplitude (quantum scale = 1e-10 m)
    double k_osc;          // m⁻¹ — standing-wave wavenumber
    double omega_osc;      // rad/s — oscillation angular frequency

    // ── Perturbation ──────────────────────────────────────────
    double delta_rho;      // kg/m³ — density perturbation

    // ── Gravitational Proxy ───────────────────────────────────
    double M_proxy;        // kg  — PN ejected mass (2 M_sun from Session 89)
    double Ug1_proxy;      // m/s² — G*M_proxy/r² (used by U_g4i term)

    // ── Cached Unique Physics (PAPER_314/315/316) ────────────
    double F_DPM_cache;            // N      — PAPER_314 lobe force
    double a_DPM_cache;            // m/s²   — PAPER_314 seed term
    double ratio_FPN_compact;      // dimless — PAPER_314 amplification ratio
    double Gamma_THz_cache;        // dimless — PAPER_315 THz amplifier
    double a_THz_cache;            // m/s²   — PAPER_315
    double r_cross_cache;          // m      — PAPER_315 VacDiff-THz crossover
    double vac_THz_ratio_cache;    // dimless — PAPER_315 dominance at PN scale
    double A_sc_cache;             // dimless — PAPER_316 Cooper-DPM amplitude
    double a_super_cache;          // m/s²   — PAPER_316

    // ── UQFF 2.0 State ────────────────────────────────────────
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

    // ── Private Methods ───────────────────────────────────────
    void initializeDefaults();
    void updateCache();
    double computeDPMTerm();
    double computeTHzTerm();
    double computeVacDiffTerm();
    double computeSuperFreqTerm();
    double computeAetherResTerm();
    double computeU_g4iTerm();
    double computeQuantumFreqTerm();
    double computeAetherFreqTerm();
    double computeFluidFreqTerm();
    double computeOscTerm(double t);
    double computeExpFreqTerm();

public:
    // ── Constructor / UQFF 2.0 Interface ─────────────────────
    NGC6302ResonanceUQFFModule();

    void setEnableLogging(bool en) { logging_enabled = en; }
    bool getLoggingEnabled() const  { return logging_enabled; }

    void setDynamicParameter(const std::string& key, double val) {
        dynamic_params[key] = val;
    }
    double getDynamicParameter(const std::string& key) const {
        auto it = dynamic_params.find(key);
        return (it != dynamic_params.end()) ? it->second
                                            : std::numeric_limits<double>::quiet_NaN();
    }

    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta) {
        addToVariable(name, -delta);
    }

    double computeG(double t);
    std::string getEquationText();
    void printParameters();
    void printVariables() { printParameters(); }  // legacy alias
    void exportState(const std::string& filename);

    template<typename OtherModule>
    double cross_validate(OtherModule& other, double t) {
        double g_this  = this->computeG(t);
        double g_other = other.computeG(t);
        double ratio   = (g_other != 0.0) ? g_this / g_other : 0.0;
        if (logging_enabled)
            std::cout << "[CROSS_VALIDATE] NGC6302_RES g=" << g_this
                      << " | other g=" << g_other
                      << " | ratio=" << ratio << "\n";
        return ratio;
    }
};

// ── Constructor ───────────────────────────────────────────────
NGC6302ResonanceUQFFModule::NGC6302ResonanceUQFFModule() {
    logging_enabled = false;
    initializeDefaults();
    updateCache();
}

void NGC6302ResonanceUQFFModule::initializeDefaults() {
    // Lobe geometry
    r           = 1.42e16;                        // m  (~1.5 ly half-span)
    V_sys       = (4.0/3.0) * PI * r*r*r;         // 1.199e49 m³
    A_area      = PI * r * r;                      // 6.333e32 m²
    rho_fluid   = 1e-21;                           // kg/m³ (lobe density)
    v_exp       = 2.68e5;                          // m/s (268 km/s, HST)

    // DPM
    I_wind      = 1e20;                            // A
    omega_1     = 1e-3;                            // rad/s
    omega_2     = -1e-3;                           // rad/s
    delta_omega = omega_1 - omega_2;               // 2e-3 rad/s
    f_DPM       = 1e12;                            // Hz (wind-aligned, 1e12 class)

    // Resonance frequencies
    f_THz       = 1e12;                            // Hz
    f_vac_diff  = 0.143;                           // Hz  (note: cancels in VacDiff)
    f_super     = 1.411e16;                        // Hz  (Cooper)
    f_aether    = 1e4;                             // Hz
    f_react     = 1e10;                            // Hz
    f_quantum   = 1.445e-17;                       // Hz
    f_Aether_eff= 1.576e-35;                       // Hz
    f_fluid     = 1.269e-14;                       // Hz
    f_exp       = 1.373e-8;                        // Hz

    // Energy / field
    E_0         = 6.381e-36;                       // J/m³
    Lambda_res  = 1.1e-52;                         // m⁻²
    B           = 1e-5;                            // T   (lobe field)
    B_crit      = 1e11;                            // T   (Meissner critical)
    f_TRZ       = 0.1;
    f_sc        = 1.0;
    V_fluid     = 1e3;                             // m³ (fluid resonance volume element)

    // Oscillation — PAPER_288 canonical A_osc=1e-10 m (quantum scale)
    // A_osc is NOT A_area; fix for stub naming collision
    A_osc       = 1e-10;                           // m  (quantum standing-wave amplitude)
    k_osc       = 2.0*PI / r;                      // m⁻¹
    omega_osc   = 2.0*PI * C_LIGHT / r;           // rad/s

    // Perturbation
    delta_rho   = 0.1 * rho_fluid;                // kg/m³

    // Gravitational proxy — M_PN = 2 M_sun (Session 89 NGC 6302 value)
    M_proxy     = 2.0 * 1.989e30;                 // 3.978e30 kg
    Ug1_proxy   = G_NEWTON * M_proxy / (r * r);   // m/s²
}

void NGC6302ResonanceUQFFModule::updateCache() {
    // Recompute geometry
    V_sys      = (4.0/3.0) * PI * r*r*r;
    A_area     = PI * r * r;
    k_osc      = 2.0*PI / r;
    omega_osc  = 2.0*PI * C_LIGHT / r;
    delta_omega= omega_1 - omega_2;
    Ug1_proxy  = G_NEWTON * M_proxy / (r * r);

    // PAPER_314 [WOLFRAM: NGC6302_RES_DPM_LOBE]
    // F_DPM = I_wind * A_area * delta_omega = 1.267e50 N
    F_DPM_cache      = I_wind * A_area * delta_omega;
    // a_DPM = F_DPM * f_DPM * E_vac_neb / (c * V_sys) = 2.497e-31 m/s²
    a_DPM_cache      = (F_DPM_cache * f_DPM * E_VAC_NEB) / (C_LIGHT * V_sys);
    // Ratio to compact DPM force (PAPER_293 F_DPM=6.284e36 N)
    ratio_FPN_compact = F_DPM_cache / 6.284e36;   // 2.017e13

    // PAPER_315 [WOLFRAM: NGC6302_RES_THz_EXPANSION]
    // Gamma_THz = VAC_RATIO * f_THz * v_exp / c = 8.939e9
    Gamma_THz_cache  = VAC_RATIO * f_THz * v_exp / C_LIGHT;
    // a_THz = Gamma_THz * a_DPM = 2.232e-21 m/s²
    a_THz_cache      = Gamma_THz_cache * a_DPM_cache;
    // VacDiff-THz crossover radius: E_0*(4π/3)*r³_cross/hbar = Gamma_THz
    // => r³_cross = 3*hbar*Gamma_THz / (4π*E_0)
    double r3_cross  = (3.0 * HBAR * Gamma_THz_cache) / (4.0 * PI * E_0);
    r_cross_cache    = std::cbrt(r3_cross);       // 3.280e3 m = 3.280 km
    // Dominance ratio at this system's PN scale
    vac_THz_ratio_cache = (E_0 * V_sys / HBAR) / Gamma_THz_cache;  // 8.118e37

    // PAPER_316 [WOLFRAM: NGC6302_RES_COOPER_SC]
    // A_sc = hbar * f_super * f_DPM / (E_vac_ISM * c) = 6.994e21
    A_sc_cache    = (HBAR * f_super * f_DPM) / (E_VAC_ISM * C_LIGHT);
    // a_super = A_sc * a_DPM = 1.747e-9 m/s²
    a_super_cache = A_sc_cache * a_DPM_cache;
}

// ── Term Computations ─────────────────────────────────────────
// PAPER_314 | WOLFRAM: NGC6302_RES_DPM_LOBE
// F_DPM = I_wind * A_area * delta_omega; a_DPM = F_DPM*f_DPM*E_vac_neb/(c*V_sys)
double NGC6302ResonanceUQFFModule::computeDPMTerm() {
    return a_DPM_cache;
}

// PAPER_315 | WOLFRAM: NGC6302_RES_THz_EXPANSION
// a_THz = Gamma_THz * a_DPM; Gamma_THz = VAC_RATIO * f_THz * v_exp / c
double NGC6302ResonanceUQFFModule::computeTHzTerm() {
    return a_THz_cache;
}

// PAPER_315 canonical VacDiff — f_vac_diff cancels; result: (E_0*V_sys/hbar)*a_DPM
// WOLFRAM: NGC6302_RES_BASE
double NGC6302ResonanceUQFFModule::computeVacDiffTerm() {
    return (E_0 * V_sys / HBAR) * a_DPM_cache;
}

// PAPER_316 | WOLFRAM: NGC6302_RES_COOPER_SC
// a_super = A_sc * a_DPM; A_sc = hbar*f_super*f_DPM/(E_vac_ISM*c) = 6.994e21
double NGC6302ResonanceUQFFModule::computeSuperFreqTerm() {
    return a_super_cache;
}

// a_aether_res = f_aether * (B/B_crit) * f_DPM * (1+f_TRZ) * a_DPM
double NGC6302ResonanceUQFFModule::computeAetherResTerm() {
    return f_aether * (B / B_crit) * f_DPM * (1.0 + f_TRZ) * a_DPM_cache;
}

// a_u_g4i = f_sc * Ug1_proxy * f_react * a_DPM / (E_vac_ISM * c)
double NGC6302ResonanceUQFFModule::computeU_g4iTerm() {
    return f_sc * Ug1_proxy * f_react * a_DPM_cache / (E_VAC_ISM * C_LIGHT);
}

// a_quantum = f_quantum * VAC_RATIO * a_DPM
double NGC6302ResonanceUQFFModule::computeQuantumFreqTerm() {
    return f_quantum * VAC_RATIO * a_DPM_cache;
}

// a_Aether = f_Aether_eff * VAC_RATIO * a_DPM
double NGC6302ResonanceUQFFModule::computeAetherFreqTerm() {
    return f_Aether_eff * VAC_RATIO * a_DPM_cache;
}

// a_fluid = (f_fluid * E_vac_neb * V_fluid * rho_fluid) / (E_vac_ISM * c)
// V_fluid = 1e3 m³ (fluid resonance volume element, distinct from V_sys)
double NGC6302ResonanceUQFFModule::computeFluidFreqTerm() {
    return (f_fluid * E_VAC_NEB * V_fluid * rho_fluid) / (E_VAC_ISM * C_LIGHT);
}

// Osc = 2*A_osc*cos(omega_osc*t) + (2π/13.8)*A_osc*cos(-omega_osc*t)
// A_osc = 1e-10 m (PAPER_288 canonical quantum-scale amplitude)
// Standing + traveling wave hybrid; 13.8 = T_universe/Gyr normalization
double NGC6302ResonanceUQFFModule::computeOscTerm(double t) {
    double standing = 2.0 * A_osc * std::cos(omega_osc * t);
    double traveling = (2.0 * PI / 13.8) * A_osc * std::cos(-omega_osc * t);
    return standing + traveling;
}

// a_exp = f_exp * VAC_RATIO * a_DPM
double NGC6302ResonanceUQFFModule::computeExpFreqTerm() {
    return f_exp * VAC_RATIO * a_DPM_cache;
}

// ── Full 11-term Resonance Co-sum ─────────────────────────────
// g_NGC6302_RES(t) = [Σ 11 terms] × (1 + f_TRZ)
// WOLFRAM: NGC6302_RES_BASE
double NGC6302ResonanceUQFFModule::computeG(double t) {
    double a_DPM      = computeDPMTerm();
    double a_THz      = computeTHzTerm();
    double a_vac_diff = computeVacDiffTerm();
    double a_super    = computeSuperFreqTerm();
    double a_aether_r = computeAetherResTerm();
    double a_u_g4i    = computeU_g4iTerm();
    double a_quantum  = computeQuantumFreqTerm();
    double a_aether_f = computeAetherFreqTerm();
    double a_fluid    = computeFluidFreqTerm();
    double a_osc      = computeOscTerm(t);
    double a_exp_f    = computeExpFreqTerm();

    double g_sum   = a_DPM + a_THz + a_vac_diff + a_super
                   + a_aether_r + a_u_g4i + a_quantum
                   + a_aether_f + a_fluid + a_osc + a_exp_f;
    double g_total = g_sum * (1.0 + f_TRZ);

    if (logging_enabled) {
        std::cout << std::scientific << std::setprecision(4)
            << "[LOG NGC6302_RES] t=" << t << " s\n"
            << "  [P314] a_DPM          = " << a_DPM      << " m/s²\n"
            << "  [P315] a_THz          = " << a_THz      << " m/s²\n"
            << "  [P315] a_vac_diff     = " << a_vac_diff << " m/s²  (VacDiff dominant)\n"
            << "  [P316] a_super        = " << a_super    << " m/s²\n"
            << "         a_aether_res   = " << a_aether_r << " m/s²\n"
            << "         a_u_g4i        = " << a_u_g4i    << " m/s²\n"
            << "         a_quantum      = " << a_quantum  << " m/s²\n"
            << "         a_aether_freq  = " << a_aether_f << " m/s²\n"
            << "         a_fluid        = " << a_fluid    << " m/s²\n"
            << "         a_osc          = " << a_osc      << " m/s²\n"
            << "         a_exp          = " << a_exp_f    << " m/s²\n"
            << "  g_total (×1+f_TRZ)   = " << g_total    << " m/s²\n";
    }
    return g_total;
}

// ── getEquationText ───────────────────────────────────────────
std::string NGC6302ResonanceUQFFModule::getEquationText() {
    return
        "g_NGC6302_RES(t) = [a_DPM + a_THz + a_vac_diff + a_super + a_aether_res\n"
        "                  + a_u_g4i + a_quantum + a_Aether + a_fluid + a_osc + a_exp] × (1+f_TRZ)\n"
        "PAPER_314 [WOLFRAM: NGC6302_RES_DPM_LOBE]:\n"
        "  F_DPM = I_wind × A_area × Δω = 1e20 × 6.333e32 × 2e-3 = 1.267e50 N\n"
        "  a_DPM = F_DPM × f_DPM × E_vac_neb / (c × V_sys) = 2.497e-31 m/s²\n"
        "  F_PN/F_compact = 1.267e50 / 6.284e36 = 2.017e13 (13-order amplification)\n"
        "  FIRST UQFF DPM force at planetary nebula lobe scale (r~1.5 ly)\n"
        "PAPER_315 [WOLFRAM: NGC6302_RES_THz_EXPANSION]:\n"
        "  Γ_THz = (E_vac_neb/E_vac_ISM) × f_THz × v_exp / c = 10×1e12×2.68e5/c = 8.939e9\n"
        "  a_THz = Γ_THz × a_DPM = 2.232e-21 m/s²\n"
        "  r_cross = [3ħΓ_THz/(4πE_0)]^(1/3) = 3.280e3 m = 3.280 km\n"
        "  (for r > r_cross VacDiff dominates; at PN scale ratio = 8.118e37)\n"
        "  Confirmation: v_exp-scaling law Γ_THz∝v_exp: NGC6302/Crab = 8.939e9/5.0e10 = 0.179\n"
        "    = v_exp_NGC/v_exp_Crab = 2.68e5/1.5e6 = 0.179 EXACT ✓\n"
        "PAPER_316 [WOLFRAM: NGC6302_RES_COOPER_SC]:\n"
        "  A_sc = ħ × f_super × f_DPM / (E_vac_ISM × c) = 6.994e21\n"
        "    → confirms PAPER_295 f_DPM=1e12 class prediction ✓\n"
        "  a_super = A_sc × a_DPM = 1.747e-9 m/s²\n"
        "  PN-scale hierarchy: a_vac_diff >> a_super >> a_THz >> a_DPM\n"
        "Stub fixes: A_area/A_osc collision; B/B_crit explicit; Ug1_proxy member;\n"
        "  VacDiff canonical (f_vac_diff cancels); std::map → typed private members.";
}

// ── printParameters ───────────────────────────────────────────
void NGC6302ResonanceUQFFModule::printParameters() {
    std::cout << std::scientific << std::setprecision(4)
        << "=== NGC6302ResonanceUQFFModule Parameters (UQFF 2.0) ===\n"
        << "r           = " << r            << " m  (~1.5 ly lobe half-span)\n"
        << "V_sys       = " << V_sys        << " m³\n"
        << "A_area      = " << A_area       << " m²\n"
        << "rho_fluid   = " << rho_fluid    << " kg/m³\n"
        << "v_exp       = " << v_exp        << " m/s  (268 km/s HST)\n"
        << "I_wind      = " << I_wind       << " A\n"
        << "delta_omega = " << delta_omega  << " rad/s\n"
        << "f_DPM       = " << f_DPM        << " Hz  (1e12 class)\n"
        << "f_THz       = " << f_THz        << " Hz\n"
        << "f_super     = " << f_super      << " Hz  (Cooper)\n"
        << "f_aether    = " << f_aether     << " Hz\n"
        << "B           = " << B            << " T\n"
        << "B_crit      = " << B_crit       << " T\n"
        << "f_TRZ       = " << f_TRZ        << "\n"
        << "A_osc       = " << A_osc        << " m   (PAPER_288 quantum amplitude)\n"
        << "M_proxy     = " << M_proxy      << " kg  (2 M_sun)\n"
        << "Ug1_proxy   = " << Ug1_proxy    << " m/s²\n"
        << "--- Cached Paper Values ---\n"
        << "F_DPM       = " << F_DPM_cache       << " N      [PAPER_314]\n"
        << "a_DPM       = " << a_DPM_cache       << " m/s²   [PAPER_314]\n"
        << "F_PN/F_cpt  = " << ratio_FPN_compact << "        [PAPER_314: 13-order amplification]\n"
        << "Gamma_THz   = " << Gamma_THz_cache   << "        [PAPER_315]\n"
        << "a_THz       = " << a_THz_cache       << " m/s²   [PAPER_315]\n"
        << "r_cross     = " << r_cross_cache     << " m = "
                            << r_cross_cache/1e3 << " km [PAPER_315]\n"
        << "VD/THz_PN   = " << vac_THz_ratio_cache << "    [PAPER_315: 38-order dominance]\n"
        << "A_sc        = " << A_sc_cache        << "        [PAPER_316]\n"
        << "a_super     = " << a_super_cache     << " m/s²   [PAPER_316]\n";
}

// ── updateVariable ────────────────────────────────────────────
void NGC6302ResonanceUQFFModule::updateVariable(const std::string& name, double value) {
    if      (name == "r")        { r = value;        updateCache(); }
    else if (name == "rho")      { rho_fluid = value; }
    else if (name == "v_exp")    { v_exp = value;    updateCache(); }
    else if (name == "I_wind")   { I_wind = value;   updateCache(); }
    else if (name == "f_DPM")    { f_DPM = value;    updateCache(); }
    else if (name == "f_THz")    { f_THz = value;    updateCache(); }
    else if (name == "f_super")  { f_super = value;  updateCache(); }
    else if (name == "f_aether") { f_aether = value; }
    else if (name == "B")        { B = value; }
    else if (name == "f_TRZ")    { f_TRZ = value; }
    else if (name == "A_osc")    { A_osc = value; }
    else if (name == "M_proxy")  { M_proxy = value;  updateCache(); }
    else {
        dynamic_params[name] = value;
        if (logging_enabled)
            std::cout << "[NGC6302_RES] dynamic param: " << name << " = " << value << "\n";
    }
}

void NGC6302ResonanceUQFFModule::addToVariable(const std::string& name, double delta) {
    if      (name == "r")     { r += delta;     updateCache(); }
    else if (name == "v_exp") { v_exp += delta; updateCache(); }
    else if (name == "f_DPM") { f_DPM += delta; updateCache(); }
    else                      { dynamic_params[name] += delta; }
}

// ── exportState ───────────────────────────────────────────────
void NGC6302ResonanceUQFFModule::exportState(const std::string& filename) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "[NGC6302_RES] exportState: cannot open " << filename << "\n";
        return;
    }
    ofs << std::scientific << std::setprecision(6)
        << "# NGC6302ResonanceUQFFModule State — Session 90\n"
        << "r = " << r << "\n"
        << "V_sys = " << V_sys << "\n"
        << "A_area = " << A_area << "\n"
        << "rho_fluid = " << rho_fluid << "\n"
        << "v_exp = " << v_exp << "\n"
        << "I_wind = " << I_wind << "\n"
        << "omega_1 = " << omega_1 << "\n"
        << "omega_2 = " << omega_2 << "\n"
        << "delta_omega = " << delta_omega << "\n"
        << "f_DPM = " << f_DPM << "\n"
        << "f_THz = " << f_THz << "\n"
        << "f_vac_diff = " << f_vac_diff << "\n"
        << "f_super = " << f_super << "\n"
        << "f_aether = " << f_aether << "\n"
        << "f_react = " << f_react << "\n"
        << "f_quantum = " << f_quantum << "\n"
        << "f_Aether_eff = " << f_Aether_eff << "\n"
        << "f_fluid = " << f_fluid << "\n"
        << "f_exp = " << f_exp << "\n"
        << "E_0 = " << E_0 << "\n"
        << "Lambda_res = " << Lambda_res << "\n"
        << "B = " << B << "\n"
        << "B_crit = " << B_crit << "\n"
        << "f_TRZ = " << f_TRZ << "\n"
        << "f_sc = " << f_sc << "\n"
        << "V_fluid = " << V_fluid << "\n"
        << "A_osc = " << A_osc << "\n"
        << "k_osc = " << k_osc << "\n"
        << "omega_osc = " << omega_osc << "\n"
        << "delta_rho = " << delta_rho << "\n"
        << "M_proxy = " << M_proxy << "\n"
        << "Ug1_proxy = " << Ug1_proxy << "\n"
        << "# Cached PAPER values\n"
        << "F_DPM_cache = " << F_DPM_cache << "\n"
        << "a_DPM_cache = " << a_DPM_cache << "\n"
        << "ratio_FPN_compact = " << ratio_FPN_compact << "\n"
        << "Gamma_THz_cache = " << Gamma_THz_cache << "\n"
        << "a_THz_cache = " << a_THz_cache << "\n"
        << "r_cross_cache = " << r_cross_cache << "\n"
        << "vac_THz_ratio_cache = " << vac_THz_ratio_cache << "\n"
        << "A_sc_cache = " << A_sc_cache << "\n"
        << "a_super_cache = " << a_super_cache << "\n";
    for (const auto& kv : dynamic_params)
        ofs << "dyn." << kv.first << " = " << kv.second << "\n";
    ofs.close();
    if (logging_enabled)
        std::cout << "[NGC6302_RES] State exported to " << filename << "\n";
}

#endif // NGC6302_RESONANCE_UQFF_MODULE_H