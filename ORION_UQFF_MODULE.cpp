// ORION_UQFF_MODULE.cpp
// Session 91 — UQFF 2.0 Full Upgrade: Orion Nebula M42/NGC 1976 (33rd C++ module)
// FIRST Trapezium OB-cluster driven HII region; r=1.18e17 m (~12.5 ly), M=2000 M_sun=3.978e33 kg
// SFR=1 M_sun/yr, v_wind=8e3 m/s (8 km/s HII front), t_age=3e5 yr, rho=1e-20 kg/m^3, z=0.00034
// PAPER_317: eta_wind=28.47; a_wind=5.424e-10 m/s^2; t_erosion=467 kyr [FIRST UQFF HII ram pressure dominance]
// PAPER_318: L_trap=7.656e31 W; a_rad=1.461e8 m/s^2; eta_rad=7.664e18 [FIRST UQFF Trapezium OB UV champagne flow]
// PAPER_319: t_cross=67730 yr; sSFR=5e-4 yr-1 (50x Lagoon); m_factor(t_age)=151 [FIRST UQFF compact-HII SFR binding phase transition]
// Watermark: Copyright - Daniel T. Murphy, Session 91 March 2026.

#ifndef ORION_UQFF_MODULE_H
#define ORION_UQFF_MODULE_H

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>

// WOLFRAM_TERM registration macros
#define WOLFRAM_TERM_ORION_BASE(val)          (val)   // Base gravity + Hz expansion [Session 91]
#define WOLFRAM_TERM_ORION_WIND_RAM(val)      (val)   // Trapezium stellar wind ram dominance [PAPER_317]
#define WOLFRAM_TERM_ORION_TRAPEZIUM_RAD(val) (val)   // Trapezium OB UV radiation dominance [PAPER_318]
#define WOLFRAM_TERM_ORION_SFR_BINDING(val)   (val)   // Compact-HII SFR binding crossover [PAPER_319]

class OrionUQFFModule {
private:
    // ── Physical constants ──────────────────────────────────────────────────
    static constexpr double G_NEWTON  = 6.6743e-11;
    static constexpr double C_LIGHT   = 2.998e8;
    static constexpr double HBAR      = 1.0546e-34;
    static constexpr double LAMBDA_CC = 1.114e-52;
    static constexpr double MU_0      = 1.2566e-6;
    static constexpr double PI        = 3.14159265358979;
    static constexpr double M_SUN     = 1.989e30;
    static constexpr double L_SUN     = 3.828e26;
    static constexpr double YEAR_TO_S = 3.15576e7;
    static constexpr double T_HUBBLE  = 13.8e9 * 3.15576e7;

    // ── Orion M42 system parameters ────────────────────────────────────────
    double M;           // 2000 M_sun = 3.978e33 kg
    double r;           // 1.18e17 m (~12.5 ly)
    double V_sys;       // (4/3)pi r^3 = 6.887e51 m^3
    double V_fluid;     // 1e3 m^3 (quantum/fluid coupler)
    double rho_fluid;   // 1e-20 kg/m^3
    double v_wind;      // 8e3 m/s (HII ionization front)
    double t_age;       // 9.467e12 s (3e5 yr)
    double SFR_yr;      // 1.0 M_sun/yr
    double M_visible;   // 0.15 M
    double M_DM;        // 0.85 M
    double delta_rho;   // 0.1 * rho_fluid

    // ── Trapezium OB cluster ────────────────────────────────────────────────
    double L_trap;      // 7.656e31 W (theta1 Ori C, 2e5 L_sun)

    // ── EM / magnetic ──────────────────────────────────────────────────────
    double B;           // 1e-5 T
    double B_crit;      // 1e11 T

    // ── Quantum / oscillatory ──────────────────────────────────────────────
    double Delta_x;     // 1e-10 m
    double A_osc;       // 1e-10 m (PAPER_288 canonical)
    double k_osc;       // 2pi/r
    double omega_osc;   // v_wind/r (PDR eigenfreq)

    // ── Cosmological ──────────────────────────────────────────────────────
    double H0_kms;      // 67.15 km/s/Mpc
    double Omega_m;     // 0.3
    double Omega_L;     // 0.7
    double z_redshift;  // 0.00034
    double Mpc_to_m;    // 3.086e22

    // ── UQFF correction factors ────────────────────────────────────────────
    double f_TRZ;       // 0.1
    double f_sc;        // 1.0

    // ── 3-tier buoyancy (Sgr A* at ~8.5 kpc) ─────────────────────────────
    double beta_i;      // 0.61
    double omega_g;     // 7.3e-16
    double U_UA;        // 1e-11
    double M_GC;        // 7.956e36 kg (Sgr A*, 4e6 M_sun)
    double r_GC;        // 2.623e20 m (~8.5 kpc)

    // ── UQFF 2.0 runtime ──────────────────────────────────────────────────
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

    // ── Cached key physics (PAPER_317/318/319) ─────────────────────────────
    double g_base_cache;
    double a_wind_t0;
    double eta_wind;
    double a_wind_tage;
    double eta_wind_tage;
    double t_erosion;
    double a_rad_cache;
    double eta_rad;
    double m_factor_tage;
    double m_factor_1myr;
    double t_consume;
    double sSFR;
    double t_cross;
    double binding_ratio_tage;

    // ── Private helpers ────────────────────────────────────────────────────
    void initializeDefaults();
    void updateCache();

    double computeHz();
    double computeMsfFactor(double t_s);
    double computeUgSum();
    double computeQuantumTerm();
    double computeEMTerm();
    double computeFluidTerm();
    double computeResonantTerm(double t);
    double computeDMTerm();
    double computeWindTerm(double t);
    double computeRadTerm();

public:
    OrionUQFFModule();

    void setEnableLogging(bool enable) { logging_enabled = enable; }
    bool getLoggingEnabled() const     { return logging_enabled; }

    void setDynamicParameter(const std::string& key, double val) {
        dynamic_params[key] = val;
    }
    double getDynamicParameter(const std::string& key, double def_val = 0.0) const {
        auto it = dynamic_params.find(key);
        return (it != dynamic_params.end()) ? it->second : def_val;
    }

    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    double computeG(double t);
    std::string getEquationText();
    void printParameters();
    void printVariables() { printParameters(); }

    void exportState(const std::string& filename = "ORION_state.txt");

    template<typename OtherModule>
    double cross_validate(OtherModule& other, double t_years = 1e6) {
        double t_s    = t_years * YEAR_TO_S;
        double g_this = this->computeG(t_s);
        double g_other = other.computeG(t_s);
        double ratio   = (std::abs(g_other) > 1e-300) ? g_this / g_other : 0.0;
        if (logging_enabled) {
            std::cout << "[ORION cross_validate] g_this=" << g_this
                      << " g_other=" << g_other << " ratio=" << ratio << "\n";
        }
        return ratio;
    }
};

// ═══════════════════════════════════════════════════════════════════════════
// IMPLEMENTATION
// ═══════════════════════════════════════════════════════════════════════════

inline OrionUQFFModule::OrionUQFFModule() {
    initializeDefaults();
    updateCache();
}

inline void OrionUQFFModule::initializeDefaults() {
    M           = 2000.0 * M_SUN;                    // 3.978e33 kg
    r           = 1.18e17;                           // m
    V_sys       = (4.0/3.0) * PI * r * r * r;        // 6.887e51 m^3
    V_fluid     = 1e3;
    rho_fluid   = 1e-20;
    v_wind      = 8e3;
    t_age       = 3e5 * YEAR_TO_S;                   // 9.467e12 s
    SFR_yr      = 1.0;
    M_visible   = 0.15 * M;
    M_DM        = 0.85 * M;
    delta_rho   = 0.1 * rho_fluid;
    L_trap      = 2e5 * L_SUN;                       // 7.656e31 W
    B           = 1e-5;
    B_crit      = 1e11;
    Delta_x     = 1e-10;
    A_osc       = 1e-10;                             // PAPER_288 canonical
    k_osc       = 2.0 * PI / r;
    omega_osc   = v_wind / r;                        // PDR eigenfreq = 6.780e-14 rad/s
    H0_kms      = 67.15;
    Omega_m     = 0.3;
    Omega_L     = 0.7;
    z_redshift  = 0.00034;
    Mpc_to_m    = 3.086e22;
    f_TRZ       = 0.1;
    f_sc        = 1.0;
    beta_i      = 0.61;
    omega_g     = 7.3e-16;
    U_UA        = 1e-11;
    M_GC        = 4e6 * M_SUN;                       // Sgr A* = 7.956e36 kg
    r_GC        = 2.623e20;                          // ~8.5 kpc
    logging_enabled = false;
}

inline void OrionUQFFModule::updateCache() {
    // PAPER_317
    g_base_cache   = G_NEWTON * M / (r * r);
    a_wind_t0      = (v_wind * v_wind) / r;           // 5.424e-10 m/s^2
    eta_wind       = a_wind_t0 / g_base_cache;        // 28.47
    a_wind_tage    = a_wind_t0 * 2.0;                 // (1 + t_age/t_age)
    eta_wind_tage  = a_wind_tage / g_base_cache;      // 56.9
    t_erosion      = r / v_wind;                      // 4.675e13 s = 467 kyr

    // PAPER_318
    double A_trap  = 4.0 * PI * r * r;
    double P_rad   = L_trap / (A_trap * C_LIGHT);
    a_rad_cache    = P_rad / rho_fluid;               // 1.461e8 m/s^2
    eta_rad        = a_rad_cache / g_base_cache;      // 7.664e18

    // PAPER_319
    double M_sun_count   = M / M_SUN;                 // 2000
    double M_sf_tage     = SFR_yr * (t_age / YEAR_TO_S) / M_sun_count; // 150
    m_factor_tage        = 1.0 + M_sf_tage;           // 151
    double M_sf_1myr     = SFR_yr * 1e6 / M_sun_count;
    m_factor_1myr        = 1.0 + M_sf_1myr;           // 501
    t_consume            = M_sun_count / SFR_yr;      // 2000 yr
    sSFR                 = SFR_yr / M_sun_count;      // 5e-4 yr^-1
    // t_cross: g_base*(1 + sSFR*t) = a_wind_t0*(1 + t/(t_age/YEAR_TO_S))
    // => (g_base*sSFR - a_wind_t0/t_age_yr)*t = a_wind_t0 - g_base
    double t_age_yr = t_age / YEAR_TO_S;
    double num      = a_wind_t0 - g_base_cache;
    double den      = g_base_cache * sSFR - a_wind_t0 / t_age_yr;
    t_cross         = (std::abs(den) > 1e-60) ? num / den : 0.0;   // ~67730 yr
    double g_sfr_tage = g_base_cache * m_factor_tage;
    binding_ratio_tage = g_sfr_tage / a_wind_tage;   // 2.654
}

inline double OrionUQFFModule::computeHz() {
    double Hz_kms = H0_kms * std::sqrt(Omega_m * std::pow(1.0 + z_redshift, 3.0) + Omega_L);
    return (Hz_kms * 1e3) / Mpc_to_m;
}

inline double OrionUQFFModule::computeMsfFactor(double t_s) {
    double t_yr      = t_s / YEAR_TO_S;
    double M_sun_cnt = M / M_SUN;
    return (SFR_yr * t_yr) / M_sun_cnt;
}

inline double OrionUQFFModule::computeUgSum() {
    double Ug1 = G_NEWTON * M / (r * r);
    double Ug4 = Ug1 * f_sc;
    return Ug1 + Ug4;
}

inline double OrionUQFFModule::computeQuantumTerm() {
    double Delta_p = HBAR / Delta_x;
    double unc     = std::sqrt(Delta_x * Delta_p);
    return (HBAR / unc) * (2.0 * PI / T_HUBBLE);
}

inline double OrionUQFFModule::computeEMTerm() {
    // MHD Alfven pressure gradient: B^2/(2*mu0*rho*r) [PAPER_318 EM companion]
    return WOLFRAM_TERM_ORION_TRAPEZIUM_RAD(
        (B * B) / (2.0 * MU_0 * rho_fluid * r)
    );
}

inline double OrionUQFFModule::computeFluidTerm() {
    // Gas self-gravity of nebular volume
    double M_fluid = rho_fluid * V_sys;
    return G_NEWTON * M_fluid / (r * r);
}

inline double OrionUQFFModule::computeResonantTerm(double t) {
    // PAPER_288 canonical: standing + traveling; omega_osc = v_wind/r (PDR mode)
    double cos_term = 2.0 * A_osc * std::cos(k_osc * 0.0) * std::cos(omega_osc * t);
    double travel   = (2.0 * PI / 13.8) * A_osc * std::cos(k_osc * 0.0 - omega_osc * t);
    return omega_osc * omega_osc * (cos_term + travel);
}

inline double OrionUQFFModule::computeDMTerm() {
    // DM acceleration + density perturbation
    double a_DM  = G_NEWTON * M_DM / (r * r);
    double pert  = delta_rho / rho_fluid;             // 0.1
    return a_DM * (1.0 + pert);
}

inline double OrionUQFFModule::computeWindTerm(double t) {
    // Ram pressure acceleration: v_wind^2/r * (1+t/t_age) [PAPER_317]
    return WOLFRAM_TERM_ORION_WIND_RAM(
        (v_wind * v_wind / r) * (1.0 + t / t_age)
    );
}

inline double OrionUQFFModule::computeRadTerm() {
    // Trapezium UV radiation acceleration [PAPER_318]
    double A_trap = 4.0 * PI * r * r;
    double P_rad  = L_trap / (A_trap * C_LIGHT);
    return WOLFRAM_TERM_ORION_TRAPEZIUM_RAD(P_rad / rho_fluid);
}

inline double OrionUQFFModule::computeG(double t) {
    double Hz        = computeHz();
    double expansion = 1.0 + Hz * t;
    double sc_corr   = 1.0 - (B / B_crit);
    double tr_factor = 1.0 + f_TRZ;
    double msf       = computeMsfFactor(t);
    double m_factor  = 1.0 + msf;

    // Base gravity — SFR+Hz+SC+TRZ [PAPER_319 SFR binding amplification]
    double g_base = WOLFRAM_TERM_ORION_BASE(
        (G_NEWTON * M * m_factor / (r * r)) * expansion * sc_corr * tr_factor
    );
    double ug_sum       = computeUgSum();
    double lambda_term  = LAMBDA_CC * C_LIGHT * C_LIGHT / 3.0;
    double quantum_term = computeQuantumTerm();
    double em_term      = computeEMTerm();
    double fluid_term   = computeFluidTerm();
    double resonant_term= computeResonantTerm(t);
    double dm_term      = computeDMTerm();
    double wind_term    = WOLFRAM_TERM_ORION_SFR_BINDING(computeWindTerm(t));
    double rad_term     = computeRadTerm();

    // 3-tier buoyancy — Sgr A* ext body at r_GC=2.623e20 m
    double ug1_base = G_NEWTON * M / (r * r);
    double tier1    =  0.5 * ug1_base;
    double tier2    = -beta_i * ug1_base * omega_g * (M   / r    ) * U_UA * std::cos(PI * t);
    double tier3    = -beta_i * ug1_base * omega_g * (M_GC / r_GC) * U_UA * std::cos(PI * t);

    double g_total = g_base + ug_sum + lambda_term + quantum_term
                   + em_term + fluid_term + resonant_term + dm_term
                   + wind_term + rad_term
                   + tier1 + tier2 + tier3;

    if (logging_enabled) {
        std::cout << std::scientific << std::setprecision(4)
                  << "[ORION LOG] t=" << t << " s\n"
                  << "  g_base="    << g_base      << "  ug_sum="  << ug_sum     << "\n"
                  << "  lambda="    << lambda_term  << "  quantum=" << quantum_term << "\n"
                  << "  em_MHD="    << em_term      << "  fluid="   << fluid_term << "\n"
                  << "  resonant="  << resonant_term<< "  dm="      << dm_term    << "\n"
                  << "  wind="      << wind_term    << "  rad="      << rad_term   << "\n"
                  << "  tier1="     << tier1        << "  tier2="   << tier2      << "\n"
                  << "  tier3="     << tier3        << "  TOTAL="   << g_total    << "\n";
    }
    return g_total;
}

inline std::string OrionUQFFModule::getEquationText() {
    return
        "g_Orion(r,t) = 12-term MUGE + 3-tier buoyancy [Session 91 UQFF 2.0]\n"
        "WOLFRAM x4: ORION_BASE / ORION_WIND_RAM [P317] / ORION_TRAPEZIUM_RAD [P318] / ORION_SFR_BINDING [P319]\n"
        "g_base = (G*M*(1+M_sf(t))/r^2)*(1+Hz*t)*(1-B/B_crit)*(1+f_TRZ)  [P319 SFR amplification]\n"
        "Ug_sum = Ug1 + Ug4 = 2*G*M/r^2\n"
        "Lambda = Lambda_cc*c^2/3\n"
        "quantum = (hbar/sqrt(hbar)) * 2pi/t_Hubble\n"
        "EM_MHD = B^2/(2*mu0*rho*r)  [Alfven pressure gradient]\n"
        "fluid = G*M_fluid/r^2  [gas self-gravity]\n"
        "resonant = omega_PDR^2 * A_osc * [standing+traveling]  [PDR eigenmode]\n"
        "DM = G*M_DM/r^2 * (1+delta_rho/rho)\n"
        "wind = v_wind^2/r * (1+t/t_age)  [P317 ram pressure; eta_wind=28.47]\n"
        "rad = L_trap/(4pi*r^2*c*rho)  [P318 Trapezium UV; eta_rad=7.664e18]\n"
        "tier1 = 0.5*Ug1  [CP3/PAPER_198 static]\n"
        "tier2 = -beta_i*Ug1*omega_g*(M/r)*U_UA*cos(pi*t)  [PAPER_198 dynamic]\n"
        "tier3 = -beta_i*Ug1*omega_g*(M_GC/r_GC)*U_UA*cos(pi*t)  [Sgr A* outer-frame]\n"
        "M_sf(t) = SFR_yr*t_yr / M_sun_count  [at t_age=150; at 1Myr=500]\n"
        "omega_PDR = v_wind/r = 6.780e-14 rad/s  [PDR eigenfreq; T_PDR=2.94 Myr]\n"
        "PAPER_317: eta_wind=28.47; t_erosion=467 kyr  [wind > gravity from birth]\n"
        "PAPER_318: eta_rad=7.664e18; champagne flow condition satisfied\n"
        "PAPER_319: t_cross=67730 yr; sSFR=5e-4 yr-1=50x Lagoon; unbound->bound phase transition";
}

inline void OrionUQFFModule::printParameters() {
    std::cout << std::scientific << std::setprecision(4)
        << "=== OrionUQFFModule Parameters [Session 91 UQFF 2.0] ===\n"
        << "System: Orion Nebula M42/NGC1976 — 33rd C++ module, FIRST Trapezium OB HII region\n"
        << "M           = " << M          << " kg (2000 M_sun)\n"
        << "r           = " << r          << " m (~12.5 ly)\n"
        << "V_sys       = " << V_sys      << " m^3\n"
        << "rho_fluid   = " << rho_fluid  << " kg/m^3\n"
        << "v_wind      = " << v_wind     << " m/s (8 km/s HII front)\n"
        << "t_age       = " << t_age      << " s (3e5 yr)\n"
        << "SFR_yr      = " << SFR_yr     << " M_sun/yr\n"
        << "L_trap      = " << L_trap     << " W (Trapezium theta1 Ori C, 2e5 L_sun)\n"
        << "B           = " << B          << " T\n"
        << "f_TRZ       = " << f_TRZ      << "\n"
        << "z           = " << z_redshift << "\n"
        << "omega_PDR   = " << omega_osc  << " rad/s  [v_wind/r, PDR eigenfreq]\n"
        << "M_GC(SgrA*) = " << M_GC       << " kg (ext body, ~8.5 kpc)\n"
        << "r_GC        = " << r_GC       << " m\n"
        << "--- Cached PAPER_317 ---\n"
        << "g_base      = " << g_base_cache   << " m/s^2\n"
        << "a_wind(t=0) = " << a_wind_t0      << " m/s^2\n"
        << "eta_wind    = " << eta_wind        << " (P_ram/P_grav)\n"
        << "a_wind(tage)= " << a_wind_tage     << " m/s^2\n"
        << "eta_wind_age= " << eta_wind_tage   << "\n"
        << "t_erosion   = " << t_erosion       << " s (" << t_erosion / YEAR_TO_S / 1e3 << " kyr)\n"
        << "--- Cached PAPER_318 ---\n"
        << "a_rad       = " << a_rad_cache     << " m/s^2\n"
        << "eta_rad     = " << eta_rad         << " (radiation/gravity)\n"
        << "--- Cached PAPER_319 ---\n"
        << "m_factor(tage)=" << m_factor_tage  << " (SFR amplification at 300 kyr)\n"
        << "m_factor(1Myr)=" << m_factor_1myr  << "\n"
        << "t_consume   = " << t_consume       << " yr (gas depletion)\n"
        << "sSFR        = " << sSFR            << " yr^-1 (50x Lagoon)\n"
        << "t_cross     = " << t_cross         << " yr (unbound->bound crossover)\n"
        << "bind_ratio  = " << binding_ratio_tage << " at t_age\n";
}

inline void OrionUQFFModule::updateVariable(const std::string& name, double value) {
    if      (name == "M")         { M = value; M_visible = 0.15*M; M_DM = 0.85*M; }
    else if (name == "r")         { r = value; V_sys = (4.0/3.0)*PI*r*r*r; k_osc = 2.0*PI/r; omega_osc = v_wind/r; }
    else if (name == "rho_fluid") { rho_fluid = value; delta_rho = 0.1*value; }
    else if (name == "v_wind")    { v_wind = value; omega_osc = v_wind/r; }
    else if (name == "t_age")     { t_age = value; }
    else if (name == "B")         { B = value; }
    else if (name == "z")         { z_redshift = value; }
    else if (name == "f_TRZ")     { f_TRZ = value; }
    else                          { dynamic_params[name] = value; }
    updateCache();
}

inline void OrionUQFFModule::addToVariable(const std::string& name, double delta) {
    if      (name == "M")         { M += delta; M_visible = 0.15*M; M_DM = 0.85*M; }
    else if (name == "r")         { r += delta; V_sys = (4.0/3.0)*PI*r*r*r; }
    else if (name == "rho_fluid") { rho_fluid += delta; delta_rho = 0.1*rho_fluid; }
    else if (name == "v_wind")    { v_wind += delta; omega_osc = v_wind/r; }
    else {
        double cur = getDynamicParameter(name, 0.0);
        dynamic_params[name] = cur + delta;
    }
    updateCache();
}

inline void OrionUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

inline void OrionUQFFModule::exportState(const std::string& filename) {
    std::ofstream f(filename);
    if (!f.is_open()) { std::cerr << "[ORION] Cannot open " << filename << "\n"; return; }
    f << std::scientific << std::setprecision(6);
    f << "# ORION_UQFF_MODULE State Export — Session 91\n";
    f << "# System: Orion Nebula M42/NGC1976 — 33rd C++ module\n";
    // System params
    f << "M," << M << "\n";
    f << "r," << r << "\n";
    f << "V_sys," << V_sys << "\n";
    f << "rho_fluid," << rho_fluid << "\n";
    f << "v_wind," << v_wind << "\n";
    f << "t_age," << t_age << "\n";
    f << "SFR_yr," << SFR_yr << "\n";
    f << "L_trap," << L_trap << "\n";
    f << "B," << B << "\n";
    f << "B_crit," << B_crit << "\n";
    f << "f_TRZ," << f_TRZ << "\n";
    f << "f_sc," << f_sc << "\n";
    f << "z_redshift," << z_redshift << "\n";
    f << "H0_kms," << H0_kms << "\n";
    f << "omega_PDR," << omega_osc << "\n";
    f << "k_PDR," << k_osc << "\n";
    f << "A_osc," << A_osc << "\n";
    f << "M_GC," << M_GC << "\n";
    f << "r_GC," << r_GC << "\n";
    f << "beta_i," << beta_i << "\n";
    f << "omega_g," << omega_g << "\n";
    f << "U_UA," << U_UA << "\n";
    f << "M_visible," << M_visible << "\n";
    f << "M_DM," << M_DM << "\n";
    f << "delta_rho," << delta_rho << "\n";
    // Cached P317
    f << "g_base_cache," << g_base_cache << "\n";
    f << "a_wind_t0," << a_wind_t0 << "\n";
    f << "eta_wind," << eta_wind << "\n";
    f << "a_wind_tage," << a_wind_tage << "\n";
    f << "eta_wind_tage," << eta_wind_tage << "\n";
    f << "t_erosion_s," << t_erosion << "\n";
    // Cached P318
    f << "a_rad_cache," << a_rad_cache << "\n";
    f << "eta_rad," << eta_rad << "\n";
    // Cached P319
    f << "m_factor_tage," << m_factor_tage << "\n";
    f << "m_factor_1myr," << m_factor_1myr << "\n";
    f << "t_consume_yr," << t_consume << "\n";
    f << "sSFR," << sSFR << "\n";
    f << "t_cross_yr," << t_cross << "\n";
    f << "binding_ratio_tage," << binding_ratio_tage << "\n";
    f.close();
    std::cout << "[ORION] State exported to " << filename << "\n";
}

#endif // ORION_UQFF_MODULE_H

// Constructor: Set all variables with Orion Nebula-specific values
OrionUQFFModule::OrionUQFFModule() {
    // Base constants (universal)
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["c"] = 3e8;                           // m/s
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Lambda"] = 1.1e-52;                  // m^-2
    variables["q"] = 1.602e-19;                     // C
    variables["pi"] = 3.141592653589793;            // pi
    variables["t_Hubble"] = 13.8e9 * 3.156e7;       // s
    variables["year_to_s"] = 3.156e7;               // s/yr

    // Orion Nebula parameters
    double M_sun_val = 1.989e30;                    // kg
    variables["M_sun"] = M_sun_val;
    variables["M"] = 2000 * M_sun_val;              // Total mass kg
    variables["M0"] = variables["M"];               // Initial mass
    variables["SFR"] = 1 * M_sun_val;               // Msun/yr
    variables["M_visible"] = 0.15 * variables["M"]; // Visible fraction est.
    variables["M_DM"] = 0.85 * variables["M"];      // Dark matter/halo
    variables["r"] = 1.18e17;                       // m (half span ~12.5 ly)

    // Hubble/cosmology
    variables["H0"] = 67.15;                        // km/s/Mpc
    variables["Mpc_to_m"] = 3.086e22;               // m/Mpc
    variables["z"] = 0.00034;                       // Redshift
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["t"] = 1e6 * variables["year_to_s"];  // Default t=1 Myr s

    // Gas/wind dynamics
    variables["rho_fluid"] = 1e-20;                 // kg/m^3 (dense gas)
    variables["V"] = 1e3;                           // m^3 (arbitrary)
    variables["v_wind"] = 8e3;                      // m/s (8 km/s)
    variables["t_age"] = 3e5 * variables["year_to_s"];  // s (~300k yr)
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];

    // EM/magnetic
    variables["B"] = 1e-5;                          // T (nebula field)
    variables["B_crit"] = 1e11;                     // T

    // Quantum terms
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;

    // Resonant/oscillatory
    variables["A"] = 1e-10;
    variables["k"] = 1e20;
    variables["omega"] = 1e15;                      // rad/s
    variables["x"] = 0.0;

    // Ug subterms
    variables["Ug1"] = 0.0;
    variables["Ug2"] = 0.0;
    variables["Ug3"] = 0.0;
    variables["Ug4"] = 0.0;

    // Scale factors
    variables["scale_macro"] = 1e-12;
    variables["f_TRZ"] = 0.1;
    variables["f_sc"] = 1.0;
}

// Update variable (set to new value)
void OrionUQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "M") {
        variables["M_visible"] = 0.15 * value;
        variables["M_DM"] = 0.85 * value;
        variables["M0"] = value;
    }
}

// Add delta to variable
void OrionUQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void OrionUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(z) in s^-1
double OrionUQFFModule::computeHz() {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + variables["z"], 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// Compute Ug sum: Ug1 = G M / r^2, Ug4 = Ug1 * f_sc, others 0
double OrionUQFFModule::computeUgSum() {
    double Ug1 = (variables["G"] * variables["M"]) / (variables["r"] * variables["r"]);
    variables["Ug1"] = Ug1;
    variables["Ug4"] = Ug1 * variables["f_sc"];
    return variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
double OrionUQFFModule::computeQuantumTerm(double t_Hubble_val) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double integral_val = variables["integral_psi"];
    return (variables["hbar"] / unc) * integral_val * (2 * variables["pi"] / t_Hubble_val);
}

// Fluid term: rho_fluid * V * g
double OrionUQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
double OrionUQFFModule::computeResonantTerm(double t) {
    double cos_term = 2 * variables["A"] * std::cos(variables["k"] * variables["x"]) * std::cos(variables["omega"] * t);
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega"] * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);
    return cos_term + exp_factor * real_exp;
}

// DM term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3)
double OrionUQFFModule::computeDMTerm() {
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3 * variables["G"] * variables["M"] / (variables["r"] * variables["r"] * variables["r"]);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Star formation factor: (SFR * t_yr) / M0
double OrionUQFFModule::computeMsfFactor(double t) {
    double t_yr = t / variables["year_to_s"];
    return (variables["SFR"] * t_yr) / variables["M0"];
}

// Stellar wind term: W_stellar = rho * v_wind^2 * (1 + t / t_age)
double OrionUQFFModule::computeW_stellar(double t) {
    return variables["rho_fluid"] * std::pow(variables["v_wind"], 2) * (1.0 + t / variables["t_age"]);
}

// Full computation: g_UQFF(r, t) = ... all terms with M_sf and + W_stellar
double OrionUQFFModule::computeG(double t) {
    variables["t"] = t;
    double Hz = computeHz();
    double expansion = 1.0 + Hz * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double tr_factor = 1.0 + variables["f_TRZ"];
    double msf_factor = computeMsfFactor(t);
    double m_factor = 1.0 + msf_factor;
    double w_stellar = computeW_stellar(t);

    // Base gravity with expansion, SC, TR, M_sf
    double g_base = (variables["G"] * variables["M"] * m_factor / (variables["r"] * variables["r"])) * expansion * sc_correction * tr_factor;

    // Ug sum
    double ug_sum = computeUgSum();

    // Cosmological
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"]) / 3.0;

    // Quantum
    double quantum_term = computeQuantumTerm(variables["t_Hubble"]);

    // EM Lorentz (magnitude v_wind B)
    double em_base = variables["q"] * variables["v_wind"] * variables["B"] / 1.673e-27;
    double em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * variables["scale_macro"];

    // Fluid
    double fluid_term = computeFluidTerm(g_base);

    // Resonant
    double resonant_term = computeResonantTerm(t);

    // DM
    double dm_term = computeDMTerm();

    // Total: Sum all + W_stellar
    return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + w_stellar;
}

// Get equation text (descriptive)
std::string OrionUQFFModule::getEquationText() {
    return "g_Orion(r, t) = (G * M(t) / r^2) * (1 + H(z) * t) * (1 - B / B_crit) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 / 3) + "
           "(hbar / sqrt(Delta_x * Delta_p)) * ∫(ψ* H ψ dV) * (2π / t_Hubble) + q (v × B) + ρ_fluid * V * g + "
           "2 A cos(k x) cos(ω t) + (2π / 13.8) A exp(i (k x - ω t)) + (M_visible + M_DM) * (δρ/ρ + 3 G M / r^3) + W_stellar\n"
           "Where M(t) = M * (1 + M_sf(t)); M_sf(t) = (SFR * t_yr) / M0; W_stellar = ρ * v_wind^2 * (1 + t / t_age)\n"
           "Special Terms:\n"
           "- Quantum: Heisenberg uncertainty for gas quantum effects.\n"
           "- Fluid: Nebular gas density-volume-gravity coupling.\n"
           "- Resonant: Oscillatory Aether waves for proplyds.\n"
           "- DM: Visible+dark mass with perturbations for halo.\n"
           "- Superconductivity: (1 - B/B_crit) for quantum fields.\n"
           "- Star Formation: M_sf(t) boosts mass via SFR=1 Msun/yr.\n"
           "- Stellar Wind: W_stellar from Trapezium erodes pillars.\n"
           "Solutions: At t=1 Myr, g_Orion ~1e-11 m/s² (W_stellar/fluid dominant; g_base ~1e-12).\n"
           "Adaptations for Orion Nebula: Trapezium radiation/winds; z=0.00034; SFR for starbirth.";
}

// Print variables
void OrionUQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}