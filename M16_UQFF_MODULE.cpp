// M16_UQFF_MODULE.cpp
// UQFF 2.0 — M16 Eagle Nebula Master Gravity Module (22nd C++ module)
// FIRST UQFF NEBULAR MODULE WITH z>0 (z=0.0015, ~5700 ly Eagle Nebula cosmological redshift).
// M16: M=1200 M_sun=2.387e33 kg, r=3.31e17 m (~35 ly half-span), z=0.0015, v_gas=1e5 m/s
// g_base = G*M/r² = 1.454e-12 m/s² (nebular scale, same order as galactic modules)
// pre_sum_Ug = 52 × g_base = 7.56e-11 m/s²
// PAPER_284: Phi_dm = (1+SFR_rate*t)*(1-E_0*(1-exp(-t/tau))); gap_mult_add = -(M_sf_frac*E_rad)
//            at t=5 Myr: M_sf_frac=4164.8; E_rad=0.2433; Phi_dm=3151.9; gap=-1013.3 (24.3% less than additive)
// PAPER_285: t_half = tau*ln(2) = 6.561e13 s = 2.079 Myr; DeltaGMax = E_0*g_base = 4.36e-13 m/s²
// PAPER_286: H(z=0.0015)=70.047 km/s/Mpc; kappa_neb=[H(z)-H(z=0)]/H(z=0)=6.71e-4; first UQFF kappa_neb
// Copyright - Daniel T. Murphy, UQFF 2.0 upgrade Session 80, March 2026.

#ifndef M16_UQFF_MODULE_H
#define M16_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <cassert>
#include <random>
#include <chrono>

// UQFF 2.0 — Wolfram Term Macros (auto-registration with Wolfram KB)
#define WOLFRAM_TERM_M16_BASE     "M16UQFF:g_total=[g_dyn(t)+Ug_sum(26)+Lambda+quantum+Lorentz+fluid+g_exp]*corr_SC; z=0.0015"
#define WOLFRAM_TERM_M16_SFR_ERO  "M16UQFF:Phi_dm=(1+SFR_rate*t)*(1-E_0*(1-Exp[-t/tau])); SFR_rate=2.639e-11/s; M_sf_frac=SFR_rate*t [PAPER_284]"
#define WOLFRAM_TERM_M16_ERODE    "M16UQFF:t_half=tau*Log[2]=6.561e13s=2.079Myr; DeltaGMax=E_0*g_base=4.36e-13 m/s^2 [PAPER_285]"
#define WOLFRAM_TERM_M16_FRIEDMANN "M16UQFF:kappa_neb=[H[z=0.0015]-H[z=0]]/H[z=0]=6.71e-4; H[z=0.0015]=70.047 km/s/Mpc [PAPER_286]"

class M16UQFFModule {
private:
    // =========================================================================
    // UQFF 2.0 Self-Expanding Framework
    // =========================================================================
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;
    mutable std::mt19937 rng;
    mutable std::uniform_real_distribution<double> dis;

    // =========================================================================
    // M16 Physics Parameters
    // =========================================================================

    // Nebula mass & geometry
    double M;               // Eagle Nebula gas mass (kg): 1200 M_sun = 2.387e33
    double r;               // Half-span radius (m): 3.31e17 (~35 ly)
    double M_Sun;           // Solar mass reference (kg): 1.989e30

    // Cosmological (PAPER_286 — FIRST nebular z>0 in UQFF)
    double z;               // Cosmological redshift: 0.0015 (~5700 ly Eagle Nebula)

    // Turbulent gas (Lorentz term)
    double v_gas;           // Turbulent gas velocity (m/s): 1e5

    // Star Formation + Erosion (PAPER_284 & PAPER_285)
    double SFR_Msun_yr;     // Star formation rate (M_sun/yr): 1.0
    double M0;              // Initial nebula mass (kg): 2.387e33 (= M)
    double SFR_rate;        // SFR_rate = SFR/(M0/M_sun) in s⁻¹ ≈ 2.639e-11 (computed in updateCache)
    double tau_erode;       // Photoevaporation e-folding time (s): 3 Myr = 9.468e13
    double E_0;             // Max photoevaporation fraction (dimensionless): 0.3

    // PAPER_285 cached diagnostics
    double t_half;          // Erosion half-time: tau*ln2 = 6.561e13 s = 2.079 Myr
    double DeltaGMax;       // Max erosion gravity amplitude: E_0*g_base (m/s²)

    // Magnetic / SC correction
    double B_field;         // Nebular magnetic field (T): 1e-5
    double B_crit;          // Magnetar critical field (T): 1e11
    double q_charge;        // Reference charge (C): 1.6e-19 (proton)
    double corr_SC;         // 1 - B_field/B_crit (near-unity for nebula)

    // Fluid (molecular cloud interface)
    double rho_fluid;       // Dense pillar molecular cloud density (kg/m³): 1e-20
    double V_fluid;         // Nebula volume (m³): (4/3)pi*r³ — computed in updateCache

    // Physical constants
    double G_grav;          // Newton's G (m³/kg/s²): 6.674e-11
    double Lambda;          // Cosmological constant (m⁻²): 1.114e-52
    double c_light;         // Speed of light (m/s): 2.998e8
    double hbar;            // Reduced Planck (J·s): 1.0546e-34
    double delta_x;         // Position uncertainty (m): 1e-10
    double delta_p;         // Momentum uncertainty (kg·m/s): hbar/delta_x
    double integral_psi;    // Wavefunction normalization: 1.0
    double t_Hubble;        // Hubble time (s): 4.352e17

    // Friedmann (PAPER_286 — first nebular z>0)
    double H0;              // Hubble constant (km/s/Mpc): 70.0
    double Omega_m;         // Matter density: 0.3
    double Omega_Lam;       // Dark energy: 0.7
    double Mpc_to_m;        // m/Mpc: 3.086e22
    double kappa_neb;       // [H(z=0.0015)-H(z=0)]/H(z=0) = 6.71e-4 (cached in updateCache)

    // No DM (nebular case: f_DM = 0, all visible gas mass)

    // 26-layer Triadic MUGE
    std::vector<double> Ug1_vec;  // Magnetic dipole (= g_base per layer)
    std::vector<double> Ug2_vec;  // Zero (negligible for nebula: Ug2/Ug3≈0)
    std::vector<double> Ug3_vec;  // Zero (negligible for nebula)
    std::vector<double> Ug4_vec;  // Vacuum concentration (= Ug1 per layer)

    // Precomputed caches
    double pre_sum_Ug;      // 26-layer Triadic sum: 52 × g_base
    double g_base_cache;    // G×M/r²

    // =========================================================================
    // Private helper implementations
    // =========================================================================
    double computeQuantumTerm(double t_H) {
        return (hbar / std::sqrt(delta_x * delta_p)) * integral_psi * (2.0 * M_PI / t_H);
    }

    double computeFluidTerm(double g_bc) {
        // Use g_bc (g_base_cache, NOT g_dyn) to avoid double-counting Phi_dm modulation
        double rho_mean = M / V_fluid;
        return (rho_fluid / rho_mean) * g_bc;
    }

    double computeUgSum() { return pre_sum_Ug; }

    double computeHz(double zz) {
        // Friedmann H(z) in s⁻¹
        double H_kms = H0 * std::sqrt(Omega_m * std::pow(1.0 + zz, 3.0) + Omega_Lam);
        return H_kms * 1.0e3 / Mpc_to_m;
    }

    // PAPER_284: Dual mass co-action product
    double computeMsfFactor(double t) {
        return SFR_rate * t;                        // dimensionless growth fraction
    }

    // PAPER_284+285: Photoevaporation saturation
    double computeErad(double t) {
        return E_0 * (1.0 - std::exp(-t / tau_erode));
    }

    void updateCache() {
        g_base_cache = G_grav * M / (r * r);
        // V_fluid = (4/3)pi*r³
        V_fluid = (4.0 / 3.0) * M_PI * r * r * r;
        // 26-layer Triadic
        for (int i = 0; i < 26; ++i) {
            Ug1_vec[i] = g_base_cache;
            Ug4_vec[i] = g_base_cache;
        }
        pre_sum_Ug = 0.0;
        for (int i = 0; i < 26; ++i)
            pre_sum_Ug += Ug1_vec[i] + Ug2_vec[i] + Ug3_vec[i] + Ug4_vec[i];
        // PAPER_284: SFR_rate = SFR_Msun_yr / (M0/M_Sun) / (3.156e7 s/yr)
        //            = 1 M_sun/yr / (1200 M_sun) / (3.156e7) = 2.639e-11 s⁻¹
        SFR_rate = SFR_Msun_yr / (M0 / M_Sun) / 3.156e7;
        // PAPER_285: cached erosion half-time and max amplitude
        t_half   = tau_erode * std::log(2.0);
        DeltaGMax = E_0 * g_base_cache;
        // PAPER_286: kappa_neb = [H(z)-H(z=0)]/H(z=0)
        double Hz0  = computeHz(0.0);
        double Hz   = computeHz(z);
        kappa_neb   = (Hz - Hz0) / Hz0;
        // SC correction
        corr_SC  = 1.0 - B_field / B_crit;
        delta_p  = hbar / delta_x;
    }

    void log(const std::string& msg) const {
        if (logging_enabled) std::cout << "[LOG][M16] " << msg << "\n";
    }

public:
    // =========================================================================
    // Constructor — UQFF 2.0 with M16 catalogue defaults
    // =========================================================================
    M16UQFFModule()
        : rng(static_cast<unsigned>(std::chrono::steady_clock::now().time_since_epoch().count()))
        , dis(0.0, 1.0)
    {
        logging_enabled = false;
        dynamic_params.clear();

        // Eagle Nebula parameters
        M_Sun        = 1.989e30;                  // kg Sun mass reference
        M            = 1200.0 * M_Sun;             // kg: 1200 M_sun = 2.387e33
        M0           = M;                          // kg initial mass = M
        r            = 3.31e17;                    // m: ~35 ly half-span

        // PAPER_286 — first nebular z>0 in UQFF
        z            = 0.0015;                     // cosmological redshift (Eagle Nebula ~5700 ly)

        // Turbulent gas (Lorentz coupling)
        v_gas        = 1.0e5;                      // m/s gas turbulent velocity

        // PAPER_284 SFR / PAPER_285 erosion
        SFR_Msun_yr  = 1.0;                        // M_sun/yr star formation rate
        tau_erode    = 3.0e6 * 3.156e7;            // s: 3 Myr photoevaporation e-folding
        E_0          = 0.3;                        // 30% max photoevaporation fraction

        // Magnetic field (pillar region ~10 μG)
        B_field      = 1.0e-5;                     // T (10 μG nebular magnetic field)
        B_crit       = 1.0e11;                     // T magnetar critical reference
        q_charge     = 1.6e-19;                    // C proton reference charge

        // Fluid — dense pillar molecular cloud interface
        rho_fluid    = 1.0e-20;                    // kg/m³ dense pillar interface

        G_grav       = 6.674e-11;
        Lambda       = 1.114e-52;
        c_light      = 2.998e8;
        hbar         = 1.0546e-34;
        delta_x      = 1.0e-10;
        integral_psi = 1.0;
        t_Hubble     = 4.352e17;

        // Friedmann (PAPER_286 — z=0.0015)
        H0        = 70.0;
        Omega_m   = 0.3;
        Omega_Lam = 0.7;
        Mpc_to_m  = 3.086e22;

        // Initialise 26-layer vectors (Ug2/Ug3 = 0 for nebula)
        Ug1_vec.resize(26, 0.0);
        Ug2_vec.resize(26, 0.0);
        Ug3_vec.resize(26, 0.0);
        Ug4_vec.resize(26, 0.0);

        // Zero-init caches
        pre_sum_Ug  = 0.0;
        g_base_cache = 0.0;
        V_fluid     = 0.0;
        SFR_rate    = 0.0;
        t_half      = 0.0;
        DeltaGMax   = 0.0;
        kappa_neb   = 0.0;
        corr_SC     = 1.0;
        delta_p     = hbar / delta_x;

        updateCache();
    }

    ~M16UQFFModule() {}

    // =========================================================================
    // Dynamic variable operations
    // =========================================================================
    void updateVariable(const std::string& name, double value) {
        dynamic_params[name] = value;
        if      (name == "M")           { M           = value; updateCache(); }
        else if (name == "r")           { r           = value; updateCache(); }
        else if (name == "z")           { z           = value; updateCache(); }
        else if (name == "v_gas")       { v_gas       = value; }
        else if (name == "SFR_Msun_yr") { SFR_Msun_yr = value; updateCache(); }
        else if (name == "M0")          { M0          = value; updateCache(); }
        else if (name == "tau_erode")   { tau_erode   = value; updateCache(); }
        else if (name == "E_0")         { E_0         = value; updateCache(); }
        else if (name == "B_field")     { B_field     = value; updateCache(); }
        else if (name == "rho_fluid")   { rho_fluid   = value; }
        else if (name == "H0")          { H0          = value; updateCache(); }
        else if (name == "Omega_m")     { Omega_m     = value; updateCache(); }
        else if (name == "Omega_Lam")   { Omega_Lam   = value; updateCache(); }
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
    // Core: Full g_UQFF(r, t) for M16 Eagle Nebula
    // =========================================================================
    double computeG(double t) {
        auto _t0 = std::chrono::high_resolution_clock::now();

        // ---- PAPER_284: Dual mass co-action product (SFR growth × erosion saturation) ----
        double M_sf_frac     = computeMsfFactor(t);               // SFR_rate * t (unitless)
        double E_rad         = computeErad(t);                    // E_0*(1-exp(-t/tau))
        double Phi_dm        = (1.0 + M_sf_frac) * (1.0 - E_rad); // MULTIPLICATIVE coupling
        double gap_mult_add  = -(M_sf_frac * E_rad);              // PAPER_284: vs additive form

        double g_dyn = g_base_cache * Phi_dm;                     // Dynamic gravity: base × Phi_dm

        // ---- 26-layer Triadic (static: uses g_base, NOT g_dyn — avoids double-counting) ----
        double Ug_sum  = computeUgSum();                          // 52 × g_base

        // ---- Cosmological Λ ----
        double Lambda_tm = (Lambda * c_light * c_light) / 3.0;

        // ---- Quantum HUP ----
        double quantum = computeQuantumTerm(t_Hubble);

        // ---- Lorentz: q*v_gas*B_field/M (turbulent gas × nebular B-field coupling) ----
        double Lorentz = (q_charge * v_gas * B_field) / M;

        // ---- Fluid buoyancy (uses g_base_cache, NOT g_dyn — independent term) ----
        double fluid = computeFluidTerm(g_base_cache);

        // ---- PAPER_286: Hubble expansion with nebular Friedmann z=0.0015 ----
        double H_z   = computeHz(z);                              // H(z=0.0015) in s⁻¹
        double g_exp = g_base_cache * H_z * t;                   // Friedmann expansion coupling

        // ---- Total sum ----
        double g_sum   = g_dyn + Ug_sum + Lambda_tm + quantum + Lorentz + fluid + g_exp;
        double g_total = g_sum * corr_SC;

        if (logging_enabled) {
            auto _t1 = std::chrono::high_resolution_clock::now();
            double _ms = std::chrono::duration<double, std::milli>(_t1 - _t0).count();
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(4)
                << "computeG: t=" << t
                << " M_sf_frac=" << M_sf_frac << " E_rad=" << E_rad
                << " Phi_dm=" << Phi_dm << " gap_mult_add=" << gap_mult_add
                << " g_dyn=" << g_dyn << " Ug_sum=" << Ug_sum
                << " Lambda=" << Lambda_tm << " quantum=" << quantum
                << " Lorentz=" << Lorentz << " fluid=" << fluid
                << " g_exp=" << g_exp << " H_z=" << H_z
                << " kappa_neb=" << kappa_neb << " corr_SC=" << corr_SC
                << " t_half=" << t_half << " DeltaGMax=" << DeltaGMax
                << " g_total=" << g_total << " elapsed=" << _ms << "ms";
            log(oss.str());
        }
        return g_total;
    }

    // =========================================================================
    // Equation text output
    // =========================================================================
    std::string getEquationText() {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(6);
        oss << "M16 Eagle Nebula UQFF 2.0 (22nd C++ module — FIRST nebular z>0):\n"
            << "  g_total(r,t) = [g_dyn(t) + Ug_sum(26) + Lambda + quantum"
               " + Lorentz + fluid + g_exp] * corr_SC\n\n"
            << "  PAPER_284 — Dual Mass Co-Action Product (SFR*Erosion MULTIPLICATIVE):\n"
            << "    g_dyn(t)   = g_base * Phi_dm(t)\n"
            << "    Phi_dm(t)  = (1 + SFR_rate*t) * (1 - E_0*(1-exp(-t/tau)))\n"
            << "    k/a SFR growth [×] erosion saturation — FIRST UQFF multiplicative coupling\n"
            << "    gap_mult_add = -(M_sf_frac * E_rad) at t [24.3% below additive at 5 Myr]\n"
            << "    SFR_rate   = " << SFR_rate << " s^-1  (1 M_sun/yr / 1200 M_sun / 3.156e7)\n"
            << "    g_base     = " << g_base_cache << " m/s^2  (G*M/r^2)\n\n"
            << "  PAPER_285 — Erosion Saturation Half-Time:\n"
            << "    t_half     = tau*ln(2) = " << t_half << " s = 2.079 Myr\n"
            << "    DeltaGMax  = E_0 * g_base = " << DeltaGMax << " m/s^2\n"
            << "    tau_erode  = " << tau_erode << " s (3 Myr)\n"
            << "    E_0        = " << E_0 << " (30% max photoevaporation)\n\n"
            << "  PAPER_286 — Nebular Friedmann Redshift (FIRST UQFF kappa_neb):\n"
            << "    z          = " << z << " (Eagle Nebula ~5700 ly, FIRST UQFF nebular z>0)\n"
            << "    H(z=0.0015)= " << (computeHz(z) * Mpc_to_m / 1.0e3) << " km/s/Mpc\n"
            << "    kappa_neb  = [H(z=0.0015)-H(z=0)]/H(z=0) = " << kappa_neb << "\n"
            << "    g_exp(t)   = g_base * H(z=0.0015) * t  [Friedmann expansion coupling]\n\n"
            << "  Ug_sum(26)   = 52 * g_base = " << pre_sum_Ug << " m/s^2  [26-layer Triadic]\n"
            << "  Lorentz      = q * v_gas * B_field / M  [turbulent gas B-field coupling]\n"
            << "    v_gas      = " << v_gas << " m/s (turbulent)\n"
            << "    B_field    = " << B_field << " T (nebular 10 uG)\n"
            << "  fluid        = (rho_fluid/rho_mean) * g_base  [pillar molecular cloud interface]\n"
            << "    rho_fluid  = " << rho_fluid << " kg/m^3 (dense pillar)\n"
            << "    V_fluid    = " << V_fluid << " m^3 ((4/3)pi*r^3)\n"
            << "  corr_SC      = 1 - B/B_crit = " << corr_SC
            << "  [B=" << B_field << " T  B_crit=" << B_crit << " T]\n"
            << "  M            = " << M << " kg (1200 M_sun Eagle Nebula gas mass)\n"
            << "  r            = " << r << " m (~35 ly half-span)\n";
        return oss.str();
    }

    // =========================================================================
    // Print variables
    // =========================================================================
    void printVariables() {
        std::cout << std::scientific << std::setprecision(4);
        std::cout << "=== M16 Eagle Nebula UQFF Module Variables (UQFF 2.0) ===" << std::endl;
        std::cout << "M              : " << M             << " kg (1200 M_sun Eagle Nebula gas mass)" << std::endl;
        std::cout << "M0             : " << M0            << " kg (initial mass = M)" << std::endl;
        std::cout << "r              : " << r             << " m (~35 ly half-span)" << std::endl;
        std::cout << "g_base         : " << g_base_cache  << " m/s^2 (G*M/r^2 = 1.454e-12 m/s^2)" << std::endl;
        std::cout << "pre_sum_Ug     : " << pre_sum_Ug    << " m/s^2 (52*g_base)" << std::endl;
        std::cout << "z              : " << z             << " (FIRST UQFF nebular z>0 [PAPER_286])" << std::endl;
        std::cout << "v_gas          : " << v_gas         << " m/s (turbulent gas)" << std::endl;
        std::cout << "[PAPER_284] SFR_Msun_yr : " << SFR_Msun_yr << " M_sun/yr" << std::endl;
        std::cout << "[PAPER_284] SFR_rate    : " << SFR_rate    << " s^-1 (2.639e-11)" << std::endl;
        std::cout << "[PAPER_285] tau_erode   : " << tau_erode   << " s (3 Myr)" << std::endl;
        std::cout << "[PAPER_285] E_0         : " << E_0         << " (30% max photoevaporation)" << std::endl;
        std::cout << "[PAPER_285] t_half      : " << t_half      << " s = 2.079 Myr" << std::endl;
        std::cout << "[PAPER_285] DeltaGMax   : " << DeltaGMax   << " m/s^2 = E_0*g_base" << std::endl;
        std::cout << "[PAPER_286] kappa_neb   : " << kappa_neb   << " = [H(z)-H(0)]/H(0)" << std::endl;
        std::cout << "B_field        : " << B_field       << " T (nebular 10 uG)" << std::endl;
        std::cout << "B_crit         : " << B_crit        << " T (magnetar critical reference)" << std::endl;
        std::cout << "corr_SC        : " << corr_SC       << " (1 - B/B_crit, near-unity for nebula)" << std::endl;
        std::cout << "rho_fluid      : " << rho_fluid     << " kg/m^3 (dense pillar molecular cloud)" << std::endl;
        std::cout << "V_fluid        : " << V_fluid       << " m^3 ((4/3)*pi*r^3)" << std::endl;
        std::cout << "H0             : " << H0            << " km/s/Mpc" << std::endl;
        std::cout << "Lambda         : " << Lambda        << " m^-2" << std::endl;
        std::cout << "logging        : " << (logging_enabled ? "ON" : "OFF") << std::endl;
        std::cout << "dynamic_params : " << dynamic_params.size() << " entries" << std::endl;
    }

    // =========================================================================
    // UQFF 2.0 Self-Expanding Framework Methods
    // =========================================================================
    void setEnableLogging(bool enabled)  { logging_enabled = enabled; }
    bool getLoggingEnabled() const       { return logging_enabled; }

    void setDynamicParameter(const std::string& name, double value) {
        dynamic_params[name] = value;
        updateCache();
    }

    double getDynamicParameter(const std::string& name) const {
        auto it = dynamic_params.find(name);
        return (it != dynamic_params.end()) ? it->second : 0.0;
    }

    void setScalingFactor(const std::string& name, double value) {
        setDynamicParameter(name, value);
    }

    double getScalingFactor(const std::string& name) const {
        return getDynamicParameter(name);
    }

    void exportState(const std::string& filename = "M16UQFF_state.txt") const {
        std::ofstream ofs(filename);
        if (!ofs.is_open()) return;
        ofs << std::scientific;
        ofs << "M "              << M              << "\n";
        ofs << "M0 "             << M0             << "\n";
        ofs << "r "              << r              << "\n";
        ofs << "M_Sun "          << M_Sun          << "\n";
        ofs << "z "              << z              << "\n";
        ofs << "v_gas "          << v_gas          << "\n";
        ofs << "SFR_Msun_yr "    << SFR_Msun_yr    << "\n";
        ofs << "SFR_rate "       << SFR_rate       << "\n";
        ofs << "tau_erode "      << tau_erode      << "\n";
        ofs << "E_0 "            << E_0            << "\n";
        ofs << "t_half "         << t_half         << "\n";
        ofs << "DeltaGMax "      << DeltaGMax      << "\n";
        ofs << "B_field "        << B_field        << "\n";
        ofs << "B_crit "         << B_crit         << "\n";
        ofs << "corr_SC "        << corr_SC        << "\n";
        ofs << "q_charge "       << q_charge       << "\n";
        ofs << "rho_fluid "      << rho_fluid      << "\n";
        ofs << "V_fluid "        << V_fluid        << "\n";
        ofs << "G_grav "         << G_grav         << "\n";
        ofs << "Lambda "         << Lambda         << "\n";
        ofs << "c_light "        << c_light        << "\n";
        ofs << "hbar "           << hbar           << "\n";
        ofs << "delta_x "        << delta_x        << "\n";
        ofs << "delta_p "        << delta_p        << "\n";
        ofs << "integral_psi "   << integral_psi   << "\n";
        ofs << "t_Hubble "       << t_Hubble       << "\n";
        ofs << "H0 "             << H0             << "\n";
        ofs << "Omega_m "        << Omega_m        << "\n";
        ofs << "Omega_Lam "      << Omega_Lam      << "\n";
        ofs << "Mpc_to_m "       << Mpc_to_m       << "\n";
        ofs << "kappa_neb "      << kappa_neb      << "\n";
        ofs << "pre_sum_Ug "     << pre_sum_Ug     << "\n";
        ofs << "g_base_cache "   << g_base_cache   << "\n";
        for (auto& kv : dynamic_params)
            ofs << "dyn_" << kv.first << " " << kv.second << "\n";
        ofs.close();
    }

    template <typename T>
    double cross_validate(T& other_module, double t = 0.0) {
        double g_this  = computeG(t);
        double g_other = other_module.computeG(t);
        double denom   = std::abs(g_this) + std::abs(g_other) + 1.0e-300;
        return std::abs(g_this - g_other) / denom;
    }
};

// UQFF 2.0 Integration Note:
// #include "M16_UQFF_MODULE.cpp" in source2.cpp or MAIN_1_CoAnQi.cpp
// M16UQFFModule mod; mod.setEnableLogging(true); double g = mod.computeG(0.0);
// std::cout << mod.getEquationText(); mod.exportState();
// PAPER_284: Phi_dm = (1+SFR_rate*t)*(1-E_0*(1-exp(-t/tau))); first UQFF multiplicative SFR×erosion product
//            gap_mult_add = -(M_sf_frac×E_rad) = -1013.3 at 5 Myr (24.3% below additive)
// PAPER_285: t_half = tau*ln2 = 6.561e13 s = 2.079 Myr; DeltaGMax = 4.36e-13 m/s²
// PAPER_286: kappa_neb = 6.71e-4; H(z=0.0015) = 70.047 km/s/Mpc; FIRST UQFF nebular kappa_neb

#endif // M16_UQFF_MODULE_H