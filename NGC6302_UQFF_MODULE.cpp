// NGC6302_UQFF_MODULE.cpp
// Session 89 — UQFF 2.0 | 31st C++ Module | FIRST Bipolar Planetary Nebula (Bug Nebula, PN) UQFF
// Copyright - Daniel T. Murphy | Analyzed: Oct 09, 2025 | UQFF 2.0 Upgrade: Session 89, Mar 17, 2026
//
// NGC 6302 ("Bug Nebula / Butterfly Nebula") — Bipolar Planetary Nebula
//   Central star: T_eff ≈ 200,000 K, L_star ≈ 5000 L_sun (one of hottest known WDs, HST Szyszka 2009)
//   Total PN gas mass: M = 2.0 M_sun = 3.978e30 kg (ejected bipolar lobes)
//   Half-lobe radius:  r ≈ 1 ly = 9.46e15 m
//   Wind speed:        v_wind = 100 km/s = 1.0e5 m/s (central star fast wind)
//   Ejection age:      t_eject = 2000 yr = 6.312e10 s (bipolar lobe age)
//   Redshift:          z = 0.00095  |  Distance ≈ 1.2 kpc
//   Ionized gas:       ρ = 1e-20 kg/m³  |  Equatorial torus B = 1e-5 T
//
// WOLFRAM_TERM: NGC6302_BASE             | g_base = G*M/r² = 2.967e-12 m/s²
// WOLFRAM_TERM: NGC6302_WIND_SHOCK       | PAPER_311: η_wind=7.127e5, a_wind(t_ej)=2.114e-6 m/s²
// WOLFRAM_TERM: NGC6302_UV_RADIATION     | PAPER_312: η_rad=1.913e20, a_rad=5.672e8 m/s²
// WOLFRAM_TERM: NGC6302_TORUS_CONFINEMENT| PAPER_313: η_B_conf=3.979e5, v_Alfven=8.921e7 m/s
//
// MANDATORY ARCHITECTURE RULES: Parameterized calculator only — no hardcoded system instances.
// Include in base: #include "NGC6302_UQFF_MODULE.h"
// Usage: NGC6302UQFFModule mod; mod.computeG(t); mod.updateVariable("v_wind", 2e5);

#ifndef NGC6302_UQFF_MODULE_H
#define NGC6302_UQFF_MODULE_H

#include <string>
#include <map>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <complex>
#include <fstream>

// ============================================================
// WOLFRAM_TERM: NGC6302_BASE — Bipolar PN gravitational anchor
// ============================================================
class NGC6302UQFFModule {
public:
    // ---- Universal Physical Constants ----
    static constexpr double G_NEWTON    = 6.6743e-11;        // m³ kg⁻¹ s⁻²
    static constexpr double C_LIGHT     = 3.0e8;             // m/s
    static constexpr double HBAR        = 1.0546e-34;        // J·s
    static constexpr double LAMBDA_CC   = 1.1e-52;           // m⁻²
    static constexpr double Q_PROTON    = 1.602e-19;         // C
    static constexpr double M_H         = 1.673e-27;         // kg
    static constexpr double MU_0        = 1.2566e-6;         // H/m (permeability of free space)
    static constexpr double L_SUN       = 3.828e26;          // W
    static constexpr double PI          = 3.141592653589793;
    static constexpr double YR_TO_S     = 3.156e7;           // s/yr
    static constexpr double T_HUBBLE_S  = 13.8e9 * 3.156e7;  // s
    static constexpr double MPC_TO_M    = 3.086e22;          // m/Mpc
    static constexpr double M_SUN       = 1.989e30;          // kg

private:
    // ---- NGC 6302 Typed Physics Members ----
    double M;           // Total PN ejected mass [kg]  = 2.0 M_sun = 3.978e30
    double M_star;      // Central WD mass [kg]        = 0.64 M_sun = 1.273e30
    double r;           // PN half-lobe radius [m]     ≈ 1 ly = 9.46e15
    double t_eject;     // Nebular ejection age [s]    = 2000 yr = 6.312e10
    double v_wind;      // Stellar wind speed [m/s]    = 1.0e5 (100 km/s)
    double T_eff_star;  // Central star T_eff [K]      = 2.0e5
    double L_star;      // Central star luminosity [W] = 5000 L_sun = 1.914e30
    double rho_fluid;   // PN ionized gas ρ [kg/m³]   = 1.0e-20
    double V_fluid;     // Reference volume [m³]       = 1.0e3
    double B;           // Equatorial torus B [T]      = 1.0e-5
    double B_crit;      // Critical B field [T]        = 1.0e11
    double z;           // Redshift                    = 9.5e-4
    double H0;          // H₀ [km/s/Mpc]              = 67.15
    double Omega_m;     //                             = 0.30
    double Omega_L;     //                             = 0.70
    double f_TRZ;       // TRZ resonance factor        = 0.10
    double f_DM;        // Dark matter fraction        = 0.85
    double f_vis;       // Visible fraction            = 0.15
    double delta_rho;   // Density perturbation        = 0.1 * rho_fluid
    double Delta_x;     // Quantum uncertainty [m]     = 1.0e-10
    double A_osc;       // Oscillatory amplitude       = 1.0e-10
    double k_osc;       // Wavenumber [1/m]            = 1.0e20
    double omega_osc;   // Frequency [rad/s]           = 1.0e15

    // ---- PAPER_311 / PAPER_312 / PAPER_313 Caches ----
    // PAPER_311: Wind Shock Gravitational Dominance
    double g_base_cache;              // 2.967e-12 m/s²
    double a_wind_t0_cache;           // v_wind²/r = 1.057e-6 m/s²
    double a_wind_teject_cache;       // 2*v_wind²/r = 2.114e-6 m/s²
    double eta_wind_cache;            // 7.127e5
    double KE_over_gwell_cache;       // v_wind²/(G M / r) = 3.564e5
    // PAPER_312: UV Radiation Pressure
    double P_rad_cache;               // 5.672e-12 Pa
    double a_rad_cache;               // 5.672e8 m/s²
    double eta_rad_cache;             // 1.913e20
    // PAPER_313: Torus Magnetic Confinement
    double P_mag_cache;               // 3.979e-5 Pa
    double P_ram_cache;               // 1.000e-10 Pa
    double eta_B_conf_cache;          // 3.979e5
    double beta_plasma_cache;         // 2.513e-6
    double v_Alfven_cache;            // 8.921e7 m/s
    double v_Alfven_over_vwind_cache; // 892.1

    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

    // ---- Cache Computation ----
    void updateCache() {
        g_base_cache = G_NEWTON * M / (r * r);

        // PAPER_311: a_wind = v_wind² * (1 + t/t_eject) / r  [properly m/s²]
        a_wind_t0_cache     = (v_wind * v_wind) / r;
        a_wind_teject_cache = 2.0 * (v_wind * v_wind) / r;
        eta_wind_cache      = a_wind_teject_cache / g_base_cache;
        KE_over_gwell_cache = (v_wind * v_wind) / (G_NEWTON * M / r);

        // PAPER_312: radiation pressure from central UV WD
        double four_pi_r2_c = 4.0 * PI * r * r * C_LIGHT;
        P_rad_cache  = L_star / four_pi_r2_c;
        a_rad_cache  = P_rad_cache / rho_fluid;
        eta_rad_cache = a_rad_cache / g_base_cache;

        // PAPER_313: equatorial torus magnetic confinement
        P_mag_cache          = (B * B) / (2.0 * MU_0);
        P_ram_cache          = rho_fluid * v_wind * v_wind;
        eta_B_conf_cache     = P_mag_cache / P_ram_cache;
        beta_plasma_cache    = P_ram_cache / P_mag_cache;
        v_Alfven_cache       = B / std::sqrt(MU_0 * rho_fluid);
        v_Alfven_over_vwind_cache = v_Alfven_cache / v_wind;
    }

    double computeHz(double z_val) const {
        double Hz_kms = H0 * std::sqrt(Omega_m * std::pow(1.0 + z_val, 3.0) + Omega_L);
        return (Hz_kms * 1e3) / MPC_TO_M;
    }

    double computeUgSum() const {
        double Ug1 = G_NEWTON * M / (r * r);
        return Ug1 + Ug1; // Ug4 = Ug1 * f_sc (f_sc=1); Ug2=Ug3=0
    }

    double computeQuantumTerm() const {
        // hbar / (m_H * Delta_x²) [m/s²]
        return HBAR / (M_H * Delta_x * Delta_x);
    }

    double computeResonantTerm(double t) const {
        double cos_term = A_osc * std::cos(k_osc * r + omega_osc * t);
        std::complex<double> exp_arg(0.0, k_osc * r - omega_osc * t);
        double exp_term = A_osc * std::exp(exp_arg).real() * (2.0 * PI / 13.8);
        return cos_term + exp_term;
    }

    // PAPER_311: Wind Shock — properly normalized to m/s² via r
    double computeWindTerm(double t) const {
        return (v_wind * v_wind) * (1.0 + t / t_eject) / r;
    }

    // PAPER_312: UV radiation pressure acceleration from central WD
    double computeRadiationTerm() const {
        return L_star / (4.0 * PI * r * r * C_LIGHT * rho_fluid);
    }

    double computeDMTerm() const {
        double pert = delta_rho / rho_fluid;
        double curv = 3.0 * G_NEWTON * M / (r * r * r);
        return M * (pert + curv);
    }

    double getVariableValue(const std::string& key) const {
        auto it = dynamic_params.find(key);
        return (it != dynamic_params.end()) ? it->second : 0.0;
    }

public:
    explicit NGC6302UQFFModule() :
        M(2.0 * M_SUN),
        M_star(0.64 * M_SUN),
        r(9.46e15),
        t_eject(2000.0 * YR_TO_S),
        v_wind(1.0e5),
        T_eff_star(2.0e5),
        L_star(5000.0 * L_SUN),
        rho_fluid(1.0e-20),
        V_fluid(1.0e3),
        B(1.0e-5),
        B_crit(1.0e11),
        z(9.5e-4),
        H0(67.15),
        Omega_m(0.30),
        Omega_L(0.70),
        f_TRZ(0.10),
        f_DM(0.85),
        f_vis(0.15),
        delta_rho(1.0e-21),
        Delta_x(1.0e-10),
        A_osc(1.0e-10),
        k_osc(1.0e20),
        omega_osc(1.0e15),
        g_base_cache(0.0),
        a_wind_t0_cache(0.0),
        a_wind_teject_cache(0.0),
        eta_wind_cache(0.0),
        KE_over_gwell_cache(0.0),
        P_rad_cache(0.0),
        a_rad_cache(0.0),
        eta_rad_cache(0.0),
        P_mag_cache(0.0),
        P_ram_cache(0.0),
        eta_B_conf_cache(0.0),
        beta_plasma_cache(0.0),
        v_Alfven_cache(0.0),
        v_Alfven_over_vwind_cache(0.0),
        logging_enabled(false)
    {
        updateCache();
    }

    void setEnableLogging(bool enabled) { logging_enabled = enabled; }

    void setDynamicParameter(const std::string& key, double value) {
        dynamic_params[key] = value;
        if (logging_enabled)
            std::cout << "[NGC6302_DYN] set " << key << "=" << value << "\n";
    }

    double getDynamicParameter(const std::string& key) const {
        auto it = dynamic_params.find(key);
        return (it != dynamic_params.end()) ? it->second : 0.0;
    }

    void updateVariable(const std::string& name, double value) {
        if      (name == "M")          { M          = value; }
        else if (name == "M_star")     { M_star     = value; }
        else if (name == "r")          { r          = value; }
        else if (name == "t_eject")    { t_eject    = value; }
        else if (name == "v_wind")     { v_wind     = value; }
        else if (name == "T_eff_star") { T_eff_star = value; }
        else if (name == "L_star")     { L_star     = value; }
        else if (name == "rho_fluid")  { rho_fluid  = value; delta_rho = 0.1 * rho_fluid; }
        else if (name == "V_fluid")    { V_fluid    = value; }
        else if (name == "B")          { B          = value; }
        else if (name == "B_crit")     { B_crit     = value; }
        else if (name == "z")          { z          = value; }
        else if (name == "H0")         { H0         = value; }
        else if (name == "f_TRZ")      { f_TRZ      = value; }
        else if (name == "f_DM")       { f_DM       = value; f_vis = 1.0 - value; }
        else if (name == "Delta_x")    { Delta_x    = value; }
        else if (name == "A_osc")      { A_osc      = value; }
        else if (name == "omega_osc")  { omega_osc  = value; }
        else { setDynamicParameter(name, value); return; }
        updateCache();
    }

    void addToVariable(const std::string& name, double delta) {
        if      (name == "M")         updateVariable(name, M         + delta);
        else if (name == "v_wind")    updateVariable(name, v_wind    + delta);
        else if (name == "L_star")    updateVariable(name, L_star    + delta);
        else if (name == "rho_fluid") updateVariable(name, rho_fluid + delta);
        else if (name == "B")         updateVariable(name, B         + delta);
        else {
            dynamic_params[name] += delta;
            if (logging_enabled)
                std::cout << "[NGC6302_ADD] " << name << " += " << delta << "\n";
        }
    }

    void subtractFromVariable(const std::string& name, double delta) {
        addToVariable(name, -delta);
    }

    // ---- Full UQFF 2.0 Pipeline ----
    // g_NGC6302(r,t) = g_base*(1+Hz*t)*(1-B/B_crit)*(1+f_TRZ)
    //                + (Ug1+Ug2+Ug3+Ug4) + Lambda*Omega_L*c²/3
    //                + hbar/(m_H*Delta_x²)  [quantum HUP]
    //                + q*v_wind*B/m_H        [EM Lorentz]
    //                + rho*V*g_base/r³       [fluid coupling, dim-consistent]
    //                + A*cos(k*r+omega*t)    [resonant oscillatory]
    //                + M*(delta_rho/rho+3GM/r³) [DM perturbation]
    //                + v_wind²*(1+t/t_eject)/r  [PAPER_311: wind shock, fixed units]
    //                + L_star/(4pi*r²*c*rho)    [PAPER_312: UV radiation]
    double computeG(double t, double z_override = -1.0) {
        double z_use    = (z_override >= 0.0) ? z_override : z;
        double Hz       = computeHz(z_use);
        double expansion = 1.0 + Hz * t;
        double sc_factor = 1.0 - (B / B_crit);
        double tr_factor = 1.0 + f_TRZ;

        double g_base    = G_NEWTON * M / (r * r) * expansion * sc_factor * tr_factor;
        double ug_sum    = computeUgSum();
        double lambda_term   = LAMBDA_CC * Omega_L * C_LIGHT * C_LIGHT / 3.0;
        double quantum_term  = computeQuantumTerm();
        double em_term       = Q_PROTON * v_wind * B / M_H;
        double fluid_term    = rho_fluid * V_fluid * g_base / (rho_fluid * r * r * r);
        double resonant_term = computeResonantTerm(t);
        double dm_term       = computeDMTerm();
        double wind_term     = computeWindTerm(t);    // PAPER_311
        double rad_term      = computeRadiationTerm(); // PAPER_312

        double dyn_sum = 0.0;
        for (const auto& kv : dynamic_params) dyn_sum += kv.second;

        double total = g_base + ug_sum + lambda_term + quantum_term + em_term
                     + fluid_term + resonant_term + dm_term + wind_term + rad_term + dyn_sum;

        if (logging_enabled) {
            std::cout << std::scientific << std::setprecision(4)
                      << "[NGC6302_G] t=" << t
                      << " g_base=" << g_base << " wind=" << wind_term
                      << " rad=" << rad_term << " total=" << total << "\n";
        }
        return total;
    }

    std::string getEquationText() const {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(4);
        oss << "=== NGC 6302 (Bug Nebula) UQFF 2.0 Full Gravity Pipeline ===\n"
            << "g_NGC6302(r,t) = g_base*(1+H(z)*t)*(1-B/B_crit)*(1+f_TRZ)\n"
            << "              + (Ug1+Ug2+Ug3+Ug4)\n"
            << "              + Lambda*Omega_L*c^2/3\n"
            << "              + hbar/(m_H*Delta_x^2)            [quantum HUP]\n"
            << "              + q*v_wind*B/m_H                  [EM Lorentz]\n"
            << "              + rho*V*g_base/r^3                [fluid coupling]\n"
            << "              + A*cos(k*r+omega*t)              [resonant]\n"
            << "              + M*(delta_rho/rho+3GM/r^3)       [DM perturbation]\n"
            << "              + v_wind^2*(1+t/t_eject)/r        [PAPER_311: wind shock]\n"
            << "              + L_star/(4pi*r^2*c*rho)          [PAPER_312: UV radiation]\n"
            << "\n--- PAPER_311: Wind Shock Gravitational Dominance ---\n"
            << "  g_base               = " << g_base_cache << " m/s^2\n"
            << "  a_wind(t=0)          = " << a_wind_t0_cache << " m/s^2\n"
            << "  a_wind(t=t_eject)    = " << a_wind_teject_cache << " m/s^2\n"
            << "  eta_wind             = " << eta_wind_cache << "  [wind/gravity]\n"
            << "  KE/grav_well         = " << KE_over_gwell_cache << "  [v^2/(GM/r)]\n"
            << "\n--- PAPER_312: Central Star UV Radiation Pressure ---\n"
            << "  L_star (5000 L_sun)  = " << L_star << " W\n"
            << "  P_rad                = " << P_rad_cache << " Pa\n"
            << "  a_rad                = " << a_rad_cache << " m/s^2\n"
            << "  eta_rad              = " << eta_rad_cache << "  [rad/gravity]\n"
            << "\n--- PAPER_313: Equatorial Torus Magnetic Confinement ---\n"
            << "  P_mag = B^2/(2*mu0)  = " << P_mag_cache << " Pa\n"
            << "  P_ram = rho*v_wind^2 = " << P_ram_cache << " Pa\n"
            << "  eta_B_conf           = " << eta_B_conf_cache << "  [P_mag/P_ram]\n"
            << "  beta_plasma          = " << beta_plasma_cache << "  [P_ram/P_mag]\n"
            << "  v_Alfven             = " << v_Alfven_cache << " m/s\n"
            << "  v_Alfven/v_wind      = " << v_Alfven_over_vwind_cache << "\n";
        return oss.str();
    }

    void printParameters() const {
        std::cout << getEquationText();
        std::cout << "\n--- Active System Parameters ---\n"
                  << std::scientific << std::setprecision(4)
                  << "  M          = " << M          << " kg (2.0 M_sun)\n"
                  << "  M_star     = " << M_star     << " kg (0.64 M_sun, T_eff=" << T_eff_star << "K)\n"
                  << "  r          = " << r          << " m (~1 ly)\n"
                  << "  t_eject    = " << t_eject    << " s (2000 yr)\n"
                  << "  v_wind     = " << v_wind     << " m/s\n"
                  << "  L_star     = " << L_star     << " W (5000 L_sun)\n"
                  << "  rho_fluid  = " << rho_fluid  << " kg/m^3\n"
                  << "  B          = " << B          << " T  |  B_crit=" << B_crit << " T\n"
                  << "  z          = " << z          << "\n"
                  << "  H0         = " << H0         << " km/s/Mpc\n"
                  << "  f_TRZ      = " << f_TRZ      << "\n"
                  << "  f_DM/f_vis = " << f_DM << "/" << f_vis << "\n";
    }

    // Legacy alias
    void printVariables() const { printParameters(); }

    void exportState(const std::string& filename) const {
        std::ofstream ofs(filename);
        if (!ofs.is_open()) { std::cerr << "[NGC6302_EXPORT] Cannot open " << filename << "\n"; return; }
        ofs << std::scientific << std::setprecision(6);
        ofs << "# NGC6302 UQFF 2.0 State Export — Session 89\n"
            << "M=" << M << "\nM_star=" << M_star << "\nr=" << r
            << "\nt_eject=" << t_eject << "\nv_wind=" << v_wind
            << "\nT_eff_star=" << T_eff_star << "\nL_star=" << L_star
            << "\nrho_fluid=" << rho_fluid << "\nB=" << B << "\nB_crit=" << B_crit
            << "\nz=" << z << "\nH0=" << H0 << "\nOmega_m=" << Omega_m
            << "\nOmega_L=" << Omega_L << "\nf_TRZ=" << f_TRZ
            << "\nf_DM=" << f_DM << "\nf_vis=" << f_vis
            << "\nDelta_x=" << Delta_x << "\nA_osc=" << A_osc
            << "\nomega_osc=" << omega_osc << "\n"
            << "# Cached physics (PAPER_311/312/313)\n"
            << "g_base=" << g_base_cache
            << "\na_wind_t0=" << a_wind_t0_cache
            << "\na_wind_teject=" << a_wind_teject_cache
            << "\neta_wind=" << eta_wind_cache
            << "\nKE_over_gwell=" << KE_over_gwell_cache
            << "\nP_rad=" << P_rad_cache
            << "\na_rad=" << a_rad_cache
            << "\neta_rad=" << eta_rad_cache
            << "\nP_mag=" << P_mag_cache
            << "\nP_ram=" << P_ram_cache
            << "\neta_B_conf=" << eta_B_conf_cache
            << "\nbeta_plasma=" << beta_plasma_cache
            << "\nv_Alfven=" << v_Alfven_cache
            << "\nv_Alfven_over_vwind=" << v_Alfven_over_vwind_cache << "\n";
        for (const auto& kv : dynamic_params)
            ofs << "dyn:" << kv.first << "=" << kv.second << "\n";
        std::cout << "[NGC6302_EXPORT] State saved to " << filename << "\n";
    }

    template<typename T>
    T cross_validate(const std::string& param, T alt_value) const {
        if (logging_enabled)
            std::cout << "[NGC6302_XVAL] param=" << param << " alt=" << alt_value << "\n";
        return alt_value;
    }
};

#endif // NGC6302_UQFF_MODULE_H
