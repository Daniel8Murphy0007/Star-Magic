/**
 * ================================================================================================
 * Header: MagnetarSGR1745_2900.h
 * 
 * Description: C++ Module for SGR 1745-2900 Magnetar Class
 *              This is the second module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on magnetar evolution and gravity
 *              equations derived from Chandra X-ray Observatory datasets, high-energy lab simulations,
 *              and UQFF refinements (dated May 11, 2025, updated for full term inclusion on October 08, 2025).
 * 
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for SGR 1745-2900 magnetar
 *          evolution, including black hole proximity (Sgr A*), magnetic energy, and outburst decay.
 *          Includes ALL terms: base gravity, cosmic expansion (H(z)), magnetic decay, BH influence,
 *          UQFF Ug components with f_sc, Lambda, quantum uncertainty, EM, fluid, oscillatory waves,
 *          DM/density perturbations, magnetic energy (effective g), and decay power (cumulative energy effective g).
 *          Supports dynamic variable updates for all parameters.
 * 
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: MagnetarSGR1745_2900 mag;
 *              Compute: double g = mag.compute_g_Magnetar(t);
 * 
 * Key Features:
 *   - UQFF canonical: B=2e10 T, B_crit=1e11 T, r=20 km, Hz=2.269e-18 s⁻¹ [CP3 SGR1745BHProximityMagEnergyCalculator].
 *   - 15-term MUGE: base+Hz+f_sc, Ug(f_TRZ), Lambda, EM(dynamic v_surf), GW, BH monopole,
 *     BH 2nd-order tidal gradient, magnetic stored energy (dynamic B), burst decay (cumulative),
 *     oscillatory burst D(t), quantum, fluid, oscillatory (Mode 1), DM+tidal.
 *   - Dynamic B(t) = B0·exp(−t/tau_B): f_sc, M_mag, and EM all evolve with spin-down.
 *   - Burst modulation D(t) = D₀·cos(ω_D·t)·e^(−t/τ_D) from MagnetarSGR1745DynamicModulationCalculator (CP3).
 *   - BH tidal gradient: a_tidal = 2·G·M_BH/r_BH³·r (2nd-order at 0.92 pc from Sgr A*).
 *   - Dynamic v_surf: (2π/P_init)·r from ATNF pulse period P_init=3.76 s.
 *   - UQFF 2.0 self-expanding framework: dynamic parameters, logging, state export.
 *   - cross_validate() utility for dual-method UQFF/MUGE cross-validation with SGR0501.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 * 
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef MAGNETAR_SGR1745_2900_H
#define MAGNETAR_SGR1745_2900_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <limits>
#include <map>
#include <fstream>

// Portable PI constant (M_PI is non-standard; not guaranteed by MSVC without _USE_MATH_DEFINES)
static constexpr double UQFF_PI = 3.14159265358979323846;

class MagnetarSGR1745_2900 {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M;               // Magnetar mass
    double r;               // Radius
    double Hz;              // Hubble parameter at z (s^-1)
    double B0;              // Initial magnetic field
    double tau_B;           // B decay timescale (s) - not used in this eq, but for consistency
    double B_crit;          // Critical B field
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double v_surf;          // Surface velocity
    double f_sc;            // Superconductive factor (computed as 1 - B/B_crit)
    double rho_vac_UA;      // UA vacuum density - not used here
    double rho_vac_SCm;     // SCm vacuum density - not used here
    double P_init;          // Initial rotation period (s)
    double tau_Omega;       // Omega decay timescale (s)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass for EM acceleration
    double f_TRZ;           // Time-reversal factor [CP3 canonical: 0.1]
    double M_BH;            // Black hole mass
    double r_BH;            // Distance to black hole
    double mu0;             // Vacuum permeability
    double L0_W;            // Initial luminosity (W)
    double tau_decay;       // Decay timescale (s)
    // Upgrade 2: D(t) burst modulation (MagnetarSGR1745DynamicModulationCalculator, CP3)
    double D0_burst;        // Burst modulation amplitude (m/s^2)
    double omega_D;         // Burst modulation angular frequency (rad/s)
    double tau_D;           // Burst modulation decay timescale (s)
    double hbar;            // Reduced Planck's constant
    double t_Hubble;        // Hubble time (s)
    double t_Hubble_gyr;    // Hubble time in Gyr
    double delta_x;         // Position uncertainty (m)
    double delta_p;         // Momentum uncertainty (kg m/s)
    double integral_psi;    // Wavefunction integral approximation
    double rho_fluid;       // Fluid density (kg/m^3)
    double A_osc;           // Oscillatory amplitude (m/s^2)
    double k_osc;           // Wave number (1/m)
    double omega_osc;       // Angular frequency (rad/s)
    double x_pos;           // Position for oscillation (m)
    double M_DM_factor;     // Dark matter mass fraction
    double delta_rho_over_rho; // Density perturbation fraction

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 = G*M/r^2
    double B;               // Current B (static in this model)

    // UQFF 2.0 Self-Expanding Framework
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

public:
    // Constructor with default UQFF values
    MagnetarSGR1745_2900() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~MagnetarSGR1745_2900() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        M = 1.4 * 1.989e30;
        r = 20e3;                          // 20 km neutron star radius [CP3 canonical]
        Hz = 2.269e-18;                    // Galactic center H(z≈0.001) [CP3 canonical] (s^-1)
        B0 = 2e10;                         // SGR 1745-2900 B field [UQFF calibration] (T)
        B = B0;                            // Static B model
        tau_B = 4000 * 3.15576e7;         // Reserved for future dynamic B model
        B_crit = 1e11;                     // UQFF framework calibration constant (T)
        Lambda = 1.114e-52;                // Cosmological constant [CP3 canonical] (m^-2)
        c_light = 2.998e8;                 // Speed of light [CP3 canonical] (m/s)
        q_charge = 1.602e-19;
        v_surf = 0.0;                      // Computed in updateCache() = (2pi/P_init)*r
        f_sc = 0.0;                        // Computed in updateCache() = 1 - B/B_crit
        f_TRZ = 0.1;                       // Time-reversal factor [CP3 canonical]
        rho_vac_UA = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        P_init = 3.76;                     // ATNF real pulse period (s) [CP3 canonical]
        tau_Omega = 10000 * 3.15576e7;
        scale_EM = 1e-12;
        proton_mass = 1.673e-27;
        M_BH = 4e6 * 1.989e30;            // Sgr A* mass [CP3 canonical]
        r_BH = 2.83e16;                    // 0.92 pc distance to Sgr A* [CP3 canonical] (m)
        mu0 = 1.2566e-6;                   // Permeability of free space [CP3 canonical] (H/m)
        L0_W = 5e28;                       // Initial burst luminosity [CP3 canonical] (W)
        tau_decay = 3.5 * 3.15576e7;       // 3.5 yr burst decay [CP3 canonical] (s)
        // Upgrade 2: D(t) burst modulation parameters [MagnetarSGR1745DynamicModulationCalculator]
        D0_burst  = 1e-3;                   // Burst modulation amplitude (m/s^2)
        omega_D   = 2.0 * UQFF_PI / 11.0;  // ~11 s burst repeat period (rad/s)
        tau_D     = 3.5 * 3.15576e7;        // Burst modulation decay timescale = tau_decay (s)
        hbar = 1.0546e-34;
        t_Hubble = 1.0 / Hz;              // Hubble time = 1/Hz [CP3 canonical] (s)
        t_Hubble_gyr = 13.8;
        delta_x = 1e-15;                  // Nuclear scale for neutron star [CP3 canonical] (m)
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
        rho_fluid = 1e-9;                 // Magnetospheric ambient density [CP3 canonical] (kg/m^3)
        A_osc = 1e10;
        k_osc = 0.0;                      // Computed in updateCache() = 2pi/r
        omega_osc = 0.0;                  // Computed in updateCache() = 2pi*c/r
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;

        logging_enabled = false;
        updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        if (r <= 0.0) {
            std::cerr << "Warning: r <= 0; cache not updated." << std::endl;
            return;
        }
        if (B_crit == 0.0) {
            std::cerr << "Warning: B_crit == 0; cache not updated." << std::endl;
            return;
        }
        ug1_base  = (G * M) / (r * r);
        f_sc      = 1.0 - (B / B_crit);                      // Superconductive suppression
        v_surf    = (2.0 * UQFF_PI / P_init) * r;           // Dynamic from ATNF P_init [CP3]
        k_osc     = 2.0 * UQFF_PI / r;                      // k = 2π/r [CP3 canonical]
        omega_osc = 2.0 * UQFF_PI * c_light / r;            // ω = 2π·c/r [CP3 canonical]
        x_pos     = r;
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M") { M = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "Hz") { Hz = newValue; }
        else if (varName == "B0") { B0 = newValue; B = newValue; }
        else if (varName == "tau_B") { tau_B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "v_surf") { v_surf = newValue; }
        else if (varName == "f_sc") { f_sc = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "P_init") { P_init = newValue; }
        else if (varName == "tau_Omega") { tau_Omega = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
        else if (varName == "proton_mass") { proton_mass = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "M_BH") { M_BH = newValue; }
        else if (varName == "r_BH") { r_BH = newValue; }
        else if (varName == "mu0") { mu0 = newValue; }
        else if (varName == "L0_W") { L0_W = newValue; }
        else if (varName == "tau_decay") { tau_decay = newValue; }
        else if (varName == "D0_burst") { D0_burst = newValue; }
        else if (varName == "omega_D")  { omega_D  = newValue; }
        else if (varName == "tau_D")    { tau_D    = newValue; }
        // Full terms
        else if (varName == "hbar") { hbar = newValue; }
        else if (varName == "t_Hubble") { t_Hubble = newValue; }
        else if (varName == "t_Hubble_gyr") { t_Hubble_gyr = newValue; }
        else if (varName == "delta_x") { delta_x = newValue; }
        else if (varName == "delta_p") { delta_p = newValue; }
        else if (varName == "integral_psi") { integral_psi = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
        else if (varName == "A_osc") { A_osc = newValue; }
        else if (varName == "k_osc") { k_osc = newValue; }
        else if (varName == "omega_osc") { omega_osc = newValue; }
        else if (varName == "x_pos") { x_pos = newValue; }
        else if (varName == "M_DM_factor") { M_DM_factor = newValue; }
        else if (varName == "delta_rho_over_rho") { delta_rho_over_rho = newValue; }
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return false;
        }
        updateCache();
        return true;
    }

    // Addition method for variables
    bool addToVariable(const std::string& varName, double delta) {
        return setVariable(varName, getVariable(varName) + delta);
    }

    // Subtraction method for variables
    bool subtractFromVariable(const std::string& varName, double delta) {
        return addToVariable(varName, -delta);
    }

    // Getter for any variable (helper for add/subtract)
    double getVariable(const std::string& varName) const {
        if (varName == "G") return G;
        else if (varName == "M") return M;
        else if (varName == "r") return r;
        else if (varName == "Hz") return Hz;
        else if (varName == "B0") return B0;
        else if (varName == "tau_B") return tau_B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "v_surf") return v_surf;
        else if (varName == "f_sc") return f_sc;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "P_init") return P_init;
        else if (varName == "tau_Omega") return tau_Omega;
        else if (varName == "scale_EM") return scale_EM;
        else if (varName == "proton_mass") return proton_mass;
        else if (varName == "f_TRZ") return f_TRZ;
        else if (varName == "M_BH") return M_BH;
        else if (varName == "r_BH") return r_BH;
        else if (varName == "mu0") return mu0;
        else if (varName == "L0_W") return L0_W;
        else if (varName == "tau_decay") return tau_decay;
        else if (varName == "D0_burst") return D0_burst;
        else if (varName == "omega_D")  return omega_D;
        else if (varName == "tau_D")    return tau_D;
        // Full terms
        else if (varName == "hbar") return hbar;
        else if (varName == "t_Hubble") return t_Hubble;
        else if (varName == "t_Hubble_gyr") return t_Hubble_gyr;
        else if (varName == "delta_x") return delta_x;
        else if (varName == "delta_p") return delta_p;
        else if (varName == "integral_psi") return integral_psi;
        else if (varName == "rho_fluid") return rho_fluid;
        else if (varName == "A_osc") return A_osc;
        else if (varName == "k_osc") return k_osc;
        else if (varName == "omega_osc") return omega_osc;
        else if (varName == "x_pos") return x_pos;
        else if (varName == "M_DM_factor") return M_DM_factor;
        else if (varName == "delta_rho_over_rho") return delta_rho_over_rho;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return std::numeric_limits<double>::quiet_NaN();  // NaN prevents silent corruption
        }
    }

    // B(t) - UPGRADE 1: dynamic exponential decay  [B(t) = B0·exp(−t/tau_B)]
    // Activates when tau_B > 0; previously static. Drives dynamic f_sc, M_mag, EM.
    double B_t(double t) const {
        return B0 * std::exp(-t / tau_B);
    }

    // Omega(t) computation
    double Omega_t(double t) const {
        return (2.0 * UQFF_PI / P_init) * std::exp(-t / tau_Omega);
    }

    // dOmega/dt computation
    double dOmega_dt(double t) const {
        double omega0 = 2.0 * UQFF_PI / P_init;
        return omega0 * (-1.0 / tau_Omega) * std::exp(-t / tau_Omega);
    }

    // Ug terms computation  [CP3: (Ug1 + Ug4) × (1 + f_TRZ), superconductive suppression]
    double compute_Ug() const {
        double Ug1 = ug1_base;
        double Ug4 = Ug1 * f_sc;          // f_sc = 1 − B/B_crit
        return (Ug1 + Ug4) * (1.0 + f_TRZ);
    }

    // Volume computation
    double compute_V() const {
        return (4.0 / 3.0) * UQFF_PI * r * r * r;
    }

    // Magnetic energy M_mag(t) (J)  — now uses dynamic B(t) [Upgrade 1]
    double compute_M_mag(double t = 0.0) const {
        double Bt = B_t(t);
        double V  = compute_V();
        return (Bt * Bt / (2.0 * mu0)) * V;
    }

    // Cumulative decay energy up to t (J)  [CP3: L₀·τ_d·(1−e^(−t/τ_d))]
    double compute_cumulative_D(double t) const {
        return L0_W * tau_decay * (1.0 - std::exp(-t / tau_decay));
    }

    // Main MUGE computation (includes ALL terms)
    double compute_g_Magnetar(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Bt = B_t(t);         // Upgrade 1: dynamic B(t) = B0*exp(-t/tau_B)
        double dOdt = dOmega_dt(t);
        double current_f_sc = 1.0 - (Bt / B_crit);  // Evolves with B(t)

        // Term 1: Base + H(z) + B corrections
        double corr_H = 1 + Hz * t;
        double corr_B = current_f_sc;
        double term1 = ug1_base * corr_H * corr_B;

        // BH term (monopole)
        double term_BH = (G * M_BH) / (r_BH * r_BH);

        // UPGRADE 3: BH 2nd-order tidal gradient [a_tidal = 2·G·M_BH/r_BH³·r]
        // Differential tidal force across NS body at 0.92 pc from Sgr A*.
        // Not negligible: at r_BH=2.83e16 m, a_tidal ~ G*M_BH*r/r_BH^3 is measurable.
        double term_tidal = 2.0 * G * M_BH * r / (r_BH * r_BH * r_BH);

        // Term 2: UQFF Ug 
        double term2 = compute_Ug();

        // Term 3: Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // Term 4: Scaled EM  [CP3: includes UA/SCm density ratio factor; v_surf=(2π/P_init)·r]
        double cross_vB = v_surf * Bt;
        double em_base  = (q_charge * cross_vB) / proton_mass;
        double corr_UA  = 1.0 + (rho_vac_UA / rho_vac_SCm);
        double term4    = em_base * corr_UA * scale_EM;

        // Term 5: GW back-reaction  [CP3: (G·M²)/(c⁴·r)·(dΩ/dt)²]
        double gw_prefactor = (G * M * M) / (std::pow(c_light, 4) * r);
        double term5 = gw_prefactor * (dOdt * dOdt);

        // Quantum uncertainty term  [CP3: delta_x=1e-15 m, t_Hubble=1/Hz]
        double sqrt_unc = std::sqrt(delta_x * delta_p);
        double term_q   = (hbar / sqrt_unc) * integral_psi * (2.0 * UQFF_PI / t_Hubble);

        // Fluid term (effective acceleration)
        double V = compute_V();
        double term_fluid = (rho_fluid * V * ug1_base) / M;

        // Oscillatory term — Mode 1 standing wave  [CP3: k=2π/r, ω=2π·c/r; Mode 2 removed]
        double term_osc = 2.0 * A_osc * std::cos(k_osc * x_pos) * std::cos(omega_osc * t);

        // DM and density perturbation term (converted to acceleration)
        double M_dm = M * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M / (r * r * r);
        double term_dm_force_like = (M + M_dm) * (pert1 + pert2);
        double term_DM = term_dm_force_like / M;

        // Magnetic energy term (effective g)  — dynamic M_mag(t) [Upgrade 1]
        double M_mag = compute_M_mag(t);
        double term_mag = M_mag / (M * r);

        // Decay term (cumulative energy effective g)
        double cum_D = compute_cumulative_D(t);
        double term_decay = cum_D / (M * r);

        // UPGRADE 2: D(t) burst modulation  [MagnetarSGR1745DynamicModulationCalculator]
        // a_D(t) = D₀ · cos(ω_D · t) · e^(−t/τ_D)
        double term_burst = D0_burst * std::cos(omega_D * t) * std::exp(-t / tau_D);

        // Total g_Magnetar (all 15 terms)
        double result = term1 + term_BH + term_tidal + term2 + term3 + term4 + term5
                      + term_q + term_fluid + term_osc + term_DM + term_mag + term_decay + term_burst;
        if (logging_enabled)
            std::cout << "[LOG] compute_g_Magnetar(t=" << t << ") = " << result
                      << "  [t1=" << term1 << " tBH=" << term_BH
                      << " ttidal=" << term_tidal << " t2=" << term2
                      << " t3=" << term3 << " t4=" << term4 << " t5=" << term5
                      << " tq=" << term_q << " tf=" << term_fluid
                      << " tosc=" << term_osc << " tdm=" << term_DM
                      << " tmag=" << term_mag << " tdecay=" << term_decay
                      << " tburst=" << term_burst << "]" << std::endl;
        return result;
    }

    // Debug/Output method (comprehensive, for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::scientific << std::setprecision(6);
        os << "=== SGR 1745-2900 MUGE Parameters ===" << std::endl;
        os << "G:              " << G              << " m^3 kg^-1 s^-2" << std::endl;
        os << "M:              " << M              << " kg (1.4 M_sun)" << std::endl;
        os << "r:              " << r              << " m" << std::endl;
        os << "Hz:             " << Hz             << " s^-1 (Galactic center H)" << std::endl;
        os << "B0/B:           " << B              << " T" << std::endl;
        os << "B_crit:         " << B_crit         << " T" << std::endl;
        os << "Lambda:         " << Lambda         << " m^-2" << std::endl;
        os << "c_light:        " << c_light        << " m/s" << std::endl;
        os << "q_charge:       " << q_charge       << " C" << std::endl;
        os << "v_surf:         " << v_surf         << " m/s (= 2pi/P_init * r)" << std::endl;
        os << "f_TRZ:          " << f_TRZ          << std::endl;
        os << "f_sc:           " << f_sc           << " (= 1 - B/B_crit)" << std::endl;
        os << "P_init:         " << P_init         << " s (ATNF pulse period)" << std::endl;
        os << "tau_Omega:      " << tau_Omega      << " s" << std::endl;
        os << "scale_EM:       " << scale_EM       << std::endl;
        os << "rho_vac_UA:     " << rho_vac_UA     << " kg/m^3" << std::endl;
        os << "rho_vac_SCm:    " << rho_vac_SCm    << " kg/m^3" << std::endl;
        os << "proton_mass:    " << proton_mass     << " kg" << std::endl;
        os << "M_BH:           " << M_BH           << " kg (Sgr A*)" << std::endl;
        os << "r_BH:           " << r_BH           << " m (0.92 pc)" << std::endl;
        os << "mu0:            " << mu0            << " H/m" << std::endl;
        os << "L0_W:           " << L0_W           << " W" << std::endl;
        os << "tau_decay:      " << tau_decay      << " s" << std::endl;
        os << "D0_burst:       " << D0_burst        << " m/s^2 (burst modulation amplitude)" << std::endl;
        os << "omega_D:        " << omega_D         << " rad/s (burst repeat ~11 s)" << std::endl;
        os << "tau_D:          " << tau_D           << " s (burst decay timescale)" << std::endl;
        os << "--- Full-Term Parameters ---" << std::endl;
        os << "hbar:           " << hbar           << " J s" << std::endl;
        os << "t_Hubble:       " << t_Hubble       << " s (= 1/Hz)" << std::endl;
        os << "delta_x:        " << delta_x        << " m (nuclear scale)" << std::endl;
        os << "delta_p:        " << delta_p        << " kg m/s" << std::endl;
        os << "rho_fluid:      " << rho_fluid      << " kg/m^3 (magnetospheric)" << std::endl;
        os << "A_osc:          " << A_osc          << " m/s^2" << std::endl;
        os << "k_osc:          " << k_osc          << " rad/m (= 2pi/r)" << std::endl;
        os << "omega_osc:      " << omega_osc      << " rad/s (= 2pi*c/r)" << std::endl;
        os << "M_DM_factor:    " << M_DM_factor    << std::endl;
        os << "delta_rho/rho:  " << delta_rho_over_rho << std::endl;
        os << "--- Cached ---" << std::endl;
        os << "ug1_base:       " << ug1_base       << " m/s^2" << std::endl;
        os << "M_mag (J):      " << compute_M_mag(0.0) << " J (at t=0)" << std::endl;
    }

    // Example computation at t=1 year (for testing)
    double exampleAtOneYear() const {
        double t_example = 1.0 * 3.15576e7;
        return compute_g_Magnetar(t_example);
    }

    // Example computation at t=5000 years (CP3 canonical reference epoch)
    double exampleAt5000Years() const {
        double t_example = 5000.0 * 3.15576e7;
        return compute_g_Magnetar(t_example);
    }

    // ── UQFF 2.0 Self-Expanding Framework ────────────────────────────────────

    void setEnableLogging(bool enabled) { logging_enabled = enabled; }
    bool getLoggingEnabled() const      { return logging_enabled; }

    void setDynamicParameter(const std::string& key, double value) {
        dynamic_params[key] = value;
        if (logging_enabled)
            std::cout << "[LOG] Dynamic param set: " << key << " = " << value << std::endl;
    }

    double getDynamicParameter(const std::string& key) const {
        auto it = dynamic_params.find(key);
        if (it != dynamic_params.end()) return it->second;
        if (logging_enabled)
            std::cout << "[LOG] Dynamic param not found: " << key << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }

    void exportState(const std::string& filename) const {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "Error: Cannot open '" << filename << "' for state export." << std::endl;
            return;
        }
        ofs << "# MagnetarSGR1745_2900 state export\n";
        ofs << std::scientific << std::setprecision(10);
        const char* names[] = {
            "G","M","r","Hz","B0","tau_B","B_crit","Lambda","c_light",
            "q_charge","f_TRZ","rho_vac_UA","rho_vac_SCm",
            "P_init","tau_Omega","scale_EM","proton_mass",
            "M_BH","r_BH","mu0","L0_W","tau_decay",
            "D0_burst","omega_D","tau_D",
            "hbar","t_Hubble","delta_x","delta_p","integral_psi",
            "rho_fluid","A_osc","k_osc","omega_osc","x_pos","t_Hubble_gyr",
            "M_DM_factor","delta_rho_over_rho"
        };
        for (const char* n : names)
            ofs << n << " = " << getVariable(n) << "\n";
        for (const auto& kv : dynamic_params)
            ofs << "dynamic." << kv.first << " = " << kv.second << "\n";
        if (logging_enabled)
            std::cout << "[LOG] State exported to: " << filename << std::endl;
    }
    // UPGRADE 4: Cross-validate with MagnetarSGR0501_4516 (dual-method UQFF/MUGE pipeline)
    // Returns fractional difference |g1745 - g0501| / g0501 at the canonical t=5000 yr epoch.
    // Values < 0.01 indicate UQFF self-consistency across independent magnetar systems.
    // Include MAGNETAR_SGR0501_4516.cpp in the same translation unit before using.
    template<typename MagSGR0501>
    double cross_validate(const MagSGR0501& sgr0501, double t_years = 5000.0) const {
        double t_s    = t_years * 3.15576e7;
        double g1745  = compute_g_Magnetar(t_s);
        double g0501  = sgr0501.compute_g_Magnetar(t_s);
        double frac   = (g0501 != 0.0) ? std::abs(g1745 - g0501) / std::abs(g0501) : 
                        std::numeric_limits<double>::quiet_NaN();
        if (logging_enabled)
            std::cout << "[CROSS-VALIDATE] t=" << t_years << " yr  "
                      << "g_SGR1745=" << g1745 << "  g_SGR0501=" << g0501
                      << "  frac_diff=" << frac << std::endl;
        return frac;
    }
};

#endif // MAGNETAR_SGR1745_2900_H