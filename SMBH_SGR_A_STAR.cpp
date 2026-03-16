/**
 * ================================================================================================
 * Header: SMBHSgrAStar.h
 * 
 * Description: C++ Module for Sagittarius A* (Sgr A*) Supermassive Black Hole Class
 *              This is the third module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on SMBH evolution and gravity
 *              equations derived from Hubble datasets, high-energy lab simulations, and UQFF
 *              refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 * 
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for Sgr A* evolution.
 *          Includes ALL terms: base gravity with mass growth M(t), cosmic expansion (H_0), magnetic decay,
 *          UQFF Ug components with f_TRZ, Lambda, quantum uncertainty, EM (with B(t)), fluid dynamics,
 *          oscillatory waves, DM/density perturbations with precession sin(30°), and GW term.
 *          Supports dynamic variable updates for all parameters.
 * 
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: SMBHSgrAStar sgrA;
 *              Compute: double g = sgrA.compute_g_SgrA(t);
 * 
 * Key Features:
 *   - Default values from UQFF document, with approximations for all terms.
 *   - Units handled: B(t) converted to T (1 G = 10^-4 T); energy terms converted to effective acceleration.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_SgrA(r, t) with every term explicitly included.
 * 
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef SMBH_SGR_A_STAR_H
#define SMBH_SGR_A_STAR_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <limits>
#include <map>
#include <fstream>

class SMBHSgrAStar {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M_initial;       // Initial SMBH mass
    double r;               // Schwarzschild radius
    double H0;              // Hubble constant (s^-1)
    double B0_G;            // Initial magnetic field (G)
    double tau_B;           // B decay timescale (s)
    double B_crit;          // Critical B field (T)
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double v_surf;          // Surface velocity (arbitrary for BH)
    double f_TRZ;           // Time-reversal factor
    double M_dot_0;         // Initial mass accretion rate factor
    double tau_acc;         // Accretion timescale (s)
    double spin_factor;     // Spin factor (0.3)
    double tau_Omega;       // Omega decay timescale (s)

    // Additional parameters for full inclusion of terms
    double hbar;            // Reduced Planck's constant
    double t_Hubble;        // Hubble time (s)
    double delta_x;         // Position uncertainty (m)
    double delta_p;         // Momentum uncertainty (kg m/s)
    double integral_psi;    // Wavefunction integral approximation
    double rho_fluid;       // Fluid density (kg/m^3)
    double A_osc;           // Oscillatory amplitude (m/s^2)
    double k_osc;           // Wave number (1/m)
    double omega_osc;       // Angular frequency (rad/s)
    double x_pos;           // Position for oscillation (m)
    double t_Hubble_gyr;    // Hubble time in Gyr
    double M_DM_factor;     // Dark matter mass fraction
    double delta_rho_over_rho; // Density perturbation fraction
    double precession_angle_deg; // Precession angle (degrees)

    // UQFF upgrade members (Terms 10-12)
    double mu0;             // Vacuum permeability (H/m)
    double L0_W;            // Volumetric inductance factor (m^3)
    double tau_decay;       // Magnetic energy decay timescale (s)
    double D0_burst;        // QPO flare amplitude (m/s^2)
    double omega_D;         // QPO angular frequency (rad/s) — 20-min Sgr A* QPO
    double tau_D;           // Burst envelope decay timescale (s)
    double M_NSC;           // Nuclear star cluster mass (kg)
    double r_NSC;           // NSC distance from BH (m) ~3 pc

    // UQFF Buoyancy parameters (CP1/CP2/CP3 pipeline)
    double beta_i;     // Buoyancy coupling constant [CP3 canonical: 0.61, PAPER_198]
    double omega_g;    // Galactic rotation rate [CP1/CP2 canonical: 7.3e-16 rad/s]
    double U_UA;       // Unit charge aether parameter [CP1/CP2 canonical: 1e-11 C]

    // UQFF 2.0 framework
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 for initial M (will recompute with M(t))

    static constexpr double UQFF_PI = 3.14159265358979323846;

public:
    // Constructor with default UQFF values
    SMBHSgrAStar() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~SMBHSgrAStar() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        M_initial = 4.3e6 * 1.989e30;
        r = 1.27e10;
        H0 = 2.268e-18;                          // CP3 canonical
        B0_G = 1e4;  // G
        tau_B = 1e6 * 3.15576e7;                 // canonical s/yr
        B_crit = 1e11;  // T
        Lambda = 1.114e-52;                      // CP3 canonical
        c_light = 2.998e8;                       // CP3 canonical
        q_charge = 1.602e-19;
        v_surf = 1e6;  // Arbitrary
        f_TRZ = 0.1;
        M_dot_0 = 0.01;
        tau_acc = 9e9 * 3.15576e7;               // canonical s/yr
        spin_factor = 0.3;
        tau_Omega = 9e9 * 3.15576e7;

        // Full terms defaults (CP3-canonical)
        hbar = 1.0546e-34;
        t_Hubble = 1.0 / H0;                     // derived from H0 (CP3 canonical)
        t_Hubble_gyr = t_Hubble / 3.15576e16;
        delta_x = 1e-15;                         // nuclear scale for compact object
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
        rho_fluid = 1e-15;                       // CP3 canonical accretion disc
        A_osc = 1e12;                            // CP3 canonical
        k_osc = 2.0 * UQFF_PI / r;              // CP3 canonical (was 1.0/r)
        omega_osc = 2.0 * UQFF_PI * (c_light / r); // CP3 canonical
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;
        precession_angle_deg = 30.0;

        // UQFF upgrade members
        mu0       = 1.2566e-6;
        L0_W      = (4.0 / 3.0) * UQFF_PI * r * r * r;
        tau_decay = 1e8 * 3.15576e7;             // 100 Myr magnetic energy decay
        D0_burst  = 1e8;                         // QPO flare amplitude (m/s^2)
        omega_D   = 2.0 * UQFF_PI / 1200.0;     // 20-min QPO (Sgr A* canonical)
        tau_D     = 3.0 * 3.15576e7;             // 3-yr burst envelope decay
        M_NSC     = 2.5e7 * 1.989e30;           // Nuclear star cluster mass
        r_NSC     = 9.26e16;                     // 3 pc in metres
        logging_enabled = false;

        // UQFF buoyancy parameters (CP3/PAPER_198/CP1 canonical)
        beta_i  = 0.61;      // PAPER_198/CP3 canonical coupling [FUBiiTaxonomy, BETA_I]
        omega_g = 7.3e-16;   // Galactic rotation rate [CP1/CP2 canonical] (rad/s)
        U_UA    = 1e-11;     // Unit charge aether [CP1/CP2 canonical] (C)

        updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        ug1_base = (G * M_initial) / (r * r);
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M_initial") { M_initial = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "H0") { H0 = newValue; }
        else if (varName == "B0_G") { B0_G = newValue; }
        else if (varName == "tau_B") { tau_B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "v_surf") { v_surf = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "M_dot_0") { M_dot_0 = newValue; }
        else if (varName == "tau_acc") { tau_acc = newValue; }
        else if (varName == "spin_factor") { spin_factor = newValue; }
        else if (varName == "tau_Omega") { tau_Omega = newValue; }
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
        else if (varName == "precession_angle_deg") { precession_angle_deg = newValue; }
        else if (varName == "mu0")         { mu0 = newValue; }
        else if (varName == "L0_W")        { L0_W = newValue; }
        else if (varName == "tau_decay")   { tau_decay = newValue; }
        else if (varName == "D0_burst")    { D0_burst = newValue; }
        else if (varName == "omega_D")     { omega_D = newValue; }
        else if (varName == "tau_D")       { tau_D = newValue; }
        else if (varName == "M_NSC")       { M_NSC = newValue; }
        else if (varName == "r_NSC")       { r_NSC = newValue; }
        else if (varName == "beta_i")      { beta_i  = newValue; }
        else if (varName == "omega_g")     { omega_g = newValue; }
        else if (varName == "U_UA")        { U_UA    = newValue; }
        else if (varName == "logging_enabled") { logging_enabled = (newValue != 0.0); }
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
        else if (varName == "M_initial") return M_initial;
        else if (varName == "r") return r;
        else if (varName == "H0") return H0;
        else if (varName == "B0_G") return B0_G;
        else if (varName == "tau_B") return tau_B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "v_surf") return v_surf;
        else if (varName == "f_TRZ") return f_TRZ;
        else if (varName == "M_dot_0") return M_dot_0;
        else if (varName == "tau_acc") return tau_acc;
        else if (varName == "spin_factor") return spin_factor;
        else if (varName == "tau_Omega") return tau_Omega;
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
        else if (varName == "precession_angle_deg") return precession_angle_deg;
        else if (varName == "mu0")         return mu0;
        else if (varName == "L0_W")        return L0_W;
        else if (varName == "tau_decay")   return tau_decay;
        else if (varName == "D0_burst")    return D0_burst;
        else if (varName == "omega_D")     return omega_D;
        else if (varName == "tau_D")       return tau_D;
        else if (varName == "M_NSC")       return M_NSC;
        else if (varName == "r_NSC")       return r_NSC;
        else if (varName == "beta_i")      return beta_i;
        else if (varName == "omega_g")     return omega_g;
        else if (varName == "U_UA")        return U_UA;
        else if (varName == "logging_enabled") return logging_enabled ? 1.0 : 0.0;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    // M(t) computation
    double M_t(double t) const {
        double M_dot = M_dot_0 * exp(-t / tau_acc);
        return M_initial * (1 + M_dot);
    }

    // B(t) in T
    double B_t(double t) const {
        double B_G = B0_G * exp(-t / tau_B);
        return B_G * 1e-4;  // G to T
    }

    // Omega(t) computation
    double Omega_t(double t) const {
        double omega0 = spin_factor * c_light / r;
        return omega0 * exp(-t / tau_Omega);
    }

    // dOmega/dt computation
    double dOmega_dt(double t) const {
        double omega0 = spin_factor * c_light / r;
        return omega0 * (-1.0 / tau_Omega) * exp(-t / tau_Omega);
    }

    // Ug terms computation
    double compute_Ug(double Mt, double Bt) const {
        double Ug1 = (G * Mt) / (r * r);
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - Bt / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ);
    }

    // Volume computation for fluid
    double compute_V() const {
        return (4.0 / 3.0) * UQFF_PI * r * r * r;
    }

    // Main MUGE computation (includes ALL terms)
    double compute_g_SgrA(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Mt = M_t(t);
        double Bt = B_t(t);
        double dOdt = dOmega_dt(t);
        double ug1_t = (G * Mt) / (r * r);

        // Term 1: Base + H0 + B corrections
        double corr_H = 1 + H0 * t;
        double corr_B = 1 - Bt / B_crit;
        double term1 = ug1_t * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ
        double term2 = compute_Ug(Mt, Bt);

        // Term 3: Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // Term 4: EM (v x B, no scaling or UA here)
        double cross_vB = v_surf * Bt;  // Magnitude
        double em_base = q_charge * cross_vB / 1.673e-27;  // Acceleration
        double term4 = em_base;

        // Term 5: GW
        double gw_prefactor = (G * Mt * Mt) / (pow(c_light, 4) * r);
        double term5 = gw_prefactor * (dOdt * dOdt);

        // Quantum uncertainty term
        double sqrt_unc = sqrt(delta_x * delta_p);
        double term_q = (hbar / sqrt_unc) * integral_psi * (2.0 * UQFF_PI / t_Hubble);

        // Fluid term (effective acceleration)
        double V = compute_V();
        double term_fluid = (rho_fluid * V * ug1_t) / Mt;

        // Oscillatory term (AGN variability; Mode 2 removed — dimensionally inconsistent)
        double term_osc1 = 2.0 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double term_osc = term_osc1;

        // DM and density perturbation term with precession (converted to acceleration)
        double M_dm = Mt * M_DM_factor;
        double sin_prec = sin(precession_angle_deg * UQFF_PI / 180.0);
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * Mt / (r * r * r);
        double term_dm_force_like = (Mt + M_dm) * (pert1 + pert2 * sin_prec);
        double term_DM = term_dm_force_like / Mt;

        // Term 10: Magnetic stored energy gradient (B²·4πr²·exp(-t/τ) / (2μ₀·M(t)))
        double U_B    = (Bt * Bt) / (2.0 * mu0);  // J/m^3
        double term10 = U_B * 4.0 * UQFF_PI * r * r / Mt * exp(-t / tau_decay);

        // Term 11: QPO burst/flare modulation D(t) = D0·cos(ω_D·t)·exp(-t/τ_D)
        //          omega_D = 2π/1200s (20-min Sgr A* QPO canonical cadence)
        double term11 = D0_burst * cos(omega_D * t) * exp(-t / tau_D);

        // Term 12: Nuclear star cluster tidal gradient (2·G·M_NSC·r·r_NSC⁻³)
        //          M_NSC = 2.5×10⁷ M_sun, r_NSC = 3 pc (9.26×10¹⁶ m)
        double term12 = 2.0 * G * M_NSC * r / (r_NSC * r_NSC * r_NSC);

        // UQFF Buoyancy Terms (CP1/CP2/CP3 pipeline — previously missing from C++ modules)
        // Ubi (CP3 canonical line 1770): Static half-gravity buoyancy — dominant term
        // Ubi = 0.5 × (G·M(t)/r²) [FUBiiTaxonomyCompactObjectCalculator, PAPER_198]
        // Uses ug1_t (time-evolving mass) for accretion-consistent buoyancy
        double term_Ubi = 0.5 * ug1_t;

        // F_UBii (PAPER_198 compact object): -β_i × Ug_i × ω_g × (M(t)/r) × [UA] × cos(π·t)
        // Reference magnitude for SMBH: ~1e36 N [CP3 PAPER_198 taxonomy; BH tier]
        double term_F_UBii = -beta_i * ug1_t * omega_g * (Mt / r) * U_UA
                             * std::cos(UQFF_PI * t);

        // Ub_i (CP1 outer-frame): buoyancy from nuclear star cluster at 3 pc
        // Ub_i = -β_i × ug1_t × ω_g × (M_NSC / r_NSC) × [UA] × cos(π·t)
        // Sgr A*'s buoyancy is driven by the NSC embedding it (M_NSC = 2.5×10⁷ M_sun at 3 pc)
        // CP1 compute_buoyancy_regime: Sgr A* is NEGATIVE buoyancy (F_U_Bi_i = -8.31e211 N)
        double term_Ub_i = -beta_i * ug1_t * omega_g * (M_NSC / r_NSC) * U_UA
                           * std::cos(UQFF_PI * t);

        // FU diagnostic (CP3): FU = -(Ug_sum + Ubi) × g_base [logged only, not added again]
        double FU_diag = -(term2 + term_Ubi) * ug1_t;

        // Total g_SgrA (15-term MUGE)
        double g_total = term1 + term2 + term3 + term4 + term5 + term_q
                       + term_fluid + term_osc + term_DM + term10 + term11 + term12
                       + term_Ubi + term_F_UBii + term_Ub_i;

        if (logging_enabled) {
            std::cout << "[SMBHSgrAStar] t=" << t
                      << " T10=" << term10 << " T11=" << term11
                      << " T12=" << term12
                      << " tUbi=" << term_Ubi << " tF_UBii=" << term_F_UBii
                      << " tUb_i=" << term_Ub_i << " FU_diag=" << FU_diag
                      << " g=" << g_total << std::endl;
        }
        return g_total;
    }

    // Comprehensive debug/output method (UQFF 2.0 standard)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::scientific << std::setprecision(6);
        os << "=== SMBHSgrAStar Parameters ===" << std::endl;
        os << "G:                  " << G             << " m^3/kg/s^2" << std::endl;
        os << "M_initial:          " << M_initial      << " kg (" << M_initial/1.989e30 << " M_sun)" << std::endl;
        os << "r (Schwarzschild):  " << r              << " m" << std::endl;
        os << "H0:                 " << H0             << " s^-1 (CP3 canonical)" << std::endl;
        os << "Lambda:             " << Lambda         << " (CP3 canonical)" << std::endl;
        os << "c_light:            " << c_light        << " m/s" << std::endl;
        os << "B0_G:               " << B0_G           << " G (-> " << B0_G*1e-4 << " T)" << std::endl;
        os << "tau_B:              " << tau_B          << " s" << std::endl;
        os << "B_crit:             " << B_crit         << " T" << std::endl;
        os << "f_TRZ:              " << f_TRZ          << std::endl;
        os << "M_dot_0:            " << M_dot_0        << std::endl;
        os << "tau_acc:            " << tau_acc        << " s" << std::endl;
        os << "spin_factor:        " << spin_factor    << std::endl;
        os << "delta_x:            " << delta_x        << " m (nuclear scale)" << std::endl;
        os << "rho_fluid:          " << rho_fluid      << " kg/m^3 (CP3 canonical)" << std::endl;
        os << "A_osc:              " << A_osc          << " m/s^2" << std::endl;
        os << "k_osc:              " << k_osc          << " m^-1 (2pi/r)" << std::endl;
        os << "omega_osc:          " << omega_osc      << " rad/s (2pi*c/r)" << std::endl;
        os << "M_DM_factor:        " << M_DM_factor    << std::endl;
        os << "precession_angle:   " << precession_angle_deg << " deg" << std::endl;
        os << "mu0:                " << mu0            << " H/m" << std::endl;
        os << "L0_W:               " << L0_W           << " m^3" << std::endl;
        os << "tau_decay:          " << tau_decay      << " s (mag energy decay)" << std::endl;
        os << "D0_burst:           " << D0_burst       << " m/s^2 (Sgr A* QPO)" << std::endl;
        os << "omega_D:            " << omega_D        << " rad/s (2pi/1200 = 20-min QPO)" << std::endl;
        os << "tau_D:              " << tau_D          << " s (burst envelope)" << std::endl;
        os << "M_NSC:              " << M_NSC          << " kg (2.5e7 M_sun NSC)" << std::endl;
        os << "r_NSC:              " << r_NSC          << " m (3 pc)" << std::endl;
        os << "ug1_base:           " << ug1_base       << " m/s^2" << std::endl;
        os << "logging_enabled:    " << (logging_enabled ? "true" : "false") << std::endl;
        if (!dynamic_params.empty()) {
            os << "Dynamic params:     " << dynamic_params.size() << " registered" << std::endl;
            for (const auto& kv : dynamic_params)
                os << "  " << kv.first << " = " << kv.second << std::endl;
        }
    }

    // Example computation at t=4.5 Gyr (for testing)
    double exampleAt4_5Gyr() const {
        double t_example = 4.5e9 * 3.15576e7;  // canonical s/yr
        return compute_g_SgrA(t_example);
    }

    // -----------------------------------------------------------------------
    // UQFF 2.0 Framework
    // -----------------------------------------------------------------------
    void setEnableLogging(bool enable) { logging_enabled = enable; }
    bool isLoggingEnabled() const { return logging_enabled; }

    void setDynamicParameter(const std::string& key, double value) {
        dynamic_params[key] = value;
        if (logging_enabled)
            std::cout << "[SMBHSgrAStar] dyn_param " << key << " = " << value << std::endl;
    }

    double getDynamicParameter(const std::string& key, double defaultVal = 0.0) const {
        auto it = dynamic_params.find(key);
        return (it != dynamic_params.end()) ? it->second : defaultVal;
    }

    void exportState(const std::string& filename = "SMBH_SGR_A_STAR_state.txt") const {
        std::ofstream f(filename);
        if (!f.is_open()) {
            std::cerr << "[SMBHSgrAStar] exportState: cannot open " << filename << std::endl;
            return;
        }
        f << std::scientific << std::setprecision(15);
        f << "G="         << G         << "\n";
        f << "M_initial=" << M_initial << "\n";
        f << "r="         << r         << "\n";
        f << "H0="        << H0        << "\n";
        f << "Lambda="    << Lambda    << "\n";
        f << "c_light="   << c_light   << "\n";
        f << "B0_G="      << B0_G      << "\n";
        f << "tau_B="     << tau_B     << "\n";
        f << "B_crit="    << B_crit    << "\n";
        f << "f_TRZ="     << f_TRZ     << "\n";
        f << "M_dot_0="   << M_dot_0   << "\n";
        f << "tau_acc="   << tau_acc   << "\n";
        f << "delta_x="   << delta_x   << "\n";
        f << "rho_fluid=" << rho_fluid << "\n";
        f << "mu0="       << mu0       << "\n";
        f << "L0_W="      << L0_W      << "\n";
        f << "tau_decay=" << tau_decay << "\n";
        f << "D0_burst="  << D0_burst  << "\n";
        f << "omega_D="   << omega_D   << "\n";
        f << "tau_D="     << tau_D     << "\n";
        f << "M_NSC="     << M_NSC     << "\n";
        f << "r_NSC="     << r_NSC     << "\n";
        for (const auto& kv : dynamic_params)
            f << "dyn_" << kv.first << "=" << kv.second << "\n";
        f.close();
        if (logging_enabled)
            std::cout << "[SMBHSgrAStar] State exported to " << filename << std::endl;
    }

    // cross_validate<OtherModule>() — header-only comparison template
    template <typename OtherModule>
    double cross_validate(double t) const {
        OtherModule other;
        double g_this  = compute_g_SgrA(t);
        double g_other = other.compute_g_SgrA(t);
        double ratio   = (g_other != 0.0) ? g_this / g_other
                                          : std::numeric_limits<double>::quiet_NaN();
        if (logging_enabled)
            std::cout << "[SMBHSgrAStar] cross_validate: g_this=" << g_this
                      << " g_other=" << g_other << " ratio=" << ratio << std::endl;
        return ratio;
    }
};

#endif // SMBH_SGR_A_STAR_H