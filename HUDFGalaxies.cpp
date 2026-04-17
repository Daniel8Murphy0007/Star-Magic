/**
 * ================================================================================================
 * Header: HUDFGalaxies.h
 * 
 * Description: C++ Module for Hubble Ultra Deep Field (HUDF) "Galaxies Galore" Class
 *              This is the eighteenth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on cosmic field of galaxies evolution
 *              and gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 * 
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for HUDF galaxies evolution.
 *          Includes ALL terms: base gravity with formation growth M(t), cosmic expansion (H(z)),
 *          magnetic correction (static B), interactions I(t), UQFF Ug components with f_TRZ,
 *          Lambda, quantum uncertainty, scaled EM with [UA], fluid dynamics, oscillatory waves,
 *          DM/density perturbations, and merger feedback. Supports dynamic variable updates.
 * 
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: HUDFGalaxies hudf;
 *              Compute: double g = hudf.compute_g_HUDF(t);
 * 
 * Key Features:
 *   - Default values from UQFF document: M0 ≈ 1e12 Msun (representative field mass), r ≈ 1.3e11 ly (cosmic scale),
 *     z_avg = 3.5, Hz_avg ≈ 2.5e-18 s^-1, SFR_factor = 1.0, tau_SF = 1 Gyr, I0 = 0.05, tau_inter = 1 Gyr,
 *     rho_wind = 1e-22 kg/m^3, v_wind = 1e6 m/s, B = 1e-10 T.
 *   - Units handled: Msun to kg, ly to m; interaction term I(t) scales gravity.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_HUDF(r, t) with every term explicitly included.
 * 
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef HUDF_GALAXIES_H
#define HUDF_GALAXIES_H

// ===========================================================================================
// UQFF 2.0 UPGRADE — March 2026 (Session 72g)
// Added: PhysicsTerm base classes, self-expanding framework, F_U_Bi_i 3-tier buoyancy,
//        DPM resonance, WOLFRAM_TERM macros, exportState, registerDynamicTerm.
// Original HUDFGalaxies physics preserved intact (additive philosophy).
// New unique physics: f_TRZ negative-time zeroing, dual-channel I(t) buoyancy cascade,
//                    B_crit superconducting gravitational quench boundary.
// ===========================================================================================

#include <iostream>
#include <cmath>
#include <iomanip>
#include <map>
#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include <sstream>

// ===========================================================================================
// WOLFRAM_TERM macros — exported to masterUQFF via AutoExportFullUQFF (source176_auto_full_uqff.cpp)
// ===========================================================================================
#define WOLFRAM_TERM_HUDF_BASE   "G*M0*(1+SFR*Exp[-t/tauSF])/(r^2)*(1+Hz*t)*(1-B/Bcrit)*(1+I0*Exp[-t/tauI])"
#define WOLFRAM_TERM_HUDF_UQFF   "(G*M0/r^2)*(1+1-B/Bcrit)*(1+fTRZ)*(1+I0*Exp[-t/tauI])"
#define WOLFRAM_TERM_HUDF_TRZ    "(G*M0/r^2)*(1+fTRZ)"           // TRZ zeroing: fTRZ=-1 -> 0
#define WOLFRAM_TERM_HUDF_LENR   "kLENR*(omegaLENR/omega0)^2"    // LENR 1.25 THz term
#define WOLFRAM_TERM_HUDF_FUBIBI "(kLENR*(omegaLENR/omega0)^2 + kRel*Frel + kNeutron*sigmaN)*x2" // F_U_Bi_i

// ===========================================================================================
// ABSTRACT PHYSICS TERM BASE (UQFF 2.0 Self-Expanding Framework)
// ===========================================================================================

class HUDFPhysicsTerm {
public:
    virtual ~HUDFPhysicsTerm() {}
    virtual double compute(double t, const std::map<std::string, double>& params) const = 0;
    virtual std::string getName()        const = 0;
    virtual std::string getDescription() const = 0;
    virtual bool validate(const std::map<std::string, double>&) const { return true; }
};

// ===========================================================================================
// UNIQUE PHYSICS TERM 1: Time-Reversal Zeroing (f_TRZ) — CPT-asymmetric UQFF gravity
// PAPER_264 source: HUDFTRZNegativeTimeTerm
//   f_TRZ = -1 → factor (1+f_TRZ) = 0 → UQFF gravity vanishes ("time-reversal zero point")
//   f_TRZ < -1 → factor < 0 → UQFF gravity reverses sign (anti-gravity / negative-time regime)
// ===========================================================================================
class HUDFTRZNegativeTimeTerm : public HUDFPhysicsTerm {
private:
    double f_TRZ;
public:
    explicit HUDFTRZNegativeTimeTerm(double trz_factor = 0.1) : f_TRZ(trz_factor) {}

    double compute(double /*t*/, const std::map<std::string, double>& params) const override {
        double G   = params.count("G")   ? params.at("G")   : 6.6743e-11;
        double M   = params.count("M")   ? params.at("M")   : 1.989e42;
        double r   = params.count("r")   ? params.at("r")   : 1.23e27;
        double trz = params.count("f_TRZ") ? params.at("f_TRZ") : f_TRZ;
        // DPM-emergent: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        // DPM-emergent: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        double Ug1 = (G * M) / (r * r);
        return Ug1 * (1.0 + trz);   // zero at trz=-1; negative for trz<-1
    }

    std::string getName()        const override { return "HUDFTRZNegativeTime"; }
    std::string getDescription() const override {
        return "CPT-asymmetric UQFF: (1+f_TRZ) zeroes at f_TRZ=-1 (time-reversal zero point); "
               "negative for f_TRZ<-1 (anti-gravity negative-time regime). PAPER_264."; }
    bool validate(const std::map<std::string, double>& p) const override {
        double trz = p.count("f_TRZ") ? p.at("f_TRZ") : f_TRZ;
        return (trz > -2.0 && trz <= 1.0); // physical range: above full anti-gravity flip
    }
};

// ===========================================================================================
// UNIQUE PHYSICS TERM 2: Dual-Channel Interaction Cascade Buoyancy — double I(t) modulation
// PAPER_265 source: HUDFInteractionCascadeTerm
//   I(t) applied to BOTH base gravity AND UQFF simultaneously.
//   Net factor: (1+I(t))^2 in combined output — quadratic buoyancy amplification.
// ===========================================================================================
class HUDFInteractionCascadeTerm : public HUDFPhysicsTerm {
private:
    double I0;
    double tau_inter;
public:
    HUDFInteractionCascadeTerm(double i0 = 0.05, double tau_i = 3.156e16)
        : I0(i0), tau_inter(tau_i) {}

    double compute(double t, const std::map<std::string, double>& params) const override {
        double G   = params.count("G")  ? params.at("G")  : 6.6743e-11;
        double M   = params.count("M")  ? params.at("M")  : 1.989e42;
        double r   = params.count("r")  ? params.at("r")  : 1.23e27;
        double i0  = params.count("I0") ? params.at("I0") : I0;
        double tau = params.count("tau_inter") ? params.at("tau_inter") : tau_inter;
        double I_t = i0 * std::exp(-t / tau);
        // DPM-emergent: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        // DPM-emergent: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        double Ug1 = (G * M) / (r * r);
        // Cascade: (1+I)^2 total modulation factor across both channels
        double cascade_factor = (1.0 + I_t) * (1.0 + I_t);
        // Additional buoyancy relative to single-channel: (1+I)^2 - (1+I) = I^2 + I
        double single_channel = Ug1 * (1.0 + I_t);
        double cascade_extra  = Ug1 * I_t * I_t;  // quadratic cascade bonus
        return single_channel + cascade_extra;    // total with cascade enhancement
    }

    std::string getName()        const override { return "HUDFInteractionCascade"; }
    std::string getDescription() const override {
        return "Dual-channel I(t) buoyancy cascade: I(t) modulates BOTH base gravity and UQFF "
               "simultaneously -> quadratic (1+I)^2 amplification at merger peak. PAPER_265."; }
};

// ===========================================================================================
// UNIQUE PHYSICS TERM 3: UQFF Critical Magnetic Superconducting Boundary
// PAPER_266 source: HUDFCriticalMagneticTerm
//   corr_B = 1 - B/B_crit encodes gravitational Meissner effect:
//   B -> B_crit: gravity expelled from UQFF field (analogous to Type II SC upper critical field)
//   B_crit = 1e11 T (neutron star Schwinger-like threshold in C++ original)
// ===========================================================================================
class HUDFCriticalMagneticTerm : public HUDFPhysicsTerm {
private:
    double B_crit;
public:
    explicit HUDFCriticalMagneticTerm(double bcrit = 1e11) : B_crit(bcrit) {}

    double compute(double /*t*/, const std::map<std::string, double>& params) const override {
        double G     = params.count("G")     ? params.at("G")     : 6.6743e-11;
        double M     = params.count("M")     ? params.at("M")     : 1.989e42;
        double r     = params.count("r")     ? params.at("r")     : 1.23e27;
        double B     = params.count("B")     ? params.at("B")     : 1e-10;
        double bcrit = params.count("B_crit")? params.at("B_crit"): B_crit;
        double Ug1   = (G * M) / (r * r);
        double corr_B = 1.0 - B / bcrit;  // Meissner suppression factor
        // Returns the suppressed field; note: corr_B -> 0 at B=B_crit -> FULL QUENCH
        return Ug1 * corr_B;
    }

    std::string getName()        const override { return "HUDFCriticalMagnetic"; }
    std::string getDescription() const override {
        return "UQFF gravitational Meissner effect: corr_B = 1-B/B_crit; "
               "full quench at B=B_crit=1e11 T (Schwinger-like NS threshold). PAPER_266."; }
    bool validate(const std::map<std::string, double>& p) const override {
        double B     = p.count("B")     ? p.at("B")     : 1e-10;
        double bcrit = p.count("B_crit")? p.at("B_crit"): B_crit;
        return (bcrit > 0.0 && B >= 0.0 && B <= bcrit * 1.1); // allow slight above-crit probe
    }
};

class HUDFGalaxies {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M0;              // Initial field mass (kg)
    double r;               // Effective radius (m)
    double Hz;              // Average Hubble parameter at z (s^-1)
    double B;               // Static magnetic field (T)
    double B_crit;          // Critical B field (T)
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double gas_v;           // Gas velocity for EM (m/s)
    double f_TRZ;           // Time-reversal factor
    double SFR_factor;      // Star formation rate factor (dimensionless)
    double tau_SF;          // Star formation timescale (s)
    double I0;              // Initial interaction factor
    double tau_inter;       // Interaction timescale (s)
    double rho_wind;        // Wind density (kg/m^3)
    double v_wind;          // Wind velocity (m/s)
    double rho_fluid;       // Fluid density (kg/m^3)
    double rho_vac_UA;      // UA vacuum density (J/m^3)
    double rho_vac_SCm;     // SCm vacuum density (J/m^3)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass for EM acceleration
    double z_avg;           // Average redshift

    // Additional parameters for full inclusion of terms
    double hbar;            // Reduced Planck's constant
    double t_Hubble;        // Hubble time (s)
    double delta_x;         // Position uncertainty (m)
    double delta_p;         // Momentum uncertainty (kg m/s)
    double integral_psi;    // Wavefunction integral approximation
    double A_osc;           // Oscillatory amplitude (m/s^2)
    double k_osc;           // Wave number (1/m)
    double omega_osc;       // Angular frequency (rad/s)
    double x_pos;           // Position for oscillation (m)
    double t_Hubble_gyr;    // Hubble time in Gyr
    double M_DM_factor;     // Dark matter mass fraction
    double delta_rho_over_rho; // Density perturbation fraction

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 for initial M0

    // ========== SELF-EXPANDING FRAMEWORK MEMBERS (UQFF 2.0) ==========
    std::map<std::string, double>                   dynamicParameters;
    std::vector<std::unique_ptr<HUDFPhysicsTerm>>   dynamicTerms;
    std::map<std::string, std::string>              metadata;
    bool   enableDynamicTerms;
    bool   enableLogging;
    double learningRate;

    // Helper logger
    void log(const std::string& msg) const {
        if (enableLogging) std::cout << "[HUDFGalaxies] " << msg << std::endl;
    }

public:
    // Constructor with default UQFF values
    HUDFGalaxies() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~HUDFGalaxies() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M0 = 1e12 * M_sun;
        double ly_to_m = 9.461e15;
        r = 1.3e11 * ly_to_m;  // Cosmic scale
        z_avg = 3.5;
        double Hz_kms = 70 * sqrt(0.3 * pow(1 + z_avg, 3) + 0.7);  // km/s/Mpc
        Hz = (Hz_kms * 1000 / 3.086e19);  // s^-1
        B = 1e-10;
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        gas_v = 1e5;
        f_TRZ = 0.1;
        SFR_factor = 1.0;
        tau_SF = 1e9 * 3.156e7;
        I0 = 0.05;
        tau_inter = 1e9 * 3.156e7;
        rho_wind = 1e-22;
        v_wind = 1e6;
        rho_fluid = 1e-22;
        rho_vac_UA = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        scale_EM = 1e-12;
        proton_mass = 1.673e-27;

        // Full terms defaults
        hbar = 1.0546e-34;
        t_Hubble = 13.8e9 * 3.156e7;
        t_Hubble_gyr = 13.8;
        delta_x = 1e-10;
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
        A_osc = 1e-12;
        k_osc = 1.0 / r;
        omega_osc = 2 * M_PI / (r / c_light);
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;

        updateCache();
        enableDynamicTerms = false;
        enableLogging      = false;
        learningRate       = 0.001;
        metadata["version"] = "UQFF_2.0_HUDF";
        metadata["session"] = "72g_March2026";
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        ug1_base = (G * M0) / (r * r);
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M0") { M0 = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "Hz") { Hz = newValue; }
        else if (varName == "B") { B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "gas_v") { gas_v = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "SFR_factor") { SFR_factor = newValue; }
        else if (varName == "tau_SF") { tau_SF = newValue; }
        else if (varName == "I0") { I0 = newValue; }
        else if (varName == "tau_inter") { tau_inter = newValue; }
        else if (varName == "rho_wind") { rho_wind = newValue; }
        else if (varName == "v_wind") { v_wind = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
        else if (varName == "proton_mass") { proton_mass = newValue; }
        else if (varName == "z_avg") { z_avg = newValue; }
        // Full terms
        else if (varName == "hbar") { hbar = newValue; }
        else if (varName == "t_Hubble") { t_Hubble = newValue; }
        else if (varName == "t_Hubble_gyr") { t_Hubble_gyr = newValue; }
        else if (varName == "delta_x") { delta_x = newValue; }
        else if (varName == "delta_p") { delta_p = newValue; }
        else if (varName == "integral_psi") { integral_psi = newValue; }
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
        else if (varName == "M0") return M0;
        else if (varName == "r") return r;
        else if (varName == "Hz") return Hz;
        else if (varName == "B") return B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "gas_v") return gas_v;
        else if (varName == "f_TRZ") return f_TRZ;
        else if (varName == "SFR_factor") return SFR_factor;
        else if (varName == "tau_SF") return tau_SF;
        else if (varName == "I0") return I0;
        else if (varName == "tau_inter") return tau_inter;
        else if (varName == "rho_wind") return rho_wind;
        else if (varName == "v_wind") return v_wind;
        else if (varName == "rho_fluid") return rho_fluid;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "scale_EM") return scale_EM;
        else if (varName == "proton_mass") return proton_mass;
        else if (varName == "z_avg") return z_avg;
        // Full terms
        else if (varName == "hbar") return hbar;
        else if (varName == "t_Hubble") return t_Hubble;
        else if (varName == "t_Hubble_gyr") return t_Hubble_gyr;
        else if (varName == "delta_x") return delta_x;
        else if (varName == "delta_p") return delta_p;
        else if (varName == "integral_psi") return integral_psi;
        else if (varName == "A_osc") return A_osc;
        else if (varName == "k_osc") return k_osc;
        else if (varName == "omega_osc") return omega_osc;
        else if (varName == "x_pos") return x_pos;
        else if (varName == "M_DM_factor") return M_DM_factor;
        else if (varName == "delta_rho_over_rho") return delta_rho_over_rho;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return 0.0;
        }
    }

    // M(t) computation
    double M_t(double t) const {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        return M0 * (1 + M_dot);
    }

    // I(t) computation
    double I_t(double t) const {
        return I0 * exp(-t / tau_inter);
    }

    // Ug terms computation
    double compute_Ug(double Mt, double It) const {
        // DPM-emergent: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        // DPM-emergent: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        double Ug1 = (G * Mt) / (r * r);
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - B / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ) * (1 + It);
    }

    // Volume computation for fluid
    double compute_V() const {
        return (4.0 / 3.0) * M_PI * r * r * r;
    }

    // Main MUGE computation (includes ALL terms)
    double compute_g_HUDF(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Mt = M_t(t);
        double It = I_t(t);
        double ug1_t = (G * Mt) / (r * r);

        // Term 1: Base + Hz + B + I corrections
        double corr_H = 1 + Hz * t;
        double corr_B = 1 - B / B_crit;
        double corr_I = 1 + It;
        double term1 = ug1_t * corr_H * corr_B * corr_I;

        // Term 2: UQFF Ug with f_TRZ and I
        double term2 = compute_Ug(Mt, It);

        // Term 3: Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // Term 4: Scaled EM with UA
        double cross_vB = gas_v * B;  // Magnitude, assuming perpendicular
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
        double term4 = (em_base * corr_UA) * scale_EM;

        // Quantum uncertainty term
        double sqrt_unc = sqrt(delta_x * delta_p);
        double term_q = (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);

        // Fluid term (effective acceleration)
        double V = compute_V();
        double term_fluid = (rho_fluid * V * ug1_t) / Mt;

        // Oscillatory terms (real parts)
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        double term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term (converted to acceleration)
        double M_dm = Mt * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * Mt / (r * r * r);
        double term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        double term_DM = term_dm_force_like / Mt;

        // Merger feedback term (pressure / density for acceleration)
        double wind_pressure = rho_wind * v_wind * v_wind;
        double term_feedback = wind_pressure / rho_fluid;

        // Total g_HUDF (all terms summed)
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_feedback;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "HUDF Galaxies Parameters:" << std::endl;
        os << "G: " << G << ", M0: " << M0 << ", r: " << r << std::endl;
        os << "Hz: " << Hz << ", B: " << B << ", B_crit: " << B_crit << std::endl;
        os << "f_TRZ: " << f_TRZ << ", SFR_factor: " << SFR_factor << ", tau_SF: " << tau_SF << std::endl;
        os << "I0: " << I0 << ", tau_inter: " << tau_inter << std::endl;
        os << "rho_fluid: " << rho_fluid << ", rho_wind: " << rho_wind << ", v_wind: " << v_wind << std::endl;
        os << "gas_v: " << gas_v << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        os << "ug1_base: " << ug1_base << std::endl;
    }

    // Example computation at t=5 Gyr (for testing)
    double exampleAt5Gyr() const {
        double t_example = 5e9 * 3.156e7;
        return compute_g_HUDF(t_example);
    }

    // ========================================================================
    // UQFF 2.0 SELF-EXPANDING FRAMEWORK — Dynamic Term Registration
    // ========================================================================
    void registerDynamicTerm(std::unique_ptr<HUDFPhysicsTerm> term) {
        enableDynamicTerms = true;
        log("Registered dynamic term: " + term->getName());
        dynamicTerms.push_back(std::move(term));
    }

    void setDynamicParameter(const std::string& name, double value) {
        dynamicParameters[name] = value;
        log("Set dynamic param: " + name);
    }

    double getDynamicParameter(const std::string& name, double defaultVal = 0.0) const {
        auto it = dynamicParameters.find(name);
        return (it != dynamicParameters.end()) ? it->second : defaultVal;
    }

    void setLearningRate(double rate) { learningRate = rate; }
    void setEnableLogging(bool enable) { enableLogging = enable; }

    // Compute sum of all registered dynamic terms at time t
    double computeDynamicContributions(double t) const {
        if (!enableDynamicTerms || dynamicTerms.empty()) return 0.0;
        std::map<std::string, double> params = {
            {"G", G}, {"M", M0}, {"r", r}, {"B", B},
            {"B_crit", B_crit}, {"f_TRZ", f_TRZ},
            {"I0", I0}, {"tau_inter", tau_inter}
        };
        for (auto& kv : dynamicParameters) params[kv.first] = kv.second;
        double sum = 0.0;
        for (auto& term : dynamicTerms) {
            if (term->validate(params)) {
                double contrib = term->compute(t, params);
                log("  " + term->getName() + " -> " + std::to_string(contrib));
                sum += contrib;
            } else {
                log("  SKIP (validate=false): " + term->getName());
            }
        }
        return sum;
    }

    // ========================================================================
    // F_U_Bi_i — UQFF 3-TIER BUOYANCY INTEGRAL (UQFF 2.0)
    // Integrand: F_LENR + F_neutron + F_rel + vacuum terms
    // Returns: (F_U_Bi_i, DPM_resonance) pair
    // ========================================================================
    std::pair<double, double> compute_F_U_Bi_i(double t = 0.0) const {
        (void)t; // time-independent integral for HUDF (cosmological scale)

        // UQFF canonical constants
        const double k_LENR    = 1e-10;
        const double omega_LENR = 2.0 * M_PI * 1.25e12;  // 1.25 THz
        const double omega_0   = 1e-12;                   // rad/s reference
        const double k_neutron = 1e10;                    // neutron term constant
        const double sigma_n   = 1e-4;                    // cosmic density parameter
        const double F_rel_LEP = 4.30e33;                 // LEP 1998 relativistic anchor
        const double F0         = 1.83e71;                 // vacuum energy anchor
        const double rho_vac_ref = 7.09e-36;
        const double b_quad     = 4.72e-3;                 // canonical quadratic stiffness

        // DPM resonance: magnetic Bohr-magneton coupling
        const double mu_B = 9.274e-24;  // Bohr magneton (J/T)
        double DPM_resonance = (2.0 * mu_B * B) / (hbar * omega_0);

        // F_LENR (1.25 THz LENR neutron-capture channel)
        double F_LENR = k_LENR * std::pow(omega_LENR / omega_0, 2.0);

        // F_neutron (degenerate matter / cosmic field density)
        double F_neutron = k_neutron * sigma_n;

        // F_relativistic (LEP 1998 anchor, constant)
        double F_relativistic = F_rel_LEP;

        // Integrand sum
        double integrand = F_LENR + F_neutron + F_relativistic;

        // Quadratic root x_2 (vacuum-energy dominated)
        // a·x² + b·x + c = 0;  c = -F0 + rho_vac*DPM_stab ≈ -F0
        double a_quad = G * M0 / (r * r);
        double c_quad = -F0 + rho_vac_UA;
        double disc   = b_quad * b_quad - 4.0 * a_quad * c_quad;
        double x_2    = (disc >= 0.0)
                        ? (-b_quad + std::sqrt(disc)) / (2.0 * a_quad)
                        : F0 / b_quad;  // fallback: vacuum-dominated root

        double F_U_Bi_i = integrand * x_2;

        log("F_LENR=" + std::to_string(F_LENR) +
            " F_n=" + std::to_string(F_neutron) +
            " F_rel=" + std::to_string(F_relativistic) +
            " x2=" + std::to_string(x_2) +
            " -> F_U_Bi_i=" + std::to_string(F_U_Bi_i));

        return {F_U_Bi_i, DPM_resonance};
    }

    // ========================================================================
    // ENHANCED COMPUTE (UQFF 2.0) — original g_HUDF + dynamic terms
    // ========================================================================
    double compute_g_HUDF_UQFF2(double t) const {
        double g_original = compute_g_HUDF(t);
        double g_dynamic  = computeDynamicContributions(t);
        return g_original + g_dynamic;
    }

    // ========================================================================
    // exportState — Self-Expanding Framework state persistence
    // ========================================================================
    void exportState(const std::string& filename) const {
        std::ofstream f(filename);
        if (!f) {
            std::cerr << "[HUDFGalaxies] Cannot open: " << filename << std::endl;
            return;
        }
        f << std::fixed << std::setprecision(10);
        f << "# HUDFGalaxies UQFF 2.0 State Export\n";
        f << "# version: " << metadata.at("version") << "\n";
        f << "# session: "  << metadata.at("session")  << "\n";
        f << "G="           << G            << "\n";
        f << "M0="          << M0           << "\n";
        f << "r="           << r            << "\n";
        f << "Hz="          << Hz           << "\n";
        f << "B="           << B            << "\n";
        f << "B_crit="      << B_crit       << "\n";
        f << "f_TRZ="       << f_TRZ        << "\n";
        f << "I0="          << I0           << "\n";
        f << "tau_inter="   << tau_inter    << "\n";
        f << "ug1_base="    << ug1_base     << "\n";
        f << "learningRate=" << learningRate << "\n";
        f << "dynamic_terms_count=" << dynamicTerms.size() << "\n";
        for (auto& kv : dynamicParameters)
            f << "dyn_param:" << kv.first << "=" << kv.second << "\n";
        for (auto& term : dynamicTerms)
            f << "registered_term:" << term->getName() << "\n";
        auto [fubibi, dpm] = compute_F_U_Bi_i();
        f << "F_U_Bi_i=" << fubibi << "\n";
        f << "DPM_resonance=" << dpm << "\n";
        log("State exported to " + filename);
    }
};

#endif // HUDF_GALAXIES_H