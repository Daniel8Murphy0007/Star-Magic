/**
 * ================================================================================================
 * Header: GalaxyNGC1792.h
 * 
 * Description: C++ Module for NGC 1792 ("The Stellar Forge") Starburst Galaxy Class
 *              This is the nineteenth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on starburst galaxy evolution and
 *              gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 * 
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for NGC 1792 evolution.
 *          Includes ALL terms: base gravity with star formation growth M(t), cosmic expansion (H(z)),
 *          magnetic correction (static B), UQFF Ug components with f_TRZ, Lambda, quantum uncertainty,
 *          scaled EM with [UA], fluid dynamics, oscillatory waves, DM/density perturbations, and
 *          supernova feedback (pressure / density for acc). Supports dynamic variable updates.
 * 
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: GalaxyNGC1792 ngc1792;
 *              Compute: double g = ngc1792.compute_g_NGC1792(t);
 * 
 * Key Features:
 *   - Default values from UQFF document: M0 = 1e10 Msun, r = 7.569e20 m (80k ly), z = 0.0095,
 *     Hz â‰ˆ 2.19e-18 s^-1, SFR_factor = 10 / 1e10 (normalized), tau_SF = 100 Myr, B = 1e-5 T,
 *     rho_wind = 1e-21 kg/m^3, v_wind = 2e6 m/s.
 *   - Units handled: Msun to kg, ly to m; feedback term as (rho * v_wind^2) / rho_fluid for acceleration.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_NGC1792(r, t) with every term explicitly included.
 * 
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef GALAXY_NGC_1792_H
#define GALAXY_NGC_1792_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <limits>
#include <map>
#include <fstream>
#include <typeinfo>

// WOLFRAM_TERM macros â€” canonical Wolfram Language representation of each MUGE term
#define WOLFRAM_TERM_NGC1792_BASE  "G*M0*(1+SFR*Exp[-t/tauSF])/(r^2)*(1+Hz*t)*(1-B/Bcrit)"
#define WOLFRAM_TERM_NGC1792_UQFF  "(G*M0/r^2)*(1+1-B/Bcrit)*(1+fTRZ)"
#define WOLFRAM_TERM_NGC1792_BUOY  "-betai*G*M/(r^2)*omegaG*(M/r)*U_UA*Cos[Pi*t]"
#define WOLFRAM_TERM_NGC1792_SN    "rhoWind*vWind^2/rhoFluid"
#define WOLFRAM_TERM_NGC1792_OSC   "2*Aosc*Cos[kosc*x]*Cos[omegaOsc*t] + (2*Pi/tHubble)*Aosc*Cos[kosc*x - omegaOsc*t]"

class GalaxyNGC1792 {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M0;              // Initial mass (kg)
    double r;               // Radius (m)
    double Hz;              // Hubble parameter at z (s^-1)
    double B;               // Static magnetic field (T)
    double B_crit;          // Critical B field (T)
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double gas_v;           // Gas velocity for EM (m/s)
    double f_TRZ;           // Time-reversal factor
    double SFR_factor;      // Star formation rate factor (dimensionless)
    double tau_SF;          // Star formation timescale (s)
    double rho_wind;        // Wind density (kg/m^3)
    double v_wind;          // Wind velocity (m/s)
    double rho_fluid;       // Fluid density (kg/m^3)
    double rho_vac_UA;      // UA vacuum density (J/m^3)
    double rho_vac_SCm;     // SCm vacuum density (J/m^3)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass for EM acceleration
    double z_gal;           // Galaxy redshift

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

    // UQFF 2.0 â€” 3-tier buoyancy parameters (CP3/PAPER_198 standard)
    double beta_i;          // Buoyancy coupling constant (calibrated 0.61)
    double omega_g;         // Gravitational angular frequency (7.3e-16 rad/s)
    double U_UA;            // [UA] vacuum coupling constant (1e-11)
    double M_Fornax;        // Fornax Cluster mass (kg) â€” outer-frame external body
    double r_Fornax;        // NGC 1792 â†’ Fornax Cluster distance (m, ~20 Mpc)
    // UQFF 2.0 â€” Self-expanding framework
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

public:
    // Constructor with default UQFF values
    GalaxyNGC1792() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~GalaxyNGC1792() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M0 = 1e10 * M_sun;
        double ly_to_m = 9.461e15;
        r = 80000.0 * ly_to_m;
        z_gal = 0.0095;
        double Hz_kms = 70 * sqrt(0.3 * pow(1 + z_gal, 3) + 0.7);  // km/s/Mpc
        Hz = (Hz_kms * 1000 / 3.086e19);  // s^-1
        B = 1e-5;
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        gas_v = 1e5;
        f_TRZ = 0.1;
        SFR_factor = 10.0 / (1e10);  // Normalized SFR
        tau_SF = 100e6 * 3.15576e7;
        rho_wind = 1e-21;
        v_wind = 2e6;
        rho_fluid = 1e-21;
        rho_vac_UA = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        scale_EM = 1e-12;
        proton_mass = 1.673e-27;

        // Full terms defaults
        hbar = 1.0546e-34;
        t_Hubble = 13.8e9 * 3.15576e7;
        t_Hubble_gyr = 13.8;
        delta_x = 1e-10;
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
        A_osc = 1e-10;
        k_osc = 1.0 / r;
        omega_osc = 2 * M_PI / (r / c_light);
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;

        // UQFF 2.0 â€” 3-tier buoyancy initialization
        beta_i   = 0.61;
        omega_g  = 7.3e-16;
        U_UA     = 1e-11;
        M_Fornax = 1.393e44;  // 7e13 M_sun Ã— 1.989e30 kg (Fornax Cluster)
        r_Fornax = 6.17e23;   // ~20 Mpc = 20 Ã— 3.086e22 m
        logging_enabled = false;

        updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        ug1_base = B * r * G * M0  /* DPM: mu_s * grad(M_s/r) */;
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
        else if (varName == "rho_wind") { rho_wind = newValue; }
        else if (varName == "v_wind") { v_wind = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
        else if (varName == "proton_mass") { proton_mass = newValue; }
        else if (varName == "z_gal") { z_gal = newValue; }
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
        // UQFF 2.0 buoyancy parameters
        else if (varName == "beta_i")   { beta_i = newValue; }
        else if (varName == "omega_g")  { omega_g = newValue; }
        else if (varName == "U_UA")     { U_UA = newValue; }
        else if (varName == "M_Fornax") { M_Fornax = newValue; }
        else if (varName == "r_Fornax") { r_Fornax = newValue; }
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
        else if (varName == "rho_wind") return rho_wind;
        else if (varName == "v_wind") return v_wind;
        else if (varName == "rho_fluid") return rho_fluid;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "scale_EM") return scale_EM;
        else if (varName == "proton_mass") return proton_mass;
        else if (varName == "z_gal") return z_gal;
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
        // UQFF 2.0 buoyancy parameters
        else if (varName == "beta_i")   return beta_i;
        else if (varName == "omega_g")  return omega_g;
        else if (varName == "U_UA")     return U_UA;
        else if (varName == "M_Fornax") return M_Fornax;
        else if (varName == "r_Fornax") return r_Fornax;
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

    // Ug terms computation
    double compute_Ug(double Mt) const {
        // DPM-seeded: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        // DPM-seeded: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        double Ug1 = B * r * G * Mt  /* DPM: mu_s * grad(M_s/r) */;
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - B / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ);
    }

    // Volume computation for fluid
    double compute_V() const {
        return (4.0 / 3.0) * M_PI * r * r * r;
    }

    // Main MUGE computation (includes ALL terms)
    double compute_g_NGC1792(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Mt = M_t(t);
        double ug1_t = B * r * G * Mt  /* DPM: mu_s * grad(M_s/r) */;

        // Term 1: Base + Hz + B corrections
        double corr_H = 1 + Hz * t;
        double corr_B = 1 - B / B_crit;
        double term1 = ug1_t * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ
        double term2 = compute_Ug(Mt);

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
        // PAPER_268 fix: use t_Hubble (seconds) not t_Hubble_gyr (dimensionless Gyr)
        double term_osc2 = (2 * M_PI / t_Hubble) * A_osc * cos(arg);
        double term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term (converted to acceleration)
        double M_dm = Mt * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * B * G * Mt  /* DPM tidal */;
        double term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        double term_DM = term_dm_force_like / Mt;

        // Supernova feedback term (pressure / density for acceleration)
        double wind_pressure = rho_wind * v_wind * v_wind;
        double term_feedback = wind_pressure / rho_fluid;

        // === UQFF 2.0: 3-Tier Buoyancy (CP3/PAPER_198) ===
        // Tier 1: Static half-gravity Ubi (CP3 PAPER_198)
        double term_Ubi    = 0.5 * ug1_t;
        // Tier 2: Dynamic compact cos modulation over M(t) â€” F_UBii (PAPER_198)
        double term_F_UBii = -beta_i * ug1_t * omega_g * (Mt / r) * U_UA * cos(M_PI * t);
        // Tier 3: Outer-frame via Fornax Cluster external body â€” Ub_i (CP1)
        double term_Ub_i   = -beta_i * ug1_t * omega_g * (M_Fornax / r_Fornax) * U_UA * cos(M_PI * t);

        // FU diagnostic (logged only â€” not summed into MUGE)
        double FU_diag = -(term2 + term_Ubi) * ug1_t;

        if (logging_enabled) {
            log("[NGC1792] term1=" + std::to_string(term1)
              + " term2=" + std::to_string(term2)
              + " term3=" + std::to_string(term3)
              + " term4=" + std::to_string(term4)
              + " term_q=" + std::to_string(term_q)
              + " term_fluid=" + std::to_string(term_fluid)
              + " term_osc=" + std::to_string(term_osc)
              + " term_DM=" + std::to_string(term_DM)
              + " term_feedback=" + std::to_string(term_feedback)
              + " term_Ubi=" + std::to_string(term_Ubi)
              + " term_F_UBii=" + std::to_string(term_F_UBii)
              + " term_Ub_i=" + std::to_string(term_Ub_i)
              + " FU_diag=" + std::to_string(FU_diag));
        }

        // 12-term MUGE total (9 original terms + 3 UQFF 2.0 buoyancy tiers)
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc
             + term_DM + term_feedback + term_Ubi + term_F_UBii + term_Ub_i;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "NGC 1792 Parameters:" << std::endl;
        os << "G: " << G << ", M0: " << M0 << ", r: " << r << std::endl;
        os << "Hz: " << Hz << ", B: " << B << ", B_crit: " << B_crit << std::endl;
        os << "f_TRZ: " << f_TRZ << ", SFR_factor: " << SFR_factor << ", tau_SF: " << tau_SF << std::endl;
        os << "rho_fluid: " << rho_fluid << ", rho_wind: " << rho_wind << ", v_wind: " << v_wind << std::endl;
        os << "gas_v: " << gas_v << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        os << "ug1_base: " << ug1_base << std::endl;
        // UQFF 2.0 buoyancy parameters
        os << "beta_i: " << beta_i << ", omega_g: " << omega_g << ", U_UA: " << U_UA << std::endl;
        os << "Ext body (Fornax Cluster): M_Fornax=" << M_Fornax << " kg, r_Fornax=" << r_Fornax << " m" << std::endl;
        os << "logging_enabled: " << (logging_enabled ? "true" : "false") << std::endl;
    }

    // Example computation at t=50 Myr (for testing)
    double exampleAt50Myr() const {
        double t_example = 50e6 * 3.15576e7;  // canonical seconds per year
        return compute_g_NGC1792(t_example);
    }

    // === UQFF 2.0 Self-Expanding Framework Methods ===
    void setEnableLogging(bool enable) { logging_enabled = enable; }
    bool getLoggingEnabled() const { return logging_enabled; }

    void setDynamicParameter(const std::string& key, double value) {
        dynamic_params[key] = value;
    }

    double getDynamicParameter(const std::string& key) const {
        auto it = dynamic_params.find(key);
        if (it != dynamic_params.end()) return it->second;
        return std::numeric_limits<double>::quiet_NaN();
    }

    void exportState(const std::string& filename = "GalaxyNGC1792_state.txt") const {
        std::ofstream f(filename);
        if (!f.is_open()) return;
        f << "# GalaxyNGC1792 UQFF 2.0 State Export\n";
        f << "G=" << G << "\nM0=" << M0 << "\nr=" << r << "\nbeta_i=" << beta_i
          << "\nomega_g=" << omega_g << "\nU_UA=" << U_UA
          << "\nM_Fornax=" << M_Fornax << "\nr_Fornax=" << r_Fornax << "\n";
        for (const auto& kv : dynamic_params)
            f << kv.first << "=" << kv.second << "\n";
    }

    // Cross-validation template (compare two GalaxyNGC1792 instances at time t)
    template<typename T>
    static double cross_validate(const T& a, const T& b, double t) {
        double ga = a.compute_g_NGC1792(t);
        double gb = b.compute_g_NGC1792(t);
        return (ga != 0.0) ? std::abs((ga - gb) / ga) : std::abs(ga - gb);
    }

private:
    void log(const std::string& msg) const {
        std::cout << "[LOG][NGC1792] " << msg << std::endl;
    }
};

#endif // GALAXY_NGC_1792_H
