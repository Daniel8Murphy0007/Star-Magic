/**
 * ================================================================================================
 * Header: RingsOfRelativity.h
 * 
 * Description: C++ Module for "Rings of Relativity" (GAL-CLUS-022058s Einstein Ring) Class
 *              This is the eighth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on gravitational lensing ring evolution
 *              and gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 * 
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for "Rings of Relativity" evolution.
 *          Includes ALL terms: base gravity (static M), cosmic expansion (H(z)), magnetic correction (static B),
 *          lensing amplification L(t) (static), UQFF Ug components with f_TRZ, Lambda, quantum uncertainty,
 *          scaled EM with [UA], fluid dynamics, oscillatory waves, DM/density perturbations, and stellar
 *          wind feedback (pressure / density for acc, added for completeness). Supports dynamic variable updates.
 * 
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: RingsOfRelativity rings;
 *              Compute: double g = rings.compute_g_Rings(t);
 * 
 * Key Features:
 *   - Default values from UQFF document: M = 1e14 Msun, r = 3.086e20 m (~10 kpc), B = 1e-5 T,
 *     H_z for z=0.5 ≈ 2.42e-18 s^-1, L = (G*M)/(c^2*r) * 0.67, etc.
 *   - Units handled: Msun to kg, kpc to m; wind term as (rho * v_wind^2) / rho_fluid for acceleration.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_Rings(r, t) with every term explicitly included.
 * 
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef RINGS_OF_RELATIVITY_H
#define RINGS_OF_RELATIVITY_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <limits>
#include <map>
#include <fstream>

class RingsOfRelativity {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M;               // Lensing mass (kg)
    double r;               // Einstein radius (m)
    double Hz;              // Hubble parameter at z (s^-1)
    double B;               // Static magnetic field (T)
    double B_crit;          // Critical B field (T)
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double gas_v;           // Gas velocity for EM (m/s)
    double f_TRZ;           // Time-reversal factor
    double L_factor;        // Lensing factor (D_LS / D_S ≈ 0.67)
    double rho_vac_UA;      // UA vacuum density (J/m^3)
    double rho_vac_SCm;     // SCm vacuum density (J/m^3)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass for EM acceleration
    double z_lens;          // Lens redshift

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
    double rho_wind;        // Wind density (kg/m^3)
    double v_wind;          // Wind velocity (m/s)

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 = G*M/r^2
    double L_t;             // Cached lensing term

    // UQFF buoyancy pipeline parameters (CP3/PAPER_198 canonical)
    double beta_i;          // Buoyancy coupling constant (PAPER_198 = 0.61)
    double omega_g;         // Effective rotation rate (rad/s)
    double U_UA;            // Unit charge aether parameter (C)
    // External gravitational body: Virgo Cluster (~16.5 Mpc from reference frame)
    double M_VC;            // Virgo Cluster mass ~1.2e15 M_sun (kg)
    double r_VC;            // Distance to Virgo Cluster ~16.5 Mpc in meters

    // UQFF 2.0 Self-Expanding Framework
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

public:
    // Constructor with default UQFF values
    RingsOfRelativity() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~RingsOfRelativity() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M = 1e14 * M_sun;
        r = 3.086e20;
        z_lens = 0.5;
        double Hz_kms = 70 * sqrt(0.3 * pow(1 + z_lens, 3) + 0.7);  // km/s/Mpc
        Hz = (Hz_kms * 1000 / 3.086e19);  // s^-1
        B = 1e-5;
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        gas_v = 1e5;
        f_TRZ = 0.1;
        L_factor = 0.67;
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
        rho_fluid = 1e-21;
        A_osc = 1e-12;  // Small for lensing scale
        k_osc = 1.0 / r;
        omega_osc = 2 * M_PI / (r / c_light);
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;
        rho_wind = 1e-21;
        v_wind = 2e6;

        // UQFF buoyancy pipeline (CP3/PAPER_198 canonical)
        beta_i = 0.61;               // PAPER_198/CP3 buoyancy coupling constant
        omega_g = 7.3e-16;           // Effective rotation rate (rad/s)
        U_UA = 1e-11;                // Unit charge aether parameter (C)
        M_VC = 2.387e45;             // Virgo Cluster ~1.2e15 M_sun (external body, kg)
        r_VC = 5.09e23;              // Distance to Virgo Cluster ~16.5 Mpc in meters

        logging_enabled = false;
        updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        ug1_base = (G * M) / (r * r);
        L_t = ((G * M) / (pow(c_light, 2) * r)) * L_factor;
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M") { M = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "Hz") { Hz = newValue; }
        else if (varName == "B") { B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "gas_v") { gas_v = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "L_factor") { L_factor = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
        else if (varName == "proton_mass") { proton_mass = newValue; }
        else if (varName == "z_lens") { z_lens = newValue; }
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
        else if (varName == "rho_wind") { rho_wind = newValue; }
        else if (varName == "v_wind") { v_wind = newValue; }
        // UQFF buoyancy parameters
        else if (varName == "beta_i") { beta_i = newValue; }
        else if (varName == "omega_g") { omega_g = newValue; }
        else if (varName == "U_UA") { U_UA = newValue; }
        else if (varName == "M_VC") { M_VC = newValue; }
        else if (varName == "r_VC") { r_VC = newValue; }
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
        else if (varName == "B") return B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "gas_v") return gas_v;
        else if (varName == "f_TRZ") return f_TRZ;
        else if (varName == "L_factor") return L_factor;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "scale_EM") return scale_EM;
        else if (varName == "proton_mass") return proton_mass;
        else if (varName == "z_lens") return z_lens;
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
        else if (varName == "rho_wind") return rho_wind;
        else if (varName == "v_wind") return v_wind;
        // UQFF buoyancy parameters
        else if (varName == "beta_i") return beta_i;
        else if (varName == "omega_g") return omega_g;
        else if (varName == "U_UA") return U_UA;
        else if (varName == "M_VC") return M_VC;
        else if (varName == "r_VC") return r_VC;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return 0.0;
        }
    }

    // Ug terms computation
    double compute_Ug(double /*Mt*/) const {  // Mt static as M
        double Ug1 = ug1_base;
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
    double compute_g_Rings(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        // Term 1: Base + Hz + B + L corrections
        double corr_H = 1 + Hz * t;
        double corr_B = 1 - B / B_crit;
        double corr_L = 1 + L_t;
        double term1 = ug1_base * corr_H * corr_B * corr_L;

        // Term 2: UQFF Ug with f_TRZ
        double term2 = compute_Ug(0);  // No Mt variation

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
        double term_fluid = (rho_fluid * V * ug1_base) / M;

        // Oscillatory terms (real parts)
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        double term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term (converted to acceleration)
        double M_dm = M * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M / (r * r * r);
        double term_dm_force_like = (M + M_dm) * (pert1 + pert2);
        double term_DM = term_dm_force_like / M;

        // Stellar wind feedback term (pressure / density for acceleration)
        double wind_pressure = rho_wind * v_wind * v_wind;
        double term_wind = wind_pressure / rho_fluid;

        // CP3/PAPER_198 Tier-1: static half-gravity Ubi (0.5 * g_base)
        // Note: M is static for galaxy cluster lensing system — use ug1_base throughout
        double term_Ubi = 0.5 * ug1_base;
        // CP3/PAPER_198 Tier-2: dynamic compact F_UBii — cosine modulation over Einstein radius
        double term_F_UBii = -beta_i * ug1_base * omega_g * (M / r) * U_UA * cos(M_PI * t);
        // CP1 Tier-3: outer-frame Ub_i via Virgo Cluster (~16.5 Mpc, nearest rich cluster)
        double term_Ub_i = -beta_i * ug1_base * omega_g * (M_VC / r_VC) * U_UA * cos(M_PI * t);
        double FU_diag = -(term2 + term_Ubi) * ug1_base;  // diagnostic only (not added to total)

        // Total g_Rings (12 terms: 9 original + 3 UQFF buoyancy tiers)
        double result = term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind
                      + term_Ubi + term_F_UBii + term_Ub_i;
        if (logging_enabled)
            std::cout << "[LOG] compute_g_Rings(t=" << t << ") = " << result
                      << "  [t1=" << term1 << " t2=" << term2 << " t3=" << term3
                      << " t4=" << term4 << " tq=" << term_q << " tf=" << term_fluid
                      << " tosc=" << term_osc << " tdm=" << term_DM << " twind=" << term_wind
                      << " tUbi=" << term_Ubi << " tF_UBii=" << term_F_UBii
                      << " tUb_i=" << term_Ub_i << " FU_diag=" << FU_diag << "]" << std::endl;
        return result;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "Rings of Relativity Parameters:" << std::endl;
        os << "G: " << G << ", M: " << M << ", r: " << r << std::endl;
        os << "Hz: " << Hz << ", B: " << B << ", B_crit: " << B_crit << std::endl;
        os << "f_TRZ: " << f_TRZ << ", L_t: " << L_t << ", L_factor: " << L_factor << std::endl;
        os << "rho_fluid: " << rho_fluid << ", rho_wind: " << rho_wind << ", v_wind: " << v_wind << std::endl;
        os << "gas_v: " << gas_v << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        os << "ug1_base: " << ug1_base << std::endl;
        os << "beta_i: " << beta_i << ", omega_g: " << omega_g << ", U_UA: " << U_UA << std::endl;
        os << "M_VC: " << M_VC << " (Virgo Cluster ext body), r_VC: " << r_VC << " (~16.5 Mpc)" << std::endl;
    }

    // Example computation at t=5 Gyr (for testing)
    double exampleAt5Gyr() const {
        double t_example = 5e9 * 3.156e7;
        return compute_g_Rings(t_example);
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

    void exportState(const std::string& filename = "RINGS_OF_RELATIVITY_state.txt") const {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cerr << "Error: Cannot open '" << filename << "' for state export." << std::endl;
            return;
        }
        ofs << "# RingsOfRelativity state export\n";
        ofs << std::scientific << std::setprecision(10);
        const char* names[] = {
            "G","M","r","Hz","B","B_crit","Lambda","c_light",
            "q_charge","gas_v","f_TRZ","L_factor","z_lens",
            "rho_vac_UA","rho_vac_SCm","scale_EM","proton_mass",
            "hbar","t_Hubble","t_Hubble_gyr","delta_x","delta_p","integral_psi",
            "rho_fluid","A_osc","k_osc","omega_osc","x_pos","M_DM_factor","delta_rho_over_rho",
            "rho_wind","v_wind",
            "beta_i","omega_g","U_UA","M_VC","r_VC"
        };
        for (const char* n : names)
            ofs << n << " = " << getVariable(n) << "\n";
        for (const auto& kv : dynamic_params)
            ofs << "dynamic." << kv.first << " = " << kv.second << "\n";
        if (logging_enabled)
            std::cout << "[LOG] State exported to: " << filename << std::endl;
    }

    // Cross-validate against PillarsOfCreation at t=5 Gyr — fractional |gRings − gPillars| / |gRings|
    template<typename PillarsOfCreationT>
    double cross_validate(const PillarsOfCreationT& pillars, double t_years = 5.0e9) const {
        double t_s     = t_years * 3.156e7;
        double g_rings = compute_g_Rings(t_s);
        double g_pil   = pillars.compute_g_Pillars(t_s);
        double frac    = (g_rings != 0.0)
                         ? std::abs(g_pil - g_rings) / std::abs(g_rings)
                         : std::numeric_limits<double>::quiet_NaN();
        if (logging_enabled)
            std::cout << "[XVAL] Rings g(" << t_years << "yr)=" << g_rings
                      << "  Pillars g=" << g_pil
                      << "  frac_diff=" << frac << std::endl;
        return frac;
    }
};

#endif // RINGS_OF_RELATIVITY_H