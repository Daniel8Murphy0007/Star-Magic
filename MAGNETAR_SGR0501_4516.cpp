/**
 * ================================================================================================
 * Header: MagnetarSGR0501_4516.h
 * 
 * Description: C++ Module for SGR 0501+4516 Magnetar Class
 *              This is the first module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on magnetar evolution and gravity
 *              equations derived from Hubble datasets, high-energy lab simulations, and UQFF
 *              refinements (dated May 08, 2025, updated for integration on October 08, 2025).
 * 
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for SGR 0501+4516 magnetar
 *          evolution. Supports dynamic variable updates (addition, subtraction, direct set) for
 *          all parameters to enable parametric studies and framework advancements.
 * 
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: MagnetarSGR0501_4516 mag;
 *              Compute: double g = mag.compute_g_Magnetar(t);
 * 
 * Key Features:
 *   - Default values from UQFF document (M = 2.785e30 kg, r = 20e3 m, etc.).
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_Magnetar(r, t) with all terms (base gravity, cosmic expansion, magnetic decay,
 *     UQFF Ug components with f_TRZ, Lambda, scaled EM, GW).
 *   - Negligible terms (quantum, fluid, DM perturbations) omitted as per derivation.
 * 
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef MAGNETAR_SGR0501_4516_H
#define MAGNETAR_SGR0501_4516_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <limits>
#include <map>
#include <fstream>

// Portable PI constant (M_PI is non-standard; not guaranteed by MSVC without _USE_MATH_DEFINES)
static constexpr double UQFF_PI = 3.14159265358979323846;

class MagnetarSGR0501_4516 {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M;               // Magnetar mass
    double r;               // Radius
    double H0;              // Hubble constant (s^-1)
    double B0;              // Initial magnetic field
    double tau_B;           // B decay timescale (s)
    double B_crit;          // Critical B field
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double v_surf;          // Surface velocity
    double f_TRZ;           // Time-reversal factor
    double rho_vac_UA;      // UA vacuum density
    double rho_vac_SCm;     // SCm vacuum density
    double P_init;          // Initial rotation period (s)
    double tau_Omega;       // Omega decay timescale (s)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass for EM acceleration

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 = G*M/r^2
    double I_ns;            // Cached neutron star moment of inertia = 0.4*M*r^2

    // UQFF 2.0 Self-Expanding Framework
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

public:
    // Constructor with default UQFF values
    MagnetarSGR0501_4516() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~MagnetarSGR0501_4516() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        M = 1.4 * 1.989e30;
        r = 20e3;
        H0 = 2.184e-18;
        B0 = 1.9e14;           // SGR 0501+4516 observed dipole surface field (T) — Swift/RXTE timing
        tau_B = 4000 * 3.156e7;
        B_crit = 4.4e13;       // Schwinger quantum critical field (T) — standard magnetar reference
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        v_surf = 1e6;
        f_TRZ = 0.1;
        rho_vac_UA = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        P_init = 5.0;
        tau_Omega = 10000 * 3.156e7;
        scale_EM = 1e-12;
        proton_mass = 1.673e-27;
        I_ns = 0.0;            // will be set by updateCache()
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
        ug1_base = (G * M) / (r * r);
        I_ns     = 0.4 * M * r * r;  // Neutron star moment of inertia (solid sphere)
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M") { M = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "H0") { H0 = newValue; }
        else if (varName == "B0") { B0 = newValue; }
        else if (varName == "tau_B") { tau_B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "v_surf") { v_surf = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "P_init") { P_init = newValue; }
        else if (varName == "tau_Omega") { tau_Omega = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
        else if (varName == "proton_mass") { proton_mass = newValue; }
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return false;
        }
        updateCache();
        return true;
    }

    // Addition method for variables
    bool addToVariable(const std::string& varName, double delta) {
        if (!setVariable(varName, getVariable(varName) + delta)) {
            return false;
        }
        updateCache();
        return true;
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
        else if (varName == "H0") return H0;
        else if (varName == "B0") return B0;
        else if (varName == "tau_B") return tau_B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "v_surf") return v_surf;
        else if (varName == "f_TRZ") return f_TRZ;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "P_init") return P_init;
        else if (varName == "tau_Omega") return tau_Omega;
        else if (varName == "scale_EM") return scale_EM;
        else if (varName == "proton_mass") return proton_mass;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return std::numeric_limits<double>::quiet_NaN();  // NaN prevents silent corruption
        }
    }

    // B(t) computation
    double B_t(double t) const {
        return B0 * exp(-t / tau_B);
    }

    // Omega(t) computation
    double Omega_t(double t) const {
        if (P_init <= 0.0) { std::cerr << "Warning: P_init <= 0." << std::endl; return 0.0; }
        return (2.0 * UQFF_PI / P_init) * exp(-t / tau_Omega);
    }

    // dOmega/dt computation
    double dOmega_dt(double t) const {
        if (P_init <= 0.0) { return 0.0; }
        double omega0 = 2.0 * UQFF_PI / P_init;
        return omega0 * (-1.0 / tau_Omega) * exp(-t / tau_Omega);
    }

    // Ug terms computation (requires t for Ug3 rotational coupling)
    double compute_Ug(double Bt, double t) const {
        // Ug1: base Newtonian gravity
        double Ug1 = ug1_base;
        // Ug2: charge-reactivity — EM acceleration from surface charge carrier in magnetar B field
        double Ug2 = (q_charge * v_surf * Bt) / (proton_mass * r);
        // Ug3: string-rotation coupling — angular momentum flux contribution to effective gravity
        double omega_t = Omega_t(t);
        double dOdt    = dOmega_dt(t);
        double Ug3     = r * omega_t * std::abs(dOdt);
        // Ug4: magnetic field suppression/reversal of base gravity (super-critical when B > B_crit)
        double Ug4 = Ug1 * (1.0 - Bt / B_crit);
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1.0 + f_TRZ);
    }

    // Main MUGE computation
    double compute_g_Magnetar(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Bt = B_t(t);
        double dOdt = dOmega_dt(t);

        // Term 1: Base + H0 + B corrections
        double corr_H = 1 + H0 * t;
        double corr_B = 1 - Bt / B_crit;
        double term1 = ug1_base * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ (Ug1+Ug2+Ug3+Ug4)
        double term2 = compute_Ug(Bt, t);

        // Term 3: Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // Term 4: Scaled EM
        double cross_vB = v_surf * Bt;  // Magnitude, assuming perpendicular
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
        double term4 = (em_base * corr_UA) * scale_EM;

        // Term 5: GW back-reaction (spin-down power → gravity correction)
        // g_GW = I_ns * Omega(t) * |dOmega/dt| / (M * c * r)
        //      = 0.4 * r * Omega(t) * |dOmega/dt| / c_light   [units: m/s^2]
        double omega_t = Omega_t(t);
        double term5   = (I_ns * omega_t * std::abs(dOdt)) / (M * c_light * r);

        double result = term1 + term2 + term3 + term4 + term5;
        if (logging_enabled)
            std::cout << "[LOG] compute_g_Magnetar(t=" << t << ") = " << result
                      << "  [t1=" << term1 << " t2=" << term2 << " t3=" << term3
                      << " t4=" << term4 << " t5=" << term5 << "]" << std::endl;
        return result;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::scientific << std::setprecision(6);
        os << "=== SGR 0501+4516 MUGE Parameters ==" << std::endl;
        os << "G:           " << G            << " m^3 kg^-1 s^-2" << std::endl;
        os << "M:           " << M            << " kg  (1.4 M_sun)" << std::endl;
        os << "r:           " << r            << " m" << std::endl;
        os << "H0:          " << H0           << " s^-1" << std::endl;
        os << "B0:          " << B0           << " T  (observed dipole field)" << std::endl;
        os << "tau_B:       " << tau_B        << " s" << std::endl;
        os << "B_crit:      " << B_crit       << " T  (Schwinger field)" << std::endl;
        os << "Lambda:      " << Lambda       << " m^-2" << std::endl;
        os << "c_light:     " << c_light      << " m/s" << std::endl;
        os << "q_charge:    " << q_charge     << " C" << std::endl;
        os << "v_surf:      " << v_surf       << " m/s" << std::endl;
        os << "f_TRZ:       " << f_TRZ        << std::endl;
        os << "rho_vac_UA:  " << rho_vac_UA   << " kg/m^3" << std::endl;
        os << "rho_vac_SCm: " << rho_vac_SCm  << " kg/m^3" << std::endl;
        os << "P_init:      " << P_init       << " s" << std::endl;
        os << "tau_Omega:   " << tau_Omega    << " s" << std::endl;
        os << "scale_EM:    " << scale_EM     << std::endl;
        os << "proton_mass: " << proton_mass  << " kg" << std::endl;
        os << "--- Cached ---" << std::endl;
        os << "ug1_base:    " << ug1_base     << " m/s^2" << std::endl;
        os << "I_ns:        " << I_ns         << " kg m^2" << std::endl;
    }

    // Example computation at t=5000 years (for testing)
    double exampleAt5000Years() const {
        double t_example = 5000 * 3.156e7;
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
        ofs << "# MagnetarSGR0501_4516 state export\n";
        ofs << std::scientific << std::setprecision(10);
        const char* names[] = {
            "G","M","r","H0","B0","tau_B","B_crit","Lambda","c_light",
            "q_charge","v_surf","f_TRZ","rho_vac_UA","rho_vac_SCm",
            "P_init","tau_Omega","scale_EM","proton_mass"
        };
        for (const char* n : names)
            ofs << n << " = " << getVariable(n) << "\n";
        for (const auto& kv : dynamic_params)
            ofs << "dynamic." << kv.first << " = " << kv.second << "\n";
        if (logging_enabled)
            std::cout << "[LOG] State exported to: " << filename << std::endl;
    }
};

#endif // MAGNETAR_SGR0501_4516_H