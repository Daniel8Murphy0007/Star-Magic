/**
 * ================================================================================================
 * Header: UQFFSource10.h
 * 
 * Description: C++ Module for UQFF Source10 Text Module (Catalogue of Equations & Variables)
 *              This is the Source10 primary text module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, aggregating all general equations, variables, and solutions
 *              from previous documents (e.g., magnetars, galaxies, nebulae). Derived from Hubble datasets,
 *              high-energy lab simulations, and UQFF refinements (dated May 09, 2025, updated October 08, 2025).
 *              Installed as the first primary text module into Source12.cpp (main architecture, e.g., Star Magic)
 *              via aliases (#include and using declarations for previous modules).
 * 
 * Purpose: Central catalogue class for UQFF core components, buoyancy forces, resonance, and long-form calculations.
 *          Supports dynamic variable updates and integration with prior modules (e.g., MagnetarSGR0501_4516).
 *          Computes F_U_Bi_i, g_UQFF(r,t), and example solutions (e.g., Eta Carinae).
 * 
 * Integration: Designed for inclusion in Source12.cpp (main). Aliases installed via #include and using namespace.
 *              Instantiate: UQFFSource10 source10;
 *              Compute: double f = source10.compute_F_U_Bi_i(t);
 * 
 * Key Features:
 *   - Aggregates all variables (e.g., F_U_Bi_i terms, g_H = 1.252e46, neutron_factor).
 *   - Long-form calculations preserved as methods (e.g., DPM_resonance).
 *   - Aliases: Includes previous headers (e.g., #include "MagnetarSGR0501_4516.h") for UQFF inheritance.
 *   - Advancement: Unifies lab (Colman-Gillespie, Sweet, Kozima) to cosmic scales.
 * 
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef UQFF_SOURCE10_H
#define UQFF_SOURCE10_H

// Installed Aliases for Integration into Source12.cpp (Main Architecture, e.g., Star Magic)
// Primary includes for UQFF modules aggregated by Source10 (first 10 representative)
#include "MagnetarSGR0501_4516.h"
#include "MagnetarSGR1745_2900.h"
#include "SMBHSgrAStar.h"
#include "StarbirthTapestry.h"
#include "Westerlund2.h"
#include "PillarsOfCreation.h"
#include "RingsOfRelativity.h"
#include "GalaxyNGC2525.h"
#include "NGC3603.h"
#include "BubbleNebula.h"
#include "GalaxyNGC1792.h"
// Extend: #include "<ModuleName>.h" for all 500+ as needed in Source12.cpp

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <random>
#include <chrono>

// UQFF 2.0 — Wolfram Term Macros (auto-registration with Wolfram KB)
#define WOLFRAM_TERM_SOURCE10_BASE "UQFFSource10:F_U_Bi_i=integrand*x2+LENR*act*exp(-t/tau)+DE+resonance*neutron+rel*(1+f_TRZ)"
#define WOLFRAM_TERM_SOURCE10_UQFF "UQFFSource10:g_UQFF(r,t)=Sum[Ug1i+Ug2i+Ug3i+Ug4i,{i,1,26}]+Lambda*c^2/3+hbar/sqrt(dx*dp)*psi*(2Pi/t_H)+3TierBuoyancy"
#define WOLFRAM_TERM_SOURCE10_DPM  "UQFFSource10:DPM_resonance=g_H*mu_B*B0/(hbar*omega0)*2.82e-56; g_H=1.252e46 cosmicOrbitalAmplifier [PAPER_270]"
#define WOLFRAM_TERM_SOURCE10_THZ  "UQFFSource10:F_thz_shock=k_thz*(omega_thz/omega_0)^2*neutron_factor*conduit_scale; THz_gate=1.2 [PAPER_271]"
#define WOLFRAM_TERM_SOURCE10_VAC  "UQFFSource10:F_vac_rep=G*delta_rho_vac*M*v; k_vac=G=6.674e-11 gravitational-vacuum-drag [PAPER_272]"

using namespace std;

// Forward declarations for UQFF core types (from thread summaries)
struct UQFFCore {
    double F_U_Bi_i;  // Buoyancy force
    double integrand; // Integral component
    double x_2;       // Position factor
    // ... (full struct from doc)
};

struct VacuumRepulsion {
    double F_vac_rep;
    double k_vac;
    double delta_rho_vac;
    // ...
};

class UQFFSource10 {
private:
    // Key Dialogue Summary Sections (captured from thread as member variables with comments)
    // 1. UQFF Core: Buoyancy F_U_Bi_i = integrand * x_2, with terms for LENR, activation, DE, resonance, neutron, rel.
    double F_U_Bi_i;        // Buoyancy force (N)
    double integrand;       // Integral term
    double x_2;             // x^2 factor
    double LENR_term;       // LENR contribution
    double activation_term; // Activation energy
    double DE_term;         // Dark energy
    double resonance_term;  // Resonance (THz)
    double neutron_term;    // Neutron factor
    double rel_term;        // Relativistic

    // 2. Vacuum Repulsion: Analogy to surface tension spike/drop; F_vac_rep = k_vac * Δρ_vac * M * v.
    double F_vac_rep;
    double k_vac;
    double delta_rho_vac;
    double M_vac;           // Mass
    double v_vac;           // Velocity

    // 3. Tail Star Formation: 26 layers Um with THz comm; F_thz_shock = k_thz * (ω_thz / ω_0)^2 * neutron_factor * conduit_scale.
    double F_thz_shock;
    double k_thz;
    double omega_thz;
    double omega_0;
    double neutron_factor;  // 1=stable, 0=unstable
    double conduit_scale;

    // 4. Conduit: H + H2O abundance → COx; F_conduit = k_conduit * (H_abundance * water_state) * neutron_factor.
    double F_conduit;
    double k_conduit;
    double H_abundance;
    double water_state;     // 1 for incompressible/stable

    // 5. Spooky Action: Quantum string/wave; F_spooky = k_spooky * (string_wave / ω_0).
    double F_spooky;
    double k_spooky;
    double string_wave;

    // 6-10. Additional (neutron_factor, water_state, push-pull balance, systems, predictions as above)

    // From Triadic Clone: Compressed UQFF eq g(r,t) = sum_{i=1 to 26} (Ug1_i + Ug2_i + Ug3_i + Ug4_i)
    vector<double> Ug1_vec; // 26 layers
    vector<double> Ug2_vec;
    vector<double> Ug3_vec;
    vector<double> Ug4_vec;
    double E_DPM;           // (hbar * c / r_i^2) * Q_i * [SCm]_i
    double R_t;             // sum cos terms for resonance

    // Catalogue Variables (all from documents, e.g., g_H = 1.252e46)
    double g_H;             // Hydrogen UQFF orbital g-factor (1.252e46) [PAPER_270]
    double mu_B;            // Bohr magneton (9.274e-24 J/T)
    double B0;              // Magnetic field (T)
    double h_planck;        // Reduced Planck hbar (1.0546e-34 J·s)
    double omega_0_base;    // Base UQFF frequency (rad/s)

    // Physics constants required for compute_g_UQFF()
    double Lambda;          // Cosmological constant (1.114e-52 m^-2)
    double c_light;         // Speed of light (2.998e8 m/s)
    double hbar;            // Reduced Planck = h_planck (J·s)
    double delta_x;         // Position uncertainty (m)
    double delta_p;         // Momentum uncertainty (kg·m/s)
    double integral_psi;    // Wavefunction normalization (dimensionless)
    double t_Hubble;        // Hubble time (4.352e17 s)

    // F_U_Bi_i time-evolution variables
    double tau_SF;          // Star-formation decay timescale (s)
    double f_TRZ;           // TRZ phase coupling factor

    // UQFF 2.0 — 3-tier buoyancy members (CP3/PAPER_198 canonical)
    double beta_i;          // Buoyancy coupling constant (0.61)
    double omega_g;         // Gravitational buoyancy frequency (rad/s)
    double U_UA;            // [UA] unit analytic factor
    double M_EtaCar;        // Eta Carinae primary mass (kg) — Tier 2 compact body
    double r_EtaCar;        // Eta Carinae distance (m)
    double M_GC;            // Galactic Center Sgr A* mass (kg) — Tier 3 outer-frame
    double r_GC;            // Galactic Center distance (m)

    // UQFF 2.0 — Self-Expanding Framework
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;
    mutable std::mt19937 rng;
    mutable std::uniform_real_distribution<double> dis;

    // Computed caches
    double DPM_resonance;   // Resonance energy density (J/m^3)

public:
    // Constructor initializes defaults from catalogue
    UQFFSource10()
        : rng(static_cast<unsigned>(std::chrono::steady_clock::now().time_since_epoch().count()))
        , dis(0.0, 1.0)
    {
        initializeCatalogue();
    }

    ~UQFFSource10() {}

    // Load key=value scaling overrides from a config file into dynamic_params
    // Format: one "key=value" per line; lines beginning with '#' are comments
    void loadConfig(const std::string& config_file) {
        if (config_file.empty()) return;
        std::ifstream file(config_file);
        if (!file.is_open()) {
            if (logging_enabled) log("loadConfig: cannot open '" + config_file + "'");
            return;
        }
        std::string line;
        int loaded = 0;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            size_t eq = line.find('=');
            if (eq == std::string::npos) continue;
            std::string key = line.substr(0, eq);
            try {
                double val = std::stod(line.substr(eq + 1));
                dynamic_params[key] = val;
                ++loaded;
            } catch (...) {}
        }
        if (logging_enabled) {
            std::ostringstream oss;
            oss << "loadConfig: loaded " << loaded << " params from '" << config_file << "'";
            log(oss.str());
        }
        updateCache();
    }

    // Initialization from document catalogue
    void initializeCatalogue() {
        // UQFF Core defaults (Eta Carinae benchmark values)
        F_U_Bi_i        = 2.11e208;    // Eta Carinae benchmark buoyancy force (N)
        integrand       = 1.56e36;
        x_2             = 1.35e172;
        LENR_term       = 1.0e12;      // Colman-Gillespie LENR contribution
        activation_term = 1.0e-10;     // Activation energy scaling
        DE_term         = 1.0e-9;      // Dark energy contribution
        resonance_term  = 1.2e3;       // THz resonance amplitude
        neutron_term    = 1.0e49;      // Kozima neutron-drop factor (NS regime)
        rel_term        = 1.0e5;       // Sweet relativistic correction

        // Vacuum Repulsion (k_vac = G — gravitational vacuum drag, PAPER_272)
        F_vac_rep       = 1.23e45;
        k_vac           = 6.674e-11;   // Newton's constant G (PAPER_272 discovery)
        delta_rho_vac   = 1.0e-9;     // kg/m^3 vacuum density gradient
        M_vac           = 1.989e30;   // 1 M_sun reference
        v_vac           = 1.0e6;      // m/s reference velocity

        // Tail Star Formation / THz shock (Colman-Gillespie 1.25 THz gate, PAPER_271)
        F_thz_shock     = 4.56e78;
        k_thz           = 1.38e-23;   // Boltzmann constant as THz coupling
        omega_thz       = 1.2e12;     // rad/s (≈1.25 THz Colman-Gillespie gate)
        omega_0         = 1.0e12;     // rad/s base (THz gate ratio = 1.2)
        neutron_factor  = 1.0;        // Kozima gate: 1=stable/open, 0=unstable/closed
        conduit_scale   = 1.0e12;     // conduit amplification scale

        // Conduit: H + H2O → COx (dual gate PAPER_271)
        F_conduit       = 3.45e67;
        k_conduit       = 8.99e9;     // Coulomb constant as conduit coupling
        H_abundance     = 0.74;       // Cosmic hydrogen mass fraction
        water_state     = 1.0;        // Fluid gate: 1=incompressible/open, 0=closed

        // Spooky Action (quantum string / DPM coherence)
        F_spooky        = 2.71e89;
        k_spooky        = 1.0546e-34; // k_spooky = hbar
        string_wave     = 5.0e14;     // Optical string frequency (rad/s)

        // Triadic 26-layer MUGE structure (canonical base values)
        Ug1_vec.resize(26, 4.645e11); // Ug1 magnetic dipole base
        Ug2_vec.resize(26, 0.0);      // Ug2 charge-reactivity (zero unless set)
        Ug3_vec.resize(26, 0.0);      // Ug3 string rotation (zero unless set)
        Ug4_vec.resize(26, 4.645e11); // Ug4 vacuum concentration (= Ug1 base)
        E_DPM           = 3.11e9;     // DPM energy density J/m^3 (Eta Carinae)
        R_t             = 1.0;        // Sum cos resonance factor (neutral)

        // Catalogue constants
        g_H             = 1.252e46;   // UQFF cosmic orbital g-factor (PAPER_270)
        mu_B            = 9.274e-24;  // Bohr magneton (J/T)
        B0              = 1.0e-4;     // Magnetic field (T)
        h_planck        = 1.0546e-34; // Reduced Planck hbar (J·s)
        omega_0_base    = 1.0e-12;    // Base UQFF frequency (rad/s)

        // Physics constants for compute_g_UQFF()
        Lambda          = 1.114e-52;  // Cosmological constant (m^-2)
        c_light         = 2.998e8;    // Speed of light (m/s)
        hbar            = 1.0546e-34; // = h_planck
        delta_x         = 1.0e-15;   // Position uncertainty (m, nuclear scale)
        delta_p         = hbar / delta_x; // Momentum uncertainty from HUP (kg·m/s)
        integral_psi    = 1.0;        // Normalized wavefunction overlap
        t_Hubble        = 4.352e17;   // Hubble time (s = 13.8e9 yr × 3.15576e7 s/yr)

        // F_U_Bi_i time-evolution
        tau_SF          = 1.0e6 * 3.15576e7; // 1 Myr in seconds (canonical)
        f_TRZ           = 0.1;        // TRZ phase coupling factor (canonical)

        // UQFF 2.0 — 3-tier buoyancy (CP3/PAPER_198 canonical)
        beta_i          = 0.61;       // Buoyancy coupling constant
        omega_g         = 7.3e-16;    // Gravitational buoyancy frequency (rad/s)
        U_UA            = 1.0e-11;    // [UA] unit analytic factor
        M_EtaCar        = 120.0 * 1.989e30; // 120 M_sun = 2.387e32 kg
        r_EtaCar        = 7.11e19;    // 7,500 ly = 7.11e19 m
        M_GC            = 7.956e36;   // Sgr A* ~4e6 M_sun (outer-frame)
        r_GC            = 2.623e20;   // ~8.5 kpc

        // UQFF 2.0 self-expanding framework
        logging_enabled = false;
        dynamic_params.clear();

        updateCache();
    }

    // Cache update (e.g., DPM_resonance)
    void updateCache() {
        // Long-form DPM_resonance calculation (from doc example for Eta Carinae)
        double step1_g_H = 1.252e46;
        double step2_muB_B0 = mu_B * B0;  // 9.274e-24 * 1e-4 = 9.274e-28
        double step3_g_muB_B0 = step1_g_H * step2_muB_B0;  // 1.16e19
        double step4_h_omega0 = h_planck * omega_0_base;  // 1.0546e-46
        DPM_resonance = step3_g_muB_B0 / step4_h_omega0;  // 1.1e65 base
        // Scaled for Q_wave ≈ 3.11e9 J/m³
        DPM_resonance *= 2.82e-56;  // Adjustment factor from doc
    }

    // Universal setter for catalogue variables
    bool setVariable(const std::string& varName, double newValue) {
        if      (varName == "F_U_Bi_i")       { F_U_Bi_i       = newValue; }
        else if (varName == "g_H")             { g_H             = newValue; }
        else if (varName == "neutron_factor")  { neutron_factor  = newValue; }
        else if (varName == "water_state")     { water_state     = newValue; }
        else if (varName == "H_abundance")     { H_abundance     = newValue; }
        else if (varName == "omega_thz")       { omega_thz       = newValue; }
        else if (varName == "omega_0")         { omega_0         = newValue; }
        else if (varName == "B0")              { B0              = newValue; }
        else if (varName == "mu_B")            { mu_B            = newValue; }
        else if (varName == "h_planck")        { h_planck        = newValue; hbar = newValue; }
        else if (varName == "omega_0_base")    { omega_0_base    = newValue; }
        else if (varName == "LENR_term")       { LENR_term       = newValue; }
        else if (varName == "activation_term") { activation_term = newValue; }
        else if (varName == "DE_term")         { DE_term         = newValue; }
        else if (varName == "resonance_term")  { resonance_term  = newValue; }
        else if (varName == "rel_term")        { rel_term        = newValue; }
        else if (varName == "f_TRZ")           { f_TRZ           = newValue; }
        else if (varName == "tau_SF")          { tau_SF          = newValue; }
        else if (varName == "Lambda")          { Lambda          = newValue; }
        else if (varName == "t_Hubble")        { t_Hubble        = newValue; }
        else if (varName == "delta_x")         { delta_x         = newValue; }
        else if (varName == "delta_p")         { delta_p         = newValue; }
        else if (varName == "beta_i")          { beta_i          = newValue; }
        else if (varName == "omega_g")         { omega_g         = newValue; }
        else if (varName == "U_UA")            { U_UA            = newValue; }
        else if (varName == "M_EtaCar")        { M_EtaCar        = newValue; }
        else if (varName == "r_EtaCar")        { r_EtaCar        = newValue; }
        else if (varName == "M_GC")            { M_GC            = newValue; }
        else if (varName == "r_GC")            { r_GC            = newValue; }
        else if (varName == "k_vac")           { k_vac           = newValue; }
        else if (varName == "conduit_scale")   { conduit_scale   = newValue; }
        else if (varName == "string_wave")     { string_wave     = newValue; }
        else if (varName == "E_DPM")           { E_DPM           = newValue; }
        else if (varName == "k_thz")           { k_thz           = newValue; }
        else if (varName == "k_conduit")       { k_conduit       = newValue; }
        else {
            std::cerr << "[SOURCE10] Error: Unknown variable '" << varName << "'.\n";
            return false;
        }
        updateCache();
        return true;
    }

    // Compute F_U_Bi_i (UQFF Core Buoyancy — Catalogue master 4-term equation)
    double compute_F_U_Bi_i(double t) {
        auto _t0 = std::chrono::high_resolution_clock::now();
        // Long-form 4-term: integrand×x2 + LENR(t) + DE+resonance + rel×TRZ
        double term1 = integrand * x_2;
        double term2 = LENR_term * activation_term * std::exp(-t / tau_SF);
        double term3 = DE_term + resonance_term * neutron_factor;
        double term4 = rel_term * (1.0 + f_TRZ);
        double result = term1 + term2 + term3 + term4;
        if (logging_enabled) {
            auto _t1 = std::chrono::high_resolution_clock::now();
            double _ms = std::chrono::duration<double, std::milli>(_t1 - _t0).count();
            std::ostringstream oss;
            oss << std::scientific << "compute_F_U_Bi_i: t=" << t
                << " term1=" << term1 << " term2=" << term2
                << " term3=" << term3 << " term4=" << term4
                << " result=" << result << " elapsed=" << _ms << "ms";
            log(oss.str());
        }
        return result;
    }

    // Batch compute F_U_Bi_i over a vector of time points — returns one result per t
    std::vector<double> batch_compute_F_U_Bi_i(const std::vector<double>& times) {
        std::vector<double> results;
        results.reserve(times.size());
        auto _t0 = std::chrono::high_resolution_clock::now();
        for (double t : times) {
            double term1 = integrand * x_2;
            double term2 = LENR_term * activation_term * std::exp(-t / tau_SF);
            double term3 = DE_term + resonance_term * neutron_factor;
            double term4 = rel_term * (1.0 + f_TRZ);
            results.push_back(term1 + term2 + term3 + term4);
        }
        if (logging_enabled) {
            auto _t1 = std::chrono::high_resolution_clock::now();
            double _ms = std::chrono::duration<double, std::milli>(_t1 - _t0).count();
            std::ostringstream oss;
            oss << "batch_compute_F_U_Bi_i: " << times.size() << " steps, elapsed=" << _ms << "ms";
            log(oss.str());
        }
        return results;
    }

    // Compute g_UQFF(r, t) — 26-layer Triadic sum + Lambda + quantum + 3-tier CP3/PAPER_198 buoyancy
    double compute_g_UQFF(double r_input, double t) {
        auto _t0 = std::chrono::high_resolution_clock::now();
        // 26-layer Triadic MUGE: g = Σ(Ug1_i + Ug2_i + Ug3_i + Ug4_i, i=1..26)
        double sum_Ug = 0.0;
        for (int i = 0; i < 26; ++i) {
            sum_Ug += Ug1_vec[i] + Ug2_vec[i] + Ug3_vec[i] + Ug4_vec[i];
        }
        // Cosmological and quantum terms
        double Lambda_term  = (Lambda * c_light * c_light) / 3.0;
        double quantum_term = (hbar / std::sqrt(delta_x * delta_p)) * integral_psi
                              * (2.0 * M_PI / t_Hubble);
        // 3-tier CP3/PAPER_198 buoyancy (Eta Carinae Tier 2 + Sgr A* Tier 3 outer-frame)
        double ug1_base   = (r_input > 0.0) ? (6.674e-11 * M_EtaCar / (r_input * r_input)) : 0.0;
        double term_Ubi   = 0.5 * ug1_base;                                          // Tier 1 static
        double term_FUBii = -beta_i * ug1_base * omega_g * (M_EtaCar / r_input)
                            * U_UA * std::cos(M_PI * t);                             // Tier 2 compact
        double term_Ub_i  = -beta_i * ug1_base * omega_g * (M_GC / r_GC)
                            * U_UA * std::cos(M_PI * t);                             // Tier 3 outer-frame
        double FU_diag    = -(sum_Ug / 26.0 + term_Ubi) * ug1_base;                // diagnostic
        double g_total    = sum_Ug + Lambda_term + quantum_term
                            + term_Ubi + term_FUBii + term_Ub_i;
        if (logging_enabled) {
            auto _t1 = std::chrono::high_resolution_clock::now();
            double _ms = std::chrono::duration<double, std::milli>(_t1 - _t0).count();
            std::ostringstream oss;
            oss << std::scientific
                << "compute_g_UQFF: r=" << r_input << " t=" << t
                << " sum_Ug=" << sum_Ug << " Lambda=" << Lambda_term
                << " quantum=" << quantum_term
                << " Ubi=" << term_Ubi << " FUBii=" << term_FUBii
                << " Ub_i=" << term_Ub_i << " FU_diag=" << FU_diag
                << " g_total=" << g_total << " elapsed=" << _ms << "ms";
            log(oss.str());
        }
        return g_total;
    }

    // Long-form Resonance Solution (DPM_resonance method)
    double compute_DPM_resonance() {
        // From doc long-form (Eta Carinae example)
        double g_H_val = g_H;
        double muB_B0 = mu_B * B0;
        double g_muB_B0 = g_H_val * muB_B0;
        double h_omega0 = h_planck * omega_0_base;
        double base = g_muB_B0 / h_omega0;
        double adjusted = base * 2.82e-56;  // Scaled to 3.11e9 J/m³
        return adjusted;
    }

    // Debug/Output: Print Catalogue Summary (UQFF 2.0)
    void printCatalogue(std::ostream& os = std::cout) const {
        os << std::scientific << std::setprecision(4);
        os << "=== UQFF Source10 Catalogue (UQFF 2.0) ===" << std::endl;
        os << "F_U_Bi_i       : " << F_U_Bi_i       << " N (Eta Carinae benchmark)" << std::endl;
        os << "g_H            : " << g_H             << " (UQFF cosmic orbital amplifier)" << std::endl;
        os << "DPM_resonance  : " << DPM_resonance   << " J/m^3" << std::endl;
        os << "E_DPM          : " << E_DPM           << " J/m^3" << std::endl;
        os << "neutron_factor : " << neutron_factor  << " (Kozima gate: 1=open)" << std::endl;
        os << "water_state    : " << water_state     << " (fluid gate: 1=open)" << std::endl;
        os << "H_abundance    : " << H_abundance     << " (cosmic H fraction)" << std::endl;
        os << "omega_thz      : " << omega_thz       << " rad/s (THz gate)" << std::endl;
        os << "omega_0        : " << omega_0         << " rad/s" << std::endl;
        os << "THz_ratio      : " << (omega_thz / omega_0) << " (Colman-Gillespie ~1.25)" << std::endl;
        os << "k_vac=G        : " << k_vac           << " m^3/kg/s^2 (gravitational vacuum drag)" << std::endl;
        os << "F_vac_rep      : " << F_vac_rep       << " N" << std::endl;
        os << "F_thz_shock    : " << F_thz_shock     << " N" << std::endl;
        os << "F_conduit      : " << F_conduit       << " N" << std::endl;
        os << "F_spooky       : " << F_spooky        << " N" << std::endl;
        os << "f_TRZ          : " << f_TRZ           << " (TRZ coupling)" << std::endl;
        os << "tau_SF         : " << tau_SF          << " s" << std::endl;
        os << "beta_i         : " << beta_i          << " (buoyancy coupling)" << std::endl;
        os << "omega_g        : " << omega_g         << " rad/s" << std::endl;
        os << "M_EtaCar       : " << M_EtaCar        << " kg (Tier 2 compact)" << std::endl;
        os << "r_EtaCar       : " << r_EtaCar        << " m" << std::endl;
        os << "M_GC (Sgr A*)  : " << M_GC           << " kg (Tier 3 outer-frame)" << std::endl;
        os << "r_GC           : " << r_GC            << " m" << std::endl;
        os << "t_Hubble       : " << t_Hubble        << " s" << std::endl;
        os << "Lambda         : " << Lambda          << " m^-2" << std::endl;
        os << "26-layer vecs  : " << Ug1_vec.size()  << " entries each" << std::endl;
        os << "logging        : " << (logging_enabled ? "ON" : "OFF") << std::endl;
        os << "dynamic_params : " << dynamic_params.size() << " entries" << std::endl;
    }

    // Example: Eta Carinae Buoyancy (from doc)
    double exampleEtaCarinae(double t) {
        return compute_F_U_Bi_i(t);  // ~2.11e208 N
    }

    // -------------------------------------------------------------------------
    // UQFF 2.0 Self-Expanding Framework Methods
    // -------------------------------------------------------------------------

    void setEnableLogging(bool enabled) { logging_enabled = enabled; }
    bool getLoggingEnabled() const { return logging_enabled; }

    void setDynamicParameter(const std::string& name, double value) {
        dynamic_params[name] = value;
        updateCache();
    }

    double getDynamicParameter(const std::string& name) const {
        auto it = dynamic_params.find(name);
        return (it != dynamic_params.end()) ? it->second : 0.0;
    }

    // Named scaling factor interface — wraps setDynamicParameter for LENR/DE/resonance scaling
    void setScalingFactor(const std::string& name, double value) {
        setDynamicParameter(name, value);
    }
    double getScalingFactor(const std::string& name) const {
        return getDynamicParameter(name);
    }

    void exportState(const std::string& filename = "UQFFSource10_state.txt") const {
        std::ofstream ofs(filename);
        if (!ofs.is_open()) return;
        ofs << std::scientific;
        ofs << "F_U_Bi_i "        << F_U_Bi_i        << "\n";
        ofs << "integrand "       << integrand       << "\n";
        ofs << "x_2 "             << x_2             << "\n";
        ofs << "LENR_term "       << LENR_term       << "\n";
        ofs << "activation_term " << activation_term << "\n";
        ofs << "DE_term "         << DE_term         << "\n";
        ofs << "resonance_term "  << resonance_term  << "\n";
        ofs << "rel_term "        << rel_term        << "\n";
        ofs << "f_TRZ "           << f_TRZ           << "\n";
        ofs << "tau_SF "          << tau_SF          << "\n";
        ofs << "F_vac_rep "       << F_vac_rep       << "\n";
        ofs << "k_vac "           << k_vac           << "\n";
        ofs << "delta_rho_vac "   << delta_rho_vac   << "\n";
        ofs << "M_vac "           << M_vac           << "\n";
        ofs << "v_vac "           << v_vac           << "\n";
        ofs << "F_thz_shock "     << F_thz_shock     << "\n";
        ofs << "omega_thz "       << omega_thz       << "\n";
        ofs << "omega_0 "         << omega_0         << "\n";
        ofs << "neutron_factor "  << neutron_factor  << "\n";
        ofs << "water_state "     << water_state     << "\n";
        ofs << "H_abundance "     << H_abundance     << "\n";
        ofs << "F_conduit "       << F_conduit       << "\n";
        ofs << "F_spooky "        << F_spooky        << "\n";
        ofs << "g_H "             << g_H             << "\n";
        ofs << "mu_B "            << mu_B            << "\n";
        ofs << "B0 "              << B0              << "\n";
        ofs << "h_planck "        << h_planck        << "\n";
        ofs << "omega_0_base "    << omega_0_base    << "\n";
        ofs << "DPM_resonance "   << DPM_resonance   << "\n";
        ofs << "E_DPM "           << E_DPM           << "\n";
        ofs << "Lambda "          << Lambda          << "\n";
        ofs << "c_light "         << c_light         << "\n";
        ofs << "t_Hubble "        << t_Hubble        << "\n";
        ofs << "delta_x "         << delta_x         << "\n";
        ofs << "delta_p "         << delta_p         << "\n";
        ofs << "beta_i "          << beta_i          << "\n";
        ofs << "omega_g "         << omega_g         << "\n";
        ofs << "U_UA "            << U_UA            << "\n";
        ofs << "M_EtaCar "        << M_EtaCar        << "\n";
        ofs << "r_EtaCar "        << r_EtaCar        << "\n";
        ofs << "M_GC "            << M_GC            << "\n";
        ofs << "r_GC "            << r_GC            << "\n";
        for (auto& kv : dynamic_params) {
            ofs << "dyn_" << kv.first << " " << kv.second << "\n";
        }
        ofs.close();
    }

    template <typename T>
    double cross_validate(T& other_module, double t_years = 50.0e6) {
        double t_s     = t_years * 3.15576e7;
        double g_this  = compute_g_UQFF(r_EtaCar, t_s);
        double g_other = other_module.compute_g_UQFF(r_EtaCar, t_s);
        double denom   = std::abs(g_this) + std::abs(g_other) + 1.0e-300;
        return std::abs(g_this - g_other) / denom;
    }

private:
    void log(const std::string& msg) const {
        if (logging_enabled) { std::cout << "[LOG][SOURCE10] " << msg << "\n"; }
    }
};

// Integration Note for Source12.cpp: Include this header as #include "UQFFSource10.h"
// Use: UQFFSource10 source10; source10.setVariable("g_H", 1.252e46);
// This installs Source10 as first primary text module, aggregating all prior for Star Magic architecture.

#endif // UQFF_SOURCE10_H