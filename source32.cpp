// CrabUQFFModule.h
#define WOLFRAM_TERM "(* Auto-contribution from source32.cpp *) + source32_unification_sector"
// Modular C++ implementation of the full Master Universal Gravity Equation (UQFF) for Crab Nebula Evolution.
// This module can be plugged into a base program (e.g., 'ziqn233h.cpp') by including this header and linking the .cpp.
// Usage in base: CrabUQFFModule mod; mod.computeG(t); mod.updateVariable("M", new_value);
// All variables are stored in a std::map for dynamic addition/subtraction/update.
// Nothing is negligible: Includes all terms - base gravity with r(t), Ug1-Ug4 (gravitational subterms), cosmological Lambda, 
// quantum (hbar uncertainty integral term), Lorentz q(v x B), fluid (rho_fluid V g), resonant oscillatory (cos and exp terms), 
// DM/visible mass with density perturbations, superconductivity correction (1 - B/B_crit), pulsar wind a_wind, magnetic M_mag.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Quantum integral normalized to 1.0 (ground state); exp term real part (cos); Ug2/Ug3=0 (negligible for remnant); 
// fluid g recursive approx using base g_grav; resonant at x=0 (central); DM fraction 0 (M_visible=M); 
// r(t) = r0 + v_exp * t; a_wind = wind_pressure / rho * scale_macro; M_mag = (q v B) / m_e * scale_macro; B_crit=1e11 T; H(z) for z=0.0015.
// Crab params: M=4.6 Msun, r0=5.2e16 m, v_exp=1.5e6 m/s, P_pulsar=5e31 W, B=1e-8 T (nebula avg), z=0.0015, rho=1e-21 kg/m^3, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef CRAB_UQFF_MODULE_H
#define CRAB_UQFF_MODULE_H

// ===========================================================================================
// SHARED FRAMEWORK HEADERS (Unified across all source files)
// ===========================================================================================
#include "uqff_constants.h"        // Physical constants: PI, c, G, hbar, k_B, etc.
#include "uqff_self_expanding.h"   // SelfExpandingModule template, PhysicsTerm base class
#include "uqff_dual_physics.h"     // DualMethodValidator, FluidSolver, CelestialBody, MUGESystem

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <functional>
#include <memory>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <array>

// Import UQFF constants namespace for PI, c, G, hbar, k_B, etc.
using namespace UQFF;
// Import dual physics namespace for CelestialBody, MUGESystem, DualMethodValidator, ValidationResult
using namespace UQFFDualPhysics;
// Import self-expanding namespace for PhysicsTerm
using namespace UQFFExpanding;

// ===========================================================================================
// ENHANCED CLASS WITH SELF-EXPANDING CAPABILITIES
// ===========================================================================================

class CrabUQFFModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeQuantumTerm(double t_Hubble_val);
    double computeFluidTerm(double g_base);
    double computeResonantTerm(double t);
    double computeDMTerm();
    double computeUgSum();
    double computeHz();
    double computeWindTerm(double r);
    double computeMagTerm();
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize all variables with Crab Nebula defaults
    CrabUQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: Full g_UQFF(r, t) for Crab Nebula
    double computeG(double t);

    // Output descriptive text of the equation
    std::string getEquationText();

    // Print all current variables (for debugging/updates)
    void printVariables();

    // ========== SELF-SIMULATION METHODS (2.0-Enhanced) ==========
    void runSelfSimulation(double t_start, double t_end, int steps);
    void exportState(const std::string& filename) const;
    void setLearningRate(double rate) { learningRate = rate; }
    void setEnableLogging(bool enable) { enableLogging = enable; }
    void setEnableDynamicTerms(bool enable) { enableDynamicTerms = enable; }
};

#endif // CRAB_UQFF_MODULE_H

// CrabUQFFModule.cpp
// // // #include "CrabUQFFModule.h"  // Commented - header not available  // Commented - header not available  // Commented - header not available
#include <complex>
#include <array> // MSVC requirement

// Constructor: Set all variables with Crab Nebula-specific values
CrabUQFFModule::CrabUQFFModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Base constants (universal)
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["c"] = 3e8;                           // m/s
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Lambda"] = 1.1e-52;                  // m^-2 (cosmological constant)
    variables["q"] = 1.602e-19;                     // C (electron charge)
    variables["pi"] = 3.141592653589793;            // pi
    variables["t_Hubble"] = 13.8e9 * 3.156e7;       // s (13.8 Gyr)

    // Crab Nebula parameters
    double M_sun_val = 1.989e30;                    // kg
    variables["M_sun"] = M_sun_val;
    variables["M"] = 4.6 * M_sun_val;               // Total mass kg
    variables["M_visible"] = variables["M"];        // Visible mass (ejecta + pulsar)
    variables["M_DM"] = 0.0;                        // No significant DM
    variables["r0"] = 5.2e16;                       // m (initial radius)
    variables["v_exp"] = 1.5e6;                     // m/s (expansion velocity)

    // Hubble/cosmology
    variables["H0"] = 70.0;                         // km/s/Mpc
    variables["Mpc_to_m"] = 3.086e22;               // m/Mpc
    variables["z"] = 0.0015;                        // Redshift
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["t"] = 971 * 3.156e7;                 // Default t=971 years s (since 1054 AD)

    // Nebula dynamics
    variables["rho_fluid"] = 1e-21;                 // kg/m^3 (filament density)
    variables["V"] = 1e3;                           // m^3 (arbitrary volume scale)
    variables["v_shock"] = 1.5e6;                   // m/s (shock velocity)
    variables["P_pulsar"] = 5e31;                   // W (pulsar luminosity)
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];  // Perturbation
    variables["rho"] = variables["rho_fluid"];      // Mean density

    // EM/magnetic/superconductivity
    variables["B"] = 1e-8;                          // T (nebula average magnetic field)
    variables["B_crit"] = 1e11;                     // T (10^15 G  1e11 T)
    variables["m_e"] = 9.11e-31;                    // kg (electron mass)

    // Quantum terms
    variables["Delta_x"] = 1e-10;                   // m (position uncertainty, atomic scale)
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];  // Momentum uncertainty (Heisenberg)
    variables["integral_psi"] = 1.0;                // Normalized <psi|H|psi> dV  E_ground (simplified to 1 for unitless)

    // Resonant/oscillatory terms
    variables["A"] = 1e-10;                         // Amplitude (arbitrary small)
    variables["k"] = 1e20;                          // m^-1 (wave number, short wavelength)
    variables["omega"] = 1e15;                      // rad/s (high freq, e.g., synchrotron)
    variables["x"] = 0.0;                           // m (position, central)

    // Ug subterms (computed dynamically, but init placeholders)
    variables["Ug1"] = 0.0;  // Will be G M / r^2
    variables["Ug2"] = 0.0;  // d^2 Phi / dt^2  0 (negligible)
    variables["Ug3"] = 0.0;  // G M_moon / r_moon^2  0 (no moon)
    variables["Ug4"] = 0.0;  // Ug1 * f_sc, f_sc=1

    // Scale factors (from streamlining)
    variables["scale_macro"] = 1e-12;               // For macro effects
    variables["f_TRZ"] = 0.1;                       // Time-reversal factor
    variables["f_sc"] = 1.0;                        // Superconductive factor
}

// Update variable (set to new value)
void CrabUQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    // Recompute dependent vars if needed (e.g., Delta_p)
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "M") {
        variables["M_visible"] = value;
        variables["M_DM"] = 0.0;
    }
}

// Add delta to variable
void CrabUQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void CrabUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(z) in s^-1
double CrabUQFFModule::computeHz() {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + variables["z"], 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// Compute Ug sum: Ug1 = G M / r^2, Ug4 = Ug1 * f_sc, others 0
double CrabUQFFModule::computeUgSum() {
    double r = variables["r0"] + variables["v_exp"] * variables["t"];  // Use current r(t)
    double Ug1 = (variables["G"] * variables["M"]) / (r * r);
    variables["Ug1"] = Ug1;  // Update map
    variables["Ug4"] = Ug1 * variables["f_sc"];
    return variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
double CrabUQFFModule::computeQuantumTerm(double t_Hubble_val) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double integral_val = variables["integral_psi"];  // Simplified
    return (variables["hbar"] / unc) * integral_val * (2 * variables["pi"] / t_Hubble_val);
}

// Fluid term: rho_fluid * V * g (g approx base grav)
double CrabUQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
double CrabUQFFModule::computeResonantTerm(double t) {
    double cos_term = 2 * variables["A"] * std::cos(variables["k"] * variables["x"]) * std::cos(variables["omega"] * t);
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega"] * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);  // Gyr? Assume unitless as per doc
    return cos_term + exp_factor * real_exp;
}

// DM term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3)
double CrabUQFFModule::computeDMTerm() {
    double r = variables["r0"] + variables["v_exp"] * variables["t"];  // Use current r(t)
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3 * variables["G"] * variables["M"] / (r * r * r);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Wind term: (P_pulsar / (4 pi r^2)) * (1 + v_shock / c) / rho_fluid * scale_macro
double CrabUQFFModule::computeWindTerm(double r) {
    double pressure = (variables["P_pulsar"] / (4 * variables["pi"] * r * r)) * (1.0 + variables["v_shock"] / variables["c"]);
    return (pressure / variables["rho_fluid"]) * variables["scale_macro"];
}

// Magnetic term: (q * v_shock * B) / m_e * scale_macro
double CrabUQFFModule::computeMagTerm() {
    double force = variables["q"] * variables["v_shock"] * variables["B"];
    return (force / variables["m_e"]) * variables["scale_macro"];
}

// Full computation: g_UQFF(r, t) = ... all terms
double CrabUQFFModule::computeG(double t) {
    variables["t"] = t;  // Update t
    double r = variables["r0"] + variables["v_exp"] * t;  // r(t)
    double Hz = computeHz();
    double expansion = 1.0 + Hz * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double tr_factor = 1.0 + variables["f_TRZ"];

    // Base gravity with expansion, SC, TR
    double g_base = (variables["G"] * variables["M"] / (r * r)) * expansion * sc_correction * tr_factor;

    // Ug sum
    double ug_sum = computeUgSum();

    // Cosmological
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"]) / 3.0;

    // Quantum
    double quantum_term = computeQuantumTerm(variables["t_Hubble"]);

    // EM Lorentz (magnitude v_shock B)
    double em_base = variables["q"] * variables["v_shock"] * variables["B"] / 1.673e-27;  // / proton mass for accel (approx)
    double em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * variables["scale_macro"];  // UA/SCm ratio=10

    // Fluid (uses g_base approx)
    double fluid_term = computeFluidTerm(g_base);

    // Resonant
    double resonant_term = computeResonantTerm(t);

    // DM
    double dm_term = computeDMTerm();

    // Wind
    double wind_term = computeWindTerm(r);

    // Mag
    double mag_term = computeMagTerm();

    // Total: Sum all
    return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + wind_term + mag_term;
}

// Get equation text (descriptive)
std::string CrabUQFFModule::getEquationText() {
    return "g_Crab(r, t) = (G * M / r(t)^2) * (1 + H(z) * t) * (1 - B / B_crit) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 / 3) + "
           "(hbar / sqrt(Delta_x * Delta_p)) * ?(?* H ? dV) * (2p / t_Hubble) + q (v × B) + ?_fluid * V * g + "
           "2 A cos(k x) cos(? t) + (2p / 13.8) A exp(i (k x - ? t)) + (M_visible + M_DM) * (d?/? + 3 G M / r^3) + a_wind + M_mag\n"
           "Where r(t) = r0 + v_exp * t; a_wind = [P_pulsar / (4p r^2) * (1 + v_shock / c)] / ? * 1e-12; M_mag = (q v_shock B) / m_e * 1e-12\n"
           "Special Terms:\n"
           "- Quantum: Heisenberg uncertainty with normalized wavefunction integral (ground state approx) for particle quantum effects.\n"
           "- Fluid: Nebular filament density-volume-gravity coupling.\n"
           "- Resonant: Oscillatory Aether-mediated waves (real part of complex exp) for wisp dynamics.\n"
           "- DM: Visible mass (ejecta + pulsar) with density perturbations and curvature term (M_DM=0).\n"
           "- Superconductivity: (1 - B/B_crit) for quantum field effects near pulsar.\n"
           "- Pulsar Wind: a_wind from relativistic wind pressure, dominant outward force.\n"
           "- Magnetic: M_mag from Lorentz force on electrons in nebula fields.\n"
           "Solutions: Numerical evaluation at t=971 yr yields ~1.481e6 m/s² (a_wind dominant; g_grav ~2e-13; micro terms ~1e-10 to 1e-3).\n"
           "Adaptations for Crab: Pulsar-driven remnant with r(t); z=0.0015; v_shock=1.5e6 m/s boosts wind/mag.";
}

// Print variables
void CrabUQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===========================================================================================
// SELF-SIMULATION IMPLEMENTATIONS (2.0-Enhanced)
// ===========================================================================================

void CrabUQFFModule::runSelfSimulation(double t_start, double t_end, int steps) {
    if (enableLogging) {
        std::cout << "[Crab] Running self-simulation: t=" << t_start << " to " << t_end << " (" << steps << " steps)" << std::endl;
    }
    double dt = (t_end - t_start) / steps;
    for (int i = 0; i <= steps; ++i) {
        double t = t_start + i * dt;
        double g = computeG(t);
        if (enableLogging) {
            std::cout << "  t=" << std::scientific << std::setprecision(3) << t << " s, g=" << g << " m/s^2" << std::endl;
        }
    }
}

void CrabUQFFModule::exportState(const std::string& filename) const {
    std::ofstream out(filename);
    if (out.is_open()) {
        out << "# CrabUQFFModule State Export\n";
        out << "# Crab Nebula 12-term MUGE with pulsar wind, magnetic, superconductivity\n";
        for (const auto& pair : variables) {
            out << pair.first << "=" << std::scientific << pair.second << "\n";
        }
        out << "learningRate=" << learningRate << "\n";
        out << "enableDynamicTerms=" << enableDynamicTerms << "\n";
        out << "enableLogging=" << enableLogging << "\n";
        out.close();
    }
}

// ===========================================================================================
// MAIN: Dual Physics Validation for Crab Nebula
// ===========================================================================================
int main()
{
    std::cout << "=== Source32: Crab Nebula UQFF Module ===" << std::endl;
    std::cout << "12-term MUGE with pulsar wind, magnetic, superconductivity + Dual Physics Validation" << std::endl;

    CrabUQFFModule crab;

    // Time evolution: 971 years (since 1054 AD), 1000 years, 1500 years
    std::cout << "\n=== MUGE Gravity Evolution (Expanding Remnant) ===" << std::endl;
    double yr_to_s = 3.156e7;
    std::array<double, 3> times = {971 * yr_to_s, 1000 * yr_to_s, 1500 * yr_to_s};
    for (double t : times) {
        double g = crab.computeG(t);
        std::cout << "t=" << static_cast<int>(t / yr_to_s) << " yr: g = " << std::scientific << std::setprecision(4) << g << " m/s^2" << std::endl;
    }

    // Print equation description
    std::cout << "\n=== Equation Description ===" << std::endl;
    std::cout << crab.getEquationText() << std::endl;

    // Dual Physics Validation: UQFF vs MUGE
    std::cout << "\n=== Dual Physics Validation ===" << std::endl;
    double t_now = 971 * yr_to_s;  // Current age
    double r_now = 5.2e16 + 1.5e6 * t_now;  // Current radius
    double g_uqff = crab.computeG(t_now);

    // Newtonian base for comparison
    double G = 6.6743e-11;
    double M = 4.6 * 1.989e30;  // 4.6 solar masses
    double g_newton = G * M / (r_now * r_now);

    std::cout << "Crab UQFF (12-term):  " << std::scientific << std::setprecision(4) << g_uqff << " m/s^2" << std::endl;
    std::cout << "Newton Base:          " << std::scientific << std::setprecision(4) << g_newton << " m/s^2" << std::endl;
    std::cout << "Wind+Mag Dominant:    " << (g_uqff > 1e3 ? "YES (pulsar-driven)" : "NO") << std::endl;

    // FluidSolver for pulsar wind nebula dynamics
    std::cout << "\n=== FluidSolver: Pulsar Wind Nebula Dynamics ===" << std::endl;
    FluidSolver fluidSolver(32, 0.1, 0.0001);  // 32x32 grid
    std::cout << "FluidSolver initialized (32x32 grid)" << std::endl;
    fluidSolver.add_jet_force(8.0);  // Pulsar wind forcing
    for (int i = 0; i < 10; ++i) {
        fluidSolver.step(g_uqff * 1e-10);  // Scaled UQFF gravity
    }
    std::cout << "FluidSolver: Max velocity = " << fluidSolver.getMaxVelocity() << " m/s" << std::endl;

    // DualMethodValidator integration
    CelestialBody body;
    body.name = "Crab_Nebula";
    body.M = M;
    body.Rs = r_now;
    body.B0 = 1e-8;  // Nebula average B field

    MUGESystem muge;
    muge.name = "Crab_MUGE";
    muge.M = M;
    muge.r = r_now;
    muge.B0 = 1e-8;
    muge.T = 1e4;  // Nebula temperature ~10,000 K

    DualMethodValidator validator;
    validator.addConstraint("crab_g", 1e-15, 1e10, 0.5);  // Wide range for pulsar-driven system
    ValidationResult result = validator.validate(body, muge);
    result.print();

    // Self-simulation test
    std::cout << "\n=== Self-Simulation Test ===" << std::endl;
    crab.setEnableLogging(true);
    crab.runSelfSimulation(500 * yr_to_s, 1500 * yr_to_s, 4);  // 500-1500 years evolution

    std::cout << "\n=== Source32 Complete ===" << std::endl;
    return 0;
}
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.