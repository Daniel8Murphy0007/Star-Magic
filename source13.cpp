// source13.cpp - UQFF Magnetar SGR 1745-2900 MUGE Calculator
// Version: 2.0-Enhanced (Module-Aligned)
// Framework: Self-Expanding UQFF with complete MUGE (12 physics terms)
// Refactored: Added UQFFModule13 wrapper, standard interface, main()
// Enhanced: January 22, 2026 - Added FluidSolver, DualMethodValidator, DualPhysicsMethods

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <iomanip>
#include <random>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "uqff_constants.h"
#include "uqff_self_expanding.h"
#include "uqff_dual_physics.h"

using namespace UQFF;
// Note: UQFFDualPhysics and UQFFExpanding used locally to avoid
// conflicts with source13's own struct definitions

// ============================================================================
// UQFFConfig13 - Singleton Configuration for Magnetar Physics
// ============================================================================
class UQFFConfig13 {
public:
    static UQFFConfig13& getInstance() {
        static UQFFConfig13 instance;
        return instance;
    }
    
    // Magnetar SGR 1745-2900 parameters
    double M_magnetar = 1.4 * SUN_MASS_KG;  // 1.4 solar masses
    double r_magnetar = 1e4;                 // 10 km radius
    double B0 = 2e10;                        // Initial magnetic field (T)
    double B_crit = 1e11;                    // Critical B field (T)
    double P_init = 3.76;                    // Initial rotation period (s)
    double tau_Omega = 10000 * 3.15576e7;    // Omega decay timescale (s)
    double tau_decay = 3.5 * 365.25 * 24 * 3600;  // 3.5 years decay
    double L0_W = 5e28;                      // Initial luminosity (W)
    
    // Black hole (Sgr A*) parameters
    double M_BH = 4e6 * SUN_MASS_KG;         // 4 million solar masses
    double r_BH = 2.83e16;                   // Distance to Sgr A* (m)
    
    // Cosmological parameters
    double Hz = 2.269e-18;                   // Hubble parameter at z (s^-1)
    double Lambda = 1.1e-52;                 // Cosmological constant
    double t_Hubble = 13.8e9 * 3.15576e7;    // Hubble time (s)
    
    // Quantum/EM parameters
    double delta_x = 1e-10;                  // Position uncertainty (m)
    double q_charge = 1.602e-19;             // Proton charge (C)
    double proton_mass = 1.673e-27;          // Proton mass (kg)
    double v_surf = 1e6;                     // Surface velocity (m/s)
    double scale_EM = 1e-12;                 // EM scaling factor
    
    // Fluid/DM parameters
    double rho_fluid = 1e17;                 // Fluid density (kg/m³)
    double M_DM_factor = 0.1;                // Dark matter mass fraction
    double delta_rho_over_rho = 1e-5;        // Density perturbation
    
    // Oscillatory parameters
    double A_osc = 1e10;                     // Oscillatory amplitude (m/s²)
    
private:
    UQFFConfig13() = default;
    UQFFConfig13(const UQFFConfig13&) = delete;
    UQFFConfig13& operator=(const UQFFConfig13&) = delete;
};

// ============================================================================
// MagnetarSGR1745_2900 - Core Physics Class (preserved from original)
// ============================================================================
class MagnetarSGR1745_2900 {
private:
    double M, r, B, B_crit, f_sc;
    double Hz, Lambda, P_init, tau_Omega;
    double M_BH, r_BH, L0_W, tau_decay;
    double delta_x, delta_p, integral_psi;
    double rho_fluid, A_osc, k_osc, omega_osc, x_pos;
    double M_DM_factor, delta_rho_over_rho;
    double q_charge, v_surf, scale_EM, proton_mass;
    double t_Hubble, t_Hubble_gyr;
    double ug1_base;  // Cached G*M/r²

public:
    MagnetarSGR1745_2900() {
        auto& cfg = UQFFConfig13::getInstance();
        M = cfg.M_magnetar;
        r = cfg.r_magnetar;
        B = cfg.B0;
        B_crit = cfg.B_crit;
        f_sc = 1.0 - (B / B_crit);
        Hz = cfg.Hz;
        Lambda = cfg.Lambda;
        P_init = cfg.P_init;
        tau_Omega = cfg.tau_Omega;
        M_BH = cfg.M_BH;
        r_BH = cfg.r_BH;
        L0_W = cfg.L0_W;
        tau_decay = cfg.tau_decay;
        delta_x = cfg.delta_x;
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
        rho_fluid = cfg.rho_fluid;
        A_osc = cfg.A_osc;
        k_osc = 1.0 / r;
        omega_osc = 2.0 * PI / P_init;
        x_pos = r;
        M_DM_factor = cfg.M_DM_factor;
        delta_rho_over_rho = cfg.delta_rho_over_rho;
        q_charge = cfg.q_charge;
        v_surf = cfg.v_surf;
        scale_EM = cfg.scale_EM;
        proton_mass = cfg.proton_mass;
        t_Hubble = cfg.t_Hubble;
        t_Hubble_gyr = 13.8;
        
        updateCache();
    }
    
    void updateCache() {
        ug1_base = (G * M) / (r * r);
        f_sc = 1.0 - (B / B_crit);
    }
    
    // Getters
    double getMass() const { return M; }
    double getRadius() const { return r; }
    double getMagneticField() const { return B; }
    double getUg1Base() const { return ug1_base; }
    
    // Setters
    void setMass(double mass) { M = mass; updateCache(); }
    void setRadius(double radius) { r = radius; k_osc = 1.0/r; x_pos = r; updateCache(); }
    void setMagneticField(double field) { B = field; updateCache(); }
    
    // Omega(t) - spin frequency decay
    double Omega_t(double t) const {
        return (2.0 * PI / P_init) * std::exp(-t / tau_Omega);
    }
    
    // dOmega/dt
    double dOmega_dt(double t) const {
        double omega0 = 2.0 * PI / P_init;
        return omega0 * (-1.0 / tau_Omega) * std::exp(-t / tau_Omega);
    }
    
    // Ug terms (UQFF gravity components)
    double compute_Ug() const {
        double Ug1 = ug1_base;
        double Ug2 = 0.0;  // Charge term (not active for this magnetar)
        double Ug3 = 0.0;  // Rotation term (absorbed in Omega)
        double Ug4 = ug1_base * f_sc;  // Superconductive correction
        return Ug1 + Ug2 + Ug3 + Ug4;
    }
    
    // Volume
    double compute_V() const {
        return (4.0 / 3.0) * PI * r * r * r;
    }
    
    // Magnetic energy (J)
    double compute_M_mag() const {
        double V = compute_V();
        return (B * B / (2.0 * mu_0)) * V;
    }
    
    // Cumulative decay energy (J)
    double compute_cumulative_D(double t) const {
        return L0_W * tau_decay * (1.0 - std::exp(-t / tau_decay));
    }
    
    // Main MUGE computation - 12 physics terms
    double compute_g_Magnetar(double t) const {
        if (t < 0) return 0.0;
        
        double dOdt = dOmega_dt(t);
        double current_f_sc = 1.0 - (B / B_crit);
        
        // Term 1: Base gravity + H(z) + B corrections
        double corr_H = 1.0 + Hz * t;
        double corr_B = current_f_sc;
        double term1 = ug1_base * corr_H * corr_B;
        
        // Term 2: Black hole influence (Sgr A*)
        double term_BH = (G * M_BH) / (r_BH * r_BH);
        
        // Term 3: UQFF Ug sum
        double term2 = compute_Ug();
        
        // Term 4: Cosmological constant Λ
        double term3 = (Lambda * c * c) / 3.0;
        
        // Term 5: EM acceleration (v × B)
        double cross_vB = v_surf * B;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double term4 = em_base * scale_EM;
        
        // Term 6: Gravitational wave term
        double gw_prefactor = (G * M * M) / (std::pow(c, 4) * r);
        double term5 = gw_prefactor * (dOdt * dOdt);
        
        // Term 7: Quantum uncertainty
        double sqrt_unc = std::sqrt(delta_x * delta_p);
        double term_q = (hbar / sqrt_unc) * integral_psi * (2.0 * PI / t_Hubble);
        
        // Term 8: Fluid dynamics
        double V = compute_V();
        double term_fluid = (rho_fluid * V * ug1_base) / M;
        
        // Term 9: Oscillatory waves
        double term_osc1 = 2.0 * A_osc * std::cos(k_osc * x_pos) * std::cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2.0 * PI / t_Hubble_gyr) * A_osc * std::cos(arg);
        double term_osc = term_osc1 + term_osc2;
        
        // Term 10: Dark matter + density perturbation
        double M_dm = M * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3.0 * G * M / (r * r * r);
        double term_dm_force = (M + M_dm) * (pert1 + pert2);
        double term_DM = term_dm_force / M;
        
        // Term 11: Magnetic energy effective g
        double M_mag = compute_M_mag();
        double term_mag = M_mag / (M * r);
        
        // Term 12: Cumulative decay energy effective g
        double cum_D = compute_cumulative_D(t);
        double term_decay = cum_D / (M * r);
        
        // Total MUGE
        return term1 + term_BH + term2 + term3 + term4 + term5 + 
               term_q + term_fluid + term_osc + term_DM + term_mag + term_decay;
    }
    
    // Individual term breakdown for analysis
    void computeTermBreakdown(double t, double& term1, double& term_BH, double& term_Ug,
                              double& term_Lambda, double& term_EM, double& term_GW,
                              double& term_quantum, double& term_fluid, double& term_osc,
                              double& term_DM, double& term_mag, double& term_decay) const {
        double dOdt = dOmega_dt(t);
        double current_f_sc = 1.0 - (B / B_crit);
        
        term1 = ug1_base * (1.0 + Hz * t) * current_f_sc;
        term_BH = (G * M_BH) / (r_BH * r_BH);
        term_Ug = compute_Ug();
        term_Lambda = (Lambda * c * c) / 3.0;
        term_EM = (q_charge * v_surf * B / proton_mass) * scale_EM;
        term_GW = (G * M * M) / (std::pow(c, 4) * r) * (dOdt * dOdt);
        
        double sqrt_unc = std::sqrt(delta_x * delta_p);
        term_quantum = (hbar / sqrt_unc) * integral_psi * (2.0 * PI / t_Hubble);
        
        double V = compute_V();
        term_fluid = (rho_fluid * V * ug1_base) / M;
        
        term_osc = 2.0 * A_osc * std::cos(k_osc * x_pos) * std::cos(omega_osc * t) +
                   (2.0 * PI / t_Hubble_gyr) * A_osc * std::cos(k_osc * x_pos - omega_osc * t);
        
        double M_dm = M * M_DM_factor;
        term_DM = (M + M_dm) * (delta_rho_over_rho + 3.0 * G * M / (r * r * r)) / M;
        
        term_mag = compute_M_mag() / (M * r);
        term_decay = compute_cumulative_D(t) / (M * r);
    }
};

// ============================================================================
// UQFFModule13 - Main Module Wrapper
// ============================================================================
class UQFFModule13 {
private:
    MagnetarSGR1745_2900 magnetar;
    std::map<std::string, double> dynamicParams;
    bool enableLogging = true;
    double learningRate = 0.01;
    int updateCounter = 0;

public:
    UQFFModule13() = default;
    
    // === Standard Interface Methods ===
    
    void exportState(const std::string& filename) {
        std::ofstream out(filename);
        if (!out.is_open()) {
            std::cerr << "[UQFFModule13] Failed to export state to " << filename << std::endl;
            return;
        }
        out << "# UQFFModule13 State Export\n";
        out << "# Magnetar SGR 1745-2900 MUGE Calculator\n";
        out << "magnetar_mass=" << magnetar.getMass() << "\n";
        out << "magnetar_radius=" << magnetar.getRadius() << "\n";
        out << "magnetic_field=" << magnetar.getMagneticField() << "\n";
        out << "ug1_base=" << magnetar.getUg1Base() << "\n";
        out << "learning_rate=" << learningRate << "\n";
        out << "update_counter=" << updateCounter << "\n";
        out << "logging_enabled=" << (enableLogging ? "true" : "false") << "\n";
        out << "\n# Dynamic Parameters\n";
        for (const auto& [key, val] : dynamicParams) {
            out << "param." << key << "=" << std::scientific << val << "\n";
        }
        if (enableLogging) {
            std::cout << "[UQFFModule13] State exported to: " << filename << std::endl;
        }
    }
    
    void importState(const std::string& filename) {
        std::ifstream in(filename);
        if (!in.is_open()) {
            std::cerr << "[UQFFModule13] Failed to import state from " << filename << std::endl;
            return;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            size_t eq = line.find('=');
            if (eq == std::string::npos) continue;
            std::string key = line.substr(0, eq);
            std::string val = line.substr(eq + 1);
            if (key == "learning_rate") learningRate = std::stod(val);
            else if (key == "update_counter") updateCounter = std::stoi(val);
            else if (key == "logging_enabled") enableLogging = (val == "true");
            else if (key.substr(0, 6) == "param.") {
                dynamicParams[key.substr(6)] = std::stod(val);
            }
        }
        if (enableLogging) {
            std::cout << "[UQFFModule13] State imported from: " << filename << std::endl;
        }
    }
    
    void setDynamicParameter(const std::string& name, double value) {
        dynamicParams[name] = value;
        updateCounter++;
        if (enableLogging) {
            std::cout << "[UQFFModule13] Set param " << name << " = " << value << std::endl;
        }
    }
    
    double getDynamicParameter(const std::string& name) const {
        auto it = dynamicParams.find(name);
        return (it != dynamicParams.end()) ? it->second : 0.0;
    }
    
    void setEnableLogging(bool enable) { enableLogging = enable; }
    void setLearningRate(double rate) { learningRate = rate; }
    
    void printInfo() const {
        std::cout << "=== UQFFModule13 Info (source13.cpp) ===" << std::endl;
        std::cout << "Version: 2.0-Enhanced" << std::endl;
        std::cout << "Framework: Self-Expanding UQFF" << std::endl;
        std::cout << "Features: Magnetar SGR 1745-2900 MUGE (12 terms)" << std::endl;
        std::cout << "Magnetar Mass: " << std::scientific << magnetar.getMass() << " kg" << std::endl;
        std::cout << "Magnetar Radius: " << magnetar.getRadius() << " m" << std::endl;
        std::cout << "Magnetic Field: " << magnetar.getMagneticField() << " T" << std::endl;
        std::cout << "Base g (Ug1): " << magnetar.getUg1Base() << " m/s²" << std::endl;
        std::cout << "Dynamic Parameters: " << dynamicParams.size() << std::endl;
        std::cout << "Learning Rate: " << learningRate << std::endl;
        std::cout << "Update Counter: " << updateCounter << std::endl;
        std::cout << "Logging: " << (enableLogging ? "Enabled" : "Disabled") << std::endl;
    }

    // Self-simulation capability (aligned with self-expanding framework)
    template<typename Func>
    void runSimulation(double t_start, double t_end, int steps, Func computeFunc) {
        if (enableLogging) {
            std::cout << "[UQFFModule13] Running simulation: t=" << t_start << " to " << t_end << " (" << steps << " steps)" << std::endl;
        }
        double dt = (t_end - t_start) / steps;
        for (int i = 0; i <= steps; ++i) {
            double t = t_start + i * dt;
            double result = computeFunc(t);
            if (enableLogging) {
                std::cout << "[UQFFModule13] t=" << t << ": g=" << result << std::endl;
            }
        }
    }
    
    // === Physics Computations ===
    
    MagnetarSGR1745_2900& getMagnetar() { return magnetar; }
    const MagnetarSGR1745_2900& getMagnetar() const { return magnetar; }
    
    double compute_g_Magnetar(double t) const {
        return magnetar.compute_g_Magnetar(t);
    }
    
    void computeTimeEvolution(double t_start, double t_end, int steps) {
        if (enableLogging) {
            std::cout << "\n--- MUGE Time Evolution ---" << std::endl;
            std::cout << "Time range: " << t_start << " to " << t_end << " s" << std::endl;
        }
        
        double dt = (t_end - t_start) / steps;
        std::vector<double> times, g_values;
        
        for (int i = 0; i <= steps; ++i) {
            double t = t_start + i * dt;
            double g = magnetar.compute_g_Magnetar(t);
            times.push_back(t);
            g_values.push_back(g);
        }
        
        // Statistics
        double g_min = *std::min_element(g_values.begin(), g_values.end());
        double g_max = *std::max_element(g_values.begin(), g_values.end());
        double g_sum = std::accumulate(g_values.begin(), g_values.end(), 0.0);
        double g_mean = g_sum / g_values.size();
        
        std::cout << "g_min: " << std::scientific << g_min << " m/s²" << std::endl;
        std::cout << "g_max: " << g_max << " m/s²" << std::endl;
        std::cout << "g_mean: " << g_mean << " m/s²" << std::endl;
    }
    
    void printTermBreakdown(double t) const {
        double term1, term_BH, term_Ug, term_Lambda, term_EM, term_GW;
        double term_quantum, term_fluid, term_osc, term_DM, term_mag, term_decay;
        
        magnetar.computeTermBreakdown(t, term1, term_BH, term_Ug, term_Lambda, term_EM, term_GW,
                                       term_quantum, term_fluid, term_osc, term_DM, term_mag, term_decay);
        
        double total = term1 + term_BH + term_Ug + term_Lambda + term_EM + term_GW +
                       term_quantum + term_fluid + term_osc + term_DM + term_mag + term_decay;
        
        std::cout << "\n--- MUGE Term Breakdown at t=" << std::scientific << t << " s ---" << std::endl;
        std::cout << "  1. Base+H(z)+B:    " << term1 << " m/s²" << std::endl;
        std::cout << "  2. BH (Sgr A*):    " << term_BH << " m/s²" << std::endl;
        std::cout << "  3. UQFF Ug sum:    " << term_Ug << " m/s²" << std::endl;
        std::cout << "  4. Lambda (Λ):     " << term_Lambda << " m/s²" << std::endl;
        std::cout << "  5. EM (v×B):       " << term_EM << " m/s²" << std::endl;
        std::cout << "  6. GW (dΩ/dt)²:    " << term_GW << " m/s²" << std::endl;
        std::cout << "  7. Quantum:        " << term_quantum << " m/s²" << std::endl;
        std::cout << "  8. Fluid:          " << term_fluid << " m/s²" << std::endl;
        std::cout << "  9. Oscillatory:    " << term_osc << " m/s²" << std::endl;
        std::cout << " 10. DM+pert:        " << term_DM << " m/s²" << std::endl;
        std::cout << " 11. Magnetic E:     " << term_mag << " m/s²" << std::endl;
        std::cout << " 12. Decay E:        " << term_decay << " m/s²" << std::endl;
        std::cout << "     TOTAL g:        " << total << " m/s²" << std::endl;
    }
};

// ============================================================================
// Main Entry Point
// ============================================================================
int main() {
    std::cout << "=======================================" << std::endl;
    std::cout << "  Source13 Magnetar SGR 1745-2900 MUGE" << std::endl;
    std::cout << "  Version: 2.0-Enhanced (Module-Aligned)" << std::endl;
    std::cout << "=======================================" << std::endl;
    
    UQFFModule13 module;
    module.printInfo();
    
    // Time points for analysis
    double t_0 = 0.0;                              // Initial
    double t_1yr = 1.0 * 365.25 * 24 * 3600;       // 1 year
    double t_3_5yr = 3.5 * 365.25 * 24 * 3600;     // 3.5 years (decay timescale)
    double t_10yr = 10.0 * 365.25 * 24 * 3600;     // 10 years
    
    std::cout << "\n--- Computing g_Magnetar at key times ---" << std::endl;
    
    std::cout << "\nt = 0 (initial):" << std::endl;
    double g0 = module.compute_g_Magnetar(t_0);
    std::cout << "  g_Magnetar = " << std::scientific << g0 << " m/s²" << std::endl;
    module.printTermBreakdown(t_0);
    
    std::cout << "\nt = 1 year:" << std::endl;
    double g1 = module.compute_g_Magnetar(t_1yr);
    std::cout << "  g_Magnetar = " << std::scientific << g1 << " m/s²" << std::endl;
    
    std::cout << "\nt = 3.5 years (decay timescale):" << std::endl;
    double g3_5 = module.compute_g_Magnetar(t_3_5yr);
    std::cout << "  g_Magnetar = " << std::scientific << g3_5 << " m/s²" << std::endl;
    
    std::cout << "\nt = 10 years:" << std::endl;
    double g10 = module.compute_g_Magnetar(t_10yr);
    std::cout << "  g_Magnetar = " << std::scientific << g10 << " m/s²" << std::endl;
    module.printTermBreakdown(t_10yr);
    
    // Time evolution
    module.computeTimeEvolution(0, t_10yr, 100);
    
    // Compare with Newtonian
    auto& cfg = UQFFConfig13::getInstance();
    double g_newton = G * cfg.M_magnetar / (cfg.r_magnetar * cfg.r_magnetar);
    std::cout << "\n--- Comparison ---" << std::endl;
    std::cout << "g_Newton (GM/r²): " << std::scientific << g_newton << " m/s²" << std::endl;
    std::cout << "g_MUGE (t=0):     " << g0 << " m/s²" << std::endl;
    std::cout << "Enhancement:      " << (g0 / g_newton) << "x" << std::endl;
    
    // Export state
    module.exportState("source13_state.txt");
    
    // ==================== DUAL PHYSICS VALIDATION ====================
    std::cout << "\n=== Dual Physics Validation (FluidSolver + DualMethodValidator) ===" << std::endl;
    
    using namespace UQFFDualPhysics;
    
    // Initialize FluidSolver
    FluidSolver fluidSolver(32, 0.1, 0.0001);
    std::cout << "FluidSolver initialized (32x32 grid)" << std::endl;
    
    // Run fluid simulation with magnetar gravity
    fluidSolver.add_jet_force(10.0);
    for (int step = 0; step < 10; ++step) {
        fluidSolver.step(g0 * 1e-10);
    }
    std::cout << "FluidSolver: Max velocity = " << fluidSolver.getMaxVelocity() << " m/s" << std::endl;
    
    // Initialize DualMethodValidator
    DualMethodValidator validator("source13_dual_physics.log");
    validator.addConstraint("SGR1745-2900", 1e10, 1e14, 15.0);
    
    // Validate magnetar using dual physics methods
    UQFFDualPhysics::CelestialBody magnetar_body("SGR1745-2900", cfg.M_magnetar, cfg.r_magnetar, cfg.B0);
    magnetar_body.Pcore = 1.0 / cfg.P_init;
    
    UQFFDualPhysics::MUGESystem magnetar_muge("SGR1745-2900", cfg.M_magnetar, cfg.r_magnetar);
    magnetar_muge.B0 = cfg.B0;
    magnetar_muge.Lambda = cfg.Lambda;
    magnetar_muge.H0 = cfg.Hz * MPC_TO_M / 1000.0;
    
    auto result = validator.validate(magnetar_body, magnetar_muge, t_0, 0.0);
    result.print();
    
    std::cout << "Dual Physics: IMPLEMENTED" << std::endl;
    // ================================================================
    
    std::cout << "\n=== Source13 Complete ===" << std::endl;
    return 0;
}
