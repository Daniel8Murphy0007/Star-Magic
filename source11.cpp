// source11.cpp - UQFF Unified Field Calculator with Full Physics Terms
// Version: 2.0-Enhanced (Module-Aligned)
// Framework: Self-Expanding UQFF with compute_FU, Ug1-4, Ubi, Um
// Refactored: Removed graphics dependencies, added standard interface
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
// conflicts with source11's own struct definitions

// ============================================================================
// UQFFConfig11 - Singleton Configuration
// ============================================================================
class UQFFConfig11 {
public:
    static UQFFConfig11& getInstance() {
        static UQFFConfig11 instance;
        return instance;
    }
    
    // UQFF Physics Parameters
    double Omega_g = 7.3e-16;       // Galactic rotation frequency
    double Mbh = 8.15e36;           // Black hole mass (Sgr A*)
    double dg = 2.55e20;            // Galactic distance scale
    double v_SCm = 0.99 * c;        // SCm velocity (near light speed)
    double rho_A = 1e-23;           // Aether density
    double rho_sw = 8e-21;          // Solar wind density
    double v_sw = 5e5;              // Solar wind velocity
    double QA = 1e-10;              // Aether charge
    double kappa = 0.0005;          // Decay constant
    double alpha = 0.001;           // Time decay factor
    double gamma_decay = 0.00005;   // Magnetic decay
    double delta_sw = 0.01;         // Solar wind modifier
    double epsilon_sw = 0.001;      // Wind epsilon
    double delta_def = 0.01;        // Defect factor
    double HSCm = 1.0;              // SCm field strength
    double UUA = 1.0;               // UA coupling
    double eta = 1e-22;             // Metric perturbation
    double k1 = 1.5, k2 = 1.2, k3 = 1.8, k4 = 2.0;  // Coupling constants
    double beta_i = 0.6;            // Buoyancy coupling
    double rho_v = 6e-27;           // Vacuum density
    double C_concentration = 1.0;   // Concentration factor
    double f_feedback = 0.1;        // Feedback coefficient
    double num_strings = 1e9;       // String count
    double Ts00 = 1.27e3 + 1.11e7;  // Stress-energy T00
    
private:
    UQFFConfig11() = default;
    UQFFConfig11(const UQFFConfig11&) = delete;
    UQFFConfig11& operator=(const UQFFConfig11&) = delete;
};

// ============================================================================
// CelestialBody - Data structure matching bodies.csv
// ============================================================================
struct CelestialBody {
    std::string name;
    double Ms;          // Mass (kg)
    double Rs;          // Radius (m)
    double Rb;          // Boundary radius (m)
    double Ts_surface;  // Surface temperature (K)
    double omega_s;     // Surface angular velocity (rad/s)
    double Bs_avg;      // Average magnetic field (T)
    double SCm_density; // SCm density (kg/m³)
    double QUA;         // UA charge
    double Pcore;       // Core pressure factor
    double PSCm;        // SCm pressure factor
    double omega_c;     // Core angular velocity (rad/s)
};

// ============================================================================
// UQFFModule11 - Main Computation Module
// ============================================================================
class UQFFModule11 {
private:
    std::vector<CelestialBody> bodies;
    std::map<std::string, double> dynamicParams;
    bool enableLogging = true;
    double learningRate = 0.01;
    int updateCounter = 0;
    std::mt19937 rng;
    
    // Metric tensor (Minkowski background)
    std::vector<std::vector<double>> g_mu_nu = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, -1.0, 0.0, 0.0},
        {0.0, 0.0, -1.0, 0.0},
        {0.0, 0.0, 0.0, -1.0}
    };

public:
    UQFFModule11() : rng(std::random_device{}()) {}
    
    // === Standard Interface Methods ===
    
    void exportState(const std::string& filename) {
        std::ofstream out(filename);
        if (!out.is_open()) {
            std::cerr << "[UQFFModule11] Failed to export state to " << filename << std::endl;
            return;
        }
        out << "# UQFFModule11 State Export\n";
        out << "bodies_count=" << bodies.size() << "\n";
        out << "learning_rate=" << learningRate << "\n";
        out << "update_counter=" << updateCounter << "\n";
        out << "logging_enabled=" << (enableLogging ? "true" : "false") << "\n";
        out << "\n# Dynamic Parameters\n";
        for (const auto& [key, val] : dynamicParams) {
            out << "param." << key << "=" << std::scientific << val << "\n";
        }
        out << "\n# Loaded Bodies\n";
        for (const auto& b : bodies) {
            out << "body=" << b.name << ",Ms=" << b.Ms << ",Rs=" << b.Rs << "\n";
        }
        if (enableLogging) {
            std::cout << "[UQFFModule11] State exported to: " << filename << std::endl;
        }
    }
    
    void importState(const std::string& filename) {
        std::ifstream in(filename);
        if (!in.is_open()) {
            std::cerr << "[UQFFModule11] Failed to import state from " << filename << std::endl;
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
            std::cout << "[UQFFModule11] State imported from: " << filename << std::endl;
        }
    }
    
    void setDynamicParameter(const std::string& name, double value) {
        dynamicParams[name] = value;
        updateCounter++;
        if (enableLogging) {
            std::cout << "[UQFFModule11] Set param " << name << " = " << value << std::endl;
        }
    }
    
    double getDynamicParameter(const std::string& name) const {
        auto it = dynamicParams.find(name);
        return (it != dynamicParams.end()) ? it->second : 0.0;
    }
    
    void setEnableLogging(bool enable) { enableLogging = enable; }
    void setLearningRate(double rate) { learningRate = rate; }
    
    void printInfo() const {
        std::cout << "=== UQFFModule11 Info (source11.cpp) ===" << std::endl;
        std::cout << "Version: 2.0-Enhanced" << std::endl;
        std::cout << "Framework: Self-Expanding UQFF" << std::endl;
        std::cout << "Features: compute_FU, Ug1-4, Ubi, Um, A_mu_nu" << std::endl;
        std::cout << "Bodies Loaded: " << bodies.size() << std::endl;
        std::cout << "Dynamic Parameters: " << dynamicParams.size() << std::endl;
        std::cout << "Learning Rate: " << learningRate << std::endl;
        std::cout << "Update Counter: " << updateCounter << std::endl;
        std::cout << "Logging: " << (enableLogging ? "Enabled" : "Disabled") << std::endl;
    }

    // Self-simulation capability (aligned with self-expanding framework)
    template<typename Func>
    void runSimulation(double t_start, double t_end, int steps, Func computeFunc) {
        if (enableLogging) {
            std::cout << "[UQFFModule11] Running simulation: t=" << t_start << " to " << t_end << " (" << steps << " steps)" << std::endl;
        }
        double dt = (t_end - t_start) / steps;
        for (int i = 0; i <= steps; ++i) {
            double t = t_start + i * dt;
            double result = computeFunc(t);
            if (enableLogging) {
                std::cout << "[UQFFModule11] t=" << t << ": g=" << result << std::endl;
            }
        }
    }
    
    // === Data Loading ===
    
    bool loadBodies(const std::string& filename) {
        std::ifstream in(filename);
        if (!in.is_open()) {
            std::cerr << "[UQFFModule11] Failed to open: " << filename << std::endl;
            return false;
        }
        bodies.clear();
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);
            CelestialBody body;
            std::string token;
            
            std::getline(ss, body.name, ',');
            if (body.name == "name") continue;  // Skip header
            
            std::getline(ss, token, ','); body.Ms = std::stod(token);
            std::getline(ss, token, ','); body.Rs = std::stod(token);
            std::getline(ss, token, ','); body.Rb = std::stod(token);
            std::getline(ss, token, ','); body.Ts_surface = std::stod(token);
            std::getline(ss, token, ','); body.omega_s = std::stod(token);
            std::getline(ss, token, ','); body.Bs_avg = std::stod(token);
            std::getline(ss, token, ','); body.SCm_density = std::stod(token);
            std::getline(ss, token, ','); body.QUA = std::stod(token);
            std::getline(ss, token, ','); body.Pcore = std::stod(token);
            std::getline(ss, token, ','); body.PSCm = std::stod(token);
            std::getline(ss, token, ','); body.omega_c = std::stod(token);
            
            bodies.push_back(body);
        }
        if (enableLogging) {
            std::cout << "[UQFFModule11] Loaded " << bodies.size() << " bodies from " << filename << std::endl;
        }
        return !bodies.empty();
    }
    
    const std::vector<CelestialBody>& getBodies() const { return bodies; }
    
    // === UQFF Physics Computations ===
    
    // Step function for boundary conditions
    double step_function(double r, double Rb) const {
        return (r > Rb) ? 1.0 : 0.0;
    }
    
    // Reaction energy from SCm-Aether interaction
    double compute_Ereact(double t, double rho_SCm) const {
        auto& cfg = UQFFConfig11::getInstance();
        if (cfg.rho_A <= 0.0) return 0.0;
        return (rho_SCm * cfg.v_SCm * cfg.v_SCm / cfg.rho_A) * std::exp(-cfg.kappa * t);
    }
    
    // Magnetic moment from field and radius
    double compute_mu_s(double t, double Bs, double omega_c, double Rs, double SCm_contrib = 1e3) const {
        double Bs_t = Bs + 0.4 * std::sin(omega_c * t) + SCm_contrib;
        return Bs_t * std::pow(Rs, 3);
    }
    
    // Gravitational gradient
    double compute_grad_Ms_r(double Ms, double Rs) const {
        if (Rs <= 0.0) return 0.0;
        return G * Ms / (Rs * Rs);
    }
    
    // Jet magnetic field
    double compute_Bj(double t, double omega_c, double SCm_contrib = 1e3) const {
        return 1e-3 + 0.4 * std::sin(omega_c * t) + SCm_contrib;
    }
    
    // Time-dependent angular velocity
    double compute_omega_s_t(double t, double omega_s, double omega_c) const {
        return omega_s - 0.4e-6 * std::sin(omega_c * t);
    }
    
    // Jet magnetic moment
    double compute_mu_j(double t, double omega_c, double Rs, double SCm_contrib = 1e3) const {
        double Bj = compute_Bj(t, omega_c, SCm_contrib);
        return Bj * std::pow(Rs, 3);
    }
    
    // Ug1: Magnetic dipole contribution
    double compute_Ug1(const CelestialBody& body, double r, double t, double tn) const {
        if (r <= 0.0) return 0.0;
        auto& cfg = UQFFConfig11::getInstance();
        double mu_s = compute_mu_s(t, body.Bs_avg, body.omega_c, body.Rs);
        double grad_Ms_r = compute_grad_Ms_r(body.Ms, body.Rs);
        double defect = 1.0 + cfg.delta_def * std::sin(0.001 * t);
        return cfg.k1 * mu_s * grad_Ms_r * std::exp(-cfg.alpha * t) * std::cos(PI * tn) * defect;
    }
    
    // Ug2: Charge-reactivity contribution
    double compute_Ug2(const CelestialBody& body, double r, double t, double tn) const {
        if (r <= 0.0) return 0.0;
        auto& cfg = UQFFConfig11::getInstance();
        double Ereact = compute_Ereact(t, body.SCm_density);
        double S = step_function(r, body.Rb);
        double wind_mod = 1.0 + cfg.delta_sw * cfg.v_sw;
        return cfg.k2 * (cfg.QA + body.QUA) * body.Ms / (r * r) * S * wind_mod * cfg.HSCm * Ereact;
    }
    
    // Ug3: String rotation contribution
    double compute_Ug3(const CelestialBody& body, double r, double t, double tn, double theta) const {
        auto& cfg = UQFFConfig11::getInstance();
        double Ereact = compute_Ereact(t, body.SCm_density);
        double omega_s_t = compute_omega_s_t(t, body.omega_s, body.omega_c);
        double Bj = compute_Bj(t, body.omega_c);
        return cfg.k3 * Bj * std::cos(omega_s_t * t * PI) * body.Pcore * Ereact;
    }
    
    // Ug4: Vacuum concentration contribution
    double compute_Ug4(double t, double tn) const {
        auto& cfg = UQFFConfig11::getInstance();
        if (cfg.dg <= 0.0) return 0.0;
        double decay = std::exp(-cfg.alpha * t);
        double cycle = std::cos(PI * tn);
        return cfg.k4 * cfg.rho_v * cfg.C_concentration * cfg.Mbh / cfg.dg * decay * cycle * (1 + cfg.f_feedback);
    }
    
    // Um: Magnetism contribution from string ensemble
    double compute_Um(const CelestialBody& body, double t, double tn, double rj, double phi_hat = 1.0) const {
        if (rj <= 0.0) return 0.0;
        auto& cfg = UQFFConfig11::getInstance();
        double Ereact = compute_Ereact(t, body.SCm_density);
        double mu_j = compute_mu_j(t, body.omega_c, body.Rs);
        double decay = 1.0 - std::exp(-cfg.gamma_decay * t * std::cos(PI * tn));
        double single = mu_j / rj * decay * phi_hat;
        return single * cfg.num_strings * body.PSCm * Ereact;
    }
    
    // Ubi: Buoyancy force (negative = attractive, positive = repulsive)
    double compute_Ubi(double Ugi, double tn) const {
        auto& cfg = UQFFConfig11::getInstance();
        if (cfg.dg <= 0.0) return 0.0;
        double wind_mod = 1.0 + cfg.epsilon_sw * cfg.rho_sw;
        return -cfg.beta_i * Ugi * cfg.Omega_g * cfg.Mbh / cfg.dg * wind_mod * cfg.UUA * std::cos(PI * tn);
    }
    
    // Metric perturbation tensor A_mu_nu
    std::vector<std::vector<double>> compute_A_mu_nu(double tn) const {
        auto& cfg = UQFFConfig11::getInstance();
        std::vector<std::vector<double>> A = g_mu_nu;
        double mod = cfg.eta * cfg.Ts00 * std::cos(PI * tn);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                A[i][j] += mod;
            }
        }
        return A;
    }
    
    // Complete Unified Field F_U
    double compute_FU(const CelestialBody& body, double r, double t, double tn, double theta) const {
        // Compute gravity components
        double Ug1 = compute_Ug1(body, r, t, tn);
        double Ug2 = compute_Ug2(body, r, t, tn);
        double Ug3 = compute_Ug3(body, r, t, tn, theta);
        double Ug4 = compute_Ug4(t, tn);
        double sum_Ugi = Ug1 + Ug2 + Ug3 + Ug4;
        
        // Compute buoyancy forces
        double Ubi1 = compute_Ubi(Ug1, tn);
        double Ubi2 = compute_Ubi(Ug2, tn);
        double Ubi3 = compute_Ubi(Ug3, tn);
        double Ubi4 = compute_Ubi(Ug4, tn);
        double sum_Ubi = Ubi1 + Ubi2 + Ubi3 + Ubi4;
        
        // Magnetic contribution
        double Um = compute_Um(body, t, tn, body.Rb);
        
        // Metric tensor scalar (trace of perturbation)
        auto A = compute_A_mu_nu(tn);
        double A_scalar = A[0][0] + A[1][1] + A[2][2] + A[3][3];
        
        return sum_Ugi + sum_Ubi + Um + A_scalar;
    }
    
    // Newtonian gravity for comparison
    double compute_g_newton(const CelestialBody& body, double r) const {
        if (r <= 0.0) return 0.0;
        return G * body.Ms / (r * r);
    }
    
    // Compute all fields for a body
    void computeAllFields(const CelestialBody& body, double t, double tn, double theta) {
        double r = body.Rs * 10.0;  // Evaluate at 10 radii
        
        double Ug1 = compute_Ug1(body, r, t, tn);
        double Ug2 = compute_Ug2(body, r, t, tn);
        double Ug3 = compute_Ug3(body, r, t, tn, theta);
        double Ug4 = compute_Ug4(t, tn);
        double Um = compute_Um(body, t, tn, body.Rb);
        double Ubi_total = compute_Ubi(Ug1 + Ug2 + Ug3 + Ug4, tn);
        double FU = compute_FU(body, r, t, tn, theta);
        double g_N = compute_g_newton(body, r);
        
        std::cout << "\nSystem: " << body.name << std::endl;
        std::cout << "  Mass: " << std::scientific << body.Ms << " kg" << std::endl;
        std::cout << "  Radius: " << body.Rs << " m" << std::endl;
        std::cout << "  Eval distance: " << r << " m (10 Rs)" << std::endl;
        std::cout << "  Ug1 (magnetic): " << Ug1 << std::endl;
        std::cout << "  Ug2 (charge): " << Ug2 << std::endl;
        std::cout << "  Ug3 (rotation): " << Ug3 << std::endl;
        std::cout << "  Ug4 (vacuum): " << Ug4 << std::endl;
        std::cout << "  Um (magnetism): " << Um << std::endl;
        std::cout << "  Ubi (buoyancy): " << Ubi_total << std::endl;
        std::cout << "  F_U (unified): " << FU << std::endl;
        std::cout << "  g_Newton: " << g_N << " m/s²" << std::endl;
        
        // Buoyancy interpretation
        if (Ubi_total < 0) {
            std::cout << "  [ATTRACTIVE buoyancy - enhanced gravity]" << std::endl;
        } else if (Ubi_total > 0) {
            std::cout << "  [REPULSIVE buoyancy - reduced gravity]" << std::endl;
        }
    }
};

// ============================================================================
// Main Entry Point
// ============================================================================
int main() {
    std::cout << "=======================================" << std::endl;
    std::cout << "  Source11 UQFF Unified Field Calculator" << std::endl;
    std::cout << "  Version: 2.0-Enhanced (Module-Aligned)" << std::endl;
    std::cout << "=======================================" << std::endl;
    
    UQFFModule11 module;
    module.printInfo();
    
    // Load celestial bodies
    if (!module.loadBodies("bodies.csv")) {
        std::cerr << "Failed to load bodies.csv" << std::endl;
        return 1;
    }
    
    std::cout << "\n--- Available Bodies ---" << std::endl;
    for (const auto& body : module.getBodies()) {
        std::cout << "  " << body.name << " (M=" << std::scientific << body.Ms << " kg)" << std::endl;
    }
    
    // Computation parameters
    double t = 1.0;       // Time (s)
    double tn = 0.5;      // Normalized time cycle
    double theta = 0.0;   // Angular position
    
    std::cout << "\n--- Computing UQFF for all bodies ---" << std::endl;
    std::cout << "Time: t=" << t << " s, tn=" << tn << ", theta=" << theta << std::endl;
    
    std::vector<double> FU_values;
    for (const auto& body : module.getBodies()) {
        module.computeAllFields(body, t, tn, theta);
        double r = body.Rs * 10.0;
        FU_values.push_back(module.compute_FU(body, r, t, tn, theta));
    }
    
    // Statistics
    if (!FU_values.empty()) {
        double min_FU = *std::min_element(FU_values.begin(), FU_values.end());
        double max_FU = *std::max_element(FU_values.begin(), FU_values.end());
        double sum_FU = std::accumulate(FU_values.begin(), FU_values.end(), 0.0);
        double mean_FU = sum_FU / FU_values.size();
        
        std::cout << "\n--- F_U Summary ---" << std::endl;
        std::cout << "Min: " << std::scientific << min_FU << std::endl;
        std::cout << "Max: " << max_FU << std::endl;
        std::cout << "Mean: " << mean_FU << std::endl;
    }
    
    // Export state
    module.exportState("source11_state.txt");
    
    // ==================== DUAL PHYSICS VALIDATION ====================
    std::cout << "\n=== Dual Physics Validation (FluidSolver + DualMethodValidator) ===" << std::endl;
    
    using namespace UQFFDualPhysics;
    
    // Initialize FluidSolver
    FluidSolver fluidSolver(32, 0.1, 0.0001);
    std::cout << "FluidSolver initialized (32x32 grid)" << std::endl;
    
    // Run fluid simulation
    double test_g = FU_values.empty() ? 9.81 : std::abs(FU_values[0]);
    fluidSolver.add_jet_force(10.0);
    for (int step = 0; step < 10; ++step) {
        fluidSolver.step(test_g * 1e-10);
    }
    std::cout << "FluidSolver: Max velocity = " << fluidSolver.getMaxVelocity() << " m/s" << std::endl;
    
    // Initialize DualMethodValidator
    DualMethodValidator validator("source11_dual_physics.log");
    
    // Validate first body
    if (!module.getBodies().empty()) {
        const auto& body0 = module.getBodies()[0];
        
        UQFFDualPhysics::CelestialBody uqff_body(body0.name, body0.Ms, body0.Rs, body0.Bs_avg);
        uqff_body.SCm_density = body0.SCm_density;
        uqff_body.QUA = body0.QUA;
        uqff_body.omega_c = body0.omega_c;
        
        UQFFDualPhysics::MUGESystem muge(body0.name, body0.Ms, body0.Rb);
        muge.B0 = body0.Bs_avg;
        muge.omega = body0.omega_s;
        muge.T = body0.Ts_surface;
        
        auto result = validator.validate(uqff_body, muge, t, theta);
        result.print();
    }
    
    std::cout << "Dual Physics: IMPLEMENTED" << std::endl;
    // ================================================================
    
    std::cout << "\n=== Source11 Complete ===" << std::endl;
    return 0;
}
