// source9.cpp - UQFF Unified Field Calculator with CoAnQi Namespace Organization
// 2.0-Enhanced with Self-Expanding Framework
// Refactored for module alignment with source2/4/5/200/6/7
// Enhanced: January 22, 2026 - Added FluidSolver, DualMethodValidator, DualPhysicsMethods

#include "uqff_constants.h"
#include "uqff_self_expanding.h"
#include "uqff_dual_physics.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <memory>
#include <chrono>

using namespace UQFF;
// Note: UQFFDualPhysics and UQFFExpanding used locally to avoid
// conflicts with source9's own struct definitions

// =============================================================================
// UQFFConfig9 Singleton - Physics Parameters
// =============================================================================
struct UQFFConfig9 {
    // SCm and Aether parameters
    double v_SCm = 0.99 * c;
    double rho_A = 1e-23;
    double rho_sw = 8e-21;
    double v_sw = 5e5;
    double QA = 1e-10;
    double Qs = 0.0;
    double kappa = 0.0005;
    double alpha = 0.001;
    double gamma_param = 0.00005;
    double delta_sw = 0.01;
    double epsilon_sw = 0.001;
    double delta_def = 0.01;
    double HSCm = 1.0;
    double UUA = 1.0;
    double eta = 1e-22;
    
    // Coupling constants
    double k1 = 1.5;
    double k2 = 1.2;
    double k3 = 1.8;
    double k4 = 2.0;
    
    // Buoyancy and feedback
    double beta_i = 0.6;
    double rho_v = 6e-27;
    double C_concentration = 1.0;
    double f_feedback = 0.1;
    double num_strings = 1e9;
    
    // Black hole parameters
    double Omega_g = 7.3e-16;
    double Mbh = 8.15e36;
    double dg = 2.55e20;
    
    // Stress-energy
    double Ts00 = 1.27e3 + 1.11e7;
    
    // Metric tensor (Minkowski)
    std::vector<std::vector<double>> g_mu_nu = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, -1.0, 0.0, 0.0},
        {0.0, 0.0, -1.0, 0.0},
        {0.0, 0.0, 0.0, -1.0}
    };
    
    static UQFFConfig9& instance() {
        static UQFFConfig9 inst;
        return inst;
    }
private:
    UQFFConfig9() = default;
};

// =============================================================================
// CelestialBody Structure
// =============================================================================
struct CelestialBody {
    std::string name;
    double Ms;           // Mass (kg)
    double Rs;           // Radius (m)
    double Rb;           // Boundary radius (m)
    double Ts_surface;   // Surface temperature (K)
    double omega_s;      // Angular velocity (rad/s)
    double Bs_avg;       // Average magnetic field (T)
    double SCm_density;  // SCm density
    double QUA;          // Charge parameter
    double Pcore;        // Core pressure
    double PSCm;         // SCm pressure
    double omega_c;      // Cyclic frequency (rad/s)
};

// =============================================================================
// MUGESystem Structure
// =============================================================================
struct MUGESystem {
    std::string name;
    double I;            // Moment of inertia
    double A;            // Area
    double omega1;       // Angular frequency 1
    double omega2;       // Angular frequency 2
    double Vsys;         // System volume
    double vexp;         // Expansion velocity
    double t;            // Time
    double z;            // Redshift
    double ffluid;       // Fluid frequency
    double M;            // Mass
    double r;            // Radius
    double B;            // Magnetic field
    double Bcrit;        // Critical magnetic field
    double rho_fluid;    // Fluid density
    double g_local;      // Local gravity
    double M_DM;         // Dark matter mass
    double delta_rho_rho;// Density perturbation
};

// =============================================================================
// ResonanceParams Structure
// =============================================================================
struct ResonanceParams {
    double fDPM = 1e12;
    double fTHz = 1e12;
    double Evac_neb = 7.09e-36;
    double Evac_ISM = 7.09e-37;
    double Delta_Evac = 6.381e-36;
    double Fsuper = 6.287e-19;
    double UA_SCM = 10;
    double omega_i = 1e-8;
    double k4_res = 1.0;
    double freact = 1e10;
    double fquantum = 1.445e-17;
    double fAether = 1.576e-35;
    double fosc = 4.57e14;
    double fTRZ = 0.1;
    double c_res = 3e8;
};

// =============================================================================
// UQFFModule9 - Self-Expanding Module Class (CoAnQi Namespace Organization)
// =============================================================================
class UQFFModule9 {
private:
    std::string version = "2.0-Enhanced";
    bool enable_logging = true;
    double learning_rate = 0.01;
    int update_counter = 0;
    std::map<std::string, double> dynamic_parameters;
    std::map<std::string, std::string> metadata;
    
public:
    UQFFModule9() {
        metadata["module"] = "source9";
        metadata["framework"] = "Self-Expanding UQFF";
        metadata["features"] = "CoAnQi,Physics,MUGE,Namespace";
    }
    
    // =========================================================================
    // Helper Functions
    // =========================================================================
    double step_function(double r, double Rb) const {
        return (r > Rb) ? 1.0 : 0.0;
    }
    
    double compute_Ereact(double t, double rho_SCm, double v_SCm, double rho_A, double kappa) const {
        if (rho_A <= 0.0) throw std::runtime_error("Invalid rho_A value");
        return (rho_SCm * v_SCm * v_SCm / rho_A) * std::exp(-kappa * t);
    }
    
    double compute_mu_s(double t, double Bs, double omega_c, double Rs, double SCm_contrib = 1e3) const {
        double Bs_t = Bs + 0.4 * std::sin(omega_c * t) + SCm_contrib;
        return Bs_t * std::pow(Rs, 3);
    }
    
    double compute_grad_Ms_r(double Ms, double Rs) const {
        if (Rs <= 0.0) throw std::runtime_error("Invalid Rs value");
        return G * Ms / (Rs * Rs);
    }
    
    double compute_Bj(double t, double omega_c, double SCm_contrib = 1e3) const {
        return 1e-3 + 0.4 * std::sin(omega_c * t) + SCm_contrib;
    }
    
    double compute_omega_s_t(double t, double omega_s, double omega_c) const {
        return omega_s - 0.4e-6 * std::sin(omega_c * t);
    }
    
    double compute_mu_j(double t, double omega_c, double Rs, double SCm_contrib = 1e3) const {
        double Bj = compute_Bj(t, omega_c, SCm_contrib);
        return Bj * std::pow(Rs, 3);
    }
    
    // =========================================================================
    // UQFF Core Physics Computations (CoAnQi::Physics)
    // =========================================================================
    double compute_Ug1(const CelestialBody& body, double r, double t, double tn) const {
        auto& cfg = UQFFConfig9::instance();
        if (r <= 0.0) throw std::runtime_error("Invalid r value");
        double mu_s = compute_mu_s(t, body.Bs_avg, body.omega_c, body.Rs);
        double grad_Ms_r = compute_grad_Ms_r(body.Ms, body.Rs);
        double defect = 1.0 + cfg.delta_def * std::sin(0.001 * t);
        return cfg.k1 * mu_s * grad_Ms_r * std::exp(-cfg.alpha * t) * std::cos(PI * tn) * defect;
    }
    
    double compute_Ug2(const CelestialBody& body, double r, double t, double tn) const {
        auto& cfg = UQFFConfig9::instance();
        if (r <= 0.0) throw std::runtime_error("Invalid r value");
        double Ereact = compute_Ereact(t, body.SCm_density, cfg.v_SCm, cfg.rho_A, cfg.kappa);
        double S = step_function(r, body.Rb);
        double wind_mod = 1.0 + cfg.delta_sw * cfg.v_sw;
        return cfg.k2 * (cfg.QA + body.QUA) * body.Ms / (r * r) * S * wind_mod * cfg.HSCm * Ereact;
    }
    
    double compute_Ug3(const CelestialBody& body, double r, double t, double tn, double theta) const {
        auto& cfg = UQFFConfig9::instance();
        double Ereact = compute_Ereact(t, body.SCm_density, cfg.v_SCm, cfg.rho_A, cfg.kappa);
        double omega_s_t = compute_omega_s_t(t, body.omega_s, body.omega_c);
        double Bj = compute_Bj(t, body.omega_c);
        return cfg.k3 * Bj * std::cos(omega_s_t * t * PI) * body.Pcore * Ereact;
    }
    
    double compute_Ug4(double t, double tn) const {
        auto& cfg = UQFFConfig9::instance();
        if (cfg.dg <= 0.0) throw std::runtime_error("Invalid dg value");
        double decay = std::exp(-cfg.alpha * t);
        double cycle = std::cos(PI * tn);
        return cfg.k4 * cfg.rho_v * cfg.C_concentration * cfg.Mbh / cfg.dg * decay * cycle * (1 + cfg.f_feedback);
    }
    
    double compute_Ubi(double Ugi, double tn) const {
        auto& cfg = UQFFConfig9::instance();
        if (cfg.dg <= 0.0) throw std::runtime_error("Invalid dg value");
        double wind_mod = 1.0 + cfg.epsilon_sw * cfg.rho_sw;
        return -cfg.beta_i * Ugi * cfg.Omega_g * cfg.Mbh / cfg.dg * wind_mod * cfg.UUA * std::cos(PI * tn);
    }
    
    double compute_Um(const CelestialBody& body, double t, double tn, double rj) const {
        auto& cfg = UQFFConfig9::instance();
        if (rj <= 0.0) throw std::runtime_error("Invalid rj value");
        double Ereact = compute_Ereact(t, body.SCm_density, cfg.v_SCm, cfg.rho_A, cfg.kappa);
        double mu_j = compute_mu_j(t, body.omega_c, body.Rs);
        double decay = 1.0 - std::exp(-cfg.gamma_param * t * std::cos(PI * tn));
        double single = mu_j / rj * decay;
        return single * cfg.num_strings * body.PSCm * Ereact;
    }
    
    std::vector<std::vector<double>> compute_A_mu_nu(double tn) const {
        auto& cfg = UQFFConfig9::instance();
        std::vector<std::vector<double>> A = cfg.g_mu_nu;
        double mod = cfg.eta * cfg.Ts00 * std::cos(PI * tn);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                A[i][j] += mod;
            }
        }
        return A;
    }
    
    double compute_FU(const CelestialBody& body, double r, double t, double tn, double theta) {
        if (enable_logging) {
            std::cout << "[UQFFModule9] Computing FU for " << body.name << std::endl;
        }
        
        double Ug1 = compute_Ug1(body, r, t, tn);
        double Ug2 = compute_Ug2(body, r, t, tn);
        double Ug3 = compute_Ug3(body, r, t, tn, theta);
        double Ug4 = compute_Ug4(t, tn);
        double sum_Ugi = Ug1 + Ug2 + Ug3 + Ug4;
        
        double Ubi1 = compute_Ubi(Ug1, tn);
        double Ubi2 = compute_Ubi(Ug2, tn);
        double Ubi3 = compute_Ubi(Ug3, tn);
        double Ubi4 = compute_Ubi(Ug4, tn);
        double sum_Ubi = Ubi1 + Ubi2 + Ubi3 + Ubi4;
        
        double Um = compute_Um(body, t, tn, body.Rb);
        
        auto A = compute_A_mu_nu(tn);
        double A_scalar = A[0][0] + A[1][1] + A[2][2] + A[3][3];
        
        update_counter++;
        return sum_Ugi + sum_Ubi + Um + A_scalar;
    }
    
    // =========================================================================
    // MUGE Computations (CoAnQi::MUGE)
    // =========================================================================
    double compute_compressed_base(const MUGESystem& sys) const {
        if (sys.r <= 0.0) throw std::runtime_error("Invalid r value");
        return G * sys.M / (sys.r * sys.r);
    }
    
    double compute_compressed_expansion(const MUGESystem& sys, double H0 = 2.269e-18) const {
        return std::exp(H0 * sys.t);
    }
    
    double compute_compressed_super_adj(const MUGESystem& sys) const {
        if (sys.Bcrit <= 0.0) return 1.0;
        return 1.0 - sys.B / sys.Bcrit;
    }
    
    double compute_compressed_MUGE(const MUGESystem& sys) {
        if (enable_logging) {
            std::cout << "[UQFFModule9] Computing compressed MUGE for " << sys.name << std::endl;
        }
        double base = compute_compressed_base(sys);
        double expansion = compute_compressed_expansion(sys);
        double super_adj = compute_compressed_super_adj(sys);
        return base * expansion * super_adj;
    }
    
    double compute_aDPM(const MUGESystem& sys, const ResonanceParams& res) const {
        double FDPM = sys.I * sys.A * (sys.omega1 - sys.omega2);
        return FDPM * res.fDPM * res.Evac_neb * res.c_res * sys.Vsys;
    }
    
    double compute_resonance_MUGE(const MUGESystem& sys, const ResonanceParams& res) {
        if (enable_logging) {
            std::cout << "[UQFFModule9] Computing resonance MUGE for " << sys.name << std::endl;
        }
        double aDPM = compute_aDPM(sys, res);
        double afluid = res.freact * sys.ffluid * sys.Vsys;
        return aDPM + afluid;
    }
    
    // =========================================================================
    // Self-Expanding Framework Methods
    // =========================================================================
    void exportState(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out.is_open()) {
            std::cerr << "[UQFFModule9] Failed to open file for export: " << filename << std::endl;
            return;
        }
        out << "# UQFFModule9 State Export\n";
        out << "version=" << version << "\n";
        out << "enable_logging=" << (enable_logging ? "true" : "false") << "\n";
        out << "learning_rate=" << learning_rate << "\n";
        out << "update_counter=" << update_counter << "\n";
        out << "# Dynamic Parameters\n";
        for (const auto& kv : dynamic_parameters) {
            out << "param:" << kv.first << "=" << kv.second << "\n";
        }
        out << "# Metadata\n";
        for (const auto& kv : metadata) {
            out << "meta:" << kv.first << "=" << kv.second << "\n";
        }
        out.close();
        if (enable_logging) {
            std::cout << "[UQFFModule9] State exported to: " << filename << std::endl;
        }
    }
    
    void importState(const std::string& filename) {
        std::ifstream in(filename);
        if (!in.is_open()) {
            std::cerr << "[UQFFModule9] Failed to open file for import: " << filename << std::endl;
            return;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            size_t eq = line.find('=');
            if (eq == std::string::npos) continue;
            std::string key = line.substr(0, eq);
            std::string val = line.substr(eq + 1);
            if (key == "learning_rate") learning_rate = std::stod(val);
            else if (key == "enable_logging") enable_logging = (val == "true");
            else if (key.substr(0, 6) == "param:") {
                dynamic_parameters[key.substr(6)] = std::stod(val);
            }
        }
        in.close();
        if (enable_logging) {
            std::cout << "[UQFFModule9] State imported from: " << filename << std::endl;
        }
    }
    
    void setDynamicParameter(const std::string& name, double value) {
        dynamic_parameters[name] = value;
        if (enable_logging) {
            std::cout << "[UQFFModule9] Set parameter " << name << " = " << value << std::endl;
        }
    }
    
    double getDynamicParameter(const std::string& name, double default_val = 0.0) const {
        auto it = dynamic_parameters.find(name);
        return (it != dynamic_parameters.end()) ? it->second : default_val;
    }
    
    void setEnableLogging(bool enable) {
        enable_logging = enable;
    }
    
    void printInfo() const {
        std::cout << "\n=== UQFFModule9 Info (source9.cpp) ===" << std::endl;
        std::cout << "Version: " << version << std::endl;
        std::cout << "Framework: Self-Expanding UQFF" << std::endl;
        std::cout << "Features: CoAnQi Namespace, Physics, MUGE" << std::endl;
        std::cout << "Dynamic Parameters: " << dynamic_parameters.size() << std::endl;
        std::cout << "Learning Rate: " << learning_rate << std::endl;
        std::cout << "Update Counter: " << update_counter << std::endl;
        std::cout << "Logging: " << (enable_logging ? "Enabled" : "Disabled") << std::endl;
    }

    // Self-simulation capability (aligned with self-expanding framework)
    template<typename Func>
    void runSimulation(double t_start, double t_end, int steps, Func computeFunc) {
        if (enable_logging) {
            std::cout << "[UQFFModule9] Running simulation: t=" << t_start << " to " << t_end << " (" << steps << " steps)" << std::endl;
        }
        double dt = (t_end - t_start) / steps;
        for (int i = 0; i <= steps; ++i) {
            double t = t_start + i * dt;
            double result = computeFunc(t);
            if (enable_logging) {
                std::cout << "[UQFFModule9] t=" << t << ": g=" << result << std::endl;
            }
        }
    }
    
    // =========================================================================
    // JSON Parameter Output (CoAnQi::Physics::output_json_params)
    // =========================================================================
    void output_json_params(const CelestialBody& body) const {
        auto& cfg = UQFFConfig9::instance();
        std::cout << "{" << std::endl;
        std::cout << "  \"name\": \"" << body.name << "\"," << std::endl;
        std::cout << "  \"SCm_density\": " << body.SCm_density << "," << std::endl;
        std::cout << "  \"UA\": " << body.QUA << "," << std::endl;
        std::cout << "  \"Qs\": " << cfg.Qs << std::endl;
        std::cout << "}" << std::endl;
    }
};

// =============================================================================
// CSV Body Loading
// =============================================================================
std::vector<CelestialBody> load_bodies_csv(const std::string& filename) {
    std::vector<CelestialBody> bodies;
    std::ifstream in(filename);
    if (!in.is_open()) {
        std::cerr << "Warning: Could not open " << filename << std::endl;
        return bodies;
    }
    std::string line;
    bool skipped_header = false;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (!skipped_header) { skipped_header = true; continue; }
        try {
            std::stringstream ss(line);
            CelestialBody body;
            std::string token;
            std::getline(ss, body.name, ',');
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
        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to parse line: " << line << std::endl;
        }
    }
    return bodies;
}

// =============================================================================
// Default Bodies
// =============================================================================
std::vector<CelestialBody> get_default_bodies() {
    std::vector<CelestialBody> bodies;
    bodies.push_back({"Sun", 1.989e30, 6.96e8, 1.496e13, 5778.0, 2.5e-6, 1e-4, 1e15, 1e-11, 1.0, 1.0, 
                      2 * PI / (11.0 * 365.25 * 24 * 3600)});
    bodies.push_back({"Earth", 5.972e24, 6.371e6, 1e7, 288.0, 7.292e-5, 3e-5, 1e12, 1e-12, 1e-3, 1e-3, 
                      2 * PI / (365.25 * 24 * 3600)});
    bodies.push_back({"Jupiter", 1.898e27, 6.9911e7, 1e8, 165.0, 1.76e-4, 4e-4, 1e13, 1e-11, 1e-3, 1e-3, 
                      2 * PI / (11.86 * 365.25 * 24 * 3600)});
    bodies.push_back({"Neptune", 1.024e26, 2.4622e7, 5e7, 72.0, 1.08e-4, 1e-4, 1e11, 1e-13, 1e-3, 1e-3, 
                      2 * PI / (164.8 * 365.25 * 24 * 3600)});
    return bodies;
}

// =============================================================================
// Default MUGE Systems
// =============================================================================
std::vector<MUGESystem> get_default_muge_systems() {
    std::vector<MUGESystem> systems;
    systems.push_back({"SGR 1745-2900", 1e21, 3.142e8, 1e-3, -1e-3, 4.189e12, 1e3, 3.799e10, 0.0009, 
                       1.269e-14, 2.984e30, 1e4, 1e10, 1e11, 1e-15, 10.0, 0.0, 1e-5});
    systems.push_back({"Sagittarius A*", 1e23, 2.813e30, 1e-5, -1e-5, 3.552e45, 5e6, 3.786e14, 0.0009, 
                       3.465e-8, 8.155e36, 1e12, 1e-5, 1e-4, 1e-20, 1e-5, 1e37, 1e-3});
    systems.push_back({"Pillars of Creation", 1e21, 2.813e32, 1e-3, -1e-3, 3.552e48, 2e3, 3.156e13, 0.0, 
                       8.457e-14, 1.989e32, 9.46e15, 1e-4, 1e-3, 1e-21, 1e-8, 0.0, 1e-5});
    systems.push_back({"Westerlund 2", 1e22, 1e35, 1e-4, -1e-4, 1e53, 1e4, 3.156e13, 0.0, 
                       1e-12, 1.989e35, 3.086e17, 1e-4, 1e-3, 1e-21, 1e-8, 1e35, 1e-4});
    return systems;
}

// =============================================================================
// Main Function
// =============================================================================
int main() {
    std::cout << "=== Source9 UQFF Unified Field Calculator (CoAnQi) ===" << std::endl;
    std::cout << "Version: 2.0-Enhanced (Module-Aligned)\n" << std::endl;
    
    UQFFModule9 module;
    module.printInfo();
    
    // Load or use default bodies
    std::vector<CelestialBody> bodies = load_bodies_csv("bodies.csv");
    if (bodies.empty()) {
        std::cout << "\nUsing default celestial bodies..." << std::endl;
        bodies = get_default_bodies();
    }
    
    // Compute UQFF for celestial bodies
    std::cout << "\n--- Computing Unified Field Strength for " << bodies.size() << " bodies ---" << std::endl;
    double t = 0.0, tn = 0.0, theta = 0.0;
    std::vector<double> fu_values;
    
    for (const auto& body : bodies) {
        double r = body.Rb;
        double FU = module.compute_FU(body, r, t, tn, theta);
        fu_values.push_back(FU);
        std::cout << "Body: " << body.name << std::endl;
        std::cout << "  r = " << r << " m" << std::endl;
        std::cout << "  FU = " << FU << std::endl;
        std::cout << "  Ug1 = " << module.compute_Ug1(body, r, t, tn) << std::endl;
        std::cout << "  Ug2 = " << module.compute_Ug2(body, r, t, tn) << std::endl;
        std::cout << "  Ug3 = " << module.compute_Ug3(body, r, t, tn, theta) << std::endl;
        std::cout << "  Ug4 = " << module.compute_Ug4(t, tn) << std::endl;
        std::cout << "  Um  = " << module.compute_Um(body, t, tn, body.Rb) << std::endl;
        module.output_json_params(body);
        std::cout << std::endl;
    }
    
    // MUGE Systems
    std::cout << "--- Computing MUGE for astrophysical systems ---" << std::endl;
    std::vector<MUGESystem> muge_systems = get_default_muge_systems();
    ResonanceParams res;
    std::vector<double> compressed_values, resonance_values;
    
    for (const auto& sys : muge_systems) {
        double compressed_g = module.compute_compressed_MUGE(const_cast<MUGESystem&>(sys));
        double resonance_g = module.compute_resonance_MUGE(const_cast<MUGESystem&>(sys), res);
        compressed_values.push_back(compressed_g);
        resonance_values.push_back(resonance_g);
        std::cout << "System: " << sys.name << std::endl;
        std::cout << "  Compressed MUGE g = " << compressed_g << " m/s^2" << std::endl;
        std::cout << "  Resonance MUGE g  = " << resonance_g << " m/s^2" << std::endl;
        std::cout << std::endl;
    }
    
    // Export state
    module.exportState("source9_state.txt");
    
    // Summary statistics
    if (!fu_values.empty()) {
        double min_fu = *std::min_element(fu_values.begin(), fu_values.end());
        double max_fu = *std::max_element(fu_values.begin(), fu_values.end());
        double sum_fu = 0.0;
        for (double v : fu_values) sum_fu += v;
        double mean_fu = sum_fu / fu_values.size();
        std::cout << "FU summary - Min: " << min_fu << ", Max: " << max_fu << ", Mean: " << mean_fu << std::endl;
    }
    
    if (!compressed_values.empty()) {
        double min_g = *std::min_element(compressed_values.begin(), compressed_values.end());
        double max_g = *std::max_element(compressed_values.begin(), compressed_values.end());
        double sum_g = 0.0;
        for (double v : compressed_values) sum_g += v;
        double mean_g = sum_g / compressed_values.size();
        std::cout << "Compressed MUGE summary - Min: " << min_g << ", Max: " << max_g << ", Mean: " << mean_g << std::endl;
    }
    
    // ==================== DUAL PHYSICS VALIDATION ====================
    std::cout << "\n=== Dual Physics Validation (FluidSolver + DualMethodValidator) ===" << std::endl;
    
    using namespace UQFFDualPhysics;
    
    // Initialize FluidSolver for Navier-Stokes simulation
    FluidSolver fluidSolver(32, 0.1, 0.0001);
    std::cout << "FluidSolver initialized (32x32 grid)" << std::endl;
    
    // Run fluid simulation with computed MUGE gravity
    double test_g = compressed_values.empty() ? 9.81 : compressed_values[0];
    fluidSolver.add_jet_force(10.0);
    for (int step = 0; step < 10; ++step) {
        fluidSolver.step(test_g * 1e-10);
    }
    std::cout << "FluidSolver simulation: Max velocity = " << fluidSolver.getMaxVelocity() << " m/s" << std::endl;
    
    // Initialize DualMethodValidator
    DualMethodValidator validator("source9_dual_physics.log");
    
    // Validate first MUGE system using dual physics methods
    if (!muge_systems.empty() && !bodies.empty()) {
        const auto& body0 = bodies[0];
        const auto& muge0 = muge_systems[0];
        
        // Create UQFFDualPhysics structures
        UQFFDualPhysics::CelestialBody uqff_body(body0.name, body0.Ms, body0.Rs, body0.Bs_avg);
        uqff_body.SCm_density = body0.SCm_density;
        uqff_body.QUA = body0.QUA;
        uqff_body.omega_c = body0.omega_c;
        
        UQFFDualPhysics::MUGESystem dual_muge(muge0.name, muge0.M, muge0.r);
        dual_muge.B0 = muge0.B;  // Use B from local struct
        
        auto result = validator.validate(uqff_body, dual_muge, t, theta);
        result.print();
    }
    
    std::cout << "Dual Physics: IMPLEMENTED" << std::endl;
    // ================================================================
    
    std::cout << "\n=== Source9 Complete ===" << std::endl;
    return 0;
}
