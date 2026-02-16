/**
 * ================================================================================================
 * UQFF Dual Physics Framework Header
 * 
 * Description: Shared infrastructure for dual-method physics validation across all UQFF modules.
 *              Includes FluidSolver (Navier-Stokes), DualMethodValidator (UQFF vs MUGE),
 *              and dual physics computation methods.
 * 
 * Features:
 *   - FluidSolver: Jos Stam's Stable Fluids (Navier-Stokes simulation)
 *   - DualMethodValidator: Cross-validation between UQFF and MUGE methods
 *   - DualPhysicsMethods: Simultaneous UQFF + MUGE gravity calculations
 *   - CelestialBody and MUGESystem data structures
 *   - Comprehensive validation with tolerance-based convergence checking
 * 
 * Usage:
 *   #include "uqff_dual_physics.h"
 *   using namespace UQFFDualPhysics;
 * 
 * Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
 * Created: January 22, 2026
 * ================================================================================================
 */

#ifndef UQFF_DUAL_PHYSICS_H
#define UQFF_DUAL_PHYSICS_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <fstream>
#include <sstream>
#include <chrono>
#include <ctime>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace UQFFDualPhysics {

// ===========================================================================================
// PHYSICAL CONSTANTS
// ===========================================================================================

constexpr double G_CONST = 6.6743e-11;      // Gravitational constant m³/kg/s²
constexpr double C_LIGHT = 2.998e8;          // Speed of light m/s
constexpr double HBAR = 1.055e-34;           // Reduced Planck constant J·s
constexpr double K_BOLTZ = 1.381e-23;        // Boltzmann constant J/K
constexpr double M_SUN = 1.989e30;           // Solar mass kg
constexpr double MPC_TO_M = 3.086e22;        // Megaparsec to meters
constexpr double KPC_TO_M = 3.086e19;        // Kiloparsec to meters
constexpr double PC_TO_M = 3.086e16;         // Parsec to meters

// ===========================================================================================
// CELESTIAL BODY DATA STRUCTURE
// ===========================================================================================

struct CelestialBody {
    std::string name = "Unknown";
    double M = 1.989e30;                     // Mass (kg)
    double Rs = 1e6;                         // Characteristic radius (m)
    double B0 = 1e10;                        // Magnetic field (T)
    double SCm_density = 1e15;               // SCm density
    double QUA = 1e-11;                      // Quantum UA factor
    double Pcore = 1e-3;                     // Core period
    double PSCm = 1.0;                       // SCm period
    double omega_c = 1e-6;                   // Angular frequency
    double Z = 0.0;                          // Charge number
    double C_R = 1.0;                        // Reactivity coefficient
    
    CelestialBody() = default;
    CelestialBody(const std::string& n, double mass, double radius, double b0 = 1e10)
        : name(n), M(mass), Rs(radius), B0(b0) {}
};

// ===========================================================================================
// MUGE SYSTEM DATA STRUCTURE
// ===========================================================================================

struct MUGESystem {
    std::string name = "Unknown";
    double M = 1.989e30;                     // Mass (kg)
    double r = 1e6;                          // Distance (m)
    double H0 = 70.0;                        // Hubble constant (km/s/Mpc)
    double Lambda = 1.11e-52;                // Cosmological constant (m⁻²)
    double B0 = 1e10;                        // Magnetic field (T)
    double omega = 1e-6;                     // Angular frequency (rad/s)
    double rho_dm = 1e-21;                   // Dark matter density (kg/m³)
    double v_flow = 1e3;                     // Flow velocity (m/s)
    double T = 1e6;                          // Temperature (K)
    
    MUGESystem() = default;
    MUGESystem(const std::string& n, double mass, double radius)
        : name(n), M(mass), r(radius) {}
};

// ===========================================================================================
// RESONANCE PARAMETERS
// ===========================================================================================

struct ResonanceParams {
    double aDPM = 1e-12;
    double aTHz = 1e12;
    double Avac_diff = 1e-15;
    double aSuperFreq = 1e-8;
    double aAetherRes = 1e-20;
    double Ug4i = 1e-10;
    double aQuantumFreq = 1e15;
    double aAetherFreq = 1e-18;
    double aFluidFreq = 1e6;
    double Osc_term = 1e-5;
    double aExpFreq = 1e-16;
    double fTRZ = 1e9;
    double wormhole_metric = 1e-30;
};

// Default resonance parameters
inline ResonanceParams getDefaultResonanceParams() {
    return ResonanceParams{};
}

// ===========================================================================================
// VALIDATION RESULT STRUCTURE
// ===========================================================================================

struct ValidationResult {
    double uqff_value = 0.0;
    double muge_compressed_value = 0.0;
    double muge_resonance_value = 0.0;
    double uqff_muge_diff = 0.0;
    double compressed_resonance_diff = 0.0;
    bool convergence_achieved = false;
    std::string discrepancy_analysis;
    
    void print() const {
        std::cout << "=== Dual Physics Validation Result ===" << std::endl;
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "  UQFF Value:           " << uqff_value << " m/s²" << std::endl;
        std::cout << "  MUGE Compressed:      " << muge_compressed_value << " m/s²" << std::endl;
        std::cout << "  MUGE Resonance:       " << muge_resonance_value << " m/s²" << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "  UQFF-MUGE Diff:       " << uqff_muge_diff << "%" << std::endl;
        std::cout << "  Compressed-Res Diff:  " << compressed_resonance_diff << "%" << std::endl;
        std::cout << "  Convergence:          " << (convergence_achieved ? "PASS" : "FAIL") << std::endl;
        std::cout << "  Analysis:             " << discrepancy_analysis << std::endl;
    }
};

// ===========================================================================================
// PHYSICS CONSTRAINTS
// ===========================================================================================

struct PhysicsConstraints {
    double min_gravity = 1e-20;
    double max_gravity = 1e20;
    double tolerance = 25.0;               // Percentage tolerance
};

// ===========================================================================================
// FLUIDSOLVER CLASS - JOS STAM'S STABLE FLUIDS
// Navier-Stokes solver for fluid dynamics simulations
// ===========================================================================================

class FluidSolver {
private:
    int N;                                  // Grid size
    double dt;                              // Time step
    double visc;                            // Viscosity
    
    inline int IX(int i, int j) const { return (i) + (N + 2) * (j); }
    
public:
    std::vector<double> u, v, u_prev, v_prev, dens, dens_prev;
    
    FluidSolver(int gridSize = 32, double timeStep = 0.1, double viscosity = 0.0001)
        : N(gridSize), dt(timeStep), visc(viscosity) {
        int size = (N + 2) * (N + 2);
        u.resize(size, 0.0);
        v.resize(size, 0.0);
        u_prev.resize(size, 0.0);
        v_prev.resize(size, 0.0);
        dens.resize(size, 0.0);
        dens_prev.resize(size, 0.0);
    }
    
    void add_source(std::vector<double>& x, std::vector<double>& s) {
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] += dt * s[i];
        }
    }
    
    void diffuse(int b, std::vector<double>& x, std::vector<double>& x0, double diff) {
        double a = dt * diff * N * N;
        for (int k = 0; k < 20; ++k) {
            for (int i = 1; i <= N; ++i) {
                for (int j = 1; j <= N; ++j) {
                    x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                                       x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
                                  (1 + 4 * a);
                }
            }
            set_bnd(b, x);
        }
    }
    
    void advect(int b, std::vector<double>& d, std::vector<double>& d0) {
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= N; ++j) {
                double x = i - dt * N * u[IX(i, j)];
                double y = j - dt * N * v[IX(i, j)];
                if (x < 0.5) x = 0.5;
                if (x > N + 0.5) x = N + 0.5;
                if (y < 0.5) y = 0.5;
                if (y > N + 0.5) y = N + 0.5;
                int i0 = static_cast<int>(x), i1 = i0 + 1;
                int j0 = static_cast<int>(y), j1 = j0 + 1;
                double s1 = x - i0, s0 = 1 - s1;
                double t1 = y - j0, t0 = 1 - t1;
                d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                              s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
            }
        }
        set_bnd(b, d);
    }
    
    void project(std::vector<double>& vel_u, std::vector<double>& vel_v, 
                 std::vector<double>& p, std::vector<double>& div) {
        double h = 1.0 / N;
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= N; ++j) {
                div[IX(i, j)] = -0.5 * h * (vel_u[IX(i + 1, j)] - vel_u[IX(i - 1, j)] + 
                                            vel_v[IX(i, j + 1)] - vel_v[IX(i, j - 1)]);
                p[IX(i, j)] = 0;
            }
        }
        set_bnd(0, div);
        set_bnd(0, p);
        for (int k = 0; k < 20; ++k) {
            for (int i = 1; i <= N; ++i) {
                for (int j = 1; j <= N; ++j) {
                    p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] +
                                   p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4;
                }
            }
            set_bnd(0, p);
        }
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= N; ++j) {
                vel_u[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
                vel_v[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
            }
        }
        set_bnd(1, vel_u);
        set_bnd(2, vel_v);
    }
    
    void set_bnd(int b, std::vector<double>& x) {
        for (int i = 1; i <= N; ++i) {
            x[IX(0, i)] = (b == 1) ? -x[IX(1, i)] : x[IX(1, i)];
            x[IX(N + 1, i)] = (b == 1) ? -x[IX(N, i)] : x[IX(N, i)];
            x[IX(i, 0)] = (b == 2) ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, N + 1)] = (b == 2) ? -x[IX(i, N)] : x[IX(i, N)];
        }
        x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
        x[IX(0, N + 1)] = 0.5 * (x[IX(1, N + 1)] + x[IX(0, N)]);
        x[IX(N + 1, 0)] = 0.5 * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
        x[IX(N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
    }
    
    void step(double uqff_g = 0.0) {
        // Add UQFF gravity as body force in v direction
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= N; ++j) {
                v[IX(i, j)] += dt * uqff_g;
            }
        }
        
        diffuse(1, u_prev, u, visc);
        diffuse(2, v_prev, v, visc);
        project(u_prev, v_prev, u, v);
        advect(1, u, u_prev);
        advect(2, v, v_prev);
        project(u, v, u_prev, v_prev);
    }
    
    void add_jet_force(double force) {
        // Add force in the center as a jet
        for (int i = N / 4; i <= 3 * N / 4; ++i) {
            v[IX(i, N / 2)] += force;
        }
    }
    
    double getMaxVelocity() const {
        double max_vel = 0.0;
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= N; ++j) {
                double vel = std::sqrt(u[IX(i, j)] * u[IX(i, j)] + v[IX(i, j)] * v[IX(i, j)]);
                if (vel > max_vel) max_vel = vel;
            }
        }
        return max_vel;
    }
    
    double getAverageVelocity() const {
        double sum_vel = 0.0;
        int count = 0;
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= N; ++j) {
                sum_vel += std::sqrt(u[IX(i, j)] * u[IX(i, j)] + v[IX(i, j)] * v[IX(i, j)]);
                count++;
            }
        }
        return sum_vel / count;
    }
    
    void print_velocity_field() {
        std::cout << "Velocity field (magnitude):" << std::endl;
        for (int j = N; j >= 1; --j) {
            for (int i = 1; i <= N; ++i) {
                double mag = std::sqrt(u[IX(i, j)] * u[IX(i, j)] + v[IX(i, j)] * v[IX(i, j)]);
                char sym = (mag > 1.0) ? '#' : (mag > 0.5) ? '+' : (mag > 0.1) ? '.' : ' ';
                std::cout << sym;
            }
            std::cout << std::endl;
        }
    }
    
    void reset() {
        std::fill(u.begin(), u.end(), 0.0);
        std::fill(v.begin(), v.end(), 0.0);
        std::fill(u_prev.begin(), u_prev.end(), 0.0);
        std::fill(v_prev.begin(), v_prev.end(), 0.0);
        std::fill(dens.begin(), dens.end(), 0.0);
        std::fill(dens_prev.begin(), dens_prev.end(), 0.0);
    }
};

// ===========================================================================================
// DUAL PHYSICS METHODS - UQFF AND MUGE CALCULATIONS
// ===========================================================================================

class DualPhysicsMethods {
public:
    // UQFF Core Components
    static double compute_Ug1(const CelestialBody& body, double r) {
        // Magnetic dipole contribution
        double mu_0 = 4.0 * M_PI * 1e-7;
        return (mu_0 * body.B0 * body.B0 * body.Rs * body.Rs * body.Rs) / (4.0 * M_PI * r * r * r * r);
    }
    
    static double compute_Ug2(const CelestialBody& body, double r) {
        // Charge-reactivity contribution
        double k_e = 8.99e9;
        double Q = body.Z * 1.602e-19;
        return (k_e * Q * body.C_R) / (r * r);
    }
    
    static double compute_Ug3(const CelestialBody& body, double r, double t) {
        // String rotation contribution
        double omega = 2.0 * M_PI / body.Pcore;
        return body.SCm_density * omega * omega * r * std::sin(omega * t);
    }
    
    static double compute_Ug4(const CelestialBody& body, double r) {
        // Vacuum concentration contribution
        double rho_vac = 7.09e-36;
        return rho_vac * body.QUA * G_CONST * body.M / (r * r);
    }
    
    static double compute_Ubi(const CelestialBody& body, double r) {
        // Buoyancy force (UQFF characteristic)
        double rho_env = 1e-20;
        double rho_body = body.M / (4.0 / 3.0 * M_PI * body.Rs * body.Rs * body.Rs);
        return (rho_env - rho_body) * G_CONST * body.M / (r * r);
    }
    
    static double compute_Um(const CelestialBody& body, double r) {
        // Magnetism contribution
        double mu_0 = 4.0 * M_PI * 1e-7;
        return (mu_0 * body.B0 * body.M) / (4.0 * M_PI * r * r * r);
    }
    
    static double compute_UQFF(const CelestialBody& body, double r, double t, double theta = 0.0) {
        // Complete UQFF unified field
        double Ug1 = compute_Ug1(body, r);
        double Ug2 = compute_Ug2(body, r);
        double Ug3 = compute_Ug3(body, r, t);
        double Ug4 = compute_Ug4(body, r);
        double Ubi = compute_Ubi(body, r);
        double Um = compute_Um(body, r);
        
        // Angular modulation
        double angular_factor = 1.0 + 0.1 * std::cos(theta);
        
        return angular_factor * (Ug1 + Ug2 + Ug3 + Ug4 + Ubi + Um);
    }
    
    // MUGE Core Components
    static double compute_Newtonian(const MUGESystem& sys) {
        // Standard Newtonian gravity
        return G_CONST * sys.M / (sys.r * sys.r);
    }
    
    static double compute_HubbleExpansion(const MUGESystem& sys) {
        // Hubble expansion correction
        double H0_SI = sys.H0 * 1000.0 / MPC_TO_M;
        return -H0_SI * H0_SI * sys.r;
    }
    
    static double compute_MagneticSuppression(const MUGESystem& sys) {
        // Magnetic field pressure gradient contribution to effective gravity
        // P_mag = B²/(2μ₀), then g_B = (P_mag × r²) / M = magnetic pressure effect
        // This converts magnetic energy density to acceleration: [J/m³] × [m³/kg] = [m/s²]
        double mu_0 = 4.0 * M_PI * 1e-7;
        double P_mag = (sys.B0 * sys.B0) / (2.0 * mu_0);  // Magnetic pressure (Pa)
        double r3_over_M = (sys.r * sys.r * sys.r) / sys.M;  // Volume/Mass
        return -P_mag * r3_over_M / sys.r;  // Converts to acceleration, negative = suppression
    }
    
    static double compute_CosmologicalConstant(const MUGESystem& sys) {
        // Cosmological constant Λ term
        return (C_LIGHT * C_LIGHT * sys.Lambda * sys.r) / 3.0;
    }
    
    static double compute_QuantumCorrection(const MUGESystem& sys) {
        // Quantum correction (ℏ-based)
        return (HBAR * sys.omega) / (sys.M * sys.r);
    }
    
    static double compute_FluidDynamics(const MUGESystem& sys) {
        // Navier-Stokes based correction
        return sys.v_flow * sys.v_flow / sys.r;
    }
    
    static double compute_DarkMatterPerturbation(const MUGESystem& sys) {
        // Dark matter density perturbation
        return 4.0 * M_PI * G_CONST * sys.rho_dm * sys.r / 3.0;
    }
    
    static double compute_MUGE_compressed(const MUGESystem& sys) {
        // Complete MUGE compressed formula
        double g_N = compute_Newtonian(sys);
        double g_H = compute_HubbleExpansion(sys);
        double g_B = compute_MagneticSuppression(sys);
        double g_Lambda = compute_CosmologicalConstant(sys);
        double g_Q = compute_QuantumCorrection(sys);
        double g_F = compute_FluidDynamics(sys);
        double g_DM = compute_DarkMatterPerturbation(sys);
        
        return g_N + g_H + g_B + g_Lambda + g_Q + g_F + g_DM;
    }
    
    static double compute_MUGE_resonance(const MUGESystem& sys, const ResonanceParams& params) {
        // MUGE resonance formula with all resonance modes
        double g_base = compute_MUGE_compressed(sys);
        
        // Resonance corrections
        double res_DPM = params.aDPM * std::sin(params.aTHz * 1e-12);
        double res_vac = params.Avac_diff * sys.r;
        double res_super = params.aSuperFreq * sys.omega;
        double res_aether = params.aAetherRes * G_CONST * sys.M / (sys.r * sys.r);
        double res_quantum = params.aQuantumFreq * HBAR / (sys.M * sys.r);
        double res_fluid = params.aFluidFreq * sys.v_flow / sys.r;
        double res_osc = params.Osc_term * std::sin(params.fTRZ * 1e-9);
        double res_wormhole = params.wormhole_metric * C_LIGHT * C_LIGHT / sys.r;
        
        return g_base + res_DPM + res_vac + res_super + res_aether + 
               res_quantum + res_fluid + res_osc + res_wormhole;
    }
};

// ===========================================================================================
// DUAL METHOD VALIDATOR CLASS
// Cross-validation between UQFF (buoyancy-based) and MUGE (Newtonian+corrections)
// ===========================================================================================

class DualMethodValidator {
private:
    std::map<std::string, PhysicsConstraints> constraints;
    std::string log_path;
    bool logging_enabled;
    
    void initialize_constraints() {
        constraints["SGR1745"] = {1e9, 1e13, 10.0};           // Magnetar
        constraints["SagA"] = {1e8, 1e12, 15.0};              // SMBH
        constraints["SagittariusA"] = {1e8, 1e12, 15.0};      // SMBH (alternate name)
        constraints["Tapestry"] = {1e3, 1e7, 20.0};           // Nebula
        constraints["Westerlund2"] = {1e2, 1e6, 20.0};        // Cluster
        constraints["Pillars"] = {1e1, 1e5, 25.0};            // Molecular cloud
        constraints["Rings"] = {1e10, 1e14, 10.0};            // Gravitational lens
        constraints["StudentGuide"] = {1e-20, 1e-16, 50.0};   // Cosmological
        constraints["Sun"] = {200.0, 300.0, 5.0};             // Solar
        constraints["Earth"] = {9.0, 11.0, 2.0};              // Terrestrial
        constraints["Jupiter"] = {20.0, 30.0, 5.0};           // Gas giant
        constraints["NeutronStar"] = {1e11, 1e13, 15.0};      // Neutron star
        constraints["BlackHole"] = {1e12, 1e16, 20.0};        // Black hole
        constraints["Generic"] = {1e-30, 1e30, 50.0};         // Generic fallback
    }
    
public:
    DualMethodValidator(const std::string& logPath = "dual_physics_validation.log")
        : log_path(logPath), logging_enabled(true) {
        initialize_constraints();
    }
    
    void setLogging(bool enabled) { logging_enabled = enabled; }
    
    void addConstraint(const std::string& name, double min_g, double max_g, double tol) {
        constraints[name] = {min_g, max_g, tol};
    }
    
    ValidationResult validate(const CelestialBody& body, const MUGESystem& muge,
                              double t = 0.0, double theta = 0.0) {
        ValidationResult result;
        
        // Calculate UQFF gravity
        result.uqff_value = DualPhysicsMethods::compute_UQFF(body, muge.r, t, theta);
        
        // Calculate MUGE compressed gravity
        result.muge_compressed_value = DualPhysicsMethods::compute_MUGE_compressed(muge);
        
        // Calculate MUGE resonance gravity
        ResonanceParams params = getDefaultResonanceParams();
        result.muge_resonance_value = DualPhysicsMethods::compute_MUGE_resonance(muge, params);
        
        // Compute percentage differences
        if (std::abs(result.uqff_value) > 1e-50) {
            result.uqff_muge_diff = std::abs(
                (result.uqff_value - result.muge_compressed_value) / result.uqff_value * 100.0
            );
        } else {
            result.uqff_muge_diff = 0.0;
        }
        
        if (std::abs(result.muge_compressed_value) > 1e-50) {
            result.compressed_resonance_diff = std::abs(
                (result.muge_compressed_value - result.muge_resonance_value) / 
                result.muge_compressed_value * 100.0
            );
        } else {
            result.compressed_resonance_diff = 0.0;
        }
        
        // Check convergence against constraints
        std::string system_name = muge.name;
        auto it = constraints.find(system_name);
        if (it == constraints.end()) {
            it = constraints.find("Generic");
        }
        
        if (it != constraints.end()) {
            const auto& c = it->second;
            double g_test = std::abs(result.muge_compressed_value);
            bool within_range = (g_test >= c.min_gravity && g_test <= c.max_gravity);
            bool uqff_converged = (result.uqff_muge_diff < c.tolerance);
            bool muge_converged = (result.compressed_resonance_diff < c.tolerance);
            
            result.convergence_achieved = within_range && uqff_converged && muge_converged;
            
            if (!result.convergence_achieved) {
                std::ostringstream oss;
                if (!within_range) {
                    oss << "RANGE_VIOLATION: |g|=" << g_test 
                        << " outside [" << c.min_gravity << ", " << c.max_gravity << "]; ";
                }
                if (!uqff_converged) {
                    oss << "UQFF_MUGE_DIVERGENCE: " << result.uqff_muge_diff 
                        << "% > " << c.tolerance << "%; ";
                }
                if (!muge_converged) {
                    oss << "MUGE_VARIANT_DIVERGENCE: " << result.compressed_resonance_diff 
                        << "% > " << c.tolerance << "%; ";
                }
                result.discrepancy_analysis = oss.str();
            } else {
                result.discrepancy_analysis = "CONVERGED: All methods agree within tolerance";
            }
        } else {
            result.convergence_achieved = true;
            result.discrepancy_analysis = "NO_CONSTRAINTS: Using default tolerance";
        }
        
        // Log if enabled
        if (logging_enabled) {
            log_validation(system_name, result);
        }
        
        return result;
    }
    
    ValidationResult validate_system(const std::string& name, double M, double r,
                                     double B0 = 1e10, double t = 0.0) {
        CelestialBody body(name, M, r * 0.01, B0);  // Rs = r/100 estimate
        MUGESystem muge(name, M, r);
        muge.B0 = B0;
        return validate(body, muge, t, 0.0);
    }
    
    void log_validation(const std::string& name, const ValidationResult& result) {
        std::ofstream log_file(log_path, std::ios::app);
        if (!log_file.is_open()) return;
        
        auto now = std::chrono::system_clock::now();
        auto time = std::chrono::system_clock::to_time_t(now);
        
        log_file << "=== " << name << " Validation [" << std::ctime(&time) << "] ===" << std::endl;
        log_file << std::scientific << std::setprecision(6);
        log_file << "  UQFF:               " << result.uqff_value << " m/s²" << std::endl;
        log_file << "  MUGE Compressed:    " << result.muge_compressed_value << " m/s²" << std::endl;
        log_file << "  MUGE Resonance:     " << result.muge_resonance_value << " m/s²" << std::endl;
        log_file << std::fixed << std::setprecision(2);
        log_file << "  UQFF-MUGE Diff:     " << result.uqff_muge_diff << "%" << std::endl;
        log_file << "  Comp-Res Diff:      " << result.compressed_resonance_diff << "%" << std::endl;
        log_file << "  Convergence:        " << (result.convergence_achieved ? "PASS" : "FAIL") << std::endl;
        log_file << "  Analysis:           " << result.discrepancy_analysis << std::endl;
        log_file << std::endl;
        log_file.close();
    }
    
    void run_full_validation_suite() {
        std::cout << "=== Running Full Dual Physics Validation Suite ===" << std::endl;
        
        // Test standard astrophysical systems
        std::vector<std::tuple<std::string, double, double, double>> systems = {
            {"Sun", M_SUN, 6.96e8, 1e-4},
            {"Earth", 5.972e24, 6.371e6, 5e-5},
            {"Jupiter", 1.898e27, 6.991e7, 4.2e-4},
            {"NeutronStar", 2.8 * M_SUN, 1e4, 1e8},
            {"SGR1745", 1.4 * M_SUN, 1e4, 1e11}
        };
        
        int passed = 0, failed = 0;
        for (const auto& sys : systems) {
            std::string name = std::get<0>(sys);
            double M = std::get<1>(sys);
            double r = std::get<2>(sys) * 10;  // Test at 10x radius
            double B0 = std::get<3>(sys);
            
            auto result = validate_system(name, M, r, B0);
            
            std::cout << "  " << name << ": " 
                      << (result.convergence_achieved ? "PASS" : "FAIL") << std::endl;
            
            if (result.convergence_achieved) passed++;
            else failed++;
        }
        
        std::cout << "=== Validation Complete: " << passed << "/" << (passed + failed) 
                  << " passed ===" << std::endl;
    }
};

} // namespace UQFFDualPhysics

#endif // UQFF_DUAL_PHYSICS_H
