/*
 * SIMULTANEOUS 7-LAYER SOLVER v1.5 - MINIMAL BUOYANCY TEST
 *
 * Strategy: Take working v1.0 solver (which plateaued) and add ONLY Ubi term to Layer 7
 * Purpose: Verify that buoyancy term actually improves convergence
 * Approach: Minimal change - only modify Layer 7 equation
 * 
 * Layer 7 OLD (v1.0): E_pair = 2*E_single + E_DPM + E_coulomb
 * Layer 7 NEW (v1.5): E_pair = 2*E_single + E_DPM + E_coulomb - Ubi_term
 *
 * If this fixes convergence, then we know:
 * 1. Buoyancy IS the missing physics
 * 2. E_DPM doesn't need to be dynamic
 * 3. Layer 6-7 over-constraint is correct design
 *
 * Date: May 24, 2026
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <iomanip>
#include <cassert>
#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace UQFF {
namespace Solver {
namespace V1_5 {

// ============================================================================
// PHYSICAL CONSTANTS
// ============================================================================

constexpr double M_E = 9.1093837e-31;
constexpr double HBAR = 1.054571817e-34;
constexpr double C = 299792458.0;
constexpr double G = 6.67430e-11;
constexpr double A0 = 0.529e-10;
constexpr double E_RY = 13.6;
constexpr double ALPHA_FINE = 1.0 / 137.036;
constexpr double RHO_VAC_SCM = 7.09e-37;
constexpr double PAIR_COST_EV = 1.022e6;

// ============================================================================
// BUOYANCY PARAMETER (MINIMAL)
// ============================================================================

constexpr double BETA_I = 0.603;                 // Universal coupling constant
constexpr double UBI_COUPLING = 1e-12;           // Ubi strength relative to energy

// ============================================================================
// LAYER STATE
// ============================================================================

struct LayerState {
    double r_s;
    double g_quantum;
    double v_orb;
    double E_single;
    double psi_norm;
    double E_DPM;
    double E_pair;
    double E_neutrino;
    double Ubi;        // NEW: Track buoyancy force for diagnostics
};

struct ProblemParameters {
    int Z;
    int n;
    double M_nucleus;
    double separation_m;
    double time_s;
    double E_neutrino_in;
};

// ============================================================================
// MINIMAL BUOYANCY CALCULATION
// ============================================================================

/*
 * Calculate minimal Ubi based on energy scales
 * Simple model: Ubi ∝ |E_single| to provide restoring force
 */
inline double calculate_Ubi_minimal(const LayerState& state, int Z, double t_norm = 0.5)
{
    // Ubi is proportional to single-particle energy (which scales with Z²)
    double energy_scale = std::abs(state.E_single);
    
    // Add Z-dependence to buoyancy strength
    double Z_factor = std::pow(Z, 1.0);  // Linear in Z
    
    // Time oscillation
    double oscillation = std::cos(M_PI * t_norm);
    
    // Result: small positive force that balances some of the energy imbalance
    double Ubi = BETA_I * UBI_COUPLING * energy_scale * Z_factor * oscillation;
    
    return Ubi;
}

// ============================================================================
// RESIDUAL CALCULATION (v1.0 with Ubi term added to Layer 7)
// ============================================================================

class Simultaneous7LayerSolverV1_5 {
public:
    Simultaneous7LayerSolverV1_5() = default;
    
    std::vector<double> calculate_residuals(
        LayerState& state,
        const ProblemParameters& params)
    {
        std::vector<double> R(7, 0.0);
        double t_norm = std::fmod(params.time_s, 1.0);
        
        // Layer 1: Buoyancy crossing (unchanged from v1)
        double r_s_target = 2.0 * A0 * ALPHA_FINE * params.Z / (params.n * params.n);
        R[0] = state.r_s - r_s_target;
        
        // Layer 2: Quantum gravity (unchanged)
        double g_classical = G * params.M_nucleus / (state.r_s * state.r_s);
        double quantum_correction = (ALPHA_FINE * params.Z / params.n) * 
                                    ((HBAR * HBAR) / (M_E * C * C * state.r_s * state.r_s));
        double g_quantum_target = g_classical * (1.0 - quantum_correction);
        R[1] = state.g_quantum - g_quantum_target;
        
        // Layer 3: Orbital mechanics (unchanged)
        double v_orb_target = C * ALPHA_FINE * params.Z / params.n;
        R[2] = state.v_orb - v_orb_target;
        
        // Layer 4: Single-particle energy (unchanged)
        double E_single_target = -E_RY * (params.Z * params.Z) / (params.n * params.n);
        double fine_structure = (ALPHA_FINE * ALPHA_FINE * params.Z * params.Z * params.Z * params.Z) / 
                               (params.n * params.n * params.n);
        E_single_target -= fine_structure;
        R[3] = state.E_single - E_single_target;
        
        // Layer 5: Superposition (unchanged)
        double d_spooky = 2.0 * state.r_s;
        double l_coherence = HBAR / (M_E * C);
        double overlap = std::exp(-d_spooky / l_coherence);
        double psi_norm_target = 1.0 / std::sqrt(1.0 + overlap);
        R[4] = state.psi_norm - psi_norm_target;
        
        // Layer 6: E_DPM (unchanged - STILL FIXED)
        double E_DPM_target = PAIR_COST_EV;
        R[5] = state.E_DPM - E_DPM_target;
        
        // Layer 7: NOW WITH UBI TERM ADDED (v1.5 change)
        // Calculate Ubi for this state
        state.Ubi = calculate_Ubi_minimal(state, params.Z, t_norm);
        
        // E_pair target: include Ubi term as energy correction
        // Physical meaning: Ubi provides damping/stabilization of pair energy
        double E_pair_target = 2.0 * state.E_single + 
                              state.E_DPM + 
                              state.E_neutrino;
        
        // Electron-electron repulsion
        double E_coulomb = (8.99e9 * 1.602e-19 * 1.602e-19) / (d_spooky * 1.602e-19);
        E_pair_target += E_coulomb;
        
        // V1.5 CHANGE: Subtract Ubi term from expected E_pair
        // This allows the solver to find E_pair that balances against buoyancy
        E_pair_target -= state.Ubi;  // NEW: buoyancy correction
        
        R[6] = state.E_pair - E_pair_target;
        
        return R;
    }
    
    std::vector<std::vector<double>> calculate_jacobian(
        const LayerState& state,
        const ProblemParameters& params)
    {
        std::vector<std::vector<double>> J(7, std::vector<double>(7, 0.0));
        
        std::vector<double> X_vals = {state.r_s, state.g_quantum, state.v_orb,
                                      state.E_single, state.psi_norm, 
                                      state.E_DPM, state.E_pair};
        
        std::vector<double> dX(7);
        for (int j = 0; j < 7; ++j) {
            double abs_val = std::abs(X_vals[j]);
            dX[j] = (abs_val > 1.0) ? std::sqrt(1e-8) * abs_val : std::sqrt(1e-8);
            dX[j] = std::max(dX[j], 1e-10);
        }
        
        LayerState state_copy = state;
        auto R0 = calculate_residuals(state_copy, params);
        
        for (int j = 0; j < 7; ++j) {
            LayerState state_pert = state;
            
            if (j == 0) state_pert.r_s += dX[j];
            else if (j == 1) state_pert.g_quantum += dX[j];
            else if (j == 2) state_pert.v_orb += dX[j];
            else if (j == 3) state_pert.E_single += dX[j];
            else if (j == 4) state_pert.psi_norm += dX[j];
            else if (j == 5) state_pert.E_DPM += dX[j];
            else if (j == 6) state_pert.E_pair += dX[j];
            
            auto R_pert = calculate_residuals(state_pert, params);
            
            for (int i = 0; i < 7; ++i) {
                double dR = R_pert[i] - R0[i];
                J[i][j] = dR / dX[j];
                
                if (!std::isfinite(J[i][j])) {
                    J[i][j] = 0.0;
                }
            }
        }
        
        return J;
    }
    
    std::vector<double> solve_linear_system(
        std::vector<std::vector<double>>& J,
        const std::vector<double>& R)
    {
        int n = 7;
        std::vector<std::vector<double>> A = J;
        std::vector<double> b = R;
        
        for (int col = 0; col < n; ++col) {
            int pivot_row = col;
            double max_val = std::abs(A[col][col]);
            
            for (int row = col + 1; row < n; ++row) {
                if (std::abs(A[row][col]) > max_val) {
                    max_val = std::abs(A[row][col]);
                    pivot_row = row;
                }
            }
            
            if (max_val < 1e-14) {
                std::vector<double> nan_vec(n, std::numeric_limits<double>::quiet_NaN());
                return nan_vec;
            }
            
            std::swap(A[col], A[pivot_row]);
            std::swap(b[col], b[pivot_row]);
            
            for (int row = col + 1; row < n; ++row) {
                double factor = A[row][col] / A[col][col];
                for (int j = col; j < n; ++j) {
                    A[row][j] -= factor * A[col][j];
                }
                b[row] -= factor * b[col];
            }
        }
        
        std::vector<double> x(n);
        for (int i = n - 1; i >= 0; --i) {
            double sum = b[i];
            for (int j = i + 1; j < n; ++j) {
                sum -= A[i][j] * x[j];
            }
            x[i] = sum / A[i][i];
        }
        
        return x;
    }
    
    double perform_line_search(
        const LayerState& state,
        const std::vector<double>& delta_X,
        const ProblemParameters& params)
    {
        LayerState state_current = state;
        auto R_current = calculate_residuals(state_current, params);
        double norm_current = 0.0;
        for (double r : R_current) norm_current += r * r;
        norm_current = std::sqrt(norm_current);
        
        double alpha = 1.0;
        double threshold = 0.99;
        
        for (int iter = 0; iter < 15; ++iter) {
            LayerState state_candidate = state;
            
            state_candidate.r_s += alpha * delta_X[0];
            state_candidate.g_quantum += alpha * delta_X[1];
            state_candidate.v_orb += alpha * delta_X[2];
            state_candidate.E_single += alpha * delta_X[3];
            state_candidate.psi_norm += alpha * delta_X[4];
            state_candidate.E_DPM += alpha * delta_X[5];
            state_candidate.E_pair += alpha * delta_X[6];
            
            state_candidate.r_s = std::max(1e-15, std::min(1e-5, state_candidate.r_s));
            state_candidate.g_quantum = std::max(0.0, std::min(1e20, state_candidate.g_quantum));
            state_candidate.v_orb = std::max(0.0, std::min(0.99*C, state_candidate.v_orb));
            
            auto R_candidate = calculate_residuals(state_candidate, params);
            double norm_candidate = 0.0;
            for (double r : R_candidate) norm_candidate += r * r;
            norm_candidate = std::sqrt(norm_candidate);
            
            if (norm_candidate < threshold * norm_current || alpha <= 1e-8) {
                return alpha;
            }
            
            alpha *= 0.5;
        }
        
        return 1e-8;
    }
    
    bool solve(LayerState& state, const ProblemParameters& params, int max_iterations = 100)
    {
        std::cout << "\n=== SIMULTANEOUS 7-LAYER SOLVER v1.5 (MINIMAL BUOYANCY TEST) ===\n";
        std::cout << "Z=" << params.Z << ", n=" << params.n << "\n";
        std::cout << "Change: Layer 7 E_pair now includes Ubi damping term\n";
        std::cout << std::string(70, '-') << "\n";
        
        for (int iteration = 0; iteration < max_iterations; ++iteration) {
            auto R = calculate_residuals(state, params);
            
            double norm_R = 0.0;
            for (double r : R) norm_R += r * r;
            norm_R = std::sqrt(norm_R);
            
            std::cout << "Iter " << std::setw(3) << iteration << ": ||R|| = " 
                     << std::scientific << std::setprecision(6) << norm_R << "\n";
            
            if (norm_R < 1e-12) {
                std::cout << "FULL CONVERGENCE at iteration " << iteration << "\n";
                std::cout << "Ubi = " << std::scientific << state.Ubi << " (buoyancy active)\n";
                return true;
            }
            
            auto J = calculate_jacobian(state, params);
            std::vector<double> neg_R(7);
            for (int i = 0; i < 7; ++i) neg_R[i] = -R[i];
            
            auto delta_X = solve_linear_system(J, neg_R);
            
            if (!std::isfinite(delta_X[0])) {
                std::cout << "Singular Jacobian\n";
                return false;
            }
            
            double alpha = perform_line_search(state, delta_X, params);
            
            state.r_s += alpha * delta_X[0];
            state.g_quantum += alpha * delta_X[1];
            state.v_orb += alpha * delta_X[2];
            state.E_single += alpha * delta_X[3];
            state.psi_norm += alpha * delta_X[4];
            state.E_DPM += alpha * delta_X[5];
            state.E_pair += alpha * delta_X[6];
        }
        
        std::cout << "Max iterations reached\n";
        std::cout << "Ubi = " << std::scientific << state.Ubi << " (buoyancy term)\n";
        return false;
    }
};

}}} // namespace

// ============================================================================
// MAIN: TEST v1.5
// ============================================================================

int main()
{
    using namespace UQFF::Solver::V1_5;
    
    std::cout << "\n" << std::string(70, '#') << "\n";
    std::cout << "SOLVER v1.5 TEST: Does Ubi term fix H/He/Ne convergence?\n";
    std::cout << std::string(70, '#') << "\n";
    
    std::vector<int> Z_vals = {1, 2, 10, 54};
    std::vector<std::string> names = {"Hydrogen", "Helium", "Neon", "Xenon"};
    
    for (size_t idx = 0; idx < Z_vals.size(); ++idx) {
        int Z = Z_vals[idx];
        std::string name = names[idx];
        
        std::cout << "\n### " << name << " (Z=" << Z << ") ###\n";
        
        ProblemParameters params;
        params.Z = Z;
        params.n = 1;
        params.M_nucleus = 1.673e-27;
        params.separation_m = 2.0 * A0;
        params.time_s = 0.0;
        params.E_neutrino_in = 0.511;
        
        LayerState state;
        state.r_s = A0 / Z;
        state.g_quantum = 1e13;
        state.v_orb = C * ALPHA_FINE * Z;
        state.E_single = -E_RY * Z * Z;
        state.psi_norm = 0.7;
        state.E_DPM = PAIR_COST_EV * 0.9;
        state.E_pair = 2.0 * state.E_single + state.E_DPM + 1.0;
        state.E_neutrino = 0.511;
        state.Ubi = 0.0;
        
        Simultaneous7LayerSolverV1_5 solver;
        bool converged = solver.solve(state, params);
        
        std::cout << "Result: " << (converged ? "CONVERGED ✓" : "PLATEAU/PARTIAL");
        std::cout << " | Final ||R|| analysis depends on Ubi term\n";
    }
    
    std::cout << "\n" << std::string(70, '#') << "\n";
    std::cout << "HYPOTHESIS TEST:\n";
    std::cout << "If v1.5 converges where v1.0 plateaued:\n";
    std::cout << "  → Buoyancy IS the missing physics\n";
    std::cout << "  → E_DPM doesn't need to be dynamic (v1 design correct)\n";
    std::cout << "If v1.5 still plateaus:\n";
    std::cout << "  → Need complete architectural redesign (v2.0 path)\n";
    std::cout << std::string(70, '#') << "\n\n";
    
    return 0;
}
