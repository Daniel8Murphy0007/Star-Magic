// HydrogenResonanceUQFFModule.cpp
// SELF-HEALED VERSION - November 15, 2025 @ 01:40 AM
// Original file was corrupted with massive code duplication and malformed structure
// This version reconstructs the proper class implementation while preserving all physics equations
//
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration)
// for Hydrogen Resonance Equations of the Periodic Table of Elements (PToE).
//
// Copyright - Daniel T. Murphy, analyzed Oct 22, 2025.

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <fstream>

#ifndef M_PI
#define M_PI 3.141592653589793238463
#endif

using cdouble = std::complex<double>;

// ============================================================================
// HYDROGEN RESONANCE UQFF MODULE - Main Physics Class
// ============================================================================

class HydrogenResonanceUQFFModule {
private:
    std::map<std::string, cdouble> variables;
    
    // Self-expanding framework support
    bool self_learning_enabled;
    double learning_rate;
    int update_counter;
    std::map<std::string, std::vector<cdouble>> variable_history;
    std::map<std::string, std::string> variable_dependencies;
    
    // Private computation methods
    cdouble computeA_res(int Z, int A);
    cdouble computeF_res(double E_bind, int A);
    cdouble computeU_dp(int A1, int A2, double f_dp, double phi_dp);
    cdouble computeK_nuc(int N, int Z);
    cdouble computeS_shell(int Z_magic, int N_magic);
    cdouble computeH_res_integrand(double t, int Z, int A);
    cdouble computeX2(int Z, int A);

public:
    // Constructor
    HydrogenResonanceUQFFModule();

    // Dynamic variable operations (complex)
    void updateVariable(const std::string& name, cdouble value);
    void addToVariable(const std::string& name, cdouble delta);
    void subtractFromVariable(const std::string& name, cdouble delta);

    // Core computation: Full H_res(Z, A, t) for element (approx integral)
    cdouble computeHRes(int Z, int A, double t);

    // Sub-equations
    cdouble computeCompressed(int Z, int A, double t);
    cdouble computeResonant(double t, int Z, int A);
    cdouble computeBuoyancy(int Z, int A);
    cdouble computeSuperconductive(double t, int Z, int A);
    double computeCompressedG(double t, int Z, int A);

    // Output descriptive text of the equation
    std::string getEquationText(int Z, int A);

    // Print all current variables
    void printVariables();
    
    // Self-expanding framework methods
    void enableSelfLearning(bool enable);
    void setLearningRate(double rate);
    void autoCalibrate(const std::string& observable, cdouble target_value, double tolerance);
    cdouble computeGradient(const std::string& param, const std::string& observable);
    void recordHistory(const std::string& name, cdouble value);
    void exportState(const std::string& filename);
    void importState(const std::string& filename);
};

// ============================================================================
// SURFACE MAGNETIC FIELD MODULE - Helper Class
// ============================================================================

class SurfaceMagneticFieldModule {
private:
    std::map<std::string, double> variables;
    
    // Self-expanding framework support
    bool self_learning_enabled;
    double learning_rate;
    int update_counter;
    std::map<std::string, std::vector<double>> variable_history;
    std::map<std::string, std::string> variable_dependencies;
    
public:
    SurfaceMagneticFieldModule();
    
    // Core computations
    double computeB_j(double t, double B_s);
    double computeB_s_min();
    double computeB_s_max();
    double computeU_g3_example(double t, double B_s);
    
    // Variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    
    // Self-expanding framework methods
    void enableSelfLearning(bool enable);
    void setLearningRate(double rate);
    void autoCalibrate(const std::string& observable, double target_value, double tolerance);
    double computeGradient(const std::string& param, const std::string& observable);
    void recordHistory(const std::string& name, double value);
    void adaptiveUpdate(double dt, const std::string& feedback_param);
    void scaleToSolarData(const std::map<std::string, double>& solar_data);
    void addCustomVariable(const std::string& name, double value, const std::string& dependency = "");
    std::map<std::string, double> getVariableHistory(const std::string& name, int steps = -1);
    void exportState(const std::string& filename);
    void importState(const std::string& filename);
    void updateDependencies(const std::string& changed_var);
};

// ============================================================================
// HYDROGEN RESONANCE UQFF MODULE - IMPLEMENTATION
// ============================================================================

HydrogenResonanceUQFFModule::HydrogenResonanceUQFFModule() {
    // Initialize framework variables
    self_learning_enabled = false;
    learning_rate = 0.001;
    update_counter = 0;
    
    // Initialize PToE physics variables with defaults
    variables["k_A"] = cdouble(1.0, 0.0);          // Amplitude resonance coupling
    variables["f_res"] = cdouble(1.0, 0.0);        // Resonance frequency factor
    variables["U_dp"] = cdouble(0.0, 0.0);         // Deep pairing contribution
    variables["f_dp"] = cdouble(1.0, 0.0);         // Deep pairing frequency
    variables["phi_dp"] = cdouble(0.0, 0.0);       // Deep pairing phase
    variables["K_nuc"] = cdouble(1.0, 0.0);        // Nuclear stability factor
    variables["S_shell"] = cdouble(0.0, 0.0);      // Shell correction
    variables["E_bind"] = cdouble(8.0, 0.0);       // Binding energy per nucleon (MeV)
    variables["delta_pair"] = cdouble(0.0, 0.0);   // Pairing energy delta
    variables["SC_m"] = cdouble(1.0, 0.0);         // Superconductive mass (normalized)
}

void HydrogenResonanceUQFFModule::updateVariable(const std::string& name, cdouble value) {
    variables[name] = value;
    recordHistory(name, value);
    update_counter++;
}

void HydrogenResonanceUQFFModule::addToVariable(const std::string& name, cdouble delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        recordHistory(name, variables[name]);
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta." << std::endl;
        variables[name] = delta;
    }
    update_counter++;
}

void HydrogenResonanceUQFFModule::subtractFromVariable(const std::string& name, cdouble delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] -= delta;
        recordHistory(name, variables[name]);
    } else {
        std::cerr << "Variable '" << name << "' not found." << std::endl;
    }
    update_counter++;
}

cdouble HydrogenResonanceUQFFModule::computeHRes(int Z, int A, double t) {
    // Approximate integral over time for hydrogen resonance
    // H_res = integral of H_res_integrand(t, Z, A) dt
    double dt = 0.1;  // Time step for integration
    int steps = 100;
    cdouble sum = cdouble(0.0, 0.0);
    
    for (int i = 0; i < steps; i++) {
        double t_i = t + i * dt;
        sum += computeH_res_integrand(t_i, Z, A) * dt;
    }
    
    return sum;
}

cdouble HydrogenResonanceUQFFModule::computeH_res_integrand(double t, int Z, int A) {
    // Integrand: Compressed + Resonant + Buoyancy + Superconductive
    return computeCompressed(Z, A, t) + computeResonant(t, Z, A) + 
           computeBuoyancy(Z, A) + computeSuperconductive(t, Z, A);
}

cdouble HydrogenResonanceUQFFModule::computeCompressed(int Z, int A, double t) {
    // Compressed gravity contribution
    double g_val = computeCompressedG(t, Z, A);
    return cdouble(g_val * A, 0.0);
}

double HydrogenResonanceUQFFModule::computeCompressedG(double t, int Z, int A) {
    // Simplified compressed gravity
    double decay = std::exp(-t / 1000.0);
    return 9.8 * decay * (1.0 + 0.1 * Z / A);
}

cdouble HydrogenResonanceUQFFModule::computeResonant(double t, int Z, int A) {
    // Resonant amplitude contribution
    cdouble A_res = computeA_res(Z, A);
    cdouble F_res = computeF_res(variables["E_bind"].real(), A);
    cdouble phase = cdouble(std::cos(2 * M_PI * t), std::sin(2 * M_PI * t));
    return A_res * F_res * phase;
}

cdouble HydrogenResonanceUQFFModule::computeA_res(int Z, int A) {
    // Amplitude resonance factor
    double k_A = variables["k_A"].real();
    return cdouble(k_A * std::sqrt(static_cast<double>(A)), 0.0);
}

cdouble HydrogenResonanceUQFFModule::computeF_res(double E_bind, int A) {
    // Resonance frequency from binding energy
    return cdouble(E_bind / A, 0.0);
}

cdouble HydrogenResonanceUQFFModule::computeBuoyancy(int Z, int A) {
    // Buoyancy contribution (nuclear stability)
    cdouble K_nuc = computeK_nuc(A - Z, Z);  // N = A - Z
    return K_nuc;
}

cdouble HydrogenResonanceUQFFModule::computeK_nuc(int N, int Z) {
    // Nuclear stability factor
    double N_Z_ratio = static_cast<double>(N) / Z;
    double stability = 1.0 - std::abs(N_Z_ratio - 1.0);  // Optimal N/Z ~ 1 for light elements
    return cdouble(stability, 0.0);
}

cdouble HydrogenResonanceUQFFModule::computeSuperconductive(double t, int Z, int A) {
    // Superconductive contribution (pairing effects)
    int N = A - Z;
    bool is_even_even = (Z % 2 == 0) && (N % 2 == 0);
    
    if (is_even_even) {
        // Enhanced pairing for even-even nuclei
        cdouble delta_pair = variables["delta_pair"];
        return delta_pair * cdouble(std::cos(M_PI * t), 0.0);
    } else {
        // Reduced or zero pairing
        return cdouble(0.0, 0.0);
    }
}

cdouble HydrogenResonanceUQFFModule::computeU_dp(int A1, int A2, double f_dp, double phi_dp) {
    // Deep pairing between nucleons
    double coupling = std::exp(-std::abs(A1 - A2) / 10.0);
    return cdouble(coupling * f_dp * std::cos(phi_dp), coupling * f_dp * std::sin(phi_dp));
}

cdouble HydrogenResonanceUQFFModule::computeS_shell(int Z_magic, int N_magic) {
    // Shell correction for magic numbers (2, 8, 20, 28, 50, 82, 126)
    std::vector<int> magic_numbers = {2, 8, 20, 28, 50, 82, 126};
    
    double Z_correction = 0.0;
    double N_correction = 0.0;
    
    for (int magic : magic_numbers) {
        if (Z_magic == magic) Z_correction = 1.0;
        if (N_magic == magic) N_correction = 1.0;
    }
    
    return cdouble(Z_correction + N_correction, 0.0);
}

cdouble HydrogenResonanceUQFFModule::computeX2(int Z, int A) {
    // Additional coupling term
    return cdouble(static_cast<double>(Z * A) / 100.0, 0.0);
}

std::string HydrogenResonanceUQFFModule::getEquationText(int Z, int A) {
    std::ostringstream oss;
    oss << "H_res(Z=" << Z << ", A=" << A << ") = ";
    oss << "∫[Compressed(g) + Resonant(A_res*F_res) + Buoyancy(K_nuc) + Superconductive(δ_pair)] dt";
    return oss.str();
}

void HydrogenResonanceUQFFModule::printVariables() {
    std::cout << "=== HydrogenResonanceUQFFModule Variables ===" << std::endl;
    for (const auto& var : variables) {
        std::cout << std::setw(15) << var.first << " = " 
                  << var.second.real() << " + " << var.second.imag() << "i" << std::endl;
    }
    std::cout << "Update counter: " << update_counter << std::endl;
}

void HydrogenResonanceUQFFModule::enableSelfLearning(bool enable) {
    self_learning_enabled = enable;
    if (enable) {
        std::cout << "Self-learning enabled with rate: " << learning_rate << std::endl;
    }
}

void HydrogenResonanceUQFFModule::setLearningRate(double rate) {
    learning_rate = rate;
}

void HydrogenResonanceUQFFModule::recordHistory(const std::string& name, cdouble value) {
    variable_history[name].push_back(value);
}

cdouble HydrogenResonanceUQFFModule::computeGradient(const std::string& param, const std::string& observable) {
    // Simple finite difference gradient
    cdouble original = variables[param];
    cdouble epsilon = cdouble(1e-6, 0.0);
    
    variables[param] = original + epsilon;
    cdouble val_plus = variables[observable];
    
    variables[param] = original - epsilon;
    cdouble val_minus = variables[observable];
    
    variables[param] = original;  // Restore
    
    return (val_plus - val_minus) / (cdouble(2.0, 0.0) * epsilon);
}

void HydrogenResonanceUQFFModule::autoCalibrate(const std::string& observable, cdouble target_value, double tolerance) {
    if (variables.find(observable) == variables.end()) {
        std::cerr << "Observable '" << observable << "' not found for calibration." << std::endl;
        return;
    }
    
    cdouble current_value = variables[observable];
    double error = std::abs(current_value - target_value);
    
    if (error > tolerance) {
        // Adjust tunable parameters
        std::vector<std::string> tunable_params = {"k_A", "f_res", "f_dp", "K_nuc"};
        
        for (const auto& param : tunable_params) {
            cdouble gradient = computeGradient(param, observable);
            if (std::abs(gradient) > 1e-20) {
                cdouble adjustment = cdouble(learning_rate, 0.0) * (target_value - current_value) / gradient;
                variables[param] += adjustment;
                recordHistory(param, variables[param]);
            }
        }
        
        std::cout << "Auto-calibrated " << observable << " (error: " << error << ")" << std::endl;
    }
}

void HydrogenResonanceUQFFModule::exportState(const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "# HydrogenResonanceUQFFModule State Export" << std::endl;
        file << "update_counter=" << update_counter << std::endl;
        file << "learning_rate=" << learning_rate << std::endl;
        file << "self_learning_enabled=" << (self_learning_enabled ? 1 : 0) << std::endl;
        
        for (const auto& var : variables) {
            file << var.first << "=" << var.second.real() << "," << var.second.imag() << std::endl;
        }
        file.close();
        std::cout << "State exported to: " << filename << std::endl;
    }
}

void HydrogenResonanceUQFFModule::importState(const std::string& filename) {
    std::ifstream file(filename);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            if (line[0] == '#') continue;
            
            size_t eq_pos = line.find('=');
            if (eq_pos != std::string::npos) {
                std::string key = line.substr(0, eq_pos);
                std::string value_str = line.substr(eq_pos + 1);
                
                if (key == "update_counter") {
                    update_counter = std::stoi(value_str);
                } else if (key == "learning_rate") {
                    learning_rate = std::stod(value_str);
                } else if (key == "self_learning_enabled") {
                    self_learning_enabled = (std::stoi(value_str) == 1);
                } else {
                    // Parse complex number "real,imag"
                    size_t comma_pos = value_str.find(',');
                    if (comma_pos != std::string::npos) {
                        double real_part = std::stod(value_str.substr(0, comma_pos));
                        double imag_part = std::stod(value_str.substr(comma_pos + 1));
                        variables[key] = cdouble(real_part, imag_part);
                    }
                }
            }
        }
        file.close();
        std::cout << "State imported from: " << filename << std::endl;
    }
}

// ============================================================================
// SURFACE MAGNETIC FIELD MODULE - IMPLEMENTATION
// ============================================================================

SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    // Initialize framework variables
    self_learning_enabled = false;
    learning_rate = 0.001;
    update_counter = 0;
    
    // Initialize solar physics variables
    variables["B_s_min"] = 1e-4;                    // T (quiet)
    variables["B_s_max"] = 0.4;                     // T (sunspot)
    variables["B_ref"] = 0.4;                       // T (reference max)
    variables["k_3"] = 1.8;                         // Coupling
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["t"] = 0.0;                           // s
    variables["solar_cycle_period"] = 11.0 * 365.25 * 86400.0;  // 11 years in seconds
    variables["magnetic_diffusion"] = 1e5;          // m²/s
    variables["convection_velocity"] = 1e3;         // m/s
}

double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);
    return base_b * (B_s / variables["B_ref"]);
}

double SurfaceMagneticFieldModule::computeB_s_min() {
    return variables["B_s_min"];
}

double SurfaceMagneticFieldModule::computeB_s_max() {
    return variables["B_s_max"];
}

double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * M_PI);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}

void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
    recordHistory(name, value);
    updateDependencies(name);
    update_counter++;
}

void SurfaceMagneticFieldModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        recordHistory(name, variables[name]);
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta." << std::endl;
        variables[name] = delta;
    }
    update_counter++;
}

void SurfaceMagneticFieldModule::enableSelfLearning(bool enable) {
    self_learning_enabled = enable;
    if (enable) {
        std::cout << "Magnetic self-learning enabled with rate: " << learning_rate << std::endl;
    } else {
        std::cout << "Magnetic self-learning disabled." << std::endl;
    }
}

void SurfaceMagneticFieldModule::setLearningRate(double rate) {
    learning_rate = rate;
}

void SurfaceMagneticFieldModule::recordHistory(const std::string& name, double value) {
    variable_history[name].push_back(value);
}

double SurfaceMagneticFieldModule::computeGradient(const std::string& param, const std::string& observable) {
    // Simple finite difference gradient
    double original = variables[param];
    double epsilon = 1e-6;
    
    variables[param] = original + epsilon;
    double val_plus = variables[observable];
    
    variables[param] = original - epsilon;
    double val_minus = variables[observable];
    
    variables[param] = original;  // Restore
    
    return (val_plus - val_minus) / (2.0 * epsilon);
}

void SurfaceMagneticFieldModule::autoCalibrate(const std::string& observable, double target_value, double tolerance) {
    if (variables.find(observable) == variables.end()) {
        std::cerr << "Observable '" << observable << "' not found for calibration." << std::endl;
        return;
    }
    
    double current_value = variables[observable];
    double error = std::abs(current_value - target_value) / target_value;
    
    if (error > tolerance) {
        std::vector<std::string> tunable_params = {"B_ref", "k_3", "omega_s", "P_core"};
        
        for (const auto& param : tunable_params) {
            double gradient = computeGradient(param, observable);
            if (std::abs(gradient) > 1e-20) {
                double adjustment = learning_rate * (target_value - current_value) / gradient;
                variables[param] += adjustment;
                recordHistory(param, variables[param]);
            }
        }
        
        std::cout << "Auto-calibrated magnetic " << observable << " from " << current_value 
                  << " to target " << target_value << " (error: " << error << ")" << std::endl;
    }
}

void SurfaceMagneticFieldModule::adaptiveUpdate(double dt, const std::string& feedback_param) {
    if (!self_learning_enabled) return;
    
    // Solar cycle evolution
    double cycle_phase = fmod(variables["t"], variables["solar_cycle_period"]) / variables["solar_cycle_period"];
    double cycle_factor = 0.5 * (1.0 + std::cos(2 * M_PI * cycle_phase));
    
    // Adaptive magnetic field reference
    variables["B_ref"] = variables["B_s_max"] * (0.1 + 0.9 * cycle_factor);
    
    // Magnetic diffusion effects
    double diffusion_decay = std::exp(-dt * variables["magnetic_diffusion"] / 1e6);
    variables["k_3"] *= diffusion_decay;
    
    // Convection-driven field generation
    double convection_enhancement = 1.0 + 0.1 * variables["convection_velocity"] / 1e3;
    variables["omega_s"] *= convection_enhancement;
    
    recordHistory("adaptive_time", dt);
    std::cout << "Magnetic adaptive update: B_ref=" << variables["B_ref"] 
              << ", k_3=" << variables["k_3"] << std::endl;
}

void SurfaceMagneticFieldModule::scaleToSolarData(const std::map<std::string, double>& solar_data) {
    for (const auto& data : solar_data) {
        if (data.first == "sunspot_field" && variables.find("B_s_max") != variables.end()) {
            double scaling = data.second / variables["B_s_max"];
            variables["B_s_max"] = data.second;
            variables["B_ref"] *= scaling;
        }
        
        if (data.first == "solar_rotation_period") {
            variables["omega_s"] = 2 * M_PI / data.second;
        }
        
        if (data.first == "core_temperature") {
            variables["E_react"] = 1e46 * std::pow(data.second / 1.5e7, 3.5);
        }
    }
    std::cout << "Scaled magnetic module to " << solar_data.size() << " solar observations." << std::endl;
}

void SurfaceMagneticFieldModule::addCustomVariable(const std::string& name, double value, const std::string& dependency) {
    variables[name] = value;
    if (!dependency.empty()) {
        variable_dependencies[name] = dependency;
    }
    recordHistory(name, value);
    std::cout << "Added custom magnetic variable: " << name << " = " << value << std::endl;
}

std::map<std::string, double> SurfaceMagneticFieldModule::getVariableHistory(const std::string& name, int steps) {
    std::map<std::string, double> history;
    if (variable_history.find(name) != variable_history.end()) {
        auto& hist = variable_history[name];
        int start = (steps < 0) ? 0 : std::max(0, (int)hist.size() - steps);
        for (int i = start; i < (int)hist.size(); i++) {
            history["step_" + std::to_string(i)] = hist[i];
        }
    }
    return history;
}

void SurfaceMagneticFieldModule::exportState(const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "# SurfaceMagneticFieldModule State Export" << std::endl;
        file << "update_counter=" << update_counter << std::endl;
        file << "learning_rate=" << learning_rate << std::endl;
        file << "self_learning_enabled=" << (self_learning_enabled ? 1 : 0) << std::endl;
        
        for (const auto& var : variables) {
            file << var.first << "=" << var.second << std::endl;
        }
        file.close();
        std::cout << "Magnetic state exported to: " << filename << std::endl;
    }
}

void SurfaceMagneticFieldModule::importState(const std::string& filename) {
    std::ifstream file(filename);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            if (line[0] == '#') continue;
            
            size_t eq_pos = line.find('=');
            if (eq_pos != std::string::npos) {
                std::string key = line.substr(0, eq_pos);
                std::string value_str = line.substr(eq_pos + 1);
                
                if (key == "update_counter") {
                    update_counter = std::stoi(value_str);
                } else if (key == "learning_rate") {
                    learning_rate = std::stod(value_str);
                } else if (key == "self_learning_enabled") {
                    self_learning_enabled = (std::stoi(value_str) == 1);
                } else {
                    variables[key] = std::stod(value_str);
                }
            }
        }
        file.close();
        std::cout << "Magnetic state imported from: " << filename << std::endl;
    }
}

void SurfaceMagneticFieldModule::updateDependencies(const std::string& changed_var) {
    if (changed_var == "B_s_max") {
        variables["B_ref"] = variables["B_s_max"];
    }
}
