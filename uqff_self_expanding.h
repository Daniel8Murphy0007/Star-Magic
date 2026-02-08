/**
 * ================================================================================================
 * UQFF Self-Expanding Framework Header
 * 
 * Description: Shared infrastructure for self-updating, self-expanding, and self-simulating
 *              physics modules. All UQFF modules should include this header.
 * 
 * Features:
 *   - PhysicsTerm base class for dynamic term registration
 *   - Built-in derived terms (DynamicVacuumTerm, QuantumCouplingTerm, DarkMatterHaloTerm, etc.)
 *   - Auto-expanding parameter maps
 *   - Self-simulation with time evolution
 *   - Auto-optimization with learning rate
 *   - Metadata tracking and state persistence
 * 
 * Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
 * Enhanced: January 2026 - Full self-expanding capabilities
 * ================================================================================================
 */

#ifndef UQFF_SELF_EXPANDING_H
#define UQFF_SELF_EXPANDING_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace UQFFExpanding {

// ===========================================================================================
// PhysicsTerm: Abstract base class for all dynamic physics terms
// ===========================================================================================

class PhysicsTerm {
public:
    virtual ~PhysicsTerm() = default;
    
    // Core computation - MUST be implemented by derived classes
    virtual double compute(double t, const std::map<std::string, double>& params) const = 0;
    
    // Identification
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
    
    // Validation - returns true if term can be computed with given params
    virtual bool validate(const std::map<std::string, double>& params) const {
        (void)params;  // Suppress unused warning
        return true;
    }
    
    // Self-optimization - adjust internal parameters based on feedback
    virtual void optimize(double learningRate, double error) {
        (void)learningRate;
        (void)error;
    }
    
    // State persistence
    virtual std::string serializeState() const { return ""; }
    virtual void deserializeState(const std::string& state) { (void)state; }
};

// ===========================================================================================
// Built-in Dynamic Physics Terms
// ===========================================================================================

class DynamicVacuumTerm : public PhysicsTerm {
private:
    double amplitude;
    double frequency;
    double phase;
    
public:
    DynamicVacuumTerm(double amp = 1e-10, double freq = 1e-15, double ph = 0.0)
        : amplitude(amp), frequency(freq), phase(ph) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        double rho_vac = params.count("rho_vac_UA") ? params.at("rho_vac_UA") : 7.09e-36;
        return amplitude * rho_vac * std::sin(frequency * t + phase);
    }
    
    std::string getName() const override { return "DynamicVacuum"; }
    std::string getDescription() const override { return "Time-varying vacuum energy fluctuations"; }
    
    void optimize(double learningRate, double error) override {
        amplitude *= (1.0 - learningRate * error);
    }
    
    std::string serializeState() const override {
        std::ostringstream oss;
        oss << amplitude << "," << frequency << "," << phase;
        return oss.str();
    }
};

class QuantumCouplingTerm : public PhysicsTerm {
private:
    double coupling_strength;
    double decay_rate;
    
public:
    QuantumCouplingTerm(double strength = 1e-40, double decay = 1e-6)
        : coupling_strength(strength), decay_rate(decay) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        double hbar = params.count("hbar") ? params.at("hbar") : 1.0546e-34;
        double M = params.count("M") ? params.at("M") : 1.989e30;
        double r = params.count("r") ? params.at("r") : 1e4;
        return coupling_strength * (hbar * hbar) / (M * r * r) * std::cos(t * decay_rate);
    }
    
    std::string getName() const override { return "QuantumCoupling"; }
    std::string getDescription() const override { return "Non-local quantum gravitational coupling"; }
    
    void optimize(double learningRate, double error) override {
        coupling_strength *= (1.0 - learningRate * error * 0.1);
    }
};

class DarkMatterHaloTerm : public PhysicsTerm {
private:
    double halo_scale;
    double concentration;
    
public:
    DarkMatterHaloTerm(double scale = 1e20, double conc = 10.0)
        : halo_scale(scale), concentration(conc) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        (void)t;  // DM halo is typically static
        double r = params.count("r") ? params.at("r") : 1e4;
        double M = params.count("M") ? params.at("M") : 1.989e30;
        double G = params.count("G") ? params.at("G") : 6.6743e-11;
        // NFW-like profile contribution
        double x = r / halo_scale;
        double nfw_factor = std::log(1.0 + concentration * x) / (concentration * x);
        return G * M * nfw_factor / (r * r);
    }
    
    std::string getName() const override { return "DarkMatterHalo"; }
    std::string getDescription() const override { return "NFW dark matter halo contribution"; }
};

class MagneticReconnectionTerm : public PhysicsTerm {
private:
    double reconnection_rate;
    double B_initial;
    
public:
    MagneticReconnectionTerm(double rate = 0.1, double B0 = 1e-6)
        : reconnection_rate(rate), B_initial(B0) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        double B = params.count("B") ? params.at("B") : B_initial;
        double rho = params.count("rho_fluid") ? params.at("rho_fluid") : 1e-20;
        double v_A = B / std::sqrt(4.0 * M_PI * 1e-7 * rho);  // Alfven velocity
        return reconnection_rate * v_A * v_A / params.at("r") * std::exp(-t / 1e15);
    }
    
    std::string getName() const override { return "MagneticReconnection"; }
    std::string getDescription() const override { return "Magnetic reconnection energy release"; }
    
    bool validate(const std::map<std::string, double>& params) const override {
        return params.count("r") > 0;
    }
};

class TidalStrippingTerm : public PhysicsTerm {
private:
    double stripping_efficiency;
    
public:
    TidalStrippingTerm(double eff = 0.01) : stripping_efficiency(eff) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        double M = params.count("M") ? params.at("M") : 1.989e30;
        double r = params.count("r") ? params.at("r") : 1e4;
        double G = params.count("G") ? params.at("G") : 6.6743e-11;
        double tau = params.count("tau_tidal") ? params.at("tau_tidal") : 1e15;
        return stripping_efficiency * G * M / (r * r) * std::exp(-t / tau);
    }
    
    std::string getName() const override { return "TidalStripping"; }
    std::string getDescription() const override { return "Tidal stripping mass loss effect"; }
};

class RadiativeCoolingTerm : public PhysicsTerm {
private:
    double cooling_function;
    
public:
    RadiativeCoolingTerm(double cool = 1e-23) : cooling_function(cool) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        (void)t;
        double rho = params.count("rho_fluid") ? params.at("rho_fluid") : 1e-20;
        double T = params.count("temperature") ? params.at("temperature") : 1e6;
        // Cooling rate ~ n^2 * Lambda(T)
        double n = rho / 1.673e-27;  // number density
        return cooling_function * n * n * std::sqrt(T) / rho;
    }
    
    std::string getName() const override { return "RadiativeCooling"; }
    std::string getDescription() const override { return "Radiative cooling feedback"; }
};

// ===========================================================================================
// SelfExpandingModule: Template mixin for self-expanding capabilities
// ===========================================================================================

template<typename ConfigType>
class SelfExpandingModule {
protected:
    // Dynamic term registry
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    
    // Auto-expanding parameter map (accepts ANY key)
    std::map<std::string, double> params;
    
    // Metadata tracking
    std::map<std::string, std::string> metadata;
    
    // Self-optimization
    double learningRate;
    bool enableDynamicTerms;
    bool enableLogging;
    bool enableAutoOptimize;
    
    // Simulation history
    std::vector<std::pair<double, double>> simulationHistory;  // (time, value)
    
    std::string moduleName;
    
public:
    SelfExpandingModule(const std::string& name = "SelfExpandingModule")
        : learningRate(0.001),
          enableDynamicTerms(true),
          enableLogging(false),
          enableAutoOptimize(false),
          moduleName(name)
    {
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-SelfExpanding";
        metadata["created"] = getCurrentTimestamp();
    }
    
    virtual ~SelfExpandingModule() = default;
    
    // ==================== SELF-EXPANDING CAPABILITIES ====================
    
    // Register a new dynamic physics term at runtime
    void registerDynamicTerm(std::unique_ptr<PhysicsTerm> term) {
        if (term && term->validate(params)) {
            if (enableLogging) {
                std::cout << "[" << moduleName << "] Registered dynamic term: " 
                          << term->getName() << " - " << term->getDescription() << "\n";
            }
            metadata["last_term_added"] = term->getName();
            dynamicTerms.push_back(std::move(term));
        }
    }
    
    // Remove a dynamic term by name
    bool removeDynamicTerm(const std::string& name) {
        auto it = std::remove_if(dynamicTerms.begin(), dynamicTerms.end(),
            [&name](const std::unique_ptr<PhysicsTerm>& term) {
                return term->getName() == name;
            });
        if (it != dynamicTerms.end()) {
            dynamicTerms.erase(it, dynamicTerms.end());
            return true;
        }
        return false;
    }
    
    // Compute contribution from all dynamic terms
    double computeDynamicTerms(double t) const {
        if (!enableDynamicTerms) return 0.0;
        
        double total = 0.0;
        for (const auto& term : dynamicTerms) {
            if (term->validate(params)) {
                total += term->compute(t, params);
            }
        }
        return total;
    }
    
    // List all registered dynamic terms
    std::vector<std::string> listDynamicTerms() const {
        std::vector<std::string> names;
        for (const auto& term : dynamicTerms) {
            names.push_back(term->getName() + ": " + term->getDescription());
        }
        return names;
    }
    
    // ==================== AUTO-EXPANDING PARAMETERS ====================
    
    // Set ANY parameter (auto-creates if not exists)
    void setDynamicParameter(const std::string& key, double value) {
        bool isNew = (params.find(key) == params.end());
        params[key] = value;
        if (enableLogging) {
            std::cout << "[" << moduleName << "] " << (isNew ? "Created" : "Updated") 
                      << " parameter: " << key << " = " << value << "\n";
        }
        if (isNew) {
            metadata["last_param_created"] = key;
        }
    }
    
    // Get parameter (returns 0 if not found, with optional logging)
    double getDynamicParameter(const std::string& key) const {
        auto it = params.find(key);
        if (it != params.end()) {
            return it->second;
        }
        if (enableLogging) {
            std::cout << "[" << moduleName << "] Warning: Parameter '" << key << "' not found\n";
        }
        return 0.0;
    }
    
    // Check if parameter exists
    bool hasParameter(const std::string& key) const {
        return params.find(key) != params.end();
    }
    
    // Get all parameter names
    std::vector<std::string> listParameters() const {
        std::vector<std::string> keys;
        for (const auto& p : params) {
            keys.push_back(p.first);
        }
        return keys;
    }
    
    // ==================== SELF-SIMULATION ====================
    
    // Run time evolution simulation
    virtual void runSimulation(double t_start, double t_end, int steps,
                               std::function<double(double)> computeFunc) {
        simulationHistory.clear();
        double dt = (t_end - t_start) / steps;
        
        if (enableLogging) {
            std::cout << "[" << moduleName << "] Running simulation: t=" << t_start 
                      << " to " << t_end << " (" << steps << " steps)\n";
        }
        
        for (int i = 0; i <= steps; ++i) {
            double t = t_start + i * dt;
            double value = computeFunc(t);
            simulationHistory.push_back({t, value});
            
            // Auto-optimize if enabled
            if (enableAutoOptimize && i > 0) {
                double prev_value = simulationHistory[i-1].second;
                double error = std::abs(value - prev_value) / (std::abs(prev_value) + 1e-30);
                optimizeDynamicTerms(error);
            }
        }
        
        metadata["last_simulation"] = getCurrentTimestamp();
        metadata["simulation_steps"] = std::to_string(steps);
    }
    
    // Get simulation results
    const std::vector<std::pair<double, double>>& getSimulationHistory() const {
        return simulationHistory;
    }
    
    // Export simulation to file
    void exportSimulation(const std::string& filename) const {
        std::ofstream out(filename);
        if (out.is_open()) {
            out << "# " << moduleName << " Simulation Results\n";
            out << "# Time, Value\n";
            for (const auto& point : simulationHistory) {
                out << std::scientific << std::setprecision(6) 
                    << point.first << "," << point.second << "\n";
            }
            out.close();
            if (enableLogging) {
                std::cout << "[" << moduleName << "] Simulation exported to " << filename << "\n";
            }
        }
    }
    
    // ==================== SELF-OPTIMIZATION ====================
    
    void setLearningRate(double rate) {
        learningRate = rate;
        if (enableLogging) {
            std::cout << "[" << moduleName << "] Learning rate set to " << rate << "\n";
        }
    }
    
    double getLearningRate() const { return learningRate; }
    
    void setEnableAutoOptimize(bool enable) {
        enableAutoOptimize = enable;
        if (enableLogging) {
            std::cout << "[" << moduleName << "] Auto-optimize " 
                      << (enable ? "enabled" : "disabled") << "\n";
        }
    }
    
    void optimizeDynamicTerms(double error) {
        for (auto& term : dynamicTerms) {
            term->optimize(learningRate, error);
        }
    }
    
    // ==================== STATE PERSISTENCE ====================
    
    void exportState(const std::string& filename) const {
        std::ofstream out(filename);
        if (out.is_open()) {
            out << "# " << moduleName << " State Export\n";
            out << "# Version: " << metadata.at("version") << "\n";
            out << "# Timestamp: " << getCurrentTimestamp() << "\n";
            out << "\n[PARAMETERS]\n";
            for (const auto& p : params) {
                out << p.first << "=" << std::scientific << std::setprecision(10) << p.second << "\n";
            }
            out << "\n[METADATA]\n";
            for (const auto& m : metadata) {
                out << m.first << "=" << m.second << "\n";
            }
            out << "\n[DYNAMIC_TERMS]\n";
            for (const auto& term : dynamicTerms) {
                out << term->getName() << ":" << term->serializeState() << "\n";
            }
            out << "\n[SETTINGS]\n";
            out << "learningRate=" << learningRate << "\n";
            out << "enableDynamicTerms=" << enableDynamicTerms << "\n";
            out << "enableLogging=" << enableLogging << "\n";
            out << "enableAutoOptimize=" << enableAutoOptimize << "\n";
            out.close();
            
            if (enableLogging) {
                std::cout << "[" << moduleName << "] State exported to " << filename << "\n";
            }
        }
    }
    
    void importState(const std::string& filename) {
        std::ifstream in(filename);
        if (in.is_open()) {
            std::string line;
            std::string section;
            while (std::getline(in, line)) {
                if (line.empty() || line[0] == '#') continue;
                if (line[0] == '[') {
                    section = line;
                    continue;
                }
                size_t eq = line.find('=');
                if (eq != std::string::npos) {
                    std::string key = line.substr(0, eq);
                    std::string val = line.substr(eq + 1);
                    
                    if (section == "[PARAMETERS]") {
                        params[key] = std::stod(val);
                    } else if (section == "[METADATA]") {
                        metadata[key] = val;
                    } else if (section == "[SETTINGS]") {
                        if (key == "learningRate") learningRate = std::stod(val);
                        else if (key == "enableDynamicTerms") enableDynamicTerms = (val == "1");
                        else if (key == "enableLogging") enableLogging = (val == "1");
                        else if (key == "enableAutoOptimize") enableAutoOptimize = (val == "1");
                    }
                }
            }
            in.close();
            if (enableLogging) {
                std::cout << "[" << moduleName << "] State imported from " << filename << "\n";
            }
        }
    }
    
    // ==================== UTILITY ====================
    
    void setEnableLogging(bool enable) {
        enableLogging = enable;
        if (enable) {
            std::cout << "[" << moduleName << "] Logging enabled\n";
        }
    }
    
    void setEnableDynamicTerms(bool enable) {
        enableDynamicTerms = enable;
    }
    
    void setMetadata(const std::string& key, const std::string& value) {
        metadata[key] = value;
    }
    
    std::string getMetadata(const std::string& key) const {
        auto it = metadata.find(key);
        return (it != metadata.end()) ? it->second : "";
    }
    
    void printExpandedInfo() const {
        std::cout << "========================================\n";
        std::cout << "Self-Expanding Module: " << moduleName << "\n";
        std::cout << "Version: " << getMetadata("version") << "\n";
        std::cout << "========================================\n";
        std::cout << "Dynamic Terms Registered: " << dynamicTerms.size() << "\n";
        for (const auto& term : dynamicTerms) {
            std::cout << "  - " << term->getName() << ": " << term->getDescription() << "\n";
        }
        std::cout << "Parameters: " << params.size() << "\n";
        std::cout << "Learning Rate: " << learningRate << "\n";
        std::cout << "Auto-Optimize: " << (enableAutoOptimize ? "ON" : "OFF") << "\n";
        std::cout << "Dynamic Terms: " << (enableDynamicTerms ? "ON" : "OFF") << "\n";
        std::cout << "========================================\n";
    }
    
protected:
    static std::string getCurrentTimestamp() {
        auto now = std::chrono::system_clock::now();
        auto time = std::chrono::system_clock::to_time_t(now);
        std::string ts = std::ctime(&time);
        ts.pop_back();  // Remove newline
        return ts;
    }
};

}  // namespace UQFFExpanding

#endif // UQFF_SELF_EXPANDING_H
