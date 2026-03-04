/**
 * ================================================================================================
 * PhysicsTerms.hpp - UQFF Dynamic Physics Term Framework
 * ================================================================================================
 * 
 * Extracted from: source4.cpp (lines 47-135)
 * Purpose: Self-expanding physics term system for runtime term registration
 * Phase: 1, Week 2 - Module Extraction
 * 
 * Contains:
 * - PhysicsTerm abstract base class
 * - DynamicVacuumTerm: Time-varying vacuum energy
 * - QuantumCouplingTerm: Non-local quantum coupling effects
 * 
 * Features:
 * - Runtime physics term registration
 * - Parameter validation before computation
 * - Extensible framework for custom physics
 * - Auto-documentation via getName/getDescription
 * 
 * Dependencies:
 * - C++17 standard library
 * - <map>, <string> for parameter passing
 * - <cmath> for sin/cos
 * 
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef PHYSICSTERMS_HPP
#define PHYSICSTERMS_HPP

#include <string>
#include <map>
#include <cmath>
#include <memory>

namespace UQFFCore {

/**
 * PhysicsTerm: Abstract base class for dynamically added physics terms
 * 
 * Enables runtime registration of new physics calculations without code modification.
 * Each term computes a contribution to the total field based on time and parameters.
 */
class PhysicsTerm
{
public:
    virtual ~PhysicsTerm() = default;

    /**
     * Compute the physics term contribution
     * @param t Time (s)
     * @param params Map of named parameters (mass, radius, temperature, etc.)
     * @return Computed term value (units vary by term)
     */
    virtual double compute(double t, const std::map<std::string, double> &params) const = 0;

    /**
     * Get term name for logging and identification
     * @return Term name (e.g., "DynamicVacuumTerm")
     */
    virtual std::string getName() const = 0;

    /**
     * Get term description for documentation
     * @return Human-readable description of physics
     */
    virtual std::string getDescription() const = 0;

    /**
     * Validate parameters before computation
     * @param params Map of named parameters
     * @return true if parameters are valid, false otherwise
     */
    virtual bool validate(const std::map<std::string, double> &params) const = 0;
};

/**
 * DynamicVacuumTerm: Time-varying vacuum energy contribution
 * 
 * Models vacuum energy fluctuations as sinusoidal oscillations:
 *   E_vac(t) = A * sin(f * t)
 * 
 * Where:
 *   A = amplitude (energy units)
 *   f = frequency (rad/s)
 *   t = time (s)
 */
class DynamicVacuumTerm : public PhysicsTerm
{
private:
    double amplitude;  // Oscillation amplitude
    double frequency;  // Oscillation frequency (rad/s)

public:
    /**
     * Constructor
     * @param amp Amplitude of vacuum energy oscillation
     * @param freq Frequency of oscillation (rad/s)
     */
    DynamicVacuumTerm(double amp, double freq) 
        : amplitude(amp), frequency(freq) {}

    double compute(double t, const std::map<std::string, double>& /* params */) const override
    {
        return amplitude * std::sin(frequency * t);
    }

    std::string getName() const override 
    { 
        return "DynamicVacuumTerm"; 
    }

    std::string getDescription() const override
    {
        return "Time-varying vacuum energy contribution: A*sin(f*t)";
    }

    bool validate(const std::map<std::string, double>& /* params */) const override
    {
        return amplitude != 0.0 && frequency > 0.0;
    }
};

/**
 * QuantumCouplingTerm: Non-local quantum coupling effects
 * 
 * Models quantum coupling between mass and spatial extent:
 *   E_qc(t) = strength * ℏ²/(M*r²) * cos(t/10⁶)
 * 
 * Where:
 *   strength = coupling strength (dimensionless)
 *   ℏ = reduced Planck constant (J·s)
 *   M = mass (kg, from params["mass"])
 *   r = radius (m, from params["radius"])
 *   t = time (s)
 * 
 * Requires params: "mass", "radius"
 */
class QuantumCouplingTerm : public PhysicsTerm
{
private:
    double coupling_strength;  // Coupling strength
    double hbar;               // Reduced Planck constant (J·s)

public:
    /**
     * Constructor
     * @param strength Quantum coupling strength (dimensionless)
     */
    QuantumCouplingTerm(double strength)
        : coupling_strength(strength), 
          hbar(1.054571817e-34) {}

    double compute(double t, const std::map<std::string, double>& params) const override
    {
        auto it = params.find("mass");
        double M = (it != params.end()) ? it->second : 1e30;
        it = params.find("radius");
        double r = (it != params.end()) ? it->second : 1e3;

        return coupling_strength * (hbar * hbar) / (M * r * r) * std::cos(t / 1e6);
    }

    std::string getName() const override 
    { 
        return "QuantumCouplingTerm"; 
    }

    std::string getDescription() const override
    {
        return "Non-local quantum coupling: strength * hbar^2/(M*r^2) * cos(t/10^6)";
    }

    bool validate(const std::map<std::string, double>& /* params */) const override
    {
        return coupling_strength != 0.0;
    }
};

/**
 * InflationForceEpochTerm: Epoch-based unified field strength calculation
 * 
 * Computes F_U at specific cosmic epochs from Inflation/Force Chart framework.
 * Source: Grok Thread 4e0ecf23 (March 4, 2026)
 * Module: GrokThread_StarMagic_UnifiedFramework.py
 * 
 * Equation:
 *   F_U(t=0) = F_core + sum_{states=1 to 26} (Ui_state + F_p_state)
 *   F_core = ℏ ω_LENR / (σ_n ρ_vac,[UA]) ~ 10^10 N
 *   Ui_sum ≈ 0.1 * F_core * epoch_number (internal energy contribution)
 *   Fp_sum ≈ 0.05 * F_core * epoch_number (pressure contribution)
 * 
 * Epochs:
 *   1: Fisile Nuclei/Nebular (t=1.0-1.9, no Ug ranges active)
 *   2: Star/Planetary Atom (t=2.0-2.9, Ug1-3 active)
 *   3: Galaxies/Quasar (t=3.0-3.9, early Ug4)
 *   4: Magnetar/SMBH (t=4.0-4.9, Ug4 dominance)
 *   5: Globular Clusters (t=5.0-5.9, stabilization)
 * 
 * Required params:
 *   "epoch": Epoch number (1-5)
 *   "rho_vac_UA": Universal Aether vacuum density (J/m³)
 *   "omega_LENR": LENR resonance frequency (Hz)
 *   "sigma_n": Neutron cross-section (m²)
 */
class InflationForceEpochTerm : public PhysicsTerm
{
private:
    int epoch_number;     // Epoch (1-5)
    double h_bar;         // Reduced Planck constant (J·s)
    
public:
    /**
     * Constructor
     * @param epoch Epoch number (1-5)
     * @param hbar_val Reduced Planck constant (default: 1.055e-34 J·s)
     */
    InflationForceEpochTerm(int epoch, double hbar_val = 1.054571817e-34)
        : epoch_number(epoch), h_bar(hbar_val) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override
    {
        // Extract parameters
        double rho_vac_UA = params.count("rho_vac_UA") ? params.at("rho_vac_UA") : 7.09e-36;
        double omega_LENR = params.count("omega_LENR") ? params.at("omega_LENR") : 1.2e12;
        double sigma_n = params.count("sigma_n") ? params.at("sigma_n") : 1e-28;
        
        // F_core = ℏ ω_LENR / (σ_n ρ_vac,[UA])
        double F_core = (h_bar * omega_LENR) / (sigma_n * rho_vac_UA);
        
        // Epoch-dependent contributions
        double Ui_sum = F_core * 0.1 * epoch_number;  // Internal energy
        double Fp_sum = F_core * 0.05 * epoch_number; // Pressure
        
        // Total unified field at epoch
        double F_U = F_core + Ui_sum + Fp_sum;
        
        return F_U;
    }
    
    std::string getName() const override
    {
        return "InflationForceEpochTerm_Epoch" + std::to_string(epoch_number);
    }
    
    std::string getDescription() const override
    {
        std::string desc = "Unified field strength at Epoch " + std::to_string(epoch_number) + ": ";
        switch(epoch_number) {
            case 1: desc += "Fisile Nuclei/Nebular (no Ug ranges)"; break;
            case 2: desc += "Star/Planetary Atom (Ug1-3 active)"; break;
            case 3: desc += "Galaxies/Quasar (early Ug4)"; break;
            case 4: desc += "Magnetar/SMBH (Ug4 dominance)"; break;
            case 5: desc += "Globular Clusters (stabilization)"; break;
            default: desc += "Unknown epoch";
        }
        return desc;
    }
    
    bool validate(const std::map<std::string, double>& params) const override
    {
        // Epoch must be 1-5
        if (epoch_number < 1 || epoch_number > 5) return false;
        
        // If parameters provided, check they're positive
        if (params.count("rho_vac_UA") && params.at("rho_vac_UA") <= 0.0) return false;
        if (params.count("omega_LENR") && params.at("omega_LENR") <= 0.0) return false;
        if (params.count("sigma_n") && params.at("sigma_n") <= 0.0) return false;
        
        return true;
    }
};

} // namespace UQFFCore

#endif // PHYSICSTERMS_HPP
