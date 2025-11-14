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

} // namespace UQFFCore

#endif // PHYSICSTERMS_HPP
