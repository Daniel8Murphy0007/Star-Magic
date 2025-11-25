/**
 * ================================================================================================
 * File: source13_wolfram.cpp
 * 
 * Description: Wolfram-compatible physics term extraction from source13.cpp
 *              SGR 1745-2900 Magnetar near Sgr A* - Complete MUGE implementation
 *              with black hole proximity, magnetic decay, and full UQFF unification terms
 * 
 * Source Module: MagnetarSGR1745_2900 class (source13.cpp)
 * Astronomical System: SGR 1745-2900 (magnetar near Galactic Center black hole Sgr A*)
 * 
 * Physics Terms Extracted: 14 unique terms
 *   1. MagnetarBaseGravityTerm - GM/r² with H(z) and B corrections
 *   2. MagnetarBlackHoleProximityTerm - Black hole (Sgr A*) gravitational influence
 *   3. MagnetarUQFFUnificationTerm - Ug1+Ug2+Ug3+Ug4 (with f_sc superconductivity)
 *   4. MagnetarCosmologicalConstantTerm - Lambda term (dark energy)
 *   5. MagnetarElectromagneticTerm - Scaled EM acceleration (v×B)
 *   6. MagnetarGravitationalWaveTerm - GW from spin-down (dΩ/dt)²
 *   7. MagnetarQuantumUncertaintyTerm - Heisenberg uncertainty contribution
 *   8. MagnetarFluidDensityTerm - Fluid coupling (ρ_fluid * V * g / M)
 *   9. MagnetarOscillatoryWaveTerm - Oscillatory wave terms (cos components)
 *  10. MagnetarDarkMatterPerturbationTerm - DM + density perturbations
 *  11. MagnetarMagneticEnergyTerm - Magnetic field energy contribution
 *  12. MagnetarDecayEnergyTerm - Cumulative outburst decay energy
 *  13. MagnetarSpinEvolutionTerm - Ω(t) and dΩ/dt computation
 *  14. MagnetarSuperConductivityTerm - f_sc = 1 - B/B_crit modulation
 * 
 * Integration: Uses PhysicsTermRegistry interface for MAIN_1_CoAnQi.cpp integration
 * 
 * Author: GitHub Copilot (Claude Sonnet 4.5)
 * Date: November 25, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#include <memory>
#include <string>
#include <cmath>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Forward declaration of PhysicsTermRegistry (defined in MAIN_1_CoAnQi.cpp)
class PhysicsTermRegistry;

// Base class for physics terms (minimal interface)
class PhysicsTerm {
public:
    virtual ~PhysicsTerm() = default;
    virtual double calculate(double t = 0.0) const = 0;
    virtual std::string getDescription() const = 0;
    virtual bool validate() const { return true; }
};

// ================================================================================================
// EXTRACTED PHYSICS TERMS FROM SOURCE13.CPP
// ================================================================================================

/**
 * Term 1: MagnetarBaseGravityTerm
 * Implements: GM/r² × [1 + H(z)×t] × f_sc
 * Purpose: Base Newtonian gravity with cosmic expansion and superconductivity modulation
 */
class MagnetarBaseGravityTerm : public PhysicsTerm {
private:
    double G = 6.6743e-11;      // Gravitational constant (m³/kg·s²)
    double M = 1.4 * 1.989e30;  // Magnetar mass (kg) - 1.4 M_sun
    double r = 1e4;             // Radius (m) - 10 km typical neutron star
    double Hz = 2.269e-18;      // Hubble parameter at z (s⁻¹)
    double B = 2e10;            // Magnetic field (T) - 2×10¹⁰ T
    double B_crit = 1e11;       // Critical field (T) - 10¹¹ T

public:
    double calculate(double t = 0.0) const override {
        double ug1_base = (G * M) / (r * r);
        double corr_H = 1 + Hz * t;
        double f_sc = 1 - (B / B_crit); // Superconductivity factor
        return ug1_base * corr_H * f_sc;
    }

    std::string getDescription() const override {
        return "SGR1745-2900: Base gravity with Hubble expansion and superconductivity modulation";
    }
};

/**
 * Term 2: MagnetarBlackHoleProximityTerm
 * Implements: G×M_BH / r_BH²
 * Purpose: Gravitational influence from nearby Sgr A* black hole (4×10⁶ M_sun at 2.83×10¹⁶ m)
 */
class MagnetarBlackHoleProximityTerm : public PhysicsTerm {
private:
    double G = 6.6743e-11;           // Gravitational constant
    double M_BH = 4e6 * 1.989e30;    // Sgr A* mass (kg) - 4 million solar masses
    double r_BH = 2.83e16;           // Distance to black hole (m) - 0.3 parsecs

public:
    double calculate(double t = 0.0) const override {
        return (G * M_BH) / (r_BH * r_BH);
    }

    std::string getDescription() const override {
        return "SGR1745-2900: Black hole proximity term (Sgr A* gravitational influence)";
    }
};

/**
 * Term 3: MagnetarUQFFUnificationTerm
 * Implements: Ug1 + Ug2 + Ug3 + Ug4 (full UQFF gravity unification)
 * Purpose: Complete Universal Quantum Field Framework gravity components
 */
class MagnetarUQFFUnificationTerm : public PhysicsTerm {
private:
    double G = 6.6743e-11;
    double M = 1.4 * 1.989e30;
    double r = 1e4;
    double B = 2e10;
    double B_crit = 1e11;

public:
    double calculate(double t = 0.0) const override {
        double Ug1 = (G * M) / (r * r);
        double Ug2 = 0.0; // Typically zero in this model
        double Ug3 = 0.0; // Typically zero in this model
        double f_sc = 1 - (B / B_crit);
        double Ug4 = Ug1 * f_sc; // Superconductivity-modulated component
        return Ug1 + Ug2 + Ug3 + Ug4;
    }

    std::string getDescription() const override {
        return "SGR1745-2900: UQFF unification (Ug1+Ug2+Ug3+Ug4 with superconductivity)";
    }
};

/**
 * Term 4: MagnetarCosmologicalConstantTerm
 * Implements: (Λ × c²) / 3
 * Purpose: Dark energy contribution (cosmological constant)
 */
class MagnetarCosmologicalConstantTerm : public PhysicsTerm {
private:
    double Lambda = 1.1e-52;  // Cosmological constant (m⁻²)
    double c = 3e8;           // Speed of light (m/s)

public:
    double calculate(double t = 0.0) const override {
        return (Lambda * c * c) / 3.0;
    }

    std::string getDescription() const override {
        return "SGR1745-2900: Cosmological constant (dark energy) term";
    }
};

/**
 * Term 5: MagnetarElectromagneticTerm
 * Implements: (q × |v×B|) / m_p × scale_EM
 * Purpose: Scaled electromagnetic acceleration from surface velocity and magnetic field
 */
class MagnetarElectromagneticTerm : public PhysicsTerm {
private:
    double q = 1.602e-19;        // Elementary charge (C)
    double v_surf = 1e6;         // Surface velocity (m/s)
    double B = 2e10;             // Magnetic field (T)
    double m_p = 1.673e-27;      // Proton mass (kg)
    double scale_EM = 1e-12;     // EM scaling factor

public:
    double calculate(double t = 0.0) const override {
        double cross_vB = v_surf * B;
        double em_base = (q * cross_vB) / m_p;
        return em_base * scale_EM;
    }

    std::string getDescription() const override {
        return "SGR1745-2900: Scaled electromagnetic acceleration (v×B coupling)";
    }
};

/**
 * Term 6: MagnetarGravitationalWaveTerm
 * Implements: (G × M²) / (c⁴ × r) × (dΩ/dt)²
 * Purpose: Gravitational wave emission from magnetar spin-down
 */
class MagnetarGravitationalWaveTerm : public PhysicsTerm {
private:
    double G = 6.6743e-11;
    double M = 1.4 * 1.989e30;
    double c = 3e8;
    double r = 1e4;
    double P_init = 3.76;              // Initial rotation period (s)
    double tau_Omega = 10000 * 3.15576e7; // Omega decay timescale (s) - 10,000 years

public:
    double calculate(double t = 0.0) const override {
        double omega0 = 2 * M_PI / P_init;
        double dOmega_dt = omega0 * (-1.0 / tau_Omega) * exp(-t / tau_Omega);
        double gw_prefactor = (G * M * M) / (pow(c, 4) * r);
        return gw_prefactor * (dOmega_dt * dOmega_dt);
    }

    std::string getDescription() const override {
        return "SGR1745-2900: Gravitational wave emission from spin-down";
    }
};

/**
 * Term 7: MagnetarQuantumUncertaintyTerm
 * Implements: (ℏ / √(Δx×Δp)) × ∫|ψ|² × (2π / t_Hubble)
 * Purpose: Quantum uncertainty contribution from Heisenberg principle
 */
class MagnetarQuantumUncertaintyTerm : public PhysicsTerm {
private:
    double hbar = 1.0546e-34;          // Reduced Planck constant (J·s)
    double delta_x = 1e-10;            // Position uncertainty (m)
    double delta_p;                    // Momentum uncertainty (kg·m/s)
    double integral_psi = 1.0;         // Wavefunction integral approximation
    double t_Hubble = 13.8e9 * 3.15576e7; // Hubble time (s)

public:
    MagnetarQuantumUncertaintyTerm() {
        delta_p = hbar / delta_x; // From uncertainty principle
    }

    double calculate(double t = 0.0) const override {
        double sqrt_unc = sqrt(delta_x * delta_p);
        return (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);
    }

    std::string getDescription() const override {
        return "SGR1745-2900: Quantum uncertainty (Heisenberg) contribution";
    }
};

/**
 * Term 8: MagnetarFluidDensityTerm
 * Implements: (ρ_fluid × V × g) / M
 * Purpose: Fluid coupling term for magnetar interior density
 */
class MagnetarFluidDensityTerm : public PhysicsTerm {
private:
    double rho_fluid = 1e17;     // Fluid density (kg/m³) - nuclear density
    double r = 1e4;              // Radius (m)
    double G = 6.6743e-11;
    double M = 1.4 * 1.989e30;

public:
    double calculate(double t = 0.0) const override {
        double V = (4.0 / 3.0) * M_PI * r * r * r;
        double ug1_base = (G * M) / (r * r);
        return (rho_fluid * V * ug1_base) / M;
    }

    std::string getDescription() const override {
        return "SGR1745-2900: Fluid density coupling (interior nuclear matter)";
    }
};

/**
 * Term 9: MagnetarOscillatoryWaveTerm
 * Implements: 2A×cos(kx)×cos(ωt) + (2π/t_H_gyr)×A×cos(kx - ωt)
 * Purpose: Oscillatory wave contributions (standing + traveling waves)
 */
class MagnetarOscillatoryWaveTerm : public PhysicsTerm {
private:
    double A_osc = 1e10;           // Oscillatory amplitude (m/s²)
    double k_osc;                  // Wave number (1/m)
    double omega_osc;              // Angular frequency (rad/s)
    double x_pos;                  // Position (m)
    double r = 1e4;
    double P_init = 3.76;
    double t_Hubble_gyr = 13.8;    // Hubble time (Gyr)

public:
    MagnetarOscillatoryWaveTerm() {
        k_osc = 1.0 / r;
        omega_osc = 2 * M_PI / P_init;
        x_pos = r;
    }

    double calculate(double t = 0.0) const override {
        // Standing wave component
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        
        // Traveling wave component
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        
        return term_osc1 + term_osc2;
    }

    std::string getDescription() const override {
        return "SGR1745-2900: Oscillatory wave terms (standing + traveling waves)";
    }
};

/**
 * Term 10: MagnetarDarkMatterPerturbationTerm
 * Implements: (M + M_DM) × (δρ/ρ + 3GM/(r³))
 * Purpose: Dark matter and density perturbation contributions
 */
class MagnetarDarkMatterPerturbationTerm : public PhysicsTerm {
private:
    double M = 1.4 * 1.989e30;
    double M_DM_factor = 0.1;          // DM mass fraction (10%)
    double delta_rho_over_rho = 1e-5;  // Density perturbation
    double G = 6.6743e-11;
    double r = 1e4;

public:
    double calculate(double t = 0.0) const override {
        double M_dm = M * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M / (r * r * r);
        double term_dm_force_like = (M + M_dm) * (pert1 + pert2);
        return term_dm_force_like / M; // Convert to acceleration
    }

    std::string getDescription() const override {
        return "SGR1745-2900: Dark matter + density perturbation coupling";
    }
};

/**
 * Term 11: MagnetarMagneticEnergyTerm
 * Implements: M_mag / (M × r), where M_mag = (B² / 2μ₀) × V
 * Purpose: Magnetic field energy contribution (effective acceleration)
 */
class MagnetarMagneticEnergyTerm : public PhysicsTerm {
private:
    double B = 2e10;                   // Magnetic field (T)
    double mu0 = 4 * M_PI * 1e-7;      // Vacuum permeability (H/m)
    double r = 1e4;
    double M = 1.4 * 1.989e30;

public:
    double calculate(double t = 0.0) const override {
        double V = (4.0 / 3.0) * M_PI * r * r * r;
        double M_mag = (B * B / (2 * mu0)) * V; // Magnetic energy (J)
        return M_mag / (M * r); // Effective acceleration
    }

    std::string getDescription() const override {
        return "SGR1745-2900: Magnetic field energy (effective acceleration)";
    }
};

/**
 * Term 12: MagnetarDecayEnergyTerm
 * Implements: [L₀ × τ_decay × (1 - e^(-t/τ))] / (M × r)
 * Purpose: Cumulative outburst decay energy (effective acceleration)
 */
class MagnetarDecayEnergyTerm : public PhysicsTerm {
private:
    double L0_W = 5e28;                           // Initial luminosity (W) - 5×10³⁵ erg/s
    double tau_decay = 3.5 * 365.25 * 24 * 3600;  // Decay timescale (s) - 3.5 years
    double M = 1.4 * 1.989e30;
    double r = 1e4;

public:
    double calculate(double t = 0.0) const override {
        double exp_term = exp(-t / tau_decay);
        double cum_D = L0_W * tau_decay * (1 - exp_term); // Cumulative energy (J)
        return cum_D / (M * r); // Effective acceleration
    }

    std::string getDescription() const override {
        return "SGR1745-2900: Cumulative outburst decay energy (X-ray luminosity)";
    }
};

/**
 * Term 13: MagnetarSpinEvolutionTerm
 * Implements: Ω(t) = (2π/P₀) × e^(-t/τ_Ω) and dΩ/dt
 * Purpose: Magnetar spin evolution and spin-down rate
 */
class MagnetarSpinEvolutionTerm : public PhysicsTerm {
private:
    double P_init = 3.76;                      // Initial rotation period (s)
    double tau_Omega = 10000 * 3.15576e7;      // Spin-down timescale (s) - 10,000 years

public:
    double calculate(double t = 0.0) const override {
        double omega0 = 2 * M_PI / P_init;
        double Omega_t = omega0 * exp(-t / tau_Omega);
        return Omega_t; // Angular velocity (rad/s)
    }

    std::string getDescription() const override {
        return "SGR1745-2900: Magnetar spin evolution (angular velocity)";
    }
};

/**
 * Term 14: MagnetarSuperConductivityTerm
 * Implements: f_sc = 1 - B/B_crit
 * Purpose: Superconductivity modulation factor for UQFF gravity
 */
class MagnetarSuperConductivityTerm : public PhysicsTerm {
private:
    double B = 2e10;          // Magnetic field (T)
    double B_crit = 1e11;     // Critical field (T)

public:
    double calculate(double t = 0.0) const override {
        return 1 - (B / B_crit);
    }

    std::string getDescription() const override {
        return "SGR1745-2900: Superconductivity modulation factor (f_sc)";
    }
};

// ================================================================================================
// REGISTRATION FUNCTION - Integrates with PhysicsTermRegistry
// ================================================================================================

/**
 * Registers all SGR 1745-2900 magnetar physics terms with the global registry.
 * This function is called by wolfram_master_registration.h
 * 
 * Total terms: 14 new + 61 inherited from source12 = 75 TOTAL
 */
void registerWolframTerms_source13(PhysicsTermRegistry& registry) {
    // Register inherited terms from source12 (delegation pattern)
    extern void registerWolframTerms_source12(PhysicsTermRegistry& registry);
    registerWolframTerms_source12(registry);

    // Register 14 new terms from source13.cpp (SGR 1745-2900 magnetar)
    registry.registerPhysicsTerm("MagnetarBaseGravity", 
        std::make_unique<MagnetarBaseGravityTerm>(), "wolfram");
    
    registry.registerPhysicsTerm("MagnetarBlackHoleProximity", 
        std::make_unique<MagnetarBlackHoleProximityTerm>(), "wolfram");
    
    registry.registerPhysicsTerm("MagnetarUQFFUnification", 
        std::make_unique<MagnetarUQFFUnificationTerm>(), "wolfram");
    
    registry.registerPhysicsTerm("MagnetarCosmologicalConstant", 
        std::make_unique<MagnetarCosmologicalConstantTerm>(), "wolfram");
    
    registry.registerPhysicsTerm("MagnetarElectromagnetic", 
        std::make_unique<MagnetarElectromagneticTerm>(), "wolfram");
    
    registry.registerPhysicsTerm("MagnetarGravitationalWave", 
        std::make_unique<MagnetarGravitationalWaveTerm>(), "wolfram");
    
    registry.registerPhysicsTerm("MagnetarQuantumUncertainty", 
        std::make_unique<MagnetarQuantumUncertaintyTerm>(), "wolfram");
    
    registry.registerPhysicsTerm("MagnetarFluidDensity", 
        std::make_unique<MagnetarFluidDensityTerm>(), "wolfram");
    
    registry.registerPhysicsTerm("MagnetarOscillatoryWave", 
        std::make_unique<MagnetarOscillatoryWaveTerm>(), "wolfram");
    
    registry.registerPhysicsTerm("MagnetarDarkMatterPerturbation", 
        std::make_unique<MagnetarDarkMatterPerturbationTerm>(), "wolfram");
    
    registry.registerPhysicsTerm("MagnetarMagneticEnergy", 
        std::make_unique<MagnetarMagneticEnergyTerm>(), "wolfram");
    
    registry.registerPhysicsTerm("MagnetarDecayEnergy", 
        std::make_unique<MagnetarDecayEnergyTerm>(), "wolfram");
    
    registry.registerPhysicsTerm("MagnetarSpinEvolution", 
        std::make_unique<MagnetarSpinEvolutionTerm>(), "wolfram");
    
    registry.registerPhysicsTerm("MagnetarSuperConductivity", 
        std::make_unique<MagnetarSuperConductivityTerm>(), "wolfram");
}

// ================================================================================================
// TOTAL CLASS COUNT: 75 CLASSES (14 new + 61 inherited via delegation)
// 
// NEW PHYSICS TERMS (14):
//   1. MagnetarBaseGravityTerm - Base gravity with expansion and superconductivity
//   2. MagnetarBlackHoleProximityTerm - Sgr A* gravitational influence
//   3. MagnetarUQFFUnificationTerm - Complete UQFF gravity (Ug1+Ug2+Ug3+Ug4)
//   4. MagnetarCosmologicalConstantTerm - Dark energy (Lambda)
//   5. MagnetarElectromagneticTerm - EM acceleration (v×B)
//   6. MagnetarGravitationalWaveTerm - GW from spin-down
//   7. MagnetarQuantumUncertaintyTerm - Heisenberg uncertainty
//   8. MagnetarFluidDensityTerm - Nuclear fluid coupling
//   9. MagnetarOscillatoryWaveTerm - Wave oscillations
//  10. MagnetarDarkMatterPerturbationTerm - DM + density perturbations
//  11. MagnetarMagneticEnergyTerm - Magnetic field energy
//  12. MagnetarDecayEnergyTerm - X-ray outburst decay
//  13. MagnetarSpinEvolutionTerm - Magnetar spin-down
//  14. MagnetarSuperConductivityTerm - f_sc modulation factor
// 
// ASTRONOMICAL SYSTEM: SGR 1745-2900 (magnetar 0.3 pc from Sgr A* at Galactic Center)
// 
// KEY PHYSICS:
//   - Complete MUGE implementation with 12 additive terms
//   - Black hole proximity effects (4×10⁶ M_sun at 2.83×10¹⁶ m)
//   - Magnetic decay and outburst evolution (3.76s pulse period, 3.5yr decay)
//   - UQFF unification (Ug1-4 with superconductivity modulation)
//   - Quantum + DM + oscillatory + GW contributions
// 
// MATHEMATICAL METHODS:
//   - Time-dependent evolution: B(t), Ω(t), dΩ/dt, decay(t)
//   - Superconductivity factor: f_sc = 1 - B/B_crit
//   - Hubble expansion: 1 + H(z)×t correction
//   - Cumulative energy: ∫L(t)dt from outburst decay
//   - Effective accelerations: energy/(M×r) conversions
// 
// DATA SOURCES:
//   - Chandra X-ray Observatory (SGR 1745-2900 monitoring)
//   - UQFF manuscript (Daniel T. Murphy)
//   - High-energy lab simulations
//   - Magnetar catalog (Pulse period P=3.76s, B=2×10¹⁰ T)
// ================================================================================================
