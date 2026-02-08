// Wolfram-Enhanced MUGE Compressed Component Classes from source5.cpp
// Generated: November 29, 2025 (Phase 2)
// Source: UQFF Module 5 - Compressed MUGE Terms
// Total Classes: 9 (Modular breakdown of compute_compressed_MUGE)
// Status: PHASE 2 COMPLETE - Compressed MUGE modularization

#include <cmath>
#include <string>
#include <map>
#include <memory>
#include <stdexcept>

// ============================================================================
// MUGE COMPRESSED COMPONENT CLASSES (9 Classes)
// Modular breakdown of compute_compressed_MUGE (Source5.cpp lines 558-640)
// ============================================================================

// CLASS 22: Source5CompressedBaseTerm - Newtonian gravitational base
class Source5CompressedBaseTerm : public PhysicsTerm {
private:
    static constexpr double G = 6.67430e-11;  // Gravitational constant
    
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        auto it_M = params.find("M");
        auto it_r = params.find("r");
        
        double M = (it_M != params.end()) ? it_M->second : 1.989e30;  // Solar mass default
        double r = (it_r != params.end()) ? it_r->second : 1.496e11;  // 1 AU default
        
        if (r == 0.0) throw std::runtime_error("Division by zero in r");
        return G * M / (r * r);  // m/s²
    }
    
    std::string getName() const override { return "Source5_CompressedBase"; }
    
    std::string getDescription() const override {
        return "Compressed MUGE base: G*M/r^2 (Newtonian gravitational acceleration)";
    }
    
    bool validate(const std::map<std::string, double>& params) const override {
        auto it_M = params.find("M");
        auto it_r = params.find("r");
        return (it_M != params.end() && it_r != params.end() && it_M->second > 0 && it_r->second > 0);
    }
};

// CLASS 23: Source5ExpansionTerm - Hubble expansion modulation
class Source5ExpansionTerm : public PhysicsTerm {
private:
    static constexpr double H0 = 2.269e-18;  // Hubble constant (1/s)
    
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        auto it_t = params.find("t_age");
        double t_sys = (it_t != params.end()) ? it_t->second : t;
        
        double H_tz = H0 * t_sys;
        return 1.0 + H_tz;  // Dimensionless expansion factor
    }
    
    std::string getName() const override { return "Source5_Expansion"; }
    
    std::string getDescription() const override {
        return "Hubble expansion factor: 1 + H0*t where H0 = 2.269e-18 s^-1";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 24: Source5SuperAdjustmentTerm - Superconductive magnetic field adjustment
class Source5SuperAdjustmentTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        auto it_B = params.find("B");
        auto it_Bcrit = params.find("Bcrit");
        
        double B = (it_B != params.end()) ? it_B->second : 1e10;
        double Bcrit = (it_Bcrit != params.end()) ? it_Bcrit->second : 1e11;
        
        if (Bcrit == 0.0) throw std::runtime_error("Division by zero in Bcrit");
        return 1.0 - B / Bcrit;  // Dimensionless
    }
    
    std::string getName() const override { return "Source5_SuperconductiveAdjustment"; }
    
    std::string getDescription() const override {
        return "Superconductive magnetic adjustment: 1 - B/Bcrit (suppression factor)";
    }
    
    bool validate(const std::map<std::string, double>& params) const override {
        auto it_Bcrit = params.find("Bcrit");
        return (it_Bcrit != params.end() && it_Bcrit->second > 0);
    }
};

// CLASS 25: Source5EnvelopeTerm - Envelope modulation (placeholder)
class Source5EnvelopeTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        return 1.0;  // Neutral envelope (future extension point)
    }
    
    std::string getName() const override { return "Source5_Envelope"; }
    
    std::string getDescription() const override {
        return "Envelope modulation factor (currently neutral = 1.0, future stellar envelopes)";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 26: Source5UgSumTerm - Sum of Ug1-3 components (placeholder)
class Source5UgSumTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        return 0.0;  // Simplified (could sum Ug1-3 if needed)
    }
    
    std::string getName() const override { return "Source5_UgSum"; }
    
    std::string getDescription() const override {
        return "Sum of Ug1-3 components (simplified to 0, can aggregate Ug terms)";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 27: Source5CosmologicalTerm - Cosmological constant contribution
class Source5CosmologicalTerm : public PhysicsTerm {
private:
    static constexpr double Lambda = 1.1e-52;  // Cosmological constant (m^-2)
    static constexpr double c = 3.0e8;         // Speed of light (m/s)
    
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        return Lambda * c * c / 3.0;  // m/s²
    }
    
    std::string getName() const override { return "Source5_Cosmological"; }
    
    std::string getDescription() const override {
        return "Cosmological constant: Lambda*c^2/3 where Lambda = 1.1e-52 m^-2 (dark energy)";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 28: Source5QuantumTerm - Quantum uncertainty contribution
class Source5QuantumTerm : public PhysicsTerm {
private:
    static constexpr double hbar = 1.0546e-34;         // Reduced Planck constant (J·s)
    static constexpr double Delta_x_p = 1e-68;         // Uncertainty product (J·m)
    static constexpr double integral_psi = 2.176e-18;  // Wavefunction integral (m^-1)
    static constexpr double tHubble = 4.35e17;         // Hubble time (s)
    static constexpr double PI = 3.141592653589793;
    
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        if (Delta_x_p == 0.0) throw std::runtime_error("Division by zero in Delta_x_p");
        return (hbar / Delta_x_p) * integral_psi * (2.0 * PI / tHubble);  // m/s²
    }
    
    std::string getName() const override { return "Source5_Quantum"; }
    
    std::string getDescription() const override {
        return "Quantum uncertainty: (hbar/Delta_xp)*integral_psi*(2*PI/tHubble) (quantum correction)";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 29: Source5FluidTerm - Fluid dynamics contribution (Navier-Stokes coupling)
class Source5FluidTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        auto it_rho = params.find("rho_fluid");
        auto it_Vsys = params.find("Vsys");
        auto it_g = params.find("g_local");
        
        double rho_fluid = (it_rho != params.end()) ? it_rho->second : 1e-15;
        double Vsys = (it_Vsys != params.end()) ? it_Vsys->second : 4.189e12;
        double g_local = (it_g != params.end()) ? it_g->second : 10.0;
        
        return rho_fluid * Vsys * g_local;  // kg·m/s² (force)
    }
    
    std::string getName() const override { return "Source5_Fluid"; }
    
    std::string getDescription() const override {
        return "Fluid dynamics: rho_fluid * Vsys * g_local (Navier-Stokes coupling, units: N)";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 30: Source5PerturbationTerm - Dark matter + density fluctuations
class Source5PerturbationTerm : public PhysicsTerm {
private:
    static constexpr double G = 6.67430e-11;  // Gravitational constant
    
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        auto it_M = params.find("M");
        auto it_M_DM = params.find("M_DM");
        auto it_r = params.find("r");
        auto it_delta = params.find("delta_rho_rho");
        
        double M = (it_M != params.end()) ? it_M->second : 1.989e30;
        double M_DM = (it_M_DM != params.end()) ? it_M_DM->second : 5e30;
        double r = (it_r != params.end()) ? it_r->second : 1.496e11;
        double delta_rho_rho = (it_delta != params.end()) ? it_delta->second : 1e-5;
        
        if (r == 0.0) throw std::runtime_error("Division by zero in r^3");
        
        return (M + M_DM) * (delta_rho_rho + 3.0 * G * M / (r * r * r));  // m/s²
    }
    
    std::string getName() const override { return "Source5_Perturbation"; }
    
    std::string getDescription() const override {
        return "Perturbation: (M+M_DM)*(delta_rho/rho + 3GM/r^3) (dark matter + density fluctuations)";
    }
    
    bool validate(const std::map<std::string, double>& params) const override {
        auto it_r = params.find("r");
        return (it_r != params.end() && it_r->second > 0);
    }
};

// ============================================================================
// COMBINED COMPRESSED MUGE TERM (Wrapper)
// ============================================================================

class Source5CompressedMUGETerm : public PhysicsTerm {
private:
    std::unique_ptr<Source5CompressedBaseTerm> base_term;
    std::unique_ptr<Source5ExpansionTerm> expansion_term;
    std::unique_ptr<Source5SuperAdjustmentTerm> super_adj_term;
    std::unique_ptr<Source5EnvelopeTerm> envelope_term;
    std::unique_ptr<Source5UgSumTerm> ug_sum_term;
    std::unique_ptr<Source5CosmologicalTerm> cosm_term;
    std::unique_ptr<Source5QuantumTerm> quantum_term;
    std::unique_ptr<Source5FluidTerm> fluid_term;
    std::unique_ptr<Source5PerturbationTerm> pert_term;
    
public:
    Source5CompressedMUGETerm() 
        : base_term(std::make_unique<Source5CompressedBaseTerm>()),
          expansion_term(std::make_unique<Source5ExpansionTerm>()),
          super_adj_term(std::make_unique<Source5SuperAdjustmentTerm>()),
          envelope_term(std::make_unique<Source5EnvelopeTerm>()),
          ug_sum_term(std::make_unique<Source5UgSumTerm>()),
          cosm_term(std::make_unique<Source5CosmologicalTerm>()),
          quantum_term(std::make_unique<Source5QuantumTerm>()),
          fluid_term(std::make_unique<Source5FluidTerm>()),
          pert_term(std::make_unique<Source5PerturbationTerm>()) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Compute base components
        double base = base_term->compute(t, params);
        double expansion = expansion_term->compute(t, params);
        double super_adj = super_adj_term->compute(t, params);
        double env = envelope_term->compute(t, params);
        
        // Adjusted base = base * expansion * super_adj * env
        double adjusted_base = base * expansion * super_adj * env;
        
        // Add remaining terms
        double Ug_sum = ug_sum_term->compute(t, params);
        double cosm = cosm_term->compute(t, params);
        double quantum = quantum_term->compute(t, params);
        double fluid = fluid_term->compute(t, params);
        double perturbation = pert_term->compute(t, params);
        
        return adjusted_base + Ug_sum + cosm + quantum + fluid + perturbation;
    }
    
    std::string getName() const override { return "Source5_CompressedMUGE_Full"; }
    
    std::string getDescription() const override {
        return "Full Compressed MUGE: base*expansion*super_adj*env + Ug_sum + cosm + quantum + fluid + perturbation";
    }
    
    bool validate(const std::map<std::string, double>& params) const override {
        return base_term->validate(params) && super_adj_term->validate(params);
    }
};

// ============================================================================
// REGISTRATION FUNCTION
// ============================================================================

void registerCompressedTerms_source5(PhysicsTermRegistry& registry) {
    // Individual Component Classes (9)
    registry.registerPhysicsTerm("Source5_CompressedBase", std::make_unique<Source5CompressedBaseTerm>(), "wolfram_compressed");
    registry.registerPhysicsTerm("Source5_Expansion", std::make_unique<Source5ExpansionTerm>(), "wolfram_compressed");
    registry.registerPhysicsTerm("Source5_SuperconductiveAdjustment", std::make_unique<Source5SuperAdjustmentTerm>(), "wolfram_compressed");
    registry.registerPhysicsTerm("Source5_Envelope", std::make_unique<Source5EnvelopeTerm>(), "wolfram_compressed");
    registry.registerPhysicsTerm("Source5_UgSum", std::make_unique<Source5UgSumTerm>(), "wolfram_compressed");
    registry.registerPhysicsTerm("Source5_Cosmological", std::make_unique<Source5CosmologicalTerm>(), "wolfram_compressed");
    registry.registerPhysicsTerm("Source5_Quantum", std::make_unique<Source5QuantumTerm>(), "wolfram_compressed");
    registry.registerPhysicsTerm("Source5_Fluid", std::make_unique<Source5FluidTerm>(), "wolfram_compressed");
    registry.registerPhysicsTerm("Source5_Perturbation", std::make_unique<Source5PerturbationTerm>(), "wolfram_compressed");
    
    // Combined Wrapper (1)
    registry.registerPhysicsTerm("Source5_CompressedMUGE_Full", std::make_unique<Source5CompressedMUGETerm>(), "wolfram_compressed");
}

// ============================================================================
// TOTAL: 10 CLASSES (9 components + 1 wrapper)
// - CompressedBase: G*M/r^2
// - Expansion: 1 + H0*t
// - SuperconductiveAdjustment: 1 - B/Bcrit
// - Envelope: Neutral placeholder
// - UgSum: Ug1-3 aggregation placeholder
// - Cosmological: Lambda*c^2/3
// - Quantum: Quantum uncertainty correction
// - Fluid: Navier-Stokes coupling
// - Perturbation: Dark matter + density fluctuations
// - CompressedMUGE_Full: Complete equation wrapper
// ============================================================================
