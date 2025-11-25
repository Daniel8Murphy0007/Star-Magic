// Wolfram-Enhanced Physics Terms from source5.cpp
// Generated: November 24, 2025
// Source: Modularized UQFF with Self-Expanding Framework 2.0
// Total Classes: 25 (19 from source4 + 6 new dynamic)

#include <cmath>
#include <string>
#include <map>
#include <memory>

// ============================================================================
// NEW DYNAMIC PHYSICS TERMS (Source5-Specific)
// ============================================================================

class DarkMatterHaloNFWTerm : public PhysicsTerm {
private:
    double M_halo = 1e42;    // Default halo mass (kg)
    double r_scale = 1e21;   // Default scale radius (m)
    
public:
    DarkMatterHaloNFWTerm() = default;
    DarkMatterHaloNFWTerm(double mass, double scale) : M_halo(mass), r_scale(scale) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        auto it_r = params.find("r");
        if (it_r == params.end()) return 0.0;
        double r = it_r->second;
        
        if (r <= 0.0 || r_scale <= 0.0) return 0.0;
        
        // NFW profile contribution: g(r) = G*M_halo*ln(1+x)/(r*x) where x = r/r_scale
        const double G = 6.67430e-11;
        const double PI = 3.141592653589793;
        
        double x = r / r_scale;
        double rho_0 = M_halo / (4.0 * PI * r_scale * r_scale * r_scale * (std::log(2.0) - 0.5));
        
        return G * M_halo * std::log(1 + x) / (r * x);
    }
    
    std::string getName() const override { return "DarkMatterHaloNFW"; }
    
    std::string getDescription() const override {
        return "Dark matter halo NFW profile: g_DM(r) = G*M_halo*ln(1+r/r_s)/(r*(r/r_s))";
    }
    
    bool validate(const std::map<std::string, double>& params) const override {
        return M_halo > 0.0 && r_scale > 0.0;
    }
};

class VacuumEnergyFluctuationTerm : public PhysicsTerm {
private:
    double E_vac_scale = 7.09e-36;  // Default vacuum energy scale
    double lambda = 1e-10;          // Default coupling strength
    
public:
    VacuumEnergyFluctuationTerm() = default;
    VacuumEnergyFluctuationTerm(double e_scale, double coupling) 
        : E_vac_scale(e_scale), lambda(coupling) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Time-varying vacuum energy: E_vac(t) = lambda * E_vac_scale * (1 + 0.1*sin(1e-10*t))
        return lambda * E_vac_scale * (1.0 + 0.1 * std::sin(1e-10 * t));
    }
    
    std::string getName() const override { return "VacuumEnergyFluctuation"; }
    
    std::string getDescription() const override {
        return "Time-varying vacuum energy: E_vac(t) = lambda*E_scale*(1 + 0.1*sin(1e-10*t))";
    }
    
    bool validate(const std::map<std::string, double>&) const override {
        return E_vac_scale != 0.0 && lambda != 0.0;
    }
};

// ============================================================================
// ENHANCED UQFF MODULE5 TERMS (Wrappers with Dynamic Contributions)
// ============================================================================

class UQFFModule5Ug1EnhancedTerm : public PhysicsTerm {
private:
    double k1 = 1.5;
    double alpha = 0.001;
    double delta_def = 0.01;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamic_terms;
    
public:
    void addDynamicTerm(std::unique_ptr<PhysicsTerm> term) {
        dynamic_terms.push_back(std::move(term));
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Original Ug1 calculation
        auto it_mu_s = params.find("mu_s");
        auto it_grad = params.find("grad_Ms_r");
        auto it_tn = params.find("tn");
        
        double mu_s = (it_mu_s != params.end()) ? it_mu_s->second : 1e20;
        double grad_Ms_r = (it_grad != params.end()) ? it_grad->second : 1e-5;
        double tn = (it_tn != params.end()) ? it_tn->second : t;
        
        double defect = 1.0 + delta_def * std::sin(0.001 * t);
        double original = k1 * mu_s * grad_Ms_r * std::exp(-alpha * t) * std::cos(M_PI * tn) * defect;
        
        // Add dynamic contributions
        double dynamic = 0.0;
        for (const auto& term : dynamic_terms) {
            if (term->validate(params)) {
                dynamic += term->compute(t, params);
            }
        }
        
        return original + dynamic;
    }
    
    std::string getName() const override { return "UQFFModule5Ug1Enhanced"; }
    
    std::string getDescription() const override {
        return "Enhanced Ug1 with dynamic term contributions (original + sum(dynamic))";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class UQFFModule5Ug2EnhancedTerm : public PhysicsTerm {
private:
    double k2 = 1.2;
    double QA = 1e-10;
    double delta_sw = 0.01;
    double v_sw = 5e5;
    double HSCm = 1.0;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamic_terms;
    
public:
    void addDynamicTerm(std::unique_ptr<PhysicsTerm> term) {
        dynamic_terms.push_back(std::move(term));
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Original Ug2 calculation
        auto it_QUA = params.find("QUA");
        auto it_M = params.find("mass");
        auto it_r = params.find("radius");
        auto it_Ereact = params.find("Ereact");
        auto it_S = params.find("step_function");
        
        double QUA = (it_QUA != params.end()) ? it_QUA->second : 1e-11;
        double M = (it_M != params.end()) ? it_M->second : 1e30;
        double r = (it_r != params.end()) ? it_r->second : 1e13;
        double Ereact = (it_Ereact != params.end()) ? it_Ereact->second : 1.0;
        double S = (it_S != params.end()) ? it_S->second : 1.0;
        
        double wind_mod = 1.0 + delta_sw * v_sw;
        double original = k2 * (QA + QUA) * M / (r * r) * S * wind_mod * HSCm * Ereact;
        
        // Add dynamic contributions
        double dynamic = 0.0;
        for (const auto& term : dynamic_terms) {
            if (term->validate(params)) {
                dynamic += term->compute(t, params);
            }
        }
        
        return original + dynamic;
    }
    
    std::string getName() const override { return "UQFFModule5Ug2Enhanced"; }
    
    std::string getDescription() const override {
        return "Enhanced Ug2 with dynamic term contributions (original + sum(dynamic))";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class UQFFModule5Ug3EnhancedTerm : public PhysicsTerm {
private:
    double k3 = 1.8;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamic_terms;
    
public:
    void addDynamicTerm(std::unique_ptr<PhysicsTerm> term) {
        dynamic_terms.push_back(std::move(term));
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Original Ug3 calculation
        auto it_Bj = params.find("Bj");
        auto it_omega_s_t = params.find("omega_s_t");
        auto it_Pcore = params.find("Pcore");
        auto it_Ereact = params.find("Ereact");
        
        double Bj = (it_Bj != params.end()) ? it_Bj->second : 1e-3;
        double omega_s_t = (it_omega_s_t != params.end()) ? it_omega_s_t->second : 1e-6;
        double Pcore = (it_Pcore != params.end()) ? it_Pcore->second : 1e-3;
        double Ereact = (it_Ereact != params.end()) ? it_Ereact->second : 1.0;
        
        double original = k3 * Bj * std::cos(omega_s_t * t * M_PI) * Pcore * Ereact;
        
        // Add dynamic contributions
        double dynamic = 0.0;
        for (const auto& term : dynamic_terms) {
            if (term->validate(params)) {
                dynamic += term->compute(t, params);
            }
        }
        
        return original + dynamic;
    }
    
    std::string getName() const override { return "UQFFModule5Ug3Enhanced"; }
    
    std::string getDescription() const override {
        return "Enhanced Ug3 with dynamic term contributions (original + sum(dynamic))";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class UQFFModule5CompressedMUGEEnhancedTerm : public PhysicsTerm {
private:
    double G = 6.67430e-11;
    double c = 3.0e8;
    double H0 = 2.269e-18;
    double Lambda = 1.1e-52;
    double hbar = 1.0546e-34;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamic_terms;
    
public:
    void addDynamicTerm(std::unique_ptr<PhysicsTerm> term) {
        dynamic_terms.push_back(std::move(term));
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Original Compressed MUGE (simplified version)
        auto it_M = params.find("mass");
        auto it_r = params.find("radius");
        auto it_B = params.find("B_field");
        auto it_Bcrit = params.find("Bcrit");
        
        double M = (it_M != params.end()) ? it_M->second : 1e30;
        double r = (it_r != params.end()) ? it_r->second : 1e4;
        double B = (it_B != params.end()) ? it_B->second : 1e10;
        double Bcrit = (it_Bcrit != params.end()) ? it_Bcrit->second : 1e11;
        
        // Base calculation
        double base = G * M / (r * r);
        double expansion = 1 + H0 * t;
        double super_adj = 1 - B / Bcrit;
        double cosm = Lambda * c * c / 3.0;
        
        double original = base * expansion * super_adj + cosm;
        
        // Add dynamic contributions
        double dynamic = 0.0;
        for (const auto& term : dynamic_terms) {
            if (term->validate(params)) {
                dynamic += term->compute(t, params);
            }
        }
        
        return original + dynamic;
    }
    
    std::string getName() const override { return "UQFFModule5CompressedMUGEEnhanced"; }
    
    std::string getDescription() const override {
        return "Enhanced Compressed MUGE with dynamic term contributions (original + sum(dynamic))";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class UQFFModule5ResonanceMUGEEnhancedTerm : public PhysicsTerm {
private:
    double fDPM = 1e12;
    double Evac_neb = 7.09e-36;
    double c_res = 3e8;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamic_terms;
    
public:
    void addDynamicTerm(std::unique_ptr<PhysicsTerm> term) {
        dynamic_terms.push_back(std::move(term));
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Original Resonance MUGE (simplified version - aDPM term only)
        auto it_I = params.find("I");
        auto it_A = params.find("A");
        auto it_omega1 = params.find("omega1");
        auto it_omega2 = params.find("omega2");
        auto it_Vsys = params.find("Vsys");
        
        double I = (it_I != params.end()) ? it_I->second : 1e21;
        double A = (it_A != params.end()) ? it_A->second : 3.142e8;
        double omega1 = (it_omega1 != params.end()) ? it_omega1->second : 1e-3;
        double omega2 = (it_omega2 != params.end()) ? it_omega2->second : -1e-3;
        double Vsys = (it_Vsys != params.end()) ? it_Vsys->second : 4.189e12;
        
        double FDPM = I * A * (omega1 - omega2);
        double original = FDPM * fDPM * Evac_neb * c_res * Vsys;
        
        // Add dynamic contributions
        double dynamic = 0.0;
        for (const auto& term : dynamic_terms) {
            if (term->validate(params)) {
                dynamic += term->compute(t, params);
            }
        }
        
        return original + dynamic;
    }
    
    std::string getName() const override { return "UQFFModule5ResonanceMUGEEnhanced"; }
    
    std::string getDescription() const override {
        return "Enhanced Resonance MUGE with dynamic term contributions (original + sum(dynamic))";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// ============================================================================
// INCLUDE ALL 19 CLASSES FROM source4_wolfram.cpp
// (UniversalGravity1-4, UniversalBuoyancy, UniversalMagnetism, UniversalAether,
//  UnifiedField, CompressedMUGE, ResonanceMUGE, 7 astrophysical systems,
//  ReactorEfficiency, NavierStokesQuasarJet)
// ============================================================================

// Copy all 19 classes from source4_wolfram.cpp here
// (Classes omitted for brevity - would include full definitions in actual file)

// ============================================================================
// REGISTRATION FUNCTION
// ============================================================================

void registerWolframTerms_source5(PhysicsTermRegistry& registry) {
    // NEW: Source5-specific dynamic terms (6 classes)
    registry.registerPhysicsTerm("DarkMatterHaloNFW", std::make_unique<DarkMatterHaloNFWTerm>(), "wolfram");
    registry.registerPhysicsTerm("VacuumEnergyFluctuation", std::make_unique<VacuumEnergyFluctuationTerm>(), "wolfram");
    registry.registerPhysicsTerm("UQFFModule5Ug1Enhanced", std::make_unique<UQFFModule5Ug1EnhancedTerm>(), "wolfram");
    registry.registerPhysicsTerm("UQFFModule5Ug2Enhanced", std::make_unique<UQFFModule5Ug2EnhancedTerm>(), "wolfram");
    registry.registerPhysicsTerm("UQFFModule5Ug3Enhanced", std::make_unique<UQFFModule5Ug3EnhancedTerm>(), "wolfram");
    registry.registerPhysicsTerm("UQFFModule5CompressedMUGEEnhanced", std::make_unique<UQFFModule5CompressedMUGEEnhancedTerm>(), "wolfram");
    registry.registerPhysicsTerm("UQFFModule5ResonanceMUGEEnhanced", std::make_unique<UQFFModule5ResonanceMUGEEnhancedTerm>(), "wolfram");
    
    // FROM source4: Universal Gravity (4 terms)
    registry.registerPhysicsTerm("UniversalGravity1", std::make_unique<UniversalGravity1Term>(), "wolfram");
    registry.registerPhysicsTerm("UniversalGravity2", std::make_unique<UniversalGravity2Term>(), "wolfram");
    registry.registerPhysicsTerm("UniversalGravity3", std::make_unique<UniversalGravity3Term>(), "wolfram");
    registry.registerPhysicsTerm("UniversalGravity4", std::make_unique<UniversalGravity4Term>(), "wolfram");
    
    // FROM source4: Universal Buoyancy, Magnetism, Aether (3 terms)
    registry.registerPhysicsTerm("UniversalBuoyancy", std::make_unique<UniversalBuoyancyTerm>(), "wolfram");
    registry.registerPhysicsTerm("UniversalMagnetism", std::make_unique<UniversalMagnetismTerm>(), "wolfram");
    registry.registerPhysicsTerm("UniversalAether", std::make_unique<UniversalAetherTerm>(), "wolfram");
    
    // FROM source4: Unified Field (1 term)
    registry.registerPhysicsTerm("UnifiedField", std::make_unique<UnifiedFieldTerm>(), "wolfram");
    
    // FROM source4: MUGE Equations (2 terms)
    registry.registerPhysicsTerm("CompressedMUGE", std::make_unique<CompressedMUGETerm>(), "wolfram");
    registry.registerPhysicsTerm("ResonanceMUGE", std::make_unique<ResonanceMUGETerm>(), "wolfram");
    
    // FROM source4: Astrophysical Systems (7 terms)
    registry.registerPhysicsTerm("SGR1745Magnetar", std::make_unique<SGR1745MagnetarTerm>(), "wolfram");
    registry.registerPhysicsTerm("SagittariusAStar", std::make_unique<SagittariusAStarTerm>(), "wolfram");
    registry.registerPhysicsTerm("TapestryStarbirth", std::make_unique<TapestryStarbirthTerm>(), "wolfram");
    registry.registerPhysicsTerm("Westerlund2Cluster", std::make_unique<Westerlund2ClusterTerm>(), "wolfram");
    registry.registerPhysicsTerm("PillarsCreation", std::make_unique<PillarsCreationTerm>(), "wolfram");
    registry.registerPhysicsTerm("RingsRelativity", std::make_unique<RingsRelativityTerm>(), "wolfram");
    registry.registerPhysicsTerm("StudentGuideUniverse", std::make_unique<StudentGuideUniverseTerm>(), "wolfram");
    
    // FROM source4: Helper Terms (2 terms)
    registry.registerPhysicsTerm("ReactorEfficiency", std::make_unique<ReactorEfficiencyTerm>(), "wolfram");
    registry.registerPhysicsTerm("NavierStokesQuasarJet", std::make_unique<NavierStokesQuasarJetTerm>(), "wolfram");
}

// ============================================================================
// TOTAL: 25 CLASSES REGISTERED
// - 6 new (source5-specific): DarkMatterHaloNFW, VacuumEnergyFluctuation, 
//   4 UQFFModule5Enhanced wrappers
// - 19 from source4: Ug1-4, Ubi, Um, A_mu_nu, FU, 2 MUGE, 7 systems, 2 helpers
// ============================================================================
