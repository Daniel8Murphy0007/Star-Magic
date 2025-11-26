// Wolfram-Enhanced Physics Terms from source8.cpp
// Generated: November 24, 2025
// Source: MEGA SCIENTIFIC CALCULATOR UI (100+ libraries, infrastructure only)
// Total Classes: 51 (10 new computational infrastructure + 41 from source7)

#include <cmath>
#include <string>
#include <map>
#include <memory>
#include <vector>

// ============================================================================
// NEW SOURCE8-SPECIFIC CLASSES (10 COMPUTATIONAL INFRASTRUCTURE ADDITIONS)
// ============================================================================

class DimensionalAnalysisTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Returns unit consistency score: 1.0 if all units match, 0.0 if mismatch
        auto it_mass1 = params.find("mass_exp1");
        auto it_mass2 = params.find("mass_exp2");
        auto it_length1 = params.find("length_exp1");
        auto it_length2 = params.find("length_exp2");
        auto it_time1 = params.find("time_exp1");
        auto it_time2 = params.find("time_exp2");
        
        if (it_mass1 == params.end() || it_mass2 == params.end() ||
            it_length1 == params.end() || it_length2 == params.end() ||
            it_time1 == params.end() || it_time2 == params.end()) {
            return 0.0; // Missing unit data
        }
        
        bool match = (std::abs(it_mass1->second - it_mass2->second) < 1e-9 &&
                      std::abs(it_length1->second - it_length2->second) < 1e-9 &&
                      std::abs(it_time1->second - it_time2->second) < 1e-9);
        
        return match ? 1.0 : 0.0;
    }
    
    std::string getName() const override { return "DimensionalAnalysis"; }
    
    std::string getDescription() const override {
        return "Unit consistency check: M^a L^b T^c matching for dimensional analysis";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class QAOAOptimizationTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Returns QAOA expectation value: <C> = Σ H_i
        auto it_layers = params.find("qaoa_layers");
        auto it_beta = params.find("beta_angle");
        auto it_gamma = params.find("gamma_angle");
        
        int layers = (it_layers != params.end()) ? static_cast<int>(it_layers->second) : 1;
        double beta = (it_beta != params.end()) ? it_beta->second : 0.5;
        double gamma = (it_gamma != params.end()) ? it_gamma->second : 1.0;
        
        // Simplified: <C> = layers * cos(beta) * sin(gamma)
        return layers * std::cos(beta) * std::sin(gamma);
    }
    
    std::string getName() const override { return "QAOAOptimization"; }
    
    std::string getDescription() const override {
        return "QAOA quantum optimization: <C> = layers * cos(beta) * sin(gamma)";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class CategoryFunctorTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Returns functor transformation score (category theory mapping)
        auto it_objects = params.find("category_objects");
        auto it_morphisms = params.find("category_morphisms");
        
        double objects = (it_objects != params.end()) ? it_objects->second : 1.0;
        double morphisms = (it_morphisms != params.end()) ? it_morphisms->second : 0.0;
        
        return objects + morphisms * 0.5; // Complexity of functor
    }
    
    std::string getName() const override { return "CategoryFunctor"; }
    
    std::string getDescription() const override {
        return "Category theory functor: complexity = objects + morphisms*0.5";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class LLVMJITCompilerTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Returns JIT compilation speedup factor
        auto it_opcodes = params.find("llvm_opcodes");
        auto it_opt_level = params.find("optimization_level");
        
        double opcodes = (it_opcodes != params.end()) ? it_opcodes->second : 100.0;
        double opt = (it_opt_level != params.end()) ? it_opt_level->second : 2.0; // O0-O3
        
        // Speedup = (1 + opt*0.3) * sqrt(opcodes) (more opcodes = better JIT payoff)
        return (1.0 + opt * 0.3) * std::sqrt(opcodes);
    }
    
    std::string getName() const override { return "LLVMJITCompiler"; }
    
    std::string getDescription() const override {
        return "LLVM JIT speedup: (1 + opt_level*0.3) * sqrt(opcodes)";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class FederatedLearningTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Returns federated model accuracy after N rounds
        auto it_clients = params.find("fl_clients");
        auto it_rounds = params.find("fl_rounds");
        auto it_local_epochs = params.find("local_epochs");
        
        double clients = (it_clients != params.end()) ? it_clients->second : 10.0;
        double rounds = (it_rounds != params.end()) ? it_rounds->second : 5.0;
        double epochs = (it_local_epochs != params.end()) ? it_local_epochs->second : 3.0;
        
        // Accuracy = 0.5 + 0.4*(1 - exp(-rounds/10)) * sqrt(clients*epochs/100)
        return 0.5 + 0.4 * (1.0 - std::exp(-rounds / 10.0)) * std::sqrt(clients * epochs / 100.0);
    }
    
    std::string getName() const override { return "FederatedLearning"; }
    
    std::string getDescription() const override {
        return "Federated learning accuracy: 0.5 + 0.4*(1-exp(-rounds/10))*sqrt(clients*epochs/100)";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class NeuralSymbolicEvalTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Returns neural-symbolic hybrid prediction error
        auto it_sym_complexity = params.find("symbolic_complexity");
        auto it_nn_layers = params.find("neural_layers");
        auto it_data_size = params.find("training_samples");
        
        double sym = (it_sym_complexity != params.end()) ? it_sym_complexity->second : 10.0;
        double layers = (it_nn_layers != params.end()) ? it_nn_layers->second : 3.0;
        double data = (it_data_size != params.end()) ? it_data_size->second : 1000.0;
        
        // Error = sym/(layers*sqrt(data)) (more layers + data = lower error)
        return sym / (layers * std::sqrt(data));
    }
    
    std::string getName() const override { return "NeuralSymbolicEval"; }
    
    std::string getDescription() const override {
        return "Neural-symbolic hybrid error: symbolic_complexity / (layers * sqrt(data))";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class NeuromorphicAcceleratorTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Returns neuromorphic hardware speedup vs CPU
        auto it_neurons = params.find("neuromorphic_neurons");
        auto it_spikes_per_sec = params.find("spikes_per_second");
        
        double neurons = (it_neurons != params.end()) ? it_neurons->second : 1000.0;
        double spikes = (it_spikes_per_sec != params.end()) ? it_spikes_per_sec->second : 1e6;
        
        // Speedup = (neurons * spikes) / 1e9 (normalized)
        return (neurons * spikes) / 1e9;
    }
    
    std::string getName() const override { return "NeuromorphicAccelerator"; }
    
    std::string getDescription() const override {
        return "Neuromorphic speedup: (neurons * spikes_per_sec) / 1e9";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class BlockchainECDSATerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Returns ECDSA signature verification time (ms)
        auto it_signatures = params.find("ecdsa_signatures");
        auto it_curve_bits = params.find("secp256k1_bits");
        
        double sigs = (it_signatures != params.end()) ? it_signatures->second : 1.0;
        double bits = (it_curve_bits != params.end()) ? it_curve_bits->second : 256.0;
        
        // Time = sigs * (bits/128) * 0.5 ms baseline
        return sigs * (bits / 128.0) * 0.5;
    }
    
    std::string getName() const override { return "BlockchainECDSA"; }
    
    std::string getDescription() const override {
        return "ECDSA verification time: signatures * (curve_bits/128) * 0.5 ms";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class OperationalTransformTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Returns OT convergence time for collaborative editing (ms)
        auto it_clients = params.find("ot_clients");
        auto it_ops = params.find("ot_operations");
        auto it_latency = params.find("network_latency_ms");
        
        double clients = (it_clients != params.end()) ? it_clients->second : 2.0;
        double ops = (it_ops != params.end()) ? it_ops->second : 10.0;
        double latency = (it_latency != params.end()) ? it_latency->second : 50.0;
        
        // Convergence = clients * ops * latency / 100 (simplified)
        return clients * ops * latency / 100.0;
    }
    
    std::string getName() const override { return "OperationalTransform"; }
    
    std::string getDescription() const override {
        return "OT convergence time: clients * operations * latency / 100 ms";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class MPIDistributedTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Returns MPI parallel efficiency (0-1)
        auto it_procs = params.find("mpi_processes");
        auto it_comm_overhead = params.find("communication_overhead_pct");
        
        double procs = (it_procs != params.end()) ? it_procs->second : 1.0;
        double overhead = (it_comm_overhead != params.end()) ? it_comm_overhead->second : 10.0; // %
        
        // Efficiency = 1 / (1 + overhead/100) / log2(procs)
        return 1.0 / (1.0 + overhead / 100.0) / std::log2(std::max(2.0, procs));
    }
    
    std::string getName() const override { return "MPIDistributed"; }
    
    std::string getDescription() const override {
        return "MPI parallel efficiency: 1 / ((1 + overhead/100) * log2(procs))";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// ============================================================================
// INCLUDE ALL 41 CLASSES FROM source7_wolfram.cpp
// (2 infrastructure from source7 + 39 from source6/5/4)
// ============================================================================

// YAMLConfigLoader, ArchiveMediaManager (source7)
// OpenGLRender, VulkanRender, MeshLoaderOBJ, ... (14 from source6)
// DarkMatterHaloNFW, VacuumEnergy, 4 Enhanced (7 from source5)
// Ug1-4, systems, MUGE, helpers (18 from source4)
// (Full class definitions would be copied here - omitted for brevity)

// ============================================================================
// REGISTRATION FUNCTION
// ============================================================================

void registerWolframTerms_source8(PhysicsTermRegistry& registry) {
    // NEW: Source8-specific computational infrastructure (10 classes)
    registry.registerPhysicsTerm("DimensionalAnalysis", std::make_unique<DimensionalAnalysisTerm>(), "wolfram");
    registry.registerPhysicsTerm("QAOAOptimization", std::make_unique<QAOAOptimizationTerm>(), "wolfram");
    registry.registerPhysicsTerm("CategoryFunctor", std::make_unique<CategoryFunctorTerm>(), "wolfram");
    registry.registerPhysicsTerm("LLVMJITCompiler", std::make_unique<LLVMJITCompilerTerm>(), "wolfram");
    registry.registerPhysicsTerm("FederatedLearning", std::make_unique<FederatedLearningTerm>(), "wolfram");
    registry.registerPhysicsTerm("NeuralSymbolicEval", std::make_unique<NeuralSymbolicEvalTerm>(), "wolfram");
    registry.registerPhysicsTerm("NeuromorphicAccelerator", std::make_unique<NeuromorphicAcceleratorTerm>(), "wolfram");
    registry.registerPhysicsTerm("BlockchainECDSA", std::make_unique<BlockchainECDSATerm>(), "wolfram");
    registry.registerPhysicsTerm("OperationalTransform", std::make_unique<OperationalTransformTerm>(), "wolfram");
    registry.registerPhysicsTerm("MPIDistributed", std::make_unique<MPIDistributedTerm>(), "wolfram");
    
    // FROM source7: Infrastructure (2 classes)
    registry.registerPhysicsTerm("YAMLConfigLoader", std::make_unique<YAMLConfigLoaderTerm>(), "wolfram");
    registry.registerPhysicsTerm("ArchiveMediaManager", std::make_unique<ArchiveMediaManagerTerm>(), "wolfram");
    
    // FROM source6: Graphics infrastructure (14 classes)
    registry.registerPhysicsTerm("OpenGLRender", std::make_unique<OpenGLRenderTerm>(), "wolfram");
    registry.registerPhysicsTerm("VulkanRender", std::make_unique<VulkanRenderTerm>(), "wolfram");
    registry.registerPhysicsTerm("MeshLoaderOBJ", std::make_unique<MeshLoaderOBJTerm>(), "wolfram");
    registry.registerPhysicsTerm("ProceduralLandscape", std::make_unique<ProceduralLandscapeTerm>(), "wolfram");
    registry.registerPhysicsTerm("MeshExtrude", std::make_unique<MeshExtrudeTerm>(), "wolfram");
    registry.registerPhysicsTerm("MeshBoolean", std::make_unique<MeshBooleanTerm>(), "wolfram");
    registry.registerPhysicsTerm("TextureLoader", std::make_unique<TextureLoaderTerm>(), "wolfram");
    registry.registerPhysicsTerm("ShaderCompile", std::make_unique<ShaderCompileTerm>(), "wolfram");
    registry.registerPhysicsTerm("CameraViewMatrix", std::make_unique<CameraViewMatrixTerm>(), "wolfram");
    registry.registerPhysicsTerm("BoneAnimation", std::make_unique<BoneAnimationTerm>(), "wolfram");
    registry.registerPhysicsTerm("LaTeXRender", std::make_unique<LaTeXRenderTerm>(), "wolfram");
    registry.registerPhysicsTerm("MultiViewport", std::make_unique<MultiViewportTerm>(), "wolfram");
    registry.registerPhysicsTerm("SimulationEntityUpdate", std::make_unique<SimulationEntityUpdateTerm>(), "wolfram");
    registry.registerPhysicsTerm("ToolPathExecution", std::make_unique<ToolPathExecutionTerm>(), "wolfram");
    
    // FROM source5: Dynamic framework (7 classes)
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
// TOTAL: 51 CLASSES REGISTERED
// - 10 new (source8-specific): DimensionalAnalysis, QAOA, CategoryFunctor, LLVM, FL, NeuralSymbolic, Neuromorphic, Blockchain, OT, MPI
// - 2 from source7: YAMLConfigLoader, ArchiveMediaManager
// - 14 from source6: Graphics infrastructure
// - 7 from source5: Dynamic framework
// - 18 from source4: Core UQFF physics
// ============================================================================
