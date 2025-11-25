// Wolfram-Enhanced Physics Terms from source7.cpp
// Generated: November 24, 2025
// Source: YAML Support + Archive System + Same UQFF Physics as source6
// Total Classes: 41 (39 from source6 + 2 new: YAML loader + Archive manager)

#include <cmath>
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <filesystem>

// ============================================================================
// NEW SOURCE7-SPECIFIC CLASSES (2 ADDITIONS)
// ============================================================================

class YAMLConfigLoaderTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Returns YAML file complexity (number of key-value pairs loaded)
        auto it_keys = params.find("num_keys");
        auto it_nested = params.find("nested_levels");
        
        double keys = (it_keys != params.end()) ? it_keys->second : 10.0;
        double nested = (it_nested != params.end()) ? it_nested->second : 2.0;
        
        return keys * (1.0 + nested * 0.5); // complexity increases with nesting
    }
    
    std::string getName() const override { return "YAMLConfigLoader"; }
    
    std::string getDescription() const override {
        return "YAML config complexity: complexity = num_keys * (1 + nested_levels*0.5)";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

class ArchiveMediaManagerTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Returns total archived media size in MB
        auto it_img = params.find("image_size_mb");
        auto it_vid = params.find("video_size_mb");
        auto it_plugin = params.find("plugin_size_mb");
        
        double img = (it_img != params.end()) ? it_img->second : 0.0;
        double vid = (it_vid != params.end()) ? it_vid->second : 0.0;
        double plugin = (it_plugin != params.end()) ? it_plugin->second : 0.0;
        
        return img + vid + plugin; // total MB
    }
    
    std::string getName() const override { return "ArchiveMediaManager"; }
    
    std::string getDescription() const override {
        return "Archive media total size: total_MB = image + video + plugin sizes";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// ============================================================================
// INCLUDE ALL 39 CLASSES FROM source6_wolfram.cpp
// (14 graphics + 25 physics from source5/4)
// ============================================================================

// OpenGLRender, VulkanRender, MeshLoaderOBJ, ProceduralLandscape, MeshExtrude
// MeshBoolean, TextureLoader, ShaderCompile, CameraViewMatrix, BoneAnimation
// LaTeXRender, MultiViewport, SimulationEntityUpdate, ToolPathExecution
// DarkMatterHaloNFW, VacuumEnergyFluctuation, 4 UQFFModule5Enhanced wrappers
// UniversalGravity1-4, Buoyancy, Magnetism, Aether, UnifiedField
// CompressedMUGE, ResonanceMUGE, 7 systems, ReactorEfficiency, NavierStokes
// (Full class definitions would be copied here - omitted for brevity)

// ============================================================================
// REGISTRATION FUNCTION
// ============================================================================

void registerWolframTerms_source7(PhysicsTermRegistry& registry) {
    // NEW: Source7-specific additions (2 classes)
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
// TOTAL: 41 CLASSES REGISTERED
// - 2 new (source7-specific): YAMLConfigLoader, ArchiveMediaManager
// - 14 from source6 (graphics): OpenGL, Vulkan, Mesh, Landscape, etc.
// - 7 from source5 (dynamic): DarkMatterHalo, VacuumEnergy, Enhanced wrappers
// - 18 from source4 (core physics): Ug1-4, systems, MUGE, helpers
// ============================================================================
