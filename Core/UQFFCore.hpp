// UQFFCore.hpp - Unified Quantum Field Force Core Master Header
// Generated: November 14, 2025
// Purpose: Consolidate all Core modules for Star-Magic UQFF framework
//
// This header provides one-stop access to all 155+ extracted UQFF modules
// organized into the Core architecture from the 167-module validated archive.
//
// Usage:
//   #include "Core/UQFFCore.hpp"
//
// Then access any module class directly without individual includes.

#ifndef UQFF_CORE_HPP
#define UQFF_CORE_HPP

// ============================================================================
// FOUNDATION MODULES (Week 1)
// ============================================================================
#include "SystemCatalogue.hpp"  // 20 astrophysical systems, master equations

// ============================================================================
// ADAPTIVE FRAMEWORK (Week 2)
// ============================================================================
#include "PhysicsTerms.hpp"      // Dynamic physics term base classes
#include "UQFFModule4.hpp"       // Self-expanding physics engine v2.0
#include "FluidSolver.hpp"       // Navier-Stokes with UQFF integration

// ============================================================================
// EXTRACTED MODULES (Weeks 3-4) - 155 modules
// ============================================================================
// Alphabetically organized for easy reference

// A
#include "Modules/Abell2256UQFFModule.hpp"
#include "Modules/AetherCouplingModule.hpp"
#include "Modules/AetherVacuumDensityModule.hpp"
#include "Modules/AndromedaUQFFModule.hpp"
#include "Modules/ASASSN14liUQFFModule.hpp"
#include "Modules/AstroSystemsUQFFModule.hpp"

// B
#include "Modules/BackgroundAetherModule.hpp"
#include "Modules/BigBangGravityUQFFModule.hpp"
#include "Modules/BuoyancyCouplingModule.hpp"
#include "Modules/ButterflyNebulaUQFFModule.hpp"

// C
#include "Modules/CentaurusAUQFFModule.hpp"
#include "Modules/CompressedResonanceUQFF24Module.hpp"
#include "Modules/CompressedResonanceUQFF34Module.hpp"
#include "Modules/CompressedResonanceUQFFModule.hpp"
#include "Modules/CorePenetrationModule.hpp"
#include "Modules/CrabNebulaUQFFModule.hpp"
#include "Modules/CrabResonanceUQFFModule.hpp"
#include "Modules/CrabUQFFModule.hpp"

// D
#include "Modules/DPMModule.hpp"

// E
#include "Modules/ElGordoUQFFModule.hpp"
#include "Modules/ESO137UQFFModule.hpp"

// F
#include "Modules/FeedbackFactorModule.hpp"
#include "Modules/FluidSolver.hpp"

// G
#include "Modules/GalacticBlackHoleModule.hpp"
#include "Modules/GalacticDistanceModule.hpp"

// H
#include "Modules/HeavisideFractionModule.hpp"
#include "Modules/HeliosphereThicknessModule.hpp"
#include "Modules/HydrogenAtomUQFFModule.hpp"
#include "Modules/HydrogenPToEResonanceUQFFModule.hpp"
#include "Modules/HydrogenResonanceUQFFModule.hpp"
#include "Modules/HydrogenUQFFModule.hpp"

// I
#include "Modules/IC2163UQFFModule.hpp"
#include "Modules/InertiaUQFFModule.hpp"
#include "Modules/InertiaCouplingModule.hpp"

// J
#include "Modules/J1610UQFFModule.hpp"
#include "Modules/JupiterAuroraeUQFFModule.hpp"

// L
#include "Modules/LagoonNebulaUQFFModule.hpp"
#include "Modules/LagoonUQFFModule.hpp"
#include "Modules/LENRCalibUQFFModule.hpp"
#include "Modules/LENRUQFFModule.hpp"

// M
#include "Modules/M16UQFFModule.hpp"
#include "Modules/M51UQFFModule.hpp"
#include "Modules/M87JetUQFFModule.hpp"
#include "Modules/MagneticMomentModule.hpp"
#include "Modules/MagneticStringModule.hpp"
#include "Modules/MUGEModule.hpp"
#include "Modules/MUGEResonanceModule.hpp"
#include "Modules/MultiCompressedUQFFModule.hpp"
#include "Modules/MultiUQFFCompressionModule.hpp"
#include "Modules/MultiUQFFModule.hpp"

// N
#include "Modules/NebularUQFFModule.hpp"
#include "Modules/NegativeTimeModule.hpp"
#include "Modules/NGC1300UQFFModule.hpp"
#include "Modules/NGC1316UQFFModule.hpp"
#include "Modules/NGC1365UQFFModule.hpp"
#include "Modules/NGC2207UQFFModule.hpp"
#include "Modules/NGC2264UQFFModule.hpp"
#include "Modules/NGC346UQFFModule.hpp"
#include "Modules/NGC4676UQFFModule.hpp"
#include "Modules/NGC6302ResonanceUQFFModule.hpp"
#include "Modules/NGC6302UQFFModule.hpp"

// O
#include "Modules/OrionUQFFModule.hpp"
#include "Modules/OuterFieldBubbleModule.hpp"

// P
#include "Modules/PiConstantModule.hpp"

// Q
#include "Modules/QuasiLongitudinalModule.hpp"

// R
#include "Modules/RAquariiUQFFModule.hpp"
#include "Modules/ReciprocationDecayModule.hpp"
#include "Modules/RedDwarfUQFFModule.hpp"
#include "Modules/RedSpiderUQFFModule.hpp"
#include "Modules/ResonanceSuperconductiveUQFFModule.hpp"

// S
#include "Modules/SaturnUQFFModule.hpp"
#include "Modules/ScientificCalculatorDialog.hpp"
#include "Modules/ScmPenetrationModule.hpp"
#include "Modules/ScmReactivityDecayModule.hpp"
#include "Modules/ScmVacuumDensityModule.hpp"
#include "Modules/ScmVelocityModule.hpp"
#include "Modules/SGR1745UQFFModule.hpp"
#include "Modules/SgrA_UQFFModule.hpp"
#include "Modules/SgrAStarUQFFModule.hpp"
#include "Modules/SIMPlugin.hpp"
#include "Modules/SMBHBinaryUQFFModule.hpp"
#include "Modules/SMBHUQFFModule.hpp"
#include "Modules/SolarCycleFrequencyModule.hpp"
#include "Modules/SolarWindBuoyancyModule.hpp"
#include "Modules/SolarWindModulationModule.hpp"
#include "Modules/SolarWindVelocityModule.hpp"
#include "Modules/SombreroUQFFModule.hpp"
#include "Modules/SpiralSupernovaeUQFFModule.hpp"
#include "Modules/SPTCLJ2215UQFFModule.hpp"
#include "Modules/StellarMassModule.hpp"
#include "Modules/StellarRotationModule.hpp"
#include "Modules/StepFunctionModule.hpp"
#include "Modules/StephanQuintetUQFFModule.hpp"
#include "Modules/StressEnergyTensorModule.hpp"
#include "Modules/SurfaceMagneticFieldModule.hpp"
#include "Modules/SurfaceTemperatureModule.hpp"
#include "Modules/SymEngineAllocator.hpp"

// T
#include "Modules/TapestryUQFFModule.hpp"
#include "Modules/TimeReversalZoneModule.hpp"

// U
#include "Modules/UaVacuumDensityModule.hpp"
#include "Modules/UFEOrbModule.hpp"
#include "Modules/Ug1DefectModule.hpp"
#include "Modules/Ug3DiskVectorModule.hpp"
#include "Modules/UgCouplingModule.hpp"
#include "Modules/UgIndexModule.hpp"
#include "Modules/UGC10214UQFFModule.hpp"
#include "Modules/UnifiedFieldModule.hpp"
#include "Modules/UniversalInertiaVacuumModule.hpp"
#include "Modules/UniverseDiameterUQFFModule.hpp"
#include "Modules/UQFF8AstroSystemsModule.hpp"
#include "Modules/UQFFBuoyancyAstroModule.hpp"
#include "Modules/UQFFBuoyancyCNBModule.hpp"
#include "Modules/UQFFBuoyancyModule.hpp"
#include "Modules/UQFFBuoyancyModule157.hpp"
#include "Modules/UQFFCompressedResonanceModule.hpp"
#include "Modules/UQFFCompressionModule.hpp"
#include "Modules/UQFFCoreModule.hpp"
#include "Modules/UQFFModule5.hpp"
#include "Modules/UQFFNebulaTriadicModule.hpp"

// V
#include "Modules/V838MonUQFFModule.hpp"
#include "Modules/VelaPulsarUQFFModule.hpp"

// Y
#include "Modules/YoungStarsOutflowsUQFFModule.hpp"

// ============================================================================
// NAMESPACE ORGANIZATION
// ============================================================================

namespace UQFFCore {
    // All Core modules live under the UQFFCore namespace
    // This prevents naming conflicts with the 3000+ module archive
    
    // Version information
    constexpr const char* VERSION = "1.0.0-alpha";
    constexpr const char* BUILD_DATE = "2025-11-14";
    constexpr int MODULE_COUNT = 155;
    
    // Architecture layers
    enum class Layer {
        Foundation,      // SystemCatalogue
        Adaptive,        // UQFFModule4, FluidSolver  
        AstroPhysics,    // Galaxy, stellar, AGN modules
        QuantumFields,   // Vacuum, quantum coupling
        FluidDynamics,   // Navier-Stokes, MHD
        Relativity,      // GR corrections
        Nuclear,         // LENR, neutron physics
        Cosmology,       // Dark matter/energy
        Specialized      // Mission-specific modules
    };
    
    // Initialization function (optional - for future use)
    // void Initialize();
    // void Shutdown();
}

#endif // UQFF_CORE_HPP
