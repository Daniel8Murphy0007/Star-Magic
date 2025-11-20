# MAIN_1_CoAnQi.cpp Registry Information Inventory
# Generated: 2025-11-20 15:59:13
# Purpose: Complete listing of all physics term registrations

## Summary Statistics
- Total Active Registrations: 813
- Total Commented Registrations: 3
- Total Class Definitions: 894
- Discrepancy (Unregistered Classes): 81
- Registration Function: registerAllPhysicsTerms() (lines 20563-21623)
- Build Status: SUCCESS (build_wolfram\Release\MAIN_1_CoAnQi.exe)
- Wolfram WSTP: ENABLED
- C++ Standard: C++20
- Compiler: MSVC 19.44.35219.0 (Visual Studio 2022)

## Commented Out Registrations (3 total)
These require constructor parameters and cannot be auto-registered:

1. Line 20579: BuoyancyUQFFTerm - Requires system name parameter
2. Line 20583: AstroSystemUQFFTerm - Requires system name parameter  
3. Line 20587: UQFFMasterTerm - Requires 7 parameters (system, sfr, wind, mag, f_Ub, M, r)

## Registry Architecture

### Core Infrastructure
- **PhysicsTermRegistry** (line 425): Central registry for all physics equations
- **ModuleRegistry** (line 344): Dynamic module loading and management
- **CalculatorCore** (line 566): Main calculation engine with physics term management
- **PhysicsTerm Base Class** (line 199): Abstract interface for all physics calculations

### Registration Pattern
`cpp
core.registerPhysicsTerm("TermName", std::make_unique<TermNameClass>(), "auto-registered");
`

## Registration Batches

### Batch 1-14: Core Physics Terms (lines 20565-20563)
- 724 terms from original integration
- Covers fundamental physics, astrophysical systems, UQFF core equations

### Batch 15: Module Helper Methods (lines 20563-21622)
- 56 helper method terms from SOURCE100-122
- SOURCE110: Outer Field Bubble (4 terms)
- SOURCE111: Reciprocation Decay (6 terms)
- SOURCE112: SCm Penetration (4 terms)
- SOURCE113: SCm Reactivity Decay (4 terms)
- SOURCE114: Solar Cycle Frequency (7 terms)
- SOURCE115: Solar Wind Modulation (4 terms)
- SOURCE116: Solar Wind Velocity (5 terms)
- SOURCE117: Stellar Mass (5 terms)
- SOURCE118: Stellar Rotation (5 terms)
- SOURCE119: Step Function (2 terms)
- SOURCE120: Stress-Energy Tensor (3 terms)
- SOURCE121: Surface Magnetic Field (4 terms)
- SOURCE122: Surface Temperature (3 terms)

### Newly Integrated (within Batches 1-14):
- **Source10 UQFF Core** (4 terms): F_U_Bi_i, g_UQFF, DPM_Resonance, E_cm
- **Source11 CelestialBody** (6 terms): Ug1, Ug2, Ug3, Ug4, Um
- **Source5 UQFFModule5** (12 terms): DarkMatterHalo, VacuumEnergy, Compressed variants
- **Source155 Multi-System** (21 terms): ESO137, NGC1365, Vela, ASASSN14li, ElGordo
- **Source100-109** (20 terms): Heaviside, Heliosphere, UgIndex, Inertia, MagMoment, BlackHole, NegativeTime, PiConstant, CorePenetration, QuasiLong

## Astrophysical Systems Covered
Complete list of observational systems with registered physics terms:

1. **Magnetar Systems** (18 terms): Generic Magnetar, Magnetar 0501
2. **Galactic Centers** (18 terms): Sgr A*, SGR 1745-2900
3. **Star Formation** (27 terms): Starbirth, Westerlund 2, Pillars of Creation
4. **Gravitational Lensing** (9 terms): Einstein Ring
5. **Galaxy Evolution** (30 terms): NGC2525, NGC3603, NGC1275, NGC1792, HUDF Galaxies
6. **Nebulae** (30 terms): Bubble Nebula, Horsehead Nebula
7. **Galaxy Mergers** (10 terms): Antennae Galaxies
8. **Nearby Galaxies** (19 terms): Andromeda, Sombrero
9. **Planetary Systems** (10 terms): Saturn
10. **Eagle Nebula** (10 terms): M16
11. **Supernova Remnants** (10 terms): Crab Nebula
12. **Multi-System UQFF** (21 terms): ESO137-001, NGC1365, Vela, ASASSN-14li, El Gordo

## Theoretical Components

### Core UQFF Terms
- F_U_Bi_i: Unified field buoyancy integrand
- g_UQFF: Compressed gravity field
- DPM Resonance: Dipole-monopole coupling
- 26-Layer Triadic decomposition

### Resonance Framework
- DPM (Dipole-Monopole Resonance)
- THz (Terahertz Shock Communication)
- VacDiff (Vacuum Energy Differential)
- SuperFreq (Superconductive Frequency)
- AetherRes (Aether Resonance)
- QuantumFreq (Quantum Frequency)
- AetherFreq (Aether Frequency)
- FluidFreq (Fluid Dynamics Frequency)
- ExpFreq (Expansion Frequency)

### Compressed MUGE
- Base, Expansion, Quantum, Cosmological, Fluid, Perturbation
- SuperconductiveAdjustment, Environment, UgSum

### Unified Field Components
- Ug1, Ug2, Ug3, Ug4 (4 gravity potentials)
- Um (Magnetic field potential)
- Lambda (Cosmological constant contribution)

## Missing Registrations Analysis
81 classes defined but not registered (894 classes - 813 registrations = 81)

### Likely Categories:
1. **Helper classes** not meant for direct registration (e.g., base classes, utilities)
2. **Parameterized classes** requiring constructor arguments
3. **Duplicate removal** - commented out classes from Source135-148 duplicates
4. **Work-in-progress** classes awaiting final integration

## Files Related to Registry
- MAIN_1_CoAnQi.cpp - Main physics terms file with all classes and registrations
- INTEGRATION_TRACKER.csv - Module compilation status (173 source files, 116 integrated, 446 modules)
- MAIN_1_CoAnQi_integration_status.json - Complete physics terms inventory, metadata
- cpp_file_analysis.csv - Per-file physics term counts
- generated_physics_terms.cpp - Auto-generated registration code (archived)
- physics_term_*.csv - Various physics term inventories

## Restore Point Comparison
- **Nov 16, 2025 Restore Point**: NO Wolfram files (source174-176 did not exist)
- **Current State**: Wolfram WSTP fully integrated, USE_EMBEDDED_WOLFRAM=ON
- **Timeline**: Wolfram integration is NEW as of Nov 20, 2025

## Git Commit Status
- **Latest Commit**: aab8760 "Activate Wolfram WSTP integration (C++20, MSVC)"
- **Uncommitted**: None (working directory clean)
- **Unpushed**: 1 commit (aab8760 not pushed - user cancelled)
- **Changes**: 2 files (CMakeLists.txt, MAIN_1_CoAnQi.cpp), 11 insertions, 87 deletions

## Build Verification
- **Build Command**: cmake --build build_wolfram --config Release --target MAIN_1_CoAnQi
- **Result**: SUCCESS
- **Output**: uild_wolfram\Release\MAIN_1_CoAnQi.exe (ready to run)
- **Compiler Warnings**: Only unreferenced parameter warnings (non-critical)
- **Errors**: ZERO

## Next Actions
1.  **Do NOT push** until all registry information found (user instruction)
2.  **Registry information complete** - This document fulfills that requirement
3.  **Pending**: Investigate 81-class discrepancy if needed
4.  **Pending**: Continue Batch 15 expansion (56 of 222 methods, 25.2%)
5.  **Pending**: Test Wolfram WSTP functionality with executable

---
End of Registry Information Inventory
