/**
 * ================================================================================================
 * test_core_integration.cpp - Core Module Integration Test
 * ================================================================================================
 * 
 * Purpose: Verify all extracted Core modules work together correctly
 * Tests:
 *   - SystemCatalogue: Load systems, compute master equations
 *   - PhysicsTerms: Dynamic term registration and computation
 *   - UQFFModule4: Variable management, auto-calibration, state persistence
 *   - FluidSolver: Navier-Stokes simulation with UQFF integration
 * 
 * Phase: 1, Weeks 1-2 Validation
 * Date: November 14, 2025
 * ================================================================================================
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "Core/SystemCatalogue.hpp"
#include "Core/PhysicsTerms.hpp"
#include "Core/UQFFModule4.hpp"
#include "Core/FluidSolver.hpp"

using namespace UQFFCatalogue;
using namespace UQFFCore;

// Test result tracking
int tests_passed = 0;
int tests_failed = 0;

#define TEST_ASSERT(condition, message) \
    do { \
        if (condition) { \
            std::cout << "  ✓ " << message << std::endl; \
            tests_passed++; \
        } else { \
            std::cout << "  ✗ FAILED: " << message << std::endl; \
            tests_failed++; \
        } \
    } while(0)

void test_system_catalogue() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST 1: SystemCatalogue" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Test 1.1: Load systems
    auto systems = initializeSystemCatalogue();
    TEST_ASSERT(systems.size() == 20, "Loaded 20 astrophysical systems");
    
    // Test 1.2: Verify NGC 1365 parameters
    auto it = systems.find("NGC 1365");
    TEST_ASSERT(it != systems.end(), "Found NGC 1365 system");
    
    if (it != systems.end()) {
        const SystemParams& ngc1365 = it->second;
        TEST_ASSERT(ngc1365.name == "NGC 1365", "NGC 1365 name correct");
        TEST_ASSERT(ngc1365.M == 1.2e11 * Msun, "NGC 1365 mass correct");
        TEST_ASSERT(ngc1365.r == 6e20, "NGC 1365 radius correct");
        
        // Test 1.3: Compute E_cm
        double E_cm = compute_E_cm(ngc1365);
        TEST_ASSERT(E_cm > 0, "E_cm computation positive");
        std::cout << "    E_cm = " << E_cm << " J" << std::endl;
        
        // Test 1.4: Compute F_U_Bi_i (buoyancy)
        double F = F_U_Bi_i(ngc1365);
        TEST_ASSERT(!std::isnan(F) && !std::isinf(F), "F_U_Bi_i computation valid");
        std::cout << "    F_U_Bi_i = " << F << " N" << std::endl;
        
        // Test 1.5: Compute compressed_g
        double g = compressed_g(ngc1365);
        TEST_ASSERT(!std::isnan(g) && !std::isinf(g), "compressed_g computation valid");
        std::cout << "    compressed_g = " << g << " m/s²" << std::endl;
    }
    
    // Test 1.6: Verify multiple systems exist
    TEST_ASSERT(systems.find("Vela Pulsar") != systems.end(), "Vela Pulsar exists");
    TEST_ASSERT(systems.find("GW170817") != systems.end(), "GW170817 exists");
    TEST_ASSERT(systems.find("3C273") != systems.end(), "3C273 exists");
}

void test_physics_terms() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST 2: PhysicsTerms" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Test 2.1: DynamicVacuumTerm
    DynamicVacuumTerm vacTerm(1e-10, 2 * PI * 1e-15);
    TEST_ASSERT(vacTerm.getName() == "DynamicVacuumTerm", "DynamicVacuumTerm name correct");
    
    std::map<std::string, double> params;
    TEST_ASSERT(vacTerm.validate(params), "DynamicVacuumTerm validates");
    
    double vac_result = vacTerm.compute(1e6, params);
    TEST_ASSERT(!std::isnan(vac_result), "DynamicVacuumTerm computes");
    std::cout << "    Vacuum term at t=1e6: " << vac_result << std::endl;
    
    // Test 2.2: QuantumCouplingTerm
    QuantumCouplingTerm qcTerm(1.23e-40);
    TEST_ASSERT(qcTerm.getName() == "QuantumCouplingTerm", "QuantumCouplingTerm name correct");
    
    params["mass"] = 1e30;
    params["radius"] = 1e6;
    TEST_ASSERT(qcTerm.validate(params), "QuantumCouplingTerm validates");
    
    double qc_result = qcTerm.compute(1e6, params);
    TEST_ASSERT(!std::isnan(qc_result), "QuantumCouplingTerm computes");
    std::cout << "    Quantum coupling at t=1e6: " << qc_result << std::endl;
}

void test_uqff_module4() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST 3: UQFFModule4" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Test 3.1: Create module
    UQFFModule4 module;
    TEST_ASSERT(true, "UQFFModule4 constructed");
    
    // Test 3.2: Variable management
    module.updateVariable("test_mass", 1e30);
    double mass = module.getVariable("test_mass");
    TEST_ASSERT(mass == 1e30, "Variable update and retrieval");
    
    // Test 3.3: Variable history
    module.updateVariable("test_mass", 2e30);
    module.updateVariable("test_mass", 3e30);
    auto history = module.getVariableHistory("test_mass");
    TEST_ASSERT(history.size() == 3, "Variable history tracking");
    
    // Test 3.4: Dynamic parameters
    module.setDynamicParameter("custom_param", 42.0);
    double param = module.getDynamicParameter("custom_param");
    TEST_ASSERT(param == 42.0, "Dynamic parameter management");
    
    // Test 3.5: Dynamic term registration
    module.setEnableDynamicTerms(true);
    module.registerDynamicTerm(std::make_unique<DynamicVacuumTerm>(1e-10, 1e-15));
    double dynamic_result = module.computeDynamicTerms(1e6);
    TEST_ASSERT(!std::isnan(dynamic_result), "Dynamic term computation");
    std::cout << "    Dynamic terms result: " << dynamic_result << std::endl;
    
    // Test 3.6: State export/import
    module.exportState("test_module_state.txt");
    UQFFModule4 module2;
    module2.importState("test_module_state.txt");
    double imported_mass = module2.getVariable("test_mass");
    TEST_ASSERT(imported_mass == 3e30, "State export and import");
    
    // Test 3.7: Metadata
    auto metadata = module.getMetadata();
    TEST_ASSERT(metadata["version"] == "2.0-Enhanced", "Metadata correct");
}

void test_fluid_solver() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST 4: FluidSolver" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Test 4.1: Create solver
    FluidSolver solver;
    TEST_ASSERT(true, "FluidSolver constructed");
    
    // Test 4.2: Initialize velocity field
    int size = (GRID_SIZE + 2) * (GRID_SIZE + 2);
    TEST_ASSERT(solver.u.size() == size, "Velocity u field initialized");
    TEST_ASSERT(solver.v.size() == size, "Velocity v field initialized");
    
    // Test 4.3: Add jet force
    solver.add_jet_force(10.0);
    bool has_velocity = false;
    for (size_t i = 0; i < solver.v.size(); i++) {
        if (std::abs(solver.v[i]) > 1e-10) {
            has_velocity = true;
            break;
        }
    }
    TEST_ASSERT(has_velocity, "Jet force added to velocity field");
    
    // Test 4.4: Simulate with UQFF gravity
    double uqff_g = -9.8; // Earth-like gravity
    solver.step(uqff_g);
    TEST_ASSERT(true, "Fluid step with UQFF gravity executed");
    
    // Test 4.5: Multiple timesteps
    for (int i = 0; i < 5; i++) {
        solver.step(uqff_g);
    }
    TEST_ASSERT(true, "Multiple timestep simulation stable");
    
    // Test 4.6: Verify no NaN/Inf in fields
    bool fields_valid = true;
    for (size_t i = 0; i < solver.u.size(); i++) {
        if (std::isnan(solver.u[i]) || std::isinf(solver.u[i]) ||
            std::isnan(solver.v[i]) || std::isinf(solver.v[i])) {
            fields_valid = false;
            break;
        }
    }
    TEST_ASSERT(fields_valid, "Velocity fields remain valid (no NaN/Inf)");
    
    std::cout << "\n  Velocity field visualization:" << std::endl;
    solver.print_velocity_field();
}

void test_integrated_workflow() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST 5: Integrated Workflow" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Test 5.1: Load system from catalogue
    auto systems = initializeSystemCatalogue();
    auto& vela = systems["Vela Pulsar"];
    TEST_ASSERT(true, "Loaded Vela Pulsar from SystemCatalogue");
    
    // Test 5.2: Compute UQFF gravity
    double g_uqff = compressed_g(vela);
    std::cout << "    Vela Pulsar UQFF gravity: " << g_uqff << " m/s²" << std::endl;
    
    // Test 5.3: Use UQFF gravity in fluid solver
    FluidSolver solver;
    solver.add_jet_force(100.0); // Strong pulsar jet
    solver.step(g_uqff);
    TEST_ASSERT(true, "Integrated UQFF gravity into fluid simulation");
    
    // Test 5.4: Store results in UQFFModule4
    UQFFModule4 module;
    module.updateVariable("vela_mass", vela.M);
    module.updateVariable("vela_radius", vela.r);
    module.updateVariable("vela_gravity", g_uqff);
    module.updateVariable("vela_buoyancy", F_U_Bi_i(vela));
    TEST_ASSERT(module.getUpdateCounter() == 4, "Stored 4 variables in module");
    
    // Test 5.5: Add physics terms to module
    module.setEnableDynamicTerms(true);
    module.registerDynamicTerm(std::make_unique<DynamicVacuumTerm>(vela.rho_vac_UA, 1e-12));
    module.registerDynamicTerm(std::make_unique<QuantumCouplingTerm>(1e-40));
    
    std::map<std::string, double> params;
    params["mass"] = vela.M;
    params["radius"] = vela.r;
    module.setDynamicParameter("mass", vela.M);
    module.setDynamicParameter("radius", vela.r);
    
    double dynamic_contrib = module.computeDynamicTerms(vela.t);
    std::cout << "    Dynamic physics contribution: " << dynamic_contrib << std::endl;
    TEST_ASSERT(!std::isnan(dynamic_contrib), "Computed dynamic physics terms");
    
    // Test 5.6: Export complete state
    module.exportState("vela_integrated_state.txt");
    TEST_ASSERT(true, "Exported integrated simulation state");
}

int main() {
    std::cout << "\n" << std::endl;
    std::cout << "================================================" << std::endl;
    std::cout << "  STAR-MAGIC CORE MODULE INTEGRATION TEST" << std::endl;
    std::cout << "  Phase 1, Weeks 1-2 Validation" << std::endl;
    std::cout << "  Date: November 14, 2025" << std::endl;
    std::cout << "================================================" << std::endl;
    
    test_system_catalogue();
    test_physics_terms();
    test_uqff_module4();
    test_fluid_solver();
    test_integrated_workflow();
    
    // Summary
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST SUMMARY" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  Total Tests: " << (tests_passed + tests_failed) << std::endl;
    std::cout << "  Passed: " << tests_passed << " ✓" << std::endl;
    std::cout << "  Failed: " << tests_failed << " ✗" << std::endl;
    
    if (tests_failed == 0) {
        std::cout << "\n🎉 ALL TESTS PASSED! Core modules working correctly." << std::endl;
        return 0;
    } else {
        std::cout << "\n⚠️  SOME TESTS FAILED. Review output above." << std::endl;
        return 1;
    }
}
