// Test suite for Dimensional Analysis System
// Tests SI units, UQFF derived quantities, dimension checking
// Part of calculator_advanced framework extracted from Grok thread Iterations #30-32

#include "../include/dimensional_analysis.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

// Test helper macro
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        std::cerr << "TEST FAILED: " << message << std::endl; \
        return false; \
    }

// Test Case 1: Base Quantity Construction
bool test_base_quantities() {
    std::cout << "Running test_base_quantities..." << std::endl;
    
    DimensionalSystem system;
    
    // Test SI base quantities
    auto length = system.getBaseQuantity("length");
    TEST_ASSERT(length.dimensions[1] == 1, "Length should have L dimension");
    
    auto mass = system.getBaseQuantity("mass");
    TEST_ASSERT(mass.dimensions[0] == 1, "Mass should have M dimension");
    
    auto time = system.getBaseQuantity("time");
    TEST_ASSERT(time.dimensions[2] == 1, "Time should have T dimension");
    
    std::cout << "test_base_quantities PASSED" << std::endl;
    return true;
}

// Test Case 2: Derived Quantities (Force)
bool test_force_dimensions() {
    std::cout << "Running test_force_dimensions..." << std::endl;
    
    DimensionalSystem system;
    
    // Force = M * L * T^-2
    auto force = system.getBaseQuantity("force");
    TEST_ASSERT(force.dimensions[0] == 1, "Force should have M^1");
    TEST_ASSERT(force.dimensions[1] == 1, "Force should have L^1");
    TEST_ASSERT(force.dimensions[2] == -2, "Force should have T^-2");
    
    std::cout << "test_force_dimensions PASSED" << std::endl;
    return true;
}

// Test Case 3: Energy Dimensions
bool test_energy_dimensions() {
    std::cout << "Running test_energy_dimensions..." << std::endl;
    
    DimensionalSystem system;
    
    // Energy = M * L^2 * T^-2
    auto energy = system.getBaseQuantity("energy");
    TEST_ASSERT(energy.dimensions[0] == 1, "Energy should have M^1");
    TEST_ASSERT(energy.dimensions[1] == 2, "Energy should have L^2");
    TEST_ASSERT(energy.dimensions[2] == -2, "Energy should have T^-2");
    
    std::cout << "test_energy_dimensions PASSED" << std::endl;
    return true;
}

// Test Case 4: Addition Compatibility Check
bool test_addition_compatibility() {
    std::cout << "Running test_addition_compatibility..." << std::endl;
    
    // Same dimensions - should be compatible
    PhysicalQuantity length1{10.0, {0, 1, 0, 0, 0, 0, 0}};  // 10 meters
    PhysicalQuantity length2{5.0, {0, 1, 0, 0, 0, 0, 0}};   // 5 meters
    
    TEST_ASSERT(length1.isCompatibleWith(length2), "Same dimensions should be compatible");
    
    auto sum = length1 + length2;
    TEST_ASSERT(std::abs(sum.value - 15.0) < 1e-10, "10m + 5m should equal 15m");
    
    // Different dimensions - should throw exception
    PhysicalQuantity time1{2.0, {0, 0, 1, 0, 0, 0, 0}};  // 2 seconds
    bool threw_exception = false;
    try {
        auto invalid = length1 + time1;  // Should throw
    } catch (const std::exception& e) {
        threw_exception = true;
        std::cout << "  Expected exception: " << e.what() << std::endl;
    }
    TEST_ASSERT(threw_exception, "Adding incompatible dimensions should throw exception");
    
    std::cout << "test_addition_compatibility PASSED" << std::endl;
    return true;
}

// Test Case 5: Multiplication
bool test_multiplication() {
    std::cout << "Running test_multiplication..." << std::endl;
    
    // Force = mass * acceleration
    // mass: M^1, acceleration: L^1 * T^-2
    PhysicalQuantity mass{2.0, {1, 0, 0, 0, 0, 0, 0}};          // 2 kg
    PhysicalQuantity accel{10.0, {0, 1, -2, 0, 0, 0, 0}};       // 10 m/s^2
    
    auto force = mass * accel;
    std::cout << "  Force value: " << force.value << " N" << std::endl;
    TEST_ASSERT(std::abs(force.value - 20.0) < 1e-10, "2 kg * 10 m/s^2 = 20 N");
    TEST_ASSERT(force.dimensions[0] == 1, "Force M dimension");
    TEST_ASSERT(force.dimensions[1] == 1, "Force L dimension");
    TEST_ASSERT(force.dimensions[2] == -2, "Force T dimension");
    
    std::cout << "test_multiplication PASSED" << std::endl;
    return true;
}

// Test Case 6: Division
bool test_division() {
    std::cout << "Running test_division..." << std::endl;
    
    // Velocity = distance / time
    PhysicalQuantity distance{100.0, {0, 1, 0, 0, 0, 0, 0}};  // 100 meters
    PhysicalQuantity time{10.0, {0, 0, 1, 0, 0, 0, 0}};        // 10 seconds
    
    auto velocity = distance / time;
    std::cout << "  Velocity value: " << velocity.value << " m/s" << std::endl;
    TEST_ASSERT(std::abs(velocity.value - 10.0) < 1e-10, "100m / 10s = 10 m/s");
    TEST_ASSERT(velocity.dimensions[1] == 1, "Velocity L dimension");
    TEST_ASSERT(velocity.dimensions[2] == -1, "Velocity T^-1 dimension");
    
    std::cout << "test_division PASSED" << std::endl;
    return true;
}

// Test Case 7: Power Calculation (Energy/Time)
bool test_power_dimensions() {
    std::cout << "Running test_power_dimensions..." << std::endl;
    
    DimensionalSystem system;
    
    // Power = Energy / Time = M * L^2 * T^-3
    auto power = system.getBaseQuantity("power");
    TEST_ASSERT(power.dimensions[0] == 1, "Power should have M^1");
    TEST_ASSERT(power.dimensions[1] == 2, "Power should have L^2");
    TEST_ASSERT(power.dimensions[2] == -3, "Power should have T^-3");
    
    // Verify by calculation
    PhysicalQuantity energy{1000.0, {1, 2, -2, 0, 0, 0, 0}};  // 1000 J
    PhysicalQuantity time{10.0, {0, 0, 1, 0, 0, 0, 0}};        // 10 s
    auto calcPower = energy / time;
    
    TEST_ASSERT(std::abs(calcPower.value - 100.0) < 1e-10, "1000J / 10s = 100 W");
    TEST_ASSERT(calcPower.dimensions[0] == 1, "Power M dimension");
    TEST_ASSERT(calcPower.dimensions[1] == 2, "Power L dimension");
    TEST_ASSERT(calcPower.dimensions[2] == -3, "Power T dimension");
    
    std::cout << "test_power_dimensions PASSED" << std::endl;
    return true;
}

// Test Case 8: Frequency Dimensions
bool test_frequency_dimensions() {
    std::cout << "Running test_frequency_dimensions..." << std::endl;
    
    DimensionalSystem system;
    
    // Frequency = T^-1
    auto frequency = system.getBaseQuantity("frequency");
    TEST_ASSERT(frequency.dimensions[2] == -1, "Frequency should have T^-1");
    
    std::cout << "test_frequency_dimensions PASSED" << std::endl;
    return true;
}

// Test Case 9: UQFF Force Dimensions (F_U_Bi_i)
bool test_uqff_force() {
    std::cout << "Running test_uqff_force..." << std::endl;
    
    // F_U_Bi_i = F_rel * (E_cm/E_LEP) * Q_wave * g
    // Where g = G*M/r^2 has dimensions M^1 * L^1 * T^-2
    // Result should be force dimensions
    
    PhysicalQuantity F_rel{4.30e33, {1, 1, -2, 0, 0, 0, 0}};  // Force
    PhysicalQuantity dimensionless{1.0, {0, 0, 0, 0, 0, 0, 0}};  // (E_cm/E_LEP) and Q_wave
    PhysicalQuantity g{1e-10, {0, 1, -2, 0, 0, 0, 0}};  // Acceleration
    
    auto F_U_Bi_i = F_rel * dimensionless * dimensionless * g;
    
    // Check dimensions (should still be force-like after multiplying by acceleration gives extra L/T^2)
    std::cout << "  F_U_Bi_i dimensions: [" 
              << F_U_Bi_i.dimensions[0] << "," 
              << F_U_Bi_i.dimensions[1] << "," 
              << F_U_Bi_i.dimensions[2] << "]" << std::endl;
    
    TEST_ASSERT(F_U_Bi_i.dimensions[2] == -4, "Time dimension should be T^-4");
    
    std::cout << "test_uqff_force PASSED" << std::endl;
    return true;
}

// Test Case 10: Dimensional Validation
bool test_dimensional_validation() {
    std::cout << "Running test_dimensional_validation..." << std::endl;
    
    DimensionalSystem system;
    
    // Test that system validates dimensions correctly
    PhysicalQuantity valid{10.0, {0, 1, 0, 0, 0, 0, 0}};
    bool isValid = system.validateDimensions(valid);
    TEST_ASSERT(isValid, "Valid dimensions should pass validation");
    
    std::cout << "test_dimensional_validation PASSED" << std::endl;
    return true;
}

// Test Case 11: Gravitational Acceleration Dimensions
bool test_gravitational_dimensions() {
    std::cout << "Running test_gravitational_dimensions..." << std::endl;
    
    // g = G * M / r^2
    // G: M^-1 * L^3 * T^-2
    // M: M^1
    // r^2: L^2
    // Result: L * T^-2 (acceleration)
    
    PhysicalQuantity G{6.67e-11, {-1, 3, -2, 0, 0, 0, 0}};
    PhysicalQuantity M{1.989e30, {1, 0, 0, 0, 0, 0, 0}};
    PhysicalQuantity r{1e11, {0, 1, 0, 0, 0, 0, 0}};
    PhysicalQuantity r2 = r * r;
    
    auto g = (G * M) / r2;
    
    std::cout << "  g dimensions: [" 
              << g.dimensions[0] << "," 
              << g.dimensions[1] << "," 
              << g.dimensions[2] << "]" << std::endl;
    
    TEST_ASSERT(g.dimensions[0] == 0, "g should have no mass dimension");
    TEST_ASSERT(g.dimensions[1] == 1, "g should have L^1");
    TEST_ASSERT(g.dimensions[2] == -2, "g should have T^-2");
    
    std::cout << "test_gravitational_dimensions PASSED" << std::endl;
    return true;
}

// Main test runner
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Dimensional Analysis Test Suite" << std::endl;
    std::cout << "calculator_advanced framework" << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::vector<bool (*)(void)> tests = {
        test_base_quantities,
        test_force_dimensions,
        test_energy_dimensions,
        test_addition_compatibility,
        test_multiplication,
        test_division,
        test_power_dimensions,
        test_frequency_dimensions,
        test_uqff_force,
        test_dimensional_validation,
        test_gravitational_dimensions
    };
    
    int passed = 0;
    int failed = 0;
    
    for (auto test : tests) {
        try {
            if (test()) {
                passed++;
            } else {
                failed++;
            }
        } catch (const std::exception& e) {
            std::cerr << "EXCEPTION: " << e.what() << std::endl;
            failed++;
        }
        std::cout << std::endl;
    }
    
    std::cout << "========================================" << std::endl;
    std::cout << "Results: " << passed << " passed, " << failed << " failed" << std::endl;
    std::cout << "========================================" << std::endl;
    
    return (failed == 0) ? 0 : 1;
}
