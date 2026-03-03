/**
 * Dimensional Analysis System
 * Thread b6d9bc22 Priority 2 - Iteration #31 UQFF Integration
 * 
 * Automated unit checking for physics equations
 * Ensures dimensional consistency across UQFF framework
 */

#pragma once

#include <string>
#include <map>
#include <vector>
#include <array>

namespace CalculatorAdvanced {

/**
 * @brief SI base units: [M, L, T, I, Θ, N, J]
 * M = Mass (kg)
 * L = Length (m)
 * T = Time (s)
 * I = Electric current (A)
 * Θ = Temperature (K)
 * N = Amount of substance (mol)
 * J = Luminous intensity (cd)
 */
using DimensionVector = std::array<int, 7>;

/**
 * @brief Physical quantity with dimensions
 */
struct PhysicalQuantity {
    std::string symbol;             // e.g., "F_U_Bi_i", "M_BH"
    std::string description;        // Human-readable name
    DimensionVector dimensions;     // [M, L, T, I, Θ, N, J] exponents
    std::string si_unit;            // e.g., "N" (Newton), "kg/m^3"
    
    /**
     * @brief Common UQFF quantities
     */
    static PhysicalQuantity force();           // [M L T^-2] = N
    static PhysicalQuantity mass();            // [M] = kg
    static PhysicalQuantity length();          // [L] = m
    static PhysicalQuantity time();            // [T] = s
    static PhysicalQuantity energy();          // [M L^2 T^-2] = J
    static PhysicalQuantity energy_density();  // [M L^-1 T^-2] = J/m³
    static PhysicalQuantity acceleration();    // [L T^-2] = m/s²
    static PhysicalQuantity angular_freq();    // [T^-1] = rad/s
    static PhysicalQuantity magnetic_field();  // [M T^-2 I^-1] = T (Tesla)
    static PhysicalQuantity electric_field();  // [M L T^-3 I^-1] = V/m
};

/**
 * @brief Dimensional analysis engine
 */
class DimensionalSystem {
public:
    DimensionalSystem();
    
    /**
     * @brief Register UQFF quantity
     * @param quantity Physical quantity with dimensions
     */
    void registerQuantity(const PhysicalQuantity& quantity);
    
    /**
     * @brief Check if equation is dimensionally consistent
     * @param lhs Left-hand side quantity name
     * @param rhs Right-hand side quantity name
     * @return true if [lhs] == [rhs]
     */
    bool checkConsistency(const std::string& lhs, const std::string& rhs) const;
    
    /**
     * @brief Compute dimensions of expression
     * @param expression Symbolic expression (e.g., "M * v^2")
     * @param var_dimensions Map of variable → dimensions
     * @return Resulting dimensions
     */
    DimensionVector computeDimensions(
        const std::string& expression,
        const std::map<std::string, DimensionVector>& var_dimensions
    ) const;
    
    /**
     * @brief Get SI unit string for dimensions
     * @param dims Dimension vector
     * @return SI unit (e.g., "m/s^2", "kg*m/s^2" = "N")
     */
    std::string getSIUnit(const DimensionVector& dims) const;
    
    /**
     * @brief Check if dimensions are dimensionless [0,0,0,0,0,0,0]
     */
    bool isDimensionless(const DimensionVector& dims) const;
    
    /**
     * @brief List all registered UQFF quantities
     */
    std::vector<PhysicalQuantity> listQuantities() const;
    
private:
    std::map<std::string, PhysicalQuantity> quantities_;
    
    // Pre-register common UQFF quantities
    void registerUQFFQuantities();
};

/**
 * @brief Dimensional operators
 */
DimensionVector operator+(const DimensionVector& a, const DimensionVector& b);  // Multiply units
DimensionVector operator-(const DimensionVector& a, const DimensionVector& b);  // Divide units
DimensionVector operator*(const DimensionVector& dims, int exponent);           // Raise to power
bool operator==(const DimensionVector& a, const DimensionVector& b);
bool operator!=(const DimensionVector& a, const DimensionVector& b);

/**
 * @brief Example usage:
 * 
 * DimensionalSystem dims;
 * 
 * // Register UQFF quantities
 * dims.registerQuantity(PhysicalQuantity::force());  // F_U_Bi_i
 * dims.registerQuantity(PhysicalQuantity::mass());   // M_BH
 * 
 * // Check equation: F = M * a
 * auto F_dims = PhysicalQuantity::force().dimensions;       // [1, 1, -2, 0, 0, 0, 0]
 * auto M_dims = PhysicalQuantity::mass().dimensions;        // [1, 0,  0, 0, 0, 0, 0]
 * auto a_dims = PhysicalQuantity::acceleration().dimensions; // [0, 1, -2, 0, 0, 0, 0]
 * auto Ma_dims = M_dims + a_dims;  // [1, 1, -2, 0, 0, 0, 0]
 * 
 * assert(F_dims == Ma_dims);  // ✅ Dimensionally consistent
 * 
 * // Check UQFF equation: compressed_g = G*M/r^2
 * dims.checkConsistency("compressed_g", "G*M/r^2");  // ✅ true
 * 
 * // Verify energy density: [SCm] = J/m³
 * auto SCm_dims = PhysicalQuantity::energy_density().dimensions;
 * std::string unit = dims.getSIUnit(SCm_dims);  // "J/m³"
 */

} // namespace CalculatorAdvanced
