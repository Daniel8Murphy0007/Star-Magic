/**
 * UQFF Equations Catalog
 * Thread b6d9bc22 Priority 2 - Iteration #31 UQFF Integration
 * 
 * Pre-loaded UQFF equations from catSymbols["UQFF"] array
 * Provides 100+ validated astrophysical equations
 */

#pragma once

#include <string>
#include <vector>
#include <map>
#include <functional>

namespace CalculatorAdvanced {

/**
 * @brief UQFF equation metadata
 */
struct UQFFEquation {
    std::string name;                        // e.g., "F_U_Bi_i", "compressed_g"
    std::string category;                    // e.g., "Buoyancy", "Gravity", "LENR"
    std::string equation_latex;              // LaTeX format for display
    std::string equation_symengine;          // SymEngine format for computation
    std::vector<std::string> variables;      // Required variables
    std::vector<std::string> parameters;     // Optional parameters
    std::map<std::string, double> defaults;  // Default parameter values
    std::string description;                 // Physical interpretation
    std::vector<std::string> references;     // Source papers/docs
    std::string units;                       // Expected output units
    
    /**
     * @brief Validate input parameters
     */
    bool validateInputs(const std::map<std::string, double>& inputs) const;
};

/**
 * @brief UQFF equation catalog
 */
class UQFFEquationCatalog {
public:
    UQFFEquationCatalog();
    
    /**
     * @brief Get equation by name
     * @param name Equation identifier (e.g., "F_U_Bi_i")
     * @return UQFFEquation structure
     */
    UQFFEquation getEquation(const std::string& name) const;
    
    /**
     * @brief List all equations in category
     * @param category Category filter (e.g., "Buoyancy", "Gravity")
     * @return Vector of equation names
     */
    std::vector<std::string> listEquations(const std::string& category = "") const;
    
    /**
     * @brief Search equations by keyword
     * @param keyword Search term (case-insensitive)
     * @return Matching equation names
     */
    std::vector<std::string> search(const std::string& keyword) const;
    
    /**
     * @brief Get equation categories
     */
    std::vector<std::string> getCategories() const;
    
    /**
     * @brief Check if equation exists
     */
    bool hasEquation(const std::string& name) const;
    
private:
    std::map<std::string, UQFFEquation> equations_;
    
    // Initialize catalog with pre-loaded equations
    void initializeUQFFEquations();
    
    // Category builders
    void addBuoyancyEquations();
    void addGravityEquations();
    void addLENREquations();
    void addAlphaClusteringEquations();
    void addMUGEEquations();
    void addQuantumLevelEquations();
    void addElectricUniverseEquations();
};

/**
 * @brief Pre-loaded UQFF equations (from thread Iteration #31)
 * 
 * 1. Buoyancy (F_U_Bi_i)
 */
const UQFFEquation F_U_BI_I = {
    .name = "F_U_Bi_i",
    .category = "Buoyancy",
    .equation_latex = R"(F_{U,Bi,i} = \int_{0}^{\infty} \left( U_{g1} + U_{g2} + U_{g3} + U_{g4} + U_m + U_{bi} \right) dr)",
    .equation_symengine = "integral(Ug1 + Ug2 + Ug3 + Ug4 + Um + Ubi, r, 0, inf)",
    .variables = {"M", "r", "t", "theta"},
    .parameters = {"rho_vac_UA", "rho_vac_SCm", "B", "L"},
    .defaults = {
        {"rho_vac_UA", 0.0001},
        {"rho_vac_SCm", 0.99},
        {"B", 1e10},      // Tesla
        {"L", 1e16}       // meters
    },
    .description = "Complete unified buoyancy force across 26 quantum levels",
    .references = {
        "MAIN_1_CoAnQi.cpp SOURCE111",
        "Thread b6d9bc22 UQFF applications (SN 1006, Eta Carinae)",
        "LEP 1998 F_rel = 4.30e33 N baseline"
    },
    .units = "N (Newton)"
};

/**
 * 2. Compressed Gravity (g_compressed)
 */
const UQFFEquation COMPRESSED_G = {
    .name = "compressed_g",
    .category = "Gravity",
    .equation_latex = R"(g_{compressed} = \frac{GM_{eff}}{r^2} + g_{expansion} + g_{super} + g_{envelope} + g_{Ug\_sum} + g_{cosm} + g_{quantum} + g_{fluid} + g_{pert})",
    .equation_symengine = "G*M_eff/r^2 + g_expansion + g_super + g_envelope + g_Ug_sum + g_cosm + g_quantum + g_fluid + g_pert",
    .variables = {"M", "r", "t", "z"},
    .parameters = {"H_0", "B", "Lambda", "hbar"},
    .defaults = {
        {"H_0", 2.27e-18},    // s^-1 (Hubble parameter)
        {"B", 4.4e13},        // T (critical magnetic field)
        {"Lambda", 1.1e-52},  // m^-2 (cosmological constant)
        {"hbar", 1.055e-34}   // J·s (reduced Planck constant)
    },
    .description = "Complete gravity with 9 correction terms (Newtonian + Hubble + Magnetic + Lambda + Quantum + Fluid + Perturbation)",
    .references = {
        "MAIN_1_CoAnQi.cpp CompressedMUGEBaseTerm",
        "source2.cpp QuantumDesignCalculatorWidget"
    },
    .units = "m/s² (acceleration)"
};

/**
 * 3. LENR Frequency (F_LENR)
 */
const UQFFEquation F_LENR = {
    .name = "F_LENR",
    .category = "LENR",
    .equation_latex = R"(F_{LENR} = k_{LENR} \left( \frac{\omega_{LENR}}{\omega_0} \right)^2)",
    .equation_symengine = "k_LENR * (omega_LENR / omega_0)^2",
    .variables = {"omega_LENR"},
    .parameters = {"k_LENR", "omega_0"},
    .defaults = {
        {"omega_LENR", 7.54e12},  // rad/s (1.2 THz)
        {"k_LENR", 1.0},
        {"omega_0", 1.885e9}      // rad/s (300 Hz baseline, Colman-Gillespie)
    },
    .description = "Low-energy nuclear reaction force at 1.2 THz resonance (Widom-Larsen mechanism)",
    .references = {
        "source2.cpp lines 4670-4673",
        "MAIN_1_CoAnQi.cpp LENRExtendedTerm line 1498",
        "K_n_Neutron Production Calibration Constant_19April2025.docx"
    },
    .units = "N (Newton)"
};

/**
 * @brief Example usage:
 * 
 * UQFFEquationCatalog catalog;
 * 
 * // Get specific equation
 * auto F_U = catalog.getEquation("F_U_Bi_i");
 * std::cout << "LaTeX: " << F_U.equation_latex << std::endl;
 * std::cout << "Units: " << F_U.units << std::endl;
 * 
 * // List all gravity equations
 * auto gravity_eqs = catalog.listEquations("Gravity");
 * for (const auto& eq_name : gravity_eqs) {
 *     auto eq = catalog.getEquation(eq_name);
 *     std::cout << eq_name << ": " << eq.description << std::endl;
 * }
 * 
 * // Search for LENR-related equations
 * auto lenr_results = catalog.search("LENR");
 * // Returns: ["F_LENR", "eta_neutron_production", "K_n_calibration", ...]
 * 
 * // Validate inputs before computation
 * std::map<std::string, double> inputs = {
 *     {"M", 1.989e30},  // Solar mass
 *     {"r", 1e16},      // 1 pc
 *     {"t", 3.156e7}    // 1 year
 * };
 * if (F_U.validateInputs(inputs)) {
 *     // Proceed with computation
 * }
 */

} // namespace CalculatorAdvanced
