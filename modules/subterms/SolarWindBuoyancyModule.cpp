// SolarWindBuoyancyModule.cpp
#include "SolarWindBuoyancyModule.h"

// Constructor: Set framework defaults
SolarWindBuoyancyModule::SolarWindBuoyancyModule() {
    // Universal constants
    variables["epsilon_sw"] = 0.001;                // Buoyancy modulation (unitless)
    variables["rho_vac_sw"] = 8e-21;                // J/m^3 (solar wind energy density)
    variables["beta_1"] = 0.6;                      // From buoyancy coupling
    variables["U_g1"] = 1.39e26;                    // J/m^3 (Ug1 example)
    variables["Omega_g"] = 7.3e-16;                 // rad/s
    variables["M_bh"] = 8.15e36;                    // kg
    variables["d_g"] = 2.55e20;                     // m
    variables["E_react"] = 1.0;                     // Normalized
    variables["U_UA"] = 1.0;                        // Universal Aether factor
    variables["t_n"] = 0.0;                         // s
    variables["pi"] = 3.141592653589793;

    // Derived defaults
    variables["modulation_factor"] = computeModulationFactor();
}

// Update variable
void SolarWindBuoyancyModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "epsilon_sw" || name == "rho_vac_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SolarWindBuoyancyModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "epsilon_sw" || name == "rho_vac_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SolarWindBuoyancyModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ε_sw (fixed 0.001)
double SolarWindBuoyancyModule::computeEpsilon_sw() {
    return variables["epsilon_sw"];
}

// Compute modulation factor 1 + ε_sw * ρ_vac,sw
double SolarWindBuoyancyModule::computeModulationFactor() {
    return 1.0 + variables["epsilon_sw"] * variables["rho_vac_sw"];
}

// Compute example U_b1 with modulation
double SolarWindBuoyancyModule::computeU_b1() {
    double beta_1 = variables["beta_1"];
    double U_g1 = variables["U_g1"];
    double Omega_g = variables["Omega_g"];
    double M_bh_over_d_g = variables["M_bh"] / variables["d_g"];
    double E_react = variables["E_react"];
    double mod_factor = computeModulationFactor();
    double U_UA = variables["U_UA"];
    double cos_term = std::cos(variables["pi"] * variables["t_n"]);
    return -beta_1 * U_g1 * Omega_g * M_bh_over_d_g * E_react * mod_factor * U_UA * cos_term;
}

// Equation text
std::string SolarWindBuoyancyModule::getEquationText() {
    return "Modulation Factor = 1 + ε_sw * ρ_vac,sw\n"
           "Where ε_sw = 0.001 (unitless); ρ_vac,sw = 8e-21 J/m³.\n"
           "In U_bi: ... * (1 + ε_sw * ρ_vac,sw) * ... ≈1 (negligible correction ~8e-24).\n"
           "U_b1 = -β_1 U_g1 Ω_g (M_bh / d_g) * modulation * U_UA * cos(π t_n)\n"
           "≈ -1.94e27 J/m³ (at t_n=0, Sun params; modulation ≈1).\n"
           "Role: Minor solar wind density effect on buoyancy; stabilizes heliosphere/nebulae.\n"
           "UQFF: Scales counterforce to Ug; negligible but flexible for variations.";
}

// Print variables
void SolarWindBuoyancyModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "SolarWindBuoyancyModule.h"
// int main() {
//     SolarWindBuoyancyModule mod;
//     double mod_factor = mod.computeModulationFactor();
//     std::cout << "Modulation Factor = " << mod_factor << std::endl;
//     double u_b1 = mod.computeU_b1();
//     std::cout << "U_b1 = " << u_b1 << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("epsilon_sw", 0.002);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o sw_mod_test sw_mod_test.cpp SolarWindBuoyancyModule.cpp -lm
// Sample Output: Modulation ≈1.000000000000008; U_b1 ≈ -1.94e27 J/m³ (unchanged).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Coupling Constant  of Ugi.docx
File

Encode this attachment using the template.

Thoughts
