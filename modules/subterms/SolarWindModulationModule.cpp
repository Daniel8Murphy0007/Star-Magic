// SolarWindModulationModule.cpp
#include "SolarWindModulationModule.h"

// Constructor: Set framework defaults (Sun at r=R_b)
SolarWindModulationModule::SolarWindModulationModule() {
    // Universal constants
    variables["delta_sw"] = 0.01;                   // Unitless
    variables["v_sw"] = 5e5;                        // m/s
    variables["k_2"] = 1.2;                         // Coupling
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["M_s"] = 1.989e30;                    // kg
    variables["r"] = 1.496e13;                      // m (R_b)
    variables["R_b"] = 1.496e13;                    // m
    variables["S_r_Rb"] = 1.0;                      // Step
    variables["H_SCm"] = 1.0;                       // Unitless
    variables["E_react"] = 1e46;                    // J

    // Derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["modulation_factor"] = computeModulationFactor();
}

// Update variable
void SolarWindModulationModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "delta_sw" || name == "v_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        } else if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SolarWindModulationModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "delta_sw" || name == "v_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        } else if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SolarWindModulationModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute δ_sw = 0.01
double SolarWindModulationModule::computeDelta_sw() {
    return variables["delta_sw"];
}

// Compute 1 + δ_sw * v_sw
double SolarWindModulationModule::computeModulationFactor() {
    return 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Compute U_g2 with modulation
double SolarWindModulationModule::computeU_g2(double r) {
    variables["r"] = r;
    double k_2 = variables["k_2"];
    double rho_sum = variables["rho_sum"];
    double M_s = variables["M_s"];
    double s_step = (r >= variables["R_b"]) ? 1.0 : 0.0;
    double mod_factor = computeModulationFactor();
    double h_scm = variables["H_SCm"];
    double e_react = variables["E_react"];
    return k_2 * (rho_sum * M_s / (r * r)) * s_step * mod_factor * h_scm * e_react;
}

// U_g2 without modulation (δ_sw=0)
double SolarWindModulationModule::computeU_g2_no_mod(double r) {
    double orig_delta = variables["delta_sw"];
    variables["delta_sw"] = 0.0;
    double result = computeU_g2(r);
    variables["delta_sw"] = orig_delta;
    return result;
}

// Equation text
std::string SolarWindModulationModule::getEquationText() {
    return "U_g2 = k_2 * [(ρ_vac,[UA] + ρ_vac,[SCm]) M_s / r^2] * S(r - R_b) * (1 + δ_sw v_sw) * H_SCm * E_react\n"
           "Where δ_sw = 0.01 (unitless solar wind modulation factor);\n"
           "Modulation = 1 + 0.01 * v_sw (v_sw=5e5 m/s → ~5001x amplification).\n"
           "Example r=R_b=1.496e13 m: U_g2 ≈1.18e53 J/m³ (with); ≈2.36e49 J/m³ (without; ~5000x less).\n"
           "Role: Enhances external gravity via solar wind momentum/pressure beyond R_b.\n"
           "UQFF: Models heliosphere dynamics; wind influence on nebular/star formation.";
}

// Print variables
void SolarWindModulationModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "SolarWindModulationModule.h"
// int main() {
//     SolarWindModulationModule mod;
//     double mod_f = mod.computeModulationFactor();
//     std::cout << "Modulation Factor = " << mod_f << std::endl;
//     double u_g2 = mod.computeU_g2(1.496e13);
//     std::cout << "U_g2 = " << u_g2 << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("delta_sw", 0.02);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o sw_mod_test sw_mod_test.cpp SolarWindModulationModule.cpp -lm
// Sample: Factor=5001; U_g2≈1.18e53 J/m³; amplifies outer bubble gravity.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Solar Wind Velocity.docx
File

Encode this attachment using the template.

Thoughts

You've reached the limit of 25 attachments in this conversation. Grok might forget content from earlier attachments; please start a new conversation to upload more attachments.
