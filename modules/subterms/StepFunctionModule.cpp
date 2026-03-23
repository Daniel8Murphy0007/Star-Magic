// StepFunctionModule.cpp
#include "StepFunctionModule.h"

// Constructor: Set framework defaults (Sun at r=R_b)
StepFunctionModule::StepFunctionModule() {
    // Universal constants
    variables["R_b"] = 1.496e13;                    // m (100 AU)
    variables["k_2"] = 1.2;                         // Coupling
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["M_s"] = 1.989e30;                    // kg
    variables["r"] = 1.496e13;                      // m (default = R_b)
    variables["delta_sw"] = 0.01;                   // Unitless
    variables["v_sw"] = 5e5;                        // m/s
    variables["H_SCm"] = 1.0;                       // Unitless
    variables["E_react"] = 1e46;                    // J

    // Derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Update variable
void StepFunctionModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        } else if (name == "delta_sw" || name == "v_sw") {
            variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void StepFunctionModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        } else if (name == "delta_sw" || name == "v_sw") {
            variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void StepFunctionModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute S(r - R_b): 1 if r > R_b, 0 otherwise (treat = as 1 per examples)
double StepFunctionModule::computeS_r_Rb(double r) {
    return (r >= variables["R_b"]) ? 1.0 : 0.0;
}

// Compute U_g2 with S(r - R_b)
double StepFunctionModule::computeU_g2(double r) {
    variables["r"] = r;
    double k_2 = variables["k_2"];
    double rho_sum = variables["rho_sum"];
    double M_s = variables["M_s"];
    double s_step = computeS_r_Rb(r);
    double swirl_factor = variables["swirl_factor"];
    double h_scm = variables["H_SCm"];
    double e_react = variables["E_react"];
    return k_2 * (rho_sum * M_s / (r * r)) * s_step * swirl_factor * h_scm * e_react;
}

// Equation text
std::string StepFunctionModule::getEquationText() {
    return "U_g2 = k_2 * [(ρ_vac,[UA] + ρ_vac,[SCm]) M_s / r^2] * S(r - R_b) * (1 + δ_sw v_sw) * H_SCm * E_react\n"
           "Where S(r - R_b) = 1 (r > R_b), 0 otherwise (Heaviside step; =1 at boundary).\n"
           "Defines outer bubble activation beyond R_b=1.496e13 m (100 AU).\n"
           "Example r=1.496e13 m: S=1, U_g2 ≈1.18e53 J/m³;\n"
           "r=1e11 m: S=0, U_g2=0; r=1e14 m: S=1, U_g2≈1.18e51 J/m³.\n"
           "Role: Sharp transition internal/external gravity; heliopause-like boundary.\n"
           "UQFF: Separates regimes for heliodynamics/nebular formation.";
}

// Print variables
void StepFunctionModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "StepFunctionModule.h"
// int main() {
//     StepFunctionModule mod;
//     double s = mod.computeS_r_Rb(1.5e13);
//     std::cout << "S(1.5e13 - R_b) = " << s << std::endl;
//     double u_g2 = mod.computeU_g2(1e11);  // Inside
//     std::cout << "U_g2 (inside) = " << u_g2 << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("R_b", 2e13);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o step_test step_test.cpp StepFunctionModule.cpp -lm
// Sample: S> R_b=1; U_g2=0 inside; activates outer bubble.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Stress Energy Tensor.docx
File

Encode this attachment using the template.

Thoughts

You've reached the limit of 25 attachments in this conversation. Grok might forget content from earlier attachments; please start a new conversation to upload more attachments.
