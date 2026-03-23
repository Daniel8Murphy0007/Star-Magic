// SurfaceTemperatureModule.cpp
#include "SurfaceTemperatureModule.h"

// Constructor: Set framework defaults (Sun)
SurfaceTemperatureModule::SurfaceTemperatureModule() {
    // Universal constants
    variables["T_s"] = 5778.0;                      // K (Sun effective)
    variables["T_s_ref"] = 5778.0;                  // K (reference)
    variables["k_3"] = 1.8;                         // Coupling
    variables["B_ref"] = 1e3;                       // Base T (string)
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
}

// Update variable
void SurfaceTemperatureModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SurfaceTemperatureModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SurfaceTemperatureModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute T_s (K)
double SurfaceTemperatureModule::computeT_s() {
    return variables["T_s"];
}

// Hypothetical B_j scaled by T_s / T_s_ref
double SurfaceTemperatureModule::computeB_j_hypothetical(double t, double T_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Cycle
    return base_b * (T_s / variables["T_s_ref"]);
}

// U_g3 example with scaled B_j
double SurfaceTemperatureModule::computeU_g3_example(double t, double T_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j_hypothetical(t, T_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}

// Equation text
std::string SurfaceTemperatureModule::getEquationText() {
    return "B_j ≈ (10^3 + 0.4 sin(ω_s t)) * (T_s / T_s,ref) T (hypothetical);\n"
           "U_g3 = k_3 * ∑ B_j * cos(ω_s t π) * P_core * E_react\n"
           "Where T_s = 5778 K (Sun effective photosphere; °C=5505).\n"
           "T_s,ref=5778 K; scales string fields by temperature.\n"
           "Example t=0, T_s=5778 K: B_j≈1e3 T, U_g3≈1.8e49 J/m³;\n"
           "T_s=10000 K: B_j≈1730 T, U_g3≈3.11e49 J/m³ (+73%).\n"
           "Role: Thermal baseline for magnetic strength; variability in U_g3/disks.\n"
           "UQFF: Temperature-dependent fields; extensible for radiation/formation.";
}

// Print variables
void SurfaceTemperatureModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "SurfaceTemperatureModule.h"
// int main() {
//     SurfaceTemperatureModule mod;
//     double t_s = mod.computeT_s();
//     std::cout << "T_s = " << t_s << " K\n";
//     double u_g3 = mod.computeU_g3_example(0.0, 10000.0);
//     std::cout << "U_g3 (T_s=10000 K) = " << u_g3 << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("T_s", 6000.0);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o temp_test temp_test.cpp SurfaceTemperatureModule.cpp -lm
// Sample: T_s=5778 K; U_g3 (hot star)≈3.11e49 J/m³; thermal scaling.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Time Reversal Zone Factor.docx
File

Encode this attachment using the template.

Thoughts

You've reached the limit of 25 attachments in this conversation. Grok might forget content from earlier attachments; please start a new conversation to upload more attachments.
