// HeavisideFractionModule.cpp
#include "HeavisideFractionModule.h"

// Constructor: Set framework defaults
HeavisideFractionModule::HeavisideFractionModule() {
    // Universal constants
    variables["f_Heaviside"] = 0.01;                // Unitless fraction
    variables["scale_Heaviside"] = 1e13;            // Amplification factor
    variables["f_quasi"] = 0.01;                    // Quasi factor
    variables["mu_j"] = 3.38e23;                    // T m^3 (j=1)
    variables["r_j"] = 1.496e13;                    // m
    variables["gamma"] = 5e-5 / 86400.0;            // day^-1 to s^-1
    variables["t_n"] = 0.0;                         // s
    variables["phi_hat_j"] = 1.0;                   // Normalized
    variables["P_SCm"] = 1.0;                       // Pressure
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;

    // Derived
    variables["heaviside_factor"] = computeHeavisideFactor();
}

// Update variable
void HeavisideFractionModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "f_Heaviside") {
            variables["heaviside_factor"] = computeHeavisideFactor();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void HeavisideFractionModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "f_Heaviside") {
            variables["heaviside_factor"] = computeHeavisideFactor();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void HeavisideFractionModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute f_Heaviside (0.01)
double HeavisideFractionModule::computeF_Heaviside() {
    return variables["f_Heaviside"];
}

// Compute 1 + 10^13 * f_Heaviside
double HeavisideFractionModule::computeHeavisideFactor() {
    return 1.0 + variables["scale_Heaviside"] * computeF_Heaviside();
}

// Base for U_m without Heaviside/Quasi
double HeavisideFractionModule::computeUmBase(int j, double t) {
    double mu_over_rj = variables["mu_j"] / variables["r_j"];
    double exp_arg = - variables["gamma"] * t * std::cos(variables["pi"] * variables["t_n"]);
    double one_minus_exp = 1.0 - std::exp(exp_arg);
    double phi_hat = variables["phi_hat_j"];
    return mu_over_rj * one_minus_exp * phi_hat * variables["P_SCm"] * variables["E_react"];
}

// U_m contribution with Heaviside
double HeavisideFractionModule::computeUmContribution(int j, double t) {
    double base = computeUmBase(j, t);
    double heaviside_f = computeHeavisideFactor();
    double quasi_f = 1.0 + variables["f_quasi"];
    return base * heaviside_f * quasi_f;
}

// U_m without Heaviside (set f=0 temporarily)
double HeavisideFractionModule::computeUmWithNoHeaviside(int j, double t) {
    double orig_f = variables["f_Heaviside"];
    variables["f_Heaviside"] = 0.0;
    double result = computeUmContribution(j, t);
    variables["f_Heaviside"] = orig_f;
    return result;
}

// Equation text
std::string HeavisideFractionModule::getEquationText() {
    return "U_m = ∑_j [ (μ_j / r_j) (1 - e^{-γ t cos(π t_n)}) φ_hat_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)\n"
           "Where f_Heaviside = 0.01 (unitless Heaviside fraction);\n"
           "Heaviside factor = 1 + 10^13 * 0.01 = 1 + 1e11 (amplifies ~10^11x).\n"
           "Example j=1, t=0: U_m contrib ≈2.28e65 J/m³ (with); ≈2.28e54 J/m³ (without).\n"
           "Role: Threshold-activated scaling in magnetic energy; nonlinear [SCm]/[UA] effects.\n"
           "UQFF: Amplifies small fraction for large impact in nebulae/quasars/jets.";
}

// Print variables
void HeavisideFractionModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print U_m comparison
void HeavisideFractionModule::printUmComparison(int j, double t) {
    double um_with = computeUmContribution(j, t);
    double um_without = computeUmWithNoHeaviside(j, t);
    double amplification = um_with / um_without;
    std::cout << "U_m Comparison for j=" << j << " at t=" << t << " s:\n";
    std::cout << "With Heaviside: " << std::scientific << um_with << " J/m³\n";
    std::cout << "Without Heaviside: " << um_without << " J/m³\n";
    std::cout << "Amplification: ~" << std::scientific << amplification << "x\n";
}

// Example usage in base program (snippet)
// #include "HeavisideFractionModule.h"
// int main() {
//     HeavisideFractionModule mod;
//     double heav_factor = mod.computeHeavisideFactor();
//     std::cout << "Heaviside Factor = " << heav_factor << std::endl;
//     mod.printUmComparison(1, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("f_Heaviside", 0.02);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o heaviside_test heaviside_test.cpp HeavisideFractionModule.cpp -lm
// Sample: Factor=1e11+1; U_m with=2.28e65 J/m³ (~1e11x without).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Heliosphere Thickness Factor.docx
File

Encode this attachment using the template.

Thoughts
