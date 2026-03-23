// SolarCycleFrequencyModule.cpp
#include "SolarCycleFrequencyModule.h"

// Constructor: Set framework defaults
SolarCycleFrequencyModule::SolarCycleFrequencyModule() {
    // Universal constants
    variables["pi"] = 3.141592653589793;
    variables["period"] = 3.96e8;                   // s (~12.55 years)
    variables["base_mu"] = 3.38e20;                 // T·m³
    variables["B_j"] = 1e3;                         // Base T
    variables["t"] = 0.0;                           // s

    // Derived
    variables["omega_c"] = computeOmega_c();
}

// Update variable
void SolarCycleFrequencyModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "period") {
            variables["omega_c"] = computeOmega_c();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SolarCycleFrequencyModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "period") {
            variables["omega_c"] = computeOmega_c();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SolarCycleFrequencyModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ω_c = 2π / period
double SolarCycleFrequencyModule::computeOmega_c() {
    return 2.0 * variables["pi"] / variables["period"];
}

// Compute sin(ω_c t)
double SolarCycleFrequencyModule::computeSinOmegaCT(double t) {
    variables["t"] = t;
    return std::sin(variables["omega_c"] * t);
}

// Example μ_j = (10^3 + 0.4 sin(ω_c t)) * 3.38e20
double SolarCycleFrequencyModule::computeMuJExample(double t) {
    double sin_omega = computeSinOmegaCT(t);
    double b_j = variables["B_j"] + 0.4 * sin_omega;
    return b_j * variables["base_mu"];
}

// Equation text
std::string SolarCycleFrequencyModule::getEquationText() {
    return "ω_c = 2π / 3.96e8 s⁻¹ ≈1.59e-8 rad/s (period ~12.55 yr, near 11-yr solar cycle);\n"
           "In U_m: μ_j = (10^3 + 0.4 sin(ω_c t)) * 3.38e20 T·m³ (cyclic magnetic variation).\n"
           "In U_g3: ... cos(ω_s t π) ... (ω_s Sun rotation, but ω_c for cycle).\n"
           "Example t=0: sin=0 → μ_j=3.38e23 T·m³;\n"
           "t=3.14e7 s (~1 yr): sin≈0.477 → μ_j≈3.381e23 T·m³ (+0.019%).\n"
           "Role: Models solar cycle periodicity; magnetic activity in strings/fields.\n"
           "UQFF: Cyclic effects in jets/nebulae/formation; near 11-yr Hale cycle.";
}

// Print variables
void SolarCycleFrequencyModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "SolarCycleFrequencyModule.h"
// int main() {
//     SolarCycleFrequencyModule mod;
//     double omega = mod.computeOmega_c();
//     std::cout << "ω_c ≈ " << omega << " rad/s\n";
//     double mu = mod.computeMuJExample(0.0);
//     std::cout << "μ_j (t=0) = " << mu << " T·m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("period", 3.8e8);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o solar_cycle_test solar_cycle_test.cpp SolarCycleFrequencyModule.cpp -lm
// Sample: ω_c≈1.59e-8 rad/s; μ_j≈3.38e23 T·m³; periodic variation.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Solar Wind Modulation Factor.docx
File

Encode this attachment using the template.

Thoughts

You've reached the limit of 25 attachments in this conversation. Grok might forget content from earlier attachments; please start a new conversation to upload more attachments.
