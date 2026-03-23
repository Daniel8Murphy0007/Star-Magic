// PiConstantModule.cpp
#include "PiConstantModule.h"

// Constructor: Set framework defaults
PiConstantModule::PiConstantModule() {
    // Mathematical constants
    variables["pi"] = 3.141592653589793;            // Unitless
    variables["t_n"] = 0.0;                         // Negative time factor
    variables["t"] = 0.0;                           // Time
    variables["period"] = 3.96e8;                   // s (example solar cycle)
    variables["omega_c"] = 2.0 * variables["pi"] / variables["period"];  // rad/s
    variables["base_mu"] = 3.38e20;                 // T·m^3
    variables["B_j"] = 1e3;                         // Base T
}

// Update variable
void PiConstantModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "period") {
            variables["omega_c"] = 2.0 * variables["pi"] / value;
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void PiConstantModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "period") {
            variables["omega_c"] = 2.0 * variables["pi"] / variables[name];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void PiConstantModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute π
double PiConstantModule::computePi() {
    return variables["pi"];
}

// Compute cos(π t_n)
double PiConstantModule::computeCosPiTn(double t_n) {
    variables["t_n"] = t_n;
    return std::cos(computePi() * t_n);
}

// Compute sin(ω_c t)
double PiConstantModule::computeSinOmegaCT(double t) {
    variables["t"] = t;
    return std::sin(variables["omega_c"] * t);
}

// Example μ_j = (10^3 + 0.4 sin(ω_c t)) * 3.38e20 T·m^3
double PiConstantModule::computeMuJExample(double t) {
    double sin_omega = computeSinOmegaCT(t);
    double b_j = variables["B_j"] + 0.4 * sin_omega;
    return b_j * variables["base_mu"];
}

// Example cos(π t_n) in U_g1
double PiConstantModule::computeUg1CosTerm(double t_n) {
    return computeCosPiTn(t_n);
}

// Equation text
std::string PiConstantModule::getEquationText() {
    return "π ≈ 3.141592653589793 (unitless mathematical constant).\n"
           "Role: Defines periodicity in oscillations; C=2π r; trig args (sin/cos with 2π cycle).\n"
           "In U_m: μ_j = (10^3 + 0.4 sin(ω_c t)) * 3.38e20; ω_c = 2π / period.\n"
           "In U_g1: ... cos(π t_n) ... (time-reversal oscillations).\n"
           "Example t=0, t_n=0: sin(ω_c t)=0 → μ_j=3.38e23 T·m^3; cos(π t_n)=1.\n"
           "UQFF: Ensures cyclic/TRZ dynamics; solar cycles, rotations in nebulae/quasars.";
}

// Print variables
void PiConstantModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "PiConstantModule.h"
// int main() {
//     PiConstantModule mod;
//     double pi_val = mod.computePi();
//     std::cout << "π ≈ " << pi_val << std::endl;
//     double mu = mod.computeMuJExample(0.0);
//     std::cout << "μ_j (t=0) = " << mu << " T·m^3\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("t_n", 1.0);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o pi_test pi_test.cpp PiConstantModule.cpp -lm
// Sample: π=3.14159; μ_j≈3.38e23 T·m^3; cos(π*0)=1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Planetary Penetration factor.docx
File

Encode this attachment using the template.

Thoughts
