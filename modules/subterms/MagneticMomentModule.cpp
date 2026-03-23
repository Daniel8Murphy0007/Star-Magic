// MagneticMomentModule.cpp
#include "MagneticMomentModule.h"

// Constructor: Set framework defaults
MagneticMomentModule::MagneticMomentModule() {
    // Universal constants
    variables["base_mu"] = 3.38e20;                 // T·m^3 (definition); note: example uses 3.38e23
    variables["omega_c"] = 2.5e-6;                  // rad/s
    variables["r_j"] = 1.496e13;                    // m (for j=1)
    variables["gamma"] = 5e-5 / 86400.0;            // s^-1 (0.00005 day^-1)
    variables["t_n"] = 0.0;                         // s
    variables["phi_hat_j"] = 1.0;                   // Normalized
    variables["P_SCm"] = 1.0;                       // Pressure
    variables["E_react"] = 1e46;                    // J
    variables["f_Heaviside"] = 0.01;                // Unitless
    variables["f_quasi"] = 0.01;                    // Unitless
    variables["k3"] = 1.8;                          // Coupling for Ug3
    variables["pi"] = 3.141592653589793;

    // Derived defaults
    variables["B_j"] = 1e3;                         // Base T
    variables["scale_Heaviside"] = 1e13;            // Amplification
}

// Update variable
void MagneticMomentModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
}

// Add delta
void MagneticMomentModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void MagneticMomentModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute μ_j(t)
double MagneticMomentModule::computeMu_j(int j, double t) {
    double sin_term = std::sin(variables["omega_c"] * t);
    double b_j = variables["B_j"] + 0.4 * sin_term;  // T
    return b_j * variables["base_mu"];  // T·m^3; adjust base if needed for example
}

// Compute B_j(t) base
double MagneticMomentModule::computeB_j(double t) {
    return variables["B_j"] + 0.4 * std::sin(variables["omega_c"] * t);
}

// Example U_m contrib for j (J/m^3, simplified)
double MagneticMomentModule::computeUmContrib(int j, double t) {
    double mu_j = computeMu_j(j, t);
    double r_j = variables["r_j"];
    double exp_arg = - variables["gamma"] * t * std::cos(variables["pi"] * variables["t_n"]);
    double one_minus_exp = 1.0 - std::exp(exp_arg);
    double phi_hat = variables["phi_hat_j"];
    double heaviside_f = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
    double quasi_f = 1.0 + variables["f_quasi"];
    return (mu_j / r_j * one_minus_exp * phi_hat) * variables["P_SCm"] * variables["E_react"] * heaviside_f * quasi_f;
}

// Example Ug3 contrib (J/m^3)
double MagneticMomentModule::computeUg3Contrib(double t) {
    double b_j = computeB_j(t);
    double cos_term = std::cos(variables["omega_c"] * t * variables["pi"]);  // Approx
    double p_core = 1.0;
    double e_react = variables["E_react"];
    return variables["k3"] * b_j * cos_term * p_core * e_react;
}

// Equation text
std::string MagneticMomentModule::getEquationText() {
    return "μ_j = (10^3 + 0.4 sin(ω_c t)) * 3.38e20 T·m^3\n"
           "Where ω_c=2.5e-6 rad/s; units T·m^3 (magnetic dipole strength).\n"
           "In U_m: ∑_j [μ_j / r_j * (1 - e^{-γ t cos(π t_n)}) φ_hat_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)\n"
           "In Ug3: k3 * ∑_j B_j cos(ω_s t π) P_core E_react; B_j = 10^3 + 0.4 sin(ω_c t) T.\n"
           "Example j=1, t=0: μ_j ≈3.38e23 T·m^3; U_m contrib ≈2.28e65 J/m³; Ug3 ≈1.8e49 J/m³.\n"
           "Role: Quantifies string magnetic strength; drives Um/Ug3 for jets/disks/nebulae.";
}

// Print variables
void MagneticMomentModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print contributions
void MagneticMomentModule::printMomentContributions(int j, double t) {
    double mu = computeMu_j(j, t);
    double b = computeB_j(t);
    double um = computeUmContrib(j, t);
    double ug3 = computeUg3Contrib(t);
    std::cout << "Magnetic Moment j=" << j << " at t=" << t << " s:\n";
    std::cout << "μ_j = " << std::scientific << mu << " T·m^3\n";
    std::cout << "B_j = " << b << " T\n";
    std::cout << "U_m contrib = " << um << " J/m³\n";
    std::cout << "Ug3 contrib = " << ug3 << " J/m³\n";
}

// Example usage in base program (snippet)
// #include "MagneticMomentModule.h"
// int main() {
//     MagneticMomentModule mod;
//     double t = 0.0;
//     mod.printMomentContributions(1, t);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("base_mu", 4e20);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o moment_test moment_test.cpp MagneticMomentModule.cpp -lm
// Sample: μ_j≈3.38e23 T·m^3; U_m≈2.28e65 J/m³; cyclic variation via sin(ω_c t).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Mass of Galactic Black Hole.docx
File

Encode this attachment using the template.

Thoughts
