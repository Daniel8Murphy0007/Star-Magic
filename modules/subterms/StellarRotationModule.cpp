// StellarRotationModule.cpp
#include "StellarRotationModule.h"

// Constructor: Set framework defaults (Sun at t=0)
StellarRotationModule::StellarRotationModule() {
    // Universal constants
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["k_3"] = 1.8;                         // Coupling U_g3
    variables["B_j"] = 1e3;                         // T
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["lambda_i"] = 1.0;                    // Unitless U_i
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["f_TRZ"] = 0.1;                       // Unitless
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
    variables["t_n"] = 0.0;                         // s

    // Derived
    variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
    variables["day_to_s"] = 86400.0;                // s/day
}

// Update variable
void StellarRotationModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "rho_vac_SCm" || name == "rho_vac_UA") {
            variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void StellarRotationModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "rho_vac_SCm" || name == "rho_vac_UA") {
            variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void StellarRotationModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ω_s (rad/s)
double StellarRotationModule::computeOmega_s() {
    return variables["omega_s"];
}

// ω_s(t), simplified as constant (no t dep in example)
double StellarRotationModule::computeOmega_s_t(double t) {
    variables["t"] = t;
    return computeOmega_s();  // Constant for Sun
}

// Rotation period in days
double StellarRotationModule::computePeriod_days() {
    double period_s = 2.0 * M_PI / computeOmega_s();
    return period_s / variables["day_to_s"];
}

// U_g3 example
double StellarRotationModule::computeU_g3(double t) {
    double k_3 = variables["k_3"];
    double b_j = variables["B_j"];
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}

// U_i example
double StellarRotationModule::computeU_i(double t, double t_n) {
    double lambda_i = variables["lambda_i"];
    double rho_sc = variables["rho_vac_SCm"];
    double rho_ua = variables["rho_vac_UA"];
    double omega_s_t = computeOmega_s_t(t);
    double cos_pi_tn = std::cos(variables["pi"] * t_n);
    double trz_factor = 1.0 + variables["f_TRZ"];
    return lambda_i * rho_sc * rho_ua * omega_s_t * cos_pi_tn * trz_factor;
}

// Equation text
std::string StellarRotationModule::getEquationText() {
    return "U_g3 = k_3 * ∑ B_j * cos(ω_s(t) t π) * P_core * E_react\n"
           "U_i = λ_i * ρ_vac,[SCm] * ρ_vac,[UA] * ω_s(t) * cos(π t_n) * (1 + f_TRZ)\n"
           "Where ω_s = 2.5e-6 rad/s (~29-day Sun equatorial rotation);\n"
           "Scales rotational oscillations/inertia.\n"
           "Example t=0, t_n=0: U_g3 ≈1.8e49 J/m³; U_i ≈1.38e-47 J/m³.\n"
           "Role: Introduces spin in disk gravity/inertia; stellar/planetary dynamics.\n"
           "UQFF: Rotational effects in nebulae/disks/formation/mergers.";
}

// Print variables
void StellarRotationModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "StellarRotationModule.h"
// int main() {
//     StellarRotationModule mod;
//     double omega = mod.computeOmega_s();
//     std::cout << "ω_s = " << omega << " rad/s (~" << mod.computePeriod_days() << " days)\n";
//     double u_g3 = mod.computeU_g3(0.0);
//     std::cout << "U_g3 = " << u_g3 << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("omega_s", 3e-6);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o rotation_test rotation_test.cpp StellarRotationModule.cpp -lm
// Sample: ω_s=2.5e-6 rad/s (~29 days); U_g3≈1.8e49 J/m³; scales rotation.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Step Function.docx
File

Encode this attachment using the template.

Thoughts

You've reached the limit of 25 attachments in this conversation. Grok might forget content from earlier attachments; please start a new conversation to upload more attachments.
