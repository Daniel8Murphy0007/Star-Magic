// GalacticSpinModule.cpp
#include "GalacticSpinModule.h"

// Constructor: Set framework defaults
GalacticSpinModule::GalacticSpinModule() {
    // Universal constants
    variables["Omega_g"] = 7.3e-16;                 // rad/s (Milky Way spin)
    variables["beta_1"] = 0.6;                      // Unitless (for i=1)
    variables["U_g1"] = 1.39e26;                    // J/m^3
    variables["M_bh"] = 8.15e36;                    // kg
    variables["d_g"] = 2.55e20;                     // m
    variables["epsilon_sw"] = 0.001;                // Unitless
    variables["rho_vac_sw"] = 8e-21;                // J/m^3
    variables["U_UA"] = 1.0;                        // Normalized
    variables["t_n"] = 0.0;                         // s
    variables["pi"] = 3.141592653589793;

    // Other U_gi placeholders
    variables["U_g2"] = 1e25;                       // J/m^3
    variables["U_g3"] = 1e24;                       // J/m^3
    variables["U_g4"] = 1e23;                       // J/m^3
}

// Update variable
void GalacticSpinModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void GalacticSpinModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void GalacticSpinModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute Ω_g (rad/s)
double GalacticSpinModule::computeOmega_g() {
    return variables["Omega_g"];
}

// Compute M_bh / d_g (kg/m)
double GalacticSpinModule::computeMbhOverDg() {
    return variables["M_bh"] / variables["d_g"];
}

// Compute U_bi for i (J/m^3)
double GalacticSpinModule::computeU_bi(int i) {
    std::string ug_key = "U_g" + std::to_string(i);
    double beta_i = (i == 1) ? variables["beta_1"] : 0.6;  // Assume uniform
    double U_gi = variables[ug_key];
    double Omega_g = computeOmega_g();
    double mbh_over_dg = computeMbhOverDg();
    double swirl_factor = 1.0 + variables["epsilon_sw"] * variables["rho_vac_sw"];
    double U_UA = variables["U_UA"];
    double cos_term = std::cos(variables["pi"] * variables["t_n"]);
    return -beta_i * U_gi * Omega_g * mbh_over_dg * swirl_factor * U_UA * cos_term;
}

// Equation text
std::string GalacticSpinModule::getEquationText() {
    return "U_bi = -β_i * U_gi * Ω_g * (M_bh / d_g) * (1 + ε_sw * ρ_vac,sw) * U_UA * cos(π t_n)\n"
           "Where Ω_g = 7.3e-16 rad/s (Milky Way spin rate);\n"
           "M_bh / d_g ≈3.20e16 kg/m; Ω_g * (M_bh / d_g) ≈2.34e1 kg/(m s).\n"
           "Example i=1, t_n=0: U_b1 ≈ -1.94e27 J/m³.\n"
           "Role: Introduces galactic rotation into buoyancy; stabilizes clouds/disks/nebulae.\n"
           "UQFF: Modulates counterforce to Ug; subtle over cosmic time, key for mergers.";
}

// Print variables
void GalacticSpinModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "GalacticSpinModule.h"
// int main() {
//     GalacticSpinModule mod;
//     double omega = mod.computeOmega_g();
//     std::cout << "Ω_g = " << omega << " rad/s\n";
//     double u_b1 = mod.computeU_bi(1);
//     std::cout << "U_b1 = " << u_b1 << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("Omega_g", 8e-16);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o spin_test spin_test.cpp GalacticSpinModule.cpp -lm
// Sample: Ω_g=7.3e-16 rad/s; U_b1 ≈ -1.94e27 J/m³; scales rotation in buoyancy.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Heaviside component fraction.docx
File

Encode this attachment using the template.

Thoughts
