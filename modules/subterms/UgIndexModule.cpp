// UgIndexModule.cpp
#include "UgIndexModule.h"

// Constructor: Set defaults for Sun at t=0
UgIndexModule::UgIndexModule() {
    // Coupling constants (unitless)
    k_values = {1.5, 1.2, 1.8, 1.0};               // k1 to k4

    // U_gi defaults (J/m^3, Sun t=0)
    variables["U_g1"] = 1.39e26;                    // Internal Dipole
    variables["U_g2"] = 1.18e53;                    // Outer Field Bubble
    variables["U_g3"] = 1.8e49;                     // Magnetic Strings Disk
    variables["U_g4"] = 2.50e-20;                   // Star-Black Hole Interactions

    // Shared params (placeholders)
    variables["t_n"] = 0.0;                         // s
    variables["pi"] = 3.141592653589793;
}

// Update variable
void UgIndexModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
}

// Add delta
void UgIndexModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void UgIndexModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Get range of i
int UgIndexModule::getIndexRange() {
    return 4;  // i=1 to 4
}

// Compute U_gi
double UgIndexModule::computeU_gi(int i) {
    std::string key = "U_g" + std::to_string(i);
    if (variables.find(key) != variables.end()) {
        return variables[key];
    }
    std::cerr << "U_g" << i << " not found. Returning 0." << std::endl;
    return 0.0;
}

// Compute k_i (1-based)
double UgIndexModule::computeK_i(int i) {
    if (i < 1 || i > 4) {
        std::cerr << "Invalid i: " << i << ". Using k1." << std::endl;
        return k_values[0];
    }
    return k_values[i-1];
}

// Compute k_i * U_gi
double UgIndexModule::computeKUgi(int i) {
    return computeK_i(i) * computeU_gi(i);
}

// Sum over i_min to i_max
double UgIndexModule::computeSumKUgi(int i_min, int i_max) {
    double sum = 0.0;
    for (int i = i_min; i <= i_max; ++i) {
        sum += computeKUgi(i);
    }
    return sum;
}

// Equation text
std::string UgIndexModule::getEquationText() {
    return "F_U = ∑_{i=1}^4 [k_i * U_gi(r,t,M_s,ω_s,T_s,B_s,ρ_vac,[SCm],ρ_vac,[UA],t_n) - β_i * ... ] + other terms\n"
           "i (dimensionless integer): Labels Ug ranges; i=1: Internal Dipole, i=2: Outer Bubble,\n"
           "i=3: Magnetic Disk, i=4: Star-BH.\n"
           "Discretizes gravity for summation; enables scale-specific modeling.\n"
           "Example Sun t=0: ∑ k_i U_gi ≈1.42e53 J/m³ (Ug2 dominant).\n"
           "Role: Structures Ug contributions; extensible for more ranges.";
}

// Print variables
void UgIndexModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
    std::cout << "k_i: k1=1.5, k2=1.2, k3=1.8, k4=1.0\n";
}

// Print breakdown
void UgIndexModule::printIndexBreakdown() {
    std::cout << "Ug Index Breakdown (i=1 to 4):\n";
    for (int i = 1; i <= 4; ++i) {
        double ugi = computeU_gi(i);
        double ki = computeK_i(i);
        double kugi = computeKUgi(i);
        std::string label;
        switch(i) {
            case 1: label = "Internal Dipole"; break;
            case 2: label = "Outer Field Bubble"; break;
            case 3: label = "Magnetic Strings Disk"; break;
            case 4: label = "Star-Black Hole"; break;
            default: label = "Unknown";
        }
        std::cout << "i=" << i << " (" << label << "): U_g" << i << "=" << std::scientific << ugi
                  << ", k" << i << "=" << ki << ", k_i U_gi=" << kugi << " J/m³\n";
    }
    std::cout << "Sum ∑ k_i U_gi = " << std::scientific << computeSumKUgi() << " J/m³\n";
}

// Example usage in base program (snippet)
// #include "UgIndexModule.h"
// int main() {
//     UgIndexModule mod;
//     double sum = mod.computeSumKUgi();
//     std::cout << "∑ k_i U_gi = " << sum << " J/m³\n";
//     mod.printIndexBreakdown();
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("U_g3", 2e49);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ug_index_test ug_index_test.cpp UgIndexModule.cpp -lm
// Sample: Sum ≈1.42e53 J/m³; i structures 4 Ug ranges.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Inertia coupling constant.docx
File

Encode this attachment using the template.

Thoughts
