// FeedbackFactorModule.cpp
#include "FeedbackFactorModule.h"

// Constructor: Set framework defaults
FeedbackFactorModule::FeedbackFactorModule() {
    // Universal constants
    variables["f_feedback"] = 0.1;                  // Unitless, for ΔM_BH=1 dex
    variables["delta_M_BH_dex"] = 1.0;              // 1 dex = factor 10
    variables["M_bh_initial"] = 8.15e36;            // kg (Sgr A*)
    variables["k_4"] = 1.0;                         // Coupling for Ug4
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["d_g"] = 2.55e20;                     // m
    variables["alpha"] = 0.001 / 86400.0;           // day^-1 to s^-1
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
    variables["t_n"] = 0.0;                         // s
}

// Update variable
void FeedbackFactorModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    if (name == "delta_M_BH_dex") {
        // Recalculate M_bh_final if dex changes
        computeM_bh_final();
    }
}

// Add delta
void FeedbackFactorModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void FeedbackFactorModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute f_feedback (fixed 0.1 for 1 dex)
double FeedbackFactorModule::computeF_feedback() {
    return variables["f_feedback"];
}

// Compute ΔM_BH in dex
double FeedbackFactorModule::computeDeltaM_BH() {
    return variables["delta_M_BH_dex"];
}

// Compute M_bh_final = M_bh_initial * 10^{ΔM_BH_dex}
double FeedbackFactorModule::computeM_bh_final() {
    double factor = std::pow(10.0, computeDeltaM_BH());
    double initial = variables["M_bh_initial"];
    variables["M_bh_final"] = initial * factor;
    return variables["M_bh_final"];
}

// Compute U_g4 with feedback
double FeedbackFactorModule::computeU_g4(double t, double t_n) {
    double k_4 = variables["k_4"];
    double rho_vac_SCm = variables["rho_vac_SCm"];
    double M_bh = computeM_bh_final();  // Use final mass for feedback scenario
    double d_g = variables["d_g"];
    double alpha = variables["alpha"];
    double pi = variables["pi"];
    double f_feedback = computeF_feedback();
    double exp_term = std::exp( - alpha * t );
    double cos_term = std::cos( pi * t_n );
    double feedback_factor = 1.0 + f_feedback;
    return k_4 * (rho_vac_SCm * M_bh / d_g) * exp_term * cos_term * feedback_factor;
}

// Compute U_g4 without feedback (f_feedback=0)
double FeedbackFactorModule::computeU_g4_no_feedback(double t, double t_n) {
    double orig_f = variables["f_feedback"];
    variables["f_feedback"] = 0.0;
    double result = computeU_g4(t, t_n);
    variables["f_feedback"] = orig_f;  // Restore
    return result;
}

// Equation text
std::string FeedbackFactorModule::getEquationText() {
    return "U_g4 = k_4 * (ρ_vac,[SCm] M_bh / d_g) * e^{-α t} * cos(π t_n) * (1 + f_feedback)\n"
           "Where f_feedback = 0.1 (unitless, for ΔM_BH = 1 dex = 10x mass increase);\n"
           "ΔM_BH =1 dex → M_bh_final = 10 * M_bh_initial (8.15e36 kg → 8.15e37 kg).\n"
           "Example t=0, t_n=0: U_g4 ≈2.75e-20 J/m³ (with); ≈2.50e-20 J/m³ (without; +10%).\n"
           "Role: Scales feedback in star-BH interactions; regulates AGN, mergers, star formation.\n"
           "UQFF: Enhances energy density for 10x M_BH; resolves final parsec, quasar jets.";
}

// Print variables
void FeedbackFactorModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print U_g4 comparison
void FeedbackFactorModule::printU_g4_comparison(double t, double t_n) {
    double u_with = computeU_g4(t, t_n);
    double u_without = computeU_g4_no_feedback(t, t_n);
    double delta_percent = ((u_with - u_without) / u_without) * 100.0;
    std::cout << "U_g4 Comparison at t=" << t << " s, t_n=" << t_n << " s:\n";
    std::cout << "With feedback: " << std::scientific << u_with << " J/m³\n";
    std::cout << "Without feedback: " << std::scientific << u_without << " J/m³\n";
    std::cout << "Difference: +" << std::fixed << std::setprecision(1) << delta_percent << "%\n";
}

// Example usage in base program (snippet)
// #include "FeedbackFactorModule.h"
// int main() {
//     FeedbackFactorModule mod;
//     double m_final = mod.computeM_bh_final();
//     std::cout << "M_bh_final = " << m_final << " kg\n";
//     mod.printU_g4_comparison(0.0, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("delta_M_BH_dex", 2.0);  // 100x mass
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o feedback_test feedback_test.cpp FeedbackFactorModule.cpp -lm
// Sample: M_bh_final=8.15e37 kg; U_g4 with=2.75e-20 J/m³ (+10% vs without).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

FU.docx
File

Encode this attachment using the template.

Thoughts
