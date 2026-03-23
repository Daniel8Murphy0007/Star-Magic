// ReciprocationDecayModule.cpp
#include "ReciprocationDecayModule.h"

// Constructor: Set framework defaults
ReciprocationDecayModule::ReciprocationDecayModule() {
    // Universal constants
    variables["gamma_day"] = 0.00005;               // day⁻¹
    variables["day_to_s"] = 86400.0;                // s/day
    variables["t_n"] = 0.0;                         // days
    variables["t_day"] = 0.0;                       // days
    variables["pi"] = 3.141592653589793;
    variables["mu_over_rj"] = 2.26e10;              // T m²
    variables["P_SCm"] = 1.0;                       // Normalized
    variables["E_react"] = 1e46;                    // J
    variables["heaviside_f"] = 1e11 + 1.0;          // 1 + 10^13 * 0.01
    variables["quasi_f"] = 1.01;                    // 1 + 0.01

    // Derived
    variables["gamma_s"] = computeGamma_s();
}

// Update variable
void ReciprocationDecayModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "gamma_day") {
            variables["gamma_s"] = computeGamma_s();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void ReciprocationDecayModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "gamma_day") {
            variables["gamma_s"] = computeGamma_s();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void ReciprocationDecayModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute γ (day⁻¹)
double ReciprocationDecayModule::computeGamma_day() {
    return variables["gamma_day"];
}

// Compute γ in s⁻¹
double ReciprocationDecayModule::computeGamma_s() {
    return computeGamma_day() / variables["day_to_s"];
}

// Compute cos(π t_n)
double ReciprocationDecayModule::computeCosPiTn(double t_n) {
    variables["t_n"] = t_n;
    return std::cos(variables["pi"] * t_n);
}

// Compute exp(-γ t cos(π t_n)) (t in days)
double ReciprocationDecayModule::computeExpTerm(double t_day, double t_n) {
    variables["t_day"] = t_day;
    double cos_pi_tn = computeCosPiTn(t_n);
    double arg = - computeGamma_day() * t_day * cos_pi_tn;
    return std::exp(arg);
}

// Compute 1 - exp(-γ t cos(π t_n))
double ReciprocationDecayModule::computeOneMinusExp(double t_day, double t_n) {
    return 1.0 - computeExpTerm(t_day, t_n);
}

// Simplified U_m example (J/m³)
double ReciprocationDecayModule::computeUmExample(double t_day, double t_n, double mu_over_rj) {
    double one_minus_exp = computeOneMinusExp(t_day, t_n);
    double phi_hat = 1.0;
    double p_scm = variables["P_SCm"];
    double e_react = variables["E_react"];
    double heaviside_f = variables["heaviside_f"];
    double quasi_f = variables["quasi_f"];
    return (mu_over_rj * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
}

// Equation text
std::string ReciprocationDecayModule::getEquationText() {
    return "γ = 0.00005 day⁻¹ (~5.8e-10 s⁻¹; timescale ~55 years);\n"
           "In U_m: ... (1 - exp(-γ t cos(π t_n))) ... (t days, reciprocating decay/growth).\n"
           "Negative cos(π t_n): exp(+γ t) >1 (growth, negentropic TRZ).\n"
           "Example t=1000 days, t_n=0: 1-exp ≈0.049, U_m ≈1.12e66 J/m³.\n"
           "UQFF: Slow decay for magnetic strings; cyclic via cos(π t_n) in jets/nebulae/mergers.";
}

// Print variables
void ReciprocationDecayModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print effects
void ReciprocationDecayModule::printDecayEffects(double t_day, double t_n) {
    double cos_pi = computeCosPiTn(t_n);
    double exp_val = computeExpTerm(t_day, t_n);
    double one_minus = computeOneMinusExp(t_day, t_n);
    double um_ex = computeUmExample(t_day, t_n);
    std::cout << "Decay Effects at t=" << t_day << " days, t_n=" << t_n << ":\n";
    std::cout << "cos(π t_n) = " << cos_pi << "\n";
    std::cout << "exp(-γ t cos(π t_n)) = " << exp_val << "\n";
    std::cout << "1 - exp(...) = " << one_minus << "\n";
    std::cout << "U_m example contrib = " << um_ex << " J/m³\n";
}

// Example usage in base program (snippet)
// #include "ReciprocationDecayModule.h"
// int main() {
//     ReciprocationDecayModule mod;
//     double gamma = mod.computeGamma_day();
//     std::cout << "γ = " << gamma << " day⁻¹\n";
//     mod.printDecayEffects(1000.0, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("gamma_day", 0.0001);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o decay_test decay_test.cpp ReciprocationDecayModule.cpp -lm
// Sample: γ=5e-5 day⁻¹; t=1000 days: 1-exp≈0.049; U_m≈1.12e66 J/m³.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

SCm Penetration Factor.docx
File

Encode this attachment using the template.

Thoughts

You've reached the limit of 25 attachments in this conversation. Grok might forget content from earlier attachments; please start a new conversation to upload more attachments.
