// NegativeTimeModule.cpp
#include "NegativeTimeModule.h"

// Constructor: Set framework defaults
NegativeTimeModule::NegativeTimeModule() {
    // Universal constants
    variables["t_0"] = 0.0;                         // Reference time (s/days)
    variables["t"] = 0.0;                           // Current time
    variables["gamma"] = 5e-5;                      // day^-1 (example)
    variables["pi"] = 3.141592653589793;
    variables["mu_over_rj"] = 2.26e10;              // T m^2 (example)
    variables["P_SCm"] = 1.0;                       // Normalized
    variables["E_react"] = 1e46;                    // J
    variables["heaviside_f"] = 1e11 + 1.0;          // 1 + 10^13 * 0.01
    variables["quasi_f"] = 1.01;                    // 1 + 0.01
}

// Update variable
void NegativeTimeModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void NegativeTimeModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void NegativeTimeModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute t_n = t - t_0
double NegativeTimeModule::computeT_n(double t) {
    variables["t"] = t;
    return t - variables["t_0"];
}

// Compute cos(π t_n)
double NegativeTimeModule::computeCosPiTn(double t) {
    double t_n = computeT_n(t);
    return std::cos(variables["pi"] * t_n);
}

// Compute exp(-γ t cos(π t_n))
double NegativeTimeModule::computeExpTerm(double gamma, double t) {
    double cos_pi_tn = computeCosPiTn(t);
    double arg = - gamma * t * cos_pi_tn;
    return std::exp(arg);
}

// Compute 1 - exp(-γ t cos(π t_n))
double NegativeTimeModule::computeOneMinusExp(double gamma, double t) {
    return 1.0 - computeExpTerm(gamma, t);
}

// Simplified U_m example contrib
double NegativeTimeModule::computeUmExample(double t, double mu_over_rj) {
    double gamma = variables["gamma"];
    double one_minus_exp = computeOneMinusExp(gamma, t);
    double phi_hat = 1.0;
    double p_scm = variables["P_SCm"];
    double e_react = variables["E_react"];
    double heaviside_f = variables["heaviside_f"];
    double quasi_f = variables["quasi_f"];
    return (mu_over_rj * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
}

// Equation text
std::string NegativeTimeModule::getEquationText() {
    return "t_n = t - t_0 (s/days, allows t_n < 0 for time-reversal);\n"
           "Used in: cos(π t_n) for oscillations; exp(-γ t cos(π t_n)) for decay/growth.\n"
           "In U_m: ... (1 - exp(-γ t cos(π t_n))) ...;\n"
           "Negative t_n: e.g., t_n=-1 → cos(-π)=-1 → exp(γ t) >1 (growth, negentropic).\n"
           "Example t=1000 days, γ=5e-5 day^-1, t_0=0: 1-exp ≈0.049, U_m ≈1.12e66 J/m³.\n"
           "t_n=-1000: same (cos even); t_n=-1: 1-exp ≈ -0.051 (growth phase).\n"
           "Role: Models cyclic/TRZ dynamics; forward/reverse time in nebulae/mergers/jets.";
}

// Print variables
void NegativeTimeModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print effects for positive/negative t_n
void NegativeTimeModule::printTnEffects(double t, double gamma) {
    double t_n_pos = computeT_n(t);  // Positive example
    double cos_pos = computeCosPiTn(t);
    double exp_pos = computeExpTerm(gamma, t);
    double one_minus_pos = computeOneMinusExp(gamma, t);
    double um_pos = computeUmExample(t);

    // Negative t_n: adjust t_0 to make t_n negative
    double orig_t0 = variables["t_0"];
    variables["t_0"] = t + 1.0;  // t_n = t - (t+1) = -1
    double t_n_neg = computeT_n(t);
    double cos_neg = computeCosPiTn(t);
    double exp_neg = computeExpTerm(gamma, t);
    double one_minus_neg = computeOneMinusExp(gamma, t);
    double um_neg = computeUmExample(t);

    variables["t_0"] = orig_t0;  // Restore

    std::cout << "t_n Effects at t=" << t << " (γ=" << gamma << "):\n";
    std::cout << "Positive t_n (" << t_n_pos << "): cos(π t_n)=" << cos_pos << ", 1-exp=" << one_minus_pos << ", U_m≈" << um_pos << " J/m³\n";
    std::cout << "Negative t_n (" << t_n_neg << "): cos(π t_n)=" << cos_neg << ", 1-exp=" << one_minus_neg << ", U_m≈" << um_neg << " J/m³\n";
}

// Example usage in base program (snippet)
// #include "NegativeTimeModule.h"
// int main() {
//     NegativeTimeModule mod;
//     double t = 1000.0;  // days
//     mod.printTnEffects(t);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("t_0", 500.0);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o tn_test tn_test.cpp NegativeTimeModule.cpp -lm
// Sample: Positive: 1-exp≈0.049; Negative t_n=-1: 1-exp≈-0.051 (growth); U_m scales accordingly.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

PI mathematical constant.docx
File

Encode this attachment using the template.

Thoughts
