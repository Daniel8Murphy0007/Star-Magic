// Source165.cpp - UQFFBuoyancyModule
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration)
// for Buoyancy Equations across M74, Eagle Nebula (M16), M84, Centaurus A, Supernova Survey.
// All variables stored in std::map for dynamic addition/subtraction/update, using complex<double>.
// Nothing is negligible: Includes all terms - base force, momentum, gravity, vacuum stability,
// LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino.
// Approximations: Integral approximated as integrand * x2 (quadratic root); imag parts small and not fully scaled;
// LENR dominant due to low ω_0; x2 from quadratic solver approx; F_rel from 1998 LEP.
// Multi-system params: M74 M=7.17e41 kg r=9.46e20 m; M16 M=1e36 kg r=2.36e17 m;
// M84 M=1.46e45 kg r=3.09e22 m; Centaurus A M=4e41 kg r=3.09e21 m; Supernova Survey M=1e30 kg r=1e10 m.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 22, 2025.

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <sstream>
#include <vector>
#include <fstream>

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

using cdouble = std::complex<double>;

// ============================================================================
// HEADER SECTION
// ============================================================================

class UQFFBuoyancyModule
{
private:
    std::map<std::string, cdouble> variables;
    cdouble computeIntegrand(double t, const std::string &system);
    cdouble computeDPM_resonance(const std::string &system);
    cdouble computeX2(const std::string &system);
    cdouble computeQuadraticRoot(cdouble a, cdouble b, cdouble c);
    cdouble computeLENRTerm(const std::string &system);
    double computeG(double t, const std::string &system);
    void setSystemParams(const std::string &system);

public:
    // Constructor: Initialize all variables with multi-system defaults
    UQFFBuoyancyModule();

    // Dynamic variable operations (complex)
    void updateVariable(const std::string &name, cdouble value);
    void addToVariable(const std::string &name, cdouble delta);
    void subtractFromVariable(const std::string &name, cdouble delta);

    // Core computation: Full F_U_Bi_i(r, t) for system (approx integral)
    cdouble computeFBi(const std::string &system, double t);

    // Sub-equations
    cdouble computeCompressed(const std::string &system, double t); // Integrand
    cdouble computeResonant(const std::string &system);
    cdouble computeBuoyancy(const std::string &system);
    cdouble computeSuperconductive(const std::string &system, double t);
    double computeCompressedG(const std::string &system, double t); // g(r,t)

    // Output descriptive text of the equation
    std::string getEquationText(const std::string &system);

    // Print all current variables (for debugging/updates)
    void printVariables();

    // Export state for integration
    void exportState(const std::string &filename);
};

// ============================================================================
// IMPLEMENTATION SECTION
// ============================================================================

// Constructor: Set all variables with multi-system defaults
UQFFBuoyancyModule::UQFFBuoyancyModule()
{
    cdouble zero = {0.0, 0.0};
    cdouble i_small = {0.0, 1e-37};
    (void)zero;    // Suppress unused warning
    (void)i_small; // Suppress unused warning

    // Base constants (universal)
    variables["G"] = {6.6743e-11, 0.0};
    variables["c"] = {3e8, 0.0};
    variables["hbar"] = {1.0546e-34, 0.0};
    variables["q"] = {1.6e-19, 0.0};
    variables["pi"] = {M_PI, 0.0};
    variables["m_e"] = {9.11e-31, 0.0};
    variables["mu_B"] = {9.274e-24, 0.0};
    variables["g_Lande"] = {2.0, 0.0};
    variables["k_B"] = {1.38e-23, 0.0};
    variables["mu0"] = {4 * M_PI * 1e-7, 0.0};

    // Shared params
    variables["F_rel"] = {4.30e33, 0.0}; // Relativistic coherence from LEP 1998
    variables["F0"] = {1.83e71, 0.0};
    variables["V"] = {1e-3, 0.0};         // Default particle velocity
    variables["theta"] = {M_PI / 4, 0.0}; // 45 deg
    variables["phi"] = {M_PI / 4, 0.0};
    variables["omega_act"] = {2 * M_PI * 300, 0.0};
    variables["k_act"] = {1e-6, 0.0};
    variables["k_DE"] = {1e-30, 0.0};
    variables["k_neutron"] = {1e10, 0.0};
    variables["sigma_n"] = {1e-4, 0.0};
    variables["k_rel"] = {1e-10, 0.0};
    variables["E_cm_astro"] = {1.24e24, 0.0};
    variables["E_cm"] = {3.0264e-8, 0.0}; // 189 GeV in J
    variables["F_neutrino"] = {9.07e-43, 1e-43};
    variables["k_LENR"] = {1e-10, 0.0};
    variables["omega_LENR"] = {2 * M_PI * 1.25e12, 0.0};
    variables["rho_vac_UA"] = {7.09e-36, 1e-37};
    variables["DPM_momentum"] = {0.93, 0.05};
    variables["DPM_gravity"] = {1.0, 0.1};
    variables["DPM_stability"] = {0.01, 0.001};
    variables["beta_i"] = {0.6, 0.0};
    variables["V_infl_UA"] = {1e-6, 1e-7};
    variables["rho_vac_A"] = {1e-30, 1e-31};
    variables["a_universal"] = {1e12, 1e11};
    variables["lambda_i"] = {1.0, 0.0};
    variables["rho_vac_SCm"] = {7.09e-37, 1e-38};
    variables["omega_s"] = {2.5e-6, 1e-7};
    variables["f_TRZ"] = {0.1, 0.0};
    variables["t_scale"] = {1e16, 0.0};

    // Quadratic approx
    variables["x2"] = {-1.35e172, 0.0};
}

// Set system-specific params
void UQFFBuoyancyModule::setSystemParams(const std::string &system)
{
    if (system == "M74")
    {
        variables["M"] = {7.17e41, 0.0};
        variables["r"] = {9.46e20, 0.0};
        variables["L_X"] = {1e35, 0.0};
        variables["B0"] = {1e-9, 0.0};
        variables["rho_gas"] = {1e-24, 0.0};
        variables["t"] = {1e16, 0.0};
        variables["omega0"] = {1e-15, 0.0};
    }
    else if (system == "M16")
    {
        variables["M"] = {1e36, 0.0};
        variables["r"] = {2.36e17, 0.0};
        variables["L_X"] = {1e32, 0.0};
        variables["B0"] = {1e-5, 0.0};
        variables["rho_gas"] = {1e-21, 0.0};
        variables["t"] = {1e13, 0.0};
        variables["omega0"] = {1e-12, 0.0};
    }
    else if (system == "M84")
    {
        variables["M"] = {1.46e45, 0.0};
        variables["r"] = {3.09e22, 0.0};
        variables["L_X"] = {1e38, 0.0};
        variables["B0"] = {1e-10, 0.0};
        variables["rho_gas"] = {1e-24, 0.0};
        variables["t"] = {2.21e16, 0.0};
        variables["omega0"] = {1e-15, 0.0};
    }
    else if (system == "CentaurusA")
    {
        variables["M"] = {4e41, 0.0};
        variables["r"] = {3.09e21, 0.0};
        variables["L_X"] = {1e35, 0.0};
        variables["B0"] = {1e-5, 0.0};
        variables["rho_gas"] = {1e-24, 0.0};
        variables["t"] = {3.156e14, 0.0};
        variables["omega0"] = {1e-15, 0.0};
    }
    else if (system == "SupernovaSurvey")
    {
        variables["M"] = {1e30, 0.0};
        variables["r"] = {1e10, 0.0};
        variables["L_X"] = {1e40, 0.0};
        variables["B0"] = {1e-6, 0.0};
        variables["rho_gas"] = {1e-20, 0.0};
        variables["t"] = {1e7, 0.0};
        variables["omega0"] = {1e-12, 0.0};
    }
}

// Update variable (set to new complex value)
void UQFFBuoyancyModule::updateVariable(const std::string &name, cdouble value)
{
    variables[name] = value;
}

// Add delta (complex) to variable
void UQFFBuoyancyModule::addToVariable(const std::string &name, cdouble delta)
{
    if (variables.find(name) != variables.end())
    {
        variables[name] += delta;
    }
    else
    {
        variables[name] = delta;
    }
}

// Subtract delta (complex)
void UQFFBuoyancyModule::subtractFromVariable(const std::string &name, cdouble delta)
{
    addToVariable(name, -delta);
}

// Compute DPM_resonance (Zeeman splitting)
cdouble UQFFBuoyancyModule::computeDPM_resonance(const std::string & /* system */)
{
    cdouble g = variables["g_Lande"];
    cdouble muB = variables["mu_B"];
    cdouble B = variables["B0"];
    cdouble hbar = variables["hbar"];
    cdouble omega0 = variables["omega0"];
    return (g * muB * B / (hbar * omega0));
}

// Compute LENR term (nuclear resonance)
cdouble UQFFBuoyancyModule::computeLENRTerm(const std::string & /* system */)
{
    cdouble k = variables["k_LENR"];
    cdouble omegaL = variables["omega_LENR"];
    cdouble omega0 = variables["omega0"];
    return k * pow(omegaL / omega0, 2.0);
}

// Compute integrand for F_U_Bi_i (11-term force equation)
cdouble UQFFBuoyancyModule::computeIntegrand(double t_user, const std::string &system)
{
    setSystemParams(system);
    variables["t"] = {t_user, 0.0};
    double cos_theta = cos(variables["theta"].real());
    double sin_theta = sin(variables["theta"].real());
    double cos_act = cos(variables["omega_act"].real() * t_user + variables["phi"].real());

    cdouble term_base = -variables["F0"];
    cdouble term_mom = (variables["m_e"] * pow(variables["c"], 2) / pow(variables["r"], 2)) * variables["DPM_momentum"] * cos_theta;
    cdouble term_grav = (variables["G"] * variables["M"] / pow(variables["r"], 2)) * variables["DPM_gravity"];
    cdouble term_vac = variables["rho_vac_UA"] * variables["DPM_stability"];
    cdouble term_LENR = computeLENRTerm(system);
    cdouble term_act = variables["k_act"] * cos_act;
    cdouble term_DE = variables["k_DE"] * variables["L_X"];
    cdouble term_res = cdouble(2.0, 0.0) * variables["q"] * variables["B0"] * variables["V"] * sin_theta * computeDPM_resonance(system);
    cdouble term_neut = variables["k_neutron"] * variables["sigma_n"];
    cdouble term_rel = variables["k_rel"] * pow(variables["E_cm_astro"] / variables["E_cm"], 2.0);
    cdouble term_neutrino = variables["F_neutrino"];

    return term_base + term_mom + term_grav + term_vac + term_LENR + term_act + term_DE + term_res + term_neut + term_rel + term_neutrino;
}

// Approx x2 (quadratic root)
cdouble UQFFBuoyancyModule::computeX2(const std::string & /* system */)
{
    return variables["x2"];
}

// Quadratic root helper
cdouble UQFFBuoyancyModule::computeQuadraticRoot(cdouble a, cdouble b, cdouble c)
{
    cdouble disc = sqrt(b * b - cdouble(4.0, 0.0) * a * c);
    return (-b - disc) / (cdouble(2.0, 0.0) * a); // Negative root approx
}

// Full F_U_Bi_i approx as integrand * x2
cdouble UQFFBuoyancyModule::computeFBi(const std::string &system, double t)
{
    cdouble integ = computeIntegrand(t, system);
    cdouble x2_val = computeX2(system);
    return integ * x2_val;
}

// Compressed (integrand)
cdouble UQFFBuoyancyModule::computeCompressed(const std::string &system, double t)
{
    return computeIntegrand(t, system);
}

// Resonant DPM
cdouble UQFFBuoyancyModule::computeResonant(const std::string &system)
{
    return computeDPM_resonance(system);
}

// Buoyancy Ub1 (inflation-driven)
cdouble UQFFBuoyancyModule::computeBuoyancy(const std::string & /* system */)
{
    cdouble beta = variables["beta_i"];
    cdouble V = variables["V_infl_UA"];
    cdouble rho = variables["rho_vac_A"];
    cdouble a = variables["a_universal"];
    return beta * V * rho * a;
}

// Superconductive Ui (time-dependent)
cdouble UQFFBuoyancyModule::computeSuperconductive(const std::string & /* system */, double t)
{
    double tn = t / variables["t_scale"].real();
    cdouble lambda = variables["lambda_i"];
    cdouble rho_sc = variables["rho_vac_SCm"];
    cdouble rho_ua = variables["rho_vac_UA"];
    cdouble omega_s = variables["omega_s"];
    double cos_term = cos(M_PI * tn);
    cdouble f_trz = variables["f_TRZ"];
    return lambda * (rho_sc / rho_ua * omega_s * cos_term * (1.0 + f_trz.real()));
}

// Compressed g(r,t) with thermal and curvature corrections
double UQFFBuoyancyModule::computeCompressedG(const std::string &system, double /* t */)
{
    setSystemParams(system);
    double G_val = variables["G"].real();
    double M_val = variables["M"].real();
    double rho = variables["rho_gas"].real();
    double r_val = variables["r"].real();
    double kB = variables["k_B"].real();
    double T_val = 1e7; // Generic temperature (K)
    double m_e_val = variables["m_e"].real();
    double c_val = variables["c"].real();
    double dpm_curv = 1e-22; // DPM curvature correction

    double term1 = -(G_val * M_val * rho) / r_val;
    double term2 = -(kB * T_val * rho) / (m_e_val * c_val * c_val);
    double term3 = dpm_curv * pow(c_val, 4) / (G_val * r_val * r_val);

    return term1 + term2 + term3;
}

// Compressed g(r,t) wrapper
double UQFFBuoyancyModule::computeG(double t, const std::string &system)
{
    return computeCompressedG(system, t);
}

// Get equation text (descriptive LaTeX)
std::string UQFFBuoyancyModule::getEquationText(const std::string &system)
{
    std::ostringstream oss;
    oss << "F_U_{Bi_i} = \\int_0^{x_2} \\left[ -F_0 + \\left( \\frac{m_e c^2}{r^2} \\right) DPM_{momentum} \\cos\\theta + ";
    oss << "\\left( \\frac{G M}{r^2} \\right) DPM_{gravity} + \\rho_{vac,[UA]} DPM_{stability} + ";
    oss << "k_{LENR} \\left( \\frac{\\omega_{LENR}}{\\omega_0} \\right)^2 + k_{act} \\cos(\\omega_{act} t + \\phi) + ";
    oss << "k_{DE} L_X + 2 q B_0 V \\sin\\theta DPM_{resonance} + k_{neutron} \\sigma_n + ";
    oss << "k_{rel} \\left( \\frac{E_{cm,astro}}{E_{cm}} \\right)^2 + F_{neutrino} \\right] dx\n";
    oss << "Approximated as (integrand * x2) for " << system << "\n";
    oss << "Buoyancy: U_{bi} = \\beta_i V_{infl} \\rho_{vac,A} a_{universal}\n";
    oss << "Superconductivity: U_i = \\lambda_i (\\rho_{SC}/\\rho_{UA}) \\omega_s \\cos(\\pi t_n) (1 + f_{TRZ})";
    return oss.str();
}

// Print variables (complex)
void UQFFBuoyancyModule::printVariables()
{
    std::cout << "Current Variables:\n";
    for (const auto &pair : variables)
    {
        std::cout << pair.first << " = " << std::scientific << std::setprecision(6)
                  << pair.second.real() << " + i*" << pair.second.imag() << std::endl;
    }
}

// Export state for integration
void UQFFBuoyancyModule::exportState(const std::string &filename)
{
    std::ofstream outfile(filename);
    outfile << "=== UQFFBuoyancyModule State Export ===\n";
    for (const auto &pair : variables)
    {
        outfile << pair.first << " = " << std::scientific << std::setprecision(10)
                << pair.second.real() << " + i*" << pair.second.imag() << "\n";
    }
    outfile.close();
    std::cout << "State exported to " << filename << std::endl;
}

// ============================================================================
// MAIN FUNCTION - DEMONSTRATION
// ============================================================================

int main()
{
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "=== Source165 - UQFFBuoyancyModule ===\n\n";

    UQFFBuoyancyModule module;

    // Test all 5 systems
    std::vector<std::string> systems = {"M74", "M16", "M84", "CentaurusA", "SupernovaSurvey"};
    double t_test = 1e12; // Test time (s)

    std::cout << "=== Multi-System UQFF Buoyancy Calculations ===\n";
    for (const auto &sys : systems)
    {
        cdouble F = module.computeFBi(sys, t_test);
        std::cout << "[" << sys << "] F_U_Bi_i(t=" << t_test << ") = "
                  << F.real() << " + i*" << F.imag() << " N\n";
    }

    std::cout << "\n=== Buoyancy Component Test (M74) ===\n";
    cdouble U_bi = module.computeBuoyancy("M74");
    std::cout << "U_{bi} (inflation buoyancy) = " << U_bi.real() << " + i*" << U_bi.imag() << "\n";

    std::cout << "\n=== Superconductivity Test (M84) ===\n";
    cdouble U_sc = module.computeSuperconductive("M84", t_test);
    std::cout << "U_i (superconductivity) = " << U_sc.real() << " + i*" << U_sc.imag() << "\n";

    std::cout << "\n=== Compressed Gravity Test (CentaurusA) ===\n";
    double g_comp = module.computeCompressedG("CentaurusA", t_test);
    std::cout << "g(r,t) compressed = " << g_comp << " m/s^2\n";

    std::cout << "\n=== Resonance Test (M16) ===\n";
    cdouble DPM_res = module.computeResonant("M16");
    std::cout << "DPM resonance = " << DPM_res.real() << " + i*" << DPM_res.imag() << "\n";

    std::cout << "\n=== Integrand Components (SupernovaSurvey) ===\n";
    cdouble integrand = module.computeCompressed("SupernovaSurvey", t_test);
    std::cout << "Integrand (11 terms) = " << integrand.real() << " + i*" << integrand.imag() << "\n";

    // Export state
    module.exportState("source165_state.txt");

    std::cout << "\n=== Module Ready for Integration ===\n";
    std::cout << "Physics covered:\n";
    std::cout << "  - Inflation-driven buoyancy (β_i × V_infl × ρ_vac × a_universal)\n";
    std::cout << "  - Time-dependent superconductivity (λ × ω_s × cos(πt_n))\n";
    std::cout << "  - Neutron scattering (k_neutron × σ_n)\n";
    std::cout << "  - Relativistic energy ratios (k_rel × (E_astro/E_cm)²)\n";
    std::cout << "  - DPM resonance (Zeeman splitting)\n";
    std::cout << "  - 11-term unified force integrand\n";

    return 0;
}
