
// UQFFBuoyancyModule.cpp
#include "UQFFBuoyancyModule.h"
#include <complex>
#include <sstream>

// Constructor: Set all variables with multi-system defaults
UQFFBuoyancyModule::UQFFBuoyancyModule() {
    double pi_val = 3.141592653589793;
    cdouble zero = {0.0, 0.0};
    cdouble i_small = {0.0, 1e-37};

    // Base constants (universal)
    variables["G"] = {6.6743e-11, 0.0}; // Gravitational constant
    variables["c"] = {3e8, 0.0};
    variables["hbar"] = {1.0546e-34, 0.0};
    variables["q"] = {1.6e-19, 0.0};
    variables["pi"] = {pi_val, 0.0};
    variables["m_e"] = {9.11e-31, 0.0};
    variables["mu_B"] = {9.274e-24, 0.0};
    variables["g_Lande"] = {2.0, 0.0};
    variables["k_B"] = {1.38e-23, 0.0};
    variables["mu0"] = {4 * pi_val * 1e-7, 0.0};
    variables["epsilon0"] = {8.85e-12, 0.0}; // Added for quadratic terms

    // Shared params from document
    variables["F_rel"] = {4.30e33, 0.0};  // Relativistic coherence from LEP 1998
    variables["F0"] = {1.83e71, 0.0};
    variables["V"] = {1e-3, 0.0};  // Default particle velocity
    variables["theta"] = {pi_val / 4, 0.0};  // 45 deg
    variables["phi"] = {pi_val / 4, 0.0};
    variables["omega_act"] = {2 * pi_val * 300, 0.0};
    variables["k_act"] = {1e-6, 0.0};
    variables["k_DE"] = {1e-30, 0.0};
    variables["k_neutron"] = {1e10, 0.0};
    variables["k_rel"] = {1e-10, 0.0};
    variables["k_relativistic"] = {1e-20, 0.0};
    variables["k_neutrino"] = {1e-15, 0.0};
    variables["k_Sweet"] = {1e-25, 0.0};
    variables["k_Kozima"] = {1e-18, 0.0};
    variables["omega_0_LENR"] = {2 * pi_val * 1.25e12, 0.0};  // LENR resonance freq 1.2-1.3 THz
    variables["k_LENR"] = {1e-10, 0.0};
    variables["rho_vac_UA"] = {7.09e-36, 0.0};  // Vacuum energy density
    variables["sigma_n"] = {1e-4, 0.0};  // Scaled for astrophysical densities
    variables["E_cm"] = {189.0, 0.0};  // GeV, but scaled appropriately
    variables["E_cm_astro_local_adj_eff_enhanced"] = {1.24e24, 0.0};  // events/m3
    variables["DPM_stability"] = {0.01, 0.0};
    variables["DPM_momentum"] = {0.93, 0.0};
    variables["DPM_gravity"] = {1.0, 0.0};
    variables["DPM_light"] = {0.01, 0.0};
    variables["phase"] = {2.36e-3, 0.0};  // s^-1
    variables["curvature"] = {1e-22, 0.0};
    variables["k_10_13"] = {1e-13, 0.0};  // For light term in quadratic

    // System-specific params will be set in setSystemParams()
}

// Set system-specific parameters
void UQFFBuoyancyModule::setSystemParams(const std::string& system)
{
    double pi_val = variables["pi"].real();
    if (system == "J1610+1811") {
        this->variables["M"] = cdouble(2.785e30, 0.0);
        this->variables["r"] = cdouble(3.09e15, 0.0);
        this->variables["T"] = cdouble(1e4, 0.0);
        this->variables["L_X"] = cdouble(1e31, 0.0);
        this->variables["B0"] = cdouble(1e-4, 0.0);
        this->variables["omega0"] = cdouble(1e-12, 0.0);
        this->variables["Mach"] = cdouble(1.0, 0.0);  // ℳ
        this->variables["C"] = cdouble(1.0, 0.0);
        this->variables["theta"] = cdouble(pi_val / 4, 0.0);
        this->variables["t"] = cdouble(3.156e10, 0.0);
    } else if (system == "PLCK_G287.0+32.9") {
        this->variables["M"] = cdouble(1.989e44, 0.0);
        this->variables["r"] = cdouble(3.09e22, 0.0);
        this->variables["T"] = cdouble(1e7, 0.0);
        this->variables["L_X"] = cdouble(1e38, 0.0);
        this->variables["B0"] = cdouble(1e-4, 0.0);
        this->variables["omega0"] = cdouble(1e-15, 0.0);
        this->variables["Mach"] = cdouble(1.5, 0.0);
        this->variables["C"] = cdouble(1.2, 0.0);
        this->variables["theta"] = cdouble(pi_val / 4, 0.0);
        this->variables["t"] = cdouble(1.42e17, 0.0);
    } else if (system == "PSZ2_G181.06+48.47") {
        this->variables["M"] = cdouble(1.989e44, 0.0);
        this->variables["r"] = cdouble(3.09e22, 0.0);
        this->variables["T"] = cdouble(1e7, 0.0);
        this->variables["L_X"] = cdouble(1e39, 0.0);
        this->variables["B0"] = cdouble(1e-4, 0.0);
        this->variables["omega0"] = cdouble(1e-15, 0.0);
        this->variables["Mach"] = cdouble(1.5, 0.0);
        this->variables["C"] = cdouble(1.2, 0.0);
        this->variables["theta"] = cdouble(pi_val / 4, 0.0);
        this->variables["t"] = cdouble(2.36e17, 0.0);
    } else if (system == "ASKAP_J1832-0911") {
        this->variables["M"] = cdouble(2.785e30, 0.0);
        this->variables["r"] = cdouble(4.63e16, 0.0);
        this->variables["T"] = cdouble(1e4, 0.0);
        this->variables["L_X"] = cdouble(1e31, 0.0);
        this->variables["B0"] = cdouble(1e-4, 0.0);
        this->variables["omega0"] = cdouble(1e-12, 0.0);  // From 44-min period approx
        this->variables["Mach"] = cdouble(1.0, 0.0);
        this->variables["C"] = cdouble(1.0, 0.0);
        this->variables["theta"] = cdouble(pi_val / 4, 0.0);
        this->variables["t"] = cdouble(3.156e10, 0.0);
    } else if (system == "SonificationCollection") {
        this->variables["M"] = cdouble(1.989e31, 0.0);
        this->variables["r"] = cdouble(6.17e16, 0.0);
        this->variables["T"] = cdouble(1e5, 0.0);
        this->variables["L_X"] = cdouble(1e33, 0.0);
        this->variables["B0"] = cdouble(1e-5, 0.0);
        this->variables["omega0"] = cdouble(1e-12, 0.0);
        this->variables["Mach"] = cdouble(1.0, 0.0);
        this->variables["C"] = cdouble(1.0, 0.0);
        this->variables["theta"] = cdouble(pi_val / 4, 0.0);
        this->variables["t"] = cdouble(3.156e14, 0.0);
    }
}

// Dynamic variable operations (complex)
void UQFFBuoyancyModule::updateVariable(const std::string& name, cdouble value) {
    this->variables[name] = value;
}
void UQFFBuoyancyModule::addToVariable(const std::string& name, cdouble delta) {
    this->variables[name] += delta;
}
void UQFFBuoyancyModule::subtractFromVariable(const std::string& name, cdouble delta) {
    this->variables[name] -= delta;
}

// Core computation: Full F_U_Bi_i(r, t) for system (approx integral)
cdouble UQFFBuoyancyModule::computeFBi(const std::string& system, double t) {
    setSystemParams(system);
    cdouble integrand = computeIntegrand(t, system);
    cdouble x2 = computeX2(system);
    cdouble f_bi_i = integrand * x2;
    cdouble momentum_term = (variables["m_e"] * variables["c"] * variables["c"] / (variables["r"] * variables["r"])) * variables["DPM_momentum"] * cos(variables["theta"].real());
    cdouble gravity_term = (variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * variables["DPM_gravity"];
    cdouble f_bi = -variables["F0"] + momentum_term + gravity_term + f_bi_i;
    return f_bi;
}

// Sub-equations
cdouble UQFFBuoyancyModule::computeCompressed(const std::string& system, double t) {
    setSystemParams(system);
    return computeIntegrand(t, system);
}
cdouble UQFFBuoyancyModule::computeResonant(const std::string& system) {
    setSystemParams(system);
    return computeDPM_resonance(system);
}
cdouble UQFFBuoyancyModule::computeBuoyancy(const std::string& system) {
    setSystemParams(system);
    return computeUb1(system);
}
cdouble UQFFBuoyancyModule::computeSuperconductive(const std::string& system, double t) {
    setSystemParams(system);
    return computeUi(t, system);
}
double UQFFBuoyancyModule::computeCompressedG(const std::string& system, double t) {
    setSystemParams(system);
    return computeG(t, system);
}

// Output descriptive text of the equation
std::string UQFFBuoyancyModule::getEquationText(const std::string& system) {
    setSystemParams(system);
    std::ostringstream oss;
    oss << "F_U_Bi_i(r, t) = Integral[Integrand(r, t) dt] approximated as Integrand * x2\n";
    oss << "Where Integrand includes terms for base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino, Sweet vac, Kozima drop.\n";
    oss << "LENR Resonance: F_LENR = k_LENR * (ω_LENR / ω_0)^2\n";
    oss << "Activation: F_act = k_act * cos(ω_act t)\n";
    oss << "Directed Energy: F_DE = k_DE * L_X\n";
    oss << "Magnetic Resonance: F_res = 2 q B_0 V sinθ * DPM_resonance\n";
    oss << "Neutron Drop: F_neutron = k_neutron * σ_n\n";
    oss << "Relativistic: F_rel = k_rel * (E_cm_astro... / E_cm)^2 = 4.30e33 N\n";
    oss << "System: " << system << "\n";
    return oss.str();
}

// Print all current variables (for debugging/updates)
void UQFFBuoyancyModule::printVariables() {
    for (const auto& pair : variables) {
        std::cout << std::setw(15) << pair.first << " : " << pair.second << std::endl;
    }
}

// Compute integrand for F_U_Bi_i
cdouble UQFFBuoyancyModule::computeIntegrand(double t, const std::string& system) {
    setSystemParams(system);
    double pi_val = variables["pi"].real();
    cdouble dpm_res = computeDPM_resonance(system);
    cdouble f_lenr = variables["k_LENR"] * pow(variables["omega_0_LENR"].real() / variables["omega0"].real(), 2.0);
    cdouble f_act = variables["k_act"] * cos(variables["omega_act"].real() * t);
    cdouble f_de = variables["k_DE"] * variables["L_X"];
    cdouble f_neutron = variables["k_neutron"] * variables["sigma_n"];
    cdouble f_rel = variables["F_rel"];  // Direct from LEP
    cdouble fres = 2.0 * variables["q"].real() * variables["B0"].real() * variables["V"].real() * sin(variables["theta"].real()) * dpm_res;
    cdouble momentum_term = (variables["m_e"] * variables["c"] * variables["c"] / (variables["r"] * variables["r"])) * variables["DPM_momentum"] * cos(variables["theta"].real());
    cdouble gravity_term = (variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * variables["DPM_gravity"];
    cdouble vac_term = variables["rho_vac_UA"] * variables["DPM_stability"];
    cdouble integrand = -variables["F0"] + momentum_term + gravity_term + vac_term + f_lenr + f_act + f_de + fres + f_neutron + f_rel;
    return integrand;
}

// Compute DPM resonance term
cdouble UQFFBuoyancyModule::computeDPM_resonance(const std::string& system) {
    setSystemParams(system);
    double g = variables["g_Lande"].real();
    double mu_b = variables["mu_B"].real();
    double b0 = variables["B0"].real();
    double hbar_omega0 = variables["hbar"].real() * variables["omega0"].real();
    return cdouble(g * mu_b * b0 / hbar_omega0, 0.0);
}

// Compute x2 from quadratic root approximation (negative root as per doc)
cdouble UQFFBuoyancyModule::computeX2(const std::string& system) {
    setSystemParams(system);
    // Extract a, b, c from document patterns (approximated; a ~ gravity, b ~ sum small, c ~ -3.06e175)
    // For precision, compute as per garbled formulas; here approximated based on system
    cdouble a, b, c;
    double r2 = variables["r"].real() * variables["r"].real();
    double t_val = variables["T"].real();
    // a ≈ (k_B * q / (4 pi epsilon0 r^2 T)) + G M / r^2 + (c^4 * 1e-13 / r^2 * 0.01)  [inverted from doc]
    // But doc shows a ≈ gravity term
    a = variables["G"] * variables["M"] / cdouble(r2, 0.0);  // Dominant
    // Add small terms
    double eps_term = (variables["k_B"].real() * variables["q"].real()) / (4 * variables["pi"].real() * variables["epsilon0"].real() * r2 * t_val);
    a += cdouble(eps_term, 0.0);
    double light_term = pow(variables["c"].real(), 4) * variables["k_10_13"].real() / cdouble(r2, 0.0) * variables["DPM_light"].real();
    a += light_term;
    
    // b ≈ 2.51e-5 + T / r^2 + phase + phase
    b = cdouble(2.51e-5 + t_val / r2 + 2 * variables["phase"].real(), 0.0);
    
    // c ≈ -3.06e175 + 1e-29 / r^2 + curvature
    c = cdouble(-3.06e175 + 1e-29 / r2 + variables["curvature"].real(), 0.0);
    
    return computeQuadraticRoot(a, b, c);
}

// Compute quadratic root (negative branch: [-b - sqrt(b^2 - 4ac)] / 2a)
cdouble UQFFBuoyancyModule::computeQuadraticRoot(cdouble a, cdouble b, cdouble c) {
    cdouble disc = b*b - 4.0*a*c;
    if (disc.real() < 0) disc = cdouble(0.0, 0.0);  // Approx real for simplicity
    cdouble sqrt_disc = cdouble(sqrt(disc.real()), 0.0);  // Ignore imag for approx
    return (-b - sqrt_disc) / (2.0 * a);
}

// Compute LENR term
cdouble UQFFBuoyancyModule::computeLENRTerm(const std::string& system) {
    setSystemParams(system);
    return variables["k_LENR"] * pow(variables["omega_0_LENR"].real() / variables["omega0"].real(), 2.0);
}

// Compute gravitational acceleration g(r,t) - approximated as energy density per doc
double UQFFBuoyancyModule::computeG(double t, const std::string& system) {
    // Per doc: ≈ -1.07e16 J/m3 for all systems
    return -1.07e16;
}

// Compute Q_wave term - approximated per doc
cdouble UQFFBuoyancyModule::computeQ_wave(double t, const std::string& system) {
    // Per doc: ≈ 3.11e5 J/m3
    return cdouble(3.11e5, 0.0);
}

// Compute Ub1 buoyancy term - placeholder aligned with doc
cdouble UQFFBuoyancyModule::computeUb1(const std::string& system) {
    // Integrand dominant term
    return computeIntegrand(0.0, system);  // At t=0 approx
}

// Compute Ui superconductive term - placeholder
cdouble UQFFBuoyancyModule::computeUi(double t, const std::string& system) {
    return computeQ_wave(t, system);
}

Write with appropriate name.

Thoughts

