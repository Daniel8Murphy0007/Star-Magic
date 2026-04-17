#include "../include/uqff_equations.h"
#include <cmath>
#include <stdexcept>

UQFFEquationCatalog::UQFFEquationCatalog() {
    initializeEquations();
}

void UQFFEquationCatalog::initializeEquations() {
    // From Grok thread - Universal Buoyancy Force
    addEquation({
        "F_U_Bi_i",
        "Universal Buoyancy Force (Relativistic)",
        "F_U_Bi_i = F_rel * (E_cm / E_LEP) * Q_wave * g(r,t)",
        {"F_rel", "E_cm", "E_LEP", "Q_wave", "r", "t"},
        [](const std::map<std::string, double>& params) {
            double F_rel = params.at("F_rel");          // 4.30e33 N (LEP baseline)
            double E_cm = params.at("E_cm");            // Center-of-mass energy
            double E_LEP = params.at("E_LEP");          // LEP energy (200 GeV)
            double Q_wave = params.at("Q_wave");        // Wave factor (~1e12)
            double r = params.at("r");
            double t = params.at("t");
            
            // Gravity term (simplified)
            double G = 6.67430e-11;
            double M = params.count("M") ? params.at("M") : 1.989e30;
            double B_field = params.count("B") ? params.at("B") : 1e-4;
            double mu_s = B_field * r * r * r;  // magnetic dipole moment
            double g = mu_s * G * M / (r * r);  // DPM-emergent: mu_s * grad(M_s/r)
            
            return F_rel * (E_cm / E_LEP) * Q_wave * g;
        }
    });
    
    // Universal Magnetism (from thread Iteration #32)
    addEquation({
        "Um",
        "Universal Magnetism (Time-dependent)",
        "Um(t,r,n) = mu_j(t) / r * (1 - exp(-gamma*t*cos(pi*t*n))) * phi * P_SCm * E_react(t) * (1 + 1e13*f_Heav) * (1 + f_quasi)",
        {"t", "r", "n", "gamma", "phi", "P_SCm", "f_Heav", "f_quasi"},
        [](const std::map<std::string, double>& params) {
            double t = params.at("t");
            double r = params.at("r");
            double n = params.count("n") ? params.at("n") : 0;
            double gamma = params.count("gamma") ? params.at("gamma") : 5e-5; // day^-1
            double phi = params.count("phi") ? params.at("phi") : 1.0;
            double P_SCm = params.count("P_SCm") ? params.at("P_SCm") : 1.0;
            double f_Heav = params.count("f_Heav") ? params.at("f_Heav") : 0.01;
            double f_quasi = params.count("f_quasi") ? params.at("f_quasi") : 0.01;
            
            // Magnetic moment μ_j(t) = (1000 + 0.4*sin(omega_c*t)) * 3.38e20 T pm³
            double omega_c = 1.585e-8; // rad/s
            double mu_j = (1000 + 0.4 * std::sin(omega_c * t)) * 3.38e20;
            
            // Reactivity energy
            double E_react = 1e46 * std::exp(-0.0005 * t);
            
            // Full equation
            double term1 = mu_j / r;
            double term2 = 1.0 - std::exp(-gamma * t * std::cos(M_PI * t * n));
            double term3 = phi * P_SCm * E_react;
            double term4 = (1.0 + 1e13 * f_Heav) * (1.0 + f_quasi);
            
            return term1 * term2 * term3 * term4;
        }
    });
    
    // Master Universal Gravity Equation (MUGE) - Hydrogen Atom
    addEquation({
        "g_MUGE_H",
        "MUGE for Hydrogen Atom",
        "g = G*m_eff*m_p/r^2 + sum(G*M_Z/r_Z^2 * (1 + f_sc) * exp(H_0*t/c))",
        {"r", "t", "m_eff", "f_sc"},
        [](const std::map<std::string, double>& params) {
            double r = params.at("r");
            double t = params.count("t") ? params.at("t") : 0;
            double m_eff = params.count("m_eff") ? params.at("m_eff") : 1.67e-27; // proton mass
            double f_sc = params.count("f_sc") ? params.at("f_sc") : 0.1;
            
            const double G = 6.67430e-11;
            const double m_p = 1.67262e-27;
            const double H_0 = 2.2e-18; // Hubble constant in s^-1
            const double c = 299792458;
            
            double g_base = G * m_eff * m_p / (r * r);  // DPM: atomic-scale mass gradient
            double z_correction = (1.0 + f_sc) * std::exp(H_0 * t / c);
            
            return g_base * z_correction;
        }
    });
    
    // Magnetar Gravity (from MUGE thread)
    addEquation({
        "g_Magnetar",
        "MUGE for Magnetar (Field decay)",
        "g = (G*M/r^2) * (1 + H_0*t) * (1 - B(t)/B_crit) + (Ug1 + Ug2 + Ug3 + Ug4)",
        {"r", "t", "M", "B_t"},
        [](const std::map<std::string, double>& params) {
            double r = params.at("r");
            double t = params.count("t") ? params.at("t") : 0;
            double M = params.count("M") ? params.at("M") : 3.0e30; // 1.5 M_sun
            double B_t = params.count("B_t") ? params.at("B_t") : 1e10; // Tesla
            
            const double G = 6.67430e-11;
            const double H_0 = 2.2e-18;
            const double B_crit = 4.4e13; // Tesla
            
            double B_sys = params.count("B_t") ? params.at("B_t") : 1e10;
            double mu_s = B_sys * r * r * r;  // magnetar dipole moment
            double g_base = mu_s * (G * M / (r * r)) * (1.0 + H_0 * t) * (1.0 - B_t / B_crit);  // DPM-emergent
            
            // Simplified Ug terms (would need full parameters)
            double Ug_sum = 0.0; // Placeholder
            
            return g_base + Ug_sum;
        }
    });
    
    // Sgr A* Evolution (from MUGE thread)
    addEquation({
        "g_SgrA",
        "MUGE for Sagittarius A* SMBH",
        "g = (G*M(t)/r^2) * (1 + H_0*t) * (1 - B(t)/B_crit) + Ug_terms",
        {"r", "t"},
        [](const std::map<std::string, double>& params) {
            double r = params.at("r");
            double t = params.count("t") ? params.at("t") : 0;
            
            const double G = 6.67430e-11;
            const double M_0 = 4.3e6 * 1.989e30; // 4.3 million solar masses
            const double H_0 = 2.2e-18;
            const double t_age = 9e9 * 365.25 * 86400; // 9 Gyr in seconds
            
            // Mass growth: M(t) = M_0 * (1 + M_dot * exp(-t/t_age))
            double M_dot = 0.001; // growth rate
            double M_t = M_0 * (1.0 + M_dot * std::exp(-t / t_age));
            
            // DPM-emergent (B negligible for SMBH; mu_s ~ B*r^3)
            double B_smbh = 1e-4;  // SMBH accretion disk B [T]
            double mu_s = B_smbh * r * r * r;
            double g = mu_s * (G * M_t / (r * r)) * (1.0 + H_0 * t);  // DPM-emergent
            
            return g;
        }
    });
    
    // Clustering Probability (from thread Batch 20)
    addEquation({
        "P_alpha",
        "Alpha Clustering Probability (Buoyancy-stabilized)",
        "P_alpha = 1 - exp(-|F_U_Bi_i| / E_th)",
        {"F_U_Bi_i", "E_th"},
        [](const std::map<std::string, double>& params) {
            double F = params.at("F_U_Bi_i");
            double E_th = params.count("E_th") ? params.at("E_th") : 8e-13; // 5 MeV in Joules
            
            return 1.0 - std::exp(-std::abs(F) / E_th);
        }
    });
    
    // EU Ratio (Electric Universe vs Gravity)
    addEquation({
        "R_EU",
        "Electric Universe Dominance Ratio",
        "R = (q * Um * rho_vac * v / r) / (G * M * m / r^2)",
        {"q", "Um", "rho_vac", "v", "r", "M", "m"},
        [](const std::map<std::string, double>& params) {
            double q = params.at("q"); // charge
            double Um = params.at("Um"); // magnetic field
            double rho_vac = params.at("rho_vac"); // vacuum density
            double v = params.at("v"); // velocity
            double r = params.at("r"); // distance
            double M = params.at("M"); // mass 1
            double m = params.at("m"); // mass 2
            
            const double G = 6.67430e-11;
            
            double F_EM = q * (Um * rho_vac / r) * v;
            double F_g = G * M * m / (r * r);
            
            return F_EM / F_g;
        }
    });
    
    // Gyro Torque Nullification
    addEquation({
        "tau_gyro",
        "Gyroscopic Torque (Inertial Stabilization)",
        "tau = I * omega * alpha",
        {"I", "omega", "alpha"},
        [](const std::map<std::string, double>& params) {
            double I = params.at("I");     // moment of inertia
            double omega = params.at("omega"); // angular velocity
            double alpha = params.at("alpha"); // precession rate
            
            return I * omega * alpha;
        }
    });
    
    // Compressed Gravity (26-layer framework from source115.cpp)
    addEquation{
        "g_compressed",
        "26-Layer Compressed Gravity",
        "g(r,t) = sum(i=1 to 26) [Ug1_i + Ug2_i + Ug3_i + Ug4_i] * Q_i * [UA]_i * [SCm]_i",
        {"r", "t", "n"},
        [](const std::map<std::string, double>& params) {
            double r = params.at("r");
            double t = params.count("t") ? params.at("t") : 0;
            int n = params.count("n") ? static_cast<int>(params.at("n")) : 26;
            
            double g_total = 0.0;
            
            // Simplified 26-layer summation
            for (int i = 1; i <= n; ++i) {
                // Each layer contributes with quantum state factors
                double Q_i = std::pow(0.99, i);      // Quantum decoherence
                double UA_i = 7.09e-36 * std::pow(0.1, i); // Aether density
                double SCm_i = 7.09e-37 * std::pow(0.1, i); // SCm density
                
                // Simplified Ug terms (would need full calculation)
                double Ug_sum = 1e-10 * std::exp(-0.1 * i);
                
                g_total += Ug_sum * Q_i * UA_i * SCm_i;
            }
            
            return g_total;
        }
    });
    
    // Widom-Larsen LENR Neutron Rate (from thread K_n document)
    addEquation({
        "eta_LENR",
        "Widom-Larsen Neutron Production Rate",
        "eta = k_eta * exp(-[SSq]^n/26) * exp(-pi - t) * Um / rho_vac_UA",
        {"k_eta", "SSq", "n", "t", "Um", "rho_vac_UA"},
        [](const std::map<std::string, double>& params) {
            double k_eta = params.count("k_eta") ? params.at("k_eta") : 1e-113; // calibration
            double SSq = params.count("SSq") ? params.at("SSq") : 0.57;
            double n = params.count("n") ? params.at("n") : 0;
            double t = params.at("t");
            double Um = params.at("Um");
            double rho_vac_UA = params.count("rho_vac_UA") ? params.at("rho_vac_UA") : 7.09e-36;
            
            double exp_term1 = std::exp(-std::pow(SSq, n) / 26.0);
            double exp_term2 = std::exp(-M_PI - t);
            
            return k_eta * exp_term1 * exp_term2 * Um / rho_vac_UA;
        }
    });
}

void UQFFEquationCatalog::addEquation(const UQFFEquation& eq) {
    equations_[eq.name] = eq;
}

const UQFFEquation* UQFFEquationCatalog::getEquation(const std::string& name) const {
    auto it = equations_.find(name);
    return (it != equations_.end()) ? &it->second : nullptr;
}

std::vector<std::string> UQFFEquationCatalog::listEquations() const {
    std::vector<std::string> names;
    for (const auto& [name, eq] : equations_) {
        names.push_back(name);
    }
    return names;
}

std::vector<std::string> UQFFEquationCatalog::searchEquations(const std::string& keyword) const {
    std::vector<std::string> results;
    
    for (const auto& [name, eq] : equations_) {
        if (name.find(keyword) != std::string::npos ||
            eq.description.find(keyword) != std::string::npos ||
            eq.latex.find(keyword) != std::string::npos) {
            results.push_back(name);
        }
    }
    
    return results;
}

double UQFFEquationCatalog::calculate(
    const std::string& equationName,
    const std::map<std::string, double>& parameters) const
{
    const UQFFEquation* eq = getEquation(equationName);
    if (!eq) {
        throw std::runtime_error("Equation not found: " + equationName);
    }
    
    // Validate all required parameters are present
    for (const auto& paramName : eq->parameters) {
        if (parameters.find(paramName) == parameters.end()) {
            throw std::runtime_error("Missing parameter: " + paramName);
        }
    }
    
    return eq->calculator(parameters);
}
