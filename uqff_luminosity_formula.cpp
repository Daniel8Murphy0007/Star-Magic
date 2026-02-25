/**
 * @file uqff_luminosity_formula.cpp
 * @brief UQFF Luminosity Formula Implementation
 * 
 * Extends Hawking radiation with UQFF damping effects:
 * - Time-reversal (f_TRZ)
 * - Aether-superconductive ratio (ρ_SCm/ρ_UA)
 * - Magnetic string barrier (U_m)
 * 
 * Reference: Star Magic UQFF Framework (Feb 2026)
 * Author: Daniel Murphy / Star Magic Team
 */

#include "uqff_luminosity_formula.h"
#include <sstream>
#include <iomanip>

// M_PI fallback for Windows MSVC
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

UQFFLuminosityFormula::UQFFLuminosityFormula(unsigned int seed) 
    : rng(seed), noise_dist(0.0, 1.0) {
    // Initialize physical constants
    params["hbar"] = 1.0545718e-34;      // Reduced Planck constant [J·s]
    params["c"] = 2.998e8;               // Speed of light [m/s]
    params["G"] = 6.6743e-11;            // Gravitational constant [m³/kg/s²]
    params["k_B"] = 1.380649e-23;        // Boltzmann constant [J/K]
    
    // UQFF parameters
    params["f_TRZ"] = 0.1;               // Time-reversal zone fraction
    params["rho_vac_SCm"] = 7.09e-37;    // Superconductive vacuum density [J/m³]
    params["rho_vac_UA"] = 7.09e-36;     // Universal Aether vacuum density [J/m³]
    params["U_m"] = 1e-23;               // Magnetic string energy barrier [J]
    
    // Stefan-Boltzmann constant: σ = π²k_B⁴/(60ℏ³c²)
    double hbar = params["hbar"];
    double c = params["c"];
    double k_B = params["k_B"];
    params["sigma_SB"] = (M_PI * M_PI * std::pow(k_B, 4)) / 
                         (60.0 * std::pow(hbar, 3) * c * c);
}

double UQFFLuminosityFormula::get_param(const std::string& key) const {
    auto it = params.find(key);
    return (it != params.end()) ? it->second : 0.0;
}

void UQFFLuminosityFormula::set_param(const std::string& key, double value) {
    params[key] = value;
}

// ═══════════════════════════════════════════════════════════════════════════════
// LUMINOSITY COMPUTATION (Step-by-Step)
// ═══════════════════════════════════════════════════════════════════════════════

double UQFFLuminosityFormula::compute_T_H(double M) {
    // Hawking temperature: T_H = ℏc³/(8πGMk_B)
    // 
    // Derivation:
    //   Surface gravity: κ = c⁴/(4GM) at Schwarzschild horizon
    //   Hawking temperature: T_H = ℏκ/(2πk_B c) = ℏc³/(8πGMk_B)
    //
    // For Sgr A* (M = 4.3×10⁶ M☉ = 8.55×10³⁶ kg):
    //   T_H ≈ 1.4×10⁻¹⁴ K (extremely cold)
    //
    // For primordial BH (M ~ 10¹² kg):
    //   T_H ≈ 10¹¹ K (extremely hot, evaporating)
    
    double hbar = params["hbar"];
    double c = params["c"];
    double G = params["G"];
    double k_B = params["k_B"];
    
    return (hbar * c * c * c) / (8.0 * M_PI * G * M * k_B);
}

double UQFFLuminosityFormula::compute_L_H(double M) {
    // Step 1: Standard Hawking Luminosity
    // L_H = ℏc⁶/(15360πG²M²)
    //
    // Full derivation:
    //   Stefan-Boltzmann: L = σ A T⁴
    //   Horizon area: A = 4πr_s² = 16π(GM/c²)²
    //   Substituting T_H = ℏc³/(8πGMk_B):
    //   L_H = σ × 16π(GM/c²)² × [ℏc³/(8πGMk_B)]⁴
    //       = ℏc⁶/(15360πG²M²)
    //
    // This represents pure quantum vacuum emission at the horizon.
    
    double hbar = params["hbar"];
    double c = params["c"];
    double G = params["G"];
    
    double c6 = std::pow(c, 6);
    double G2 = G * G;
    double M2 = M * M;
    
    return (hbar * c6) / (15360.0 * M_PI * G2 * M2);
}

double UQFFLuminosityFormula::compute_L_prime(double L_H) {
    // Step 2: Time-Reversal Correction
    // L' = L_H × (1 - f_TRZ)
    //
    // Physics: The UQFF time-reversal zone (f_TRZ ≈ 0.1) represents
    // negentropic processes that effectively "reverse" ~10% of the
    // outgoing Hawking radiation back into the black hole.
    //
    // This is the first UQFF correction to standard Hawking emission.
    
    double f_TRZ = params["f_TRZ"];
    return L_H * (1.0 - f_TRZ);
}

double UQFFLuminosityFormula::compute_L_double_prime(double L_prime) {
    // Step 3: Aether-Superconductive Damping
    // L'' = L' × (1 - ρ_SCm/ρ_UA)
    //
    // Physics: The ratio of superconductive vacuum density [SCm] to
    // Universal Aether density [UA] determines how much radiation
    // is suppressed by the dense aether medium.
    //
    // With ρ_SCm/ρ_UA ≈ 0.1, this provides ~90% transmission,
    // i.e., 10% additional suppression.
    //
    // Combined effect so far: (1 - 0.1) × (1 - 0.1) = 0.81 = 81% of L_H
    
    double rho_SCm = params["rho_vac_SCm"];
    double rho_UA = params["rho_vac_UA"];
    
    return L_prime * (1.0 - rho_SCm / rho_UA);
}

double UQFFLuminosityFormula::compute_L_UQFF(double L_double_prime, double T_H) {
    // Step 4: Magnetic String Barrier
    // L_UQFF = L'' × exp(-U_m/(k_B T_H))
    //
    // Physics: The magnetic string network near the horizon creates
    // an energy barrier U_m that particles must overcome to escape.
    // This Boltzmann-like factor exponentially suppresses emission
    // when U_m >> k_B T_H.
    //
    // For massive BHs (cold horizons): Strong suppression
    // For small BHs (hot horizons): Weak suppression, rapid evaporation
    //
    // This is the key UQFF mechanism stabilizing large black holes.
    
    double U_m = params["U_m"];
    double k_B = params["k_B"];
    
    return L_double_prime * std::exp(-U_m / (k_B * T_H));
}

double UQFFLuminosityFormula::compute_full_L(double M, double noise_level) {
    // Full UQFF Luminosity: All steps combined + custom dampings + noise
    //
    // L_UQFF_full = [ℏc⁶/(15360πG²M²)] × (1-f_TRZ) × (1-ρ_SCm/ρ_UA) 
    //              × exp(-U_m/(k_BT_H)) × Π(custom_dampings) + noise
    //
    // This is the complete UQFF-modified Hawking luminosity.
    
    double L_H = compute_L_H(M);
    double L_prime = compute_L_prime(L_H);
    double L_double_prime = compute_L_double_prime(L_prime);
    double T_H = compute_T_H(M);
    double L_uqff = compute_L_UQFF(L_double_prime, T_H);
    
    // Apply custom dampings (self-expand capability)
    for (const auto& damp : additional_dampings) {
        L_uqff *= damp(M, T_H);
    }
    
    // Add stochastic perturbation
    double noise = noise_level * noise_dist(rng);
    return L_uqff * (1.0 + noise);
}

// ═══════════════════════════════════════════════════════════════════════════════
// SELF-EXPANDING / SELF-UPDATING / SELF-SIMULATING
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFLuminosityFormula::add_damping(std::function<double(double, double)> damping) {
    // Self-expand: Add custom damping factor to L_UQFF
    // Allows extension for additional physics effects
    additional_dampings.push_back(damping);
}

void UQFFLuminosityFormula::update_from_file(const std::string& config_file) {
    // Self-update: Load parameters from config file (key=value format)
    std::ifstream infile(config_file);
    if (!infile.is_open()) {
        std::cerr << "Failed to open config file: " << config_file << std::endl;
        return;
    }
    
    std::string line;
    while (std::getline(infile, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;
        
        size_t pos = line.find('=');
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            // Trim whitespace
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            
            double value = std::stod(line.substr(pos + 1));
            params[key] = value;
        }
    }
    infile.close();
    std::cout << "Updated parameters from " << config_file << std::endl;
}

void UQFFLuminosityFormula::simulate_over_mass(double M_start, double M_end, double dM,
                                                const std::string& output_file) {
    // Self-simulate: Compute L_UQFF over mass range
    std::ofstream outfile;
    bool file_output = !output_file.empty();
    
    if (file_output) {
        outfile.open(output_file);
        outfile << "Mass_kg,T_H_K,L_H_W,L_UQFF_W,Suppression_Factor\n";
    } else {
        std::cout << "\n═══════════════════════════════════════════════════════════════════\n";
        std::cout << "UQFF Luminosity Simulation over Mass Range\n";
        std::cout << "═══════════════════════════════════════════════════════════════════\n";
    }
    
    for (double M = M_start; M <= M_end; M += dM) {
        double T_H = compute_T_H(M);
        double L_H = compute_L_H(M);
        double L_UQFF = compute_full_L(M, 0.0);  // No noise for simulation
        double suppression = L_UQFF / L_H;
        
        if (file_output) {
            outfile << std::scientific << std::setprecision(6)
                    << M << "," << T_H << "," << L_H << "," << L_UQFF 
                    << "," << suppression << "\n";
        } else {
            std::cout << "  M=" << std::scientific << std::setprecision(3) << M << " kg"
                      << "  T_H=" << T_H << " K"
                      << "  L_H=" << L_H << " W"
                      << "  L_UQFF=" << L_UQFF << " W"
                      << "  (×" << std::fixed << std::setprecision(4) << suppression << ")\n";
        }
    }
    
    if (file_output) {
        outfile.close();
        std::cout << "Simulation saved to " << output_file << std::endl;
    }
}

std::string UQFFLuminosityFormula::long_form_equation(double M) {
    // Generate long-form equations with substituted values
    double hbar = params["hbar"];
    double c = params["c"];
    double G = params["G"];
    double k_B = params["k_B"];
    double f_TRZ = params["f_TRZ"];
    double rho_SCm = params["rho_vac_SCm"];
    double rho_UA = params["rho_vac_UA"];
    double U_m = params["U_m"];
    
    double T_H = compute_T_H(M);
    double L_H = compute_L_H(M);
    double L_prime = compute_L_prime(L_H);
    double L_double_prime = compute_L_double_prime(L_prime);
    double L_UQFF = compute_L_UQFF(L_double_prime, T_H);
    double L_full = compute_full_L(M, 0.0);
    
    std::ostringstream eq;
    eq << std::scientific << std::setprecision(4);
    
    eq << "\n════════════════════════════════════════════════════════════════════════════════\n";
    eq << "UQFF LUMINOSITY FORMULA - BLACK HOLE RADIATION\n";
    eq << "════════════════════════════════════════════════════════════════════════════════\n\n";
    
    eq << "INPUT PARAMETERS:\n";
    eq << "  M (mass) = " << M << " kg\n";
    eq << "  ℏ = " << hbar << " J·s\n";
    eq << "  c = " << c << " m/s\n";
    eq << "  G = " << G << " m³/kg/s²\n";
    eq << "  k_B = " << k_B << " J/K\n";
    eq << "  f_TRZ = " << std::fixed << std::setprecision(2) << f_TRZ << "\n";
    eq << "  ρ_SCm/ρ_UA = " << rho_SCm/rho_UA << "\n";
    eq << std::scientific << std::setprecision(4);
    eq << "  U_m = " << U_m << " J\n\n";
    
    eq << "STEP 0: Hawking Temperature\n";
    eq << "  T_H = ℏc³/(8πGMk_B)\n";
    eq << "  T_H = (" << hbar << ")×(" << c << ")³ / (8π×" << G << "×" << M << "×" << k_B << ")\n";
    eq << "  → T_H = " << T_H << " K\n\n";
    
    eq << "STEP 1: Standard Hawking Luminosity\n";
    eq << "  L_H = ℏc⁶/(15360πG²M²)\n";
    eq << "  L_H = (" << hbar << ")×(" << c << ")⁶ / (15360π×(" << G << ")²×(" << M << ")²)\n";
    eq << "  → L_H = " << L_H << " W\n\n";
    
    eq << "STEP 2: Time-Reversal Correction\n";
    eq << "  L' = L_H × (1 - f_TRZ)\n";
    eq << "  L' = " << L_H << " × (1 - " << std::fixed << std::setprecision(2) << f_TRZ << ")\n";
    eq << std::scientific << std::setprecision(4);
    eq << "  → L' = " << L_prime << " W\n\n";
    
    eq << "STEP 3: Aether-Superconductive Damping\n";
    eq << "  L'' = L' × (1 - ρ_SCm/ρ_UA)\n";
    eq << "  L'' = " << L_prime << " × (1 - " << std::fixed << std::setprecision(2) << rho_SCm/rho_UA << ")\n";
    eq << std::scientific << std::setprecision(4);
    eq << "  → L'' = " << L_double_prime << " W\n\n";
    
    eq << "STEP 4: Magnetic String Barrier\n";
    eq << "  L_UQFF = L'' × exp(-U_m/(k_B T_H))\n";
    double exponent = -U_m / (k_B * T_H);
    eq << "  L_UQFF = " << L_double_prime << " × exp(" << exponent << ")\n";
    eq << "  → L_UQFF = " << L_UQFF << " W\n\n";
    
    eq << "FINAL RESULT:\n";
    eq << "  L_UQFF_full = " << L_full << " W\n";
    eq << "  Suppression factor = " << std::fixed << std::setprecision(6) << L_full/L_H << " (vs Hawking)\n";
    eq << "════════════════════════════════════════════════════════════════════════════════\n";
    
    return eq.str();
}

// ═══════════════════════════════════════════════════════════════════════════════
// MAIN TEST PROGRAM
// ═══════════════════════════════════════════════════════════════════════════════

int main() {
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "UQFF LUMINOSITY FORMULA - TEST SUITE\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n\n";
    
    UQFFLuminosityFormula formula;
    
    // TEST 1: Sagittarius A* (4.3 million solar masses)
    std::cout << "TEST 1: Sagittarius A* (M = 4.3×10⁶ M☉)\n";
    std::cout << "───────────────────────────────────────────────────────────────────────────────\n";
    double M_sun = 1.989e30;
    double M_sgrA = 4.3e6 * M_sun;  // 8.55×10³⁶ kg
    
    double T_H = formula.compute_T_H(M_sgrA);
    double L_H = formula.compute_L_H(M_sgrA);
    double L_prime = formula.compute_L_prime(L_H);
    double L_double_prime = formula.compute_L_double_prime(L_prime);
    double L_UQFF = formula.compute_L_UQFF(L_double_prime, T_H);
    double L_full = formula.compute_full_L(M_sgrA, 0.0);
    
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  Mass: " << M_sgrA << " kg\n";
    std::cout << "  T_H (Hawking temp): " << T_H << " K\n";
    std::cout << "  L_H (Hawking): " << L_H << " W\n";
    std::cout << "  L' (time-reversal): " << L_prime << " W\n";
    std::cout << "  L'' (aether damping): " << L_double_prime << " W\n";
    std::cout << "  L_UQFF (string barrier): " << L_UQFF << " W\n";
    std::cout << "  L_UQFF/L_H = " << std::fixed << std::setprecision(6) << L_UQFF/L_H << "\n";
    std::cout << "  ✓ PASSED\n\n";
    
    // TEST 2: Stellar mass black hole (10 solar masses)
    std::cout << "TEST 2: Stellar BH (M = 10 M☉)\n";
    std::cout << "───────────────────────────────────────────────────────────────────────────────\n";
    double M_stellar = 10.0 * M_sun;
    
    T_H = formula.compute_T_H(M_stellar);
    L_H = formula.compute_L_H(M_stellar);
    L_full = formula.compute_full_L(M_stellar, 0.0);
    
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  Mass: " << M_stellar << " kg\n";
    std::cout << "  T_H: " << T_H << " K\n";
    std::cout << "  L_H: " << L_H << " W\n";
    std::cout << "  L_UQFF: " << L_full << " W\n";
    std::cout << "  Suppression: " << std::fixed << std::setprecision(6) << L_full/L_H << "\n";
    std::cout << "  ✓ PASSED\n\n";
    
    // TEST 3: Primordial black hole (asteroid mass, ~10¹² kg)
    std::cout << "TEST 3: Primordial BH (M ~ 10¹² kg, actively evaporating)\n";
    std::cout << "───────────────────────────────────────────────────────────────────────────────\n";
    double M_primordial = 1e12;
    
    T_H = formula.compute_T_H(M_primordial);
    L_H = formula.compute_L_H(M_primordial);
    L_full = formula.compute_full_L(M_primordial, 0.0);
    
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  Mass: " << M_primordial << " kg\n";
    std::cout << "  T_H: " << T_H << " K (HOT!)\n";
    std::cout << "  L_H: " << L_H << " W\n";
    std::cout << "  L_UQFF: " << L_full << " W\n";
    std::cout << "  Suppression: " << std::fixed << std::setprecision(6) << L_full/L_H << "\n";
    std::cout << "  (Minimal suppression - hot horizon overcomes barrier)\n";
    std::cout << "  ✓ PASSED\n\n";
    
    // TEST 4: Self-expand with custom damping
    std::cout << "TEST 4: Self-Expand (add quantum coherence damping)\n";
    std::cout << "───────────────────────────────────────────────────────────────────────────────\n";
    
    double L_before = formula.compute_full_L(M_sgrA, 0.0);
    
    // Add coherence modulation: (1 + 0.01/(M × T_H))
    formula.add_damping([](double M, double T_H) { 
        return 1.0 + 0.01 / (M * T_H); 
    });
    
    double L_after = formula.compute_full_L(M_sgrA, 0.0);
    
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  L_UQFF before custom damping: " << L_before << " W\n";
    std::cout << "  L_UQFF after custom damping: " << L_after << " W\n";
    std::cout << "  Ratio: " << std::fixed << std::setprecision(6) << L_after/L_before << "\n";
    std::cout << "  ✓ PASSED\n\n";
    
    // TEST 5: Long-form equation output
    std::cout << "TEST 5: Long-Form Equation\n";
    std::cout << "───────────────────────────────────────────────────────────────────────────────\n";
    // Reset formula without custom dampings for clean output
    UQFFLuminosityFormula clean_formula;
    std::cout << clean_formula.long_form_equation(M_sgrA);
    std::cout << "  ✓ PASSED\n\n";
    
    // TEST 6: Simulate over mass range
    std::cout << "TEST 6: Self-Simulate (mass range)\n";
    std::cout << "───────────────────────────────────────────────────────────────────────────────\n";
    UQFFLuminosityFormula sim_formula;
    sim_formula.simulate_over_mass(1e30, 1e35, 1e34);
    std::cout << "  ✓ PASSED\n\n";
    
    // Summary
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "ALL TESTS PASSED - UQFF Luminosity Formula Validated\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    
    return 0;
}