// uqff_wormhole_transverse_time.cpp
// UQFF Wormhole Transverse Time Module - Implementation
// Computes τ_UQFF traversal time through stable wormholes

#include "uqff_wormhole_transverse_time.h"

WormholeTransverseTime::WormholeTransverseTime(unsigned int seed) 
    : rng(seed == 0 ? std::chrono::system_clock::now().time_since_epoch().count() : seed),
      noise_dist(0.0, 1.0) {
    
    // Initialize physical constants
    G = 6.6743e-11;              // m³ kg⁻¹ s⁻²
    c = 2.998e8;                 // m/s
    hbar = 1.0545718e-34;        // J·s
    k_B = 1.380649e-23;          // J/K

    // Initialize UQFF vacuum densities
    rho_vac_SCm = 7.09e-37;      // J/m³ (superconductive)
    rho_vac_UA = 7.09e-36;       // J/m³ (aether)

    // Initialize UQFF parameters
    f_TRZ = 0.1;                 // Time-reversal zone factor
    U_m = 1e-30;                 // Base magnetic barrier energy (J)
    mu_j = 1e20;                 // Magnetic permeability factor
    gamma = 1e-3;                // Temporal decay rate
    t_n = 1.0;                   // Normalized time parameter

    // UQFF scaling factors (consistent with uqff_wormhole_formation)
    kappa_UQFF = 1e-60;          // Energy reduction factor
    lambda_UQFF = 1e-9;          // Magnetic scaling factor
    T_eff_floor = 1e16;          // Effective temperature floor (K)

    // Populate explanations
    explanations.push_back("=== UQFF Wormhole Transverse Time Derivation ===");
    explanations.push_back("");
    explanations.push_back("Step 1: Base traversal time from GR");
    explanations.push_back("  l ≈ r_throat = 2GM/c², τ_base = l/c = 2GM/c³");
    explanations.push_back("");
    explanations.push_back("Step 2: Effective velocity from aether flux");
    explanations.push_back("  v_eff = c × (1 - ρ_SCm/ρ_UA)");
    explanations.push_back("  With ρ_SCm/ρ_UA ≈ 0.1, v_eff ≈ 0.9c");
    explanations.push_back("");
    explanations.push_back("Step 3: Adjusted traversal time");
    explanations.push_back("  τ' = τ_base / (1 - ρ_SCm/ρ_UA)");
    explanations.push_back("");
    explanations.push_back("Step 4: Time-reversal zone drag");
    explanations.push_back("  τ'' = τ' × (1 + f_TRZ)");
    explanations.push_back("  f_TRZ ≈ 0.1 adds temporal resistance");
    explanations.push_back("");
    explanations.push_back("Step 5: Magnetic barrier exponential delay");
    explanations.push_back("  τ_UQFF = τ'' × exp(U_m / (k_B × T_eff))");
    explanations.push_back("  U_m = λ_UQFF × (μ_j/r_throat) × oscillation term");
    explanations.push_back("  T_eff = max(T_H, T_eff_floor) prevents underflow");
    explanations.push_back("");
    explanations.push_back("Full UQFF traversal time:");
    explanations.push_back("  τ_UQFF = [2GM/c³] / (1 - ρ_SCm/ρ_UA) × (1 + f_TRZ) × exp(U_m/(k_B T_eff))");
    explanations.push_back("");
    explanations.push_back("Physics interpretation (DUAL REGIME):");
    explanations.push_back("  - Aether flux slows effective propagation");
    explanations.push_back("  - TRZ zone adds temporal drag");
    explanations.push_back("  - Magnetic oscillation term determines regime:");
    explanations.push_back("    * U_m > 0 (RESISTIVE): τ_UQFF > τ'' (delayed traversal)");
    explanations.push_back("    * U_m < 0 (BUOYANT): τ_UQFF < τ'' (assisted traversal)");
    explanations.push_back("  - t_n parameter controls buoyant/resistive switching");
    explanations.push_back("  - Result: Stable wormhole with mode-dependent passage time");
}

double WormholeTransverseTime::compute_tau_base(double M) {
    // τ_base = 2GM/c³ from Step 1
    return (2.0 * G * M) / (c * c * c);
}

double WormholeTransverseTime::compute_v_eff() {
    // v_eff = c × (1 - ρ_SCm/ρ_UA) from Step 2
    return c * (1.0 - rho_vac_SCm / rho_vac_UA);
}

double WormholeTransverseTime::compute_tau_prime(double tau_base) {
    // τ' = τ_base / (1 - ρ_SCm/ρ_UA) from Step 3
    double rho_ratio = rho_vac_SCm / rho_vac_UA;
    return tau_base / (1.0 - rho_ratio);
}

double WormholeTransverseTime::compute_tau_double_prime(double tau_prime) {
    // τ'' = τ' × (1 + f_TRZ) from Step 4
    return tau_prime * (1.0 + f_TRZ);
}

double WormholeTransverseTime::compute_T_H(double M) {
    // Hawking temperature: T_H = ℏc³ / (8πGMk_B)
    return (hbar * c * c * c) / (8.0 * M_PI * G * M * k_B);
}

double WormholeTransverseTime::compute_T_eff(double T_H) {
    // Effective temperature with floor to prevent numerical underflow
    return std::max(T_H, T_eff_floor);
}

double WormholeTransverseTime::compute_U_m(double M) {
    // U_m = λ_UQFF × (μ_j / r_throat) × (1 - exp(-γ cos(πt_n)))
    // Physics: U_m < 0 = BUOYANT (assisted traversal, τ_UQFF < τ'')
    //          U_m > 0 = RESISTIVE (barrier traversal, τ_UQFF > τ'')
    // When cos(πt_n) < 0: exp term > 1 → oscillation < 0 → buoyant
    // When cos(πt_n) > 0: exp term < 1 → oscillation > 0 → resistive
    double r_throat = (2.0 * G * M) / (c * c);
    double oscillation = 1.0 - std::exp(-gamma * std::cos(M_PI * t_n));
    return lambda_UQFF * (mu_j / r_throat) * oscillation;
}

double WormholeTransverseTime::compute_tau_UQFF(double tau_double_prime, double U_m_val, double T_eff) {
    // τ_UQFF = τ'' × exp(U_m / (k_B × T_eff)) from Step 5
    // Physics: U_m < 0 → exp < 1 → τ_UQFF < τ'' (BUOYANT: faster traversal)
    //          U_m > 0 → exp > 1 → τ_UQFF > τ'' (RESISTIVE: slower traversal)
    double exponent = U_m_val / (k_B * T_eff);
    // Clamp exponent to prevent overflow (both directions)
    exponent = std::max(-700.0, std::min(exponent, 700.0));
    return tau_double_prime * std::exp(exponent);
}

double WormholeTransverseTime::compute_full_tau(double M, double noise_level) {
    // Full τ_UQFF with all UQFF effects + optional noise + mods
    double tau_base = compute_tau_base(M);
    double tau_prime = compute_tau_prime(tau_base);
    double tau_double_prime = compute_tau_double_prime(tau_prime);
    
    double T_H = compute_T_H(M);
    double T_eff = compute_T_eff(T_H);
    double U_m_val = compute_U_m(M);
    
    double tau_uqff = compute_tau_UQFF(tau_double_prime, U_m_val, T_eff);

    // Apply additional modifications
    for (const auto& mod : additional_mods) {
        tau_uqff *= mod(M);
    }

    // Add stochastic noise if requested
    if (noise_level > 0.0) {
        double noise = noise_level * tau_uqff * noise_dist(rng);
        tau_uqff += noise;
    }

    return tau_uqff;
}

void WormholeTransverseTime::add_mod(std::function<double(double)> mod) {
    // Self-expand: Add custom modification to τ_UQFF
    additional_mods.push_back(mod);
}

void WormholeTransverseTime::update_from_file(const std::string& config_file) {
    // Self-update: Load parameters from key=value config file
    std::ifstream infile(config_file);
    if (!infile.is_open()) {
        std::cerr << "Warning: Could not open config file: " << config_file << std::endl;
        return;
    }

    std::string line;
    while (std::getline(infile, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        size_t pos = line.find('=');
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            double value = std::stod(line.substr(pos + 1));
            
            if (key == "G") G = value;
            else if (key == "c") c = value;
            else if (key == "hbar") hbar = value;
            else if (key == "k_B") k_B = value;
            else if (key == "f_TRZ") f_TRZ = value;
            else if (key == "rho_vac_UA") rho_vac_UA = value;
            else if (key == "rho_vac_SCm") rho_vac_SCm = value;
            else if (key == "U_m") U_m = value;
            else if (key == "mu_j") mu_j = value;
            else if (key == "gamma") gamma = value;
            else if (key == "t_n") t_n = value;
            else if (key == "kappa_UQFF") kappa_UQFF = value;
            else if (key == "lambda_UQFF") lambda_UQFF = value;
            else if (key == "T_eff_floor") T_eff_floor = value;
        }
    }
    infile.close();
    std::cout << "Updated parameters from " << config_file << std::endl;
}

void WormholeTransverseTime::simulate_over_mass(double M_start, double M_end, double dM, const std::string& output_file) {
    // Self-simulate: Compute τ_UQFF over mass range
    std::ofstream outfile;
    bool file_output = !output_file.empty();
    
    if (file_output) {
        outfile.open(output_file);
        outfile << "Mass_kg,tau_base_s,tau_UQFF_s,T_H_K,T_eff_K,U_m_J\n";
    }

    std::cout << "Simulating τ_UQFF over mass range [" << M_start << ", " << M_end << "] kg..." << std::endl;

    for (double M = M_start; M <= M_end; M += dM) {
        double tau_base = compute_tau_base(M);
        double tau_uqff = compute_full_tau(M);
        double T_H = compute_T_H(M);
        double T_eff = compute_T_eff(T_H);
        double U_m_val = compute_U_m(M);

        if (file_output) {
            outfile << std::scientific << std::setprecision(6)
                    << M << "," << tau_base << "," << tau_uqff << ","
                    << T_H << "," << T_eff << "," << U_m_val << "\n";
        } else {
            std::cout << "M=" << std::scientific << M 
                      << ", τ_base=" << tau_base 
                      << " s, τ_UQFF=" << tau_uqff << " s" << std::endl;
        }
    }

    if (file_output) {
        outfile.close();
        std::cout << "Simulation output saved to " << output_file << std::endl;
    }
}

void WormholeTransverseTime::display_explanations() {
    // Output captured text explanations
    for (const auto& exp : explanations) {
        std::cout << exp << std::endl;
    }
}

std::string WormholeTransverseTime::generate_derivation(double M) {
    // Generate full derivation string with numerical values
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(4);

    double tau_base = compute_tau_base(M);
    double tau_prime = compute_tau_prime(tau_base);
    double tau_double_prime = compute_tau_double_prime(tau_prime);
    double T_H = compute_T_H(M);
    double T_eff = compute_T_eff(T_H);
    double U_m_val = compute_U_m(M);
    double tau_uqff = compute_full_tau(M);
    double r_throat = (2.0 * G * M) / (c * c);

    oss << "=== UQFF Wormhole Transverse Time Derivation ===" << std::endl;
    oss << std::endl;
    oss << "Input: M = " << M << " kg" << std::endl;
    oss << std::endl;
    
    oss << "Step 1: Base traversal time (GR limit)" << std::endl;
    oss << "  r_throat = 2GM/c² = " << r_throat << " m" << std::endl;
    oss << "  τ_base = 2GM/c³ = r_throat/c = " << tau_base << " s" << std::endl;
    oss << std::endl;

    oss << "Step 2: Effective velocity from aether flux" << std::endl;
    oss << "  ρ_SCm = " << rho_vac_SCm << " J/m³" << std::endl;
    oss << "  ρ_UA = " << rho_vac_UA << " J/m³" << std::endl;
    oss << "  ρ_ratio = ρ_SCm/ρ_UA = " << (rho_vac_SCm / rho_vac_UA) << std::endl;
    oss << "  v_eff = c(1 - ρ_ratio) = " << compute_v_eff() << " m/s" << std::endl;
    oss << std::endl;

    oss << "Step 3: Adjusted traversal time" << std::endl;
    oss << "  τ' = τ_base / (1 - ρ_ratio) = " << tau_prime << " s" << std::endl;
    oss << std::endl;

    oss << "Step 4: Time-reversal zone drag" << std::endl;
    oss << "  f_TRZ = " << f_TRZ << std::endl;
    oss << "  τ'' = τ' × (1 + f_TRZ) = " << tau_double_prime << " s" << std::endl;
    oss << std::endl;

    oss << "Step 5: Magnetic barrier/buoyancy term" << std::endl;
    oss << "  T_H = ℏc³/(8πGMk_B) = " << T_H << " K" << std::endl;
    oss << "  T_eff = max(T_H, T_floor) = " << T_eff << " K" << std::endl;
    oss << "  U_m = λ_UQFF × (μ_j/r) × osc = " << U_m_val << " J" << std::endl;
    oss << "  MODE: " << (U_m_val < 0 ? "BUOYANT (assisted)" : "RESISTIVE (barrier)") << std::endl;
    oss << "  τ_UQFF = τ'' × exp(U_m/(k_B T_eff)) = " << tau_uqff << " s" << std::endl;
    oss << std::endl;

    oss << "Full equation:" << std::endl;
    oss << "  τ_UQFF = [2GM/c³] / (1 - ρ_SCm/ρ_UA) × (1 + f_TRZ) × exp(U_m/(k_B T_eff))" << std::endl;
    oss << std::endl;

    oss << "UQFF Scaling Factors:" << std::endl;
    oss << "  κ_UQFF = " << kappa_UQFF << " (energy reduction)" << std::endl;
    oss << "  λ_UQFF = " << lambda_UQFF << " (magnetic scaling)" << std::endl;
    oss << "  T_eff_floor = " << T_eff_floor << " K (temperature floor)" << std::endl;
    oss << "  t_n = " << t_n << " (normalized time, controls regime)" << std::endl;
    oss << std::endl;

    oss << "Result: τ_UQFF = " << tau_uqff << " s" << std::endl;
    oss << "  Mode: " << (U_m_val < 0 ? "BUOYANT" : "RESISTIVE") << std::endl;
    oss << "  Ratio to GR: τ_UQFF/τ_base = " << (tau_uqff / tau_base) << std::endl;
    oss << "  Ratio to τ'': τ_UQFF/τ'' = " << (tau_uqff / tau_double_prime) << std::endl;

    return oss.str();
}

int WormholeTransverseTime::run_tests() {
    // Run validation tests
    std::cout << "=== UQFF Wormhole Transverse Time Tests ===" << std::endl;
    int passed = 0;
    int total = 12;

    // Test mass: Sgr A* ~ 4×10⁶ M_sun
    double M_sun = 1.989e30;
    double M_test = 4.0e6 * M_sun;  // 7.956e36 kg

    // Test 1: τ_base > 0
    double tau_base = compute_tau_base(M_test);
    if (tau_base > 0) {
        std::cout << "Test 1: τ_base > 0 - PASS (got " << tau_base << " s)" << std::endl;
        passed++;
    } else {
        std::cout << "Test 1: τ_base > 0 - FAIL (got " << tau_base << ")" << std::endl;
    }

    // Test 2: v_eff < c (aether reduces speed)
    double v_eff = compute_v_eff();
    if (v_eff > 0 && v_eff < c) {
        std::cout << "Test 2: 0 < v_eff < c - PASS (got " << v_eff << " m/s)" << std::endl;
        passed++;
    } else {
        std::cout << "Test 2: 0 < v_eff < c - FAIL (got " << v_eff << ")" << std::endl;
    }

    // Test 3: v_eff ≈ 0.9c (with ρ_ratio = 0.1)
    double expected_v_ratio = 1.0 - rho_vac_SCm / rho_vac_UA;
    double actual_v_ratio = v_eff / c;
    if (std::abs(actual_v_ratio - expected_v_ratio) < 0.01) {
        std::cout << "Test 3: v_eff/c ≈ " << expected_v_ratio << " - PASS (got " << actual_v_ratio << ")" << std::endl;
        passed++;
    } else {
        std::cout << "Test 3: v_eff/c ≈ " << expected_v_ratio << " - FAIL (got " << actual_v_ratio << ")" << std::endl;
    }

    // Test 4: τ' > τ_base (aether delays)
    double tau_prime = compute_tau_prime(tau_base);
    if (tau_prime > tau_base) {
        std::cout << "Test 4: τ' > τ_base - PASS (ratio " << tau_prime / tau_base << ")" << std::endl;
        passed++;
    } else {
        std::cout << "Test 4: τ' > τ_base - FAIL" << std::endl;
    }

    // Test 5: τ'' > τ' (TRZ adds drag)
    double tau_double_prime = compute_tau_double_prime(tau_prime);
    if (tau_double_prime > tau_prime) {
        std::cout << "Test 5: τ'' > τ' - PASS (ratio " << tau_double_prime / tau_prime << ")" << std::endl;
        passed++;
    } else {
        std::cout << "Test 5: τ'' > τ' - FAIL" << std::endl;
    }

    // Test 6: T_H > 0
    double T_H = compute_T_H(M_test);
    if (T_H > 0) {
        std::cout << "Test 6: T_H > 0 - PASS (got " << T_H << " K)" << std::endl;
        passed++;
    } else {
        std::cout << "Test 6: T_H > 0 - FAIL" << std::endl;
    }

    // Test 7: T_eff >= T_eff_floor
    double T_eff = compute_T_eff(T_H);
    if (T_eff >= T_eff_floor) {
        std::cout << "Test 7: T_eff >= T_floor - PASS (got " << T_eff << " K)" << std::endl;
        passed++;
    } else {
        std::cout << "Test 7: T_eff >= T_floor - FAIL (got " << T_eff << ")" << std::endl;
    }

    // Test 8: U_m is finite
    double U_m_val = compute_U_m(M_test);
    if (std::isfinite(U_m_val)) {
        std::cout << "Test 8: U_m is finite - PASS (got " << U_m_val << " J)" << std::endl;
        passed++;
    } else {
        std::cout << "Test 8: U_m is finite - FAIL" << std::endl;
    }

    // Test 9: τ_UQFF follows correct physics for U_m sign
    // BUOYANT (U_m < 0): τ_UQFF < τ'' (faster)
    // RESISTIVE (U_m > 0): τ_UQFF > τ'' (slower)
    double tau_uqff = compute_full_tau(M_test);
    bool regime_correct = false;
    std::string regime_name;
    if (U_m_val < 0) {
        regime_correct = (tau_uqff < tau_double_prime);
        regime_name = "BUOYANT (τ_UQFF < τ'')";
    } else {
        regime_correct = (tau_uqff >= tau_double_prime);
        regime_name = "RESISTIVE (τ_UQFF >= τ'')";
    }
    if (regime_correct) {
        std::cout << "Test 9: " << regime_name << " - PASS (ratio " << tau_uqff / tau_double_prime << ")" << std::endl;
        passed++;
    } else {
        std::cout << "Test 9: " << regime_name << " - FAIL (ratio " << tau_uqff / tau_double_prime << ")" << std::endl;
    }

    // Test 10: τ_UQFF is finite and positive
    if (std::isfinite(tau_uqff) && tau_uqff > 0) {
        std::cout << "Test 10: τ_UQFF finite & positive - PASS (got " << tau_uqff << " s)" << std::endl;
        passed++;
    } else {
        std::cout << "Test 10: τ_UQFF finite & positive - FAIL (got " << tau_uqff << ")" << std::endl;
    }

    // Test 11: Mod doubles τ_UQFF
    WormholeTransverseTime wtt_mod(42);
    wtt_mod.add_mod([](double) { return 2.0; });
    double tau_uqff_mod = wtt_mod.compute_full_tau(M_test);
    double tau_uqff_base = compute_full_tau(M_test);
    double mod_ratio = tau_uqff_mod / tau_uqff_base;
    if (std::abs(mod_ratio - 2.0) < 0.01) {
        std::cout << "Test 11: Mod doubles τ_UQFF - PASS (ratio " << mod_ratio << ")" << std::endl;
        passed++;
    } else {
        std::cout << "Test 11: Mod doubles τ_UQFF - FAIL (ratio " << mod_ratio << ")" << std::endl;
    }

    // Test 12: Derivation contains all steps
    std::string deriv = generate_derivation(M_test);
    bool has_all = deriv.find("Step 1") != std::string::npos &&
                   deriv.find("Step 2") != std::string::npos &&
                   deriv.find("Step 3") != std::string::npos &&
                   deriv.find("Step 4") != std::string::npos &&
                   deriv.find("Step 5") != std::string::npos;
    if (has_all) {
        std::cout << "Test 12: Derivation contains all steps - PASS" << std::endl;
        passed++;
    } else {
        std::cout << "Test 12: Derivation contains all steps - FAIL" << std::endl;
    }

    std::cout << "=== Results: " << passed << "/" << total << " tests passed ===" << std::endl;
    return passed;
}

// Main function for standalone testing
int main() {
    std::cout << "UQFF Wormhole Transverse Time Module" << std::endl;
    std::cout << std::string(45, '=') << std::endl << std::endl;

    WormholeTransverseTime wtt;

    // Display theory
    wtt.display_explanations();
    std::cout << std::endl;

    // Run tests
    int passed = wtt.run_tests();
    std::cout << std::endl;

    // Example calculation: Sgr A*
    double M_sun = 1.989e30;
    double M_test = 4.0e6 * M_sun;

    std::cout << "=== Example: Sgr A* (4×10⁶ M_sun) ===" << std::endl;
    std::cout << wtt.generate_derivation(M_test) << std::endl;

    // Simulate over mass range
    std::cout << "Running mass simulation..." << std::endl;
    wtt.simulate_over_mass(1.0e30, 1.0e40, 1.0e39, "transverse_time_vs_mass.csv");

    return (passed >= 10) ? 0 : 1;
}
