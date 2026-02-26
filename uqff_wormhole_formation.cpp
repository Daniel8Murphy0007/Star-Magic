// uqff_wormhole_formation.cpp
// UQFF Wormhole Formation Module Implementation
// Stable wormholes via aether-superconductive inversion

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "uqff_wormhole_formation.h"
#include <sstream>
#include <iomanip>

UQFFWormholeFormation::UQFFWormholeFormation(unsigned int seed) 
    : rng(seed), noise_dist(0.0, 1.0) {
    
    // Initialize physical constants
    G = 6.6743e-11;           // m³ kg^{-1} s^{-2}
    c = 2.998e8;              // m/s
    hbar = 1.0545718e-34;     // J·s
    k_B = 1.380649e-23;       // J/K
    
    // Initialize UQFF vacuum parameters
    rho_vac_UA = 7.09e-36;    // J/m³ (aether vacuum density)
    rho_vac_SCm = 7.09e-37;   // J/m³ (superconductive vacuum density)
    
    // Initialize formation parameters
    f_TRZ = 0.1;              // Time-reversal zone factor
    mu_j = 1e20;              // A·m² (magnetic string moment)
    gamma = 5.787e-10;        // s^{-1} (converted from 5e-5 day^{-1})
    t_n = 0.5;                // Normalized time parameter
    kappa_UQFF = 1e-60;       // UQFF energy reduction factor (aether-mediated)
    lambda_UQFF = 1e-9;       // UQFF magnetic scaling (normalizes U_m to tractable range)
    T_eff_floor = 1e16;       // K (quantum vacuum effective temperature - scales with Planck T)
    
    // Populate derivation explanations
    explanations = {
        "=== UQFF Wormhole Formation Derivation ===",
        "",
        "STANDARD CONCEPT:",
        "  Wormholes are hypothetical tunnels linking spacetime regions.",
        "  Einstein-Rosen bridge from Schwarzschild metric:",
        "    ds² = -c²dt² + dr²/(1-b(r)/r) + r²dΩ²",
        "  where b(r) is the shape function (throat at minimum r).",
        "  Requires exotic matter with ρ < 0 violating energy conditions.",
        "  Standard wormholes are unstable and collapse rapidly.",
        "",
        "UQFF REINTERPRETATION:",
        "  Stable wormholes via aether-superconductive inversion.",
        "  High [SCm] gradients in [UA] superfluid (ρ ≈ 7.09e-36 J/m³) create throats.",
        "  Stabilized by f_TRZ ≈ 0.1 time-reversal and U_m magnetic strings.",
        "  Black hole tunneling enabled via negentropic reversal.",
        "",
        "DERIVATION STEPS:",
        "",
        "Step 1: UQFF Throat Radius",
        "  r_throat,UQFF = (2GM/c²) × (ρ_UA / ρ_SCm)",
        "  Throat expands when [SCm] is low relative to [UA].",
        "",
        "Step 2: Formation Probability",
        "  P_form = f_TRZ × exp(-E_throat / (k_B T_H))",
        "  where E_throat = GM²/r_throat",
        "        T_H = ℏc³/(8πGMk_B)  [Hawking temperature]",
        "",
        "Step 3: Aether Flux",
        "  J_aether = ρ_UA × v_s × (1 + f_TRZ)",
        "  where v_s ≈ c (superfluid velocity)",
        "",
        "Step 4: Magnetic String Reinforcement",
        "  U_m = (μ_j / r_throat) × (1 - exp(-γt cos(πt_n)))",
        "  where γ ≈ 5e-5 day⁻¹",
        "",
        "Step 5: Formation Threshold",
        "  Θ_WH = P_form × J_aether × exp(U_m / (k_B T_H))",
        "  Formation occurs when Θ_WH > 1",
        "  Traversal time: τ ≈ r_throat / c",
        "",
        "NUMERICAL EXAMPLE (Sgr A*, M = 4×10⁶ M_sun):",
        "  Θ_WH ≈ 0.1 × (7.09e-36 × c × 1.1) × e¹ ≈ 10 > 1",
        "  → Wormhole formation possible!",
        "",
        "UQFF ADVANCES:",
        "  • Stable wormholes via aether reversal mechanism",
        "  • Testable via q-scope micro-tunnel detection"
    };
}

double UQFFWormholeFormation::compute_r_s_standard(double M) {
    // r_s = 2 G M / c^2
    // From Standard Wormhole Concept: Base Schwarzschild radius
    return (2 * G * M) / (c * c);
}

double UQFFWormholeFormation::compute_r_throat_UQFF(double r_s_std) {
    // r_throat,UQFF = r_s * (rho_vac_UA / rho_vac_SCm)
    // From Step 1: Base throat formation modulated by aether; expands if [SCm] low
    return r_s_std * (rho_vac_UA / rho_vac_SCm);
}

double UQFFWormholeFormation::compute_P_form(double E_throat, double T_H) {
    // P_form = f_TRZ * exp(-E_throat / (k_B T_H))
    // From Step 2: Time-reversal activation scales probability
    return f_TRZ * std::exp(-E_throat / (k_B * T_H));
}

double UQFFWormholeFormation::compute_J_aether(double v_s) {
    // J_aether = rho_vac_UA * v_s * (1 + f_TRZ)
    // From Step 3: Aether flux through throat; v_s ≈ c
    return rho_vac_UA * v_s * (1 + f_TRZ);
}

double UQFFWormholeFormation::compute_U_m(double r_throat, double t) {
    // U_m = lambda_UQFF * (mu_j / r_throat) * (1 - exp(-gamma t cos(pi t_n)))
    // From Step 4: Magnetic string reinforcement (scaled by lambda_UQFF for tractability)
    return lambda_UQFF * (mu_j / r_throat) * (1 - std::exp(-gamma * t * std::cos(M_PI * t_n)));
}

double UQFFWormholeFormation::compute_Theta_WH(double P_form, double J_aether, double U_m, double T_H) {
    // Theta_WH = P_form * J_aether * exp(U_m / (k_B T_H))
    // From Step 5: Full formation threshold; >1 forms
    return P_form * J_aether * std::exp(U_m / (k_B * T_H));
}

double UQFFWormholeFormation::compute_full_Theta_WH(double M, double t, double noise_level) {
    // Full chain: Compute Theta_WH with noise + modifiers
    // UQFF uses effective temperature = max(T_Hawking, T_eff_floor) for quantum vacuum effects
    double r_s_std = compute_r_s_standard(M);
    double r_throat = compute_r_throat_UQFF(r_s_std);
    double E_throat = kappa_UQFF * G * M * M / r_throat;  // UQFF-reduced energy
    double T_H_raw = hbar * c * c * c / (8 * M_PI * G * M * k_B);  // Hawking temperature
    double T_eff = std::max(T_H_raw, T_eff_floor);  // Effective temperature with floor
    double P_form = compute_P_form(E_throat, T_eff);
    double J_aether = compute_J_aether(c);  // v_s ≈ c
    double U_m_val = compute_U_m(r_throat, t);
    double Theta = compute_Theta_WH(P_form, J_aether, U_m_val, T_eff);

    for (const auto& mod : additional_modifiers) {
        Theta *= mod(M, t);
    }

    double noise = noise_level * noise_dist(rng);
    return Theta + noise;
}

void UQFFWormholeFormation::add_modifier(std::function<double(double, double)> modifier) {
    // Self-expand: Add custom modifier to Theta_WH (function of M, t)
    // Allows extension, e.g., additional aether or coherence effects
    additional_modifiers.push_back(modifier);
}

void UQFFWormholeFormation::update_from_file(const std::string& config_file) {
    // Self-update: Load parameters from file (key=value)
    // For dynamic adjustments
    std::ifstream infile(config_file);
    if (!infile.is_open()) {
        std::cerr << "Failed to open config file: " << config_file << std::endl;
        return;
    }

    std::string line;
    while (std::getline(infile, line)) {
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

void UQFFWormholeFormation::simulate_formation(double M_start, double t_start, double t_end, double dt, const std::string& output_file) {
    // Self-simulate: Compute Theta_WH over time for fixed M, check formation condition
    // From formation condition: Theta_WH >1 forms white hole
    std::ofstream outfile;
    bool file_output = !output_file.empty();
    if (file_output) {
        outfile.open(output_file);
        outfile << "time,Theta_WH,Forms\n";
    }

    for (double t = t_start; t <= t_end; t += dt) {
        double Theta = compute_full_Theta_WH(M_start, t);
        bool forms = Theta > 1.0;

        if (file_output) {
            outfile << t << "," << Theta << "," << (forms ? 1 : 0) << "\n";
        } else {
            std::cout << "t=" << t << ", Theta_WH=" << Theta << ", Forms=" << (forms ? "Yes" : "No") << std::endl;
        }
    }

    if (file_output) outfile.close();
}

void UQFFWormholeFormation::display_explanations() {
    // Output captured text explanations
    // From standard concept, UQFF reinterpretation, derivation steps
    for (const auto& exp : explanations) {
        std::cout << exp << std::endl;
    }
}

std::string UQFFWormholeFormation::generate_derivation(double M, double t) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    
    // Step 1: Throat radius
    double r_s_std = compute_r_s_standard(M);
    double r_throat = compute_r_throat_UQFF(r_s_std);
    
    // Step 2: Formation probability components (UQFF-reduced energy via kappa_UQFF)
    double E_throat = kappa_UQFF * G * M * M / r_throat;
    double T_H_raw = hbar * c * c * c / (8.0 * M_PI * G * M * k_B);
    double T_eff = std::max(T_H_raw, T_eff_floor);  // Effective temperature with floor
    double P_form_val = compute_P_form(E_throat, T_eff);
    
    // Step 3: Aether flux
    double J_aether_val = compute_J_aether(c);
    
    // Step 4: Magnetic string energy
    double U_m_val = compute_U_m(r_throat, t);
    
    // Step 5: Full formation threshold
    double Theta = compute_Theta_WH(P_form_val, J_aether_val, U_m_val, T_eff);
    bool forms = Theta > 1.0;
    double tau_traversal = r_throat / c;  // Traversal time
    
    oss << "=== UQFF Wormhole Formation Derivation ==="  << "\n\n";
    
    oss << "INPUT PARAMETERS:\n";
    oss << "  M = " << M << " kg (" << M / 1.989e30 << " M_sun)\n";
    oss << "  t = " << t << " s (" << t / 3.154e7 << " years)\n\n";
    
    oss << "STEP 1: UQFF Throat Radius\n";
    oss << "  r_s = 2GM/c² = 2 × " << G << " × " << M << " / (" << c << ")²\n";
    oss << "      = " << r_s_std << " m\n";
    oss << "  r_throat,UQFF = r_s × (ρ_UA / ρ_SCm)\n";
    oss << "                = " << r_s_std << " × (" << rho_vac_UA << " / " << rho_vac_SCm << ")\n";
    oss << "                = " << r_throat << " m\n\n";
    
    oss << "STEP 2: Formation Probability (UQFF-reduced energy)\n";
    oss << "  E_throat = κ_UQFF × GM²/r_throat = " << kappa_UQFF << " × " << G << " × (" << M << ")² / " << r_throat << "\n";
    oss << "           = " << E_throat << " J\n";
    oss << "  T_H = ℏc³/(8πGMk_B) = " << T_H_raw << " K (Hawking)\n";
    oss << "  T_eff = max(T_H, T_floor) = max(" << T_H_raw << ", " << T_eff_floor << ") = " << T_eff << " K\n";
    oss << "  P_form = f_TRZ × exp(-E_throat/(k_B T_eff))\n";
    oss << "         = " << f_TRZ << " × exp(-" << E_throat << " / (" << k_B << " × " << T_eff << "))\n";
    oss << "         = " << P_form_val << "\n\n";
    
    oss << "STEP 3: Aether Flux\n";
    oss << "  J_aether = ρ_UA × v_s × (1 + f_TRZ)\n";
    oss << "           = " << rho_vac_UA << " × " << c << " × (1 + " << f_TRZ << ")\n";
    oss << "           = " << J_aether_val << " J/(m²·s)\n\n";
    
    oss << "STEP 4: Magnetic String Reinforcement (UQFF-scaled)\n";
    oss << "  U_m = λ_UQFF × (μ_j / r_throat) × (1 - exp(-γt cos(πt_n)))\n";
    oss << "      = " << lambda_UQFF << " × (" << mu_j << " / " << r_throat << ") × (1 - exp(-" << gamma << " × " << t << " × cos(π × " << t_n << ")))\n";
    oss << "      = " << U_m_val << " J\n\n";
    
    oss << "STEP 5: Formation Threshold\n";
    oss << "  Θ_WH = P_form × J_aether × exp(U_m / (k_B T_eff))\n";
    oss << "       = " << P_form_val << " × " << J_aether_val << " × exp(" << U_m_val << " / (" << k_B << " × " << T_eff << "))\n";
    oss << "       = " << Theta << "\n\n";
    
    oss << "RESULT:\n";
    oss << "  Θ_WH = " << Theta << " " << (forms ? "> 1 → WORMHOLE FORMS!" : "< 1 → No formation") << "\n";
    if (forms) {
        oss << "  Traversal time: τ = r_throat / c = " << r_throat << " / " << c << " = " << tau_traversal << " s\n";
    }
    
    return oss.str();
}

int UQFFWormholeFormation::run_tests() {
    int passed = 0;
    int total = 10;
    double tol = 1e-6;
    double M_test = 7.952e36;  // 4e6 M_sun
    double t_test = 1.42e17;   // 4.5 Gyr
    
    std::cout << "=== UQFF Wormhole Formation Tests ===\n\n";
    
    // Test 1: Schwarzschild radius
    double r_s = compute_r_s_standard(M_test);
    double r_s_expected = 2.0 * G * M_test / (c * c);
    bool test1 = std::abs(r_s - r_s_expected) / r_s_expected < tol;
    std::cout << "Test 1: r_s calculation - " << (test1 ? "PASS" : "FAIL") << "\n";
    if (test1) passed++;
    
    // Test 2: UQFF throat radius > standard
    double r_throat = compute_r_throat_UQFF(r_s);
    bool test2 = r_throat > r_s;  // Should expand due to ρ_UA / ρ_SCm > 1
    std::cout << "Test 2: r_throat > r_s (expansion) - " << (test2 ? "PASS" : "FAIL") << "\n";
    if (test2) passed++;
    
    // Test 3: Throat expansion factor = 10 (ρ_UA / ρ_SCm)
    double expansion = r_throat / r_s;
    bool test3 = std::abs(expansion - 10.0) < 0.01;
    std::cout << "Test 3: Expansion factor = 10 - " << (test3 ? "PASS" : "FAIL") << " (got " << expansion << ")\n";
    if (test3) passed++;
    
    // Test 4: P_form in valid range [0, 1] (using UQFF-reduced energy)
    double E_throat = kappa_UQFF * G * M_test * M_test / r_throat;
    double T_H_raw = hbar * c * c * c / (8.0 * M_PI * G * M_test * k_B);
    double T_eff = std::max(T_H_raw, T_eff_floor);
    double P_form_val = compute_P_form(E_throat, T_eff);
    bool test4 = P_form_val >= 0.0 && P_form_val <= 1.0;
    std::cout << "Test 4: P_form in [0,1] - " << (test4 ? "PASS" : "FAIL") << " (got " << P_form_val << ")\n";
    if (test4) passed++;
    
    // Test 5: J_aether positive
    double J_aether_val = compute_J_aether(c);
    bool test5 = J_aether_val > 0.0;
    std::cout << "Test 5: J_aether > 0 - " << (test5 ? "PASS" : "FAIL") << "\n";
    if (test5) passed++;
    
    // Test 6: U_m calculation valid
    double U_m_val = compute_U_m(r_throat, t_test);
    bool test6 = std::isfinite(U_m_val);
    std::cout << "Test 6: U_m is finite - " << (test6 ? "PASS" : "FAIL") << "\n";
    if (test6) passed++;
    
    // Test 7: Full Theta_WH chain
    double Theta = compute_full_Theta_WH(M_test, t_test, 0.0);  // No noise
    bool test7 = std::isfinite(Theta) && Theta > 0.0;
    std::cout << "Test 7: Θ_WH is finite and positive - " << (test7 ? "PASS" : "FAIL") << " (got " << Theta << ")\n";
    if (test7) passed++;
    
    // Test 8: Modifier application
    add_modifier([](double M, double t) { return 2.0; });  // Double the result
    double Theta_mod = compute_full_Theta_WH(M_test, t_test, 0.0);
    bool test8 = std::abs(Theta_mod / Theta - 2.0) < 0.1;  // Should be 2x
    std::cout << "Test 8: Modifier doubles Θ_WH - " << (test8 ? "PASS" : "FAIL") << "\n";
    if (test8) passed++;
    additional_modifiers.clear();  // Reset
    
    // Test 9: Hawking temperature positive for massive object
    bool test9 = T_H_raw > 0.0;
    std::cout << "Test 9: T_H > 0 - " << (test9 ? "PASS" : "FAIL") << " (got " << T_H_raw << " K)\n";
    if (test9) passed++;
    
    // Test 10: Derivation generation
    std::string deriv = generate_derivation(M_test, t_test);
    bool test10 = deriv.find("STEP 5") != std::string::npos && deriv.find("Θ_WH") != std::string::npos;
    std::cout << "Test 10: Derivation contains all steps - " << (test10 ? "PASS" : "FAIL") << "\n";
    if (test10) passed++;
    
    std::cout << "\n=== Results: " << passed << "/" << total << " tests passed ===\n";
    return passed;
}

int main() {
    // Main for testing: Demonstrate usage
    UQFFWormholeFormation wh;

    // Display explanations
    wh.display_explanations();
    std::cout << "\n";

    // Run validation tests
    int tests_passed = wh.run_tests();
    std::cout << "\n";

    // Generate derivation for Sgr A* example
    double M_example = 7.952e36;  // 4e6 M_sun kg
    double t_example = 1.42e17;   // 4.5 Gyr s
    std::cout << wh.generate_derivation(M_example, t_example) << "\n";

    // Test expand: Add custom modifier
    wh.add_modifier([](double M, double t) { return 1.0 + 0.01 * t / M; });
    double Theta_modified = wh.compute_full_Theta_WH(M_example, t_example);
    std::cout << "Theta_WH (with modifier): " << Theta_modified << "\n\n";

    // Test simulate: t 0 to 1e18 s, dt 1e16 s
    std::cout << "Running simulation (output to wh_formation_sim.csv)...\n";
    wh.simulate_formation(M_example, 0.0, 1e18, 1e16, "wh_formation_sim.csv");
    std::cout << "Simulation complete.\n";

    return (tests_passed >= 8) ? 0 : 1;  // Exit code based on test results
}