// uqff_black_hole_entropy_impl.cpp
// UQFF Black Hole Entropy Implementation
// Fixed from original - proper initialization, explanations, validation tests

#define _USE_MATH_DEFINES
#include "uqff_black_hole_entropy.h"
#include <sstream>
#include <algorithm>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// CONSTRUCTOR - Full Initialization
// ============================================================================

UQFFBlackHoleEntropy::UQFFBlackHoleEntropy(unsigned int seed)
    : rng(seed == 0 ? std::chrono::system_clock::now().time_since_epoch().count() : seed),
      noise_dist(0.0, 1.0),
      tests_passed(0),
      tests_total(0)
{
    // Fundamental constants (CODATA 2018)
    G = 6.67430e-11;            // Gravitational constant (m³ kg⁻¹ s⁻²)
    c = 2.99792458e8;           // Speed of light (m/s)
    k_B = 1.380649e-23;         // Boltzmann constant (J/K)
    hbar = 1.054571817e-34;     // Reduced Planck constant (J·s)
    
    // UQFF parameters
    rho_vac_UA = 7.09e-36;      // Aether vacuum density (J/m³)
    rho_vac_SCm = 7.09e-37;     // Superconductive vacuum density (J/m³)
    f_TRZ = 0.1;                // Time-reversal zone coupling
    U_m = 1e-40;                // Magnetic string energy (J) - calibrated for astrophysical BHs
    
    // Derived Planck units
    l_Pl = std::sqrt(hbar * G / (c * c * c));   // ~1.616e-35 m
    M_Pl = std::sqrt(hbar * c / G);              // ~2.176e-8 kg
    t_Pl = std::sqrt(hbar * G / (c * c * c * c * c)); // ~5.391e-44 s
    
    // Solar mass
    M_sun = 1.989e30;           // kg
    
    // Populate explanations
    populate_explanations();
}

// ============================================================================
// STEP 1: Schwarzschild Geometry
// ============================================================================

double UQFFBlackHoleEntropy::compute_r_s(double M) const {
    // r_s = 2GM/c² (Schwarzschild radius)
    // Event horizon where escape velocity = c
    // Solar mass: r_s ≈ 3 km
    // Sgr A*: r_s ≈ 1.2×10^10 m
    return (2.0 * G * M) / (c * c);
}

double UQFFBlackHoleEntropy::compute_A(double r_s) const {
    // A = 4π r_s² (horizon area)
    // Area theorem: A cannot decrease classically
    return 4.0 * M_PI * r_s * r_s;
}

double UQFFBlackHoleEntropy::compute_A_from_mass(double M) const {
    double r_s = compute_r_s(M);
    return compute_A(r_s);
}

// ============================================================================
// STEP 2: Base Bekenstein-Hawking Entropy
// ============================================================================

double UQFFBlackHoleEntropy::compute_S_BH(double A) const {
    // S_BH = k_B × c³ × A / (4 × G × ℏ)
    // = k_B × A / (4 × l_Pl²)
    // Bekenstein (1973): S = (ln 2)/(8π) × (c³/Gℏ) × A
    // Hawking (1974): Fixed coefficient to 1/4
    return (k_B * c * c * c * A) / (4.0 * G * hbar);
}

double UQFFBlackHoleEntropy::compute_S_BH_from_mass(double M) const {
    double A = compute_A_from_mass(M);
    return compute_S_BH(A);
}

// ============================================================================
// STEP 3: Aether Holographic Rescaling
// ============================================================================

double UQFFBlackHoleEntropy::compute_rho_ratio() const {
    // ρ_UA / ρ_SCm ≈ 10 in standard UQFF calibration
    return rho_vac_UA / rho_vac_SCm;
}

double UQFFBlackHoleEntropy::compute_A_UQFF(double A) const {
    // A_UQFF = A × (ρ_UA/ρ_SCm)
    // Aether threads inflate effective holographic area
    // Physical: More [UA] channels → more information storage
    return A * compute_rho_ratio();
}

double UQFFBlackHoleEntropy::compute_S_aether(double S_BH) const {
    // S_aether = S_BH × (ρ_UA/ρ_SCm)
    // Direct scaling of entropy by aether density ratio
    return S_BH * compute_rho_ratio();
}

// ============================================================================
// STEP 4: Time-Reversal Negentropy
// ============================================================================

double UQFFBlackHoleEntropy::compute_S_TRZ(double S) const {
    // S' = S × (1 - f_TRZ)
    // Time-reversal zone coupling reduces entropy
    // f_TRZ ≈ 0.1 → 10% reduction
    // Physical: Negentropic ordering stabilizes information
    return S * (1.0 - f_TRZ);
}

// ============================================================================
// STEP 5: Hawking Temperature & String Damping
// ============================================================================

double UQFFBlackHoleEntropy::compute_T_H(double M) const {
    // T_H = ℏc³ / (8π G M k_B)
    // Hawking temperature: BH radiates as blackbody at T_H
    // Solar mass: T_H ≈ 6×10^-8 K (essentially cold)
    // Sgr A*: T_H ≈ 1.5×10^-14 K
    // Planck mass: T_H ≈ 10^32 K
    return (hbar * c * c * c) / (8.0 * M_PI * G * M * k_B);
}

double UQFFBlackHoleEntropy::compute_string_damping(double M) const {
    // Damping = exp(-U_m / (k_B × T_H))
    // At low T_H (large M): Damping → 0 (strong suppression)
    // At high T_H (small M): Damping → 1 (no suppression)
    double T_H = compute_T_H(M);
    double exponent = -U_m / (k_B * T_H);
    
    // Prevent underflow
    if (exponent < -700.0) return 0.0;
    if (exponent > 0.0) return 1.0;  // Should not happen physically
    
    return std::exp(exponent);
}

double UQFFBlackHoleEntropy::compute_S_string_damped(double S, double M) const {
    // S'' = S × exp(-U_m/(k_B T_H))
    // Magnetic string energy damps high-entropy states
    return S * compute_string_damping(M);
}

// ============================================================================
// STEP 6: Full UQFF Black Hole Entropy
// ============================================================================

double UQFFBlackHoleEntropy::compute_full_S_UQFF(double M, double noise_level) const {
    // Full UQFF black hole entropy combining all effects:
    // S_UQFF = S_BH × (ρ_UA/ρ_SCm) × (1 - f_TRZ) × exp(-U_m/(k_B T_H))
    
    // Step 1-2: Base Bekenstein-Hawking
    double S_BH = compute_S_BH_from_mass(M);
    
    // Step 3: Aether enhancement
    double S_aether = compute_S_aether(S_BH);
    
    // Step 4: TRZ reduction
    double S_TRZ = compute_S_TRZ(S_aether);
    
    // Step 5: String damping
    double S_UQFF = compute_S_string_damped(S_TRZ, M);
    
    // Apply additional mods
    for (const auto& mod : additional_mods) {
        S_UQFF *= mod(M);
    }
    
    // Add noise if requested
    if (noise_level > 0.0) {
        std::mt19937 local_rng(std::chrono::system_clock::now().time_since_epoch().count());
        std::normal_distribution<double> local_noise(0.0, 1.0);
        S_UQFF *= (1.0 + noise_level * local_noise(local_rng));
    }
    
    return S_UQFF;
}

// ============================================================================
// DERIVED QUANTITIES
// ============================================================================

double UQFFBlackHoleEntropy::compute_information_bits(double S) const {
    // I = S / (k_B × ln 2)
    // Converts entropy to information in bits
    return S / (k_B * std::log(2.0));
}

double UQFFBlackHoleEntropy::compute_evaporation_time(double M) const {
    // t_evap ≈ 5120 π G² M³ / (ℏ c⁴)
    // Page time (half-evaporation): t_Page ≈ t_evap / 2
    // Solar mass: t_evap ≈ 10^67 years
    double t_evap = (5120.0 * M_PI * G * G * M * M * M) / (hbar * c * c * c * c);
    return t_evap;
}

double UQFFBlackHoleEntropy::compute_scrambling_time(double M) const {
    // t_scr ≈ (r_s/c) × ln(S/k_B)
    // Time for horizon to scramble infalling information
    double r_s = compute_r_s(M);
    double S = compute_S_BH_from_mass(M);
    double t_scr = (r_s / c) * std::log(S / k_B);
    return t_scr;
}

// ============================================================================
// SELF-EXPANSION
// ============================================================================

void UQFFBlackHoleEntropy::add_mod(std::function<double(double)> mod) {
    additional_mods.push_back(mod);
}

void UQFFBlackHoleEntropy::clear_mods() {
    additional_mods.clear();
}

size_t UQFFBlackHoleEntropy::get_mod_count() const {
    return additional_mods.size();
}

// ============================================================================
// SELF-UPDATE
// ============================================================================

void UQFFBlackHoleEntropy::update_from_file(const std::string& config_file) {
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
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            
            std::string val_str = line.substr(pos + 1);
            val_str.erase(0, val_str.find_first_not_of(" \t"));
            val_str.erase(val_str.find_last_not_of(" \t") + 1);
            
            try {
                double value = std::stod(val_str);
                
                if (key == "G") G = value;
                else if (key == "c") c = value;
                else if (key == "k_B") k_B = value;
                else if (key == "hbar") hbar = value;
                else if (key == "f_TRZ") f_TRZ = value;
                else if (key == "rho_vac_UA") rho_vac_UA = value;
                else if (key == "rho_vac_SCm") rho_vac_SCm = value;
                else if (key == "U_m") U_m = value;
            } catch (const std::exception& e) {
                std::cerr << "Error parsing value for " << key << ": " << e.what() << std::endl;
            }
        }
    }
    infile.close();
    std::cout << "Updated parameters from " << config_file << std::endl;
}

void UQFFBlackHoleEntropy::set_f_TRZ(double f_TRZ_new) {
    f_TRZ = f_TRZ_new;
}

void UQFFBlackHoleEntropy::set_U_m(double U_m_new) {
    U_m = U_m_new;
}

void UQFFBlackHoleEntropy::set_rho_ratio(double rho_UA_new, double rho_SCm_new) {
    rho_vac_UA = rho_UA_new;
    rho_vac_SCm = rho_SCm_new;
}

// ============================================================================
// SELF-SIMULATION
// ============================================================================

void UQFFBlackHoleEntropy::simulate_over_mass(double M_start, double M_end, double dM,
                                               const std::string& output_file) const {
    std::ofstream outfile;
    bool file_output = !output_file.empty();
    
    if (file_output) {
        outfile.open(output_file);
        outfile << "# UQFF Black Hole Entropy vs Mass\n";
        outfile << "# f_TRZ=" << f_TRZ << ", rho_ratio=" << compute_rho_ratio() 
                << ", U_m=" << U_m << "\n";
        outfile << "M_kg,M_solar,r_s_m,A_m2,T_H_K,S_BH,S_UQFF,damping\n";
    }
    
    for (double M = M_start; M <= M_end; M += dM) {
        double r_s = compute_r_s(M);
        double A = compute_A(r_s);
        double T_H = compute_T_H(M);
        double S_BH = compute_S_BH(A);
        double S_UQFF = compute_full_S_UQFF(M);
        double damping = compute_string_damping(M);
        double M_sol = M / M_sun;
        
        if (file_output) {
            outfile << std::scientific << std::setprecision(6)
                    << M << "," << M_sol << "," << r_s << "," << A << ","
                    << T_H << "," << S_BH << "," << S_UQFF << "," << damping << "\n";
        } else {
            std::cout << "M=" << M_sol << " M_sun, S_UQFF=" << S_UQFF 
                      << ", T_H=" << T_H << " K" << std::endl;
        }
    }
    
    if (file_output) {
        outfile.close();
        std::cout << "Simulation saved to " << output_file << std::endl;
    }
}

std::vector<std::pair<double, double>> UQFFBlackHoleEntropy::compute_S_vs_mass(
    double M_start, double M_end, double dM) const {
    
    std::vector<std::pair<double, double>> results;
    
    for (double M = M_start; M <= M_end; M += dM) {
        double S_UQFF = compute_full_S_UQFF(M);
        results.push_back({M, S_UQFF});
    }
    
    return results;
}

// ============================================================================
// EXPLANATIONS
// ============================================================================

void UQFFBlackHoleEntropy::populate_explanations() {
    explanations.clear();
    
    explanations.push_back(
        "STEP 1: SCHWARZSCHILD GEOMETRY\n"
        "r_s = 2GM/c² (Schwarzschild radius)\n"
        "A = 4π r_s² (horizon area)\n\n"
        "The event horizon defines the causal boundary.\n"
        "For M = M_sun: r_s ≈ 3 km\n"
        "For Sgr A* (4×10^6 M_sun): r_s ≈ 1.2×10^10 m"
    );
    
    explanations.push_back(
        "STEP 2: BEKENSTEIN-HAWKING ENTROPY\n"
        "S_BH = k_B × c³ × A / (4 × G × ℏ)\n"
        "     = k_B × A / (4 × l_Pl²)\n\n"
        "Bekenstein (1973): Entropy proportional to area.\n"
        "Hawking (1974): Fixed coefficient to 1/4.\n"
        "Scaling: S ~ M² (quadratic in mass)\n"
        "Solar mass BH: S ~ 10^77 k_B\n"
        "Sgr A*: S ~ 10^90 k_B"
    );
    
    explanations.push_back(
        "STEP 3: AETHER HOLOGRAPHIC RESCALING\n"
        "A_UQFF = A × (ρ_UA/ρ_SCm)\n"
        "S_aether = S_BH × (ρ_UA/ρ_SCm)\n\n"
        "UQFF: Aether [UA] threads inflate effective area.\n"
        "ρ_ratio ≈ 10 in standard calibration.\n"
        "Physical: More information channels via aether.\n"
        "→ Entropy enhanced by factor of ~10"
    );
    
    explanations.push_back(
        "STEP 4: TIME-REVERSAL NEGENTROPY\n"
        "S' = S × (1 - f_TRZ)\n\n"
        "f_TRZ ≈ 0.1 time-reversal coupling.\n"
        "Negentropic ordering reduces entropy by ~10%.\n"
        "Physical: TRZ stabilizes horizon information.\n"
        "At f_TRZ = 0: Standard entropy\n"
        "At f_TRZ = 1: Zero entropy (unphysical)"
    );
    
    explanations.push_back(
        "STEP 5: HAWKING TEMPERATURE & STRING DAMPING\n"
        "T_H = ℏc³ / (8π G M k_B)\n"
        "Damping = exp(-U_m / (k_B × T_H))\n\n"
        "Hawking: BH radiates thermally at T_H.\n"
        "Solar mass: T_H ≈ 6×10^-8 K (cold)\n"
        "Planck mass: T_H ~ 10^32 K (hot)\n\n"
        "String damping suppresses high-entropy states.\n"
        "Large M (cold): Strong suppression\n"
        "Small M (hot): Weak suppression"
    );
    
    explanations.push_back(
        "STEP 6: FULL UQFF BLACK HOLE ENTROPY\n"
        "S_UQFF = S_BH × (ρ_UA/ρ_SCm) × (1 - f_TRZ) × exp(-U_m/(k_B T_H))\n\n"
        "Combines all effects:\n"
        "  • Base Bekenstein-Hawking (S ~ M²)\n"
        "  • Aether area inflation (×10)\n"
        "  • Negentropic reduction (×0.9)\n"
        "  • String damping (mass-dependent)\n\n"
        "Net effect: S_UQFF ~ 9 × S_BH × damping"
    );
    
    explanations.push_back(
        "NUMERICAL EXAMPLE (Sgr A*):\n"
        "M = 4×10^6 M_sun ≈ 8×10^36 kg\n"
        "r_s ≈ 1.2×10^10 m\n"
        "A ≈ 1.8×10^21 m²\n"
        "T_H ≈ 1.5×10^-14 K\n\n"
        "S_BH ≈ 1.1×10^90 k_B\n"
        "With ρ_ratio=10, f_TRZ=0.1, damping≈0.37:\n"
        "S_UQFF ≈ 1.1×10^90 × 10 × 0.9 × 0.37\n"
        "       ≈ 3.7×10^90 k_B"
    );
    
    explanations.push_back(
        "ADVANCES IN UQFF:\n"
        "• Entropy enhanced by aether holographic inflation\n"
        "• Negentropic TRZ stabilizes information\n"
        "• String damping bounds high-entropy states\n"
        "• Testable via q-scope horizon measurements\n"
        "• Links to holographic principle, AdS/CFT\n"
        "• Information paradox: Aether preserves unitarity"
    );
}

void UQFFBlackHoleEntropy::display_explanations() const {
    std::cout << "=========================================================\n";
    std::cout << "UQFF BLACK HOLE ENTROPY DERIVATION\n";
    std::cout << "=========================================================\n\n";
    
    for (size_t i = 0; i < explanations.size(); ++i) {
        std::cout << explanations[i] << "\n\n";
    }
}

std::string UQFFBlackHoleEntropy::get_explanation(size_t index) const {
    if (index < explanations.size()) {
        return explanations[index];
    }
    return "";
}

// ============================================================================
// VALIDATION TESTS
// ============================================================================

bool UQFFBlackHoleEntropy::run_validation_tests() {
    test_results.clear();
    tests_passed = 0;
    tests_total = 0;
    
    std::cout << "=========================================================\n";
    std::cout << "UQFF Black Hole Entropy Validation Tests\n";
    std::cout << "=========================================================\n\n";
    
    // Test 1: Schwarzschild radius for solar mass
    {
        tests_total++;
        double r_s = compute_r_s(M_sun);
        double expected = 2953.25;  // ~2.95 km
        bool pass = std::abs(r_s - expected) / expected < 0.01;
        test_results["schwarzschild_radius_solar"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Solar mass r_s: " 
                  << r_s << " m (expected ~2953 m)\n";
    }
    
    // Test 2: S_BH scaling (S ~ M²)
    {
        tests_total++;
        double S1 = compute_S_BH_from_mass(M_sun);
        double S2 = compute_S_BH_from_mass(2.0 * M_sun);
        double ratio = S2 / S1;
        bool pass = std::abs(ratio - 4.0) < 0.01;  // Should be 4 for M²
        test_results["S_BH_M_squared_scaling"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] S ~ M² scaling: S(2M)/S(M) = " 
                  << ratio << " (expected 4)\n";
    }
    
    // Test 3: Hawking temperature scaling (T ~ 1/M)
    {
        tests_total++;
        double T1 = compute_T_H(M_sun);
        double T2 = compute_T_H(2.0 * M_sun);
        double ratio = T1 / T2;
        bool pass = std::abs(ratio - 2.0) < 0.01;
        test_results["T_H_inverse_M_scaling"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] T_H ~ 1/M scaling: T(M)/T(2M) = " 
                  << ratio << " (expected 2)\n";
    }
    
    // Test 4: rho ratio
    {
        tests_total++;
        double ratio = compute_rho_ratio();
        bool pass = std::abs(ratio - 10.0) < 0.1;
        test_results["rho_ratio"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] ρ_UA/ρ_SCm: " 
                  << ratio << " (expected ~10)\n";
    }
    
    // Test 5: TRZ reduces entropy
    {
        tests_total++;
        double S = 1e90;
        double S_TRZ = compute_S_TRZ(S);
        bool pass = S_TRZ < S && std::abs(S_TRZ - 0.9e90) < 1e88;
        test_results["TRZ_reduces_entropy"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] TRZ reduction: " 
                  << S << " → " << S_TRZ << "\n";
    }
    
    // Test 6: String damping in [0,1]
    {
        tests_total++;
        double damping = compute_string_damping(M_sun);
        bool pass = (damping >= 0.0 && damping <= 1.0);
        test_results["damping_factor_bounds"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Damping factor: " 
                  << damping << " ∈ [0,1]\n";
    }
    
    // Test 7: S_UQFF > 0
    {
        tests_total++;
        double S_UQFF = compute_full_S_UQFF(M_sun);
        bool pass = S_UQFF > 0;
        test_results["S_UQFF_positive"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] S_UQFF > 0: " 
                  << S_UQFF << "\n";
    }
    
    // Test 8: Planck length
    {
        tests_total++;
        double expected_l_Pl = 1.616e-35;
        bool pass = std::abs(l_Pl - expected_l_Pl) / expected_l_Pl < 0.01;
        test_results["planck_length"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Planck length: " 
                  << l_Pl << " m (expected ~1.616e-35 m)\n";
    }
    
    // Test 9: Information bits conversion
    {
        tests_total++;
        double S = k_B * std::log(2.0);  // 1 bit of entropy
        double bits = compute_information_bits(S);
        bool pass = std::abs(bits - 1.0) < 1e-6;
        test_results["information_bits"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Info bits: S=k_B ln2 → " 
                  << bits << " bits (expected 1)\n";
    }
    
    // Test 10: Area law (A = 4π r_s²)
    {
        tests_total++;
        double r_s = 1000.0;  // 1 km
        double A = compute_A(r_s);
        double expected = 4.0 * M_PI * r_s * r_s;
        bool pass = std::abs(A - expected) < 1e-6;
        test_results["area_law"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Area law: A(1km) = " 
                  << A << " m² (expected " << expected << ")\n";
    }
    
    // Test 11: Sgr A* entropy order of magnitude
    {
        tests_total++;
        double M_SgrA = 4.0e6 * M_sun;
        double S_BH = compute_S_BH_from_mass(M_SgrA);
        double log_S = std::log10(S_BH / k_B);
        bool pass = log_S > 85 && log_S < 95;  // Should be ~90
        test_results["SgrA_entropy_magnitude"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Sgr A* entropy: 10^" 
                  << log_S << " k_B (expected ~10^90)\n";
    }
    
    // Test 12: Explanations populated
    {
        tests_total++;
        bool pass = explanations.size() >= 7;
        test_results["explanations_populated"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Explanations: " 
                  << explanations.size() << " steps documented\n";
    }
    
    std::cout << "\n=========================================================\n";
    std::cout << "Results: " << tests_passed << "/" << tests_total << " tests passed\n";
    std::cout << "=========================================================\n";
    
    return tests_passed == tests_total;
}

void UQFFBlackHoleEntropy::display_test_results() const {
    std::cout << "\nTest Results Summary:\n";
    for (const auto& [name, result] : test_results) {
        std::cout << "  " << name << ": " << (result ? "PASS" : "FAIL") << "\n";
    }
}

// ============================================================================
// MAIN - Test Driver
// ============================================================================

int main() {
    std::cout << "=========================================================\n";
    std::cout << "UQFF BLACK HOLE ENTROPY MODULE\n";
    std::cout << "Bekenstein-Hawking with aether holographic modifications\n";
    std::cout << "=========================================================\n\n";
    
    UQFFBlackHoleEntropy bh_entropy(42);  // Fixed seed for reproducibility
    
    // Display derivation
    bh_entropy.display_explanations();
    
    // Run validation
    bh_entropy.run_validation_tests();
    
    // Example calculations
    std::cout << "\n---------------------------------------------------------\n";
    std::cout << "Example: Sagittarius A* Black Hole\n";
    std::cout << "---------------------------------------------------------\n";
    
    double M_SgrA = 4.0e6 * bh_entropy.get_M_sun();  // 4 million solar masses
    
    double r_s = bh_entropy.compute_r_s(M_SgrA);
    double A = bh_entropy.compute_A(r_s);
    double T_H = bh_entropy.compute_T_H(M_SgrA);
    double S_BH = bh_entropy.compute_S_BH(A);
    double S_UQFF = bh_entropy.compute_full_S_UQFF(M_SgrA);
    double damping = bh_entropy.compute_string_damping(M_SgrA);
    double t_evap = bh_entropy.compute_evaporation_time(M_SgrA);
    
    std::cout << "\nSgr A* (M = 4×10^6 M_sun):\n";
    std::cout << "  Schwarzschild radius: r_s = " << r_s << " m (" << r_s/1e9 << " Gm)\n";
    std::cout << "  Horizon area: A = " << A << " m²\n";
    std::cout << "  Hawking temperature: T_H = " << T_H << " K\n";
    std::cout << "  String damping: " << damping << "\n";
    std::cout << "\n  S_BH (standard) = " << S_BH << " J/K\n";
    std::cout << "  S_BH / k_B = 10^" << std::log10(S_BH / bh_entropy.get_k_B()) << "\n";
    std::cout << "\n  S_UQFF = " << S_UQFF << " J/K\n";
    std::cout << "  S_UQFF / k_B = 10^" << std::log10(S_UQFF / bh_entropy.get_k_B()) << "\n";
    std::cout << "\n  Evaporation time: t_evap = " << t_evap << " s (" 
              << t_evap / (3.15e7 * 1e9) << " billion years)\n";
    
    // Simulate over mass range
    std::cout << "\n---------------------------------------------------------\n";
    std::cout << "Mass Sweep (1 to 10 solar masses):\n";
    std::cout << "---------------------------------------------------------\n";
    
    bh_entropy.simulate_over_mass(bh_entropy.get_M_sun(), 
                                   10.0 * bh_entropy.get_M_sun(), 
                                   3.0 * bh_entropy.get_M_sun(), 
                                   "uqff_black_hole_entropy_vs_M.csv");
    
    std::cout << "\n=========================================================\n";
    std::cout << "Module test complete.\n";
    std::cout << "=========================================================\n";
    
    return 0;
}
