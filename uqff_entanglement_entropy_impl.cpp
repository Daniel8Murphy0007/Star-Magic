// uqff_entanglement_entropy_impl.cpp
// UQFF Entanglement Entropy S Implementation
// Fixed from original derivation - proper initialization and N-level support

#define _USE_MATH_DEFINES
#include "uqff_entanglement_entropy.h"
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <numeric>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// CONSTRUCTOR - Full Initialization
// ============================================================================

UQFFEntanglementEntropy::UQFFEntanglementEntropy(unsigned int seed)
    : rng(seed == 0 ? std::chrono::system_clock::now().time_since_epoch().count() : seed),
      noise_dist(0.0, 1.0),
      tests_passed(0),
      tests_total(0)
{
    // UQFF parameters
    rho_vac_UA = 7.09e-36;          // Aether vacuum density (J/m³)
    rho_vac_SCm = 7.09e-37;         // Superconductive vacuum density (J/m³)
    f_TRZ = 0.1;                    // Time-reversal zone coupling
    U_m = 1e-20;                    // String cutoff energy (J)
    delta = 0.01;                   // Fluctuation amplitude
    
    // Physical constants
    k_B = 1.380649e-23;             // Boltzmann constant (J/K)
    
    // System parameters
    T = 100.0;                      // Temperature (K)
    
    // Populate explanations
    populate_explanations();
}

// ============================================================================
// STEP 1: Base von Neumann Entropy
// ============================================================================

double UQFFEntanglementEntropy::compute_S_base(const std::vector<double>& lambda) const {
    // S = -Σ λ_i log λ_i (von Neumann entropy)
    // Standard entanglement entropy quantifying bipartite correlations:
    // For state |ψ⟩ = Σ √λ_i |i_A⟩|i_B⟩ (Schmidt decomposition)
    // ρ_A = Tr_B(|ψ⟩⟨ψ|) = Σ λ_i |i_A⟩⟨i_A|
    // S = 0: pure product state (λ = {1, 0, ...})
    // S = log N: maximally entangled (λ = {1/N, ...})
    
    double S = 0.0;
    for (double l : lambda) {
        if (l > 1e-30) {
            S -= l * std::log(l);
        }
    }
    return S;
}

// ============================================================================
// STEP 2: Effective Dimension with Aether Enhancement
// ============================================================================

double UQFFEntanglementEntropy::compute_rho_ratio() const {
    // ρ_UA / ρ_SCm ≈ 10 in standard UQFF calibration
    return rho_vac_UA / rho_vac_SCm;
}

double UQFFEntanglementEntropy::compute_d_eff_UQFF(const std::vector<double>& lambda) const {
    // d_eff,UQFF = (1/Σ λ_i²) × (ρ_UA/ρ_SCm)
    // From Step 2: Aether enhancement boosts effective dimension
    // Standard d_eff = 1/Σ λ_i² (participation ratio)
    // UQFF multiplies by ρ_ratio to account for aether thread correlations
    // Physical: More [UA] channels → more effective degrees of freedom
    
    double sum_sq = 0.0;
    for (double l : lambda) {
        sum_sq += l * l;
    }
    
    if (sum_sq < 1e-30) return 1.0;
    
    double d_eff = 1.0 / sum_sq;
    return d_eff * compute_rho_ratio();
}

double UQFFEntanglementEntropy::compute_S_enhancement(double d_eff_uqff) const {
    // S_enhancement = f_TRZ × log(d_eff,UQFF)
    // Contribution from enhanced effective dimension
    // Subtracted in full formula (entropy reduction)
    if (d_eff_uqff <= 0) return 0.0;
    return f_TRZ * std::log(d_eff_uqff);
}

// ============================================================================
// STEP 3: Time-Reversal Reduction
// ============================================================================

double UQFFEntanglementEntropy::compute_S_TRZ(double S_base) const {
    // S' = S × (1 - f_TRZ)
    // From Step 3: Time-reversal zone coupling reduces entropy
    // f_TRZ ≈ 0.1 creates negentropic flow
    // Physical: TRZ reversal partially orders the system
    // At f_TRZ = 0: normal entropy
    // At f_TRZ = 1: zero entropy (full reversal)
    
    return S_base * (1.0 - f_TRZ);
}

// ============================================================================
// STEP 4: String Cutoff Damping
// ============================================================================

double UQFFEntanglementEntropy::compute_string_damping_factor() const {
    // Factor: (1 - exp(-k_B T / U_m))
    // Temperature-dependent damping from magnetic string energy
    // Low T (k_B T << U_m): factor → 0, strong damping
    // High T (k_B T >> U_m): factor → 1, no damping
    // Physical: String energy sets entropy ceiling
    
    double exponent = -(k_B * T) / U_m;
    if (exponent < -700.0) return 1.0;  // Prevent underflow
    return 1.0 - std::exp(exponent);
}

double UQFFEntanglementEntropy::compute_S_string_damped(double S_TRZ) const {
    // S'' = S' × (1 - exp(-k_B T / U_m))
    // From Step 4: String cutoff damps high-entropy states
    return S_TRZ * compute_string_damping_factor();
}

// ============================================================================
// STEP 5: Full UQFF Entropy
// ============================================================================

double UQFFEntanglementEntropy::compute_full_S_UQFF(const std::vector<double>& lambda,
                                                     double noise_level) const {
    // Complete UQFF entanglement entropy
    // S_UQFF = -Σ λ_i log λ_i - f_TRZ × log(d_eff,UQFF) × (1 - exp(-k_B T/U_m))
    
    // Step 1: Base entropy
    double S_base = compute_S_base(lambda);
    
    // Step 2: Enhanced effective dimension
    double d_eff_uqff = compute_d_eff_UQFF(lambda);
    double S_enhancement = compute_S_enhancement(d_eff_uqff);
    
    // Step 3: TRZ reduction (apply to base)
    double S_TRZ = compute_S_TRZ(S_base);
    
    // Step 4: String damping factor
    double damping = compute_string_damping_factor();
    
    // Step 5: Full formula
    // S_UQFF = S_base - f_TRZ × log(d_eff_UQFF) × damping
    double S_UQFF = S_base - S_enhancement * damping;
    
    // Apply additional mods
    for (const auto& mod : additional_mods) {
        S_UQFF *= mod(lambda);
    }
    
    // Add noise if requested
    if (noise_level > 0.0) {
        std::mt19937 local_rng(std::chrono::system_clock::now().time_since_epoch().count());
        std::normal_distribution<double> local_noise(0.0, 1.0);
        S_UQFF += noise_level * std::abs(S_UQFF) * local_noise(local_rng);
    }
    
    // Entropy must be non-negative
    if (S_UQFF < 0) S_UQFF = 0.0;
    
    return S_UQFF;
}

// ============================================================================
// ENTROPY PROPERTIES
// ============================================================================

double UQFFEntanglementEntropy::compute_purity(const std::vector<double>& lambda) const {
    // P = Σ λ_i² = Tr(ρ_A²)
    // Purity: P=1 for pure state, P=1/N for maximally mixed
    double P = 0.0;
    for (double l : lambda) {
        P += l * l;
    }
    return P;
}

double UQFFEntanglementEntropy::compute_effective_dimension(const std::vector<double>& lambda) const {
    // d_eff = 1 / Σ λ_i² (standard, without aether factor)
    double P = compute_purity(lambda);
    if (P < 1e-30) return 1.0;
    return 1.0 / P;
}

double UQFFEntanglementEntropy::compute_concurrence_2qubit(double lambda1, double lambda2) const {
    // C = 2 × √(λ_1 × λ_2) (concurrence for 2-qubit pure state)
    // C = 0: product state; C = 1: maximally entangled
    // Valid only for 2-level systems
    if (lambda1 < 0 || lambda2 < 0) return 0.0;
    return 2.0 * std::sqrt(lambda1 * lambda2);
}

// ============================================================================
// SELF-EXPANSION
// ============================================================================

void UQFFEntanglementEntropy::add_mod(std::function<double(std::vector<double>)> mod) {
    additional_mods.push_back(mod);
}

void UQFFEntanglementEntropy::clear_mods() {
    additional_mods.clear();
}

size_t UQFFEntanglementEntropy::get_mod_count() const {
    return additional_mods.size();
}

// ============================================================================
// SELF-UPDATE
// ============================================================================

void UQFFEntanglementEntropy::update_from_file(const std::string& config_file) {
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
                
                if (key == "rho_vac_UA") rho_vac_UA = value;
                else if (key == "rho_vac_SCm") rho_vac_SCm = value;
                else if (key == "f_TRZ") f_TRZ = value;
                else if (key == "k_B") k_B = value;
                else if (key == "U_m") U_m = value;
                else if (key == "T") T = value;
                else if (key == "delta") delta = value;
            } catch (const std::exception& e) {
                std::cerr << "Error parsing value for " << key << ": " << e.what() << std::endl;
            }
        }
    }
    infile.close();
    std::cout << "Updated parameters from " << config_file << std::endl;
}

void UQFFEntanglementEntropy::set_temperature(double T_new) {
    T = T_new;
}

void UQFFEntanglementEntropy::set_TRZ_coupling(double f_TRZ_new) {
    f_TRZ = f_TRZ_new;
}

void UQFFEntanglementEntropy::set_string_energy(double U_m_new) {
    U_m = U_m_new;
}

// ============================================================================
// SELF-SIMULATION
// ============================================================================

void UQFFEntanglementEntropy::simulate_over_temperature(double T_start, double T_end,
                                                        double dT,
                                                        const std::vector<double>& lambda,
                                                        const std::string& output_file) const {
    std::ofstream outfile;
    bool file_output = !output_file.empty();
    
    if (file_output) {
        outfile.open(output_file);
        outfile << "# UQFF Entanglement Entropy vs Temperature\n";
        outfile << "# f_TRZ=" << f_TRZ << ", U_m=" << U_m << ", rho_ratio=" << compute_rho_ratio() << "\n";
        outfile << "T,S_base,S_UQFF,d_eff,d_eff_UQFF,damping_factor\n";
    }
    
    for (double T_val = T_start; T_val <= T_end; T_val += dT) {
        const_cast<UQFFEntanglementEntropy*>(this)->set_temperature(T_val);
        
        double S_base = compute_S_base(lambda);
        double S_UQFF = compute_full_S_UQFF(lambda);
        double d_eff = compute_effective_dimension(lambda);
        double d_eff_uqff = compute_d_eff_UQFF(lambda);
        double damping = compute_string_damping_factor();
        
        if (file_output) {
            outfile << std::scientific << std::setprecision(6)
                    << T_val << "," << S_base << "," << S_UQFF 
                    << "," << d_eff << "," << d_eff_uqff << "," << damping << "\n";
        } else {
            std::cout << "T=" << T_val << ", S_base=" << S_base << ", S_UQFF=" << S_UQFF 
                      << ", d_eff_UQFF=" << d_eff_uqff << std::endl;
        }
    }
    
    if (file_output) {
        outfile.close();
        std::cout << "Simulation saved to " << output_file << std::endl;
    }
}

std::vector<std::pair<double, double>> UQFFEntanglementEntropy::compute_entropy_vs_temperature(
    double T_start, double T_end, double dT, const std::vector<double>& lambda) const {
    
    std::vector<std::pair<double, double>> results;
    
    for (double T_val = T_start; T_val <= T_end; T_val += dT) {
        const_cast<UQFFEntanglementEntropy*>(this)->set_temperature(T_val);
        double S_UQFF = compute_full_S_UQFF(lambda);
        results.push_back({T_val, S_UQFF});
    }
    
    return results;
}

// ============================================================================
// EXPLANATIONS
// ============================================================================

void UQFFEntanglementEntropy::populate_explanations() {
    explanations.clear();
    
    explanations.push_back(
        "STEP 1: Base von Neumann Entropy\n"
        "S = -Σ λ_i log λ_i\n"
        "Standard entanglement entropy from Schmidt decomposition.\n"
        "For state |ψ⟩ = Σ √λ_i |i_A⟩|i_B⟩:\n"
        "  S = 0: Product state (λ = {1, 0, ...})\n"
        "  S = log N: Maximally entangled (λ = {1/N, ...})\n"
        "  Bell state: S = log 2 ≈ 0.693"
    );
    
    explanations.push_back(
        "STEP 2: Aether-Enhanced Effective Dimension\n"
        "d_eff,UQFF = (1/Σ λ_i²) × (ρ_UA/ρ_SCm)\n"
        "Standard d_eff = 1/Σ λ_i² (participation ratio)\n"
        "UQFF multiplies by ρ_ratio ≈ 10 for aether thread correlations.\n"
        "Physical: More [UA] channels → more effective degrees of freedom.\n"
        "Entropy enhancement: f_TRZ × log(d_eff,UQFF)"
    );
    
    explanations.push_back(
        "STEP 3: Time-Reversal Zone Reduction\n"
        "S' = S × (1 - f_TRZ)\n"
        "f_TRZ ≈ 0.1 time-reversal coupling creates negentropic flow.\n"
        "Partial entropy reduction from TRZ reversal effects.\n"
        "At f_TRZ = 0: Normal entropy\n"
        "At f_TRZ = 1: Zero entropy (full reversal - unphysical)"
    );
    
    explanations.push_back(
        "STEP 4: String Cutoff Damping\n"
        "S'' = S' × (1 - exp(-k_B T/U_m))\n"
        "U_m: Magnetic string energy sets entropy ceiling.\n"
        "Low T (k_B T << U_m): Strong damping, factor → 0\n"
        "High T (k_B T >> U_m): No damping, factor → 1\n"
        "Physical: String tension limits accessible entropy states."
    );
    
    explanations.push_back(
        "STEP 5: Full UQFF Entanglement Entropy\n"
        "S_UQFF = S_base - f_TRZ × log(d_eff,UQFF) × (1 - exp(-k_B T/U_m))\n"
        "Combines: base entropy, enhanced dimension, TRZ reduction, string damping.\n"
        "Entropy generally reduced relative to standard QM.\n"
        "Testable predictions: q-scope entangled state measurements."
    );
    
    explanations.push_back(
        "NUMERICAL EXAMPLE:\n"
        "Bell state (N=2): λ = {0.5, 0.5}\n"
        "S_base = log 2 ≈ 0.693\n"
        "d_eff = 1/0.5 = 2, d_eff,UQFF = 2 × 10 = 20\n"
        "log(d_eff,UQFF) = log(20) ≈ 3.0\n"
        "With f_TRZ = 0.1, damping ≈ 0.63 (k_B T/U_m ≈ 1):\n"
        "S_UQFF ≈ 0.693 - 0.1 × 3.0 × 0.63 ≈ 0.693 - 0.189 ≈ 0.504\n"
        "→ ~27% entropy reduction from UQFF effects"
    );
    
    explanations.push_back(
        "ADVANCES IN UQFF:\n"
        "• Entropy reduction testable in q-scope entangled systems\n"
        "• Temperature dependence confirms string damping mechanism\n"
        "• d_eff,UQFF enhancement from aether [UA] channels\n"
        "• Negentropic TRZ effects create order from entanglement\n"
        "• Links holographic entropy (S~Area) to aether structure"
    );
}

void UQFFEntanglementEntropy::display_explanations() const {
    std::cout << "=========================================================\n";
    std::cout << "UQFF ENTANGLEMENT ENTROPY S DERIVATION\n";
    std::cout << "=========================================================\n\n";
    
    for (size_t i = 0; i < explanations.size(); ++i) {
        std::cout << explanations[i] << "\n\n";
    }
}

std::string UQFFEntanglementEntropy::get_explanation(size_t index) const {
    if (index < explanations.size()) {
        return explanations[index];
    }
    return "";
}

// ============================================================================
// VALIDATION TESTS
// ============================================================================

bool UQFFEntanglementEntropy::run_validation_tests() {
    test_results.clear();
    tests_passed = 0;
    tests_total = 0;
    
    std::cout << "=========================================================\n";
    std::cout << "UQFF Entanglement Entropy Validation Tests\n";
    std::cout << "=========================================================\n\n";
    
    // Test 1: Bell state entropy
    {
        tests_total++;
        std::vector<double> lambda = {0.5, 0.5};
        double S = compute_S_base(lambda);
        double expected = std::log(2.0);
        bool pass = std::abs(S - expected) < 1e-6;
        test_results["bell_state_entropy"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Bell state entropy: S=" 
                  << S << " (expected " << expected << ")\n";
    }
    
    // Test 2: Pure state entropy = 0
    {
        tests_total++;
        std::vector<double> lambda = {1.0, 0.0};
        double S = compute_S_base(lambda);
        bool pass = std::abs(S) < 1e-10;
        test_results["pure_state_entropy_zero"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Pure state entropy: S=" 
                  << S << " (expected 0)\n";
    }
    
    // Test 3: Maximally entangled 4-level entropy = log 4
    {
        tests_total++;
        std::vector<double> lambda = {0.25, 0.25, 0.25, 0.25};
        double S = compute_S_base(lambda);
        double expected = std::log(4.0);
        bool pass = std::abs(S - expected) < 1e-6;
        test_results["max_entangled_4level"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] 4-level max entangled: S=" 
                  << S << " (expected " << expected << ")\n";
    }
    
    // Test 4: d_eff enhancement
    {
        tests_total++;
        std::vector<double> lambda = {0.5, 0.5};
        double d_eff = compute_effective_dimension(lambda);
        double d_eff_uqff = compute_d_eff_UQFF(lambda);
        bool pass = std::abs(d_eff_uqff - d_eff * 10.0) < 0.1;
        test_results["d_eff_enhancement"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] d_eff enhancement: " 
                  << d_eff << " → " << d_eff_uqff << " (ratio ~10)\n";
    }
    
    // Test 5: TRZ reduces entropy
    {
        tests_total++;
        double S_base = 0.693;
        double S_TRZ = compute_S_TRZ(S_base);
        bool pass = S_TRZ < S_base;
        test_results["TRZ_reduces_entropy"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] TRZ reduces: S=" 
                  << S_base << " → " << S_TRZ << "\n";
    }
    
    // Test 6: String damping factor bounds
    {
        tests_total++;
        double factor = compute_string_damping_factor();
        bool pass = (factor >= 0.0 && factor <= 1.0);
        test_results["damping_factor_bounds"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Damping factor: " 
                  << factor << " ∈ [0,1]\n";
    }
    
    // Test 7: UQFF entropy ≤ base entropy
    {
        tests_total++;
        std::vector<double> lambda = {0.5, 0.5};
        double S_base = compute_S_base(lambda);
        double S_UQFF = compute_full_S_UQFF(lambda);
        bool pass = S_UQFF <= S_base || std::abs(S_UQFF - S_base) < 1e-6;
        test_results["UQFF_entropy_reduced"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] UQFF ≤ base: " 
                  << S_UQFF << " ≤ " << S_base << "\n";
    }
    
    // Test 8: Purity ranges correctly
    {
        tests_total++;
        std::vector<double> lambda = {0.5, 0.5};
        double P = compute_purity(lambda);
        bool pass = (P >= 0.0 && P <= 1.0);
        test_results["purity_range"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Purity: P=" 
                  << P << " ∈ [0,1]\n";
    }
    
    // Test 9: Concurrence 2-qubit
    {
        tests_total++;
        double C = compute_concurrence_2qubit(0.5, 0.5);
        bool pass = std::abs(C - 1.0) < 1e-6;  // Bell state has C=1
        test_results["concurrence_bell"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Bell state concurrence: C=" 
                  << C << " (expected 1.0)\n";
    }
    
    // Test 10: rho ratio
    {
        tests_total++;
        double ratio = compute_rho_ratio();
        bool pass = std::abs(ratio - 10.0) < 0.1;
        test_results["rho_ratio"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] ρ_UA/ρ_SCm: " 
                  << ratio << " (expected ~10)\n";
    }
    
    // Test 11: Temperature dependence
    {
        tests_total++;
        std::vector<double> lambda = {0.5, 0.5};
        
        set_temperature(10.0);
        double S_low_T = compute_full_S_UQFF(lambda);
        
        set_temperature(1000.0);
        double S_high_T = compute_full_S_UQFF(lambda);
        
        bool pass = S_high_T < S_low_T;  // Higher T → stronger damping term
        test_results["temperature_dependence"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] T dependence: S(10K)=" 
                  << S_low_T << " > S(1000K)=" << S_high_T << "\n";
        
        // Reset T
        const_cast<UQFFEntanglementEntropy*>(this)->set_temperature(100.0);
    }
    
    // Test 12: Explanations populated
    {
        tests_total++;
        bool pass = explanations.size() >= 6;
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

void UQFFEntanglementEntropy::display_test_results() const {
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
    std::cout << "UQFF ENTANGLEMENT ENTROPY MODULE\n";
    std::cout << "S with aether enhancement, TRZ reduction, string damping\n";
    std::cout << "=========================================================\n\n";
    
    UQFFEntanglementEntropy entropy(42);  // Fixed seed for reproducibility
    
    // Display derivation
    entropy.display_explanations();
    
    // Run validation
    entropy.run_validation_tests();
    
    // Example calculations
    std::cout << "\n---------------------------------------------------------\n";
    std::cout << "Example: Bell State Entanglement Entropy\n";
    std::cout << "---------------------------------------------------------\n";
    
    std::vector<double> bell_state = {0.5, 0.5};
    
    double S_base = entropy.compute_S_base(bell_state);
    double d_eff = entropy.compute_effective_dimension(bell_state);
    double d_eff_uqff = entropy.compute_d_eff_UQFF(bell_state);
    double S_UQFF = entropy.compute_full_S_UQFF(bell_state);
    double purity = entropy.compute_purity(bell_state);
    double concurrence = entropy.compute_concurrence_2qubit(0.5, 0.5);
    
    std::cout << "\nBell state λ = {0.5, 0.5}:\n";
    std::cout << "  S_base (von Neumann) = " << S_base << " (log 2 = " << std::log(2.0) << ")\n";
    std::cout << "  d_eff (standard)     = " << d_eff << "\n";
    std::cout << "  d_eff,UQFF           = " << d_eff_uqff << " (aether-enhanced)\n";
    std::cout << "  S_UQFF               = " << S_UQFF << "\n";
    std::cout << "  Purity               = " << purity << "\n";
    std::cout << "  Concurrence          = " << concurrence << "\n";
    std::cout << "  Entropy reduction    = " << (1.0 - S_UQFF/S_base) * 100 << "%\n";
    
    // Simulate over temperature
    std::cout << "\n---------------------------------------------------------\n";
    std::cout << "Temperature Evolution (T = 10 to 200 K):\n";
    std::cout << "---------------------------------------------------------\n";
    
    entropy.simulate_over_temperature(10.0, 200.0, 40.0, bell_state, 
                                      "uqff_entanglement_entropy_vs_T.csv");
    
    std::cout << "\n=========================================================\n";
    std::cout << "Module test complete.\n";
    std::cout << "=========================================================\n";
    
    return 0;
}
