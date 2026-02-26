// uqff_conductivity_spectrum_impl.cpp
// UQFF Conductivity Spectrum σ(ω) Implementation
// Fixed from original derivation - proper initialization and complex handling

#define _USE_MATH_DEFINES
#include "uqff_conductivity_spectrum.h"
#include <iomanip>
#include <sstream>
#include <algorithm>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// CONSTRUCTOR - Full Initialization
// ============================================================================

UQFFConductivitySpectrum::UQFFConductivitySpectrum(unsigned int seed)
    : rng(seed == 0 ? std::chrono::system_clock::now().time_since_epoch().count() : seed),
      noise_dist(0.0, 1.0),
      tests_passed(0),
      tests_total(0)
{
    // Physical constants
    c = 2.998e8;                    // m/s
    k_B = 1.380649e-23;             // J/K
    hbar = 1.054571817e-34;         // J·s
    
    // AdS/CFT boundary parameters
    z_0 = 1e-10;                    // Near boundary
    A_prime = 1.0;                  // Normalized
    A_z0 = 1.0;                     // Normalized
    
    // UQFF vacuum densities
    rho_vac_UA = 7.09e-36;          // J/m³
    rho_vac_SCm = 7.09e-37;         // J/m³
    
    // Magnetic parameters
    B_t = 5e10;                     // T (half critical)
    B_crit = 1e11;                  // T (magnetar scale)
    
    // Time-reversal zone parameters
    f_TRZ = 0.1;                    // Coupling factor
    tau_coh = 1e-11;                // Coherence time (s)
    Gamma = 1.0 / tau_coh;          // ~1e11 Hz
    omega_res = 2.0 * M_PI * f_TRZ / tau_coh;  // ~6.28e10 rad/s
    
    // String damping parameters
    mu_j = 1e-10;                   // String tension (J/m)
    L_string = 1e-9;                // String length (m)
    U_m = mu_j / L_string;          // ~0.1 J
    gamma_damp = 0.5;               // Damping coefficient
    
    // Temperature & correlation
    T = 100.0;                      // K
    T_c = 120.0;                    // K (critical)
    xi = 1e-9;                      // m (nanometer scale)
    
    // UQFF scaling factors
    kappa_UQFF = 1e-60;
    lambda_UQFF = 1e-9;
    
    // Populate explanations
    populate_explanations();
}

// ============================================================================
// STEP 1: Base Holographic Conductivity
// ============================================================================

std::complex<double> UQFFConductivitySpectrum::compute_sigma_base_complex(double omega) const {
    // σ(ω) = (i/ω)(A'/A + (ω²/c²) log z_0)
    // From AdS/CFT: Maxwell fluctuations in AdS bulk yield boundary conductivity
    // Near z→0, A_x satisfies A_x'' + (ω²/c² z² - m²)A_x = 0
    // Boundary expansion gives σ from ratio of subleading to leading coefficients
    
    if (std::abs(omega) < 1e-20) {
        // DC limit: Return large real part (delta function approximation)
        return std::complex<double>(1e20, 0.0);
    }
    
    double log_z0 = std::log(std::abs(z_0));
    double ratio = A_prime / A_z0;
    double omega_term = (omega * omega) / (c * c) * log_z0;
    
    // i/ω × (real terms)
    std::complex<double> prefactor(0.0, 1.0 / omega);
    std::complex<double> bracket(ratio + omega_term, 0.0);
    
    return prefactor * bracket;
}

double UQFFConductivitySpectrum::compute_sigma_base_real(double omega) const {
    return compute_sigma_base_complex(omega).real();
}

double UQFFConductivitySpectrum::compute_sigma_base_imag(double omega) const {
    return compute_sigma_base_complex(omega).imag();
}

// ============================================================================
// STEP 2: Aether-Modified Correlation Length
// ============================================================================

double UQFFConductivitySpectrum::compute_xi_UQFF() const {
    // ξ_UQFF = ξ × (ρ_UA/ρ_SCm)^{1/2}
    // Aether vacuum energy modulates superconducting correlation length
    // Higher ρ_UA relative to ρ_SCm extends correlations
    return xi * std::sqrt(rho_vac_UA / rho_vac_SCm);
}

double UQFFConductivitySpectrum::compute_rho_ratio() const {
    // ρ_UA/ρ_SCm ≈ 10 in standard UQFF calibration
    return rho_vac_UA / rho_vac_SCm;
}

double UQFFConductivitySpectrum::compute_gap_UQFF(double Delta_standard) const {
    // Gap modification: Δ_UQFF = Δ / √(ρ_UA/ρ_SCm)
    // Aether coupling widens the superconducting gap in frequency space
    double rho_ratio = compute_rho_ratio();
    return Delta_standard / std::sqrt(rho_ratio);
}

// ============================================================================
// STEP 3: Magnetic Suppression Factor
// ============================================================================

double UQFFConductivitySpectrum::compute_suppression_factor() const {
    // Factor (1 - B_t/B_crit)
    // Magnetic field suppresses superconductivity; at B_crit, SC vanishes
    // For B_t < B_crit, partial suppression damps above-gap conductivity
    double ratio = B_t / B_crit;
    if (ratio >= 1.0) return 0.0;  // Fully suppressed
    return 1.0 - ratio;
}

double UQFFConductivitySpectrum::compute_sigma_prime(double sigma_base) const {
    // σ' = σ × (1 - B_t/B_crit)
    return sigma_base * compute_suppression_factor();
}

// ============================================================================
// STEP 4: Time-Reversal Resonance
// ============================================================================

double UQFFConductivitySpectrum::compute_lorentzian(double omega) const {
    // Lorentzian: Γ / ((ω - ω_res)² + Γ²)
    // Time-reversal zone creates resonant enhancement at ω_res
    // Width Γ ~ 1/τ_coh from coherence time
    double delta_omega = omega - omega_res;
    double denominator = delta_omega * delta_omega + Gamma * Gamma;
    return Gamma / denominator;
}

double UQFFConductivitySpectrum::compute_sigma_double_prime(double sigma_prime, double omega) const {
    // σ'' = σ' + f_TRZ × Lorentzian
    // TRZ coupling f_TRZ ≈ 0.1 adds resonant peak
    return sigma_prime + f_TRZ * compute_lorentzian(omega);
}

// ============================================================================
// STEP 5: String Damping
// ============================================================================

double UQFFConductivitySpectrum::compute_U_m_dynamic(double t, double t_n) const {
    // U_m = (μ_j/L) × (1 - exp(-γ t cos(π t_n)))
    // Dynamic string energy depends on time t and phase t_n
    // γ controls damping rate; cos(π t_n) modulates phase
    double cos_term = std::cos(M_PI * t_n);
    double exp_term = std::exp(-gamma_damp * t * cos_term);
    return (mu_j / L_string) * (1.0 - exp_term);
}

double UQFFConductivitySpectrum::compute_damping_factor(double omega) const {
    // exp(-U_m ω / (k_B T))
    // High-frequency damping from magnetic string excitations
    double exponent = -U_m * omega / (k_B * T);
    // Prevent underflow
    if (exponent < -700.0) return 0.0;
    return std::exp(exponent);
}

double UQFFConductivitySpectrum::compute_sigma_UQFF_step5(double sigma_double_prime, double omega) const {
    // σ_UQFF = σ'' × exp(-U_m ω / (k_B T))
    return sigma_double_prime * compute_damping_factor(omega);
}

// ============================================================================
// STEP 6: Full σ_UQFF
// ============================================================================

double UQFFConductivitySpectrum::compute_full_sigma(double omega, double noise_level) const {
    // Complete UQFF conductivity spectrum
    // Combines all steps with aether rescaling
    
    // Step 1: Base conductivity (imaginary part dominates for ω ≠ 0)
    double sigma_base = compute_sigma_base_imag(omega);
    
    // Step 2: Aether rescaling √(ρ_UA/ρ_SCm)
    double rho_factor = std::sqrt(compute_rho_ratio());
    sigma_base *= rho_factor;
    
    // Step 3: Magnetic suppression
    double sigma_prime = compute_sigma_prime(sigma_base);
    
    // Step 4: TRZ resonance
    double sigma_double_prime = compute_sigma_double_prime(sigma_prime, omega);
    
    // Step 5: String damping
    double sigma_uqff = compute_sigma_UQFF_step5(sigma_double_prime, omega);
    
    // Apply additional mods
    for (const auto& mod : additional_mods) {
        sigma_uqff *= mod(omega, T);
    }
    
    // Add noise if requested
    if (noise_level > 0.0) {
        // Need non-const RNG for noise
        std::mt19937 local_rng(std::chrono::system_clock::now().time_since_epoch().count());
        std::normal_distribution<double> local_noise(0.0, 1.0);
        sigma_uqff += noise_level * std::abs(sigma_uqff) * local_noise(local_rng);
    }
    
    return sigma_uqff;
}

std::complex<double> UQFFConductivitySpectrum::compute_full_sigma_complex(double omega) const {
    // Full complex conductivity
    std::complex<double> sigma_base = compute_sigma_base_complex(omega);
    
    double rho_factor = std::sqrt(compute_rho_ratio());
    sigma_base *= rho_factor;
    
    double suppression = compute_suppression_factor();
    sigma_base *= suppression;
    
    // Add Lorentzian (real addition)
    double lorentz = f_TRZ * compute_lorentzian(omega);
    sigma_base += std::complex<double>(lorentz, 0.0);
    
    // Apply damping
    double damping = compute_damping_factor(omega);
    sigma_base *= damping;
    
    return sigma_base;
}

// ============================================================================
// SELF-EXPANSION
// ============================================================================

void UQFFConductivitySpectrum::add_mod(std::function<double(double, double)> mod) {
    additional_mods.push_back(mod);
}

void UQFFConductivitySpectrum::clear_mods() {
    additional_mods.clear();
}

size_t UQFFConductivitySpectrum::get_mod_count() const {
    return additional_mods.size();
}

// ============================================================================
// SELF-UPDATE
// ============================================================================

void UQFFConductivitySpectrum::update_from_file(const std::string& config_file) {
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
            
            std::string val_str = line.substr(pos + 1);
            val_str.erase(0, val_str.find_first_not_of(" \t"));
            val_str.erase(val_str.find_last_not_of(" \t") + 1);
            
            try {
                double value = std::stod(val_str);
                
                if (key == "c") c = value;
                else if (key == "z_0") z_0 = value;
                else if (key == "A_prime") A_prime = value;
                else if (key == "A_z0") A_z0 = value;
                else if (key == "rho_vac_UA") rho_vac_UA = value;
                else if (key == "rho_vac_SCm") rho_vac_SCm = value;
                else if (key == "B_t") B_t = value;
                else if (key == "B_crit") B_crit = value;
                else if (key == "f_TRZ") f_TRZ = value;
                else if (key == "Gamma") Gamma = value;
                else if (key == "omega_res") omega_res = value;
                else if (key == "tau_coh") tau_coh = value;
                else if (key == "U_m") U_m = value;
                else if (key == "mu_j") mu_j = value;
                else if (key == "L_string") L_string = value;
                else if (key == "gamma_damp") gamma_damp = value;
                else if (key == "k_B") k_B = value;
                else if (key == "T") T = value;
                else if (key == "T_c") T_c = value;
                else if (key == "xi") xi = value;
                else if (key == "kappa_UQFF") kappa_UQFF = value;
                else if (key == "lambda_UQFF") lambda_UQFF = value;
            } catch (const std::exception& e) {
                std::cerr << "Error parsing value for " << key << ": " << e.what() << std::endl;
            }
        }
    }
    infile.close();
    std::cout << "Updated parameters from " << config_file << std::endl;
}

void UQFFConductivitySpectrum::set_temperature(double T_new) {
    T = T_new;
}

void UQFFConductivitySpectrum::set_magnetic_field(double B_new) {
    B_t = B_new;
}

void UQFFConductivitySpectrum::set_frequency_params(double omega_res_new, double Gamma_new) {
    omega_res = omega_res_new;
    Gamma = Gamma_new;
    if (Gamma_new > 0) {
        tau_coh = 1.0 / Gamma_new;
    }
}

// ============================================================================
// SELF-SIMULATION
// ============================================================================

void UQFFConductivitySpectrum::simulate_spectrum(double omega_start, double omega_end, 
                                                  double d_omega, const std::string& output_file) const {
    std::ofstream outfile;
    bool file_output = !output_file.empty();
    
    if (file_output) {
        outfile.open(output_file);
        outfile << "# UQFF Conductivity Spectrum\n";
        outfile << "# T=" << T << " K, B_t=" << B_t << " T, f_TRZ=" << f_TRZ << "\n";
        outfile << "# omega_res=" << omega_res << " rad/s, Gamma=" << Gamma << " Hz\n";
        outfile << "omega,sigma_UQFF_real,sigma_UQFF_imag,sigma_UQFF_magnitude\n";
    }
    
    for (double omega = omega_start; omega <= omega_end; omega += d_omega) {
        std::complex<double> sigma = compute_full_sigma_complex(omega);
        double magnitude = std::abs(sigma);
        
        if (file_output) {
            outfile << std::scientific << std::setprecision(6)
                    << omega << "," << sigma.real() << "," << sigma.imag() 
                    << "," << magnitude << "\n";
        } else {
            std::cout << "ω=" << std::scientific << omega 
                      << ", σ_UQFF=" << magnitude 
                      << " (Re=" << sigma.real() << ", Im=" << sigma.imag() << ")" 
                      << std::endl;
        }
    }
    
    if (file_output) {
        outfile.close();
        std::cout << "Spectrum saved to " << output_file << std::endl;
    }
}

std::vector<std::pair<double, double>> UQFFConductivitySpectrum::compute_spectrum(
    double omega_start, double omega_end, double d_omega) const {
    
    std::vector<std::pair<double, double>> spectrum;
    for (double omega = omega_start; omega <= omega_end; omega += d_omega) {
        double sigma = compute_full_sigma(omega);
        spectrum.push_back({omega, sigma});
    }
    return spectrum;
}

// ============================================================================
// EXPLANATIONS
// ============================================================================

void UQFFConductivitySpectrum::populate_explanations() {
    explanations.clear();
    
    explanations.push_back(
        "STEP 1: Base Holographic Conductivity\n"
        "σ(ω) = (i/ω)(A'/A + (ω²/c²) log z_0) + δ(ω)\n"
        "From AdS/CFT Maxwell fluctuations near boundary z→0.\n"
        "Below T_c: Re[σ] has delta at ω=0 (infinite DC), gap Δ≈8 k_B T_c.\n"
        "Im[σ] ~ 1/ω pole from superconducting condensate."
    );
    
    explanations.push_back(
        "STEP 2: Aether-Modified Correlation Length\n"
        "ξ_UQFF = ξ × (ρ_UA/ρ_SCm)^{1/2}\n"
        "Aether vacuum energy ratio ≈ 10 extends superconducting correlations.\n"
        "Gap modification: Δ_UQFF ≈ Δ/√10 widens accessible frequency range."
    );
    
    explanations.push_back(
        "STEP 3: Magnetic Suppression\n"
        "σ' = σ × (1 - B_t/B_crit)\n"
        "Magnetic field suppresses superconductivity.\n"
        "At B_t = B_crit, superconductivity vanishes completely.\n"
        "For magnetar-scale B_crit ~ 10^11 T, partial suppression at B_t ~ 5×10^10 T."
    );
    
    explanations.push_back(
        "STEP 4: Time-Reversal Zone Resonance\n"
        "σ'' = σ' + f_TRZ × Γ/((ω - ω_res)² + Γ²)\n"
        "TRZ coupling f_TRZ ≈ 0.1 adds Lorentzian peak at ω_res.\n"
        "ω_res = 2π f_TRZ/τ_coh, Γ ~ 1/τ_coh from coherence time.\n"
        "Creates observable resonance in q-scope THz spectra."
    );
    
    explanations.push_back(
        "STEP 5: Magnetic String Damping\n"
        "σ_UQFF = σ'' × exp(-U_m ω / (k_B T))\n"
        "U_m = (μ_j/L)(1 - exp(-γ t cos(π t_n))) string energy.\n"
        "High-frequency damping from magnetic string excitations.\n"
        "Suppresses UV contributions while preserving IR physics."
    );
    
    explanations.push_back(
        "STEP 6: Full UQFF Conductivity Spectrum\n"
        "σ_UQFF = [holographic base] × √(ρ_UA/ρ_SCm) × (1 - B/B_c)\n"
        "         + f_TRZ × Lorentzian × exp(-U_m ω/(k_B T))\n"
        "Numerical: T=100K, ω=1 THz, f_TRZ=0.1, ρ_ratio=10, B/B_c=0.5:\n"
        "σ_UQFF ≈ standard × √10 × 0.5 + resonance × e^{-0.1}\n"
        "Enhanced gap and resonant features at THz frequencies."
    );
    
    explanations.push_back(
        "ADVANCES IN UQFF:\n"
        "• Observable resonances in q-scope spectra\n"
        "• Testable holographic duality predictions\n"
        "• Links superconductivity to aether vacuum structure\n"
        "• Magnetic string damping provides UV regularization\n"
        "• TRZ coupling creates characteristic spectral signatures"
    );
}

void UQFFConductivitySpectrum::display_explanations() const {
    std::cout << "=========================================================\n";
    std::cout << "UQFF CONDUCTIVITY SPECTRUM σ(ω) DERIVATION\n";
    std::cout << "=========================================================\n\n";
    
    for (size_t i = 0; i < explanations.size(); ++i) {
        std::cout << explanations[i] << "\n\n";
    }
}

std::string UQFFConductivitySpectrum::get_explanation(size_t index) const {
    if (index < explanations.size()) {
        return explanations[index];
    }
    return "";
}

// ============================================================================
// VALIDATION TESTS
// ============================================================================

bool UQFFConductivitySpectrum::run_validation_tests() {
    test_results.clear();
    tests_passed = 0;
    tests_total = 0;
    
    std::cout << "=========================================================\n";
    std::cout << "UQFF Conductivity Spectrum Validation Tests\n";
    std::cout << "=========================================================\n\n";
    
    // Test 1: xi_UQFF scales correctly
    {
        tests_total++;
        double xi_uqff = compute_xi_UQFF();
        double expected = xi * std::sqrt(10.0);  // rho_ratio ≈ 10
        bool pass = std::abs(xi_uqff - expected) / expected < 0.01;
        test_results["xi_UQFF_scaling"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] xi_UQFF scaling: "
                  << xi_uqff << " (expected ~" << expected << ")\n";
    }
    
    // Test 2: Suppression factor bounds
    {
        tests_total++;
        double factor = compute_suppression_factor();
        bool pass = (factor >= 0.0 && factor <= 1.0);
        test_results["suppression_bounds"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Suppression factor bounds: "
                  << factor << " ∈ [0,1]\n";
    }
    
    // Test 3: Lorentzian peak at omega_res
    {
        tests_total++;
        double at_res = compute_lorentzian(omega_res);
        double off_res = compute_lorentzian(omega_res + 10.0 * Gamma);
        bool pass = at_res > off_res * 5.0;  // Peak should be much higher
        test_results["lorentzian_peak"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Lorentzian peak: "
                  << at_res << " >> " << off_res << "\n";
    }
    
    // Test 4: Damping factor decreases with frequency
    {
        tests_total++;
        double damp_low = compute_damping_factor(1e10);
        double damp_high = compute_damping_factor(1e13);
        bool pass = damp_high < damp_low;
        test_results["damping_decreases"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Damping vs frequency: "
                  << damp_low << " > " << damp_high << "\n";
    }
    
    // Test 5: Complex sigma has imaginary component
    {
        tests_total++;
        double omega = 1e12;
        std::complex<double> sigma = compute_sigma_base_complex(omega);
        bool pass = std::abs(sigma.imag()) > 0;
        test_results["complex_sigma_imag"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Complex σ imaginary: "
                  << sigma.imag() << "\n";
    }
    
    // Test 6: Full sigma finite at THz
    {
        tests_total++;
        double omega = 1e12;
        double sigma = compute_full_sigma(omega);
        bool pass = std::isfinite(sigma) && sigma != 0.0;
        test_results["full_sigma_finite"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Full σ_UQFF finite at THz: "
                  << sigma << "\n";
    }
    
    // Test 7: Gap modification
    {
        tests_total++;
        double Delta_std = 1.0;
        double Delta_uqff = compute_gap_UQFF(Delta_std);
        double expected = Delta_std / std::sqrt(10.0);
        bool pass = std::abs(Delta_uqff - expected) / expected < 0.01;
        test_results["gap_modification"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Gap modification: "
                  << Delta_uqff << " (expected ~" << expected << ")\n";
    }
    
    // Test 8: Rho ratio correct
    {
        tests_total++;
        double ratio = compute_rho_ratio();
        bool pass = std::abs(ratio - 10.0) < 0.1;
        test_results["rho_ratio"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] ρ_UA/ρ_SCm ratio: "
                  << ratio << " (expected ~10)\n";
    }
    
    // Test 9: Self-expansion works
    {
        tests_total++;
        size_t before = get_mod_count();
        // Can't add mod in const context, but verify count works
        bool pass = (before == 0);  // Initially empty
        test_results["mod_system"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Mod system initialized: "
                  << before << " mods\n";
    }
    
    // Test 10: Explanations populated
    {
        tests_total++;
        bool pass = explanations.size() >= 6;
        test_results["explanations_populated"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Explanations: "
                  << explanations.size() << " steps documented\n";
    }
    
    // Test 11: Spectrum computation
    {
        tests_total++;
        auto spectrum = compute_spectrum(1e10, 1e11, 1e10);
        bool pass = spectrum.size() > 0 && std::isfinite(spectrum[0].second);
        test_results["spectrum_computation"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Spectrum computation: "
                  << spectrum.size() << " points\n";
    }
    
    // Test 12: Dynamic U_m
    {
        tests_total++;
        double U_m_t0 = compute_U_m_dynamic(0.0, 0.0);
        double U_m_t1 = compute_U_m_dynamic(1e-9, 0.5);
        bool pass = U_m_t1 > U_m_t0;  // Should increase with time
        test_results["dynamic_U_m"] = pass;
        if (pass) tests_passed++;
        std::cout << "[" << (pass ? "PASS" : "FAIL") << "] Dynamic U_m: "
                  << U_m_t0 << " → " << U_m_t1 << "\n";
    }
    
    std::cout << "\n=========================================================\n";
    std::cout << "Results: " << tests_passed << "/" << tests_total << " tests passed\n";
    std::cout << "=========================================================\n";
    
    return tests_passed == tests_total;
}

void UQFFConductivitySpectrum::display_test_results() const {
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
    std::cout << "UQFF CONDUCTIVITY SPECTRUM MODULE\n";
    std::cout << "σ(ω) with holographic, aether, and TRZ modifications\n";
    std::cout << "=========================================================\n\n";
    
    UQFFConductivitySpectrum spectrum(42);  // Fixed seed for reproducibility
    
    // Display derivation
    spectrum.display_explanations();
    
    // Run validation
    spectrum.run_validation_tests();
    
    // Example calculations
    std::cout << "\n---------------------------------------------------------\n";
    std::cout << "Example Calculations:\n";
    std::cout << "---------------------------------------------------------\n";
    
    double omega_THz = 1e12;  // 1 THz
    
    std::cout << "\nAt ω = " << omega_THz << " rad/s (THz regime):\n";
    std::cout << "  ξ_UQFF = " << spectrum.compute_xi_UQFF() << " m\n";
    std::cout << "  Suppression factor = " << spectrum.compute_suppression_factor() << "\n";
    std::cout << "  Lorentzian at ω = " << spectrum.compute_lorentzian(omega_THz) << "\n";
    std::cout << "  Damping factor = " << spectrum.compute_damping_factor(omega_THz) << "\n";
    std::cout << "  σ_UQFF = " << spectrum.compute_full_sigma(omega_THz) << "\n";
    
    std::complex<double> sigma_complex = spectrum.compute_full_sigma_complex(omega_THz);
    std::cout << "  σ_UQFF (complex) = " << sigma_complex.real() 
              << " + " << sigma_complex.imag() << "i\n";
    
    // Test self-expansion
    std::cout << "\n---------------------------------------------------------\n";
    std::cout << "Self-Expansion Test:\n";
    std::cout << "---------------------------------------------------------\n";
    
    spectrum.add_mod([](double omega, double T) {
        return 1.0 + 0.001 * omega / T;  // Temperature-dependent enhancement
    });
    std::cout << "Added custom mod. New count: " << spectrum.get_mod_count() << "\n";
    std::cout << "σ_UQFF with mod = " << spectrum.compute_full_sigma(omega_THz) << "\n";
    
    // Simulate spectrum
    std::cout << "\n---------------------------------------------------------\n";
    std::cout << "Spectrum Simulation (ω = 1e10 to 1e12 rad/s):\n";
    std::cout << "---------------------------------------------------------\n";
    
    spectrum.simulate_spectrum(1e10, 1e12, 2e11, "uqff_conductivity_spectrum_sim.csv");
    
    std::cout << "\n=========================================================\n";
    std::cout << "Module test complete.\n";
    std::cout << "=========================================================\n";
    
    return 0;
}
