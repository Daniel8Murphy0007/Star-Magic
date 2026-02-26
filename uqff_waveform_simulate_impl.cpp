// uqff_waveform_simulate_impl.cpp
// Implementation of UQFF Waveform Simulation with Chirp Evolution
// Fixes: Member initialization, correct suppression factors, chirp frequency

#define _USE_MATH_DEFINES
#include <cmath>
#include "uqff_waveform_simulate.h"
#include <iomanip>
#include <algorithm>
#include <numeric>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ═══════════════════════════════════════════════════════════════════════════════
// CONSTRUCTOR
// ═══════════════════════════════════════════════════════════════════════════════

UQFFWaveformSimulate::UQFFWaveformSimulate(unsigned int seed) 
    : G(6.6743e-11),
      c(2.99792458e8),
      M_sun(1.989e30),
      k_B(1.380649e-23),
      f_TRZ(0.1),
      B_t(1e-16),
      B_crit(1e11),
      rho_vac_UA(7.09e-36),
      alpha_UA(6.6743e-11 / (2.99792458e8 * 2.99792458e8)),
      U_m(1.0e-20),
      beta_m(0.01),
      T(2.725),
      rng(seed == 0 ? std::chrono::system_clock::now().time_since_epoch().count() : seed),
      noise_dist(0.0, 1.0) {
    
    // Initialize explanations
    explanations.push_back("UQFF Waveform Simulation with Chirp Evolution");
    explanations.push_back("═══════════════════════════════════════════════════════════");
    explanations.push_back("");
    explanations.push_back("PHYSICAL MODEL:");
    explanations.push_back("  Standard GR: h = (4G²μM)/(c⁴ar) × cos(2ωt)");
    explanations.push_back("  UQFF: h_UQFF = h × S_TRZ × S_SCm × S_aether × S_string × (1 + β_m sin(...))");
    explanations.push_back("");
    explanations.push_back("SUPPRESSION FACTORS:");
    explanations.push_back("  S_TRZ = (1 - f_TRZ) = 0.9 for f_TRZ = 0.1");
    explanations.push_back("  S_SCm = (1 - B_t/B_crit) ≈ 1.0 for stellar BHs");
    explanations.push_back("  S_aether = exp(-α_UA ρ_UA r/c) ≈ 1.0 for short r");
    explanations.push_back("  S_string = exp(-U_m/E_bind) ~ 0.37-1.0");
    explanations.push_back("");
    explanations.push_back("CHIRP EVOLUTION:");
    explanations.push_back("  ω(t) = ω₀ × (1 - t/τ_merge)^(-3/8)");
    explanations.push_back("  τ_merge = (5/256) × (c⁵a₀⁴)/(G³M²μ)");
    explanations.push_back("");
    explanations.push_back("GW150914-LIKE SIMULATION:");
    explanations.push_back("  μ = 15 M_sun, M_tot = 65 M_sun, a₀ = 100 km");
    explanations.push_back("  Expected: ~10-20% amplitude reduction from UQFF");
}

// ═══════════════════════════════════════════════════════════════════════════════
// STANDARD WAVEFORM
// ═══════════════════════════════════════════════════════════════════════════════

double UQFFWaveformSimulate::compute_h_standard(double mu, double M_tot, double a, 
                                                 double r_observer, double omega, double t) {
    /**
     * Standard quadrupole GW strain:
     *   h = (4G²μM)/(c⁴ar) × cos(2ωt)
     * 
     * Parameters:
     *   μ = reduced mass (kg)
     *   M_tot = total mass (kg)
     *   a = orbital separation (m)
     *   r_observer = distance to observer (m)
     *   ω = orbital angular frequency (rad/s)
     *   t = time (s)
     * 
     * Returns: dimensionless strain amplitude
     */
    if (a <= 0 || r_observer <= 0) return 0.0;
    
    double G2 = G * G;
    double c4 = c * c * c * c;
    double amplitude = (4.0 * G2 * mu * M_tot) / (c4 * a * r_observer);
    return amplitude * std::cos(2.0 * omega * t);
}

// ═══════════════════════════════════════════════════════════════════════════════
// UQFF SUPPRESSION FACTORS
// ═══════════════════════════════════════════════════════════════════════════════

double UQFFWaveformSimulate::compute_h_time_reversal(double h_std) {
    /**
     * Time-reversal zone suppression:
     *   S_TRZ = (1 - f_TRZ)
     * 
     * For f_TRZ = 0.1: S_TRZ = 0.9 (10% reduction)
     * Physical: Near-horizon f_TRZ couples to wave amplitude
     */
    return h_std * (1.0 - f_TRZ);
}

double UQFFWaveformSimulate::compute_h_superconducting(double h) {
    /**
     * Superconducting horizon screening:
     *   S_SCm = (1 - B_t/B_crit)
     * 
     * For stellar BHs: B_t ~ 10^-16 T, B_crit ~ 10^11 T
     *   → S_SCm ≈ 1.0 (negligible)
     * For magnetars: B_t ~ 10^15 T → significant screening
     */
    double ratio = B_t / B_crit;
    if (ratio >= 1.0) return 0.0;  // Full screening
    return h * (1.0 - ratio);
}

double UQFFWaveformSimulate::compute_h_aether(double h, double r_observer) {
    /**
     * Aether damping (cosmological):
     *   S_aether = exp(-α_UA × ρ_UA × r / c)
     * 
     * With α_UA = G/c² and ρ_UA = 7.09e-36 J/m³:
     *   Damping length ~ c/(α_UA × ρ_UA) >> Mpc
     *   → Negligible for local sources
     */
    double damping = alpha_UA * rho_vac_UA * r_observer / c;
    return h * std::exp(-damping);
}

double UQFFWaveformSimulate::compute_h_magnetic_string(double h, double M_tot, double a) {
    /**
     * String binding interference:
     *   S_string = exp(-U_m/E_bind)
     * 
     * E_bind = GM²/a (gravitational binding energy)
     * For tight binaries: E_bind >> U_m → S_string ≈ 1
     * 
     * When U_m ~ E_bind: S_string ~ e^-1 ≈ 0.37
     */
    if (a <= 0) return 0.0;
    
    double E_bind = G * M_tot * M_tot / a;
    if (E_bind <= 0) return h;
    
    double ratio = U_m / E_bind;
    return h * std::exp(-ratio);
}

double UQFFWaveformSimulate::compute_h_interference(double h, double omega) {
    /**
     * String interference modulation:
     *   Factor = (1 + β_m × sin(U_m × ω / (k_B × T)))
     * 
     * Creates ±β_m oscillations in amplitude
     * For β_m = 0.01: ±1% modulation
     */
    double phase = U_m * omega / (k_B * T);
    return h * (1.0 + beta_m * std::sin(phase));
}

// ═══════════════════════════════════════════════════════════════════════════════
// FULL UQFF WAVEFORM
// ═══════════════════════════════════════════════════════════════════════════════

double UQFFWaveformSimulate::compute_full_h_UQFF(double mu, double M_tot, double a,
                                                  double r_observer, double omega, double t,
                                                  double noise_level) {
    /**
     * Full UQFF waveform with all corrections:
     *   h_UQFF = h_std × S_TRZ × S_SCm × S_aether × S_string × (1 + β_m sin(...))
     * 
     * Pipeline:
     *   1. Standard quadrupole h
     *   2. Time-reversal suppression
     *   3. SCm horizon screening
     *   4. Aether damping
     *   5. String interference
     *   6. Modulation
     *   7. Additional user mods
     *   8. Noise (optional)
     */
    double h = compute_h_standard(mu, M_tot, a, r_observer, omega, t);
    h = compute_h_time_reversal(h);
    h = compute_h_superconducting(h);
    h = compute_h_aether(h, r_observer);
    h = compute_h_magnetic_string(h, M_tot, a);
    h = compute_h_interference(h, omega);
    
    // Apply additional modifications
    for (const auto& mod : additional_mods) {
        h *= mod(omega, t);
    }
    
    // Add noise if requested
    if (noise_level > 0.0) {
        h += noise_level * std::abs(h) * noise_dist(rng);
    }
    
    return h;
}

// ═══════════════════════════════════════════════════════════════════════════════
// CHIRP FREQUENCY EVOLUTION
// ═══════════════════════════════════════════════════════════════════════════════

double UQFFWaveformSimulate::compute_chirp_omega(double omega_0, double t, double tau_merge) {
    /**
     * Chirp frequency evolution (Peters formula):
     *   ω(t) = ω₀ × (1 - t/τ_merge)^(-3/8)
     * 
     * As t → τ_merge: ω → ∞ (merger)
     * 
     * τ_merge = (5/256) × (c⁵a₀⁴)/(G³M²μ) [Peters 1964]
     */
    double ratio = 1.0 - t / tau_merge;
    if (ratio <= 0.0) return omega_0 * 10.0;  // Cap at merger
    return omega_0 * std::pow(ratio, -3.0/8.0);
}

// ═══════════════════════════════════════════════════════════════════════════════
// SELF-EXPANDING FUNCTIONALITY
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFWaveformSimulate::add_mod(std::function<double(double, double)> mod) {
    additional_mods.push_back(mod);
}

void UQFFWaveformSimulate::update_from_file(const std::string& config_file) {
    std::ifstream infile(config_file);
    if (!infile.is_open()) {
        std::cerr << "Warning: Could not open config file: " << config_file << std::endl;
        return;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        size_t pos = line.find('=');
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            std::string val_str = line.substr(pos + 1);
            
            // Trim whitespace
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            val_str.erase(0, val_str.find_first_not_of(" \t"));
            val_str.erase(val_str.find_last_not_of(" \t") + 1);
            
            try {
                double value = std::stod(val_str);
                
                if (key == "G") G = value;
                else if (key == "c") c = value;
                else if (key == "f_TRZ") f_TRZ = value;
                else if (key == "B_t") B_t = value;
                else if (key == "B_crit") B_crit = value;
                else if (key == "U_m") U_m = value;
                else if (key == "beta_m") beta_m = value;
                else if (key == "T") T = value;
                else if (key == "rho_vac_UA") rho_vac_UA = value;
                else {
                    std::cerr << "Unknown parameter: " << key << std::endl;
                }
            } catch (const std::exception& e) {
                std::cerr << "Parse error for " << key << ": " << e.what() << std::endl;
            }
        }
    }
    infile.close();
    std::cout << "Loaded parameters from " << config_file << std::endl;
}

// ═══════════════════════════════════════════════════════════════════════════════
// WAVEFORM SIMULATION
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFWaveformSimulate::simulate_waveform(double mu, double M_tot, double a_initial,
                                              double r_observer, double omega_start,
                                              double t_start, double t_end, double dt,
                                              const std::string& output_file) {
    /**
     * Simulate waveform with linear chirp approximation.
     * For more accurate simulation, use simulate_inspiral().
     */
    if (dt <= 0 || t_end <= t_start) {
        std::cerr << "Invalid time parameters" << std::endl;
        return;
    }

    std::ofstream outfile;
    bool file_output = !output_file.empty();
    if (file_output) {
        outfile.open(output_file);
        if (!outfile.is_open()) {
            std::cerr << "Could not open output file: " << output_file << std::endl;
            return;
        }
        outfile << "time,omega,h_standard,h_UQFF,ratio" << std::endl;
    }

    // Compute merger timescale for chirp
    double tau_merge = (5.0/256.0) * std::pow(c, 5) * std::pow(a_initial, 4) / 
                       (std::pow(G, 3) * M_tot * M_tot * mu);
    
    std::vector<double> h_values;
    
    for (double t = t_start; t <= t_end; t += dt) {
        double omega = compute_chirp_omega(omega_start, t, tau_merge);
        
        // Evolve separation (Kepler approx)
        double a = a_initial * std::pow(1.0 - t/tau_merge, 1.0/4.0);
        if (a < 3.0 * G * M_tot / (c * c)) a = 3.0 * G * M_tot / (c * c);  // ISCO
        
        double h_std = compute_h_standard(mu, M_tot, a, r_observer, omega, t);
        double h_uqff = compute_full_h_UQFF(mu, M_tot, a, r_observer, omega, t);
        double ratio = (std::abs(h_std) > 1e-50) ? h_uqff / h_std : 1.0;
        
        h_values.push_back(h_uqff);
        
        if (file_output) {
            outfile << std::scientific << std::setprecision(6)
                    << t << "," << omega << "," << h_std << "," << h_uqff << "," << ratio << std::endl;
        }
    }
    
    if (file_output) {
        outfile.close();
        std::cout << "Waveform saved to " << output_file << std::endl;
        
        // Print statistics
        double max_h, min_h, mean_h, rms_h;
        get_statistics(h_values, max_h, min_h, mean_h, rms_h);
        std::cout << "Statistics: max=" << max_h << ", min=" << min_h 
                  << ", rms=" << rms_h << std::endl;
    }
}

void UQFFWaveformSimulate::simulate_inspiral(double M1, double M2, double a_initial,
                                              double r_observer, double t_duration, double dt,
                                              const std::string& output_file) {
    /**
     * Simulate inspiral with Peters orbital evolution.
     * More accurate than linear chirp.
     */
    double M_tot = M1 + M2;
    double mu = M1 * M2 / M_tot;
    
    // Merger timescale
    double tau_merge = (5.0/256.0) * std::pow(c, 5) * std::pow(a_initial, 4) / 
                       (std::pow(G, 3) * M_tot * M_tot * mu);
    
    // Initial orbital frequency
    double omega_0 = std::sqrt(G * M_tot / std::pow(a_initial, 3));
    
    std::cout << "Inspiral simulation:" << std::endl;
    std::cout << "  M1 = " << M1/M_sun << " M_sun" << std::endl;
    std::cout << "  M2 = " << M2/M_sun << " M_sun" << std::endl;
    std::cout << "  a_initial = " << a_initial/1e3 << " km" << std::endl;
    std::cout << "  tau_merge = " << tau_merge << " s" << std::endl;
    std::cout << "  f_initial = " << omega_0/(2*M_PI) << " Hz" << std::endl;
    
    simulate_waveform(mu, M_tot, a_initial, r_observer, omega_0, 0.0, t_duration, dt, output_file);
}

// ═══════════════════════════════════════════════════════════════════════════════
// UTILITIES
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFWaveformSimulate::display_explanations() {
    std::cout << std::endl;
    for (const auto& exp : explanations) {
        std::cout << exp << std::endl;
    }
    std::cout << std::endl;
}

void UQFFWaveformSimulate::get_statistics(const std::vector<double>& h_values,
                                           double& max_val, double& min_val,
                                           double& mean_val, double& rms_val) {
    if (h_values.empty()) {
        max_val = min_val = mean_val = rms_val = 0.0;
        return;
    }
    
    max_val = *std::max_element(h_values.begin(), h_values.end());
    min_val = *std::min_element(h_values.begin(), h_values.end());
    mean_val = std::accumulate(h_values.begin(), h_values.end(), 0.0) / h_values.size();
    
    double sq_sum = 0.0;
    for (double h : h_values) {
        sq_sum += h * h;
    }
    rms_val = std::sqrt(sq_sum / h_values.size());
}

// ═══════════════════════════════════════════════════════════════════════════════
// VALIDATION TESTS
// ═══════════════════════════════════════════════════════════════════════════════

int run_validation_tests() {
    std::cout << "\n╔══════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║        UQFF WAVEFORM SIMULATION - VALIDATION TESTS               ║" << std::endl;
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n" << std::endl;

    int passed = 0;
    int failed = 0;
    
    UQFFWaveformSimulate waveform(42);  // Fixed seed for reproducibility
    
    // Physical constants
    double M_sun = 1.989e30;
    double G = 6.6743e-11;
    double c = 2.99792458e8;

    // Test 1: Standard waveform amplitude (GW150914-like at 410 Mpc)
    std::cout << "Test 1: Standard waveform amplitude...";
    {
        double M1 = 36.0 * M_sun;       // 36 M_sun
        double M2 = 29.0 * M_sun;       // 29 M_sun
        double M_tot = M1 + M2;
        double mu = M1 * M2 / M_tot;    // ~16 M_sun
        double a = 350.0e3;             // 350 km (late inspiral)
        double r = 410e6 * 3.086e16;    // 410 Mpc in meters
        double omega = 35.0 * M_PI;     // f_orb = f_GW/2
        double t = 0.0;
        
        double h = waveform.compute_h_standard(mu, M_tot, a, r, omega, t);
        
        // Expected: h ~ 1e-21 at 410 Mpc (LIGO sensitivity)
        if (std::abs(h) > 1e-23 && std::abs(h) < 1e-19) {
            std::cout << " PASSED (h = " << h << ")" << std::endl;
            passed++;
        } else {
            std::cout << " FAILED (h = " << h << ", expected ~1e-21)" << std::endl;
            failed++;
        }
    }

    // Test 2: Time-reversal suppression
    std::cout << "Test 2: Time-reversal suppression (f_TRZ = 0.1)...";
    {
        waveform.set_f_TRZ(0.1);
        double h_in = 1.0e-21;
        double h_out = waveform.compute_h_time_reversal(h_in);
        double expected = 0.9e-21;
        
        if (std::abs(h_out - expected) / expected < 0.01) {
            std::cout << " PASSED (ratio = " << h_out/h_in << ")" << std::endl;
            passed++;
        } else {
            std::cout << " FAILED (got " << h_out << ", expected " << expected << ")" << std::endl;
            failed++;
        }
    }

    // Test 3: SCm screening negligible for stellar BHs
    std::cout << "Test 3: SCm screening negligible for B_t << B_crit...";
    {
        waveform.set_B_t(1e-16);
        waveform.set_B_crit(1e11);
        double h_in = 1.0e-21;
        double h_out = waveform.compute_h_superconducting(h_in);
        
        if (std::abs(h_out - h_in) / h_in < 1e-10) {
            std::cout << " PASSED (ratio = " << h_out/h_in << ")" << std::endl;
            passed++;
        } else {
            std::cout << " FAILED (got " << h_out << ")" << std::endl;
            failed++;
        }
    }

    // Test 4: SCm screening significant for magnetars
    std::cout << "Test 4: SCm screening for magnetar fields...";
    {
        waveform.set_B_t(1e10);  // Near B_crit
        waveform.set_B_crit(1e11);
        double h_in = 1.0e-21;
        double h_out = waveform.compute_h_superconducting(h_in);
        double expected_ratio = 0.9;  // 1 - 0.1
        
        if (std::abs(h_out/h_in - expected_ratio) < 0.01) {
            std::cout << " PASSED (ratio = " << h_out/h_in << ")" << std::endl;
            passed++;
        } else {
            std::cout << " FAILED (got ratio " << h_out/h_in << ", expected " << expected_ratio << ")" << std::endl;
            failed++;
        }
        waveform.set_B_t(1e-16);  // Reset
    }

    // Test 5: Aether damping negligible at short range
    std::cout << "Test 5: Aether damping negligible at short range...";
    {
        double h_in = 1.0e-21;
        double r = 1e3;  // 1 km
        double h_out = waveform.compute_h_aether(h_in, r);
        
        if (std::abs(h_out - h_in) / h_in < 1e-20) {
            std::cout << " PASSED (ratio = " << h_out/h_in << ")" << std::endl;
            passed++;
        } else {
            std::cout << " FAILED (got " << h_out << ")" << std::endl;
            failed++;
        }
    }

    // Test 6: String interference with strong binding
    std::cout << "Test 6: String interference (strong binding)...";
    {
        waveform.set_U_m(1e-20);  // Very small U_m
        double h_in = 1.0e-21;
        double M_tot = 65 * M_sun;
        double a = 100e3;
        double E_bind = G * M_tot * M_tot / a;  // ~10^47 J
        double h_out = waveform.compute_h_magnetic_string(h_in, M_tot, a);
        
        // U_m << E_bind → ratio ≈ 1
        if (std::abs(h_out - h_in) / h_in < 1e-10) {
            std::cout << " PASSED (E_bind = " << E_bind << " J)" << std::endl;
            passed++;
        } else {
            std::cout << " FAILED (got " << h_out << ")" << std::endl;
            failed++;
        }
    }

    // Test 7: Modulation oscillations
    std::cout << "Test 7: Interference modulation...";
    {
        waveform.set_beta_m(0.01);
        double h_in = 1.0e-21;
        double omega1 = 100.0;
        double omega2 = 200.0;
        double h_out1 = waveform.compute_h_interference(h_in, omega1);
        double h_out2 = waveform.compute_h_interference(h_in, omega2);
        
        // Should differ by ±1%
        double diff = std::abs((h_out1 - h_out2) / h_in);
        if (diff < 0.03 && diff > 0.0) {
            std::cout << " PASSED (variation = " << diff*100 << "%)" << std::endl;
            passed++;
        } else {
            std::cout << " FAILED (variation = " << diff*100 << "%)" << std::endl;
            failed++;
        }
    }

    // Test 8: Full UQFF pipeline
    std::cout << "Test 8: Full UQFF pipeline...";
    {
        double mu = 15.0 * M_sun;
        double M_tot = 65.0 * M_sun;
        double a = 100.0e3;
        double r = 1e9;
        double omega = 35.0 * 2 * M_PI;
        double t = 0.0;
        
        double h_std = waveform.compute_h_standard(mu, M_tot, a, r, omega, t);
        double h_uqff = waveform.compute_full_h_UQFF(mu, M_tot, a, r, omega, t);
        double ratio = h_uqff / h_std;
        
        // With defaults: should be ~0.9 from f_TRZ
        if (ratio > 0.8 && ratio < 1.0) {
            std::cout << " PASSED (ratio = " << ratio << ")" << std::endl;
            passed++;
        } else {
            std::cout << " FAILED (ratio = " << ratio << ", expected ~0.9)" << std::endl;
            failed++;
        }
    }

    // Test 9: Chirp frequency evolution
    std::cout << "Test 9: Chirp frequency evolution...";
    {
        double omega_0 = 100.0;
        double tau_merge = 1.0;
        
        double omega_0s = waveform.compute_chirp_omega(omega_0, 0.0, tau_merge);
        double omega_half = waveform.compute_chirp_omega(omega_0, 0.5, tau_merge);
        double omega_near = waveform.compute_chirp_omega(omega_0, 0.9, tau_merge);
        
        if (omega_0s < omega_half && omega_half < omega_near) {
            std::cout << " PASSED (ω increases: " << omega_0s << " → " << omega_half << " → " << omega_near << ")" << std::endl;
            passed++;
        } else {
            std::cout << " FAILED (chirp not increasing)" << std::endl;
            failed++;
        }
    }

    // Test 10: Self-expand with custom modification
    std::cout << "Test 10: Self-expand with custom mod...";
    {
        UQFFWaveformSimulate waveform2(42);
        
        double mu = 15.0 * M_sun;
        double M_tot = 65.0 * M_sun;
        double a = 100.0e3;
        double r = 1e9;
        double omega = 35.0 * 2 * M_PI;
        double t = 0.0;
        
        double h_before = waveform2.compute_full_h_UQFF(mu, M_tot, a, r, omega, t);
        
        // Add 50% reduction
        waveform2.add_mod([](double omega, double t) { return 0.5; });
        
        double h_after = waveform2.compute_full_h_UQFF(mu, M_tot, a, r, omega, t);
        double ratio = h_after / h_before;
        
        if (std::abs(ratio - 0.5) < 0.01) {
            std::cout << " PASSED (ratio = " << ratio << ")" << std::endl;
            passed++;
        } else {
            std::cout << " FAILED (ratio = " << ratio << ", expected 0.5)" << std::endl;
            failed++;
        }
    }

    // Test 11: Waveform simulation runs
    std::cout << "Test 11: Waveform simulation...";
    {
        try {
            UQFFWaveformSimulate waveform3(42);
            waveform3.simulate_waveform(15*M_sun, 65*M_sun, 100e3, 1e9, 
                                         35*2*M_PI, 0.0, 0.01, 0.001, "");
            std::cout << " PASSED (simulation completed)" << std::endl;
            passed++;
        } catch (...) {
            std::cout << " FAILED (exception thrown)" << std::endl;
            failed++;
        }
    }

    // Test 12: GW150914-like simulation
    std::cout << "Test 12: GW150914-like example...";
    {
        UQFFWaveformSimulate waveform4(42);
        
        double M1 = 36.0 * M_sun;   // ~36 M_sun
        double M2 = 29.0 * M_sun;   // ~29 M_sun
        double M_tot = M1 + M2;
        double mu = M1 * M2 / M_tot;  // ~16 M_sun
        double a = 350.0e3;            // Late inspiral
        double r = 410e6 * 3.086e16;   // 410 Mpc in meters
        double f_gw = 35.0;            // Peak frequency
        double omega = M_PI * f_gw;    // Orbital (half GW)
        
        double h_std = waveform4.compute_h_standard(mu, M_tot, a, r, omega, 0.0);
        double h_uqff = waveform4.compute_full_h_UQFF(mu, M_tot, a, r, omega, 0.0);
        
        // GW150914 observed h ~ 1e-21
        if (std::abs(h_std) > 1e-23 && std::abs(h_std) < 1e-19) {
            double ratio = h_uqff / h_std;
            std::cout << " PASSED" << std::endl;
            std::cout << "         h_std = " << h_std << std::endl;
            std::cout << "         h_UQFF = " << h_uqff << std::endl;
            std::cout << "         ratio = " << ratio << std::endl;
            passed++;
        } else {
            std::cout << " FAILED (h_std = " << h_std << ")" << std::endl;
            failed++;
        }
    }

    // Summary
    std::cout << "\n══════════════════════════════════════════════════════════════════" << std::endl;
    std::cout << "VALIDATION SUMMARY: " << passed << "/" << (passed + failed) << " tests passed" << std::endl;
    std::cout << "══════════════════════════════════════════════════════════════════\n" << std::endl;

    return failed;
}

// ═══════════════════════════════════════════════════════════════════════════════
// MAIN
// ═══════════════════════════════════════════════════════════════════════════════

int main() {
    // Run validation tests
    int failures = run_validation_tests();
    
    if (failures == 0) {
        // Demo simulation
        std::cout << "\n═══════════════════════════════════════════════════════════════" << std::endl;
        std::cout << "DEMO: GW150914-like inspiral simulation" << std::endl;
        std::cout << "═══════════════════════════════════════════════════════════════\n" << std::endl;
        
        UQFFWaveformSimulate waveform;
        waveform.display_explanations();
        
        // GW150914 parameters
        double M_sun = 1.989e30;
        double M1 = 36.0 * M_sun;
        double M2 = 29.0 * M_sun;
        double a_initial = 500.0e3;  // 500 km
        double r_observer = 410e6 * 3.086e16;  // 410 Mpc
        
        waveform.simulate_inspiral(M1, M2, a_initial, r_observer, 
                                    0.1, 0.0001, "gw150914_uqff_sim.csv");
        
        std::cout << "\nSimulation output saved to gw150914_uqff_sim.csv" << std::endl;
    }
    
    return failures;
}
