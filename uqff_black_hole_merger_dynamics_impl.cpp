// uqff_black_hole_merger_dynamics_impl.cpp
// Implementation of UQFF Black Hole Merger Dynamics Module
// Compile: g++ -std=c++20 -O2 -Wall -o uqff_bh_merger uqff_black_hole_merger_dynamics_impl.cpp

#include "uqff_black_hole_merger_dynamics.h"
#include <iomanip>
#include <sstream>
#include <cassert>

// ═══════════════════════════════════════════════════════════════════════════════
// CONSTRUCTOR
// ═══════════════════════════════════════════════════════════════════════════════

UQFFBlackHoleMergerDynamics::UQFFBlackHoleMergerDynamics(unsigned int seed) 
    : rng(seed), noise_dist(0.0, 1.0) {
    
    // Physical constants (CODATA 2018)
    G = 6.67430e-11;        // m³ kg⁻¹ s⁻²
    c = 2.99792458e8;       // m/s
    M_sun = 1.989e30;       // kg
    
    // UQFF parameters
    rho_vac_UA = 7.09e-36;  // J/m³
    B_crit = 4.4e13;        // T (magnetar critical field)
    f_TRZ = 0.1;            // Time-reversal factor
    U_m = 1e40;             // J (string binding energy - calibrated for BH scales)
    gamma = 5e-5;           // day⁻¹
    t_n = 0.5;              // Normalized time
    B_t = 1e12;             // T (typical binary field, < B_crit)
    
    populate_explanations();
}

// ═══════════════════════════════════════════════════════════════════════════════
// POPULATE EXPLANATIONS
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFBlackHoleMergerDynamics::populate_explanations() {
    explanations.clear();
    
    explanations.push_back(
        "═══════════════════════════════════════════════════════════════════════════════\n"
        "UQFF BLACK HOLE MERGER DYNAMICS\n"
        "Gravitational Wave Emission with Aether Damping Corrections\n"
        "═══════════════════════════════════════════════════════════════════════════════");
    
    explanations.push_back(
        "STEP 1: STANDARD GW POWER (Quadrupole Formula)\n"
        "────────────────────────────────────────────────────────────────────────────\n"
        "  P_GW = (32/5) × (G⁴/c⁵) × (μ² M_tot²) / a⁵\n"
        "\n"
        "where:\n"
        "  μ = M₁M₂/(M₁+M₂) = reduced mass\n"
        "  M_tot = M₁ + M₂ = total mass\n"
        "  a = orbital separation\n"
        "\n"
        "This is the Peters-Mathews formula for circular orbits.");
    
    explanations.push_back(
        "STEP 2: AETHER DAMPING\n"
        "────────────────────────────────────────────────────────────────────────────\n"
        "  P' = P_GW × exp(-ρ_UA × a × c² / (G × M_tot))\n"
        "\n"
        "Physical interpretation:\n"
        "  The [UA] aether medium provides viscous damping that reduces GW emission.\n"
        "  Damping increases with separation a (more aether to traverse).\n"
        "  Damping decreases with mass (stronger gravity overcomes medium).");
    
    explanations.push_back(
        "STEP 3: SUPERCONDUCTIVE HORIZON SUPPRESSION\n"
        "────────────────────────────────────────────────────────────────────────────\n"
        "  P'' = P' × (1 - B_t/B_crit)\n"
        "\n"
        "Physical interpretation:\n"
        "  When B_t → B_crit, the horizon becomes superconductive and GW emission\n"
        "  is suppressed. For typical binaries B_t << B_crit, so factor ≈ 1.");
    
    explanations.push_back(
        "STEP 4: TIME-REVERSAL NEGENTROPY\n"
        "────────────────────────────────────────────────────────────────────────────\n"
        "  P''' = P'' × (1 - f_TRZ)\n"
        "\n"
        "Physical interpretation:\n"
        "  f_TRZ ≈ 0.1 represents negentropic time-reversal component.\n"
        "  Reduces energy loss rate by ~10%, extending merger time.");
    
    explanations.push_back(
        "STEP 5: STRING BINDING\n"
        "────────────────────────────────────────────────────────────────────────────\n"
        "  P_UQFF = P''' × exp(-U_m / (G M_tot²/a))\n"
        "\n"
        "Physical interpretation:\n"
        "  U_m = magnetic string binding energy threading the binary.\n"
        "  Denominator G M_tot²/a is the gravitational binding energy.\n"
        "  When U_m ~ binding energy, significant suppression occurs.");
    
    explanations.push_back(
        "STEP 6: MODIFIED MERGER TIMESCALE\n"
        "────────────────────────────────────────────────────────────────────────────\n"
        "  τ_standard = (5/256) × (c⁵/G³) × (a⁴/(μ M_tot²))\n"
        "  τ_UQFF = τ_standard × (P_GW_standard / P_GW_UQFF)\n"
        "\n"
        "Since UQFF reduces power, merger takes longer:\n"
        "  τ_UQFF > τ_standard");
    
    explanations.push_back(
        "NUMERICAL EXAMPLE: GW150914-like\n"
        "────────────────────────────────────────────────────────────────────────────\n"
        "  M₁ = 36 M_sun, M₂ = 29 M_sun, a = 10¹¹ m\n"
        "  μ = 16.0 M_sun, M_tot = 65 M_sun\n"
        "  \n"
        "  UQFF factors:\n"
        "    Aether: exp(-ρ_UA a c²/(G M_tot)) ≈ 1 (negligible)\n"
        "    Horizon: (1 - B_t/B_crit) ≈ 0.9\n"
        "    f_TRZ: (1 - 0.1) = 0.9\n"
        "    String: exp(-U_m a/(G M_tot²)) varies\n"
        "  \n"
        "  Combined: τ_UQFF ≈ 1.2-3× longer than standard");
    
    explanations.push_back(
        "Q-SCOPE TESTABILITY:\n"
        "────────────────────────────────────────────────────────────────────────────\n"
        "  • Measure LIGO/Virgo waveform phase deviations from GR templates\n"
        "  • Search for anomalous inspiral rates in binary population\n"
        "  • Compare merger rates to standard predictions\n"
        "  • Detect mass retention signatures (M_final > expected)\n"
        "═══════════════════════════════════════════════════════════════════════════════");
}

// ═══════════════════════════════════════════════════════════════════════════════
// CORE COMPUTATIONS
// ═══════════════════════════════════════════════════════════════════════════════

double UQFFBlackHoleMergerDynamics::compute_mu(double M1, double M2) {
    // μ = M₁M₂ / (M₁ + M₂)
    return (M1 * M2) / (M1 + M2);
}

double UQFFBlackHoleMergerDynamics::compute_P_GW_standard(double mu, double M_tot, double a) {
    // P_GW = (32/5) × (G⁴/c⁵) × (μ² M_tot²) / a⁵
    // Peters-Mathews formula for circular orbits
    double G4 = G * G * G * G;
    double c5 = c * c * c * c * c;
    double a5 = a * a * a * a * a;
    return (32.0 / 5.0) * (G4 / c5) * (mu * mu * M_tot * M_tot) / a5;
}

double UQFFBlackHoleMergerDynamics::compute_P_GW_prime(double P_GW, double a, double M_tot) {
    // P' = P_GW × exp(-ρ_UA × a × c² / (G × M_tot))
    // Aether damping term
    double exponent = -(rho_vac_UA * a * c * c) / (G * M_tot);
    // Prevent underflow
    if (exponent < -700) return 0.0;
    return P_GW * std::exp(exponent);
}

double UQFFBlackHoleMergerDynamics::compute_P_GW_double_prime(double P_GW_prime) {
    // P'' = P' × (1 - B_t/B_crit)
    // Superconductive horizon suppression
    double factor = 1.0 - B_t / B_crit;
    if (factor < 0) factor = 0;  // Can't have negative power
    return P_GW_prime * factor;
}

double UQFFBlackHoleMergerDynamics::compute_P_GW_triple_prime(double P_GW_double_prime) {
    // P''' = P'' × (1 - f_TRZ)
    // Time-reversal negentropy
    return P_GW_double_prime * (1.0 - f_TRZ);
}

double UQFFBlackHoleMergerDynamics::compute_P_GW_UQFF(double P_GW_triple_prime, double M_tot, double a) {
    // P_UQFF = P''' × exp(-U_m / (G M_tot²/a))
    // String binding suppression
    double binding_energy = G * M_tot * M_tot / a;
    if (binding_energy == 0) return 0.0;
    
    double exponent = -U_m / binding_energy;
    // Prevent underflow
    if (exponent < -700) return 0.0;
    return P_GW_triple_prime * std::exp(exponent);
}

double UQFFBlackHoleMergerDynamics::compute_full_P_GW_UQFF(double M1, double M2, double a, double noise_level) {
    // Complete UQFF GW power computation
    double mu = compute_mu(M1, M2);
    double M_tot = M1 + M2;
    
    double P_GW = compute_P_GW_standard(mu, M_tot, a);
    double P_prime = compute_P_GW_prime(P_GW, a, M_tot);
    double P_double_prime = compute_P_GW_double_prime(P_prime);
    double P_triple_prime = compute_P_GW_triple_prime(P_double_prime);
    double P_uqff = compute_P_GW_UQFF(P_triple_prime, M_tot, a);
    
    // Apply additional mods
    for (const auto& mod : additional_mods) {
        P_uqff *= mod(M_tot, a);
    }
    
    // Add noise if requested
    if (noise_level > 0) {
        double noise = noise_level * noise_dist(rng);
        P_uqff += noise;
    }
    
    return P_uqff;
}

double UQFFBlackHoleMergerDynamics::compute_tau_merge_standard(double a, double mu, double M_tot) {
    // τ_merge = (5/256) × (c⁵/G³) × (a⁴/(μ M_tot²))
    double c5 = c * c * c * c * c;
    double G3 = G * G * G;
    double a4 = a * a * a * a;
    return (5.0 / 256.0) * (c5 / G3) * (a4 / (mu * M_tot * M_tot));
}

double UQFFBlackHoleMergerDynamics::compute_tau_merge_UQFF(double tau_std, double P_GW_uqff, double P_GW_std) {
    // τ_UQFF = τ_standard × (P_GW_standard / P_GW_UQFF)
    if (P_GW_uqff <= 0) return std::numeric_limits<double>::infinity();
    return tau_std * (P_GW_std / P_GW_uqff);
}

// ═══════════════════════════════════════════════════════════════════════════════
// ADDITIONAL COMPUTATIONS
// ═══════════════════════════════════════════════════════════════════════════════

double UQFFBlackHoleMergerDynamics::compute_f_orbital(double M_tot, double a) {
    // Kepler's law: f_orbital = (1/2π) × sqrt(G M_tot / a³)
    double a3 = a * a * a;
    return (1.0 / (2.0 * M_PI)) * std::sqrt(G * M_tot / a3);
}

double UQFFBlackHoleMergerDynamics::compute_f_GW(double M_tot, double a) {
    // GW frequency is twice orbital frequency
    return 2.0 * compute_f_orbital(M_tot, a);
}

double UQFFBlackHoleMergerDynamics::compute_chirp_mass(double M1, double M2) {
    // M_chirp = (M₁ M₂)^(3/5) / (M₁ + M₂)^(1/5)
    double M_tot = M1 + M2;
    return std::pow(M1 * M2, 0.6) / std::pow(M_tot, 0.2);
}

double UQFFBlackHoleMergerDynamics::compute_E_radiated(double M1, double M2, double a_initial, double a_final) {
    // E_radiated = G M₁ M₂ × (1/a_final - 1/a_initial) / 2
    return G * M1 * M2 * (1.0/a_final - 1.0/a_initial) / 2.0;
}

// ═══════════════════════════════════════════════════════════════════════════════
// SELF-EXPANDING CAPABILITIES
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFBlackHoleMergerDynamics::add_mod(std::function<double(double, double)> mod) {
    additional_mods.push_back(mod);
}

void UQFFBlackHoleMergerDynamics::update_from_file(const std::string& config_file) {
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
            else if (key == "rho_vac_UA") rho_vac_UA = value;
            else if (key == "B_crit") B_crit = value;
            else if (key == "f_TRZ") f_TRZ = value;
            else if (key == "U_m") U_m = value;
            else if (key == "gamma") gamma = value;
            else if (key == "t_n") t_n = value;
            else if (key == "B_t") B_t = value;
        }
    }
    infile.close();
    std::cout << "Updated parameters from " << config_file << std::endl;
}

void UQFFBlackHoleMergerDynamics::simulate_merger(double M1, double M2, double a_start, 
                                                   double t_start, double t_end, double dt,
                                                   const std::string& output_file) {
    std::ofstream outfile;
    bool file_output = !output_file.empty();
    if (file_output) {
        outfile.open(output_file);
        outfile << "time_s,a_m,P_GW_UQFF_W,P_GW_std_W,f_GW_Hz\n";
    }

    double a_current = a_start;
    double M_tot = M1 + M2;
    double mu = compute_mu(M1, M2);
    
    for (double t = t_start; t <= t_end && a_current > 0; t += dt) {
        double P_GW_std = compute_P_GW_standard(mu, M_tot, a_current);
        double P_GW_uqff = compute_full_P_GW_UQFF(M1, M2, a_current);
        double f_GW = compute_f_GW(M_tot, a_current);
        
        // Energy loss: dE/dt = P_GW
        // Orbital energy E = -G M₁ M₂ / (2a)
        // dE/da = G M₁ M₂ / (2a²)
        // da/dt = dE/dt / (dE/da) = -P_GW × 2a² / (G M₁ M₂)
        double da_dt = -P_GW_uqff * 2.0 * a_current * a_current / (G * M1 * M2);
        
        if (file_output) {
            outfile << std::scientific << std::setprecision(6);
            outfile << t << "," << a_current << "," << P_GW_uqff << "," 
                    << P_GW_std << "," << f_GW << "\n";
        }
        
        a_current += da_dt * dt;
        
        // Stop if separation gets too small
        double r_s = 2 * G * M_tot / (c * c);
        if (a_current < 3 * r_s) break;  // ISCO approximately
    }

    if (file_output) {
        outfile.close();
        std::cout << "Merger simulation saved to: " << output_file << std::endl;
    }
}

void UQFFBlackHoleMergerDynamics::display_explanations() {
    for (const auto& exp : explanations) {
        std::cout << exp << std::endl << std::endl;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// VALIDATION TESTS
// ═══════════════════════════════════════════════════════════════════════════════

void run_validation_tests() {
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "               UQFF BLACK HOLE MERGER DYNAMICS - VALIDATION TESTS              \n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n\n";
    
    UQFFBlackHoleMergerDynamics merger(42);  // Fixed seed for reproducibility
    
    int tests_passed = 0;
    int tests_total = 0;
    
    double M_sun = merger.get_M_sun();
    double G = merger.get_G();
    double c = merger.get_c();
    
    // Test parameters (GW150914-like)
    double M1 = 36.0 * M_sun;
    double M2 = 29.0 * M_sun;
    double M_tot = M1 + M2;
    double a = 1e9;  // 10^9 m separation
    
    // ─────────────────────────────────────────────────────────────────────────────
    // TEST 1: Reduced mass computation
    // ─────────────────────────────────────────────────────────────────────────────
    tests_total++;
    double mu = merger.compute_mu(M1, M2);
    double mu_expected = M1 * M2 / M_tot;
    double mu_solar = mu / M_sun;
    bool t1 = std::abs(mu - mu_expected) / mu_expected < 1e-10;
    std::cout << "TEST 1: Reduced mass computation\n";
    std::cout << "  M₁ = 36 M_sun, M₂ = 29 M_sun\n";
    std::cout << "  μ = " << mu_solar << " M_sun\n";
    std::cout << "  Expected: " << mu_expected/M_sun << " M_sun\n";
    std::cout << "  Result: " << (t1 ? "PASS ✓" : "FAIL ✗") << "\n\n";
    if (t1) tests_passed++;
    
    // ─────────────────────────────────────────────────────────────────────────────
    // TEST 2: Chirp mass
    // ─────────────────────────────────────────────────────────────────────────────
    tests_total++;
    double M_chirp = merger.compute_chirp_mass(M1, M2);
    double M_chirp_expected = std::pow(M1 * M2, 0.6) / std::pow(M_tot, 0.2);
    bool t2 = std::abs(M_chirp - M_chirp_expected) / M_chirp_expected < 1e-10;
    std::cout << "TEST 2: Chirp mass\n";
    std::cout << "  M_chirp = " << M_chirp/M_sun << " M_sun\n";
    std::cout << "  GW150914 observed: ~28 M_sun\n";
    std::cout << "  Result: " << (t2 ? "PASS ✓" : "FAIL ✗") << "\n\n";
    if (t2) tests_passed++;
    
    // ─────────────────────────────────────────────────────────────────────────────
    // TEST 3: Standard GW power (order of magnitude)
    // ─────────────────────────────────────────────────────────────────────────────
    tests_total++;
    double P_GW_std = merger.compute_P_GW_standard(mu, M_tot, a);
    // At a=10^9 m (wide separation), P ~ 1 W; at a=350 km (detection), P ~ 10^17 W
    bool t3 = P_GW_std > 0 && P_GW_std < 1e55;
    std::cout << "TEST 3: Standard GW power (positive, finite)\n";
    std::cout << "  a = 10^9 m\n";
    std::cout << "  P_GW = " << P_GW_std << " W\n";
    std::cout << "  log₁₀(P_GW) = " << std::log10(P_GW_std) << "\n";
    std::cout << "  Expected: > 0 and < 10^55 W\n";
    std::cout << "  Result: " << (t3 ? "PASS ✓" : "FAIL ✗") << "\n\n";
    if (t3) tests_passed++;
    
    // ─────────────────────────────────────────────────────────────────────────────
    // TEST 4: P_GW scales as a^(-5)
    // ─────────────────────────────────────────────────────────────────────────────
    tests_total++;
    double a2 = 2 * a;
    double P_GW_2a = merger.compute_P_GW_standard(mu, M_tot, a2);
    double ratio_expected = std::pow(0.5, 5);  // (a/2a)^5 = 1/32
    double ratio_actual = P_GW_2a / P_GW_std;
    bool t4 = std::abs(ratio_actual - ratio_expected) / ratio_expected < 1e-6;
    std::cout << "TEST 4: P_GW ∝ a^(-5) scaling\n";
    std::cout << "  P(2a)/P(a) = " << ratio_actual << "\n";
    std::cout << "  Expected: " << ratio_expected << " = 1/32\n";
    std::cout << "  Result: " << (t4 ? "PASS ✓" : "FAIL ✗") << "\n\n";
    if (t4) tests_passed++;
    
    // ─────────────────────────────────────────────────────────────────────────────
    // TEST 5: f_TRZ reduces power by (1-f_TRZ)
    // ─────────────────────────────────────────────────────────────────────────────
    tests_total++;
    double f_TRZ = merger.get_f_TRZ();
    double P_prime = merger.compute_P_GW_prime(P_GW_std, a, M_tot);
    double P_double = merger.compute_P_GW_double_prime(P_prime);
    double P_triple = merger.compute_P_GW_triple_prime(P_double);
    // P_triple = P_double × (1 - f_TRZ)
    double f_TRZ_measured = 1.0 - P_triple / P_double;
    bool t5 = std::abs(f_TRZ_measured - f_TRZ) < 1e-10;
    std::cout << "TEST 5: f_TRZ = 0.1 reduction verification\n";
    std::cout << "  f_TRZ input: " << f_TRZ << "\n";
    std::cout << "  f_TRZ measured: " << f_TRZ_measured << "\n";
    std::cout << "  P''' / P'' = " << P_triple/P_double << "\n";
    std::cout << "  Result: " << (t5 ? "PASS ✓" : "FAIL ✗") << "\n\n";
    if (t5) tests_passed++;
    
    // ─────────────────────────────────────────────────────────────────────────────
    // TEST 6: GW frequency = 2 × orbital frequency
    // ─────────────────────────────────────────────────────────────────────────────
    tests_total++;
    double f_orb = merger.compute_f_orbital(M_tot, a);
    double f_GW = merger.compute_f_GW(M_tot, a);
    bool t6 = std::abs(f_GW - 2*f_orb) < 1e-10;
    std::cout << "TEST 6: f_GW = 2 × f_orbital\n";
    std::cout << "  f_orbital = " << f_orb << " Hz\n";
    std::cout << "  f_GW = " << f_GW << " Hz\n";
    std::cout << "  Ratio: " << f_GW/f_orb << "\n";
    std::cout << "  Result: " << (t6 ? "PASS ✓" : "FAIL ✗") << "\n\n";
    if (t6) tests_passed++;
    
    // ─────────────────────────────────────────────────────────────────────────────
    // TEST 7: Merger time τ_merge > 0
    // ─────────────────────────────────────────────────────────────────────────────
    tests_total++;
    double tau_std = merger.compute_tau_merge_standard(a, mu, M_tot);
    bool t7 = tau_std > 0;
    double tau_years = tau_std / (365.25 * 24 * 3600);
    std::cout << "TEST 7: Standard merger time > 0\n";
    std::cout << "  τ_merge = " << tau_std << " s\n";
    std::cout << "  τ_merge = " << tau_years << " years\n";
    std::cout << "  Result: " << (t7 ? "PASS ✓" : "FAIL ✗") << "\n\n";
    if (t7) tests_passed++;
    
    // ─────────────────────────────────────────────────────────────────────────────
    // TEST 8: UQFF extends merger time (τ_UQFF ≥ τ_std)
    // ─────────────────────────────────────────────────────────────────────────────
    tests_total++;
    double P_UQFF = merger.compute_full_P_GW_UQFF(M1, M2, a);
    double tau_UQFF = merger.compute_tau_merge_UQFF(tau_std, P_UQFF, P_GW_std);
    bool t8 = tau_UQFF >= tau_std;
    double tau_ratio = tau_UQFF / tau_std;
    std::cout << "TEST 8: UQFF extends merger time\n";
    std::cout << "  τ_standard = " << tau_std << " s\n";
    std::cout << "  τ_UQFF = " << tau_UQFF << " s\n";
    std::cout << "  τ_UQFF / τ_standard = " << tau_ratio << "\n";
    std::cout << "  Result: " << (t8 ? "PASS ✓" : "FAIL ✗") << "\n\n";
    if (t8) tests_passed++;
    
    // ─────────────────────────────────────────────────────────────────────────────
    // TEST 9: B_t < B_crit check (standard operation)
    // ─────────────────────────────────────────────────────────────────────────────
    tests_total++;
    double B_t = merger.get_B_t();
    double B_crit = merger.get_B_crit();
    bool t9 = B_t < B_crit;
    std::cout << "TEST 9: B_t < B_crit (non-superconductive regime)\n";
    std::cout << "  B_t = " << B_t << " T\n";
    std::cout << "  B_crit = " << B_crit << " T\n";
    std::cout << "  B_t/B_crit = " << B_t/B_crit << "\n";
    std::cout << "  Result: " << (t9 ? "PASS ✓" : "FAIL ✗") << "\n\n";
    if (t9) tests_passed++;
    
    // ─────────────────────────────────────────────────────────────────────────────
    // TEST 10: Energy radiated > 0
    // ─────────────────────────────────────────────────────────────────────────────
    tests_total++;
    double r_ISCO = 6 * G * M_tot / (c * c);  // ISCO for Schwarzschild
    double E_rad = merger.compute_E_radiated(M1, M2, a, r_ISCO);
    bool t10 = E_rad > 0;
    double E_rad_solar_c2 = E_rad / (M_sun * c * c);
    std::cout << "TEST 10: Energy radiated > 0\n";
    std::cout << "  a_initial = " << a << " m\n";
    std::cout << "  a_final = r_ISCO = " << r_ISCO << " m\n";
    std::cout << "  E_radiated = " << E_rad << " J\n";
    std::cout << "  E_radiated = " << E_rad_solar_c2 << " M_sun c²\n";
    std::cout << "  Result: " << (t10 ? "PASS ✓" : "FAIL ✗") << "\n\n";
    if (t10) tests_passed++;
    
    // ─────────────────────────────────────────────────────────────────────────────
    // TEST 11: τ_merge ∝ a⁴ scaling
    // ─────────────────────────────────────────────────────────────────────────────
    tests_total++;
    double tau_2a = merger.compute_tau_merge_standard(2*a, mu, M_tot);
    double tau_ratio_a4 = tau_2a / tau_std;
    double expected_ratio = 16.0;  // 2^4 = 16
    bool t11 = std::abs(tau_ratio_a4 - expected_ratio) / expected_ratio < 1e-6;
    std::cout << "TEST 11: τ_merge ∝ a⁴ scaling\n";
    std::cout << "  τ(2a)/τ(a) = " << tau_ratio_a4 << "\n";
    std::cout << "  Expected: " << expected_ratio << " = 2⁴\n";
    std::cout << "  Result: " << (t11 ? "PASS ✓" : "FAIL ✗") << "\n\n";
    if (t11) tests_passed++;
    
    // ─────────────────────────────────────────────────────────────────────────────
    // TEST 12: UQFF power < standard power
    // ─────────────────────────────────────────────────────────────────────────────
    tests_total++;
    bool t12 = P_UQFF < P_GW_std;
    std::cout << "TEST 12: UQFF power < Standard power\n";
    std::cout << "  P_GW_std = " << P_GW_std << " W\n";
    std::cout << "  P_GW_UQFF = " << P_UQFF << " W\n";
    std::cout << "  Ratio: " << P_UQFF/P_GW_std << "\n";
    std::cout << "  Result: " << (t12 ? "PASS ✓" : "FAIL ✗") << "\n\n";
    if (t12) tests_passed++;
    
    // ═══════════════════════════════════════════════════════════════════════════════
    // SUMMARY
    // ═══════════════════════════════════════════════════════════════════════════════
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "                              VALIDATION SUMMARY                               \n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "  Tests passed: " << tests_passed << " / " << tests_total << "\n";
    std::cout << "  Status: " << (tests_passed == tests_total ? "ALL TESTS PASSED ✓" : "SOME TESTS FAILED ✗") << "\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n\n";
}

// ═══════════════════════════════════════════════════════════════════════════════
// MAIN
// ═══════════════════════════════════════════════════════════════════════════════

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    // Display module explanations
    UQFFBlackHoleMergerDynamics merger;
    merger.display_explanations();
    
    // Run validation tests
    run_validation_tests();
    
    // Example usage: GW150914-like system
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "                    EXAMPLE: GW150914-LIKE BINARY                              \n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n\n";
    
    double M_sun = merger.get_M_sun();
    (void)merger.get_G();  // Available for future use
    (void)merger.get_c();  // Available for future use
    
    double M1 = 36.0 * M_sun;
    double M2 = 29.0 * M_sun;
    double M_tot = M1 + M2;
    double mu = merger.compute_mu(M1, M2);
    
    // At LIGO detection (a ~ 350 km)
    double a_detect = 350e3;  // 350 km
    
    double f_GW = merger.compute_f_GW(M_tot, a_detect);
    double P_std = merger.compute_P_GW_standard(mu, M_tot, a_detect);
    double P_UQFF = merger.compute_full_P_GW_UQFF(M1, M2, a_detect);
    double M_chirp = merger.compute_chirp_mass(M1, M2);
    
    std::cout << "Binary Parameters:\n";
    std::cout << "  M₁ = " << M1/M_sun << " M_sun\n";
    std::cout << "  M₂ = " << M2/M_sun << " M_sun\n";
    std::cout << "  M_chirp = " << M_chirp/M_sun << " M_sun\n";
    std::cout << "  μ = " << mu/M_sun << " M_sun\n";
    std::cout << "  a (at detection) = " << a_detect/1e3 << " km\n\n";
    
    std::cout << "Gravitational Wave Properties:\n";
    std::cout << "  f_GW = " << f_GW << " Hz\n";
    std::cout << "  P_GW (standard) = " << P_std << " W\n";
    std::cout << "  P_GW (UQFF) = " << P_UQFF << " W\n";
    std::cout << "  UQFF ratio: " << P_UQFF/P_std << "\n\n";
    
    // Run simulation
    std::cout << "Running merger simulation...\n";
    merger.simulate_merger(M1, M2, 1e8, 0, 1e6, 100, "uqff_bh_merger_sim.csv");
    
    std::cout << "\nModule ready for MAIN_1_CoAnQi.cpp integration.\n";
    
    return 0;
}
