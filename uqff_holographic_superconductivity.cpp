// uqff_holographic_superconductivity.cpp
// UQFF Holographic Superconductivity Module Implementation
// Compiles: g++ -std=c++17 -O2 uqff_holographic_superconductivity.cpp -o uqff_holographic_superconductivity.exe

#include "uqff_holographic_superconductivity.h"

// ═══════════════════════════════════════════════════════════════════════════════
// CONSTRUCTOR
// ═══════════════════════════════════════════════════════════════════════════════

UQFFHolographicSuperconductivity::UQFFHolographicSuperconductivity(unsigned int seed)
    : rng(seed == 0 ? std::chrono::system_clock::now().time_since_epoch().count() : seed),
      noise_dist(0.0, 1.0) {
    
    // Physical constants
    hbar = 1.0545718e-34;     // J·s
    c = 2.998e8;              // m/s
    G = 6.6743e-11;           // m³ kg⁻¹ s⁻²
    k_B = 1.380649e-23;       // J/K
    
    // UQFF parameters
    rho_vac_UA = 7.09e-36;    // UA vacuum density J/m³
    rho_vac_SCm = 7.09e-37;   // SCm vacuum density J/m³
    f_TRZ = 0.1;              // Time-reversal factor (10%)
    B_crit = 1e11;            // Critical field for holographic SC (T)
    mu_j = 1e15;              // Magnetic string tension J·m
    gamma = 5e-5;             // String evolution rate 1/day → convert to 1/s
    
    // Superconductor parameters (cuprate-like defaults)
    d = 3.0;                  // CFT₃ (boundary dimension)
    m_scalar = 1e-30;         // Scalar mass (kg) - determines Δ± = (d ± √(d² + 4m²L²))/2
    lambda_sc = 1e-10;        // Quartic coupling
    T_c_base = 100.0;         // Base critical temperature (K) - cuprate scale
    Delta_base = 3.5 * k_B * T_c_base;  // BCS-like gap: Δ ≈ 3.5 k_B T_c
    mu_chem = k_B * T_c_base; // Chemical potential ~ T_c
    
    // UQFF scaling
    kappa_UQFF = 1e-60;
    lambda_UQFF = 1e-9;
    T_eff_floor = 1e16;
    
    // Populate explanations
    populate_explanations();
}

// ═══════════════════════════════════════════════════════════════════════════════
// POPULATE EXPLANATIONS
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFHolographicSuperconductivity::populate_explanations() {
    explanations.clear();
    
    explanations.push_back("=== UQFF Holographic Superconductivity Derivation ===");
    explanations.push_back("");
    explanations.push_back("STANDARD HOLOGRAPHIC SUPERCONDUCTIVITY (Gubser 2008, HHH 2008):");
    explanations.push_back("  AdS/CFT models high-T_c superconductors via gravity dual");
    explanations.push_back("  Couple AdS gravity to Maxwell + charged scalar (hair):");
    explanations.push_back("");
    explanations.push_back("  ACTION (AdS₄/CFT₃ s-wave):");
    explanations.push_back("  S = ∫ √(-g) [R + d(d-1)/L² - ¼F² - |Dψ|² - m²|ψ|²] d^{d+1}x");
    explanations.push_back("");
    explanations.push_back("  PHASE TRANSITION:");
    explanations.push_back("    T > T_c: ⟨ψ⟩ = 0 (normal phase)");
    explanations.push_back("    T < T_c: ⟨ψ⟩ ≠ 0 (superconducting - scalar condenses)");
    explanations.push_back("    T_c ~ μ (chemical potential sets scale)");
    explanations.push_back("");
    explanations.push_back("  OPTICAL CONDUCTIVITY:");
    explanations.push_back("    σ(ω) = n_s/m × δ(ω) + regular (superconducting pole)");
    explanations.push_back("    Gap 2Δ appears in optical response");
    explanations.push_back("");
    explanations.push_back("UQFF INTEGRATION:");
    explanations.push_back("  • AdS bulk = [UA] superfluid aether");
    explanations.push_back("  • CFT boundary = [SCm] superconducting medium");
    explanations.push_back("  • Order parameter ψ = aether flux condensate");
    explanations.push_back("");
    explanations.push_back("Step 1: Standard Action (Reference)");
    explanations.push_back("  S = ∫√(-g) [R + d(d-1)/L² - ¼F² - |Dψ|² - m²|ψ|²] d^{d+1}x");
    explanations.push_back("  For d=3 (CFT₃): yields 2+1 dimensional superconductor");
    explanations.push_back("");
    explanations.push_back("Step 2: UQFF Holographic Scale");
    explanations.push_back("  L_UQFF = (ℏc/ρ_UA)^{1/4}");
    explanations.push_back("  Aether density sets AdS radius");
    explanations.push_back("");
    explanations.push_back("Step 3: Modified Potential");
    explanations.push_back("  V_UQFF(ψ) = m²|ψ|² + (λ/2)|ψ|⁴ × (1 - B_t/B_crit)");
    explanations.push_back("  Magnetic field suppresses quartic term → phase transition");
    explanations.push_back("");
    explanations.push_back("Step 4: Enhanced Critical Temperature");
    explanations.push_back("  T_c,UQFF = T_c × (1 + f_TRZ)");
    explanations.push_back("  Time-reversal symmetry enhances T_c by 10%");
    explanations.push_back("");
    explanations.push_back("Step 5: String-Damped Gap");
    explanations.push_back("  Δ_UQFF = Δ × exp(-U_m/(k_B×T))");
    explanations.push_back("  Where U_m = (μ_j/L)(1 - exp(-γt cos(πt_n)))");
    explanations.push_back("  Magnetic strings reduce superconducting gap");
    explanations.push_back("");
    explanations.push_back("Step 6: Full UQFF Action");
    explanations.push_back("  S_UQFF = S + ∫√(-g) U_m (1 + f_TRZ) d⁴x");
    explanations.push_back("  String energy adds to bulk action");
    explanations.push_back("");
    explanations.push_back("NUMERICAL EXAMPLE (Cuprate):");
    explanations.push_back("  T_c = 100 K, f_TRZ = 0.1 → T_c,UQFF = 110 K");
    explanations.push_back("  U_m/(k_B T) = 0.5 → Δ_UQFF = Δ × e^{-0.5} ≈ 0.6Δ");
    explanations.push_back("");
    explanations.push_back("Q-SCOPE TESTABILITY:");
    explanations.push_back("  • Measure T_c enhancement in thin films (~10% from f_TRZ)");
    explanations.push_back("  • THz spectroscopy: Gap reduction at high B_t");
    explanations.push_back("  • Boundary coherence length matches L_UQFF");
}

// ═══════════════════════════════════════════════════════════════════════════════
// CORE PHYSICS COMPUTATIONS
// ═══════════════════════════════════════════════════════════════════════════════

double UQFFHolographicSuperconductivity::compute_ds2_AdS(double L, double z, double dt, double dx, double dz) {
    // AdS metric: ds² = (L²/z²)(-dt² + dx² + dz²)
    if (z <= 0) return 0.0;
    double factor = (L * L) / (z * z);
    return factor * (-dt * dt + dx * dx + dz * dz);
}

double UQFFHolographicSuperconductivity::compute_L_UQFF() {
    // L_UQFF = (ℏc/ρ_UA)^{1/4}
    return std::pow(hbar * c / rho_vac_UA, 0.25);
}

double UQFFHolographicSuperconductivity::compute_standard_action(double R, double L, double F_sq, double D_psi_sq, double psi) {
    // S = ∫√(-g) [R + d(d-1)/L² - ¼F² - |Dψ|² - m²|ψ|²] d^{d+1}x
    // Return integrand (symbolic: √(-g) = 1 for flat slice)
    if (L <= 0) return 0.0;
    double cosmological = d * (d - 1.0) / (L * L);
    double maxwell = -0.25 * F_sq;
    double kinetic = -D_psi_sq;
    double mass_term = -m_scalar * m_scalar * psi * psi;
    return R + cosmological + maxwell + kinetic + mass_term;
}

double UQFFHolographicSuperconductivity::compute_V_UQFF(double psi, double B_t) {
    // V_UQFF(ψ) = m²|ψ|² + (λ/2)|ψ|⁴ × (1 - B_t/B_crit)
    // Magnetic field suppresses quartic coupling → weakens order
    double B_ratio = B_t / B_crit;
    B_ratio = std::min(B_ratio, 1.0);  // Cap at 1
    double quartic_factor = 1.0 - B_ratio;
    return m_scalar * m_scalar * psi * psi + 
           (lambda_sc / 2.0) * std::pow(psi, 4) * quartic_factor;
}

double UQFFHolographicSuperconductivity::compute_T_c_UQFF() {
    // T_c,UQFF = T_c × (1 + f_TRZ)
    return T_c_base * (1.0 + f_TRZ);
}

double UQFFHolographicSuperconductivity::compute_U_m(double t, double t_n) {
    // U_m = (μ_j/L)(1 - exp(-γt cos(πt_n)))
    // Magnetic string energy evolves with time
    double L_uqff = compute_L_UQFF();
    if (L_uqff <= 0) return 0.0;
    
    double cos_factor = std::cos(M_PI * t_n);
    double time_factor = 1.0 - std::exp(-gamma * t * cos_factor);
    return (mu_j / L_uqff) * time_factor;
}

double UQFFHolographicSuperconductivity::compute_Delta_UQFF(double T, double t, double t_n) {
    // Δ_UQFF = Δ × exp(-U_m/(k_B×T))
    // String energy suppresses superconducting gap
    double U_m_val = compute_U_m(t, t_n);
    if (T <= 0) return 0.0;
    
    double exponent = -U_m_val / (k_B * T);
    exponent = std::max(exponent, -700.0);  // Prevent underflow
    return Delta_base * std::exp(exponent);
}

double UQFFHolographicSuperconductivity::compute_S_UQFF(double S_standard, double U_m_val) {
    // S_UQFF = S + ∫√(-g) U_m (1 + f_TRZ) d⁴x
    // String energy adds to action with time-reversal enhancement
    return S_standard + U_m_val * (1.0 + f_TRZ);
}

double UQFFHolographicSuperconductivity::compute_full_S_UQFF(double R, double F_sq, double D_psi_sq, double psi,
                                                             double B_t, double T, double t, double t_n,
                                                             double noise_level) {
    // Full S_UQFF with all components
    double L_uqff = compute_L_UQFF();
    double S_standard = compute_standard_action(R, L_uqff, F_sq, D_psi_sq, psi);
    double V_uqff = compute_V_UQFF(psi, B_t);
    double U_m_val = compute_U_m(t, t_n);
    double T_c_uqff = compute_T_c_UQFF();
    double Delta_uqff = compute_Delta_UQFF(T, t, t_n);
    
    // Combine: S = S_standard + V_UQFF + Δ_UQFF (as additional potential) + U_m term
    double S_total = S_standard + V_uqff;
    
    // Add gap contribution (normalized)
    S_total += Delta_uqff / (k_B * T_c_uqff);
    
    // Add string energy modification
    S_total = compute_S_UQFF(S_total, U_m_val);
    
    // Apply additional mods
    for (const auto& mod : additional_mods) {
        S_total *= mod(R, F_sq, psi);
    }
    
    // Add noise if requested
    if (noise_level > 0) {
        S_total += noise_level * noise_dist(rng);
    }
    
    return S_total;
}

// ═══════════════════════════════════════════════════════════════════════════════
// OPTICAL CONDUCTIVITY (Simplified Model)
// ═══════════════════════════════════════════════════════════════════════════════

double UQFFHolographicSuperconductivity::compute_sigma_dc(double T) {
    // DC conductivity: σ(ω=0) ~ (T_c - T)^{1/2} for T < T_c
    // Returns normalized conductivity
    double T_c_uqff = compute_T_c_UQFF();
    if (T >= T_c_uqff) return 0.0;  // Normal phase, no supercurrent
    return std::sqrt((T_c_uqff - T) / T_c_uqff);
}

double UQFFHolographicSuperconductivity::compute_gap_ratio() {
    // 2Δ/(k_B T_c) ratio
    // BCS: 3.53, cuprates: typically 4-6
    double T_c_uqff = compute_T_c_UQFF();
    return 2.0 * Delta_base / (k_B * T_c_uqff);
}

// ═══════════════════════════════════════════════════════════════════════════════
// SELF-EXPANSION: ADD CUSTOM MODIFIERS
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFHolographicSuperconductivity::add_mod(std::function<double(double, double, double)> mod) {
    additional_mods.push_back(mod);
}

// ═══════════════════════════════════════════════════════════════════════════════
// SELF-UPDATE: LOAD PARAMETERS FROM FILE
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFHolographicSuperconductivity::update_from_file(const std::string& config_file) {
    std::ifstream infile(config_file);
    if (!infile.is_open()) {
        std::cerr << "Failed to open config file: " << config_file << std::endl;
        return;
    }
    
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        size_t pos = line.find('=');
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            double value = std::stod(line.substr(pos + 1));
            
            if (key == "hbar") hbar = value;
            else if (key == "c") c = value;
            else if (key == "G") G = value;
            else if (key == "k_B") k_B = value;
            else if (key == "rho_vac_UA") rho_vac_UA = value;
            else if (key == "rho_vac_SCm") rho_vac_SCm = value;
            else if (key == "f_TRZ") f_TRZ = value;
            else if (key == "B_crit") B_crit = value;
            else if (key == "mu_j") mu_j = value;
            else if (key == "gamma") gamma = value;
            else if (key == "d") d = value;
            else if (key == "m_scalar") m_scalar = value;
            else if (key == "lambda_sc") lambda_sc = value;
            else if (key == "T_c_base") T_c_base = value;
            else if (key == "Delta_base") Delta_base = value;
            else if (key == "mu_chem") mu_chem = value;
        }
    }
    infile.close();
    std::cout << "Updated parameters from " << config_file << std::endl;
    
    // Re-populate explanations
    populate_explanations();
}

// ═══════════════════════════════════════════════════════════════════════════════
// SELF-SIMULATE: SWEEP OVER TEMPERATURE
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFHolographicSuperconductivity::simulate_over_temperature(double T_start, double T_end, double dT,
                                                                  const std::string& output_file) {
    std::ofstream outfile;
    bool file_output = !output_file.empty();
    
    if (file_output) {
        outfile.open(output_file);
        outfile << "T_K,T_c_UQFF_K,Delta_UQFF_J,S_UQFF,sigma_dc,phase\n";
    }
    
    double T_c_uqff = compute_T_c_UQFF();
    
    std::cout << "\nSimulating holographic superconductivity over T range [" 
              << T_start << ", " << T_end << "] K\n";
    std::cout << "T_c,UQFF = " << T_c_uqff << " K\n\n";
    
    // Example field values for simulation
    double R = 1.0;          // Ricci scalar (normalized)
    double F_sq = 1.0;       // Maxwell field squared
    double D_psi_sq = 0.1;   // Kinetic term
    double psi = 1.0;        // Order parameter
    double B_t = 0.0;        // No external field
    double t = 1e10;         // Time (seconds, ~300 years)
    double t_n = 0.5;        // Normalized time
    
    int count = 0;
    for (double T = T_start; T <= T_end; T += dT) {
        double Delta_uqff = compute_Delta_UQFF(T, t, t_n);
        double S_uqff = compute_full_S_UQFF(R, F_sq, D_psi_sq, psi, B_t, T, t, t_n, 0.0);
        double sigma = compute_sigma_dc(T);
        std::string phase = (T < T_c_uqff) ? "SUPERCONDUCTING" : "NORMAL";
        
        count++;
        
        if (file_output) {
            outfile << std::scientific << std::setprecision(6);
            outfile << T << "," << T_c_uqff << "," << Delta_uqff << "," 
                    << S_uqff << "," << sigma << "," << phase << "\n";
        } else {
            std::cout << std::fixed << std::setprecision(2);
            std::cout << "T=" << T << " K: Δ_UQFF=" << std::scientific << Delta_uqff 
                      << " J, σ_dc=" << sigma << " [" << phase << "]\n";
        }
    }
    
    if (file_output) {
        outfile.close();
        std::cout << "Simulation output saved to " << output_file << std::endl;
    }
    
    std::cout << "\nTotal points: " << count << "\n";
}

// ═══════════════════════════════════════════════════════════════════════════════
// DISPLAY EXPLANATIONS
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFHolographicSuperconductivity::display_explanations() {
    for (const auto& exp : explanations) {
        std::cout << exp << std::endl;
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// GET DERIVATION STRING
// ═══════════════════════════════════════════════════════════════════════════════

std::string UQFFHolographicSuperconductivity::get_derivation(double psi, double B_t, double T, double t, double t_n) {
    std::ostringstream oss;
    
    double L_uqff = compute_L_UQFF();
    double T_c_uqff = compute_T_c_UQFF();
    double V_uqff = compute_V_UQFF(psi, B_t);
    double U_m_val = compute_U_m(t, t_n);
    double Delta_uqff = compute_Delta_UQFF(T, t, t_n);
    double gap_ratio = compute_gap_ratio();
    double sigma = compute_sigma_dc(T);
    
    oss << std::scientific << std::setprecision(4);
    
    oss << "\n=== UQFF Holographic Superconductivity Derivation ===\n\n";
    
    oss << "PHYSICAL CONSTANTS:\n";
    oss << "  ℏ = " << hbar << " J·s\n";
    oss << "  c = " << c << " m/s\n";
    oss << "  k_B = " << k_B << " J/K\n\n";
    
    oss << "SUPERCONDUCTOR PARAMETERS:\n";
    oss << "  d = " << d << " (CFT dimension)\n";
    oss << "  m_scalar = " << m_scalar << " kg\n";
    oss << "  λ_sc = " << lambda_sc << "\n";
    oss << "  T_c (base) = " << T_c_base << " K\n";
    oss << "  Δ (base) = " << Delta_base << " J\n\n";
    
    oss << "UQFF PARAMETERS:\n";
    oss << "  ρ_vac,[UA] = " << rho_vac_UA << " J/m³\n";
    oss << "  f_TRZ = " << f_TRZ << "\n";
    oss << "  B_crit = " << B_crit << " T\n";
    oss << "  μ_j = " << mu_j << " J·m\n";
    oss << "  γ = " << gamma << " s⁻¹\n\n";
    
    oss << "INPUTS:\n";
    oss << "  ψ = " << psi << " (order parameter)\n";
    oss << "  B_t = " << B_t << " T (applied field)\n";
    oss << "  T = " << T << " K (temperature)\n";
    oss << "  t = " << t << " s (time)\n";
    oss << "  t_n = " << t_n << " (normalized time)\n\n";
    
    oss << "Step 1: Standard Holographic Action\n";
    oss << "  S = ∫√(-g) [R + d(d-1)/L² - ¼F² - |Dψ|² - m²|ψ|²] d^{d+1}x\n";
    oss << "  For AdS₄/CFT₃ (d=3): describes 2+1 dim superconductor\n\n";
    
    oss << "Step 2: UQFF Holographic Scale\n";
    oss << "  L_UQFF = (ℏc/ρ_UA)^{1/4}\n";
    oss << "         = (" << hbar << " × " << c << " / " << rho_vac_UA << ")^{0.25}\n";
    oss << "         = " << L_uqff << " m\n\n";
    
    oss << "Step 3: Modified Potential\n";
    oss << "  V_UQFF(ψ) = m²|ψ|² + (λ/2)|ψ|⁴ × (1 - B_t/B_crit)\n";
    oss << "            = " << m_scalar*m_scalar << " × " << psi*psi << " + ";
    oss << lambda_sc/2.0 << " × " << std::pow(psi,4) << " × (1 - " << B_t << "/" << B_crit << ")\n";
    oss << "            = " << V_uqff << "\n\n";
    
    oss << "Step 4: Enhanced Critical Temperature\n";
    oss << "  T_c,UQFF = T_c × (1 + f_TRZ)\n";
    oss << "           = " << T_c_base << " × (1 + " << f_TRZ << ")\n";
    oss << "           = " << T_c_uqff << " K\n";
    oss << "  Enhancement: +" << (f_TRZ * 100) << "% from time-reversal symmetry\n\n";
    
    oss << "Step 5: Magnetic String Energy\n";
    oss << "  U_m = (μ_j/L)(1 - exp(-γt cos(πt_n)))\n";
    oss << "      = (" << mu_j << "/" << L_uqff << ")(1 - exp(-" << gamma << " × " << t << " × cos(π×" << t_n << ")))\n";
    oss << "      = " << U_m_val << " J\n\n";
    
    oss << "Step 6: String-Damped Gap\n";
    oss << "  Δ_UQFF = Δ × exp(-U_m/(k_B×T))\n";
    oss << "         = " << Delta_base << " × exp(-" << U_m_val << "/(" << k_B << "×" << T << "))\n";
    double exponent = -U_m_val / (k_B * T);
    oss << "         = " << Delta_base << " × exp(" << exponent << ")\n";
    oss << "         = " << Delta_uqff << " J\n";
    oss << "  Gap ratio 2Δ/(k_B T_c) = " << gap_ratio << " (BCS: 3.53)\n\n";
    
    oss << "════════════════════════════════════════════════════════════════════════\n";
    oss << "RESULTS:\n\n";
    
    std::string phase = (T < T_c_uqff) ? "SUPERCONDUCTING" : "NORMAL";
    oss << "  ┌─────────────────────────────────────────────────────────────────┐\n";
    oss << "  │ T_c,UQFF = " << T_c_uqff << " K                                     │\n";
    oss << "  │ Δ_UQFF = " << Delta_uqff << " J                                   │\n";
    oss << "  │ σ_dc (normalized) = " << sigma << "                              │\n";
    oss << "  │ Phase at T=" << T << " K: " << phase << "                      │\n";
    oss << "  └─────────────────────────────────────────────────────────────────┘\n\n";
    
    oss << "PHYSICAL INTERPRETATION:\n";
    oss << "  • T_c enhanced by " << (f_TRZ*100) << "% due to time-reversal symmetry\n";
    oss << "  • Gap reduced by exp(-U_m/k_B T) due to magnetic strings\n";
    oss << "  • At T=" << T << " K: " << phase << " phase\n\n";
    
    oss << "Q-SCOPE TESTABILITY:\n";
    oss << "  • Thin film T_c measurement: expect +" << (f_TRZ*100) << "% enhancement\n";
    oss << "  • THz spectroscopy: Gap at 2Δ_UQFF = " << (2*Delta_uqff) << " J\n";
    oss << "  • Coherence length should match L_UQFF ~ " << L_uqff << " m\n";
    oss << "════════════════════════════════════════════════════════════════════════\n";
    
    return oss.str();
}

// ═══════════════════════════════════════════════════════════════════════════════
// MAIN: VALIDATION TESTS
// ═══════════════════════════════════════════════════════════════════════════════

int main() {
    std::cout << "UQFF Holographic Superconductivity Module\n";
    std::cout << "==========================================\n\n";
    
    // Create test instance
    UQFFHolographicSuperconductivity holo(42);
    
    // Display explanations
    holo.display_explanations();
    std::cout << "\n";
    
    // Run validation tests
    std::cout << "=== UQFF Holographic Superconductivity Validation Tests ===\n";
    int passed = 0;
    int total = 12;
    
    // Test 1: L_UQFF > 0
    double L_uqff = holo.compute_L_UQFF();
    bool t1 = L_uqff > 0;
    std::cout << "Test 1: L_UQFF > 0 - " << (t1 ? "PASS" : "FAIL")
              << " (got " << L_uqff << " m)\n";
    if (t1) passed++;
    
    // Test 2: T_c,UQFF > T_c (enhanced by f_TRZ)
    double T_c_uqff = holo.compute_T_c_UQFF();
    double T_c_base = holo.get_T_c_base();
    bool t2 = T_c_uqff > T_c_base && T_c_uqff < 2 * T_c_base;
    std::cout << "Test 2: T_c,UQFF > T_c (enhancement) - " << (t2 ? "PASS" : "FAIL")
              << " (" << T_c_base << " K → " << T_c_uqff << " K)\n";
    if (t2) passed++;
    
    // Test 3: T_c enhancement is exactly (1 + f_TRZ)
    double expected_enhancement = 1.0 + holo.get_f_TRZ();
    double actual_enhancement = T_c_uqff / T_c_base;
    bool t3 = std::abs(actual_enhancement - expected_enhancement) < 1e-10;
    std::cout << "Test 3: T_c × (1 + f_TRZ) exactly - " << (t3 ? "PASS" : "FAIL")
              << " (ratio " << actual_enhancement << ")\n";
    if (t3) passed++;
    
    // Test 4: V_UQFF(ψ, B_t=0) contains m² and λ terms
    double V_0 = holo.compute_V_UQFF(1.0, 0.0);
    bool t4 = V_0 > 0;
    std::cout << "Test 4: V_UQFF(ψ=1, B_t=0) > 0 - " << (t4 ? "PASS" : "FAIL")
              << " (got " << V_0 << ")\n";
    if (t4) passed++;
    
    // Test 5: V_UQFF decreases with B_t (quartic suppressed)
    double V_B = holo.compute_V_UQFF(1.0, 0.5 * holo.get_B_crit());
    bool t5 = V_B <= V_0;
    std::cout << "Test 5: V_UQFF decreases with B_t - " << (t5 ? "PASS" : "FAIL")
              << " (V(0)=" << V_0 << ", V(0.5B_c)=" << V_B << ")\n";
    if (t5) passed++;
    
    // Test 6: U_m(t=0, t_n) = 0 (no evolution)
    double U_m_0 = holo.compute_U_m(0.0, 0.5);
    bool t6 = std::abs(U_m_0) < 1e-300;
    std::cout << "Test 6: U_m(t=0) = 0 - " << (t6 ? "PASS" : "FAIL")
              << " (got " << U_m_0 << ")\n";
    if (t6) passed++;
    
    // Test 7: U_m increases with time
    double U_m_1 = holo.compute_U_m(1e10, 0.5);  // 300 years
    bool t7 = U_m_1 > U_m_0;
    std::cout << "Test 7: U_m(t>0) > 0 - " << (t7 ? "PASS" : "FAIL")
              << " (got " << U_m_1 << " at t=1e10 s)\n";
    if (t7) passed++;
    
    // Test 8: Δ_UQFF(t=0) = Δ_base (no string suppression at t=0)
    double Delta_t0 = holo.compute_Delta_UQFF(100.0, 0.0, 0.5);
    double Delta_base_val = holo.get_Delta_base();
    bool t8 = std::abs(Delta_t0 - Delta_base_val) < 1e-30;
    std::cout << "Test 8: Δ_UQFF(t=0) = Δ_base - " << (t8 ? "PASS" : "FAIL")
              << " (got " << Delta_t0 << ")\n";
    if (t8) passed++;
    
    // Test 9: Δ_UQFF(t>0) < Δ_base (string suppression)
    double Delta_t1 = holo.compute_Delta_UQFF(100.0, 1e10, 0.5);
    bool t9 = Delta_t1 < Delta_base_val;
    std::cout << "Test 9: Δ_UQFF(t>0) < Δ_base (suppression) - " << (t9 ? "PASS" : "FAIL")
              << " (" << Delta_t1 << " < " << Delta_base_val << ")\n";
    if (t9) passed++;
    
    // Test 10: σ_dc > 0 for T < T_c
    double sigma_low = holo.compute_sigma_dc(50.0);  // Well below T_c
    bool t10 = sigma_low > 0;
    std::cout << "Test 10: σ_dc(T<T_c) > 0 - " << (t10 ? "PASS" : "FAIL")
              << " (σ_dc(50K) = " << sigma_low << ")\n";
    if (t10) passed++;
    
    // Test 11: σ_dc = 0 for T > T_c (normal phase)
    double sigma_high = holo.compute_sigma_dc(150.0);  // Above T_c
    bool t11 = sigma_high == 0.0;
    std::cout << "Test 11: σ_dc(T>T_c) = 0 - " << (t11 ? "PASS" : "FAIL")
              << " (σ_dc(150K) = " << sigma_high << ")\n";
    if (t11) passed++;
    
    // Test 12: Derivation contains key terms
    std::string deriv = holo.get_derivation();
    bool t12 = deriv.find("T_c,UQFF") != std::string::npos &&
               deriv.find("L_UQFF") != std::string::npos &&
               deriv.find("Δ_UQFF") != std::string::npos;
    std::cout << "Test 12: Derivation contains key terms - " << (t12 ? "PASS" : "FAIL") << "\n";
    if (t12) passed++;
    
    std::cout << "=== Results: " << passed << "/" << total << " tests passed ===\n\n";
    
    // Print full derivation
    std::cout << holo.get_derivation(1.0, 0.0, 100.0, 1e10, 0.5);
    
    // Run temperature simulation
    std::cout << "\nRunning temperature simulation...\n";
    holo.simulate_over_temperature(50.0, 150.0, 10.0, "holographic_sc_vs_T.csv");
    
    return (passed == total) ? 0 : 1;
}
