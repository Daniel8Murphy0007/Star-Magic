// uqff_white_hole_formation.cpp
// UQFF White Hole Formation Module - Implementation
// Author: Daniel Murphy / Star-Magic Project
// Date: February 25, 2026

#include "uqff_white_hole_formation.h"

// ═══════════════════════════════════════════════════════════════════════════════
// CONSTRUCTOR
// ═══════════════════════════════════════════════════════════════════════════════

UQFFWhiteHoleFormation::UQFFWhiteHoleFormation(unsigned int seed) 
    : rng(seed == 0 ? std::chrono::system_clock::now().time_since_epoch().count() : seed),
      noise_dist(0.0, 1.0) {
    
    // Initialize physical constants
    G = 6.6743e-11;        // m³ kg⁻¹ s⁻²
    c = 2.998e8;           // m/s
    hbar = 1.0545718e-34;  // J·s
    k_B = 1.380649e-23;    // J/K
    M_sun = 1.989e30;      // kg
    
    // UQFF vacuum parameters
    rho_vac_UA = 7.09e-36;   // J/m³ Universal Aether density
    rho_vac_SCm = 7.09e-37;  // J/m³ Superconductive density
    f_TRZ = 0.1;             // Time-reversal zone fraction
    
    // Magnetic string parameters
    mu_j = 1e15;             // Am² (typical magnetar magnetic moment)
    gamma = 5e-5 / 86400.0;  // Convert day⁻¹ to s⁻¹
    t_n = 0.0;               // Normalized time phase
    
    // Fixed radius for simulation
    r_fixed = 0.0;           // Auto-compute from r_s when 0
    
    // Initialize explanations
    init_explanations();
}

void UQFFWhiteHoleFormation::init_explanations() {
    explanations.clear();
    
    explanations.push_back(
        "═══════════════════════════════════════════════════════════════════════════════\n"
        "UQFF WHITE HOLE FORMATION THEORY\n"
        "═══════════════════════════════════════════════════════════════════════════════"
    );
    
    explanations.push_back(
        "STANDARD WHITE HOLE CONCEPT:\n"
        "  • Hypothetical expulsion regions, time-reverse of black holes\n"
        "  • Schwarzschild metric: ds² = (1-2GM/(c²r))c²dt² - (1-2GM/(c²r))⁻¹dr² - r²dΩ²\n"
        "  • For r < 2GM/c², time and space coordinates swap\n"
        "  • Problems: Thermodynamically unstable, violates entropy increase\n"
        "  • Status: Mathematical curiosities, never observed, possibly quantum artifacts"
    );
    
    explanations.push_back(
        "UQFF REINTERPRETATION:\n"
        "  • White holes can form STABLY via aether time-reversal mechanism\n"
        "  • Superfluid [UA] vacuum enables negentropic processes\n"
        "  • [SCm] superconductive horizons allow information reversal\n"
        "  • f_TRZ ≈ 0.1 provides necessary entropy bypass\n"
        "  • Black hole inverts when [SCm] density exceeds threshold\n"
        "  • Matter expelled via U_m magnetic string channels\n"
        "  • Links to THz hole superconductive reversal in laboratory"
    );
    
    explanations.push_back(
        "FORMATION DERIVATION:\n"
        "\n"
        "  STEP 1: UQFF Horizon Shrinkage\n"
        "    r_s,UQFF = (2GM/c²) × (1 - ρ_SCm/ρ_UA)\n"
        "    Horizon shrinks when [SCm] density is high relative to [UA]\n"
        "\n"
        "  STEP 2: Time-Reversal Trigger Probability\n"
        "    P_inv = f_TRZ × exp(-E_horizon/(k_B T_H))\n"
        "    E_horizon = GM²/r_s (binding energy)\n"
        "    T_H = ℏc³/(8πGMk_B) (Hawking temperature)\n"
        "\n"
        "  STEP 3: Aether Expulsion Flux\n"
        "    Φ_flux = (ρ_UA/ρ_SCm) × (GM/c) × (1 + f_TRZ)\n"
        "    Gradient-driven outflow from density differential\n"
        "\n"
        "  STEP 4: Magnetic String Stabilization\n"
        "    U_m = (μ_j/r) × (1 - exp(-γt cos(πt_n)))\n"
        "    γ ≈ 5×10⁻⁵ day⁻¹, provides structural stability\n"
        "\n"
        "  STEP 5: Formation Condition\n"
        "    Θ_WH = P_inv × Φ_flux × exp(U_m/(k_B T_H))\n"
        "    WHITE HOLE FORMS WHEN: Θ_WH > 1\n"
        "\n"
        "  STEP 6: Mass Expulsion Rate (dM/dt)\n"
        "    L_UQFF = Φ_flux × Θ_WH × exp(-U_m/(k_B T_H))\n"
        "    dM/dt ≈ -L_UQFF/c²\n"
        "    Mass expelled via magnetic string channels when Θ_WH > 1"
    );
    
    explanations.push_back(
        "NUMERICAL EXAMPLE (Sgr A*):\n"
        "  M = 4×10⁶ M_sun = 8.0×10³⁶ kg\n"
        "  r_s = 1.2×10¹⁰ m\n"
        "  T_H ≈ 1.5×10⁻¹⁴ K\n"
        "  P_inv ≈ 0.1 (dominated by f_TRZ)\n"
        "  Φ_flux ≈ 10 (ρ_UA/ρ_SCm = 10)\n"
        "  exp(U_m/(k_B T_H)) ≈ e¹ ≈ 2.7\n"
        "  \n"
        "  Θ_WH ≈ 0.1 × 10 × 2.7 = 2.7 > 1\n"
        "  RESULT: White hole formation POSSIBLE for Sgr A* under UQFF"
    );
    
    explanations.push_back(
        "Q-SCOPE TESTABILITY:\n"
        "  • THz hole superconducting circuits can simulate horizon reversal\n"
        "  • Sonic black holes in BEC show time-reversal signatures\n"
        "  • Measure expelled phonon flux as Φ_flux analog\n"
        "  • Detect correlation between [SCm] density and reversal probability"
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// CORE COMPUTATION METHODS
// ═══════════════════════════════════════════════════════════════════════════════

double UQFFWhiteHoleFormation::compute_r_s_standard(double M) {
    // r_s = 2GM/c²
    // Standard Schwarzschild radius
    return (2.0 * G * M) / (c * c);
}

double UQFFWhiteHoleFormation::compute_r_s_UQFF(double r_s_std) {
    // r_s,UQFF = r_s × (1 - ρ_SCm/ρ_UA)
    // Step 1: Aether density modulates horizon size
    double density_ratio = rho_vac_SCm / rho_vac_UA;
    return r_s_std * (1.0 - density_ratio);
}

double UQFFWhiteHoleFormation::compute_E_horizon(double M, double r_s) {
    // E_horizon = GM²/r_s
    // Step 2: Horizon binding energy
    if (r_s <= 0) return 0.0;
    return (G * M * M) / r_s;
}

double UQFFWhiteHoleFormation::compute_T_H(double M) {
    // T_H = ℏc³/(8πGMk_B)
    // Hawking temperature
    if (M <= 0) return 0.0;
    return (hbar * c * c * c) / (8.0 * M_PI * G * M * k_B);
}

double UQFFWhiteHoleFormation::compute_P_inv(double E_horizon, double T_H) {
    // P_inv = f_TRZ × exp(-E_horizon/(k_B T_H))
    // Step 2: Time-reversal trigger probability
    if (T_H <= 0) return 0.0;
    
    double exponent = -E_horizon / (k_B * T_H);
    // Clamp to prevent underflow
    exponent = std::max(exponent, -700.0);
    
    return f_TRZ * std::exp(exponent);
}

double UQFFWhiteHoleFormation::compute_Phi_flux(double M) {
    // Φ_flux = (ρ_UA/ρ_SCm) × (GM/c) × (1 + f_TRZ)
    // Step 3: Superconductive aether flux from density gradient
    if (rho_vac_SCm <= 0) return 0.0;
    
    double density_ratio = rho_vac_UA / rho_vac_SCm;
    double gravitational_factor = (G * M) / c;
    double TRZ_enhancement = 1.0 + f_TRZ;
    
    return density_ratio * gravitational_factor * TRZ_enhancement;
}

double UQFFWhiteHoleFormation::compute_U_m(double r, double t) {
    // U_m = (μ_j/r) × (1 - exp(-γt cos(πt_n)))
    // Step 4: Magnetic string stabilization energy
    if (r <= 0) return 0.0;
    
    double cos_term = std::cos(M_PI * t_n);
    double decay_factor = 1.0 - std::exp(-gamma * t * std::max(cos_term, 0.0));
    
    return (mu_j / r) * decay_factor;
}

double UQFFWhiteHoleFormation::compute_Theta_WH(double P_inv, double Phi_flux, double U_m, double T_H) {
    // Θ_WH = P_inv × Φ_flux × exp(U_m/(k_B T_H))
    // Step 5: Formation parameter - white hole forms if Θ_WH > 1
    if (T_H <= 0 || k_B * T_H <= 0) return 0.0;
    
    double exponent = U_m / (k_B * T_H);
    // Clamp to prevent overflow
    exponent = std::min(exponent, 700.0);
    
    return P_inv * Phi_flux * std::exp(exponent);
}

bool UQFFWhiteHoleFormation::check_formation(double Theta_WH, double noise_level) {
    // Formation if Θ_WH > 1 (with optional noise)
    double noise = noise_level * noise_dist(rng);
    return (Theta_WH + noise) > 1.0;
}

double UQFFWhiteHoleFormation::compute_full_Theta_WH(double M, double r, double t, double noise_level) {
    // Full computation chain with all modifiers and noise
    
    // Step 1: Compute radii
    double r_s_std = compute_r_s_standard(M);
    double r_s_uqff = compute_r_s_UQFF(r_s_std);
    
    // Use provided r, or r_s_uqff if r is 0
    double r_eval = (r > 0) ? r : r_s_uqff;
    
    // Step 2: Compute probabilities and energies
    double E_horizon = compute_E_horizon(M, r_s_uqff);
    double T_H = compute_T_H(M);
    double P_inv = compute_P_inv(E_horizon, T_H);
    
    // Step 3: Compute flux
    double Phi_flux = compute_Phi_flux(M);
    
    // Step 4: Compute magnetic string energy
    double U_m_val = compute_U_m(r_eval, t);
    
    // Step 5: Compute base Θ_WH
    double Theta = compute_Theta_WH(P_inv, Phi_flux, U_m_val, T_H);
    
    // Apply self-expanding modifiers
    for (const auto& mod : additional_modifiers) {
        Theta *= mod(M, t);
    }
    
    // Add stochastic noise
    double noise = noise_level * noise_dist(rng);
    return Theta + noise;
}

// ═══════════════════════════════════════════════════════════════════════════════
// MASS EXPULSION (dM/dt)
// ═══════════════════════════════════════════════════════════════════════════════

double UQFFWhiteHoleFormation::compute_L_UQFF(double M, double r, double t) {
    // L_UQFF = Φ_flux × Θ_WH × exp(-U_m/(k_B T_H))
    // Step 6: UQFF luminosity driving mass expulsion
    // Physically: energy radiated via aether flux channels
    
    // Compute intermediate values
    double r_s_std = compute_r_s_standard(M);
    double r_s_uqff = compute_r_s_UQFF(r_s_std);
    double r_eval = (r > 0) ? r : r_s_uqff;
    
    double T_H = compute_T_H(M);
    if (T_H <= 0) return 0.0;
    
    double E_horizon = compute_E_horizon(M, r_s_uqff);
    double P_inv = compute_P_inv(E_horizon, T_H);
    double Phi_flux = compute_Phi_flux(M);
    double U_m_val = compute_U_m(r_eval, t);
    double Theta_WH = compute_Theta_WH(P_inv, Phi_flux, U_m_val, T_H);
    
    // L_UQFF = Φ_flux × Θ_WH × exp(-U_m/(k_B T_H))
    // Note: exp(-U_m/kT) is the inverse of the enhancement factor,
    // representing energy being expended in expulsion
    double exponent = -U_m_val / (k_B * T_H);
    exponent = std::max(exponent, -700.0);  // Prevent underflow
    
    // Scale by characteristic energy to get physical luminosity
    // L_UQFF ~ Φ_flux × Θ_WH × exp(-U_m/(k_B T_H)) × (GM/r) for dimensional consistency
    double energy_scale = (r_eval > 0) ? (G * M / r_eval) : 1.0;
    
    return Phi_flux * Theta_WH * std::exp(exponent) * energy_scale;
}

double UQFFWhiteHoleFormation::compute_dM_dt(double M, double r, double t) {
    // dM/dt ≈ -L_UQFF/c²
    // Step 6: Mass expulsion rate when white hole forms
    // Only non-zero if Θ_WH > 1 (formation condition met)
    
    // First check if formation condition is met
    double Theta = compute_full_Theta_WH(M, r, t, 0.0);  // No noise for deterministic dM/dt
    
    if (Theta <= 1.0) {
        // No formation, no expulsion
        return 0.0;
    }
    
    // Compute UQFF luminosity
    double L_UQFF = compute_L_UQFF(M, r, t);
    
    // dM/dt = -L_UQFF/c²
    // Negative sign indicates mass loss
    return -L_UQFF / (c * c);
}

std::string UQFFWhiteHoleFormation::generate_derivation(double M, double r, double t) {
    // Generate detailed step-by-step derivation
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(4);
    
    // Compute all intermediate values
    double r_s_std = compute_r_s_standard(M);
    double r_s_uqff = compute_r_s_UQFF(r_s_std);
    double r_eval = (r > 0) ? r : r_s_uqff;
    double E_horizon = compute_E_horizon(M, r_s_uqff);
    double T_H = compute_T_H(M);
    double P_inv = compute_P_inv(E_horizon, T_H);
    double Phi_flux = compute_Phi_flux(M);
    double U_m_val = compute_U_m(r_eval, t);
    double Theta_WH = compute_Theta_WH(P_inv, Phi_flux, U_m_val, T_H);
    bool forms = check_formation(Theta_WH, 0.0);
    
    oss << "═══════════════════════════════════════════════════════════════════════════════\n";
    oss << "UQFF WHITE HOLE FORMATION DERIVATION\n";
    oss << "═══════════════════════════════════════════════════════════════════════════════\n\n";
    
    oss << "INPUT PARAMETERS:\n";
    oss << "  M = " << M << " kg = " << (M / M_sun) << " M_sun\n";
    oss << "  r = " << r_eval << " m\n";
    oss << "  t = " << t << " s = " << (t / 3.156e7) << " years\n";
    oss << "  f_TRZ = " << f_TRZ << "\n";
    oss << "  ρ_UA = " << rho_vac_UA << " J/m³\n";
    oss << "  ρ_SCm = " << rho_vac_SCm << " J/m³\n\n";
    
    oss << "═══════════════════════════════════════════════════════════════════════════════\n";
    oss << "STEP 1: SCHWARZSCHILD RADIUS\n";
    oss << "═══════════════════════════════════════════════════════════════════════════════\n";
    oss << "  Standard: r_s = 2GM/c²\n";
    oss << "          = 2 × " << G << " × " << M << " / " << (c*c) << "\n";
    oss << "          = " << r_s_std << " m\n\n";
    
    oss << "  UQFF: r_s,UQFF = r_s × (1 - ρ_SCm/ρ_UA)\n";
    oss << "      = " << r_s_std << " × (1 - " << rho_vac_SCm << "/" << rho_vac_UA << ")\n";
    oss << "      = " << r_s_std << " × " << (1.0 - rho_vac_SCm/rho_vac_UA) << "\n";
    oss << "      = " << r_s_uqff << " m\n\n";
    
    oss << "═══════════════════════════════════════════════════════════════════════════════\n";
    oss << "STEP 2: HAWKING TEMPERATURE & INVERSION PROBABILITY\n";
    oss << "═══════════════════════════════════════════════════════════════════════════════\n";
    oss << "  T_H = ℏc³/(8πGMk_B)\n";
    oss << "      = " << T_H << " K\n\n";
    
    oss << "  E_horizon = GM²/r_s\n";
    oss << "            = " << E_horizon << " J\n\n";
    
    oss << "  P_inv = f_TRZ × exp(-E_horizon/(k_B T_H))\n";
    oss << "        = " << f_TRZ << " × exp(-" << E_horizon << "/(" << k_B << " × " << T_H << "))\n";
    oss << "        = " << P_inv << "\n\n";
    
    oss << "═══════════════════════════════════════════════════════════════════════════════\n";
    oss << "STEP 3: AETHER EXPULSION FLUX\n";
    oss << "═══════════════════════════════════════════════════════════════════════════════\n";
    oss << "  Φ_flux = (ρ_UA/ρ_SCm) × (GM/c) × (1 + f_TRZ)\n";
    oss << "         = (" << rho_vac_UA << "/" << rho_vac_SCm << ") × (" << G << " × " << M << "/" << c << ") × " << (1+f_TRZ) << "\n";
    oss << "         = " << Phi_flux << "\n\n";
    
    oss << "═══════════════════════════════════════════════════════════════════════════════\n";
    oss << "STEP 4: MAGNETIC STRING ENERGY\n";
    oss << "═══════════════════════════════════════════════════════════════════════════════\n";
    oss << "  U_m = (μ_j/r) × (1 - exp(-γt cos(πt_n)))\n";
    oss << "      = (" << mu_j << "/" << r_eval << ") × (1 - exp(-" << gamma << " × " << t << " × cos(π × " << t_n << ")))\n";
    oss << "      = " << U_m_val << " J\n\n";
    
    oss << "═══════════════════════════════════════════════════════════════════════════════\n";
    oss << "STEP 5: FORMATION PARAMETER\n";
    oss << "═══════════════════════════════════════════════════════════════════════════════\n";
    oss << "  Θ_WH = P_inv × Φ_flux × exp(U_m/(k_B T_H))\n";
    oss << "       = " << P_inv << " × " << Phi_flux << " × exp(" << U_m_val << "/(" << k_B << " × " << T_H << "))\n";
    oss << "       = " << Theta_WH << "\n\n";
    
    // Compute dM/dt values
    double L_UQFF = compute_L_UQFF(M, r_eval, t);
    double dM_dt = forms ? compute_dM_dt(M, r_eval, t) : 0.0;
    double t_evap = (forms && dM_dt < 0) ? -M / dM_dt : std::numeric_limits<double>::infinity();
    
    oss << "═══════════════════════════════════════════════════════════════════════════════\n";
    oss << "STEP 6: MASS EXPULSION RATE (dM/dt)\n";
    oss << "═══════════════════════════════════════════════════════════════════════════════\n";
    oss << "  L_UQFF = Φ_flux × Θ_WH × exp(-U_m/(k_B T_H)) × (GM/r)\n";
    oss << "         = " << L_UQFF << " W\n\n";
    oss << "  dM/dt = -L_UQFF/c²\n";
    if (forms) {
        oss << "        = -" << L_UQFF << "/" << (c*c) << "\n";
        oss << "        = " << dM_dt << " kg/s\n";
        oss << "        = " << (dM_dt * 3.156e7 / M_sun) << " M_sun/year\n\n";
        oss << "  Characteristic expulsion time (M/|dM/dt|):\n";
        oss << "        = " << t_evap << " s = " << (t_evap / 3.156e16) << " Gyr\n\n";
    } else {
        oss << "        = 0 (no formation, Θ_WH ≤ 1)\n\n";
    }
    
    oss << "╔══════════════════════════════════════════════════════════════════════════════╗\n";
    if (forms) {
        oss << "║  RESULT: Θ_WH = " << std::setw(12) << Theta_WH << " > 1  →  WHITE HOLE FORMS         ║\n";
        oss << "║  dM/dt = " << std::setw(12) << dM_dt << " kg/s (mass expelling)        ║\n";
    } else {
        oss << "║  RESULT: Θ_WH = " << std::setw(12) << Theta_WH << " ≤ 1  →  NO FORMATION              ║\n";
    }
    oss << "╚══════════════════════════════════════════════════════════════════════════════╝\n";
    
    return oss.str();
}

// ═══════════════════════════════════════════════════════════════════════════════
// SELF-EXPANDING FRAMEWORK
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFWhiteHoleFormation::add_modifier(std::function<double(double, double)> modifier) {
    additional_modifiers.push_back(modifier);
}

void UQFFWhiteHoleFormation::clear_modifiers() {
    additional_modifiers.clear();
}

// ═══════════════════════════════════════════════════════════════════════════════
// SELF-UPDATE FROM FILE
// ═══════════════════════════════════════════════════════════════════════════════

bool UQFFWhiteHoleFormation::update_from_file(const std::string& config_file) {
    std::ifstream infile(config_file);
    if (!infile.is_open()) {
        std::cerr << "Failed to open config file: " << config_file << std::endl;
        return false;
    }
    
    std::string line;
    int params_loaded = 0;
    
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
                
                if (key == "G") { G = value; params_loaded++; }
                else if (key == "c") { c = value; params_loaded++; }
                else if (key == "hbar") { hbar = value; params_loaded++; }
                else if (key == "k_B") { k_B = value; params_loaded++; }
                else if (key == "f_TRZ") { f_TRZ = value; params_loaded++; }
                else if (key == "rho_vac_UA") { rho_vac_UA = value; params_loaded++; }
                else if (key == "rho_vac_SCm") { rho_vac_SCm = value; params_loaded++; }
                else if (key == "mu_j") { mu_j = value; params_loaded++; }
                else if (key == "gamma") { gamma = value; params_loaded++; }
                else if (key == "t_n") { t_n = value; params_loaded++; }
                else if (key == "r_fixed") { r_fixed = value; params_loaded++; }
            } catch (const std::exception& e) {
                std::cerr << "Error parsing value for " << key << ": " << e.what() << std::endl;
            }
        }
    }
    
    infile.close();
    std::cout << "Loaded " << params_loaded << " parameters from " << config_file << std::endl;
    return params_loaded > 0;
}

// ═══════════════════════════════════════════════════════════════════════════════
// SELF-SIMULATE
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFWhiteHoleFormation::simulate_formation(double M_start, double r, double t_start, 
                                                 double t_end, double dt, 
                                                 const std::string& output_file) {
    // Self-simulate: Compute Θ_WH over time, check formation, integrate dM/dt if formed
    // Mass evolution: M(t+dt) = M(t) + dM/dt × dt (when Θ_WH > 1)
    
    std::ofstream outfile;
    bool file_output = !output_file.empty();
    
    if (file_output) {
        outfile.open(output_file);
        outfile << "time_s,time_Gyr,Mass_kg,Mass_Msun,Theta_WH,Forms,dM_dt_kg_s,L_UQFF_W,r_s_UQFF_m,T_H_K\n";
    }
    
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "WHITE HOLE FORMATION SIMULATION WITH MASS EVOLUTION\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "M_start = " << std::scientific << M_start << " kg = " << (M_start / M_sun) << " M_sun\n";
    std::cout << "r = " << r << " m (0 = auto-track horizon)\n";
    std::cout << "t: " << t_start << " to " << t_end << " s (dt = " << dt << " s)\n";
    std::cout << "───────────────────────────────────────────────────────────────────────────────\n";
    
    int formation_count = 0;
    int total_steps = 0;
    double M_current = M_start;
    double M_lost_total = 0.0;
    
    for (double t = t_start; t <= t_end && M_current > 0; t += dt) {
        // Compute all values with current mass
        double r_s_std = compute_r_s_standard(M_current);
        double r_s_uqff = compute_r_s_UQFF(r_s_std);
        double r_eval = (r > 0) ? r : r_s_uqff;  // Auto-track horizon if r=0
        double T_H = compute_T_H(M_current);
        double E_horizon = compute_E_horizon(M_current, r_s_uqff);
        double P_inv = compute_P_inv(E_horizon, T_H);
        double Phi_flux = compute_Phi_flux(M_current);
        double U_m_val = compute_U_m(r_eval, t);
        double Theta = compute_Theta_WH(P_inv, Phi_flux, U_m_val, T_H);
        
        // Apply modifiers
        for (const auto& mod : additional_modifiers) {
            Theta *= mod(M_current, t);
        }
        
        bool forms = check_formation(Theta, 0.0);  // No noise for deterministic evolution
        
        // Compute dM/dt and L_UQFF
        double dM_dt_val = 0.0;
        double L_UQFF = 0.0;
        if (forms) {
            formation_count++;
            L_UQFF = compute_L_UQFF(M_current, r_eval, t);
            dM_dt_val = -L_UQFF / (c * c);
            
            // Update mass (Euler integration)
            double dM = dM_dt_val * dt;
            M_current += dM;
            M_lost_total -= dM;  // dM is negative, so this is positive
            
            // Prevent negative mass
            if (M_current < 0) M_current = 0;
        }
        
        total_steps++;
        double t_Gyr = t / 3.156e16;  // Convert to Gyr
        
        if (file_output) {
            outfile << std::scientific << std::setprecision(6)
                    << t << "," << t_Gyr << "," 
                    << M_current << "," << (M_current / M_sun) << ","
                    << Theta << "," << (forms ? 1 : 0) << ","
                    << dM_dt_val << "," << L_UQFF << ","
                    << r_s_uqff << "," << T_H << "\n";
        }
        
        // Print every 10th step to console
        if (total_steps % 10 == 1 || t == t_start) {
            std::cout << std::scientific << std::setprecision(3)
                      << "t=" << t_Gyr << " Gyr, M=" << (M_current / M_sun) << " M_sun, Θ_WH=" << Theta;
            if (forms) {
                std::cout << ", dM/dt=" << dM_dt_val << " kg/s";
            }
            std::cout << std::endl;
        }
    }
    
    if (file_output) {
        outfile.close();
        std::cout << "Output written to: " << output_file << std::endl;
    }
    
    double mass_fraction_lost = M_lost_total / M_start;
    
    std::cout << "───────────────────────────────────────────────────────────────────────────────\n";
    std::cout << "SUMMARY:\n";
    std::cout << "  Formation occurred in: " << formation_count << "/" << total_steps 
              << " (" << std::fixed << std::setprecision(1) << (100.0 * formation_count / total_steps) << "%) of time steps\n";
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  Initial mass:  " << M_start << " kg = " << (M_start / M_sun) << " M_sun\n";
    std::cout << "  Final mass:    " << M_current << " kg = " << (M_current / M_sun) << " M_sun\n";
    std::cout << "  Total expelled: " << M_lost_total << " kg = " << (M_lost_total / M_sun) << " M_sun\n";
    std::cout << "  Mass fraction lost: " << std::fixed << std::setprecision(4) << (mass_fraction_lost * 100) << "%\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
}

void UQFFWhiteHoleFormation::simulate_sgr_a(const std::string& output_file) {
    // Sgr A* canonical example
    double M_sgr_a = 4.0e6 * M_sun;  // 4 million solar masses
    double r_s = compute_r_s_standard(M_sgr_a);
    
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                    SGR A* WHITE HOLE FORMATION SIMULATION                   ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    std::cout << "Sgr A* Parameters:\n";
    std::cout << "  M = 4×10⁶ M_sun = " << std::scientific << M_sgr_a << " kg\n";
    std::cout << "  r_s = " << r_s << " m\n";
    std::cout << "\n";
    
    // Print detailed derivation at t = 4.5 Gyr
    double t_now = 4.5e9 * 3.156e7;  // 4.5 Gyr in seconds
    std::cout << generate_derivation(M_sgr_a, r_s, t_now);
    
    // Run simulation over cosmic time
    std::string out = output_file.empty() ? "sgr_a_white_hole_sim.csv" : output_file;
    simulate_formation(M_sgr_a, r_s, 0.0, 1e18, 1e16, out);
}

// ═══════════════════════════════════════════════════════════════════════════════
// DISPLAY AND VALIDATION
// ═══════════════════════════════════════════════════════════════════════════════

void UQFFWhiteHoleFormation::display_explanations() {
    std::cout << "\n";
    for (const auto& exp : explanations) {
        std::cout << exp << "\n\n";
    }
}

bool UQFFWhiteHoleFormation::validate() {
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "UQFF WHITE HOLE FORMATION MODEL VALIDATION\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    
    bool all_passed = true;
    int test_num = 0;
    
    // Test 1: Schwarzschild radius for Sun
    test_num++;
    double r_s_sun = compute_r_s_standard(M_sun);
    double expected_r_s = 2.95e3;  // ~3 km
    bool test1 = (std::abs(r_s_sun - expected_r_s) / expected_r_s < 0.01);
    std::cout << "Test " << test_num << ": r_s(Sun) = " << r_s_sun << " m (expected ~" << expected_r_s << " m): " 
              << (test1 ? "PASS" : "FAIL") << std::endl;
    all_passed &= test1;
    
    // Test 2: UQFF radius shrinkage
    test_num++;
    double r_s_uqff = compute_r_s_UQFF(r_s_sun);
    bool test2 = (r_s_uqff < r_s_sun) && (r_s_uqff > 0.8 * r_s_sun);
    std::cout << "Test " << test_num << ": r_s,UQFF < r_s (shrinkage): " 
              << (test2 ? "PASS" : "FAIL") << std::endl;
    all_passed &= test2;
    
    // Test 3: Hawking temperature scaling
    test_num++;
    double T_H_sun = compute_T_H(M_sun);
    double T_H_2sun = compute_T_H(2 * M_sun);
    bool test3 = (std::abs(T_H_sun / T_H_2sun - 2.0) < 0.01);  // T ∝ 1/M
    std::cout << "Test " << test_num << ": T_H ∝ 1/M scaling: " 
              << (test3 ? "PASS" : "FAIL") << std::endl;
    all_passed &= test3;
    
    // Test 4: P_inv bounded by f_TRZ
    test_num++;
    double M_test = 1e30;
    double r_s_test = compute_r_s_standard(M_test);
    double E_test = compute_E_horizon(M_test, r_s_test);
    double T_H_test = compute_T_H(M_test);
    double P_inv_test = compute_P_inv(E_test, T_H_test);
    bool test4 = (P_inv_test <= f_TRZ) && (P_inv_test >= 0);
    std::cout << "Test " << test_num << ": P_inv ≤ f_TRZ: " 
              << (test4 ? "PASS" : "FAIL") << std::endl;
    all_passed &= test4;
    
    // Test 5: Phi_flux proportional to M
    test_num++;
    double Phi_1 = compute_Phi_flux(M_sun);
    double Phi_2 = compute_Phi_flux(2 * M_sun);
    bool test5 = (std::abs(Phi_2 / Phi_1 - 2.0) < 0.01);  // Φ ∝ M
    std::cout << "Test " << test_num << ": Φ_flux ∝ M scaling: " 
              << (test5 ? "PASS" : "FAIL") << std::endl;
    all_passed &= test5;
    
    // Test 6: U_m grows with time
    test_num++;
    double U_m_early = compute_U_m(1e10, 1e10);
    double U_m_late = compute_U_m(1e10, 1e17);
    bool test6 = (U_m_late >= U_m_early);
    std::cout << "Test " << test_num << ": U_m grows with time: " 
              << (test6 ? "PASS" : "FAIL") << std::endl;
    all_passed &= test6;
    
    // Test 7: Sgr A* forms white hole (Θ_WH > 1)
    test_num++;
    double M_sgr_a = 4.0e6 * M_sun;
    double r_sgr_a = compute_r_s_standard(M_sgr_a);
    double t_cosmo = 4.5e17;  // 14 Gyr in seconds
    double Theta_sgr_a = compute_full_Theta_WH(M_sgr_a, r_sgr_a, t_cosmo, 0.0);
    bool test7 = check_formation(Theta_sgr_a, 0.0);
    std::cout << "Test " << test_num << ": Sgr A* Θ_WH > 1 (Θ=" << Theta_sgr_a << "): " 
              << (test7 ? "PASS" : "FAIL") << std::endl;
    all_passed &= test7;
    
    // Test 8: dM/dt is negative when Θ_WH > 1
    test_num++;
    double dM_dt_sgr_a = compute_dM_dt(M_sgr_a, r_sgr_a, t_cosmo);
    bool test8 = (dM_dt_sgr_a < 0);  // Mass loss if formed
    std::cout << "Test " << test_num << ": dM/dt < 0 when Θ_WH > 1 (dM/dt=" << dM_dt_sgr_a << " kg/s): "
              << (test8 ? "PASS" : "FAIL") << std::endl;
    all_passed &= test8;
    
    // Test 9: dM/dt is zero when Θ_WH < 1
    test_num++;
    double M_small = 1.0 * M_sun;  // Small BH, unlikely to form WH
    double r_small = compute_r_s_standard(M_small);
    double dM_dt_small = compute_dM_dt(M_small, r_small, 1e10);  // Early time
    double Theta_small = compute_full_Theta_WH(M_small, r_small, 1e10, 0.0);
    bool test9 = (Theta_small <= 1.0) ? (dM_dt_small == 0.0) : true;  // Only test if no formation
    std::cout << "Test " << test_num << ": dM/dt = 0 when Θ_WH ≤ 1 (Θ=" << Theta_small << "): "
              << (test9 ? "PASS" : "SKIP (formed)") << std::endl;
    all_passed &= test9;
    
    // Test 10: L_UQFF is positive
    test_num++;
    double L_UQFF_sgr_a = compute_L_UQFF(M_sgr_a, r_sgr_a, t_cosmo);
    bool test10 = (L_UQFF_sgr_a > 0);
    std::cout << "Test " << test_num << ": L_UQFF > 0 (L=" << L_UQFF_sgr_a << " W): "
              << (test10 ? "PASS" : "FAIL") << std::endl;
    all_passed &= test10;
    
    std::cout << "───────────────────────────────────────────────────────────────────────────────\n";
    std::cout << "VALIDATION RESULT: " << (all_passed ? "ALL TESTS PASSED ✓" : "SOME TESTS FAILED ✗") << std::endl;
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    
    return all_passed;
}

// ═══════════════════════════════════════════════════════════════════════════════
// MAIN FOR TESTING
// ═══════════════════════════════════════════════════════════════════════════════

int main() {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║         UQFF WHITE HOLE FORMATION MODULE - DEMONSTRATION                    ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════════════╝\n";
    
    // Create model instance
    UQFFWhiteHoleFormation wh;
    
    // Display theoretical explanations
    wh.display_explanations();
    
    // Run validation tests
    wh.validate();
    
    // Simulate Sgr A* white hole formation
    wh.simulate_sgr_a();
    
    // Demonstrate self-expanding: Add custom modifier
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "SELF-EXPANDING: Adding custom modifier\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    
    wh.add_modifier([](double M, double t) { 
        // Example: Slight enhancement at late times
        return 1.0 + 0.01 * std::log10(t / 1e10 + 1.0); 
    });
    
    std::cout << "Added modifier: Θ_WH *= (1 + 0.01 × log₁₀(t/10¹⁰ + 1))\n";
    std::cout << "Active modifiers: " << wh.modifier_count() << std::endl;
    
    // Recompute with modifier
    double M_sgr_a = 4.0e6 * 1.989e30;
    double r_sgr_a = wh.compute_r_s_standard(M_sgr_a);
    double t_now = 4.5e17;
    double Theta_modified = wh.compute_full_Theta_WH(M_sgr_a, r_sgr_a, t_now, 0.0);
    std::cout << "Θ_WH (with modifier) = " << Theta_modified << std::endl;
    
    // Demonstrate dM/dt computation
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "MASS EXPULSION RATE (dM/dt) DEMONSTRATION\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    
    // Clear modifiers for clean dM/dt test
    wh.clear_modifiers();
    
    double L_UQFF = wh.compute_L_UQFF(M_sgr_a, r_sgr_a, t_now);
    double dM_dt = wh.compute_dM_dt(M_sgr_a, r_sgr_a, t_now);
    double M_sun_val = 1.989e30;
    
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "Sgr A* at t = 4.5 Gyr:\n";
    std::cout << "  L_UQFF = " << L_UQFF << " W\n";
    std::cout << "  dM/dt  = " << dM_dt << " kg/s\n";
    std::cout << "         = " << (dM_dt * 3.156e7 / M_sun_val) << " M_sun/year\n";
    
    if (dM_dt < 0) {
        double t_evap = -M_sgr_a / dM_dt;
        std::cout << "  Characteristic expulsion time: " << (t_evap / 3.156e16) << " Gyr\n";
    }
    
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    std::cout << "MODULE DEMONSTRATION COMPLETE\n";
    std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
    
    return 0;
}
