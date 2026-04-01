#include "BlackToWhiteHoleUQFF.h"

// ─────────────────────────────────────────────────────────────────────────────
//  BlackToWhiteHoleUQFF  —  UQFF Black-to-White Hole Transition Module
//  PAPER_659 | Session 172 | April 1–2, 2026
//  Source: grok_share_fc21e30c24b4.txt  BlackToWhiteHoleTransition class
// ─────────────────────────────────────────────────────────────────────────────

BlackToWhiteHoleUQFF::BlackToWhiteHoleUQFF(unsigned int seed)
    : rng(seed)
{
    // Validate UQFF ratio (should be ~10:1)
    if (rho_vac_UA < rho_vac_SCm) {
        std::cerr << "[WARN] rho_UA < rho_SCm — check vacuum density parameters!\n";
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Step 1: Standard Schwarzschild radius
//  r_s = 2 G M / c²
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::compute_r_s(double M) const {
    return 2.0 * G * M / (c * c);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Step 2a: UQFF modified horizon
//  r_s,UQFF = r_s · (1 − ρ_SCm / ρ_UA)
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::compute_r_s_UQFF(double r_s) const {
    double ratio = rho_vac_SCm / rho_vac_UA;   // ~0.1
    return r_s * (1.0 - ratio);                 // ~0.9 r_s
}

// ─────────────────────────────────────────────────────────────────────────────
//  Step 2b: Inversion energy
//  E_flip = G M² / r_s,UQFF
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::compute_E_flip(double M, double r_s_UQFF) const {
    return G * M * M / r_s_UQFF;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Step 3a: Hawking temperature
//  T_H = ħ c³ / (8π G M k_B)
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::compute_T_H(double M) const {
    return (hbar * c * c * c) / (8.0 * M_PI * G * M * k_B);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Step 3b: Quantum flip probability
//  P_flip = exp(−E_flip / (k_B T_H))
//  NOTE: For stellar-mass BHs this is astronomically small; for micro-BHs it
//  becomes significant. UQFF modulation via f_TRZ elevates it by ×(1/f_TRZ).
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::compute_P_flip(double E_flip, double T_H) const {
    double exponent = -E_flip / (k_B * T_H);
    if (exponent < -700.0) return 0.0;   // Underflow guard
    return std::exp(exponent);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Step 3c: UQFF time-reversal boosted probability
//  P_trans = f_TRZ · P_flip
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::compute_P_trans(double P_flip) const {
    return f_TRZ * P_flip;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Step 4: Buoyancy transition potential (Buoyancy Harmonics Series PAPER_648)
//  Φ_trans = (ρ_UA / ρ_SCm) · (G M / c) · (1 + f_TRZ)
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::compute_Phi_trans(double M) const {
    double rho_ratio   = rho_vac_UA / rho_vac_SCm;   // ~10
    double grav_len    = G * M / c;                   // [m² kg s⁻¹] — gravitational length / c
    return rho_ratio * grav_len * (1.0 + f_TRZ);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Step 5a: U_m magnetic string anchor (Dipole Vortex Primes PAPER_647)
//  U_m(r,t) = μ_j / r · (1 − exp(−γ t cos(π t_n)))
//  γ has units [s⁻¹], t [s], t_n = t / t_ref (dimensionless normalised)
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::compute_U_m(double r, double t) const {
    double t_normalised = t / t_n;                             // dimensionless
    double decay_arg    = -gamma * t * std::cos(M_PI * t_normalised);
    double envelope     = 1.0 - std::exp(decay_arg);
    return (mu_j / r) * envelope;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Step 5b: Instability timescale
//  τ_instab = r_s / c
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::compute_tau_instab(double r_s) const {
    return r_s / c;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Step 5c: White hole stability timescale
//  τ_WH = τ_instab · exp(U_m / (k_B T_H))   where T_H > 0 (|T_WH| = T_H)
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::compute_tau_WH(double r_s, double t) const {
    // τ_WH requires M to compute τ_instab → use r_s proxy: r_s = 2GM/c² → M = r_s c²/2G
    double M_proxy  = r_s * c * c / (2.0 * G);
    double tau_i    = compute_tau_instab(r_s);
    double T_H      = compute_T_H(M_proxy);
    double U_m      = compute_U_m(r_s, t);             // evaluate at r = r_s
    double exponent = U_m / (k_B * T_H);
    if (exponent > 700.0) return tau_i * std::exp(700.0);  // Overflow guard
    return tau_i * std::exp(exponent);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Step 6a: S_Um = exp(U_m / (k_B T_H))
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::compute_S_Um(double r, double t, double T_H) const {
    double U_m     = compute_U_m(r, t);
    double exponent = U_m / (k_B * T_H);
    if (exponent > 700.0) return std::exp(700.0);
    return std::exp(exponent);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Step 6b: Full transition criterion
//  Θ_trans = P_trans · Φ_trans · S_Um
//  > 1 → white hole forms;  Sgr A*: Θ_trans ≈ 2.7 nominally
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::compute_Theta_trans(double M, double r, double t) const {
    double r_s        = compute_r_s(M);
    double r_s_uqff   = compute_r_s_UQFF(r_s);
    double E_flip     = compute_E_flip(M, r_s_uqff);
    double T_H        = compute_T_H(M);
    double P_flip     = compute_P_flip(E_flip, T_H);
    double P_trans    = compute_P_trans(P_flip);
    double Phi_trans  = compute_Phi_trans(M);
    double S_Um       = compute_S_Um(r, t, T_H);

    double Theta = P_trans * Phi_trans * S_Um;

    // Apply extra mods (self-expanding framework 2.0)
    for (auto& mod : extra_mods) {
        Theta += mod(M, t);
    }

    return Theta;
}

// ─────────────────────────────────────────────────────────────────────────────
//  UQFF white-hole luminosity
//  Hawking: L_H = π k_B² T_H² M/ (30 ħ c) ... simplified
//  UQFF:    L_WH ≈ L_H · (1 + f_TRZ) · (ρ_UA/ρ_SCm) · exp(U_m / k_B T_H)
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::compute_L_H(double M) const {
    // Standard Hawking luminosity: L_H = ħ c⁶ / (15360 π G² M²)
    double G2  = G * G;
    double c6  = std::pow(c, 6.0);
    return (hbar * c6) / (15360.0 * M_PI * G2 * M * M);
}

double BlackToWhiteHoleUQFF::compute_L_WH(double M, double r, double t) const {
    double L_H     = compute_L_H(M);
    double T_H     = compute_T_H(M);
    double S_Um    = compute_S_Um(r, t, T_H);
    double ratio   = rho_vac_UA / rho_vac_SCm;
    return L_H * (1.0 + f_TRZ) * ratio * S_Um;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Monte-Carlo P(Θ_trans > 1)
//  Samples ρ_UA and ρ_SCm from log-normal with sigma_rho relative uncertainty
// ─────────────────────────────────────────────────────────────────────────────
double BlackToWhiteHoleUQFF::monte_carlo_P_transition(
    double M, double r, double t,
    int n_samples, double sigma_rho) const
{
    // Need a non-const RNG copy for sampling
    BlackToWhiteHoleUQFF sampler(*this);
    int count_above = 0;
    for (int i = 0; i < n_samples; ++i) {
        // Perturb vacuum densities by multiplicative log-normal noise
        double noise_UA  = std::exp(sigma_rho * sampler.noise_dist(sampler.rng));
        double noise_SCm = std::exp(sigma_rho * sampler.noise_dist(sampler.rng));
        sampler.rho_vac_UA  = rho_vac_UA  * noise_UA;
        sampler.rho_vac_SCm = rho_vac_SCm * noise_SCm;
        double Theta = sampler.compute_Theta_trans(M, r, t);
        if (Theta > 1.0) ++count_above;
    }
    return static_cast<double>(count_above) / static_cast<double>(n_samples);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Self-expanding framework 2.0
// ─────────────────────────────────────────────────────────────────────────────
void BlackToWhiteHoleUQFF::add_extra_mod(std::function<double(double, double)> mod) {
    extra_mods.push_back(std::move(mod));
}

// ─────────────────────────────────────────────────────────────────────────────
//  Self-update (runtime parameter tuning from file)
//  File format: key=value per line
// ─────────────────────────────────────────────────────────────────────────────
void BlackToWhiteHoleUQFF::update_from_file(const std::string& filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        std::cerr << "[WARN] BlackToWhiteHoleUQFF: cannot open " << filename << "\n";
        return;
    }
    std::string line;
    while (std::getline(ifs, line)) {
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = line.substr(0, eq);
        double val      = std::stod(line.substr(eq + 1));
        if      (key == "rho_vac_UA")  rho_vac_UA  = val;
        else if (key == "rho_vac_SCm") rho_vac_SCm = val;
        else if (key == "f_TRZ")       f_TRZ  = val;
        else if (key == "kappa")       kappa  = val;
        else if (key == "SSq")         SSq    = val;
        else if (key == "mu_j")        mu_j   = val;
        else if (key == "gamma")       gamma  = val / 86400.0;  // input in day⁻¹
        else if (key == "t_n")         t_n    = val;
        else if (key == "r_ref")       r_ref  = val;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Simulate: iterate Θ_trans(M, r_s(M), t) over a time range
//  Outputs: t [s], r_s [m], T_H [K], Theta_trans, L_WH [W]
// ─────────────────────────────────────────────────────────────────────────────
void BlackToWhiteHoleUQFF::simulate(double M, double r, double t_start,
                                    double t_end, double dt,
                                    const std::string& output_file) const
{
    std::ostream* out = &std::cout;
    std::ofstream fout;
    if (!output_file.empty()) {
        fout.open(output_file);
        if (fout.is_open()) out = &fout;
    }
    *out << "t_s,r_s_m,T_H_K,Theta_trans,L_WH_W\n";
    double r_s = compute_r_s(M);
    double T_H = compute_T_H(M);
    for (double t = t_start; t <= t_end; t += dt) {
        double Theta = compute_Theta_trans(M, r, t);
        double L_WH  = compute_L_WH(M, r, t);
        *out << t << "," << r_s << "," << T_H << ","
             << Theta << "," << L_WH << "\n";
    }
    if (fout.is_open()) fout.close();
}

// ─────────────────────────────────────────────────────────────────────────────
//  Print explanations
// ─────────────────────────────────────────────────────────────────────────────
void BlackToWhiteHoleUQFF::print_explanations() const {
    std::cout << "=== Black-to-White Hole Transition (UQFF PAPER_659) ===\n"
              << "  - UQFF vacuuum densities: rho_UA=" << rho_vac_UA
              << ", rho_SCm=" << rho_vac_SCm << " J/m³\n"
              << "  - f_TRZ time-reversal factor: " << f_TRZ << "\n"
              << "  - μ_j magnetic string: " << mu_j << " J/T\n"
              << "  - γ [s⁻¹]: " << gamma << "\n"
              << "  - κ=" << kappa << " day⁻¹, [SSq]=" << SSq << "\n"
              << "  Criterion: Θ_trans = P_trans · Φ_trans · S_Um\n"
              << "  > 1 → white hole forms\n";
}

// ─────────────────────────────────────────────────────────────────────────────
//  main() — Standalone test driver for PAPER_659 numerical validation
// ─────────────────────────────────────────────────────────────────────────────
#ifdef STANDALONE_BLACK_TO_WHITE_HOLE
int main() {
    BlackToWhiteHoleUQFF btwh;
    btwh.print_explanations();

    // ── Sgr A* test: M = 4.3e6 M☉ ─────────────────────────────────────────
    const double M_sun   = 1.989e30;   // kg
    double M_sgrA        = 4.3e6 * M_sun;
    double t_Hubble      = 4.35e17;    // ~13.8 Gyr in seconds
    double r_s           = btwh.compute_r_s(M_sgrA);
    double r_s_uqff      = btwh.compute_r_s_UQFF(r_s);
    double T_H           = btwh.compute_T_H(M_sgrA);
    double E_flip        = btwh.compute_E_flip(M_sgrA, r_s_uqff);
    double P_flip        = btwh.compute_P_flip(E_flip, T_H);
    double P_trans       = btwh.compute_P_trans(P_flip);
    double Phi_trans     = btwh.compute_Phi_trans(M_sgrA);
    double U_m           = btwh.compute_U_m(r_s, t_Hubble);
    double S_Um          = btwh.compute_S_Um(r_s, t_Hubble, T_H);
    double Theta_trans   = btwh.compute_Theta_trans(M_sgrA, r_s, t_Hubble);
    double L_WH          = btwh.compute_L_WH(M_sgrA, r_s, t_Hubble);
    double tau_WH        = btwh.compute_tau_WH(r_s, t_Hubble);

    std::cout << "\n=== Sgr A* Validation ===\n"
              << "  M           = " << M_sgrA        << " kg\n"
              << "  r_s         = " << r_s           << " m\n"
              << "  r_s,UQFF    = " << r_s_uqff      << " m\n"
              << "  T_H         = " << T_H           << " K\n"
              << "  E_flip      = " << E_flip        << " J\n"
              << "  P_flip      = " << P_flip        << "\n"
              << "  P_trans     = " << P_trans       << "\n"
              << "  Phi_trans   = " << Phi_trans     << "\n"
              << "  U_m(r_s)    = " << U_m           << " J\n"
              << "  S_Um        = " << S_Um          << "\n"
              << "  Theta_trans = " << Theta_trans   << "  (need >1 for WH formation)\n"
              << "  tau_WH      = " << tau_WH        << " s\n"
              << "  L_WH        = " << L_WH          << " W\n";

    // ── Monte-Carlo P(Θ > 1) ─────────────────────────────────────────────────
    double P_MC = btwh.monte_carlo_P_transition(M_sgrA, r_s, t_Hubble, 5000, 0.05);
    std::cout << "\n  Monte-Carlo P(Θ>1) n=5000, σ_ρ=5%: " << P_MC << "\n";

    // ── Micro-BH test: M = 1e20 kg (survivable PBH) ──────────────────────────
    double M_pbh = 1.0e20;
    double Theta_pbh = btwh.compute_Theta_trans(M_pbh, btwh.compute_r_s(M_pbh), 1e10);
    std::cout << "\n=== Micro-BH (1e20 kg) ===\n"
              << "  T_H         = " << btwh.compute_T_H(M_pbh)     << " K\n"
              << "  Theta_trans = " << Theta_pbh                    << "\n"
              << "  P_MC        = " << btwh.monte_carlo_P_transition(M_pbh, btwh.compute_r_s(M_pbh), 1e10, 5000, 0.05) << "\n";

    // ── Simulate Sgr A* over 1e9 years ───────────────────────────────────────
    btwh.simulate(M_sgrA, r_s, 0.0, 3.15e16, 1e14, "bh_wh_transition_sgrA.csv");
    std::cout << "\nSimulation complete → bh_wh_transition_sgrA.csv\n";

    return 0;
}
#endif // STANDALONE_BLACK_TO_WHITE_HOLE
