// ============================================================================
// BlackHoleBounceUQFF.cpp — UQFF extension of LQG Black Hole Bounces (LQC)
// PAPER_658 | Session 172 | April 1–2, 2026
// Source: grok_share_fc21e30c24b4.txt — module 60 "Loop Quantum Gravity Bounce"
// ============================================================================
#include "BlackHoleBounceUQFF.h"

// ── Constructor ────────────────────────────────────────────────────────────
BlackHoleBounceUQFF::BlackHoleBounceUQFF(unsigned int seed)
    : rng(seed)
{
    // Derived Planck quantities
    // ρ_Pl = c^5 / (ħ G²)
    rho_Planck = std::pow(c, 5.0) / (hbar * G * G);

    // t_Pl = sqrt(ħ G / c^5)
    t_Planck = std::sqrt(hbar * G / std::pow(c, 5.0));

    // a_min = sqrt(ħ G / c^3) ≈ Planck length
    a_min = std::sqrt(hbar * G / std::pow(c, 3.0));

    // LQC critical density  ρ_c = 0.41 ρ_Pl
    rho_c = rho_c_factor * rho_Planck;
}

// ── Classical Friedmann (k=0) ──────────────────────────────────────────────
double BlackHoleBounceUQFF::compute_classical_friedmann(double rho) const
{
    // (ȧ/a)² = (8πG/3) ρ
    // Singularity: ȧ/a → ∞ as ρ → ∞
    return (8.0 * M_PI * G / 3.0) * rho;
}

// ── LQC-modified Friedmann ─────────────────────────────────────────────────
double BlackHoleBounceUQFF::compute_LQC_friedmann(double rho, double a) const
{
    // (ȧ/a)² = (8πG/3) ρ (1 − ρ/ρ_c) − k c² / a²
    // When ρ → ρ_c the bracket → 0 ⟹ ȧ = 0 ⟹ bounce (not singularity)
    double H2 = (8.0 * M_PI * G / 3.0) * rho * (1.0 - rho / rho_c);
    if (k_curv != 0 && a > 0.0) {
        H2 -= static_cast<double>(k_curv) * c * c / (a * a);
    }
    // Apply any extra terms registered at runtime (self-expanding framework)
    for (const auto& term : extra_terms) {
        H2 += term(a, rho);
    }
    return H2;
}

// ── UQFF-extended critical density ────────────────────────────────────────
double BlackHoleBounceUQFF::compute_rho_c_UQFF() const
{
    // ρ_c,UQFF = ρ_c · (1 + ρ_vac_UA / ρ_vac_SCm)
    // Vacuum Density Series (PAPER_646): ratio ρ_UA/ρ_SCm = 10 for the first two terms
    return rho_c * (1.0 + rho_vac_UA / rho_vac_SCm);
}

// ── Near-bounce scale factor ───────────────────────────────────────────────
double BlackHoleBounceUQFF::compute_scale_factor_near_bounce(double t) const
{
    // a(t) ≈ a_min · cosh(t / t_Pl)
    // Valid near t=0 (bounce moment); symmetric contraction / expansion
    return a_min * std::cosh(t / t_Planck);
}

// ── UQFF-corrected scale factor ────────────────────────────────────────────
double BlackHoleBounceUQFF::compute_scale_factor_UQFF(double t) const
{
    // UQFF correction: the Buoyancy Harmonics (PAPER_648) raise a_min by
    // the aether buoyancy factor.
    // a(t,UQFF) = a_min · cosh(t/t_Pl) · (1 + f_TRZ · ρ_UA/ρ_SCm)^{1/3}
    double base = a_min * std::cosh(t / t_Planck);
    double uqff_factor = std::cbrt(1.0 + f_TRZ * rho_vac_UA / rho_vac_SCm);
    return base * uqff_factor;
}

// ── Effective equation-of-state ────────────────────────────────────────────
double BlackHoleBounceUQFF::compute_effective_eos() const
{
    // w_eff = −1 + (1+f_TRZ)(ρ_vac_UA/ρ_vac_SCm) · κ · [SSq]
    // Captures UQFF negentropic contribution near bounce
    return -1.0 + (1.0 + f_TRZ) * (rho_vac_UA / rho_vac_SCm) * kappa * SSq;
}

// ── Density rate-of-change ─────────────────────────────────────────────────
double BlackHoleBounceUQFF::compute_rho_rate(double rho, double a) const
{
    // drho/dt = −3 H (1+w) ρ,  w ≈ 0 for dust, 1/3 for radiation
    // Use UQFF effective w near bounce
    double H2 = compute_LQC_friedmann(rho, a);
    double H  = (H2 > 0.0) ? std::sqrt(H2) : 0.0;
    double w  = compute_effective_eos();
    return -3.0 * H * (1.0 + w) * rho;
}

// ── Self-expanding: register extra term ───────────────────────────────────
void BlackHoleBounceUQFF::add_extra_term(
    std::function<double(double, double)> term)
{
    extra_terms.push_back(std::move(term));
}

// ── Self-update: load parameters from file ────────────────────────────────
void BlackHoleBounceUQFF::update_from_file(const std::string& filename)
{
    std::ifstream f(filename);
    if (!f.is_open()) {
        std::cerr << "[BlackHoleBounceUQFF] Cannot open: " << filename << "\n";
        return;
    }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto pos = line.find('=');
        if (pos == std::string::npos) continue;
        std::string key = line.substr(0, pos);
        double val = std::stod(line.substr(pos + 1));
        if (key == "G")              { G = val; }
        else if (key == "c")         { c = val; }
        else if (key == "hbar")      { hbar = val; }
        else if (key == "rho_c_factor") {
            rho_c_factor = val;
            rho_Planck = std::pow(c, 5.0) / (hbar * G * G);
            rho_c = rho_c_factor * rho_Planck;
            t_Planck = std::sqrt(hbar * G / std::pow(c, 5.0));
            a_min    = std::sqrt(hbar * G / std::pow(c, 3.0));
        }
        else if (key == "rho_vac_UA")  { rho_vac_UA  = val; }
        else if (key == "rho_vac_SCm") { rho_vac_SCm = val; }
        else if (key == "f_TRZ")       { f_TRZ  = val; }
        else if (key == "kappa")       { kappa  = val; }
        else if (key == "SSq")         { SSq    = val; }
        else if (key == "k_curv")      { k_curv = static_cast<int>(val); }
    }
    f.close();
    std::cout << "[BlackHoleBounceUQFF] Parameters updated from " << filename << "\n";
}

// ── Self-simulate ─────────────────────────────────────────────────────────
void BlackHoleBounceUQFF::simulate(
    double a0, double rho0,
    double t_start, double t_end, double dt,
    const std::string& output_file) const
{
    std::ofstream fout;
    bool to_file = !output_file.empty();
    if (to_file) {
        fout.open(output_file);
        fout << "t_s,a_m,rho_kg_m3,H2_s-2,a_UQFF_m\n";
    } else {
        std::cout << "t_s,a_m,rho_kg_m3,H2_s-2,a_UQFF_m\n";
    }

    double a   = a0;
    double rho = rho0;
    for (double t = t_start; t <= t_end; t += dt) {
        double H2    = compute_LQC_friedmann(rho, a);
        double H     = (H2 > 0.0) ? std::sqrt(H2) : 0.0;
        double aUQFF = compute_scale_factor_UQFF(t - t_start);

        if (to_file) {
            fout << t << "," << a << "," << rho << "," << H2 << "," << aUQFF << "\n";
        } else {
            std::cout << t << "," << a << "," << rho << "," << H2 << "," << aUQFF << "\n";
        }

        // Euler step: ȧ = H · a,   drho/dt from continuity
        double drho_dt = compute_rho_rate(rho, a);
        a   += H * a * dt;
        rho += drho_dt * dt;
        if (rho < 0.0) rho = 0.0;
    }
    if (to_file) fout.close();
}

// ── Print explanations ─────────────────────────────────────────────────────
void BlackHoleBounceUQFF::print_explanations() const
{
    std::cout
        << "=== BlackHoleBounceUQFF === PAPER_658 ===\n"
        << "LQG replaces singularities with smooth bounces via discrete spacetime.\n"
        << "LQC modified Friedmann: (da/dt/a)^2 = (8piG/3)rho(1-rho/rho_c) - kc^2/a^2\n"
        << "Critical density rho_c = " << rho_c << " kg/m^3  (0.41 x Planck density)\n"
        << "Planck density  rho_Pl = " << rho_Planck << " kg/m^3\n"
        << "Planck time     t_Pl   = " << t_Planck << " s\n"
        << "Min scale       a_min  = " << a_min    << " m  (Planck length)\n"
        << "Near-bounce: a(t) ~ a_min * cosh(t/t_Pl)\n"
        << "UQFF rho_c: " << compute_rho_c_UQFF() << " kg/m^3  (x11 via VDS PAPER_646)\n"
        << "Eff EOS w_eff: " << compute_effective_eos() << "\n"
        << "Buoyancy Harmonics (PAPER_648): BH aether gradient prevents singularity.\n"
        << "Dipole Vortex Primes (PAPER_647): U_m oscillator stabilises BH post-bounce.\n"
        << "\n";
}

// ── main() — test driver ───────────────────────────────────────────────────
int main()
{
    BlackHoleBounceUQFF lqc;
    lqc.print_explanations();

    // 1. Compare classical vs LQC Friedmann at high density
    double rho_test = 0.9 * lqc.get_rho_c();   // 90% of critical density
    double a_test   = 2.0 * lqc.get_a_min();
    double H2_classical = lqc.compute_classical_friedmann(rho_test);
    double H2_LQC       = lqc.compute_LQC_friedmann(rho_test, a_test);
    std::cout << "Classical H^2  (rho=0.9rho_c): " << H2_classical << " s^-2\n";
    std::cout << "LQC      H^2  (rho=0.9rho_c): " << H2_LQC       << " s^-2  (reduced)\n";

    // 2. Near-bounce scale factor at t = 0 (bounce moment)
    double t_bounce = 0.0;
    std::cout << "a(t=0) near bounce  = " << lqc.compute_scale_factor_near_bounce(t_bounce)
              << " m  (should be a_min = " << lqc.get_a_min() << " m)\n";

    // 3. UQFF corrected scale factor 1 Planck time after bounce
    double t_1tPl = lqc.get_t_Planck();
    std::cout << "a_UQFF(t=t_Pl) = " << lqc.compute_scale_factor_UQFF(t_1tPl) << " m\n";

    // 4. UQFF critical density (elevated by VDS ratio)
    std::cout << "rho_c (standard) = " << lqc.get_rho_c()      << " kg/m^3\n";
    std::cout << "rho_c (UQFF)     = " << lqc.get_rho_c_UQFF() << " kg/m^3\n";

    // 5. Self-expand: add a small aether damping term proportional to a^{-6}
    //    (stiff fluid contribution from [UA] vortex lattice)
    lqc.add_extra_term([](double a, double /*rho*/) -> double {
        const double lambda_UA = 1e-200;   // illustrative small constant
        return (a > 0.0) ? lambda_UA / (a * a * a * a * a * a) : 0.0;
    });

    // 6. Short simulation: collapse start → bounce → expansion (Planck units)
    double a0   = 1.0e-33;           // slightly above a_min, contracting universe
    double rho0 = 0.3 * lqc.get_rho_c();  // well below critical density
    double t_start = -20.0 * lqc.get_t_Planck();
    double t_end   =  20.0 * lqc.get_t_Planck();
    double dt      =   0.5 * lqc.get_t_Planck();

    std::cout << "\nSimulating LQC bounce (t = " << t_start << " to " << t_end << " s):\n";
    lqc.simulate(a0, rho0, t_start, t_end, dt, "lqc_bounce_sim.csv");
    std::cout << "Simulation written to lqc_bounce_sim.csv\n";

    return 0;
}
