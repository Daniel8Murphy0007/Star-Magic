// EighteenAstroSystemsMUGE.cpp
// Entry #105: (18 astro-systems)_cpp_09May2025
//
// =============================================================
// PHYSICS SUMMARY (Master MUGE Compressed UQFF — 18 systems)
// =============================================================
// Same MUGE compressed framework as entry #101, extended to
// 18 unique astrophysical systems.  The 26D Ug1-Ug4i extended
// Quantum State terms are included per the full UQFF expansion:
//
//   E_DPM,i = (hbar*c / r_i^2) * Q_i * [SCm]_i
//             where r_i = r/i, Q_i = i, [SCm]_i = 1e-5 * i^2  (i=1..26)
//
//   Ug1_i = E_DPM,i * (1+H(z)*t) * (1-E_rad) * cos(theta_i) * (1+f_TRZ,i)
//   Ug2_i = E_DPM,i * (1-B/B_crit) * (1+M_sf) * (1+rho_UA/rho_SCm) * Sigma_j cos(omega_j*t)
//   Ug3_i = E_DPM,i * (q*v*B/m_p) * (1-T_lock) * (1+f_TRZ,i)
//   Ug4i_i= (hbar*c / r_THz,i) * (1+f_Um,i) * (1+rho_UA/rho_SCm)
//
// Full MUGE:
//   g(r,t) = G*M/r^2 * (1+H(z)*t) * (1+M_evo) * (1-E_rad) * (1+f_TRZ)
//          + F_em + Sigma_i (Ug1_i + Ug2_i + Ug3_i + Ug4i_i)
//
// 18 Systems (10 from entry #101 + 8 new JWST deep-field targets):
//
// System                     | g_MUGE (m/s^2)
// ---------------------------+----------------
// NGC 2264                   | ~1.05e-2  (EM-dominated, SFR=0.5)
// UGC 10214 (Tadpole)        | ~1.2e-11  (tidal tail, z=0.028)
// NGC 4676 (Mice Galaxies)   | ~2.0e-11  (dual mass, merger)
// Red Spider Nebula           | ~2.5e-9   (fast winds, PN)
// NGC 3372 (Carina)          | ~1.8e-9   (star-forming complex)
// AG Carinae Nebula          | ~1.6e-9   (LBV wind)
// M42 (Orion)                | ~8.4e-10  (protostellar)
// Tarantula Nebula           | ~2.2e-9   (starburst LMC)
// NGC 2841 (spiral)          | ~5.0e-11  (Milgrom regime)
// Mystic Mountain            | ~1.3e-9   (Carina pillar)
// NGC 6217 (barred spiral)   | ~4.9e-11
// Stephan's Quintet          | ~3.0e-11
// NGC 7049 (lenticular)      | ~4.5e-11
// Carina Neb. NGC 3324       | ~1.7e-9
// M74 (grand design spiral)  | ~5.1e-11
// NGC 1672 (Seyfert barred)  | ~4.7e-11
// NGC 5866 (edge-on lenticular)| ~4.6e-11
// M82 (Cigar, starburst)     | ~3.8e-10  (SFR=10)
// =============================================================

#include "EighteenAstroSystemsMUGE.h"

// ============================================================
// System definitions (entry #105, 18 unique systems + Spirograph = 19 entries)
// ============================================================
const AstroSystem105 EighteenAstroSystemsMUGE::SYSTEMS[N_SYSTEMS] = {
    // -------- 10 systems from entry #101 --------
    { "NGC 2264 (Cone Nebula)",          1.989e33, 4.73e16, 0.0006, 0.5,  1.0e-5, 1.0e6,  6.312e13, 0.20 },
    { "UGC 10214 (Tadpole Galaxy)",      1.989e41, 1.24e21, 0.028,  1.0,  1.0e-5, 2.0e5,  3.156e16, 0.05 },
    { "NGC 4676 (Mice Galaxies)",        3.978e41, 3.0e20,  0.022,  10.0, 1.0e-4, 3.0e5,  3.156e15, 0.10 },
    { "Red Spider Nebula",               1.193e30, 1.0e16,  0.0013, 0.0,  1.0e-5, 2.0e6,  3.156e11, 0.15 },
    { "NGC 3372 (Carina Nebula)",        1.989e35, 2.0e17,  0.0025, 2.0,  1.0e-5, 1.5e6,  6.312e13, 0.20 },
    { "AG Carinae Nebula",               3.978e31, 1.0e16,  0.002,  0.0,  1.0e-5, 2.0e6,  3.156e11, 0.15 },
    { "M42 (Orion Nebula)",              3.978e33, 2.0e16,  0.0004, 0.3,  1.0e-5, 1.0e6,  6.312e13, 0.10 },
    { "Tarantula Nebula (30 Doradus)",   1.989e35, 3.0e17,  0.0005, 5.0,  1.0e-4, 2.0e6,  3.156e13, 0.25 },
    { "NGC 2841 (Spiral Galaxy)",        1.989e41, 5.0e20,  0.0031, 0.5,  1.0e-5, 1.5e5,  3.156e16, 0.05 },
    { "Mystic Mountain (Carina pillar)", 1.989e32, 1.0e16,  0.0025, 0.1,  1.0e-5, 1.0e6,  3.156e12, 0.20 },
    // -------- 8 new systems from entry #105 (JWST deep-field) --------
    { "NGC 6217 (barred spiral)",        1.989e41, 3.0e20,  0.0045, 1.0,  1.0e-5, 1.5e5,  3.156e16, 0.05 },
    { "Stephan's Quintet",               9.945e41, 1.0e21,  0.022,  10.0, 1.0e-4, 3.0e5,  3.156e15, 0.10 },
    { "NGC 7049 (lenticular)",           1.989e41, 5.0e20,  0.0067, 0.2,  1.0e-5, 1.0e5,  3.156e16, 0.05 },
    { "Carina Nebula (NGC 3324)",        1.989e35, 2.0e17,  0.0025, 2.0,  1.0e-5, 1.5e6,  6.312e13, 0.20 },
    { "M74 (NGC 628 grand-design)",      1.989e41, 5.0e20,  0.0022, 1.0,  1.0e-5, 1.5e5,  3.156e16, 0.05 },
    { "NGC 1672 (Seyfert barred)",       1.989e41, 3.0e20,  0.004,  2.0,  1.0e-5, 2.0e5,  3.156e15, 0.08 },
    { "NGC 5866 (edge-on lenticular)",   1.989e41, 3.0e20,  0.0029, 0.3,  1.0e-5, 1.0e5,  3.156e16, 0.05 },
    { "M82 (Cigar Galaxy, starburst)",   1.989e40, 2.0e20,  0.0008, 10.0, 1.0e-4, 5.0e5,  3.156e14, 0.20 },
    { "Spirograph Nebula (IC 418)",      1.193e30, 1.0e16,  0.0007, 0.0,  1.0e-5, 1.5e6,  3.156e11, 0.12 }
};

// ============================================================
// H(z): Hubble parameter at redshift z
// ============================================================
double EighteenAstroSystemsMUGE::H_z(double z) const {
    double a = 1.0 + z;
    return H_0 * std::sqrt(0.3 * a * a * a + 0.7);
}

// ============================================================
// M_evo: fractional mass growth from SFR
// ============================================================
double EighteenAstroSystemsMUGE::M_evo(double t, double SFR, double M_kg) const {
    double t_yr  = t / 3.156e7;
    double M_sun = M_kg / 1.989e30;
    return SFR * t_yr / M_sun;
}

// ============================================================
// E_rad: radiation/merger erosion factor
// ============================================================
double EighteenAstroSystemsMUGE::E_rad(double t, double E_0, double tau_erode) const {
    return E_0 * (1.0 - std::exp(-t / tau_erode));
}

// ============================================================
// Ug_26D: 26D quantum state contribution (Ug1+Ug2+Ug3+Ug4i summed over i)
// ============================================================
double EighteenAstroSystemsMUGE::Ug_26D(const AstroSystem105& sys, double t) const {
    double sum = 0.0;
    double hub   = 1.0 + H_z(sys.z) * t;
    double e_fac = 1.0 - E_rad(t, sys.E_0, sys.tau_erode);
    double m_fac = 1.0 + M_evo(t, sys.SFR, sys.M_kg);
    double T_lock = 0.0; // tidal-locking factor (0 = symmetric orbit)
    double B_ratio = sys.B_T / B_crit_magnetar; // B / B_crit (<<1 for most systems)

    for (int i = 1; i <= 26; ++i) {
        double r_i   = sys.r_m / static_cast<double>(i);
        double Q_i   = static_cast<double>(i);
        double SCm_i = 1.0e-5 * static_cast<double>(i * i);
        double theta_i = static_cast<double>(i) * M_PI / 26.0; // phase angle per layer
        double f_TRZ_i = f_TRZ * (1.0 + 0.01 * static_cast<double>(i));

        double E_DPM_i = (hbar * c_light / (r_i * r_i)) * Q_i * SCm_i;

        double Ug1 = E_DPM_i * hub * e_fac * std::cos(theta_i) * (1.0 + f_TRZ_i);
        double Ug2 = E_DPM_i * (1.0 - B_ratio) * m_fac * (1.0 + aether_ratio)
                     * std::cos(static_cast<double>(i) * 6.28e-13 * t);
        double Ug3 = E_DPM_i * (q_e * sys.v_wind * sys.B_T / m_p) * (1.0 - T_lock) * (1.0 + f_TRZ_i);
        double r_THz_i = r_i * 1.0e-3; // THz scale inner radius
        double f_Um_i  = 0.01 * static_cast<double>(i);
        double Ug4i = (hbar * c_light / (r_THz_i * r_THz_i)) * (1.0 + f_Um_i) * (1.0 + aether_ratio);

        sum += Ug1 + Ug2 + Ug3 + Ug4i;
    }
    return sum;
}

// ============================================================
// g_MUGE: Master Universal Gravity Equation (Compressed UQFF + 26D)
// ============================================================
double EighteenAstroSystemsMUGE::g_MUGE(const AstroSystem105& sys, double t) const {
    double g_grav = G * sys.M_kg / (sys.r_m * sys.r_m);

    double hub   = 1.0 + H_z(sys.z) * t;
    double m_fac = 1.0 + M_evo(t, sys.SFR, sys.M_kg);
    double e_fac = 1.0 - E_rad(t, sys.E_0, sys.tau_erode);
    double trz   = 1.0 + f_TRZ;

    double F_em = q_e * sys.v_wind * sys.B_T / m_p
                  * (1.0 + aether_ratio) * em_scale;

    double ug26 = Ug_26D(sys, t);

    return g_grav * hub * m_fac * e_fac * trz + F_em + ug26 * 1.0e-40;
}

// ============================================================
// R_resonance: Master Resonance UQFF Equation
// ============================================================
double EighteenAstroSystemsMUGE::R_resonance(const AstroSystem105& sys, double t) const {
    double g_grav  = G * sys.M_kg / (sys.r_m * sys.r_m);
    double m_fac   = 1.0 + M_evo(t, sys.SFR, sys.M_kg);
    double R_grav  = g_grav * m_fac;

    double R_mag   = q_e * sys.v_wind * sys.B_T / m_p * em_scale;

    double t_sf     = sys.tau_erode;
    double omega_g  = (t_sf > 0.0) ? (2.0 * M_PI / t_sf) : 1.0e-13;
    double omega_m  = omega_g * 100.0;

    return R_grav * std::cos(omega_g * t)
         + R_mag * std::cos(omega_m * t) * aether_ratio * (1.0 + f_TRZ);
}

// ============================================================
// computeAll: return vector of g values for all N_SYSTEMS
// ============================================================
std::vector<double> EighteenAstroSystemsMUGE::computeAll(double t) const {
    std::vector<double> results;
    for (int i = 0; i < N_SYSTEMS; ++i) {
        results.push_back(g_MUGE(SYSTEMS[i], t));
    }
    return results;
}

// ============================================================
// simulate: 5-step self-simulation
// ============================================================
void EighteenAstroSystemsMUGE::simulate() {
    std::cout << "=== EighteenAstroSystemsMUGE Simulation (PAPER_733) ===\n";
    double t_step = 3.156e14; // 10 Myr
    double t = t_step;
    for (int step = 1; step <= 5; ++step) {
        std::cout << "Step " << step << " (t=" << t/3.156e7 << " Myr):\n";
        for (int i = 0; i < N_SYSTEMS; ++i) {
            double g = g_MUGE(SYSTEMS[i], t);
            std::cout << "  " << SYSTEMS[i].name
                      << ": g = " << g << " m/s^2\n";
        }
        t += t_step;
    }
}

void EighteenAstroSystemsMUGE::self_update() { /* additive only */ }
void EighteenAstroSystemsMUGE::self_expand() { /* reserved */ }

std::string EighteenAstroSystemsMUGE::primary_equation_str() const {
    return "g(r,t) = G*M/r^2 * (1+H(z)*t) * (1+M_evo) * (1-E_rad) * (1+f_TRZ) + F_em + Sigma_26D[Ug1+Ug2+Ug3+Ug4i]";
}

std::string EighteenAstroSystemsMUGE::description() const {
    return "MUGE Compressed + 26D UQFF for 18/19 astro-systems (entry #105, PAPER_733)";
}

// ============================================================
// main() — MUGE + resonance test for all 19 systems
// ============================================================
#ifdef STANDALONE_EIGHTEEN_ASTRO
int main() {
    EighteenAstroSystemsMUGE mod;
    double t = 3.156e14; // 10 Myr

    std::cout << "=== Entry #105: 18 Astro-Systems MUGE + 26D UQFF Test ===\n";
    std::cout << "Eq: " << mod.primary_equation_str() << "\n\n";

    std::cout << "--- g_MUGE at t=10 Myr ---\n";
    for (int i = 0; i < EighteenAstroSystemsMUGE::N_SYSTEMS; ++i) {
        double g = mod.g_MUGE(EighteenAstroSystemsMUGE::SYSTEMS[i], t);
        double R = mod.R_resonance(EighteenAstroSystemsMUGE::SYSTEMS[i], t);
        std::cout << i+1 << ". " << EighteenAstroSystemsMUGE::SYSTEMS[i].name
                  << "\n   g = " << g << " m/s^2"
                  << "   R = " << R << " m/s^2\n";
    }

    std::cout << "\n";
    mod.simulate();

    return 0;
}
#endif // STANDALONE_EIGHTEEN_ASTRO
