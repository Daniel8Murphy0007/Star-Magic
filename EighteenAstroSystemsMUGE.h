// EighteenAstroSystemsMUGE.h
// Entry #105: (18 astro-systems)_cpp_09May2025
// MUGE Compressed + Resonance UQFF equations for 18 astrophysical systems:
//   (10 from entry #101) PLUS:
//  11. NGC 6217  (Spiral Galaxy)
//  12. Stephan's Quintet (NGC 7317/7318a/7318b/7319/7320)
//  13. NGC 7049  (Elliptical Galaxy)
//  14. Carina Nebula (NGC 3324 region - Cosmic Cliffs)
//  15. M74       (NGC 628, face-on spiral)
//  16. NGC 1672  (Barred Spiral Galaxy)
//  17. NGC 5866  (Lenticular Galaxy / Spindle Galaxy)
//  18. M82       (Messier 82, Cigar Galaxy, starburst)
//  Also includes: Spirograph Nebula (IC 418) as the 19th system
//  named entry calls it "18" because Mystic Mountain is counted once
//
// Source: grok_share_ba508f76c8e.txt entry #105 | Session 178
// Tag: STANDALONE_EIGHTEEN_ASTRO // PAPER_733 / CP4 #317

#ifndef EIGHTEEN_ASTRO_SYSTEMS_MUGE_H
#define EIGHTEEN_ASTRO_SYSTEMS_MUGE_H

#include <cmath>
#include <string>
#include <vector>
#include <iostream>

// Re-use the same parameter struct (no cross-include needed)
struct AstroSystem105 {
    const char* name;
    double M_kg;        // total mass (kg)
    double r_m;         // characteristic radius (m)
    double z;           // redshift
    double SFR;         // star-formation rate (M_sun/yr)
    double B_T;         // magnetic field (T)
    double v_wind;      // wind velocity (m/s)
    double tau_erode;   // erosion timescale (s)
    double E_0;         // mass-loss amplitude
};

class EighteenAstroSystemsMUGE {
public:
    // ----------------------------------------------------------
    // UQFF Universal Constants (same set as all modules)
    // ----------------------------------------------------------
    static constexpr double G        = 6.6743e-11;
    static constexpr double c        = 3.0e8;
    static constexpr double hbar     = 1.0546e-34;
    static constexpr double mu_0     = 1.2566e-6;
    static constexpr double k_B      = 1.3806e-23;
    static constexpr double M_sun    = 1.989e30;
    static constexpr double kpc      = 3.086e19;
    static constexpr double Mpc      = 3.086e22;
    static constexpr double rho_UA   = 7.09e-36;
    static constexpr double rho_SCm  = 7.09e-37;
    static constexpr double f_TRZ    = 0.1;
    static constexpr double kappa    = 5.0e-4;
    static constexpr double SSq      = 0.57;
    static constexpr double H_0      = 2.269e-18;
    static constexpr double t_H      = 4.35e17;
    static constexpr double Lambda   = 1.1e-52;
    static constexpr double q_e      = 1.602e-19;
    static constexpr double m_p      = 1.673e-27;
    static constexpr double em_scale = 1.0e-12;
    static constexpr double aether_ratio = rho_UA / rho_SCm; // 10.0

    // 19 systems (entry #105 adds 9 to the 10 in entry #101;
    // "18" excludes the Spirograph Nebula overlap in naming)
    static constexpr int N_SYSTEMS = 19;
    static const AstroSystem105 SYSTEMS[19];

    EighteenAstroSystemsMUGE() : version("Session178") {}

    // ----------------------------------------------------------
    // Core MUGE methods
    // ----------------------------------------------------------
    double H_z(double z) const;
    double M_evo(double t, double SFR, double M_kg) const;
    double E_rad(double t, double E_0, double tau_erode) const;
    double g_MUGE(const AstroSystem105& sys, double t) const;
    double R_resonance(const AstroSystem105& sys, double t) const;

    std::vector<double> computeAll(double t) const;
    void simulate();

    void self_update();
    void self_expand();
    std::string primary_equation_str() const;
    std::string description() const;

    const std::string version;
};

#endif // EIGHTEEN_ASTRO_SYSTEMS_MUGE_H
