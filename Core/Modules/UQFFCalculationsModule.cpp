// UQFFCalculationsModule.cpp
// Implementation of UQFF Core calculations for 5 astrophysical systems
// Based on "5 Astro Systems_B_11Oct2025.docx"
// Systems: M82 (Messier 82), IC418 (Spirograph Nebula), Canis Major (R136),
//          NGC6302 (Butterfly Nebula), NGC7027
// Watermark: Copyright - Daniel T. Murphy, analyzed by Grok 3, dated November 17, 2025

#include "../UQFFCalculationsModule.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <string>

// ============================================================
// UQFFCore: Constructor
// ============================================================
UQFFCore::UQFFCore(double k1, double kg3, double km, double ke)
    : k1_(k1), kg3_(kg3), km_(km), ke_(ke) {}

// ============================================================
// calculate_U_g1: DPM Force (Coulomb-analog)
// U_g1 = k1 * Sum_j [f_UA'(Z_j) / r_j^2] * f_nu_THz
// ============================================================
double UQFFCore::calculate_U_g1(const std::vector<DPMVars>& vars,
                                 GeometryType geom) {
    double total = 0.0;
    for (const auto& v : vars) {
        double f_UA_prime = (UQFFCALC_Z_MAX - v.Z) / UQFFCALC_Z_MAX;
        double f_nu = 1.0 + std::sin(UQFFCALC_PI * v.nu_THz / UQFFCALC_NU_THz);
        double r_eff = (v.r > 0) ? v.r : 1.0;

        if (geom == SPHERICAL) {
            total += k1_ * (f_UA_prime / (r_eff * r_eff)) * f_nu;
        } else {
            // Toroidal: factor of 2*PI * r for toroidal shell
            total += k1_ * (f_UA_prime / (r_eff * r_eff)) * f_nu * 2.0 * UQFFCALC_PI * r_eff;
        }
    }
    return total;
}

// ============================================================
// calculate_U_m: Universal Magnetism (Heaviside + quasi-factor)
// U_m = km * B * f_UA' * f_SCm * H(f_UA')
// ============================================================
double UQFFCore::calculate_U_m(const DPMVars& v) {
    double f_UA_prime = (UQFFCALC_Z_MAX - v.Z) / UQFFCALC_Z_MAX;
    double f_SCm = v.Z / UQFFCALC_Z_MAX;
    double H_step = (f_UA_prime >= 0) ? 1.0 : 0.0; // Heaviside
    return km_ * v.B * f_UA_prime * f_SCm * H_step;
}

// ============================================================
// calculate_U_g3: Composite force (U_i + U_m)
// U_g3 = kg3 * (U_i + U_m)
// U_i is the electrostatic interaction using R_EB barrier
// ============================================================
double UQFFCore::calculate_U_g3(const DPMVars& v) {
    double R_EB = UQFFCALC_K_R * v.Z;
    double U_i = R_EB * v.E;
    double U_m = calculate_U_m(v);
    return kg3_ * (U_i + U_m);
}

// ============================================================
// calculate_E: Electric field from DPM charge density
// E = k * rho_vac / (r^2) * f_UA' * N_quantum
// ============================================================
double UQFFCore::calculate_E(const DPMVars& v) {
    double f_UA_prime = (UQFFCALC_Z_MAX - v.Z) / UQFFCALC_Z_MAX;
    double r_eff = (v.r > 0) ? v.r : 1.0;
    return ke_ * UQFFCALC_RHO_VAC_UA / (r_eff * r_eff) * f_UA_prime * UQFFCALC_N_QUANTUM;
}

// ============================================================
// calculate_eta: Neutron production rate
// eta = K_ETA_BASE * f_UA' * f_SCm * sqrt(B / RHO_VAC_UA)
// ============================================================
double UQFFCore::calculate_eta(const DPMVars& v) {
    double f_UA_prime = (UQFFCALC_Z_MAX - v.Z) / UQFFCALC_Z_MAX;
    double f_SCm = v.Z / UQFFCALC_Z_MAX;
    return UQFFCALC_K_ETA_BASE * f_UA_prime * f_SCm *
           std::sqrt(std::fabs(v.B / UQFFCALC_RHO_VAC_UA));
}

// ============================================================
// UQFFSystem: Constructor
// ============================================================
UQFFSystem::UQFFSystem(const SystemParams& params) : params_(params) {
    // Build default DPMVars from system parameters
    DPMVars v;
    v.Z      = 26.0;  // Iron-like DPM atomic number (astrophysical approximation)
    v.nu_THz = 1.0e12;
    v.theta  = UQFFCALC_PI / 4.0;
    v.phi    = 0.0;
    v.r      = params.r;
    v.B      = params.B;
    v.E      = 0.0; // Will be computed
    vars_.push_back(v);
}

double UQFFSystem::calculate_master_force(const UQFFCore& core) {
    double U_g1 = core.calculate_U_g1(vars_, params_.geom);
    double U_g3 = core.calculate_U_g3(vars_[0]);
    return U_g1 + U_g3;
}

// ============================================================
// Factory Functions for 5 systems
// ============================================================

// M82: Messier 82 (Cigar Galaxy) — starburst, distance ~3.5 Mpc
UQFFSystem create_M82_system() {
    SystemParams p = {1.0e20, 5.0, 1e-5, 0.00067, 9.46e13, SPHERICAL, "M82 (Messier 82)"};
    return UQFFSystem(p);
}

// IC418: Spirograph Nebula — planetary nebula, distance ~2 kpc
UQFFSystem create_IC418_system() {
    SystemParams p = {1.0e16, 0.0, 1e-5, 0.00014, 3.15e12, SPHERICAL, "IC418 (Spirograph Nebula)"};
    return UQFFSystem(p);
}

// Canis Major R136: super star cluster, distance ~50 kpc
UQFFSystem create_CanisMajor_system() {
    SystemParams p = {3.0e20, 7.5, 1e-4, 0.00016, 4.73e13, SPHERICAL, "Canis Major (R136)"};
    return UQFFSystem(p);
}

// NGC6302: Butterfly Nebula — extreme planetary nebula, distance ~3.9 kpc
UQFFSystem create_NGC6302_system() {
    SystemParams p = {2.5e16, 0.0, 1e-5, 0.00034, 6.31e12, TOROIDAL, "NGC6302 (Butterfly Nebula)"};
    return UQFFSystem(p);
}

// NGC7027: Carbon-rich PN — distance ~1 kpc
UQFFSystem create_NGC7027_system() {
    SystemParams p = {9.46e15, 0.1, 1e-5, 0.001, 3.15e10, SPHERICAL, "NGC7027"};
    return UQFFSystem(p);
}
