// UQFFMultiAstroSystemsModule.cpp
// Implementation of UQFF Multi-Astrophysical Systems (11 systems)
// Based on "8 Astro Systems_B_11Oct2025.docx" (final version, renamed with DeepSearch additions)
// Simultaneous Compressed, Resonance, and Buoyancy UQFF solutions for:
// NGC4826, NGC1805, NGC6307, NGC7027, Cassini Encke, Cassini Division, Cassini Maxwell,
// ESO391-12, M57, LMC, ESO510-G13
// Watermark: Copyright - Daniel T. Murphy, analyzed by Grok 3, dated November 17, 2025

#include "../UQFFMultiAstroSystemsModule.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std::complex_literals;

// ============================================================
// UQFFMultiAstroCore: Constructor
// ============================================================
UQFFMultiAstroCore::UQFFMultiAstroCore(double k_comp, double k_res, double k_buoy)
    : k_comp_(k_comp), k_res_(k_res), k_buoy_(k_buoy) {}

// ============================================================
// F_compressed: Compressed UQFF force
// F = k_comp * [rho_vac * r^2 * H(z) * E_rad * (1 + i * B * r / c)]
// ============================================================
std::complex<double> UQFFMultiAstroCore::F_compressed(const MultiAstroParams& p) {
    double r_eff = (p.r > 0) ? p.r : 1.0;
    double h_z = hubble_correction(p.z); // (1 + z)
    double re = k_comp_ * MULTI_RHO_VAC_UA * r_eff * r_eff * h_z * MULTI_E_RAD_FACTOR;
    double im = k_comp_ * MULTI_RHO_VAC_UA * p.B * r_eff / MULTI_C * h_z;
    return std::complex<double>(re, im);
}

// ============================================================
// F_resonance: Resonance UQFF force
// F = k_res * [rho_vac * B * (1 + z) * sin(omega_THz * t) * (1 + i*sfr/c)]
// ============================================================
std::complex<double> UQFFMultiAstroCore::F_resonance(const MultiAstroParams& p) {
    double h_z = hubble_correction(p.z);
    double omega_THz = 2.0 * MULTI_PI * 1.0e12; // THz angular frequency
    double t_eff = (p.t_age > 0) ? p.t_age : 1.0;
    double phase = std::sin(omega_THz * t_eff * 1e-30); // normalized to avoid overflow
    double re = k_res_ * MULTI_RHO_VAC_UA * p.B * h_z * phase;
    double im = k_res_ * MULTI_RHO_VAC_UA * p.sfr / MULTI_C * h_z;
    return std::complex<double>(re, im);
}

// ============================================================
// F_buoyancy: Buoyancy-dominant UQFF force
// F = k_buoy * [rho_vac * r^3 / (r^2) * (1 + z)^2 * (1 + i*B*sfr)]
// Vacuum displacement: rho_vac * volume × H(z)^2
// ============================================================
std::complex<double> UQFFMultiAstroCore::F_buoyancy(const MultiAstroParams& p) {
    double r_eff = (p.r > 0) ? p.r : 1.0;
    double h_z = hubble_correction(p.z);
    double re = k_buoy_ * MULTI_RHO_VAC_UA * r_eff * h_z * h_z;
    double im = k_buoy_ * MULTI_RHO_VAC_UA * p.B * p.sfr * h_z;
    return std::complex<double>(re, im);
}

// ============================================================
// F_DPM_creation: DPM pair-creation simulation
// Rate = rho_vac * c / (hbar * r^2) * (sfr + 1) * (1 + z)
// ============================================================
double UQFFMultiAstroCore::F_DPM_creation(const MultiAstroParams& p) {
    double r_eff = (p.r > 0) ? p.r : 1.0;
    double h_z = hubble_correction(p.z);
    return MULTI_RHO_VAC_UA * MULTI_C / (MULTI_HBAR * r_eff * r_eff)
           * (p.sfr + 1.0) * h_z;
}

// ============================================================
// compute_system: Compute all three solutions + DPM creation
// ============================================================
MultiAstroResult UQFFMultiAstroCore::compute_system(const MultiAstroParams& p) {
    MultiAstroResult result;
    result.system_name   = p.name;
    result.F_compressed  = F_compressed(p);
    result.F_resonance   = F_resonance(p);
    result.F_buoyancy    = F_buoyancy(p);
    result.F_DPM_creation = F_DPM_creation(p);
    return result;
}

// ============================================================
// compute_all_systems: Batch compute all 11 systems
// ============================================================
std::vector<MultiAstroResult> UQFFMultiAstroCore::compute_all_systems() {
    std::vector<MultiAstroResult> results;
    auto systems = {
        create_NGC4826_system(), create_NGC1805_system(), create_NGC6307_system(),
        create_NGC7027_system(), create_Cassini_Encke_system(), create_Cassini_Div_system(),
        create_Cassini_Max_system(), create_ESO391_12_system(), create_M57_system(),
        create_LMC_system(), create_ESO510_G13_system()
    };
    for (auto& sys : systems) {
        results.push_back(compute_system(sys.get_params()));
    }
    return results;
}

// ============================================================
// UQFFMultiAstroSystem
// ============================================================
UQFFMultiAstroSystem::UQFFMultiAstroSystem(const MultiAstroParams& params) : params_(params) {}

MultiAstroResult UQFFMultiAstroSystem::compute(UQFFMultiAstroCore& core) const {
    return core.compute_system(params_);
}

// ============================================================
// Factory functions — Parameters from DeepSearch-validated data
// ============================================================

UQFFMultiAstroSystem create_NGC4826_system() {
    MultiAstroParams p = {3.31e20, 0.5, 1e-5, 0.0014, 3.15e17, MULTI_NGC4826, "NGC4826 (Black Eye Galaxy)"};
    return UQFFMultiAstroSystem(p);
}

UQFFMultiAstroSystem create_NGC1805_system() {
    MultiAstroParams p = {3.0e17, 0.2, 1e-5, 0.0005, 9.46e14, MULTI_NGC1805, "NGC1805 (Star Cluster)"};
    return UQFFMultiAstroSystem(p);
}

UQFFMultiAstroSystem create_NGC6307_system() {
    MultiAstroParams p = {9.46e15, 0.1, 1e-5, 0.0007, 9.46e13, MULTI_NGC6307, "NGC6307 (Lenticular Galaxy)"};
    return UQFFMultiAstroSystem(p);
}

UQFFMultiAstroSystem create_NGC7027_system() {
    MultiAstroParams p = {9.46e15, 0.1, 1e-5, 0.001, 3.15e10, MULTI_NGC7027, "NGC7027 (Planetary Nebula)"};
    return UQFFMultiAstroSystem(p);
}

UQFFMultiAstroSystem create_Cassini_Encke_system() {
    // Encke Gap: 133,590 km from Saturn center
    MultiAstroParams p = {1.3359e8, 0.0, 1e-7, 0.0, 3.156e7, MULTI_CASSINI_ENCKE, "Cassini Encke Gap"};
    return UQFFMultiAstroSystem(p);
}

UQFFMultiAstroSystem create_Cassini_Div_system() {
    // Cassini Division: 117,500-122,200 km — gap between A and B rings
    MultiAstroParams p = {1.2e8, 0.0, 1e-7, 0.0, 3.156e7, MULTI_CASSINI_DIV, "Cassini Division"};
    return UQFFMultiAstroSystem(p);  // Fixed: was UQFFAastroSystem(p) (typo corrected)
}

UQFFMultiAstroSystem create_Cassini_Max_system() {
    // Maxwell Gap: 87,500 km from Saturn center
    MultiAstroParams p = {8.75e7, 0.0, 1e-7, 0.0, 3.156e7, MULTI_CASSINI_MAX, "Cassini Maxwell Gap"};
    return UQFFMultiAstroSystem(p);
}

UQFFMultiAstroSystem create_ESO391_12_system() {
    MultiAstroParams p = {4.73e20, 0.2, 1e-5, 0.0067, 9.46e13, MULTI_ESO391_12, "ESO391-12 (Galaxy)"};
    return UQFFMultiAstroSystem(p);
}

UQFFMultiAstroSystem create_M57_system() {
    // Messier 57 / Ring Nebula, distance ~2.3 kly
    MultiAstroParams p = {1.89e16, 0.0, 1e-5, 0.0008, 1.58e11, MULTI_M57, "M57 (Ring Nebula)"};
    return UQFFMultiAstroSystem(p);
}

UQFFMultiAstroSystem create_LMC_system() {
    // Large Magellanic Cloud, distance ~160 kly
    MultiAstroParams p = {1.32e20, 0.4, 1e-5, 0.0005, 4.1e17, MULTI_LMC, "Large Magellanic Cloud"};
    return UQFFMultiAstroSystem(p);
}

UQFFMultiAstroSystem create_ESO510_G13_system() {
    // ESO510-G13: warped edge-on spiral, distance ~490 Mly
    MultiAstroParams p = {9.46e20, 1.0, 1e-5, 0.011, 9.46e13, MULTI_ESO510_G13, "ESO510-G13 (Warped Spiral)"};
    return UQFFMultiAstroSystem(p);
}
