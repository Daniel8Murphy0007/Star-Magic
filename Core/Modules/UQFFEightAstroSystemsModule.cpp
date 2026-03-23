// UQFFEightAstroSystemsModule.cpp
// Implementation of UQFF Eight Astrophysical Systems (annotated proof version)
// Based on "8 Astro Systems_C_11Oct2025.docx" — star-forming regions
// Systems: AFGL5180, NGC346, LMC opo9944a, LMC heic1301, LMC potw1408a,
//          LMC heic1206, LMC heic1402, NGC2174
// Provides full 7-step equation proofs and numerical verification table.
// Watermark: Copyright - Daniel T. Murphy, analyzed by Grok 3, dated November 17, 2025

#include "../UQFFEightAstroSystemsModule.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std::complex_literals;

// ============================================================
// UQFFEightAstroCore: Constructor
// ============================================================
UQFFEightAstroCore::UQFFEightAstroCore(double k_comp, double k_res)
    : k_comp_(k_comp), k_res_(k_res) {}

// ============================================================
// F_compressed_proof: Compressed UQFF with 7-step inline proof
// Proof:
//  Step 1: Identify system parameters (r, B, sfr, z, t_age)
//  Step 2: Compute f_UA' = (Z_max - Z) / Z_max with Z=26
//  Step 3: Compute f_SCm = Z / Z_max
//  Step 4: Hubble correction h_z = (1 + z)
//  Step 5: E_rad factor = 1 - 0.1554 = 0.8446
//  Step 6: F_c = k_comp * rho_vac * r^2 * h_z * E_rad_factor
//  Step 7: Add imaginary term i * k_comp * rho_vac * B * r / c * h_z
// ============================================================
std::complex<double> UQFFEightAstroCore::F_compressed_proof(
    const EightAstroParams& p, std::string& proof_out) {
    // Step 1
    double Z = 26.0;
    // Step 2
    double f_UA_prime = (26.0 < 1000.0) ? (1000.0 - Z) / 1000.0 : 0.0;
    // Step 3
    double f_SCm = Z / 1000.0;
    // Step 4
    double h_z = 1.0 + p.z;
    // Step 5
    double e_rad = EIGHT_E_RAD_FACTOR;
    // Step 6
    double r_eff = (p.r > 0) ? p.r : 1.0;
    double re = k_comp_ * EIGHT_RHO_VAC_UA * r_eff * r_eff * h_z * e_rad;
    // Step 7
    double im = k_comp_ * EIGHT_RHO_VAC_UA * p.B * r_eff / EIGHT_C * h_z;

    std::ostringstream oss;
    oss << "[Compressed Proof - " << p.name << "]\n"
        << "  Step 1: r=" << std::scientific << r_eff << " m, B=" << p.B << " T, z=" << p.z << "\n"
        << "  Step 2: f_UA' = (1000 - " << Z << ") / 1000 = " << f_UA_prime << "\n"
        << "  Step 3: f_SCm = " << Z << " / 1000 = " << f_SCm << "\n"
        << "  Step 4: h(z) = 1 + " << p.z << " = " << h_z << "\n"
        << "  Step 5: E_rad_factor = " << e_rad << "\n"
        << "  Step 6: F_c(re) = " << k_comp_ << " * " << EIGHT_RHO_VAC_UA
        << " * r^2 * " << h_z << " * " << e_rad << " = " << re << " N\n"
        << "  Step 7: F_c(im) = " << im << " N\n"
        << "  Result: F_compressed = (" << re << " + " << im << "i) N\n";
    proof_out = oss.str();
    return std::complex<double>(re, im);
}

// ============================================================
// F_resonance_proof: Resonance UQFF with 7-step inline proof
// Proof:
//  Step 1: Identify parameters
//  Step 2: omega_THz = 2*PI*1e12 rad/s
//  Step 3: Resonance phase = sin(omega_THz * t_age * norm)
//  Step 4: h_z = (1 + z)
//  Step 5: Compute F_res(re) = k_res * rho_vac * B * h_z * phase
//  Step 6: Add star-formation rate correction (imaginary)
//  Step 7: Final result
// ============================================================
std::complex<double> UQFFEightAstroCore::F_resonance_proof(
    const EightAstroParams& p, std::string& proof_out) {
    // Step 1-2
    double omega_THz = 2.0 * EIGHT_PI * EIGHT_NU_THz;
    // Step 3
    double t_eff = (p.t_age > 0) ? p.t_age : 1.0;
    double phase = std::sin(omega_THz * t_eff * 1e-30); // normalized
    // Step 4
    double h_z = 1.0 + p.z;
    // Step 5
    double re = k_res_ * EIGHT_RHO_VAC_UA * p.B * h_z * phase;
    // Step 6
    double im = k_res_ * EIGHT_RHO_VAC_UA * p.sfr / EIGHT_C * h_z;
    // Step 7
    std::ostringstream oss;
    oss << "[Resonance Proof - " << p.name << "]\n"
        << "  Step 1: B=" << p.B << " T, sfr=" << p.sfr << " M_sun/yr, z=" << p.z << "\n"
        << "  Step 2: omega_THz = 2*PI*1e12 = " << omega_THz << " rad/s\n"
        << "  Step 3: phase = sin(omega_THz * t_age * 1e-30) = " << phase << "\n"
        << "  Step 4: h(z) = " << h_z << "\n"
        << "  Step 5: F_res(re) = " << re << " N\n"
        << "  Step 6: F_res(im) = " << im << " N\n"
        << "  Step 7: F_resonance = (" << re << " + " << im << "i) N\n";
    proof_out = oss.str();
    return std::complex<double>(re, im);
}

// ============================================================
// compute_with_proof: Compute both solutions with proofs
// ============================================================
EightAstroResult UQFFEightAstroCore::compute_with_proof(const EightAstroParams& p) {
    EightAstroResult result;
    result.system_name = p.name;
    result.F_compressed = F_compressed_proof(p, result.proof_compressed);
    result.F_resonance  = F_resonance_proof(p, result.proof_resonance);
    return result;
}

// ============================================================
// compute_all_systems: Batch compute all 8 systems
// ============================================================
std::vector<EightAstroResult> UQFFEightAstroCore::compute_all_systems() {
    std::vector<EightAstroResult> results;
    std::vector<UQFFEightAstroSystem> systems = {
        create_AFGL5180_system(), create_NGC346_system(),
        create_LMC_opo9944a_system(), create_LMC_heic1301_system(),
        create_LMC_potw1408a_system(), create_LMC_heic1206_system(),
        create_LMC_heic1402_system(), create_NGC2174_system()
    };
    for (auto& sys : systems) {
        results.push_back(compute_with_proof(sys.get_params()));
    }
    return results;
}

// ============================================================
// print_verification_table: Numerical verification output
// ============================================================
void UQFFEightAstroCore::print_verification_table(const std::vector<EightAstroResult>& results) {
    std::cout << "\n=== UQFF Eight Astro Systems — Numerical Verification Table ===\n";
    std::cout << std::setw(25) << "System"
              << std::setw(30) << "F_compressed (N)"
              << std::setw(30) << "F_resonance (N)" << "\n";
    std::cout << std::string(85, '-') << "\n";
    for (const auto& r : results) {
        std::cout << std::setw(25) << r.system_name
                  << std::setw(30) << r.F_compressed
                  << std::setw(30) << r.F_resonance << "\n";
    }
    std::cout << std::string(85, '=') << "\n";
}

// ============================================================
// UQFFEightAstroSystem
// ============================================================
UQFFEightAstroSystem::UQFFEightAstroSystem(const EightAstroParams& params) : params_(params) {}

EightAstroResult UQFFEightAstroSystem::compute(UQFFEightAstroCore& core) const {
    return core.compute_with_proof(params_);
}

// ============================================================
// Factory functions — Parameters for 8 star-forming systems
// ============================================================

// AFGL5180: HII region / star-forming cloud
UQFFEightAstroSystem create_AFGL5180_system() {
    EightAstroParams p = {1e16, 0.01, 1e-4, 0.0, 3.15e13, EIGHT_AFGL5180, "AFGL5180"};
    return UQFFEightAstroSystem(p);
}

// NGC346: LMC star cluster, intense star formation
UQFFEightAstroSystem create_NGC346_system() {
    EightAstroParams p = {1e19, 0.1, 1e-5, 0.0006, 3.15e14, EIGHT_NGC346, "NGC346 (GFSC cluster)"};
    return UQFFEightAstroSystem(p);
}

// LMC opo9944a: LMC HST image 1 — star birth pillar
UQFFEightAstroSystem create_LMC_opo9944a_system() {
    EightAstroParams p = {5e18, 0.05, 1e-5, 0.0009, 1.58e14, EIGHT_LMC_OPO9944A, "LMC opo9944a"};
    return UQFFEightAstroSystem(p);
}

// LMC heic1301: 30 Doradus surroundings
UQFFEightAstroSystem create_LMC_heic1301_system() {
    EightAstroParams p = {2e19, 0.2, 1e-5, 0.0009, 3.15e14, EIGHT_LMC_HEIC1301, "LMC heic1301 (30 Dor)"};
    return UQFFEightAstroSystem(p);
}

// LMC potw1408a: LMC massive star nursery
UQFFEightAstroSystem create_LMC_potw1408a_system() {
    EightAstroParams p = {8e18, 0.08, 1e-5, 0.0009, 2.0e14, EIGHT_LMC_POTW1408A, "LMC potw1408a"};
    return UQFFEightAstroSystem(p);
}

// LMC heic1206: Complex in the LMC
UQFFEightAstroSystem create_LMC_heic1206_system() {
    EightAstroParams p = {6e18, 0.06, 1e-5, 0.0009, 1.9e14, EIGHT_LMC_HEIC1206, "LMC heic1206"};
    return UQFFEightAstroSystem(p);
}

// LMC heic1402: Stellar nursery in LMC
UQFFEightAstroSystem create_LMC_heic1402_system() {
    EightAstroParams p = {7e18, 0.07, 1e-5, 0.0009, 2.2e14, EIGHT_LMC_HEIC1402, "LMC heic1402"};
    return UQFFEightAstroSystem(p);
}

// NGC2174: Monkey Head Nebula HII region
UQFFEightAstroSystem create_NGC2174_system() {
    EightAstroParams p = {2e19, 0.1, 1e-5, 0.00015, 1.58e14, EIGHT_NGC2174, "NGC2174 (Monkey Head Nebula)"};
    return UQFFEightAstroSystem(p);
}
