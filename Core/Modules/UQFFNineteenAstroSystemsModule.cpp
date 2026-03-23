// UQFFNineteenAstroSystemsModule.cpp
// Implementation of UQFF 19-System 26D Polynomial Framework
// Based on "19 Astro Systems_11Oct2025.docx"
// Implements g(r,t) = Sum_{i=1}^{26} [E_DPM,i / r_i^2] * f_TRZ_i * f_Um_i * H(z) * (1 - E_rad)
// 19 systems: NGC2264, UGC10214, NGC4676, Red Spider, NGC3372, AG Carinae, M42, Tarantula,
//             NGC2841, Mystic Mountain, NGC6217, Stephan's Quintet, NGC7049, Carina NGC3324,
//             M74, NGC1672, NGC5866, M82, Spirograph IC418
// Watermark: Copyright - Daniel T. Murphy, analyzed by Grok 3, dated November 17, 2025

#include "../UQFFNineteenAstroSystemsModule.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <numeric>

using namespace std::complex_literals;

// ============================================================
// UQFFNineteenAstroCore: Constructor
// ============================================================
UQFFNineteenAstroCore::UQFFNineteenAstroCore(double k_poly, double k_comp, double k_res)
    : k_poly_(k_poly), k_comp_(k_comp), k_res_(k_res) {}

// ============================================================
// build_vars: Construct 26D DPM variable set from system parameters
// omega_i(i) = H_Z_BASE * i  (dimensional frequency scaling)
// ============================================================
NineteenDPMVars26D UQFFNineteenAstroCore::build_vars(const NineteenAstroParams& p) {
    NineteenDPMVars26D v;
    for (int i = 0; i < NINETEEN_DIM; ++i) {
        double fi = static_cast<double>(i + 1); // dimensions 1..26
        double Z_i = 26.0 + fi; // atomic-like index per dimension

        v.f_UA_prime[i] = (1000.0 - Z_i) / 1000.0;
        v.f_SCm[i]      = Z_i / 1000.0;
        v.R_EB[i]       = NINETEEN_K_R * Z_i;
        v.Q_i[i]        = 1.0 / (1.0 + fi * NINETEEN_E_RAD);
        v.theta_i[i]    = NINETEEN_PI * fi / (2.0 * NINETEEN_DIM);
        // r_i: physical radius scaled by dimension index
        v.r_i[i]        = p.r * (1.0 + fi / static_cast<double>(NINETEEN_DIM));
        // f_TRZ_i: THz resonance per dimension — Einstein Boson Bridge per state
        double d_i = v.r_i[i] * 0.01;
        v.f_TRZ_i[i]   = std::exp(-NINETEEN_NU_THz * d_i / NINETEEN_C) /
                          (1.0 + NINETEEN_NU_THz * 1e-22);
        // f_Um_i: Magnetism per dimension
        v.f_Um_i[i]    = v.f_UA_prime[i] * v.f_SCm[i] * p.B *
                          std::sin(NINETEEN_PI * fi / NINETEEN_DIM);
    }
    return v;
}

// ============================================================
// compute_g_poly: 26D gravity polynomial sum
// g = k_poly * Sum_{i=1}^{26} [E_DPM,i / r_i^2] * f_TRZ_i * f_Um_i * H(z) * (1 - E_rad)
// E_DPM,i = Q_i * rho_vac * c^2 / Z_i  (DPM energy per  dimension)
// ============================================================
std::complex<double> UQFFNineteenAstroCore::compute_g_poly(
    const NineteenAstroParams& p, const NineteenDPMVars26D& vars) {
    double h_z = 1.0 + p.z;
    double e_rad_factor = 1.0 - NINETEEN_E_RAD;
    std::complex<double> gsum(0.0, 0.0);
    for (int i = 0; i < NINETEEN_DIM; ++i) {
        double fi = static_cast<double>(i + 1);
        double Z_i = 26.0 + fi;
        double E_DPM_i = vars.Q_i[i] * NINETEEN_RHO_VAC_UA * NINETEEN_C * NINETEEN_C / Z_i;
        double r_i = (vars.r_i[i] > 0) ? vars.r_i[i] : 1.0;
        double term_re = k_poly_ * E_DPM_i / (r_i * r_i)
                       * vars.f_TRZ_i[i] * vars.f_Um_i[i] * h_z * e_rad_factor;
        // Imaginary: quantum phase per dimension omega_i * t / c
        double term_im = term_re * omega_i(i + 1) * p.t_age / (NINETEEN_C * 1e15);
        gsum += std::complex<double>(term_re, term_im);
    }
    return gsum;
}

// ============================================================
// F_compressed: Compressed UQFF (26D-aware)
// ============================================================
std::complex<double> UQFFNineteenAstroCore::F_compressed(const NineteenAstroParams& p) {
    double r_eff = (p.r > 0) ? p.r : 1.0;
    double h_z = 1.0 + p.z;
    double e_rad = 1.0 - NINETEEN_E_RAD;
    double re = k_comp_ * NINETEEN_RHO_VAC_UA * r_eff * r_eff * h_z * e_rad;
    double im = k_comp_ * NINETEEN_RHO_VAC_UA * p.B * r_eff / NINETEEN_C * h_z;
    return std::complex<double>(re, im);
}

// ============================================================
// F_resonance: Resonance UQFF (26D-aware with omega_i summation)
// ============================================================
std::complex<double> UQFFNineteenAstroCore::F_resonance(const NineteenAstroParams& p) {
    double h_z = 1.0 + p.z;
    double omega_sum = 0.0;
    for (int i = 1; i <= NINETEEN_DIM; ++i) omega_sum += omega_i(i);
    double phase = std::sin(omega_sum * p.t_age * 1e-40); // normalized
    double re = k_res_ * NINETEEN_RHO_VAC_UA * p.B * h_z * phase;
    double im = k_res_ * NINETEEN_RHO_VAC_UA * p.sfr / NINETEEN_C * h_z;
    return std::complex<double>(re, im);
}

// ============================================================
// evaluate_26D_poly: Horner's method for 26D Taylor polynomial
// P(x) = a[0] + a[1]*x + ... + a[25]*x^25
// ============================================================
double UQFFNineteenAstroCore::evaluate_26D_poly(
    const std::array<double, 26>& coeffs, double x) {
    double result = coeffs[25];
    for (int i = 24; i >= 0; --i) {
        result = result * x + coeffs[i];
    }
    return result;
}

// ============================================================
// compute_system: Full 26D computation for a single system
// ============================================================
NineteenAstroResult UQFFNineteenAstroCore::compute_system(const NineteenAstroParams& p) {
    NineteenAstroResult result;
    result.system_name = p.name;
    NineteenDPMVars26D vars = build_vars(p);
    result.g_poly       = compute_g_poly(p, vars);
    result.F_compressed = F_compressed(p);
    result.F_resonance  = F_resonance(p);
    // Build Taylor coefficients array from f_UA_prime (26D polynomial basis)
    for (int i = 0; i < NINETEEN_DIM; ++i) {
        result.taylor_coeffs[i] = vars.f_UA_prime[i] * vars.Q_i[i];
    }
    return result;
}

// ============================================================
// compute_all_systems: Batch compute all 19 systems
// ============================================================
std::vector<NineteenAstroResult> UQFFNineteenAstroCore::compute_all_systems() {
    std::vector<NineteenAstroResult> results;
    std::vector<UQFFNineteenAstroSystem> systems = {
        create_NGC2264_system(), create_UGC10214_system(), create_NGC4676_system(),
        create_RedSpider_system(), create_NGC3372_system(), create_AGCarinae_system(),
        create_M42_system(), create_Tarantula_system(), create_NGC2841_system(),
        create_MysticMountain_system(), create_NGC6217_system(), create_StephansQuintet_system(),
        create_NGC7049_system(), create_CarinaNGC3324_system(), create_M74_system(),
        create_NGC1672_system(), create_NGC5866_system(), create_M82_system(),
        create_Spirograph_system()
    };
    for (auto& sys : systems) {
        results.push_back(compute_system(sys.get_params()));
    }
    return results;
}

// ============================================================
// UQFFNineteenAstroSystem
// ============================================================
UQFFNineteenAstroSystem::UQFFNineteenAstroSystem(const NineteenAstroParams& params) : params_(params) {}

NineteenAstroResult UQFFNineteenAstroSystem::compute(UQFFNineteenAstroCore& core) const {
    return core.compute_system(params_);
}

// ============================================================
// Factory functions — 19 systems
// AstroParams: {M_0, r, sfr, B, z, t_age, type, name}
// ============================================================

UQFFNineteenAstroSystem create_NGC2264_system() {
    NineteenAstroParams p = {1.989e36, 2e19, 0.5, 1e-5, 0.0006, 3e14, N19_NGC2264, "NGC2264 (Cone Nebula/Cluster)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_UGC10214_system() {
    NineteenAstroParams p = {1.989e41, 1.3e21, 1.0, 1e-5, 0.028, 9.46e13, N19_UGC10214, "UGC10214 (Tadpole Galaxy)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_NGC4676_system() {
    NineteenAstroParams p = {3.978e41, 3e20, 10.0, 1e-4, 0.022, 9.46e13, N19_NGC4676, "NGC4676 (Mice Galaxies)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_RedSpider_system() {
    NineteenAstroParams p = {1.989e30, 1e16, 0.0, 1e-5, 0.0013, 3.15e13, N19_RED_SPIDER, "Red Spider Nebula"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_NGC3372_system() {
    NineteenAstroParams p = {1.989e35, 2e17, 2.0, 1e-5, 0.0025, 3.15e13, N19_NGC3372, "NGC3372 (Eta Carinae Nebula)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_AGCarinae_system() {
    NineteenAstroParams p = {3.978e31, 1e16, 0.0, 1e-5, 0.002, 3.15e13, N19_AG_CARINAE, "AG Carinae Nebula"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_M42_system() {
    NineteenAstroParams p = {3.978e33, 2e16, 0.3, 1e-5, 0.0004, 3.15e13, N19_M42, "M42 (Orion Nebula)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_Tarantula_system() {
    NineteenAstroParams p = {1.989e35, 3e17, 5.0, 1e-4, 0.0005, 3.15e13, N19_TARANTULA, "Tarantula Nebula (LMC)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_NGC2841_system() {
    NineteenAstroParams p = {1.989e41, 5e20, 0.5, 1e-5, 0.0031, 9.46e13, N19_NGC2841, "NGC2841 (Spiral Galaxy)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_MysticMountain_system() {
    NineteenAstroParams p = {1.989e32, 1e16, 0.1, 1e-5, 0.0025, 3.15e13, N19_MYSTIC_MOUNTAIN, "Mystic Mountain (Carina)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_NGC6217_system() {
    NineteenAstroParams p = {1.989e41, 3e20, 1.0, 1e-5, 0.0045, 9.46e13, N19_NGC6217, "NGC6217 (Barred Spiral)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_StephansQuintet_system() {
    NineteenAstroParams p = {9.945e41, 1e21, 10.0, 1e-4, 0.022, 9.46e13, N19_STEPHANS_QUINTET, "Stephan's Quintet"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_NGC7049_system() {
    NineteenAstroParams p = {1.989e41, 5e20, 0.2, 1e-5, 0.0067, 9.46e13, N19_NGC7049, "NGC7049 (Lenticular)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_CarinaNGC3324_system() {
    NineteenAstroParams p = {1.989e35, 2e17, 2.0, 1e-5, 0.0025, 3.15e13, N19_CARINA_NGC3324, "Carina NGC3324"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_M74_system() {
    NineteenAstroParams p = {1.989e41, 5e20, 1.0, 1e-5, 0.0022, 9.46e13, N19_M74, "M74 (Phantom Galaxy)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_NGC1672_system() {
    NineteenAstroParams p = {1.989e41, 3e20, 2.0, 1e-5, 0.004, 9.46e13, N19_NGC1672, "NGC1672 (Barred Spiral)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_NGC5866_system() {
    NineteenAstroParams p = {1.989e41, 3e20, 0.3, 1e-5, 0.0029, 9.46e13, N19_NGC5866, "NGC5866 (Spindle Galaxy)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_M82_system() {
    NineteenAstroParams p = {1.989e40, 2e20, 10.0, 1e-4, 0.0008, 9.46e13, N19_M82, "M82 (Cigar Galaxy)"};
    return UQFFNineteenAstroSystem(p);
}

UQFFNineteenAstroSystem create_Spirograph_system() {
    NineteenAstroParams p = {1.989e30, 1e16, 0.0, 1e-5, 0.0007, 3.15e13, N19_SPIROGRAPH_IC418, "Spirograph IC418"};
    return UQFFNineteenAstroSystem(p);
}
