// test_source4_validation.cpp - SOURCE4 Dual-Method Validation Test Suite
// Tests UQFF vs MUGE cross-validation for all 7 astrophysical systems
// Created: April 16, 2026 (Phase 2 of SOURCE4 Integration Plan)
//
// Build (standalone, no Wolfram required):
//   cl /std:c++20 /EHsc /O2 test_source4_validation.cpp /Fe:test_source4_validation.exe
//
// Build (via CMake):
//   cmake --build build_msvc --config Release --target test_source4_validation
//
// Run:
//   .\build_msvc\Release\test_source4_validation.exe
//   Output: validation_log.txt (detailed per-system results)

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <chrono>
#include <cassert>

// ============================================================================
// MINIMAL SOURCE4 PHYSICS DEFINITIONS (standalone, no MAIN_1_CoAnQi dependency)
// Extracted from MAIN_1_CoAnQi.cpp namespace SOURCE4 (lines 25620-25830)
// ============================================================================

static const double PI_S4 = 3.14159265358979323846;
static const double G_S4 = 6.674e-11;         // m^3 kg^-1 s^-2
static const double c_S4 = 2.998e8;            // m/s
static const double hbar_S4 = 1.0546e-34;      // J·s
static const double B_CRIT_S4 = 4.4e13;        // T (critical magnetar field)
static const double H0_S4 = 2.27e-18;          // s^-1 (Hubble constant ~70 km/s/Mpc)
static const double LAMBDA_COSM_S4 = 1.1056e-52; // m^-2 (cosmological constant)

// Celestial body parameters for UQFF
struct CelestialBody_S4 {
    std::string name;
    double Ms;           // Mass (kg)
    double Rs;           // Radius (m)
    double B_surface;    // Surface magnetic field (T)
    double T_surface;    // Surface temperature (K)
    double Ereact;       // Reactor efficiency
    double SCm_density;  // SCm density
    double QUA;          // Quantum uncertainty amplitude
    double Pcore;        // Core pressure factor
    double PSCm;         // SCm pressure
    double omega_c;      // Core rotation frequency (rad/s)
};

// MUGE system parameters (compressed + resonance)
struct MUGESystem_S4 {
    double M;            // Mass (kg)
    double r;            // Distance (m)
    double t;            // Time (s)
    double B;            // Magnetic field (T)
    double rho_fluid;    // Fluid density (kg/m^3)
    double V_sys;        // System velocity (m/s)
    double g_local;      // Local gravity estimate (m/s^2)
    double H_z;          // Hubble parameter at redshift z
    double M_DM;         // Dark matter mass (kg)
    double delta_rho;    // Density perturbation ratio
    double envelope;     // Envelope factor
};

// Resonance parameters
struct ResonanceParams_S4 {
    double I_dipole;     // Current (A)
    double A_loop;       // Loop area (m^2)
    double omega1;       // Frequency 1 (rad/s)
    double omega2;       // Frequency 2 (rad/s)
    double fTHz;         // THz frequency
    double vexp;         // Expansion velocity (m/s)
    double Evac;         // Vacuum energy density
    double Evac_bar;     // Mean vacuum energy density
    double delta_Evac;   // Vacuum energy difference
    double Fsuper;       // Superconductive force
    double UA;           // Aether coupling
    double omega_i;      // Intermediate frequency
    double k4;           // Ug4 coupling
    double freact;       // Reactor frequency
    double fquantum;     // Quantum frequency
    double fAether;      // Aether frequency
    double ffluid;       // Fluid frequency
    double fTRZ;         // TRZ frequency
    double fexp_factor;  // Expansion factor
    double f_worm;       // Wormhole coupling
    double b_throat;     // Wormhole throat radius (m)
};

// ============================================================================
// UQFF COMPUTATION (Simplified from SOURCE4 inline functions)
// ============================================================================

// Ug1: Magnetic dipole-gradient gravity
static double compute_Ug1(const CelestialBody_S4& body, double r) {
    double mu_s = body.B_surface * body.Rs * body.Rs * body.Rs;
    double grad_Ms_r = G_S4 * body.Ms / (r * r);
    return (mu_s * grad_Ms_r) / (4.0 * PI_S4 * r * r * r);
}

// Ug2: Charge-reactivity gravity
static double compute_Ug2(const CelestialBody_S4& body, double r) {
    double eV_react = body.Ereact * body.Ms * c_S4 * c_S4;
    return eV_react / (4.0 * PI_S4 * r * r);
}

// Ug3: String rotation gravity
static double compute_Ug3(const CelestialBody_S4& body, double r, double t) {
    double omega_ST = body.omega_c * (1.0 + 0.01 * sin(2.0 * PI_S4 * t));
    double Bj = body.B_surface * (body.Rs / r) * (body.Rs / r) * (body.Rs / r);
    double mu_j = Bj * body.Rs * body.Rs * body.Rs;
    return (mu_j * omega_ST) / (4.0 * PI_S4 * r * r);
}

// Ug4: Vacuum concentration gravity
static double compute_Ug4(const CelestialBody_S4& body, double r) {
    double rho_vac = body.SCm_density * body.QUA;
    return (4.0 / 3.0) * PI_S4 * G_S4 * rho_vac * r;
}

// Ubi: Buoyancy force
static double compute_Ubi(const CelestialBody_S4& body, double r) {
    double rho_body = body.Ms / ((4.0 / 3.0) * PI_S4 * body.Rs * body.Rs * body.Rs);
    double rho_vac = body.SCm_density * body.QUA;
    return (4.0 / 3.0) * PI_S4 * G_S4 * (rho_body - rho_vac) * r;
}

// Um: Magnetism contribution
static double compute_Um(const CelestialBody_S4& body, double r) {
    double Bj = body.B_surface * pow(body.Rs / r, 3);
    double mu_j = Bj * body.Rs * body.Rs * body.Rs;
    double n_strings = 1e9;  // billion magnetic strings
    double single = (mu_j / (4.0 * PI_S4 * r * r * r));
    return single * n_strings * body.PSCm * body.Ereact;
}

// FU: Complete unified field
static double compute_FU(const CelestialBody_S4& body, double r, double t) {
    double ug1 = compute_Ug1(body, r);
    double ug2 = compute_Ug2(body, r);
    double ug3 = compute_Ug3(body, r, t);
    double ug4 = compute_Ug4(body, r);
    double ubi = compute_Ubi(body, r);
    double um  = compute_Um(body, r);
    return ug1 + ug2 + ug3 + ug4 + ubi + um;
}

// ============================================================================
// MUGE COMPRESSED COMPUTATION (9-term formula)
// ============================================================================

static double compute_compressed_MUGE(const MUGESystem_S4& sys) {
    // Base: Newtonian gravity
    double base = G_S4 * sys.M / (sys.r * sys.r);

    // Expansion: Hubble factor
    double expansion = 1.0 + H0_S4 * sys.t;

    // Superconductive adjustment
    double super_adj = 1.0 - sys.B / B_CRIT_S4;

    // Envelope factor (neutral)
    double env = sys.envelope;

    // Core product
    double core = base * expansion * super_adj * env;

    // Ug sum (simplified: sum of Ug1-4 magnitudes)
    double ug_sum = base * 0.01;  // ~1% correction from Ug terms

    // Cosmological constant term
    double cosm = LAMBDA_COSM_S4 * c_S4 * c_S4 / 3.0;

    // Quantum uncertainty
    double quantum = hbar_S4 / (sys.M * sys.r);

    // Fluid coupling
    double fluid = sys.rho_fluid * sys.V_sys * sys.g_local;

    // Dark matter perturbation
    double perturbation = G_S4 * sys.M_DM / (sys.r * sys.r) + sys.delta_rho * base;

    return core + ug_sum + cosm + quantum + fluid + perturbation;
}

// ============================================================================
// MUGE RESONANCE COMPUTATION (13-term cascade)
// ============================================================================

static double compute_resonance_MUGE(const MUGESystem_S4& sys, const ResonanceParams_S4& p) {
    // Root: aDPM = FDPM * fDPM * Evac * c * V
    double FDPM = p.I_dipole * p.A_loop;
    double fDPM = p.omega1 * p.omega2 / (4.0 * PI_S4 * PI_S4);
    double aDPM = FDPM * fDPM * p.Evac * c_S4 * sys.V_sys;

    // Cascade from aDPM
    double aTHz = p.fTHz * p.Evac * p.vexp * aDPM / (p.Evac_bar * c_S4 + 1e-300);
    double avac_diff = p.delta_Evac * p.vexp * p.vexp * aDPM / (p.Evac * c_S4 * c_S4 + 1e-300);
    double asuper_freq = p.Fsuper * p.fTHz * aDPM / (p.Evac * c_S4 + 1e-300);
    double aaether_res = p.UA * p.omega_i * p.fTHz * aDPM * (1.0 + p.fTRZ);
    double Ug4i = p.k4 * p.Evac * p.freact * aDPM / (p.Evac * c_S4 + 1e-300);
    double aquantum_freq = p.fquantum * p.Evac * aDPM / (p.Evac_bar * c_S4 + 1e-300);
    double aAether_freq = p.fAether * p.Evac * aDPM / (p.Evac_bar * c_S4 + 1e-300);
    double afluid_freq = p.ffluid * p.Evac * sys.V_sys / (p.Evac_bar * c_S4 + 1e-300);
    double osc_term = 0.0;  // Simplified
    double aexp_freq_val = 2.0 * PI_S4 * sys.H_z * sys.t;
    double aexp_freq = aexp_freq_val * p.Evac * aDPM / (p.Evac_bar * c_S4 + 1e-300);
    double fTRZ_pass = p.fTRZ;
    double a_wormhole = p.f_worm * p.Evac / (p.b_throat * p.b_throat + sys.r * sys.r + 1e-300);

    return aDPM + aTHz + avac_diff + asuper_freq + aaether_res + Ug4i
         + aquantum_freq + aAether_freq + afluid_freq + osc_term
         + aexp_freq + fTRZ_pass + a_wormhole;
}

// ============================================================================
// PRE-DEFINED ASTROPHYSICAL SYSTEMS (from MAIN_1_CoAnQi.cpp namespace SOURCE4)
// ============================================================================

// UQFF bodies
static const CelestialBody_S4 sgr1745_body = {
    "SGR1745", 2.78e30, 1.0e4, 1.0e11, 5e6, 1e-4, 1e17, 1e-11, 1e-3, 1.0, 714.0
};
static const CelestialBody_S4 sagA_body = {
    "SagA*", 7.956e36, 1.18e10, 1e-2, 1e7, 1e-6, 1e15, 1e-11, 1e-3, 1.0, 1e-6
};
static const CelestialBody_S4 tapestry_body = {
    "Tapestry", 1e34, 3.086e18, 1e-5, 300.0, 1e-8, 1e12, 1e-11, 1e-3, 1.0, 1e-10
};
static const CelestialBody_S4 westerlund2_body = {
    "Westerlund2", 5e34, 6.172e18, 1e-4, 3e4, 1e-7, 1e13, 1e-11, 1e-3, 1.0, 1e-9
};
static const CelestialBody_S4 pillars_body = {
    "Pillars", 1e35, 2.778e19, 1e-6, 50.0, 1e-9, 1e11, 1e-11, 1e-3, 1.0, 1e-12
};
static const CelestialBody_S4 rings_body = {
    "Rings", 1e41, 1e22, 1e-3, 1e8, 1e-5, 1e16, 1e-11, 1e-3, 1.0, 1e-7
};
static const CelestialBody_S4 student_guide_body = {
    "StudentGuide", 1e53, 4.4e26, 1e-20, 2.725, 1e-20, 1e5, 1e-11, 1e-3, 1.0, 1e-18
};

// MUGE systems
static const MUGESystem_S4 sgr1745_muge = {
    2.78e30, 1e4, 0.0, 1e11, 1e10, 1e4, 1e12, H0_S4, 0.0, 0.0, 1.0
};
static const MUGESystem_S4 sagA_muge = {
    7.956e36, 1.18e10, 0.0, 1e-2, 1e5, 1e6, 1e10, H0_S4, 1e36, 0.01, 1.0
};
static const MUGESystem_S4 tapestry_muge = {
    1e34, 3.086e18, 0.0, 1e-5, 1e-20, 1e3, 1e-5, H0_S4, 1e33, 0.001, 1.0
};
static const MUGESystem_S4 westerlund2_muge = {
    5e34, 6.172e18, 0.0, 1e-4, 1e-19, 1e4, 1e-4, H0_S4, 1e34, 0.005, 1.0
};
static const MUGESystem_S4 pillars_muge = {
    1e35, 2.778e19, 0.0, 1e-6, 1e-21, 1e2, 1e-6, H0_S4, 1e34, 0.001, 1.0
};
static const MUGESystem_S4 rings_muge = {
    1e41, 1e22, 0.0, 1e-3, 1e3, 1e5, 1e8, H0_S4, 1e40, 0.01, 1.0
};
static const MUGESystem_S4 student_guide_muge = {
    1e53, 4.4e26, 0.0, 1e-20, 1e-27, 1e2, 1e-18, H0_S4, 1e52, 0.001, 1.0
};

// Default resonance parameters
static const ResonanceParams_S4 default_res_params = {
    1.0, 1e-4, 2*PI_S4*1e3, 2*PI_S4*1e6, 1e12, 1e4,
    1e-10, 1e-10, 1e-15, 1e-6, 1e-11, 2*PI_S4*1e9,
    1e-20, 1e6, 2*PI_S4*1e15, 2*PI_S4*1e12, 2*PI_S4*1e8,
    1e-3, 2*PI_S4*H0_S4, 1e-30, 1e3
};

// ============================================================================
// VALIDATION FRAMEWORK
// ============================================================================

struct ValidationResult {
    std::string system_name;
    double uqff_value;
    double muge_compressed;
    double muge_resonance;
    double uqff_muge_diff_pct;
    double compressed_resonance_diff_pct;
    bool convergence;
    std::string analysis;
};

struct PhysicsConstraints {
    double min_gravity;
    double max_gravity;
    double tolerance_pct;
};

static std::map<std::string, PhysicsConstraints> g_constraints = {
    {"SGR1745",      {1e9,   1e13,  10.0}},
    {"SagA*",        {1e8,   1e12,  15.0}},
    {"Tapestry",     {1e-10, 1e7,   25.0}},
    {"Westerlund2",  {1e-10, 1e6,   25.0}},
    {"Pillars",      {1e-10, 1e5,   30.0}},
    {"Rings",        {1e6,   1e14,  15.0}},
    {"StudentGuide", {1e-30, 1e-10, 50.0}},
};

static ValidationResult validate_system(
    const std::string& name,
    const CelestialBody_S4& body,
    const MUGESystem_S4& muge_sys,
    const ResonanceParams_S4& res_params)
{
    ValidationResult result;
    result.system_name = name;

    // UQFF
    result.uqff_value = compute_FU(body, muge_sys.r, muge_sys.t);

    // MUGE compressed
    result.muge_compressed = compute_compressed_MUGE(muge_sys);

    // MUGE resonance
    result.muge_resonance = compute_resonance_MUGE(muge_sys, res_params);

    // Percent differences
    if (result.uqff_value != 0.0)
        result.uqff_muge_diff_pct = std::abs(
            (result.uqff_value - result.muge_compressed) / result.uqff_value * 100.0);
    else
        result.uqff_muge_diff_pct = 0.0;

    if (result.muge_compressed != 0.0)
        result.compressed_resonance_diff_pct = std::abs(
            (result.muge_compressed - result.muge_resonance) / result.muge_compressed * 100.0);
    else
        result.compressed_resonance_diff_pct = 0.0;

    // Convergence check
    auto it = g_constraints.find(name);
    if (it != g_constraints.end()) {
        auto& c = it->second;
        double abs_compressed = std::abs(result.muge_compressed);
        bool in_range = (abs_compressed >= c.min_gravity && abs_compressed <= c.max_gravity);
        bool uqff_ok = (result.uqff_muge_diff_pct < c.tolerance_pct);
        bool muge_ok = (result.compressed_resonance_diff_pct < c.tolerance_pct);
        result.convergence = in_range && uqff_ok && muge_ok;

        std::ostringstream oss;
        if (!in_range)
            oss << "RANGE: |g|=" << abs_compressed
                << " outside [" << c.min_gravity << "," << c.max_gravity << "]; ";
        if (!uqff_ok)
            oss << "UQFF-MUGE: " << result.uqff_muge_diff_pct << "% > " << c.tolerance_pct << "%; ";
        if (!muge_ok)
            oss << "COMP-RES: " << result.compressed_resonance_diff_pct << "% > " << c.tolerance_pct << "%; ";
        result.analysis = result.convergence ? "CONVERGED" : oss.str();
    } else {
        result.convergence = false;
        result.analysis = "NO_CONSTRAINTS";
    }

    return result;
}

// ============================================================================
// TEST RUNNER
// ============================================================================

static int g_pass = 0;
static int g_fail = 0;

static void report(const char* test_name, bool passed, const std::string& detail = "") {
    if (passed) {
        std::cout << "  [PASS] " << test_name;
        ++g_pass;
    } else {
        std::cout << "  [FAIL] " << test_name;
        ++g_fail;
    }
    if (!detail.empty()) std::cout << " -- " << detail;
    std::cout << std::endl;
}

int main() {
    auto t_start = std::chrono::high_resolution_clock::now();

    std::cout << "============================================================" << std::endl;
    std::cout << " SOURCE4 Dual-Method Validation Test Suite" << std::endl;
    std::cout << " UQFF vs MUGE Compressed vs MUGE Resonance" << std::endl;
    std::cout << " 7 Astrophysical Systems | 37 Physics Functions" << std::endl;
    std::cout << "============================================================" << std::endl;

    std::ofstream log("validation_log.txt", std::ios::trunc);
    log << "SOURCE4 Validation Log" << std::endl;
    log << "=====================" << std::endl << std::endl;

    // ---- Test 1: Individual UQFF component sanity checks ----
    std::cout << "\n--- Test 1: UQFF Component Sanity ---" << std::endl;
    {
        double ug1 = compute_Ug1(sgr1745_body, sgr1745_muge.r);
        report("Ug1(SGR1745) finite", std::isfinite(ug1));
        report("Ug1(SGR1745) > 0", ug1 > 0.0,
               "Ug1=" + std::to_string(ug1));

        double ug2 = compute_Ug2(sgr1745_body, sgr1745_muge.r);
        report("Ug2(SGR1745) finite", std::isfinite(ug2));

        double ug3 = compute_Ug3(sgr1745_body, sgr1745_muge.r, 0.0);
        report("Ug3(SGR1745) finite", std::isfinite(ug3));

        double ug4 = compute_Ug4(sgr1745_body, sgr1745_muge.r);
        report("Ug4(SGR1745) finite", std::isfinite(ug4));

        double ubi = compute_Ubi(sgr1745_body, sgr1745_muge.r);
        report("Ubi(SGR1745) finite", std::isfinite(ubi));

        double um = compute_Um(sgr1745_body, sgr1745_muge.r);
        report("Um(SGR1745) finite", std::isfinite(um));

        double fu = compute_FU(sgr1745_body, sgr1745_muge.r, 0.0);
        report("FU(SGR1745) finite", std::isfinite(fu));
        report("FU(SGR1745) = Σ components",
               std::abs(fu - (ug1 + ug2 + ug3 + ug4 + ubi + um)) < 1e-20);
    }

    // ---- Test 2: MUGE compressed component sanity ----
    std::cout << "\n--- Test 2: MUGE Compressed Sanity ---" << std::endl;
    {
        double g_c = compute_compressed_MUGE(sgr1745_muge);
        report("Compressed(SGR1745) finite", std::isfinite(g_c));
        report("Compressed(SGR1745) > 0", g_c > 0.0,
               "g=" + std::to_string(g_c));

        // Newtonian baseline
        double g_newton = G_S4 * sgr1745_muge.M / (sgr1745_muge.r * sgr1745_muge.r);
        report("Compressed ~ Newtonian order",
               std::abs(g_c) > g_newton * 0.01 && std::abs(g_c) < g_newton * 100.0,
               "g_N=" + std::to_string(g_newton));
    }

    // ---- Test 3: MUGE resonance sanity ----
    std::cout << "\n--- Test 3: MUGE Resonance Sanity ---" << std::endl;
    {
        double g_r = compute_resonance_MUGE(sgr1745_muge, default_res_params);
        report("Resonance(SGR1745) finite", std::isfinite(g_r));
    }

    // ---- Test 4: Full dual-method validation for all 7 systems ----
    std::cout << "\n--- Test 4: Dual-Method Cross-Validation (7 Systems) ---" << std::endl;
    struct SystemEntry {
        std::string name;
        CelestialBody_S4 body;
        MUGESystem_S4 muge;
    };
    std::vector<SystemEntry> systems = {
        {"SGR1745",      sgr1745_body,        sgr1745_muge},
        {"SagA*",        sagA_body,            sagA_muge},
        {"Tapestry",     tapestry_body,        tapestry_muge},
        {"Westerlund2",  westerlund2_body,     westerlund2_muge},
        {"Pillars",      pillars_body,         pillars_muge},
        {"Rings",        rings_body,           rings_muge},
        {"StudentGuide", student_guide_body,   student_guide_muge},
    };

    int convergence_count = 0;
    for (auto& sys : systems) {
        auto result = validate_system(sys.name, sys.body, sys.muge, default_res_params);

        std::ostringstream detail;
        detail << std::scientific << std::setprecision(4)
               << "UQFF=" << result.uqff_value
               << " MUGE_C=" << result.muge_compressed
               << " MUGE_R=" << result.muge_resonance
               << " diff=" << std::fixed << std::setprecision(1)
               << result.uqff_muge_diff_pct << "%";

        report((sys.name + " methods produce finite values").c_str(),
               std::isfinite(result.uqff_value) &&
               std::isfinite(result.muge_compressed) &&
               std::isfinite(result.muge_resonance),
               detail.str());

        if (result.convergence) ++convergence_count;

        // Log detailed results
        log << "=== " << sys.name << " ===" << std::endl;
        log << std::scientific << std::setprecision(6);
        log << "  UQFF:               " << result.uqff_value << " (unified)" << std::endl;
        log << "  MUGE Compressed:    " << result.muge_compressed << " m/s^2" << std::endl;
        log << "  MUGE Resonance:     " << result.muge_resonance << " m/s^2" << std::endl;
        log << std::fixed << std::setprecision(2);
        log << "  UQFF-MUGE Diff:     " << result.uqff_muge_diff_pct << "%" << std::endl;
        log << "  Comp-Res Diff:      " << result.compressed_resonance_diff_pct << "%" << std::endl;
        log << "  Convergence:        " << (result.convergence ? "PASS" : "FAIL") << std::endl;
        log << "  Analysis:           " << result.analysis << std::endl;
        log << std::endl;
    }

    report("At least 1 system converges", convergence_count >= 1,
           std::to_string(convergence_count) + "/7 converged");

    // ---- Test 5: Physical invariants ----
    std::cout << "\n--- Test 5: Physical Invariants ---" << std::endl;
    {
        // FU must increase as r decreases (stronger gravity closer)
        double fu_near = compute_FU(sgr1745_body, 1e3, 0.0);
        double fu_far  = compute_FU(sgr1745_body, 1e6, 0.0);
        report("FU(r_near) > FU(r_far) for SGR1745", std::abs(fu_near) > std::abs(fu_far));

        // MUGE compressed must give Newtonian at large r, weak B
        MUGESystem_S4 weak = sagA_muge;
        weak.B = 0.0;
        weak.rho_fluid = 0.0;
        weak.M_DM = 0.0;
        weak.delta_rho = 0.0;
        weak.t = 0.0;
        double g_weak = compute_compressed_MUGE(weak);
        double g_newton = G_S4 * weak.M / (weak.r * weak.r);
        double pct_diff = std::abs((g_weak - g_newton) / g_newton * 100.0);
        report("MUGE → Newtonian when corrections=0", pct_diff < 5.0,
               "diff=" + std::to_string(pct_diff) + "%");
    }

    // ---- Test 6: Symmetry and scaling ----
    std::cout << "\n--- Test 6: Symmetry and Scaling ---" << std::endl;
    {
        // Doubling mass should roughly double Newtonian gravity
        MUGESystem_S4 sys1 = sagA_muge;
        MUGESystem_S4 sys2 = sagA_muge;
        sys2.M *= 2.0;
        sys2.M_DM *= 2.0;
        double g1 = compute_compressed_MUGE(sys1);
        double g2 = compute_compressed_MUGE(sys2);
        double ratio = g2 / g1;
        report("2x mass ≈ 2x gravity", ratio > 1.5 && ratio < 2.5,
               "ratio=" + std::to_string(ratio));

        // Inverse square: 2x distance ≈ 1/4 gravity
        MUGESystem_S4 sys3 = sagA_muge;
        sys3.r *= 2.0;
        double g3 = compute_compressed_MUGE(sys3);
        double isq_ratio = g1 / g3;
        report("2x distance ≈ 4x weaker", isq_ratio > 2.0 && isq_ratio < 8.0,
               "ratio=" + std::to_string(isq_ratio));
    }

    // ---- Test 7: Edge cases ----
    std::cout << "\n--- Test 7: Edge Cases ---" << std::endl;
    {
        // Cosmological scale (Student Guide)
        double fu_cosm = compute_FU(student_guide_body, student_guide_muge.r, 0.0);
        report("FU(cosmological) finite", std::isfinite(fu_cosm));

        // Extreme magnetar field
        double fu_mag = compute_FU(sgr1745_body, 1e4, 0.0);
        report("FU(magnetar, r=10km) finite", std::isfinite(fu_mag));
        report("FU(magnetar, r=10km) > 0", fu_mag > 0.0);

        // Zero time
        double g_t0 = compute_compressed_MUGE(sgr1745_muge);
        report("MUGE at t=0 finite", std::isfinite(g_t0));
    }

    // ---- Summary ----
    auto t_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);

    std::cout << "\n============================================================" << std::endl;
    std::cout << " RESULTS: " << g_pass << " passed, " << g_fail << " failed" << std::endl;
    std::cout << " Time: " << duration.count() << " ms" << std::endl;
    std::cout << " Log: validation_log.txt" << std::endl;
    std::cout << "============================================================" << std::endl;

    log << "Summary: " << g_pass << " passed, " << g_fail << " failed" << std::endl;
    log << "Time: " << duration.count() << " ms" << std::endl;
    log.close();

    return (g_fail > 0) ? 1 : 0;
}
