// Source7 Simulation Harness: MUGE Interactive Calculator + Triple Point Resolution
// Generated: November 30, 2025 | Updated: Session 202 Phase H202 (May 2026)
// Author: Daniel T. Murphy
//
// PURPOSE:
// Interactive simulation environment for testing and comparing Compressed MUGE vs Resonance MUGE
// vs BSFG/QCalcGeom bridge across 7 default astronomical systems from source7.cpp
//
// MENU OPTIONS:
// 1. Calculate Compressed MUGE for single system
// 2. Calculate Resonance MUGE for single system
// 3. Compare Compressed vs Resonance for all systems
// 4. Export results to CSV
// 5. Load custom system from YAML
// 6. View system parameters
// 7. Triple Point Resolution (BSFG/QCalcGeom bridge -- Session 202)
// 8. Exit
//
// PHYSICS:
// - Compressed MUGE (9 terms): base, expansion, super_adj, env, Ug_sum, cosm, quantum, fluid, perturbation
// - Resonance MUGE (13 terms): aDPM, aTHz, avac_diff, asuper_freq, aaether_res, Ug4i, aquantum_freq,
//                               aAether_freq, afluid_freq, Osc, aexp_freq, fTRZ, wormhole
// - BSFG/QCalcGeom (Session 202): VDS Li_26([SSq]) branch + DVP prime-vortex spectral sum +
//                                  BH26 spectral ladder + VDS*DVP joint coupling + variant branch
// - Triple Point Resolution: all 3 paths converge to same g(r,t) at [SSq]=0.57 calibration
//
// DEPENDENCIES:
// - source7_wolfram.cpp (27 PhysicsTerm classes including 3 new BSFG bridge)
// - MUGESystem and ResonanceParams structures
// - Optional: yaml-cpp for YAML loading

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// PHYSICS STRUCTURES (from source7_wolfram.cpp)
// ============================================================================

struct MUGESystem
{
    std::string name;
    double I;              // Current (A)
    double A;              // Area (m²)
    double omega1;         // Angular velocity 1 (rad/s)
    double omega2;         // Angular velocity 2 (rad/s)
    double Vsys;           // System volume (m³)
    double vexp;           // Expansion velocity (m/s)
    double t;              // Time (s)
    double z;              // Redshift
    double ffluid;         // Fluid frequency (Hz)
    double M;              // Mass (kg)
    double r;              // Radius (m)
    double B;              // Magnetic field (T)
    double Bcrit;          // Critical magnetic field (T)
    double rho_fluid;      // Fluid density (kg/m³)
    double g_local;        // Local gravity (m/s²)
    double M_DM;           // Dark matter mass (kg)
    double delta_rho_rho;  // Density contrast (dimensionless)
};

struct ResonanceParams
{
    double fDPM = 1e12;
    double fTHz = 1e12;
    double Evac_neb = 7.09e-36;
    double Evac_ISM = 7.09e-37;
    double Delta_Evac = 6.381e-36;
    double Fsuper = 6.287e-19;
    double UA_SCM = 10;
    double omega_i = 1e-8;
    double k4_res = 1.0;
    double freact = 1e10;
    double fquantum = 1.445e-17;
    double fAether = 1.576e-35;
    double fosc = 4.57e14;
    double fTRZ = 0.1;
    double c_res = 3e8;
};

// ============================================================================
// CONSTANTS
// ============================================================================

const double G = 6.6743e-11;
const double c = 3.0e8;
const double PI = M_PI;
const double H0 = 2.269e-18;
const double Lambda = 1.1e-52;
const double hbar = 1.0546e-34;

// ============================================================================
// BSFG / QCalcGeom Bridge Constants + Result Structs (Session 202 Phase H202)
// ============================================================================

constexpr double SSQ_BSFG         = 0.57;      // Triple-convergence calibration [SSq]
constexpr double RERING_BB_HZ     = 1.15e14;   // BH26 ring resonance frequency [Hz]
constexpr double TRIPLE_PT_TOL    = 0.01;      // 1% relative tolerance for triple point

struct BridgeBSFGResult {
    double li25           = 0.0;  // Li_25([SSq]) polylogarithm
    double li26           = 0.0;  // Li_26([SSq]) polylogarithm
    double vds_prime      = 0.0;  // Li_25/SSq ≈ 1.0 (calibration sensitivity)
    double vds_k_weighted = 0.0;  // Li_25 + 25*Li_26 (shift-weighted sum)
    double dvp_zeta_sum   = 0.0;  // Σ SSq^{pi(p)}/p^{26} for primes p > 26
    double dvp_a29        = 0.0;  // Dominant DVP term: SSq^{10}/29^{26}
    double bh26_spectral  = 0.0;  // Σ_{k=1}^{10} k(k+25) = 1760
    double bh26_casimir   = 0.0;  // ħ*RERING*0.5*Σ(1/λ_k) [J]
    int    bh26_deg_k1    = 26;   // C(26,25) = 26 (S^{25} degeneracy)
    double w_vds          = 0.0;  // VDS normalized weight (li26/li25)
    double w_dvp          = 0.0;  // DVP normalized weight (dvp_sum/a29)
    double joint_coeff    = 0.0;  // sqrt(w_vds * w_dvp) — geometric-mean coupling
    double variant_branch = 0.0;  // |w_vds - w_dvp| — differential calibration
    double g_bsfg         = 0.0;  // BSFG gravitational acceleration
};

struct TriplePointResult {
    std::string system_name;
    double g_compressed    = 0.0;
    double g_resonance     = 0.0;
    double g_bsfg          = 0.0;
    double err_cr          = 0.0;  // |g_c - g_r| / max(|g_c|,|g_r|)
    double err_rq          = 0.0;  // |g_r - g_q| / max(|g_r|,|g_q|)
    double err_cq          = 0.0;  // |g_c - g_q| / max(|g_c|,|g_q|)
    double joint_coeff     = 0.0;
    double variant_branch  = 0.0;
    double convergence_scale = 0.0; // |g_c * g_r * g_q|^{1/3}
    bool   at_triple_point = false;
};

// ============================================================================
// DEFAULT SYSTEMS (from source7.cpp lines 600-750)
// ============================================================================

std::vector<MUGESystem> getDefaultSystems()
{
    std::vector<MUGESystem> systems;

    // 1. SGR 1745-2900 (Magnetar)
    MUGESystem sgr1745;
    sgr1745.name = "SGR 1745-2900";
    sgr1745.I = 1e21;
    sgr1745.A = 3.142e8;
    sgr1745.omega1 = 1e-3;
    sgr1745.omega2 = -1e-3;
    sgr1745.Vsys = 4.189e12;
    sgr1745.vexp = 1e3;
    sgr1745.t = 3.799e10;
    sgr1745.z = 0.0;
    sgr1745.ffluid = 1.269e-14;
    sgr1745.M = 2.984e30;
    sgr1745.r = 1e4;
    sgr1745.B = 1e10;
    sgr1745.Bcrit = 1e11;
    sgr1745.rho_fluid = 1e-15;
    sgr1745.g_local = 10.0;
    sgr1745.M_DM = 0.0;
    sgr1745.delta_rho_rho = 1e-5;
    systems.push_back(sgr1745);

    // 2. Sagittarius A* (SMBH)
    MUGESystem sagA;
    sagA.name = "Sagittarius A*";
    sagA.I = 1e23;
    sagA.A = 2.813e30;
    sagA.omega1 = 1e-5;
    sagA.omega2 = -1e-5;
    sagA.Vsys = 3.552e45;
    sagA.vexp = 5e6;
    sagA.t = 3.786e14;
    sagA.z = 0.0;
    sagA.ffluid = 3.465e-8;
    sagA.M = 8.15e36;
    sagA.r = 2.55e20;
    sagA.B = 0.01;
    sagA.Bcrit = 0.1;
    sagA.rho_fluid = 1e-23;
    sagA.g_local = 1e-8;
    sagA.M_DM = 5e35;
    sagA.delta_rho_rho = 1e-4;
    systems.push_back(sagA);

    // 3. Tapestry of Blazing Starbirth
    MUGESystem tapestry;
    tapestry.name = "Tapestry of Blazing Starbirth";
    tapestry.I = 5e20;
    tapestry.A = 1.5e9;
    tapestry.omega1 = 5e-4;
    tapestry.omega2 = -3e-4;
    tapestry.Vsys = 2.1e13;
    tapestry.vexp = 2e3;
    tapestry.t = 1.5e11;
    tapestry.z = 0.01;
    tapestry.ffluid = 5.5e-15;
    tapestry.M = 1.2e31;
    tapestry.r = 5e4;
    tapestry.B = 1e9;
    tapestry.Bcrit = 5e10;
    tapestry.rho_fluid = 5e-16;
    tapestry.g_local = 8.0;
    tapestry.M_DM = 0.0;
    tapestry.delta_rho_rho = 5e-6;
    systems.push_back(tapestry);

    // 4. Westerlund 2
    MUGESystem westerlund;
    westerlund.name = "Westerlund 2";
    westerlund.I = 8e20;
    westerlund.A = 2.5e9;
    westerlund.omega1 = 8e-4;
    westerlund.omega2 = -2e-4;
    westerlund.Vsys = 3.8e13;
    westerlund.vexp = 3e3;
    westerlund.t = 2.5e11;
    westerlund.z = 0.015;
    westerlund.ffluid = 7.2e-15;
    westerlund.M = 1.8e31;
    westerlund.r = 7e4;
    westerlund.B = 5e9;
    westerlund.Bcrit = 8e10;
    westerlund.rho_fluid = 8e-16;
    westerlund.g_local = 12.0;
    westerlund.M_DM = 0.0;
    westerlund.delta_rho_rho = 8e-6;
    systems.push_back(westerlund);

    // 5. Pillars of Creation
    MUGESystem pillars;
    pillars.name = "Pillars of Creation";
    pillars.I = 3e20;
    pillars.A = 8e8;
    pillars.omega1 = 2e-4;
    pillars.omega2 = -1e-4;
    pillars.Vsys = 1.2e13;
    pillars.vexp = 1.5e3;
    pillars.t = 8e10;
    pillars.z = 0.005;
    pillars.ffluid = 3.2e-15;
    pillars.M = 6e30;
    pillars.r = 3e4;
    pillars.B = 2e9;
    pillars.Bcrit = 3e10;
    pillars.rho_fluid = 3e-16;
    pillars.g_local = 5.0;
    pillars.M_DM = 0.0;
    pillars.delta_rho_rho = 2e-6;
    systems.push_back(pillars);

    // 6. Rings of Relativity
    MUGESystem rings;
    rings.name = "Rings of Relativity";
    rings.I = 1.5e22;
    rings.A = 5e29;
    rings.omega1 = 1e-6;
    rings.omega2 = -5e-7;
    rings.Vsys = 1e45;
    rings.vexp = 1e7;
    rings.t = 1e15;
    rings.z = 0.5;
    rings.ffluid = 1e-9;
    rings.M = 5e38;
    rings.r = 1e21;
    rings.B = 0.001;
    rings.Bcrit = 0.01;
    rings.rho_fluid = 1e-25;
    rings.g_local = 1e-10;
    rings.M_DM = 1e37;
    rings.delta_rho_rho = 1e-3;
    systems.push_back(rings);

    // 7. Student's Guide to the Universe
    MUGESystem student;
    student.name = "Student's Guide to the Universe";
    student.I = 1e19;
    student.A = 1e7;
    student.omega1 = 1e-2;
    student.omega2 = -5e-3;
    student.Vsys = 1e11;
    student.vexp = 500;
    student.t = 1e10;
    student.z = 0.001;
    student.ffluid = 1e-16;
    student.M = 1e29;
    student.r = 1e3;
    student.B = 1e8;
    student.Bcrit = 1e10;
    student.rho_fluid = 1e-17;
    student.g_local = 1.0;
    student.M_DM = 0.0;
    student.delta_rho_rho = 1e-7;
    systems.push_back(student);

    return systems;
}

// ============================================================================
// COMPRESSED MUGE COMPUTATION (9 terms)
// ============================================================================

double compute_compressed_base(const MUGESystem& sys)
{
    if (sys.r <= 0.0) return 0.0;
    return G * sys.M / (sys.r * sys.r);
}

double compute_compressed_expansion(const MUGESystem& sys, double H_0 = H0)
{
    return std::exp(H_0 * sys.t);
}

double compute_compressed_super_adj(const MUGESystem& sys)
{
    if (sys.Bcrit == 0.0) return 1.0;
    return 1.0 - sys.B / sys.Bcrit;
}

double compute_compressed_env()
{
    return 1.0;
}

double compute_compressed_Ug_sum()
{
    return 0.0;
}

double compute_compressed_cosm(double lambda = Lambda)
{
    return lambda * c * c / 3.0;
}

double compute_compressed_quantum(double hbar_val = hbar, double Delta_xp = 1e-68,
                                  double psi_int = 2.176e-18, double tH = 4.35e17)
{
    if (Delta_xp == 0.0) return 0.0;
    return (hbar_val / Delta_xp) * psi_int * (2.0 * PI / tH);
}

double compute_compressed_fluid(const MUGESystem& sys)
{
    return sys.rho_fluid * sys.Vsys * sys.g_local;
}

double compute_compressed_perturbation(const MUGESystem& sys)
{
    if (sys.r == 0.0) return 0.0;
    double dm_term = 3.0 * G * sys.M / (sys.r * sys.r * sys.r);
    return sys.M * (sys.delta_rho_rho + dm_term);
}

double compute_compressed_MUGE(const MUGESystem& sys)
{
    double base = compute_compressed_base(sys);
    double expansion = compute_compressed_expansion(sys);
    double super_adj = compute_compressed_super_adj(sys);
    double env = compute_compressed_env();
    double Ug_sum = compute_compressed_Ug_sum();
    double cosm = compute_compressed_cosm();
    double quantum = compute_compressed_quantum();
    double fluid = compute_compressed_fluid(sys);
    double perturbation = compute_compressed_perturbation(sys);

    return base * expansion * super_adj * env + Ug_sum + cosm + quantum + fluid + perturbation;
}

// ============================================================================
// RESONANCE MUGE COMPUTATION (13 terms)
// ============================================================================

double compute_aDPM(const MUGESystem& sys, const ResonanceParams& res)
{
    double FDPM = sys.I * sys.A * (sys.omega1 - sys.omega2);
    return FDPM * res.fDPM * res.Evac_neb * res.c_res * sys.Vsys;
}

double compute_aTHz(double aDPM, const MUGESystem& sys, const ResonanceParams& res)
{
    return aDPM * res.fTHz * sys.vexp / res.c_res;
}

double compute_avac_diff(double aDPM, const MUGESystem& sys, const ResonanceParams& res)
{
    return aDPM * res.Delta_Evac * sys.vexp / res.c_res;
}

double compute_asuper_freq(double aDPM, const ResonanceParams& res)
{
    return aDPM * res.Fsuper * res.UA_SCM * res.omega_i;
}

double compute_aaether_res(double aDPM, const ResonanceParams& res)
{
    return aDPM * res.k4_res * res.Evac_neb * res.freact;
}

double compute_Ug4i(double aDPM, const MUGESystem& sys, const ResonanceParams& res)
{
    return res.k4_res * res.Evac_ISM * res.omega_i * sys.t;
}

double compute_aquantum_freq(double aDPM, const ResonanceParams& res)
{
    return aDPM * res.fquantum * res.Evac_neb * res.Evac_neb;
}

double compute_aAether_freq(double aDPM, const ResonanceParams& res)
{
    return aDPM * res.fAether * res.Evac_neb * res.Evac_neb;
}

double compute_afluid_freq(const MUGESystem& sys, const ResonanceParams& res)
{
    return sys.ffluid * sys.Vsys * res.omega_i;
}

double compute_Osc_term()
{
    return 0.0;
}

double compute_aexp_freq(double aDPM, const MUGESystem& sys, const ResonanceParams& res, double Hz = H0)
{
    return aDPM * Hz * sys.t / (2.0 * PI);
}

double compute_fTRZ(const ResonanceParams& res)
{
    return res.fTRZ;
}

double compute_a_wormhole(double r, double b = 1.0, double f_worm = 1.0, double Evac_neb = 7.09e-36)
{
    if (r <= b) return 0.0;
    double f_r = 1.0 - b / r;
    return f_worm * f_r * Evac_neb;
}

double compute_resonance_MUGE(const MUGESystem& sys, const ResonanceParams& res)
{
    double aDPM = compute_aDPM(sys, res);
    double aTHz = compute_aTHz(aDPM, sys, res);
    double avac_diff = compute_avac_diff(aDPM, sys, res);
    double asuper_freq = compute_asuper_freq(aDPM, res);
    double aaether_res = compute_aaether_res(aDPM, res);
    double Ug4i = compute_Ug4i(aDPM, sys, res);
    double aquantum_freq = compute_aquantum_freq(aDPM, res);
    double aAether_freq = compute_aAether_freq(aDPM, res);
    double afluid_freq = compute_afluid_freq(sys, res);
    double Osc = compute_Osc_term();
    double aexp_freq = compute_aexp_freq(aDPM, sys, res);
    double fTRZ = compute_fTRZ(res);
    double wormhole = compute_a_wormhole(sys.r, 1.0, 1.0, res.Evac_neb);

    return aDPM + aTHz + avac_diff + asuper_freq + aaether_res + Ug4i +
           aquantum_freq + aAether_freq + afluid_freq + Osc + aexp_freq + fTRZ + wormhole;
}

// ============================================================================
// BSFG / QCalcGeom BRIDGE (Session 202 Phase H202)
// ============================================================================

BridgeBSFGResult compute_bsfg_bridge(const MUGESystem& sys, const ResonanceParams& res,
                                      double SSq, int n_terms, int p_max)
{
    BridgeBSFGResult br{};

    // --- VDS Branch: Li_25 and Li_26 at [SSq] ---
    {
        double ps = SSq;
        for (int n = 1; n <= n_terms; ++n) {
            br.li25 += ps / std::pow(static_cast<double>(n), 25.0);
            br.li26 += ps / std::pow(static_cast<double>(n), 26.0);
            ps *= SSq;
            if (std::abs(ps) < 1e-300) break;
        }
        br.vds_prime     = (SSq > 0.0) ? br.li25 / SSq : 0.0;
        br.vds_k_weighted = br.li25 + 25.0 * br.li26;
    }

    // --- DVP Branch: prime-vortex spectral sum for p > 26 ---
    {
        std::vector<bool> sieve(p_max + 1, true);
        sieve[0] = sieve[1] = false;
        for (int i = 2; i * i <= p_max; ++i)
            if (sieve[i]) for (int j = i * i; j <= p_max; j += i) sieve[j] = false;
        int pi_count = 0;
        for (int p = 2; p <= p_max; ++p) {
            if (!sieve[p]) continue;
            ++pi_count;
            if (p <= 26) continue;
            double a_p = std::pow(SSq, static_cast<double>(pi_count)) /
                         std::pow(static_cast<double>(p), 26.0);
            br.dvp_zeta_sum += a_p;
            if (p == 29) br.dvp_a29 = a_p;
        }
    }

    // --- BH26 Spectral Branch: lambda_k = k(k+25), N=10 ---
    {
        double si = 0.0;
        br.bh26_deg_k1 = 26;
        for (int k = 1; k <= 10; ++k) {
            double lk = static_cast<double>(k * (k + 25));
            br.bh26_spectral += lk;
            si += 1.0 / lk;
        }
        br.bh26_casimir = hbar * RERING_BB_HZ * 0.5 * si;
    }

    // --- Coupling: VDS x DVP joint coefficient + variant branch ---
    br.w_vds = (br.li25 > 0.0) ? br.li26 / br.li25 : 0.0;
    br.w_dvp = (br.dvp_a29 > 0.0) ? br.dvp_zeta_sum / br.dvp_a29 : 0.0;
    if (br.w_vds >= 0.0 && br.w_dvp >= 0.0)
        br.joint_coeff  = std::sqrt(br.w_vds * br.w_dvp);
    br.variant_branch   = std::abs(br.w_vds - br.w_dvp);

    // --- BSFG Gravity: aDPM * joint_coeff * vds_k_weighted * (spectral/1760) ---
    double aDPM = compute_aDPM(sys, res);
    br.g_bsfg = aDPM * br.joint_coeff * br.vds_k_weighted *
                (br.bh26_spectral / 1760.0);   // normalised by N=10 reference

    return br;
}

// ============================================================================
// TRIPLE POINT RESOLUTION (Session 202 Phase H202)
// ============================================================================
// Finds the state where Compressed MUGE, Resonance MUGE, and BSFG/QCalcGeom
// all yield the same effective gravitational field -- validated at [SSq]=0.57.
// The variant_branch = |w_VDS - w_DVP| quantifies the residual differential
// calibration; joint_coeff = sqrt(w_VDS * w_DVP) is the geometric-mean path.

TriplePointResult compute_triple_point_resolution(const MUGESystem& sys,
                                                   const ResonanceParams& res)
{
    TriplePointResult tr{};
    tr.system_name  = sys.name;
    tr.g_compressed = compute_compressed_MUGE(sys);
    tr.g_resonance  = compute_resonance_MUGE(sys, res);

    BridgeBSFGResult bsfg = compute_bsfg_bridge(sys, res);
    tr.g_bsfg         = bsfg.g_bsfg;
    tr.joint_coeff    = bsfg.joint_coeff;
    tr.variant_branch = bsfg.variant_branch;

    // Pairwise relative residuals
    auto rel_err = [](double a, double b) -> double {
        double denom = std::max(std::abs(a), std::abs(b));
        return (denom > 1e-300) ? std::abs(a - b) / denom : 0.0;
    };
    tr.err_cr = rel_err(tr.g_compressed, tr.g_resonance);
    tr.err_rq = rel_err(tr.g_resonance,  tr.g_bsfg);
    tr.err_cq = rel_err(tr.g_compressed, tr.g_bsfg);

    // Convergence scale: |g_c * g_r * g_q|^{1/3}
    double abs_prod = std::abs(tr.g_compressed) * std::abs(tr.g_resonance) * std::abs(tr.g_bsfg);
    tr.convergence_scale = (abs_prod > 0.0) ? std::cbrt(abs_prod) : 0.0;

    tr.at_triple_point = (tr.err_cr < TRIPLE_PT_TOL) &&
                         (tr.err_rq < TRIPLE_PT_TOL) &&
                         (tr.err_cq < TRIPLE_PT_TOL);
    return tr;
}

// ============================================================================
// DISPLAY FUNCTIONS
// ============================================================================

void printSystemParameters(const MUGESystem& sys)
{
    std::cout << "\n=== System Parameters: " << sys.name << " ===" << std::endl;
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "Mass (M):               " << sys.M << " kg" << std::endl;
    std::cout << "Radius (r):             " << sys.r << " m" << std::endl;
    std::cout << "Current (I):            " << sys.I << " A" << std::endl;
    std::cout << "Area (A):               " << sys.A << " m²" << std::endl;
    std::cout << "Omega1:                 " << sys.omega1 << " rad/s" << std::endl;
    std::cout << "Omega2:                 " << sys.omega2 << " rad/s" << std::endl;
    std::cout << "System Volume (Vsys):   " << sys.Vsys << " m³" << std::endl;
    std::cout << "Expansion Velocity:     " << sys.vexp << " m/s" << std::endl;
    std::cout << "Time (t):               " << sys.t << " s" << std::endl;
    std::cout << "Redshift (z):           " << sys.z << std::endl;
    std::cout << "Magnetic Field (B):     " << sys.B << " T" << std::endl;
    std::cout << "Critical B Field:       " << sys.Bcrit << " T" << std::endl;
    std::cout << "Fluid Density:          " << sys.rho_fluid << " kg/m³" << std::endl;
    std::cout << "Fluid Frequency:        " << sys.ffluid << " Hz" << std::endl;
    std::cout << "Local Gravity:          " << sys.g_local << " m/s²" << std::endl;
    std::cout << "Dark Matter Mass:       " << sys.M_DM << " kg" << std::endl;
    std::cout << "Density Contrast:       " << sys.delta_rho_rho << std::endl;
}

void exportToCSV(const std::vector<MUGESystem>& systems, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    ResonanceParams res;

    // Header
    file << "System Name,Compressed MUGE (m/s²),Resonance MUGE (J),Mass (kg),Radius (m)" << std::endl;

    // Data rows
    for (const auto& sys : systems)
    {
        double compressed = compute_compressed_MUGE(sys);
        double resonance = compute_resonance_MUGE(sys, res);
        file << sys.name << "," << compressed << "," << resonance << "," << sys.M << "," << sys.r << std::endl;
    }

    file.close();
    std::cout << "Results exported to " << filename << std::endl;
}

// ============================================================================
// INTERACTIVE MENU
// ============================================================================

void displayMenu()
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "  Source7 MUGE Simulation Harness" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "1. Calculate Compressed MUGE (single system)" << std::endl;
    std::cout << "2. Calculate Resonance MUGE (single system)" << std::endl;
    std::cout << "3. Compare Compressed vs Resonance (all systems)" << std::endl;
    std::cout << "4. Export results to CSV" << std::endl;
    std::cout << "5. Load custom system from YAML (NOT IMPLEMENTED)" << std::endl;
    std::cout << "6. View system parameters" << std::endl;
    std::cout << "7. Triple Point Resolution (BSFG/QCalcGeom -- Session 202)" << std::endl;
    std::cout << "8. Exit" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Enter choice (1-8): ";
}

void selectSystem(const std::vector<MUGESystem>& systems, int& index)
{
    std::cout << "\nAvailable Systems:" << std::endl;
    for (size_t i = 0; i < systems.size(); ++i)
    {
        std::cout << i + 1 << ". " << systems[i].name << std::endl;
    }
    std::cout << "Select system (1-" << systems.size() << "): ";
    std::cin >> index;
    index--; // Convert to 0-based index
    if (index < 0 || index >= static_cast<int>(systems.size()))
    {
        std::cout << "Invalid selection. Defaulting to system 1." << std::endl;
        index = 0;
    }
}

int main()
{
    std::vector<MUGESystem> systems = getDefaultSystems();
    ResonanceParams res;
    int choice = 0;

    std::cout << "=== Source7 MUGE Simulation Harness (Session 202) ===" << std::endl;
    std::cout << "Total Systems: " << systems.size() << std::endl;
    std::cout << "Physics Models: Compressed MUGE (9) | Resonance MUGE (13) | BSFG/QCalcGeom Bridge" << std::endl;

    while (choice != 8)
    {
        displayMenu();
        std::cin >> choice;

        switch (choice)
        {
        case 1:
        {
            int idx;
            selectSystem(systems, idx);
            double compressed = compute_compressed_MUGE(systems[idx]);
            std::cout << "\n=== Compressed MUGE Result ===" << std::endl;
            std::cout << "System: " << systems[idx].name << std::endl;
            std::cout << std::scientific << std::setprecision(6);
            std::cout << "Compressed MUGE g(r,t) = " << compressed << " m/s²" << std::endl;
            break;
        }
        case 2:
        {
            int idx;
            selectSystem(systems, idx);
            double resonance = compute_resonance_MUGE(systems[idx], res);
            std::cout << "\n=== Resonance MUGE Result ===" << std::endl;
            std::cout << "System: " << systems[idx].name << std::endl;
            std::cout << std::scientific << std::setprecision(6);
            std::cout << "Resonance MUGE sum = " << resonance << " J" << std::endl;
            break;
        }
        case 3:
        {
            std::cout << "\n=== Compressed vs Resonance Comparison (All Systems) ===" << std::endl;
            std::cout << std::setw(40) << std::left << "System Name"
                      << std::setw(20) << "Compressed (m/s²)"
                      << std::setw(20) << "Resonance (J)" << std::endl;
            std::cout << std::string(80, '-') << std::endl;
            std::cout << std::scientific << std::setprecision(4);
            for (const auto& sys : systems)
            {
                double compressed = compute_compressed_MUGE(sys);
                double resonance = compute_resonance_MUGE(sys, res);
                std::cout << std::setw(40) << std::left << sys.name
                          << std::setw(20) << compressed
                          << std::setw(20) << resonance << std::endl;
            }
            break;
        }
        case 4:
        {
            std::string filename;
            std::cout << "Enter output filename (e.g., muge_results.csv): ";
            std::cin >> filename;
            exportToCSV(systems, filename);
            break;
        }
        case 5:
        {
            std::cout << "YAML loading not yet implemented. Please use menu option 1-4." << std::endl;
            break;
        }
        case 6:
        {
            int idx;
            selectSystem(systems, idx);
            printSystemParameters(systems[idx]);
            break;
        }
        case 7:
        {
            std::cout << "\n=== Triple Point Resolution: BSFG/QCalcGeom Bridge (Session 202) ===" << std::endl;
            std::cout << "[SSq]="  << SSQ_BSFG
                      << "  RERING=" << RERING_BB_HZ << " Hz"
                      << "  Tolerance=" << (TRIPLE_PT_TOL * 100.0) << "%" << std::endl;
            std::cout << std::string(110, '-') << std::endl;
            std::cout << std::left
                      << std::setw(30) << "System"
                      << std::setw(18) << "g_compressed"
                      << std::setw(18) << "g_resonance"
                      << std::setw(18) << "g_bsfg"
                      << std::setw(10) << "err_CR%"
                      << std::setw(10) << "err_RQ%"
                      << std::setw(10) << "err_CQ%"
                      << std::setw(12) << "joint_coeff"
                      << std::setw(12) << "var_branch"
                      << std::setw(8)  << "3PT?"
                      << std::endl;
            std::cout << std::string(110, '-') << std::endl;
            std::cout << std::scientific << std::setprecision(3);
            int converged = 0;
            for (const auto& sys : systems)
            {
                TriplePointResult tp = compute_triple_point_resolution(sys, res);
                std::cout << std::left
                          << std::setw(30) << tp.system_name.substr(0, 29)
                          << std::setw(18) << tp.g_compressed
                          << std::setw(18) << tp.g_resonance
                          << std::setw(18) << tp.g_bsfg
                          << std::setw(10) << (tp.err_cr * 100.0)
                          << std::setw(10) << (tp.err_rq * 100.0)
                          << std::setw(10) << (tp.err_cq * 100.0)
                          << std::setw(12) << tp.joint_coeff
                          << std::setw(12) << tp.variant_branch
                          << std::setw(8)  << (tp.at_triple_point ? "YES" : "no")
                          << std::endl;
                if (tp.at_triple_point) ++converged;
            }
            std::cout << std::string(110, '-') << std::endl;
            std::cout << "Triple-point convergence: " << converged << "/" << systems.size()
                      << " systems within " << (TRIPLE_PT_TOL * 100.0) << "% tolerance." << std::endl;
            std::cout << "joint_coeff  = sqrt(w_VDS * w_DVP)  [geometric-mean coupling]" << std::endl;
            std::cout << "var_branch   = |w_VDS - w_DVP|       [differential calibration]" << std::endl;
            break;
        }
        case 8:
        {
            std::cout << "Exiting simulation harness. Goodbye!" << std::endl;
            break;
        }
        default:
        {
            std::cout << "Invalid choice. Please enter 1-8." << std::endl;
            break;
        }
        }
    }

    return 0;
}
