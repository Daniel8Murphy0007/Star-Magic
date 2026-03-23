// UQFFCassiniBuoyancyModule.cpp
// Implementation of UQFF Cassini Mission Buoyancy calculations
// Based on "Buoyancy_txt_12Oct2025.docx"
// Implements complex-number UQFF equations for Saturn ring gaps:
// Encke Gap (133,590 km), Cassini Division (117,500–122,200 km),
// Maxwell Gap (87,500 km from center)
// Watermark: Copyright - Daniel T. Murphy, analyzed by Grok 3, dated November 17, 2025

#include "../UQFFCassiniBuoyancyModule.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>

using namespace std::complex_literals;

// ============================================================
// UQFFCassiniCore: Constructor
// ============================================================
UQFFCassiniCore::UQFFCassiniCore(double k1, double ki, double km, double ke)
    : k1_(k1), ki_(ki), km_(km), ke_(ke) {}

// ============================================================
// calculate_U_Mi: Universal Magnetism with Heaviside reverse polarity
// U_Mi(t) = km * B * exp(-gamma * t) * H_rev(cos(phase * t))
// Complex component: imaginary part captures polarity reversal
// ============================================================
std::complex<double> UQFFCassiniCore::calculate_U_Mi(double t, double r, int n) {
    double r_eff = (r > 0) ? r : 1.0;
    double B_eff = CASS_B_GRADIENT / (r_eff * r_eff);
    double decay = std::exp(-CASS_GAMMA_DECAY * t);
    std::complex<double> cos_phase(std::cos(CASS_PHASE * t), 0.0);
    std::complex<double> H_rev = heaviside_reverse(cos_phase);
    // Imaginary component: n-th harmonic Landau level contribution
    double landau = static_cast<double>(n) + 0.5;
    return std::complex<double>(km_ * B_eff * decay * H_rev.real(),
                                km_ * B_eff * decay * landau * CASS_RHO_VAC_UA);
}

// ============================================================
// calculate_U_Ii: Universal Inertia (gyroscopic mimic)
// U_Ii = gyro_principle(U_Mi, omega) = U_Mi * exp(i * omega * PI)
// ============================================================
std::complex<double> UQFFCassiniCore::calculate_U_Ii(
    const std::complex<double>& U_Mi_val, double gyro_factor) {
    double omega = CASS_NU_THz * gyro_factor; // angular frequency
    return gyro_principle(U_Mi_val, omega);
}

// ============================================================
// calculate_U_Bi: Universal Buoyancy (calibration difference)
// U_Bi = ki * (rho_vac_UA - rho_vac_SCm) * r^3 / r^2
//       = ki * delta_k / r
// ============================================================
std::complex<double> UQFFCassiniCore::calculate_U_Bi(double delta_k) {
    // delta_k = calibration offset between [UA] and [SCm] vacuum densities
    return std::complex<double>(ki_ * delta_k, ki_ * delta_k * CASS_CURVATURE);
}

// ============================================================
// calculate_THz_hole: Resonant THz Einstein Boson Bridge effect
// THz_hole = exp(I * 2PI * nu * d/c) / (1 + resonance * CURVATURE)
// This models the quantum tunneling bridge through ring gaps
// ============================================================
std::complex<double> UQFFCassiniCore::calculate_THz_hole(double nu, double distance) {
    double d_eff = (distance > 0) ? distance : 1.0;
    std::complex<double> phase = CASS_I * 2.0 * CASS_PI * nu * d_eff / CASS_C;
    double denominator = 1.0 + CASS_NU_THz * CASS_CURVATURE;
    return std::exp(phase) / denominator;
}

// ============================================================
// calculate_delta_v_particle: q-Scope Particle Deceleration
// delta_v = -K_Q * B_grad / (rho_vac_UA * c^2) × 1e-12
// ============================================================
std::complex<double> UQFFCassiniCore::calculate_delta_v_particle(double B_grad) {
    double dv_real = -CASS_K_Q * B_grad / (CASS_RHO_VAC_UA * CASS_C * CASS_C) * 1e-12;
    // Imaginary: quantum uncertainty in particle velocity measurement
    double dv_imag = CASS_K_Q * B_grad * CASS_RHO_VAC_UA * 1e-12;
    return std::complex<double>(dv_real, dv_imag);
}

// ============================================================
// calculate_U_g1: DPM Force (using Cassini ring geometry)
// ============================================================
std::complex<double> UQFFCassiniCore::calculate_U_g1(
    const std::vector<CassDPMVars>& vars, CassiniGeometryType geom) {
    std::complex<double> total(0.0, 0.0);
    for (const auto& v : vars) {
        double r_eff = (v.r > 0) ? v.r : 1.0;
        double f_UA_prime_re = (CASS_Z_MAX - v.Z) / CASS_Z_MAX;
        std::complex<double> f_UA_prime(f_UA_prime_re, f_UA_prime_re * CASS_CURVATURE);
        std::complex<double> f_nu = f_nu_THz(std::complex<double>(v.nu_THz, 0.0));

        if (geom == CASS_TOROIDAL) {
            total += k1_ * (f_UA_prime / (r_eff * r_eff)) * f_nu * 2.0 * CASS_PI * r_eff;
        } else {
            total += k1_ * (f_UA_prime / (r_eff * r_eff)) * f_nu;
        }
    }
    return total;
}

// ============================================================
// calculate_U_g3: Composite (U_i + U_m)
// ============================================================
std::complex<double> UQFFCassiniCore::calculate_U_g3(const CassDPMVars& v) {
    // U_Mi from current state
    std::complex<double> U_mi = calculate_U_Mi(v.r / CASS_C, v.r, 1); // t ~ r/c
    // U_Ii gyroscopic
    std::complex<double> U_ii = calculate_U_Ii(U_mi, 1.0);
    return k1_ * (U_ii + U_mi);
}

// ============================================================
// calculate_master_force: Sum of all UQFF Cassini forces
// F_master = U_g1 + U_g3 + U_Bi + THz_hole + delta_v
// ============================================================
std::complex<double> UQFFCassiniCore::calculate_master_force(
    const CassiniParams& params, const std::vector<CassDPMVars>& vars) {

    std::complex<double> U_g1 = calculate_U_g1(vars, params.geom);

    // Build a representative DPMVar for U_g3
    CassDPMVars v0;
    v0.Z = 26.0; v0.nu_THz = CASS_NU_THz; v0.r = params.ring_r;
    v0.B = params.B_field;
    std::complex<double> U_g3 = calculate_U_g3(v0);

    double delta_k = CASS_RHO_VAC_UA - CASS_RHO_VAC_SCM;
    std::complex<double> U_Bi = calculate_U_Bi(delta_k);

    // THz hole: ring width as distance
    std::complex<double> THzH = calculate_THz_hole(CASS_NU_THz, params.ring_r * 0.01);

    std::complex<double> dv = calculate_delta_v_particle(params.B_field);

    return U_g1 + U_g3 + U_Bi + THzH + dv;
}

// ============================================================
// UQFFCassiniSystem: Constructor and factory
// ============================================================
UQFFCassiniSystem::UQFFCassiniSystem(const CassiniParams& params) : params_(params) {
    // Build default DPMVars for Cassini ring geometry
    CassDPMVars v;
    v.Z = 14.0;             // Silicon-like grain composition
    v.nu_THz = CASS_NU_THz;
    v.theta = CASS_PI / 4.0;
    v.phi = 0.0;
    v.r = params.ring_r;
    v.r_shell = params.orbital_r;
    v.B = params.B_field;
    v.f_U_Bi = {0.0, 0.0};
    v.U_Ii = {0.0, 0.0};
    v.U_Mi = {0.0, 0.0};
    vars_.push_back(v);
}

std::complex<double> UQFFCassiniSystem::calculate_master_force(const UQFFCassiniCore& core) {
    return core.calculate_master_force(params_, vars_);
}

// ============================================================
// Factory: create_Cassini_system (Encke Gap parameters)
// ============================================================
UQFFCassiniSystem create_Cassini_system() {
    CassiniParams p;
    p.orbital_r       = 1.43e12;       // m (Cassini orbital distance from Saturn center)
    p.ring_r          = 7.0e7;         // m (ring radius scale)
    p.saturn_mass     = 5.683e26;      // kg
    p.ring_mass       = 1.5e19;        // kg (total ring mass estimate)
    p.B_field         = 1.0e-7;        // T (Saturn ring magnetic field)
    p.wind_vel        = 500.0;         // m/s (zonal wind speed)
    p.rotation_period = 10.7 * 3600.0; // s (Saturn rotation ~10.7 hours)
    p.geom            = CASS_TOROIDAL; // Ring geometry
    return UQFFCassiniSystem(p);
}
