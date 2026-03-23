// UQFFBuoyancySNRModule.cpp
// Implementation of UQFF Buoyancy force for Supernova Remnant (SNR) systems
// Based on "content(20).docx"
// Systems: SN1006, Eta Carinae, Chandra Archive, Galactic Center, Kepler's SNR
// Master equation: F_U_Bi_i = F_LENR + F_act + F_DE + F_res + F_neutron + F_rel
// Watermark: Copyright - Daniel T. Murphy, analyzed by Grok 3, dated November 17, 2025

#include "../UQFFBuoyancySNRModule.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>

// ============================================================
// UQFFBuoyancyCore: Constructor
// ============================================================
UQFFBuoyancyCore::UQFFBuoyancyCore(double kl, double ka, double kd, double kr, double kn, double krel)
    : k_lenr_(kl), k_act_(ka), k_de_(kd), k_res_(kr), k_neutron_(kn), k_rel_(krel) {}

// ============================================================
// F_LENR: Low-Energy Nuclear Reaction force
// F_LENR = k_LENR * rho_vac * sin(omega_LENR * t) * exp(-r / r0)
// ============================================================
double UQFFBuoyancyCore::calc_F_LENR(const BuoyancySNRParams& p) {
    double phase = BUOY_OMEGA_LENR * p.t_age;
    double decay = std::exp(-p.r / (1.0e15)); // 1 kpc normalization
    return k_lenr_ * BUOY_F0 * BUOY_RHO_VAC_UA * std::sin(phase) * decay;
}

// ============================================================
// F_act: Activation energy force
// F_act = k_act * F0 * sin(omega_ACT * t) * B * Q_wave
// ============================================================
double UQFFBuoyancyCore::calc_F_act(const BuoyancySNRParams& p) {
    double phase = BUOY_OMEGA_ACT * p.t_age;
    return k_act_ * BUOY_F0 * std::sin(phase) * p.B * p.Q_wave;
}

// ============================================================
// F_DE: Dark Energy pressure force
// F_DE = k_DE * RHO_VAC_UA * r^2 * (1 + z)
// ============================================================
double UQFFBuoyancyCore::calc_F_DE(const BuoyancySNRParams& p) {
    double h_z = 1.0 + p.z;
    return k_de_ * BUOY_RHO_VAC_UA * p.r * p.r * h_z;
}

// ============================================================
// F_res: Resonance force
// F_res = k_res * F0 * cos(omega_ACT * t) / r^2
// ============================================================
double UQFFBuoyancyCore::calc_F_res(const BuoyancySNRParams& p) {
    double r_eff = (p.r > 0) ? p.r : 1.0;
    double phase = BUOY_OMEGA_ACT * p.t_age;
    return k_res_ * BUOY_F0 * std::cos(phase) / (r_eff * r_eff);
}

// ============================================================
// F_neutron: Neutron production rate force contribution
// F_neutron = k_neutron * eta * m_n * c^2 / r^2
// eta = eta0 * B / rho_vac^0.5
// ============================================================
double UQFFBuoyancyCore::calc_F_neutron(const BuoyancySNRParams& p) {
    double r_eff = (p.r > 0) ? p.r : 1.0;
    double eta = p.eta_base * p.B / std::sqrt(BUOY_RHO_VAC_UA);
    double m_n = 1.675e-27; // neutron mass (kg)
    double E_n = m_n * BUOY_C * BUOY_C; // rest energy
    return k_neutron_ * eta * E_n / (r_eff * r_eff);
}

// ============================================================
// F_rel: Relativistic correction force
// F_rel = k_rel * F_REL_BASE * exp(-v_exp / c) * (1 + z)^2
// ============================================================
double UQFFBuoyancyCore::calc_F_rel(const BuoyancySNRParams& p) {
    double rel_factor = std::exp(-p.v_exp / BUOY_C);
    double z_factor = (1.0 + p.z) * (1.0 + p.z);
    return k_rel_ * BUOY_F_REL_BASE * rel_factor * z_factor;
}

// ============================================================
// F_U_Bi_i: Master buoyancy force (sum of all components)
// ============================================================
double UQFFBuoyancyCore::F_U_Bi_i(const BuoyancySNRParams& p) {
    return calc_F_LENR(p)
         + calc_F_act(p)
         + calc_F_DE(p)
         + calc_F_res(p)
         + calc_F_neutron(p)
         + calc_F_rel(p);
}

// ============================================================
// solve_x2: Quadratic solver for buoyancy mode separation
// Given a*x^2 + b*x + c = 0, returns both roots
// Used to find dual buoyancy modes for SNR expansion
// ============================================================
std::pair<double, double> UQFFBuoyancyCore::solve_x2(double a, double b, double c) {
    if (std::fabs(a) < 1e-50) {
        throw std::runtime_error("solve_x2: coefficient a is too small (near zero)");
    }
    double disc = b * b - 4.0 * a * c;
    if (disc < 0.0) {
        // Return real parts of complex roots
        double re = -b / (2.0 * a);
        return {re, re}; // degenerate complex pair
    }
    return {(-b + std::sqrt(disc)) / (2.0 * a),
            (-b - std::sqrt(disc)) / (2.0 * a)};
}

// ============================================================
// UQFFBuoyancySystem: Constructor and calculation
// ============================================================
UQFFBuoyancySystem::UQFFBuoyancySystem(const BuoyancySNRParams& params) : params_(params) {}

BuoyancySNRResult UQFFBuoyancySystem::compute(const UQFFBuoyancyCore& core) {
    BuoyancySNRResult result;
    result.system_name = params_.name;
    result.F_total = core.F_U_Bi_i(params_);
    result.F_LENR = core.calc_F_LENR(params_);
    result.F_act = core.calc_F_act(params_);
    result.F_DE = core.calc_F_DE(params_);
    result.F_res = core.calc_F_res(params_);
    result.F_neutron = core.calc_F_neutron(params_);
    result.F_rel = core.calc_F_rel(params_);
    // g_rt placeholder: acceleration at r
    result.g_rt = (params_.r > 0) ? result.F_total / (1.989e30 * params_.Q_wave) : 0.0;
    return result;
}

// ============================================================
// Factory Functions for 5 SNR systems
// ============================================================

// SN 1006: Type Ia SNR, distance ~2.2 kpc
UQFFBuoyancySystem create_SN1006_system() {
    BuoyancySNRParams p = {
        1.989e31,   // mass (kg) ~10 solar masses ejecta
        6.17e16,    // r (m) ~ 2 kpc radius
        1e6,        // v_exp (m/s) expansion velocity ~1000 km/s
        1e32,       // L (W) luminosity
        1e-5,       // B (T)
        1e-12,      // eta_base
        1.0,        // z (redshift ~0)
        1.0,        // Q_wave
        BUOY_PI/4,  // theta
        3.213e10,   // t_age (s)
        "SN1006"
    };
    return UQFFBuoyancySystem(p);
}

// Eta Carinae: LBV / SNR precursor, distance ~2.3 kpc
UQFFBuoyancySystem create_EtaCarinae_system() {
    BuoyancySNRParams p = {
        3.978e32,   // mass (~200 solar masses)
        3.09e16,    // r ~1 kpc
        5e5,        // v_exp ~500 km/s wind
        5e31,       // L
        1e-4,       // B (strong field ~mG)
        1e-11,      // eta_base
        0.0,        // z
        1.2,        // Q_wave
        BUOY_PI/3,
        5.05e14,    // t_age ~16,000 yr
        "Eta Carinae"
    };
    return UQFFBuoyancySystem(p);
}

// Chandra Archive Collection: composite of X-ray SNR measurements
UQFFBuoyancySystem create_ChandraArchive_system() {
    BuoyancySNRParams p = {
        5.967e31,   // mass (~30 solar)
        3.09e19,    // r ~ 1 Mpc scale (collection)
        2e5,        // v_exp
        1e33,       // L
        1e-5,       // B
        5e-12,      // eta_base
        0.1,        // z
        0.8,        // Q_wave
        BUOY_PI/6,
        1.892e16,   // t_age
        "Chandra Archive Collection"
    };
    return UQFFBuoyancySystem(p);
}

// Galactic Center: dense SNR population near Sgr A*
UQFFBuoyancySystem create_GalacticCenter_system() {
    BuoyancySNRParams p = {
        1.989e36,   // mass (~10^6 M_sun effective)
        2.461e20,   // r ~ 8 kpc Galactic Center distance
        1e5,        // v_exp
        3e37,       // L
        1e-4,       // B (strong near GC)
        1e-10,      // eta_base
        0.0,        // z
        1.5,        // Q_wave
        BUOY_PI/2,
        1.577e17,   // t_age ~5 Gyr
        "Galactic Center"
    };
    return UQFFBuoyancySystem(p);
}

// Kepler's SNR (SN1604): distance ~6 kpc, Type Ia
UQFFBuoyancySystem create_KeplersSNR_system() {
    BuoyancySNRParams p = {
        1.989e31,   // mass ~10 M_sun ejecta
        1.852e17,   // r ~ 6 kpc
        7e3,        // v_exp ~7000 km/s
        1e31,       // L
        1e-5,       // B
        8e-13,      // eta_base
        0.0,        // z
        1.0,        // Q_wave
        BUOY_PI/4,
        1.293e10,   // t_age ~410 yr
        "Kepler's SNR (SN1604)"
    };
    return UQFFBuoyancySystem(p);
}
