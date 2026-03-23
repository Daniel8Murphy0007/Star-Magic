// UQFFMultiAstroSystemsModule.h
// Header file for Unified Quantum Force Field (UQFF) Multi-Astrophysical Systems
// Based on "8 Astro Systems_B_11Oct2025.docx" (final 11-system version)
// Implements simultaneous Compressed, Resonance, and Buoyancy solutions for 11 systems:
// NGC4826, NGC1805, NGC6307, NGC7027, Cassini Encke/Division/Maxwell, ESO391-12, M57, LMC, ESO510-G13
// Watermark: Copyright - Daniel T. Murphy, analyzed by Grok 3, dated November 17, 2025

#ifndef UQFF_MULTI_ASTRO_SYSTEMS_MODULE_H
#define UQFF_MULTI_ASTRO_SYSTEMS_MODULE_H

#include <vector>
#include <complex>
#include <cmath>
#include <string>
#include <map>

const double MULTI_PI          = 3.141592653589793;
const double MULTI_G           = 6.674e-11;    // Gravitational constant [m^3 kg^-1 s^-2]
const double MULTI_C           = 2.998e8;      // Speed of light [m/s]
const double MULTI_HBAR        = 1.055e-34;    // Reduced Planck [J s]
const double MULTI_RHO_VAC_UA  = 7.09e-36;    // Vacuum energy density [UA] (J/m^3)
const double MULTI_K_R         = 1.0;          // Electrostatic barrier constant
const double MULTI_E_RAD_FACTOR = 0.8446;      // 1 - 0.1554 = radiation energy factor
const double MULTI_H_Z_BASE    = 67.4;         // Hubble constant (km/s/Mpc)
const std::complex<double> MULTI_I(0.0, 1.0);  // Imaginary unit

// Enum for 11 multi-astro systems
enum MultiAstroSystemType {
    MULTI_NGC4826,
    MULTI_NGC1805,
    MULTI_NGC6307,
    MULTI_NGC7027,
    MULTI_CASSINI_ENCKE,
    MULTI_CASSINI_DIV,
    MULTI_CASSINI_MAX,
    MULTI_ESO391_12,
    MULTI_M57,
    MULTI_LMC,
    MULTI_ESO510_G13
};

// Structure for multi-system parameters
struct MultiAstroParams {
    double r;       // Radius or orbital distance (m)
    double sfr;     // Star formation rate (M_sun/yr)
    double B;       // Magnetic field (T)
    double z;       // Redshift
    double t_age;   // System age (s)
    MultiAstroSystemType type;
    std::string name;
};

// Structure for the 3-solution result (Compressed, Resonance, Buoyancy)
struct MultiAstroResult {
    std::complex<double> F_compressed;  // Compressed UQFF force
    std::complex<double> F_resonance;   // Resonance UQFF force
    std::complex<double> F_buoyancy;    // Buoyancy-dominant force
    double F_DPM_creation;              // DPM pair creation rate
    std::string system_name;
};

// Core calculation class for multi-astro systems
class UQFFMultiAstroCore {
public:
    UQFFMultiAstroCore(double k_comp = 1.0, double k_res = 1.0, double k_buoy = 1.0);

    // Compressed UQFF: Hubble-corrected, radiation-damped DPM force
    std::complex<double> F_compressed(const MultiAstroParams& p);

    // Resonance UQFF: omega resonance with THz hole term
    std::complex<double> F_resonance(const MultiAstroParams& p);

    // Buoyancy-dominant UQFF: vacuum energy + magnetism + redshift correction
    std::complex<double> F_buoyancy(const MultiAstroParams& p);

    // DPM pair creation simulation
    double F_DPM_creation(const MultiAstroParams& p);

    // Solve all three simultaneously for a system
    MultiAstroResult compute_system(const MultiAstroParams& p);

    // Batch compute all 11 predefined systems
    std::vector<MultiAstroResult> compute_all_systems();

private:
    double k_comp_, k_res_, k_buoy_;
    double hubble_correction(double z) const { return 1.0 + z; }
};

// System-level class holding params and calling core
class UQFFMultiAstroSystem {
public:
    UQFFMultiAstroSystem(const MultiAstroParams& params);
    MultiAstroResult compute(UQFFMultiAstroCore& core) const;
    MultiAstroParams get_params() const { return params_; }

private:
    MultiAstroParams params_;
};

// Factory functions for all 11 systems
UQFFMultiAstroSystem create_NGC4826_system();
UQFFMultiAstroSystem create_NGC1805_system();
UQFFMultiAstroSystem create_NGC6307_system();
UQFFMultiAstroSystem create_NGC7027_system();
UQFFMultiAstroSystem create_Cassini_Encke_system();
UQFFMultiAstroSystem create_Cassini_Div_system();
UQFFMultiAstroSystem create_Cassini_Max_system();
UQFFMultiAstroSystem create_ESO391_12_system();
UQFFMultiAstroSystem create_M57_system();
UQFFMultiAstroSystem create_LMC_system();
UQFFMultiAstroSystem create_ESO510_G13_system();

#endif // UQFF_MULTI_ASTRO_SYSTEMS_MODULE_H
