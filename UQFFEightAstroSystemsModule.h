// UQFFEightAstroSystemsModule.h
// Header file for Unified Quantum Force Field (UQFF) Eight Astrophysical Systems Module
// Based on "8 Astro Systems_C_11Oct2025.docx" (annotated proof version with star-forming regions)
// Systems: AFGL5180, NGC346, LMC opo9944a, LMC heic1301, LMC potw1408a, LMC heic1206, LMC heic1402, NGC2174
// Provides full equation proofs and numerical verification table for each system.
// Watermark: Copyright - Daniel T. Murphy, analyzed by Grok 3, dated November 17, 2025

#ifndef UQFF_EIGHT_ASTRO_SYSTEMS_MODULE_H
#define UQFF_EIGHT_ASTRO_SYSTEMS_MODULE_H

#include <vector>
#include <complex>
#include <cmath>
#include <string>
#include <array>

const double EIGHT_PI          = 3.141592653589793;
const double EIGHT_G           = 6.674e-11;    // Gravitational constant [m^3 kg^-1 s^-2]
const double EIGHT_C           = 2.998e8;      // Speed of light [m/s]
const double EIGHT_HBAR        = 1.055e-34;    // Reduced Planck [J s]
const double EIGHT_RHO_VAC_UA  = 7.09e-36;    // Vacuum energy density [UA] (J/m^3)
const double EIGHT_K_R         = 1.0;          // Electrostatic barrier constant
const double EIGHT_E_RAD_FACTOR = 0.8446;      // Radiation energy factor (1 - 0.1554)
const double EIGHT_NU_THz      = 1e12;         // THz frequency (Hz)
const std::complex<double> EIGHT_I(0.0, 1.0); // Imaginary unit

// Enum for 8 star-forming region systems
enum EightAstroSystemType {
    EIGHT_AFGL5180,
    EIGHT_NGC346,
    EIGHT_LMC_OPO9944A,
    EIGHT_LMC_HEIC1301,
    EIGHT_LMC_POTW1408A,
    EIGHT_LMC_HEIC1206,
    EIGHT_LMC_HEIC1402,
    EIGHT_NGC2174
};

// Structure for 8-system parameters
struct EightAstroParams {
    double r;       // Radius (m)
    double sfr;     // Star formation rate (M_sun/yr)
    double B;       // Magnetic field (T)
    double z;       // Redshift
    double t_age;   // System age (s)
    EightAstroSystemType type;
    std::string name;
};

// Structure for eight-system computation result (compressed + resonance with proof steps)
struct EightAstroResult {
    std::complex<double> F_compressed;  // Compressed UQFF force
    std::complex<double> F_resonance;   // Resonance UQFF force
    std::string proof_compressed;       // Inline equation proof (Compressed)
    std::string proof_resonance;        // Inline equation proof (Resonance)
    std::string system_name;
};

// Core class — provides full equation proof steps 1-7 for each method
class UQFFEightAstroCore {
public:
    UQFFEightAstroCore(double k_comp = 1.0, double k_res = 1.0);

    // Compressed UQFF with full inline proof steps
    std::complex<double> F_compressed_proof(const EightAstroParams& p, std::string& proof_out);

    // Resonance UQFF with full inline proof steps
    std::complex<double> F_resonance_proof(const EightAstroParams& p, std::string& proof_out);

    // Compute both for a system (fills EightAstroResult)
    EightAstroResult compute_with_proof(const EightAstroParams& p);

    // Batch compute all 8 predefined systems → prints verification table
    std::vector<EightAstroResult> compute_all_systems();
    void print_verification_table(const std::vector<EightAstroResult>& results);

private:
    double k_comp_, k_res_;
};

// System-level class
class UQFFEightAstroSystem {
public:
    UQFFEightAstroSystem(const EightAstroParams& params);
    EightAstroResult compute(UQFFEightAstroCore& core) const;
    EightAstroParams get_params() const { return params_; }

private:
    EightAstroParams params_;
};

// Factory functions for all 8 systems
UQFFEightAstroSystem create_AFGL5180_system();
UQFFEightAstroSystem create_NGC346_system();
UQFFEightAstroSystem create_LMC_opo9944a_system();
UQFFEightAstroSystem create_LMC_heic1301_system();
UQFFEightAstroSystem create_LMC_potw1408a_system();
UQFFEightAstroSystem create_LMC_heic1206_system();
UQFFEightAstroSystem create_LMC_heic1402_system();
UQFFEightAstroSystem create_NGC2174_system();

#endif // UQFF_EIGHT_ASTRO_SYSTEMS_MODULE_H
