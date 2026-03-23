// UQFFCalculationsModule.h
// Header file for Unified Quantum Force Field (UQFF) calculations
// Based on the provided UQFF framework document dated June 05, 2025
// Implements key equations for U_g1, U_g3, U_m, E, eta, and system-specific Master UQFF equations
// Systems: Messier 82 (M82), Spirograph Nebula IC418, Canis Major R136, NGC 6302 Butterfly Nebula, NGC 7027
// Watermark: Copyright - Daniel T. Murphy, analyzed by Grok 3, dated November 10, 2025

#ifndef UQFF_CALCULATIONS_MODULE_H
#define UQFF_CALCULATIONS_MODULE_H

#include <vector>
#include <complex>
#include <cmath>
#include <string>

// Constants (scaled as per document)
const double UQFFCALC_PI         = 3.141592653589793;
const double UQFFCALC_K_R        = 1.0;      // Electrostatic barrier constant
const double UQFFCALC_Z_MAX      = 1000.0;   // Max Z for f_UA' and f_SCm
const double UQFFCALC_GAMMA      = 1.0;      // Decay constant for U_m
const double UQFFCALC_MU_0       = 4 * UQFFCALC_PI * 1e-7; // Vacuum permeability
const double UQFFCALC_RHO_VAC_UA = 1e-27;    // Vacuum energy density [UA]
const double UQFFCALC_RHO_VAC_SCM= 1e-28;    // Vacuum energy density [SCm]
const double UQFFCALC_K_ETA_BASE = 2.75e8;   // Base neutron production rate constant
const double UQFFCALC_SSQ        = 1.0;      // [SSq] placeholder
const double UQFFCALC_N_QUANTUM  = 26.0;     // 26 quantum states

// Forward declarations
class UQFFSystem;

// Enum for geometries (spherical, toroidal)
enum UQFFCalcGeometryType { UQFFCALC_SPHERICAL, UQFFCALC_TOROIDAL };

// Structure for DPM variables
struct UQFFCalcDPMVars {
    double f_UA_prime;  // f_UA' = (Z_max - Z) / Z_max
    double f_SCm;       // f_SCm = Z / Z_max
    double R_EB;        // R_EB = k_R * Z
    double Z;           // Atomic number
    double nu_THz;      // THz frequency
    double nu_res;      // Resonance frequency
    double theta;       // Polar angle
    double phi;         // Azimuthal angle
    double r;           // Distance
    double r_shell;     // Shell distance
    double f_Ub;        // Buoyancy factor (calibration difference)
};

// Class for core UQFF force calculations
class UQFFCore {
public:
    UQFFCore(double k1 = 1.0, double ki = 1.0, double km = 1.0, double ke = 1.0,
             double k_eta = UQFFCALC_K_ETA_BASE);

    // U_g1 (DPM) force calculation
    double calculate_U_g1(const std::vector<UQFFCalcDPMVars>& vars,
                          UQFFCalcGeometryType geom = UQFFCALC_SPHERICAL);

    // U_g3 (U_i + U_m) force calculation
    double calculate_U_g3(const UQFFCalcDPMVars& vars);

    // Universal Magnetism U_m
    double calculate_U_m(double t, double r, int n,
                         double rho_vac_SCm = UQFFCALC_RHO_VAC_SCM,
                         double mu_j = UQFFCALC_MU_0);

    // Electric Field E
    double calculate_E(double U_m_val, double r,
                       double rho_vac_UA = UQFFCALC_RHO_VAC_UA);

    // Neutron Production Rate eta
    double calculate_eta(double n, double t, double U_m_val,
                         double rho_vac_UA = UQFFCALC_RHO_VAC_UA);

    double calculate_f_Ub(double calibration_diff) const { return calibration_diff; }

private:
    double k1_, ki_, km_, ke_, k_eta_;
    double f_nu_THz(double nu) const {
        return 1.0 + std::sin(UQFFCALC_PI * nu / 1e12);
    }
    double heaviside(double x) const { return (x >= 0) ? 1.0 : 0.0; }
    double quasi_factor(double x) const { return 1.0 + 1e-13 * x; }
};

// Class for system-specific UQFF Master Equations
class UQFFSystem {
public:
    UQFFSystem(const std::string& name, double sfr = 1.0,
               double wind_vel = 1000.0, double mag_field = 1e-4);

    double calculate_master_force(const std::vector<UQFFCalcDPMVars>& vars,
                                  const UQFFCore& core);

    void set_f_Ub_scale(double scale) { f_Ub_scale_ = scale; }
    std::string get_name() const { return name_; }

private:
    std::string name_;
    double sfr_;
    double wind_vel_;
    double mag_field_;
    double f_Ub_scale_;
};

// Factory functions for pre-defined systems
UQFFSystem create_M82_UQFFCalc_system();
UQFFSystem create_IC418_UQFFCalc_system();
UQFFSystem create_CanisMajor_UQFFCalc_system();
UQFFSystem create_NGC6302_UQFFCalc_system();
UQFFSystem create_NGC7027_UQFFCalc_system();

#endif // UQFF_CALCULATIONS_MODULE_H
