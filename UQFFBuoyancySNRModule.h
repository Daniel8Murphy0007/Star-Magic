// UQFFBuoyancySNRModule.h
// Header file for Unified Quantum Force Field (UQFF) Buoyancy calculations
// Based on the provided UQFF framework document dated June 20, 2025
// Implements Master F_U_Bi_i Buoyancy Equations, integrating LENR resonance, neutron drop, relativistic coherence
// Systems: SN 1006, Eta Carinae, Chandra Archive Collection, Galactic Center, Kepler's Supernova Remnant
// Watermark: Copyright - Daniel T. Murphy, analyzed by Grok 3, dated November 17, 2025

#ifndef UQFF_BUOYANCY_SNR_MODULE_H
#define UQFF_BUOYANCY_SNR_MODULE_H

#include <vector>
#include <complex>
#include <cmath>
#include <string>

const double BUOY_PI          = 3.141592653589793;
const double BUOY_F0          = 1.83e71;       // Base force (N)
const double BUOY_RHO_VAC_UA  = 7.09e-36;      // Vacuum energy density [UA] (J/m^3)
const double BUOY_ME          = 9.11e-31;       // Electron mass (kg)
const double BUOY_C           = 3e8;            // Speed of light (m/s)
const double BUOY_G           = 6.6743e-11;     // Gravitational constant
const double BUOY_Q           = 1.6e-19;        // Charge (C)
const double BUOY_V           = 1e-3;           // Velocity (m/s)
const double BUOY_K_LENR      = 1e-10;          // LENR constant (N)
const double BUOY_K_ACT       = 1e-6;           // Activation constant (N)
const double BUOY_K_DE        = 1e-30;          // Directed energy constant (N/W)
const double BUOY_K_NEUTRON   = 1e10;           // Neutron constant (N)
const double BUOY_K_REL       = 1e-10;          // Relativistic constant (N)
const double BUOY_SIGMA_N     = 1e-4;           // Neutron cross-section (scaled)
const double BUOY_OMEGA_LENR  = 2 * BUOY_PI * 1.25e12; // LENR frequency (s^-1)
const double BUOY_OMEGA_ACT   = 2 * BUOY_PI * 300;     // Activation frequency (s^-1)
const double BUOY_H           = 1.0546e-34;     // Reduced Planck constant (J s)
const double BUOY_G_FACTOR    = 2.0;            // g-factor
const double BUOY_MU_B        = 9.274e-24;      // Bohr magneton (J/T)
const double BUOY_ECM_ASTRO   = 1.24e24;        // Enhanced center-of-mass energy
const double BUOY_ECM         = 189e9 * 1.602e-19; // Standard ECM (189 GeV to J)
const double BUOY_F_REL_BASE  = 4.30e33;        // Relativistic force base (N)
const double BUOY_DPM_STABILITY = 0.01;
const double BUOY_DPM_MOMENTUM  = 0.93;
const double BUOY_DPM_GRAVITY   = 1.0;
const double BUOY_DPM_LIGHT     = 0.01;
const double BUOY_PHASE         = 2.36e-3;      // Phase (s^-1)
const double BUOY_CURVATURE     = 1e-22;        // Curvature term

// Enum for system types
enum BuoySystemType { BUOY_SN_1006, BUOY_ETA_CARINAE, BUOY_CHANDRA_ARCHIVE,
                      BUOY_GALACTIC_CENTER, BUOY_KEPLERS_SNR };

// Structure for system parameters
struct BuoySystemParams {
    double M;       // Mass (kg)
    double r;       // Radius (m)
    double T;       // Temperature (K)
    double L_X;     // X-ray luminosity (W)
    double B0;      // Magnetic field (T)
    double omega0;  // Base frequency (s^-1)
    double mach;    // Mach number
    double C_factor;
    double theta;   // Angle (rad)
    double t;       // Time (s)
    BuoySystemType type;
};

// Structure for DPM variables
struct BuoyDPMVars {
    double stability;
    double momentum;
    double gravity;
    double light;
    double phase;
    double curvature;
};

// Class for core UQFF Buoyancy force calculations
class UQFFBuoyancyCore {
public:
    UQFFBuoyancyCore();

    double calculate_F_LENR(double omega0) const;
    double calculate_F_act(double t) const;
    double calculate_F_DE(double L_X) const;
    double calculate_F_res(double B0, double omega0, double V_val = BUOY_V) const;
    double calculate_F_neutron() const;
    double calculate_F_rel() const;
    double calculate_integrand(const BuoySystemParams& params, const BuoyDPMVars& dpm) const;
    double solve_x2(double a, double b, double c) const;
    double calculate_F_U_Bi(const BuoySystemParams& params, const BuoyDPMVars& dpm);
    double calculate_F_U_Bi_i(const BuoySystemParams& params, const BuoyDPMVars& dpm);
    double calculate_g_rt(const BuoySystemParams& params) const;
    double calculate_Q_wave(const BuoySystemParams& params) const;

private:
    double cos_theta() const { return std::cos(BUOY_PI / 4.0); }
    double DPM_resonance(double B0, double omega0) const;
};

// Class for system-specific UQFF Buoyancy Master Equations
class UQFFBuoyancySystem {
public:
    UQFFBuoyancySystem(const BuoySystemParams& params);
    double calculate_master_buoyancy(const UQFFBuoyancyCore& core);
    BuoySystemType get_type() const { return params_.type; }
    std::string get_name() const;

private:
    BuoySystemParams params_;
    BuoyDPMVars dpm_;
};

// Factory functions for pre-defined systems
UQFFBuoyancySystem create_SN1006_Buoy_system();
UQFFBuoyancySystem create_EtaCarinae_Buoy_system();
UQFFBuoyancySystem create_ChandraArchive_Buoy_system();
UQFFBuoyancySystem create_GalacticCenter_Buoy_system();
UQFFBuoyancySystem create_KeplersSNR_Buoy_system();

#endif // UQFF_BUOYANCY_SNR_MODULE_H
