// UQFFConstants.h
// Shared physical constants for all UQFF modules.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef UQFF_CONSTANTS_H
#define UQFF_CONSTANTS_H

namespace UQFF {

// Gravitational
constexpr double G           = 6.6743e-11;    // m^3 kg^-1 s^-2
constexpr double c           = 3.0e8;         // m/s
constexpr double c2          = c * c;         // m^2/s^2

// Quantum
constexpr double hbar        = 1.0546e-34;    // J s
constexpr double h_planck    = 6.626e-34;     // J s
constexpr double l_planck    = 1.616e-35;     // m  (Planck length)
constexpr double t_planck    = 5.391e-44;     // s  (Planck time)
constexpr double m_proton    = 1.673e-27;     // kg
constexpr double m_electron  = 9.109e-31;     // kg
constexpr double q_e         = 1.602e-19;     // C  (elementary charge)

// Cosmological
constexpr double Lambda      = 1.1e-52;       // m^-2 (cosmological constant)
constexpr double H0          = 67.15;         // km/s/Mpc (Planck 2018)
constexpr double Omega_m     = 0.3;           // dimensionless
constexpr double Omega_Lambda = 0.7;          // dimensionless
constexpr double Omega_DM    = 0.268;         // DM fraction of critical density
constexpr double t_Hubble    = 13.8e9 * 3.156e7;   // s (age of universe)
constexpr double Mpc_to_m    = 3.086e22;      // m/Mpc

// Stellar
constexpr double M_sun       = 1.989e30;      // kg
constexpr double L_sun       = 3.828e26;      // W
constexpr double R_sun       = 6.957e8;       // m
constexpr double AU          = 1.496e11;      // m
constexpr double pc_to_m     = 3.086e16;      // m/pc
constexpr double year_to_s   = 3.156e7;       // s/yr
constexpr double Myr_s       = 3.156e13;      // s/Myr
constexpr double Gyr_s       = 3.156e16;      // s/Gyr

// Atomic
constexpr double a0          = 5.2918e-11;    // m (Bohr radius)
constexpr double m_H         = 1.6735e-27;    // kg (hydrogen atom mass)
constexpr double lambda_Halpha = 6.563e-7;    // m (H-alpha wavelength 656.3 nm)

// Gravitational wave
constexpr double h_strain_default = 1.0e-21;  // dimensionless (typical LIGO)
constexpr double lambda_gw_default = 1.0e13;  // m (GW wavelength)

// Magnetic critical
constexpr double B_crit_magnetar = 1.0e11;    // T (quantum critical field)

// Vacuum densities
constexpr double rho_vac_UA  = 7.09e-36;      // kg/m^3
constexpr double rho_vac_SCm = 7.09e-37;      // kg/m^3

// Mathematical
constexpr double pi          = 3.141592653589793;

} // namespace UQFF

#endif // UQFF_CONSTANTS_H
