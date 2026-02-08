// source82_wolfram.cpp
// Wolfram Language Physics Term Companions for SMBHUQFFModule (source82.cpp)
// Implements 10 PhysicsTerm classes for Virgo Cluster UQFF Integration
// Systems: Virgo Cluster (nearest large cluster, ~16.5 Mpc), M87 central galaxy, SMBH M-σ relation
// Auto-generated: February 6, 2026
// Module: SMBHUQFFModule - Master Universal Gravity Equation for SMBH M-σ relation
// Physics: Cluster dynamics, ICM, dark matter halo, galaxy velocity dispersion, X-ray emission
// Key parameters: M_cluster~1e15 M_sun, R_virial~2.2 Mpc, σ_v~700 km/s, T_ICM~2.3 keV
// Copyright - Daniel T. Murphy

#include <cmath>
#include <string>
#include <map>
#include <complex>
#include <vector>
#include <memory>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Physical constants
constexpr double G_CONST = 6.6743e-11;      // Gravitational constant (m³/kg·s²)
constexpr double M_SUN = 1.989e30;          // Solar mass (kg)
constexpr double MPC_TO_M = 3.086e22;       // Megaparsec to meters
constexpr double KPC_TO_M = 3.086e19;       // Kiloparsec to meters
constexpr double KEV_TO_J = 1.602e-16;      // keV to Joules
constexpr double K_BOLTZ = 1.381e-23;       // Boltzmann constant (J/K)
constexpr double M_PROTON = 1.673e-27;      // Proton mass (kg)
constexpr double C_LIGHT = 2.998e8;         // Speed of light (m/s)
constexpr double YEAR_TO_S = 3.156e7;       // Year to seconds

// ========================================
// CLASS 820: VirgoClusterMassTerm
// Category: gravity
// Physics: Total cluster mass M_cluster ~ 1.2e15 M_sun (gravitational + dark matter)
// Virgo Cluster is the nearest large galaxy cluster at ~16.5 Mpc
// ========================================
class VirgoClusterMassTerm
{
private:
    double G;              // Gravitational constant (m³/kg·s²)
    double M_cluster;      // Total cluster mass (kg, ~1.2e15 M_sun)
    double R_virial;       // Virial radius (m, ~2.2 Mpc = 6.79e22 m)
    double z;              // Redshift (0.0036 for Virgo center)

public:
    VirgoClusterMassTerm(double mass = 1.2e15 * M_SUN, double r_vir = 2.2 * MPC_TO_M, double redshift = 0.0036)
        : G(G_CONST), M_cluster(mass), R_virial(r_vir), z(redshift) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // Protect against r = 0 (or extremely small r) to avoid division by zero / numerical instability
        double r_eff = r;
        if (r_eff < 1e-10)
            r_eff = 1e-10;

        // Enclosed mass profile: M(<r) = M_cluster * (r/R_virial)^3 * (1 + (R_virial/r))^(-2)
        // Simplified NFW-like profile for cluster
        double x = r_eff / R_virial;
        double M_enclosed = M_cluster * (x * x * x) / (1.0 + x) / (1.0 + x);
        
        // Return gravitational acceleration: a = G * M_enclosed / r_eff²
        return (G * M_enclosed) / (r_eff * r_eff);
    }

    std::string toWolfram() const
    {
        return "VirgoClusterMass[r_, Mcluster_: 1.2*^15 * 1.989*^30, Rvirial_: 2.2 * 3.086*^22, G_: 6.6743*^-11] := "
               "Module[{x, Menclosed}, "
               "x = r / Rvirial; "
               "Menclosed = Mcluster * x^3 / (1 + x)^2; "
               "G * Menclosed / r^2]";
    }

    std::string getSignature() const { return "VirgoClusterMassTerm(r, params)"; }
    std::string getCategory() const { return "gravity"; }
    std::string getName() const { return "VirgoClusterMass"; }
    std::string getDescription() const { return "Virgo Cluster total mass gravitational acceleration: a = G·M(<r)/r² with M_cluster~1.2e15 M_sun"; }
};

// ========================================
// CLASS 821: VirgoClusterIntraclusterMediumTerm
// Category: gas_dynamics
// Physics: Intracluster medium (ICM) - hot X-ray emitting gas at T~2-4 keV
// ICM contains ~15% of cluster baryonic mass, creates X-ray emission and SZ effect
// ========================================
class VirgoClusterIntraclusterMediumTerm
{
private:
    double T_ICM;          // ICM temperature (keV, ~2.3 keV for Virgo)
    double n_e0;           // Central electron density (m⁻³, ~3e3)
    double r_c;            // Core radius (m, ~40 kpc = 1.23e21 m)
    double beta;           // Beta model parameter (~0.5 for Virgo)

public:
    VirgoClusterIntraclusterMediumTerm(double temp = 2.3, double ne_central = 3e3, double core_r = 40 * KPC_TO_M, double beta_param = 0.5)
        : T_ICM(temp), n_e0(ne_central), r_c(core_r), beta(beta_param) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // Beta model: n_e(r) = n_e0 * (1 + (r/r_c)²)^(-3β/2)
        double n_e = n_e0 * std::pow(1.0 + (r * r) / (r_c * r_c), -1.5 * beta);
        
        // Pressure: P = n_e * k_B * T (convert T from keV to Kelvin: T_K = T_keV * 1.16e7)
        double T_K = T_ICM * 1.16e7;
        double P_ICM = n_e * K_BOLTZ * T_K;
        
        return P_ICM;
    }

    std::string toWolfram() const
    {
        return "VirgoICM[r_, T_: 2.3, ne0_: 3000, rc_: 40 * 3.086*^19, beta_: 0.5] := "
               "Module[{ne, TK, kB}, "
               "kB = 1.381*^-23; "
               "ne = ne0 * (1 + (r/rc)^2)^(-1.5 * beta); "
               "TK = T * 1.16*^7; "
               "ne * kB * TK]";
    }

    std::string getSignature() const { return "VirgoClusterIntraclusterMediumTerm(r, params)"; }
    std::string getCategory() const { return "gas_dynamics"; }
    std::string getName() const { return "VirgoClusterICM"; }
    std::string getDescription() const { return "ICM beta-model pressure: P = n_e(r)·k_B·T with T~2.3 keV, β~0.5"; }
};

// ========================================
// CLASS 822: VirgoClusterGravitationalPotentialTerm
// Category: gravity
// Physics: Cluster gravitational potential well Φ(r) = -G·M(<r)/r
// Determines escape velocity and binding energy of member galaxies
// ========================================
class VirgoClusterGravitationalPotentialTerm
{
private:
    double G;              // Gravitational constant
    double M_cluster;      // Total cluster mass
    double R_virial;       // Virial radius
    double c_NFW;          // NFW concentration parameter (~4 for clusters)

public:
    VirgoClusterGravitationalPotentialTerm(double mass = 1.2e15 * M_SUN, double r_vir = 2.2 * MPC_TO_M, double concentration = 4.0)
        : G(G_CONST), M_cluster(mass), R_virial(r_vir), c_NFW(concentration) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // NFW potential: Φ(r) = -G·M_virial · ln(1 + c·x) / (x · f(c))
        // where x = r/R_virial, f(c) = ln(1+c) - c/(1+c)
        double x = r / R_virial;
        double f_c = std::log(1.0 + c_NFW) - c_NFW / (1.0 + c_NFW);
        
        // Avoid division by zero at r=0
        if (x < 1e-10) x = 1e-10;
        
        double Phi = -G * M_cluster * std::log(1.0 + c_NFW * x) / (x * f_c * R_virial);
        
        return Phi;
    }

    std::string toWolfram() const
    {
        return "VirgoClusterPotential[r_, Mcluster_: 1.2*^15 * 1.989*^30, Rvirial_: 2.2 * 3.086*^22, c_: 4] := "
               "Module[{x, fc, G}, "
               "G = 6.6743*^-11; "
               "x = Max[r / Rvirial, 10^-10]; "
               "fc = Log[1 + c] - c/(1 + c); "
               "-G * Mcluster * Log[1 + c*x] / (x * fc * Rvirial)]";
    }

    std::string getSignature() const { return "VirgoClusterGravitationalPotentialTerm(r, params)"; }
    std::string getCategory() const { return "gravity"; }
    std::string getName() const { return "VirgoClusterPotential"; }
    std::string getDescription() const { return "NFW gravitational potential: Φ(r) = -G·M·ln(1+c·x)/(x·f(c)·R_vir) with c~4"; }
};

// ========================================
// CLASS 823: VirgoClusterDarkMatterTerm
// Category: dark_matter
// Physics: Dark matter halo profile (NFW) - ~85% of cluster mass
// ρ_DM(r) = ρ_s / ((r/r_s)(1 + r/r_s)²)
// ========================================
class VirgoClusterDarkMatterTerm
{
private:
    double rho_s;          // Scale density (kg/m³, ~3e-23 for Virgo)
    double r_s;            // Scale radius (m, ~550 kpc = 1.70e22 m)
    double M_DM;           // Total dark matter mass (~1e15 M_sun)

public:
    VirgoClusterDarkMatterTerm(double scale_density = 3e-23, double scale_radius = 550 * KPC_TO_M, double dm_mass = 1e15 * M_SUN)
        : rho_s(scale_density), r_s(scale_radius), M_DM(dm_mass) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // NFW density profile: ρ(r) = ρ_s / ((r/r_s) * (1 + r/r_s)²)
        double x = r / r_s;
        if (x < 1e-10) x = 1e-10;  // Avoid singularity at r=0
        
        double rho_DM = rho_s / (x * (1.0 + x) * (1.0 + x));
        
        return rho_DM;
    }

    std::string toWolfram() const
    {
        return "VirgoClusterDarkMatter[r_, rhoS_: 3*^-23, rS_: 550 * 3.086*^19] := "
               "Module[{x}, "
               "x = Max[r/rS, 10^-10]; "
               "rhoS / (x * (1 + x)^2)]";
    }

    std::string getSignature() const { return "VirgoClusterDarkMatterTerm(r, params)"; }
    std::string getCategory() const { return "dark_matter"; }
    std::string getName() const { return "VirgoClusterDarkMatter"; }
    std::string getDescription() const { return "NFW dark matter profile: ρ(r) = ρ_s/((r/r_s)·(1+r/r_s)²) with r_s~550 kpc"; }
};

// ========================================
// CLASS 824: VirgoClusterM87JetTerm
// Category: agn
// Physics: M87 (Virgo A) central AGN relativistic jet energy injection
// Jet power ~1e44 erg/s = 1e37 W, affects ICM heating and cluster evolution
// ========================================
class VirgoClusterM87JetTerm
{
private:
    double L_jet;          // Jet luminosity (W, ~1e37 for M87)
    double theta_jet;      // Jet opening angle (rad, ~0.1)
    double v_jet;          // Jet velocity (m/s, ~0.99c)
    double r_jet_base;     // Jet base radius (m, ~0.01 pc = 3.086e14 m)

public:
    VirgoClusterM87JetTerm(double luminosity = 1e37, double angle = 0.1, double velocity = 0.99 * C_LIGHT, double base_r = 0.01 * 3.086e16)
        : L_jet(luminosity), theta_jet(angle), v_jet(velocity), r_jet_base(base_r) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // Jet energy density: u_jet(r) = L_jet / (4π·r²·v_jet·Ω_jet)
        // where Ω_jet = 2π(1 - cos(θ_jet)) is the solid angle
        double Omega_jet = 2.0 * M_PI * (1.0 - std::cos(theta_jet));
        
        // Energy density in jet cone
        double u_jet = L_jet / (4.0 * M_PI * r * r * v_jet * (Omega_jet / (4.0 * M_PI)));
        
        // Lorentz factor for relativistic correction
        double gamma = 1.0 / std::sqrt(1.0 - (v_jet * v_jet) / (C_LIGHT * C_LIGHT));
        
        return u_jet * gamma;
    }

    std::string toWolfram() const
    {
        return "VirgoM87Jet[r_, Ljet_: 10^37, theta_: 0.1, vjet_: 0.99 * 2.998*^8] := "
               "Module[{OmegaJet, uJet, gamma, c}, "
               "c = 2.998*^8; "
               "OmegaJet = 2 * Pi * (1 - Cos[theta]); "
               "uJet = Ljet / (4 * Pi * r^2 * vjet * (OmegaJet/(4*Pi))); "
               "gamma = 1/Sqrt[1 - (vjet/c)^2]; "
               "uJet * gamma]";
    }

    std::string getSignature() const { return "VirgoClusterM87JetTerm(r, params)"; }
    std::string getCategory() const { return "agn"; }
    std::string getName() const { return "VirgoM87Jet"; }
    std::string getDescription() const { return "M87 AGN jet energy density: u_jet = L_jet·γ/(4π·r²·v_jet·Ω) with L~10³⁷ W"; }
};

// ========================================
// CLASS 825: VirgoClusterTidalStrippingTerm
// Category: dynamics
// Physics: Tidal stripping of infalling galaxies by cluster potential
// Tidal radius r_t = r·(M_gal/(3·M_cluster(<r)))^(1/3)
// ========================================
class VirgoClusterTidalStrippingTerm
{
private:
    double M_cluster;      // Cluster mass (kg)
    double R_virial;       // Virial radius (m)
    double G;              // Gravitational constant

public:
    VirgoClusterTidalStrippingTerm(double mass = 1.2e15 * M_SUN, double r_vir = 2.2 * MPC_TO_M)
        : M_cluster(mass), R_virial(r_vir), G(G_CONST) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // Get galaxy mass from params, default to 1e11 M_sun (typical spiral)
        double M_gal = (params.count("M_gal") ? params.at("M_gal") : 1e11 * M_SUN);
        
        // Enclosed cluster mass at radius r (NFW-like approximation)
        double x = r / R_virial;
        double M_enclosed = M_cluster * (x * x * x) / (1.0 + x) / (1.0 + x);
        
        // Tidal radius: r_t = r * (M_gal / (3 * M_enclosed))^(1/3)
        double r_tidal = r * std::pow(M_gal / (3.0 * M_enclosed), 1.0 / 3.0);
        
        // Tidal acceleration at galaxy edge
        double a_tidal = 2.0 * G * M_enclosed * r_tidal / (r * r * r);
        
        return a_tidal;
    }

    std::string toWolfram() const
    {
        return "VirgoTidalStripping[r_, Mgal_, Mcluster_: 1.2*^15 * 1.989*^30, Rvirial_: 2.2 * 3.086*^22] := "
               "Module[{x, Menclosed, rtidal, G}, "
               "G = 6.6743*^-11; "
               "x = r / Rvirial; "
               "Menclosed = Mcluster * x^3 / (1 + x)^2; "
               "rtidal = r * (Mgal / (3 * Menclosed))^(1/3); "
               "2 * G * Menclosed * rtidal / r^3]";
    }

    std::string getSignature() const { return "VirgoClusterTidalStrippingTerm(r, params)"; }
    std::string getCategory() const { return "dynamics"; }
    std::string getName() const { return "VirgoTidalStripping"; }
    std::string getDescription() const { return "Tidal stripping acceleration: a_tidal = 2·G·M(<r)·r_t/r³ with r_t from Jacobi radius"; }
};

// ========================================
// CLASS 826: VirgoClusterVirialTerm
// Category: thermodynamics
// Physics: Virial theorem equilibrium: 2K + U = 0, σ_v² = G·M_vir/(3·R_vir)
// Velocity dispersion σ_v ~ 700 km/s for Virgo
// ========================================
class VirgoClusterVirialTerm
{
private:
    double M_virial;       // Virial mass (kg)
    double R_virial;       // Virial radius (m)
    double sigma_v;        // Velocity dispersion (m/s, ~700 km/s)
    double G;              // Gravitational constant

public:
    VirgoClusterVirialTerm(double mass = 1.2e15 * M_SUN, double r_vir = 2.2 * MPC_TO_M, double sigma = 700e3)
        : M_virial(mass), R_virial(r_vir), sigma_v(sigma), G(G_CONST) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // r parameter included for API consistency (unused for virial equilibrium)
        (void)r;
        
        // Virial relation: σ² = G·M_vir / (3·R_vir)
        // Returns virial equilibrium parameter (should be ~1 if in equilibrium)
        double sigma_virial = std::sqrt(G * M_virial / (3.0 * R_virial));
        
        // Ratio of observed to virial velocity dispersion
        double virial_ratio = sigma_v / sigma_virial;
        
        return virial_ratio;
    }

    // Helper method (extends base interface for specific virial mass calculations)
    double computeVirialMass() const
    {
        // M_vir = 3·σ²·R_vir / G
        return 3.0 * sigma_v * sigma_v * R_virial / G;
    }

    std::string toWolfram() const
    {
        return "VirgoClusterVirial[Mvir_: 1.2*^15 * 1.989*^30, Rvir_: 2.2 * 3.086*^22, sigmaV_: 700000] := "
               "Module[{sigmaVirial, G}, "
               "G = 6.6743*^-11; "
               "sigmaVirial = Sqrt[G * Mvir / (3 * Rvir)]; "
               "sigmaV / sigmaVirial]";
    }

    std::string getSignature() const { return "VirgoClusterVirialTerm(r, params)"; }
    std::string getCategory() const { return "thermodynamics"; }
    std::string getName() const { return "VirgoClusterVirial"; }
    std::string getDescription() const { return "Virial equilibrium ratio: σ_obs/σ_vir where σ_vir² = G·M/(3·R), σ~700 km/s"; }
};

// ========================================
// CLASS 827: VirgoClusterXRayLuminosityTerm
// Category: radiation
// Physics: X-ray emission from hot ICM: L_X ∝ n_e²·Λ(T)·V
// Virgo X-ray luminosity ~10^43 erg/s = 10^36 W
// ========================================
class VirgoClusterXRayLuminosityTerm
{
private:
    double n_e0;           // Central electron density (m⁻³)
    double T_ICM;          // ICM temperature (keV)
    double r_c;            // Core radius (m)
    double beta;           // Beta model parameter

public:
    VirgoClusterXRayLuminosityTerm(double ne_central = 3e3, double temp = 2.3, double core_r = 40 * KPC_TO_M, double beta_param = 0.5)
        : n_e0(ne_central), T_ICM(temp), r_c(core_r), beta(beta_param) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // Electron density profile (beta model)
        double n_e = n_e0 * std::pow(1.0 + (r * r) / (r_c * r_c), -1.5 * beta);
        
        // Cooling function Λ(T) ~ 3e-23 * sqrt(T_keV) [W·m³] (simplified)
        double Lambda_T = 3e-23 * std::sqrt(T_ICM);
        
        // X-ray emissivity: ε_X = n_e² * Λ(T) [W/m³]
        double epsilon_X = n_e * n_e * Lambda_T;
        
        return epsilon_X;
    }

    double computeTotalLuminosity(double R_max) const
    {
        // Numerically integrate emissivity over spherical volume:
        // L_X = ∫_0^{R_max} ε_X(r) 4π r² dr
        if (R_max <= 0.0)
        {
            return 0.0;
        }

        const int N = 1000; // number of radial integration steps (trapezoidal rule)
        const double dr = R_max / static_cast<double>(N);
        double integral = 0.0;

        // Empty parameter map since compute(r, params) does not currently use params
        const std::map<std::string, double> emptyParams;

        for (int i = 0; i <= N; ++i)
        {
            double r = dr * static_cast<double>(i);
            double epsilon_X = compute(r, emptyParams); // emissivity at radius r
            double volumeElement = 4.0 * M_PI * r * r;

            double integrand = epsilon_X * volumeElement;

            // Trapezoidal weighting: 0.5 at endpoints, 1.0 elsewhere
            if (i == 0 || i == N)
            {
                integral += 0.5 * integrand;
            }
            else
            {
                integral += integrand;
            }
        }

        double L_X = integral * dr;
        return L_X;
    }

    std::string toWolfram() const
    {
        return "VirgoXRay[r_, ne0_: 3000, T_: 2.3, rc_: 40 * 3.086*^19, beta_: 0.5] := "
               "Module[{ne, LambdaT}, "
               "ne = ne0 * (1 + (r/rc)^2)^(-1.5 * beta); "
               "LambdaT = 3*^-23 * Sqrt[T]; "
               "ne^2 * LambdaT]";
    }

    std::string getSignature() const { return "VirgoClusterXRayLuminosityTerm(r, params)"; }
    std::string getCategory() const { return "radiation"; }
    std::string getName() const { return "VirgoXRay"; }
    std::string getDescription() const { return "X-ray emissivity: ε_X = n_e²·Λ(T) with Λ(T)~3e-23·√T W·m³"; }
};

// ========================================
// CLASS 828: VirgoClusterVelocityDispersionTerm
// Category: kinematics
// Physics: Galaxy velocity dispersion profile σ(r)
// Central σ ~ 700 km/s, decreases with radius
// ========================================
class VirgoClusterVelocityDispersionTerm
{
private:
    double sigma_0;        // Central velocity dispersion (m/s, ~700 km/s)
    double r_sigma;        // Scale radius for dispersion (m, ~500 kpc)
    double G;              // Gravitational constant
    double M_cluster;      // Cluster mass

public:
    VirgoClusterVelocityDispersionTerm(double sigma_central = 700e3, double scale_r = 500 * KPC_TO_M, double mass = 1.2e15 * M_SUN)
        : sigma_0(sigma_central), r_sigma(scale_r), G(G_CONST), M_cluster(mass) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // Velocity dispersion profile: σ(r) = σ_0 / sqrt(1 + (r/r_σ)²)
        double sigma_r = sigma_0 / std::sqrt(1.0 + (r * r) / (r_sigma * r_sigma));
        
        return sigma_r;
    }

    double computeDynamicalMass(double r) const
    {
        // M(<r) = r·σ²(r) / G (simplified spherical Jeans equation)
        double sigma_r = sigma_0 / std::sqrt(1.0 + (r * r) / (r_sigma * r_sigma));
        return r * sigma_r * sigma_r / G;
    }

    std::string toWolfram() const
    {
        return "VirgoVelocityDispersion[r_, sigma0_: 700000, rSigma_: 500 * 3.086*^19] := "
               "sigma0 / Sqrt[1 + (r/rSigma)^2]";
    }

    std::string getSignature() const { return "VirgoClusterVelocityDispersionTerm(r, params)"; }
    std::string getCategory() const { return "kinematics"; }
    std::string getName() const { return "VirgoVelocityDispersion"; }
    std::string getDescription() const { return "Velocity dispersion profile: σ(r) = σ_0/√(1+(r/r_σ)²) with σ_0~700 km/s"; }
};

// ========================================
// CLASS 829: SMBHMSigmaRelationTerm
// Category: scaling_relations
// Physics: M-σ relation from source82 SMBH module
// M_BH = 1.9e8 * (σ/200 km/s)^4.38 M_sun (McConnell & Ma 2013)
// Links SMBH mass to host galaxy velocity dispersion
// ========================================
class SMBHMSigmaRelationTerm
{
private:
    double alpha;          // M-σ slope (~4.38)
    double M_norm;         // Normalization mass (M_sun)
    double sigma_norm;     // Normalization dispersion (m/s, 200 km/s)

public:
    SMBHMSigmaRelationTerm(double slope = 4.38, double m_0 = 1.9e8, double sig_0 = 200e3)
        : alpha(slope), M_norm(m_0), sigma_norm(sig_0) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // For M-σ relation, r is interpreted as sigma (velocity dispersion in m/s)
        // This maintains API consistency while using physically meaningful parameter
        double sigma = r;  // r represents sigma for this scaling relation
        
        // Allow override from params map if provided
        if (params.count("sigma")) sigma = params.at("sigma");
        
        // M_BH = M_norm * (σ/σ_norm)^α
        double M_BH = M_norm * std::pow(sigma / sigma_norm, alpha) * M_SUN;
        
        return M_BH;
    }

    // Helper method (extends base interface for inverse calculation)
    double computeSigmaFromMass(double M_BH) const
    {
        // Inverse: σ = σ_norm * (M_BH / M_norm)^(1/α)
        return sigma_norm * std::pow(M_BH / (M_norm * M_SUN), 1.0 / alpha);
    }

    std::string toWolfram() const
    {
        return "SMBHMSigma[sigma_, alpha_: 4.38, Mnorm_: 1.9*^8, sigmaNorm_: 200000] := "
               "Mnorm * 1.989*^30 * (sigma/sigmaNorm)^alpha";
    }

    std::string getSignature() const { return "SMBHMSigmaRelationTerm(r, params)"; }
    std::string getCategory() const { return "scaling_relations"; }
    std::string getName() const { return "SMBHMSigma"; }
    std::string getDescription() const { return "M-σ relation: M_BH = 1.9e8·(σ/200 km/s)^4.38 M_sun (McConnell & Ma 2013)"; }
};

// ========================================
// REGISTRATION FUNCTION
// ========================================

// NOTE: Registration function is commented out because PhysicsTermRegistry
// is defined in MAIN_1_CoAnQi.cpp and requires proper header inclusion.
// To integrate these terms into the main calculator:
// 1. Include this file in MAIN_1_CoAnQi.cpp
// 2. Uncomment the registration function below
// 3. Call registerWolframTerms_source82(registry) from main initialization

// Forward declaration of PhysicsTermRegistry if not available
// #include "PhysicsTermRegistry.h"

/*
void registerWolframTerms_source82(PhysicsTermRegistry& registry) {
    // Register all 10 Virgo Cluster physics terms
    registry.registerPhysicsTerm("VirgoClusterMass", std::make_unique<VirgoClusterMassTerm>(), "wolfram");
    registry.registerPhysicsTerm("VirgoClusterICM", std::make_unique<VirgoClusterIntraclusterMediumTerm>(), "wolfram");
    registry.registerPhysicsTerm("VirgoClusterPotential", std::make_unique<VirgoClusterGravitationalPotentialTerm>(), "wolfram");
    registry.registerPhysicsTerm("VirgoClusterDarkMatter", std::make_unique<VirgoClusterDarkMatterTerm>(), "wolfram");
    registry.registerPhysicsTerm("VirgoM87Jet", std::make_unique<VirgoClusterM87JetTerm>(), "wolfram");
    registry.registerPhysicsTerm("VirgoTidalStripping", std::make_unique<VirgoClusterTidalStrippingTerm>(), "wolfram");
    registry.registerPhysicsTerm("VirgoClusterVirial", std::make_unique<VirgoClusterVirialTerm>(), "wolfram");
    registry.registerPhysicsTerm("VirgoXRay", std::make_unique<VirgoClusterXRayLuminosityTerm>(), "wolfram");
    registry.registerPhysicsTerm("VirgoVelocityDispersion", std::make_unique<VirgoClusterVelocityDispersionTerm>(), "wolfram");
    registry.registerPhysicsTerm("SMBHMSigma", std::make_unique<SMBHMSigmaRelationTerm>(), "wolfram");
}
*/

// ============================================================================
// SUMMARY: 10 PHYSICS TERM CLASSES FOR VIRGO CLUSTER (820-829)
// ============================================================================
// 820: VirgoClusterMassTerm - Total cluster mass gravitational acceleration
// 821: VirgoClusterIntraclusterMediumTerm - ICM beta-model pressure profile
// 822: VirgoClusterGravitationalPotentialTerm - NFW gravitational potential
// 823: VirgoClusterDarkMatterTerm - NFW dark matter density profile
// 824: VirgoClusterM87JetTerm - M87 AGN relativistic jet energy injection
// 825: VirgoClusterTidalStrippingTerm - Tidal stripping of infalling galaxies
// 826: VirgoClusterVirialTerm - Virial equilibrium σ² = GM/(3R)
// 827: VirgoClusterXRayLuminosityTerm - X-ray emissivity from hot ICM
// 828: VirgoClusterVelocityDispersionTerm - Galaxy velocity dispersion profile
// 829: SMBHMSigmaRelationTerm - M-σ relation linking SMBH to host σ
// ============================================================================
// Key Virgo Cluster Parameters:
// - Distance: 16.5 Mpc (z = 0.0036)
// - Total mass: M_cluster ~ 1.2e15 M_sun
// - Virial radius: R_vir ~ 2.2 Mpc
// - Velocity dispersion: σ_v ~ 700 km/s
// - ICM temperature: T ~ 2.3 keV
// - Central galaxy: M87 (Virgo A) with M_BH ~ 6.5e9 M_sun
// - Number of galaxies: ~1,500-2,000 members
// ============================================================================
// Companion to source82.cpp (SMBHUQFFModule)
// Copyright - Daniel T. Murphy, 2025-2026
// ============================================================================
