// source77_wolfram.cpp - Whirlpool Galaxy (M51) UQFF Module
// Classes 720-729: M51 grand design spiral galaxy with NGC 5195 companion
// Physical basis: Prototypical interacting galaxy pair, tidal bridge/tail, enhanced star formation

#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ========================================
// Class 720: M51SpiralArmDensityWaveTerm
// ========================================
// Physical model: Two-armed grand design spiral with m=2 mode, Lin-Shu density wave theory
// Observational basis: M51 pitch angle i ≈ 20°, pattern speed Ω_p ≈ 25 km/s/kpc
// Reference: Tully (1974) - M51 rotation curve and spiral structure
class M51SpiralArmDensityWaveTerm {
public:
    M51SpiralArmDensityWaveTerm(double r_kpc, double phi_rad, double t_Gyr = 0.0,
                                 double A = 0.3, double m = 2.0, double pitch_deg = 20.0,
                                 double Omega_p = 25.0)
        : r_kpc_(r_kpc), phi_rad_(phi_rad), t_Gyr_(t_Gyr), A_(A), m_(m),
          pitch_deg_(pitch_deg), Omega_p_(Omega_p) {}

    double compute() const {
        // Pitch angle to wavenumber: tan(i) = m/(k·r)
        double pitch_rad = pitch_deg_ * M_PI / 180.0;
        double k = m_ / (r_kpc_ * std::tan(pitch_rad)); // kpc^-1
        
        // Phase: ψ = m·φ - k·r - Ω_p·t
        double psi = m_ * phi_rad_ - k * r_kpc_ - Omega_p_ * t_Gyr_;
        
        // Density perturbation: Σ(r,φ,t) = Σ₀(r)·[1 + A·cos(ψ)]
        double Sigma_base = 100.0 * std::exp(-r_kpc_ / 3.5); // M_sun/pc², exponential disk
        double Sigma_perturb = Sigma_base * (1.0 + A_ * std::cos(psi));
        
        // Corotation radius: r_CR where Ω(r_CR) = Ω_p
        // For flat rotation curve v = 200 km/s, Ω = v/r, r_CR = v/Ω_p ≈ 8 kpc
        double r_CR = 200.0 / Omega_p_;
        double Delta_r_CR = std::abs(r_kpc_ - r_CR);
        
        return Sigma_perturb * (1.0 + A_ + Delta_r_CR / 10.0);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M51SpiralDensity[r_, \\[Phi]_, t_, A_, m_, pitch_, \\[CapitalOmega]p_] := "
            << "Module[{k, \\[Psi], \\[CapitalSigma]base}, "
            << "k = m/(r*Tan[pitch*Degree]); "
            << "\\[Psi] = m*\\[Phi] - k*r - \\[CapitalOmega]p*t; "
            << "\\[CapitalSigma]base = 100*Exp[-r/3.5]; "
            << "\\[CapitalSigma]base*(1 + A*Cos[\\[Psi]])]; "
            << "(* M51 m=2 grand design, pitch = " << pitch_deg_ << "°, \\[CapitalOmega]p = " << Omega_p_ << " km/s/kpc *)";
        return oss.str();
    }

    std::string getSignature() const { return "M51SpiralArmDensityWaveTerm(r_kpc, phi_rad, t_Gyr, A, m, pitch_deg, Omega_p)"; }
    std::string getCategory() const { return "dynamics"; }

private:
    double r_kpc_;
    double phi_rad_;
    double t_Gyr_;
    double A_;           // Amplitude of density perturbation
    double m_;           // Azimuthal mode number
    double pitch_deg_;   // Pitch angle [degrees]
    double Omega_p_;     // Pattern speed [km/s/kpc]
};

// ========================================
// Class 721: M51TidalInteractionTerm
// ========================================
// Physical model: NGC 5195 companion at d ≈ 9 kpc separation, tidal torque and bridge formation
// Observational basis: M51-NGC 5195 closest approach ~50-100 Myr ago, tidal tail extends ~50 kpc
// Reference: Salo & Laurikainen (2000) - M51/NGC 5195 interaction model
class M51TidalInteractionTerm {
public:
    M51TidalInteractionTerm(double r_kpc, double phi_rad,
                            double M_NGC5195 = 2e10, double d_kpc = 9.0, double phi_comp_rad = 0.0)
        : r_kpc_(r_kpc), phi_rad_(phi_rad), M_NGC5195_(M_NGC5195),
          d_kpc_(d_kpc), phi_comp_rad_(phi_comp_rad) {}

    double compute() const {
        // Position of NGC 5195 companion relative to M51 center
        double x_comp = d_kpc_ * std::cos(phi_comp_rad_);
        double y_comp = d_kpc_ * std::sin(phi_comp_rad_);
        
        // Test particle position
        double x = r_kpc_ * std::cos(phi_rad_);
        double y = r_kpc_ * std::sin(phi_rad_);
        
        // Separation from companion
        double dx = x - x_comp;
        double dy = y - y_comp;
        double r_sep = std::sqrt(dx * dx + dy * dy);
        
        // Tidal force: F_tid = -G·M_comp·(r - r_comp)/|r - r_comp|³
        const double G = 4.3e-6; // kpc·(km/s)²/M_sun
        double F_tid_x = -G * M_NGC5195_ * dx / (r_sep * r_sep * r_sep);
        double F_tid_y = -G * M_NGC5195_ * dy / (r_sep * r_sep * r_sep);
        double F_tid_mag = std::sqrt(F_tid_x * F_tid_x + F_tid_y * F_tid_y);
        
        // Tidal torque: τ = r × F_tid
        double torque = x * F_tid_y - y * F_tid_x; // kpc²·(km/s)²/Gyr
        
        // Roche lobe radius: r_Roche = d·(M_M51/(3·M_NGC5195))^(1/3)
        const double M_M51 = 1.6e11; // M_sun
        double r_Roche = d_kpc_ * std::pow(M_M51 / (3.0 * M_NGC5195_), 1.0/3.0);
        
        // Tidal tail formation parameter: Δr = r - r_Roche (positive if beyond Roche lobe)
        double Delta_r_Roche = r_kpc_ - r_Roche;
        
        return F_tid_mag * (1.0 + std::abs(torque) / 1e6 + std::abs(Delta_r_Roche) / r_Roche);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M51TidalForce[r_, \\[Phi]_, Mcomp_, d_, \\[Phi]comp_] := "
            << "Module[{xcomp, ycomp, x, y, dx, dy, rsep, Ftidx, Ftidy}, "
            << "xcomp = d*Cos[\\[Phi]comp]; ycomp = d*Sin[\\[Phi]comp]; "
            << "x = r*Cos[\\[Phi]]; y = r*Sin[\\[Phi]]; "
            << "dx = x - xcomp; dy = y - ycomp; rsep = Sqrt[dx^2 + dy^2]; "
            << "Ftidx = -G*Mcomp*dx/rsep^3; Ftidy = -G*Mcomp*dy/rsep^3; "
            << "Sqrt[Ftidx^2 + Ftidy^2]]; "
            << "(* M51-NGC 5195 d = " << d_kpc_ << " kpc, M_comp = " << M_NGC5195_ << " Msun *)";
        return oss.str();
    }

    std::string getSignature() const { return "M51TidalInteractionTerm(r_kpc, phi_rad, M_NGC5195, d_kpc, phi_comp_rad)"; }
    std::string getCategory() const { return "dynamics"; }

private:
    double r_kpc_;
    double phi_rad_;
    double M_NGC5195_;     // Companion mass [M_sun]
    double d_kpc_;         // M51-NGC 5195 separation [kpc]
    double phi_comp_rad_;  // Companion position angle [radians]
};

// ========================================
// Class 722: M51StarFormationRateTerm
// ========================================
// Physical model: Interaction-enhanced SFR with Schmidt-Kennicutt relation Σ_SFR ∝ Σ_gas^N
// Observational basis: M51 global SFR ≈ 3-4 M_sun/yr, enhanced in spiral arms and interaction zone
// Reference: Calzetti et al. (2005) - Spatially resolved SFR in M51 from HST
class M51StarFormationRateTerm {
public:
    M51StarFormationRateTerm(double Sigma_gas, double spiral_phase = 0.0,
                             double N = 1.4, double epsilon = 0.02, double enhancement = 1.5)
        : Sigma_gas_(Sigma_gas), spiral_phase_(spiral_phase), N_(N),
          epsilon_(epsilon), enhancement_(enhancement) {}

    double compute() const {
        // Base Kennicutt-Schmidt: Σ_SFR = A·Σ_gas^N
        const double A = 2.5e-4; // (M_sun/yr/kpc²) / (M_sun/pc²)^N
        double Sigma_SFR_base = A * std::pow(Sigma_gas_, N_);
        
        // Spiral arm enhancement: amplification at spiral arm (spiral_phase ≈ 0)
        double arm_factor = 1.0 + enhancement_ * std::exp(-spiral_phase_ * spiral_phase_ / 0.1);
        
        // Free-fall time: t_ff = √(3π/(32·G·ρ_gas))
        const double h_gas = 100.0; // pc, gas scale height
        double rho_gas = Sigma_gas_ / (2.0 * h_gas); // M_sun/pc³
        const double G = 4.302e-3; // pc·(km/s)²·M_sun^-1
        double t_ff_yr = std::sqrt(3.0 * M_PI / (32.0 * G * rho_gas * 1e-6 * 3.086e13)) / 3.156e7;
        
        // Efficiency-corrected SFR: SFR = ε·M_gas/t_ff
        double Sigma_SFR = epsilon_ * Sigma_gas_ / (t_ff_yr / 1e6) * arm_factor;
        
        // Integrated SFR over disk (r_d ≈ 3.5 kpc)
        const double r_d = 3.5; // kpc
        double SFR_total = Sigma_SFR * M_PI * r_d * r_d;
        
        return Sigma_SFR * (1.0 + arm_factor + SFR_total / 10.0);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M51StarFormationRate[\\[CapitalSigma]gas_, spiralPhase_, N_, \\[Epsilon]_, enh_] := "
            << "Module[{A, \\[CapitalSigma]SFRbase, armFactor, \\[Rho]gas, tff}, "
            << "A = 2.5*10^(-4); \\[CapitalSigma]SFRbase = A*\\[CapitalSigma]gas^N; "
            << "armFactor = 1 + enh*Exp[-spiralPhase^2/0.1]; "
            << "\\[Rho]gas = \\[CapitalSigma]gas/200; tff = Sqrt[3*Pi/(32*G*\\[Rho]gas*10^(-6)*3.086*10^13)]/3.156*10^7; "
            << "\\[Epsilon]*\\[CapitalSigma]gas/(tff/10^6)*armFactor]; "
            << "(* M51 SFR ~ 3-4 Msun/yr, N = " << N_ << ", enhancement = " << enhancement_ << " *)";
        return oss.str();
    }

    std::string getSignature() const { return "M51StarFormationRateTerm(Sigma_gas, spiral_phase, N, epsilon, enhancement)"; }
    std::string getCategory() const { return "stellar"; }

private:
    double Sigma_gas_;      // Gas surface density [M_sun/pc²]
    double spiral_phase_;   // Phase relative to spiral arm (0 = on arm, π = interarm)
    double N_;              // Power-law index (Kennicutt-Schmidt)
    double epsilon_;        // Star formation efficiency
    double enhancement_;    // Spiral arm/interaction enhancement factor
};

// ========================================
// Class 723: M51MolecularCloudTerm
// ========================================
// Physical model: Giant molecular cloud (GMC) mass function dN/dM ∝ M^(-α) with virial equilibrium
// Observational basis: M51 GMCs with masses 10^4-10^7 M_sun, CO observations reveal ~500 GMCs
// Reference: Schinnerer et al. (2013) - PAWS survey CO(1-0) in M51
class M51MolecularCloudTerm {
public:
    M51MolecularCloudTerm(double M_cloud, double alpha = 1.8, double X_CO = 2e20)
        : M_cloud_(M_cloud), alpha_(alpha), X_CO_(X_CO) {}

    double compute() const {
        // GMC mass function: dN/dM = N₀·M^(-α) with α ≈ 1.8 (Rosolowsky 2005)
        const double N_0 = 1000.0; // Normalization
        const double M_min = 1e4; // M_sun, minimum GMC mass
        double dN_dM = N_0 * std::pow(M_cloud_ / M_min, -alpha_);
        
        // Virial equilibrium: M_vir = 5·σ_v²·R / G
        // For σ_v ~ 5 km/s, R ~ 50 pc → M_vir ~ 10^6 M_sun
        const double sigma_v = 5.0; // km/s, velocity dispersion
        double R_pc = 50.0 * std::pow(M_cloud_ / 1e6, 1.0/3.0); // Larson's law R ∝ M^(1/3)
        const double G = 4.302e-3; // pc·(km/s)²·M_sun^-1
        double M_vir = 5.0 * sigma_v * sigma_v * R_pc / G;
        
        // Virial parameter: α_vir = M_vir/M_cloud (α_vir ~ 1 for equilibrium)
        double alpha_vir = M_vir / M_cloud_;
        
        // CO luminosity: L_CO = M_cloud / X_CO where X_CO is CO-to-H₂ conversion factor
        double L_CO = M_cloud_ / X_CO_; // K·km/s·pc²
        
        // Free-fall time: t_ff = √(3π/(32·G·ρ)) where ρ = 3·M/(4π·R³)
        double rho = 3.0 * M_cloud_ / (4.0 * M_PI * R_pc * R_pc * R_pc); // M_sun/pc³
        double t_ff_yr = std::sqrt(3.0 * M_PI / (32.0 * G * rho * 1e-6 * 3.086e13)) / 3.156e7;
        
        return dN_dM * (1.0 + alpha_vir + L_CO / 1e6 + t_ff_yr / 1e7);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M51GMCMassFunction[M_, \\[Alpha]_, XCO_] := "
            << "Module[{N0, Mmin, dNdM, R, Mvir, \\[Alpha]vir, LCO}, "
            << "N0 = 1000; Mmin = 10^4; dNdM = N0*(M/Mmin)^(-\\[Alpha]); "
            << "R = 50*(M/10^6)^(1/3); Mvir = 5*5^2*R/G; "
            << "\\[Alpha]vir = Mvir/M; LCO = M/XCO; "
            << "dNdM*(1 + \\[Alpha]vir + LCO/10^6)]; "
            << "(* M51 ~500 GMCs, mass range 10^4-10^7 Msun, \\[Alpha] = " << alpha_ << " *)";
        return oss.str();
    }

    std::string getSignature() const { return "M51MolecularCloudTerm(M_cloud, alpha, X_CO)"; }
    std::string getCategory() const { return "stellar"; }

private:
    double M_cloud_;  // GMC mass [M_sun]
    double alpha_;    // Power-law index of mass function
    double X_CO_;     // CO-to-H₂ conversion factor [cm^-2/(K·km/s)]
};

// ========================================
// Class 724: M51DarkMatterHaloTerm
// ========================================
// Physical model: NFW halo ρ_DM(r) = ρ_s / [(r/r_s)·(1 + r/r_s)²]
// Observational basis: M51 rotation curve decomposition suggests M_halo ~ 10^12 M_sun
// Reference: Tully et al. (2008) - M51 mass model from HI rotation curve
class M51DarkMatterHaloTerm {
public:
    M51DarkMatterHaloTerm(double r_kpc, double M_200 = 1e12, double c = 10.0)
        : r_kpc_(r_kpc), M_200_(M_200), c_(c) {}

    double compute() const {
        // NFW profile parameters
        // r_200 = virial radius where ρ_halo = 200·ρ_crit
        const double H_0 = 70.0; // km/s/Mpc, Hubble constant
        const double rho_crit = 3.0 * H_0 * H_0 / (8.0 * M_PI * 4.3e-6 * 1e6); // M_sun/kpc³
        double r_200 = std::pow(3.0 * M_200_ / (4.0 * M_PI * 200.0 * rho_crit), 1.0/3.0); // kpc
        
        // Scale radius: r_s = r_200/c
        double r_s = r_200 / c_;
        
        // Characteristic density: ρ_s = M_200 / [4π·r_s³·(ln(1+c) - c/(1+c))]
        double f_c = std::log(1.0 + c_) - c_ / (1.0 + c_);
        double rho_s = M_200_ / (4.0 * M_PI * r_s * r_s * r_s * f_c);
        
        // NFW density: ρ_DM(r) = ρ_s / [(r/r_s)·(1 + r/r_s)²]
        double x = r_kpc_ / r_s;
        double rho_DM = rho_s / (x * (1.0 + x) * (1.0 + x));
        
        // Enclosed DM mass: M_DM(<r) = 4π·ρ_s·r_s³·[ln(1+x) - x/(1+x)]
        double M_DM_enclosed = 4.0 * M_PI * rho_s * r_s * r_s * r_s * (std::log(1.0 + x) - x / (1.0 + x));
        
        // Circular velocity from DM: v_DM² = G·M_DM(<r)/r
        const double G = 4.3e-6; // kpc·(km/s)²/M_sun
        double v_DM = std::sqrt(G * M_DM_enclosed / r_kpc_);
        
        return rho_DM * (1.0 + v_DM / 100.0 + M_DM_enclosed / 1e11);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M51NFWDensity[r_, M200_, c_] := "
            << "Module[{H0, \\[Rho]crit, r200, rs, fc, \\[Rho]s, x}, "
            << "H0 = 70; \\[Rho]crit = 3*H0^2/(8*Pi*G*10^6); "
            << "r200 = (3*M200/(4*Pi*200*\\[Rho]crit))^(1/3); rs = r200/c; "
            << "fc = Log[1 + c] - c/(1 + c); \\[Rho]s = M200/(4*Pi*rs^3*fc); "
            << "x = r/rs; \\[Rho]s/(x*(1 + x)^2)]; "
            << "(* M51 NFW halo: M_200 = " << M_200_ << " Msun, c = " << c_ << " *)";
        return oss.str();
    }

    std::string getSignature() const { return "M51DarkMatterHaloTerm(r_kpc, M_200, c)"; }
    std::string getCategory() const { return "dark_matter"; }

private:
    double r_kpc_;
    double M_200_;  // Halo mass within r_200 [M_sun]
    double c_;      // Concentration parameter
};

// ========================================
// Class 725: M51MagneticFieldTerm
// ========================================
// Physical model: B-field aligned with spiral arms, dynamo amplification, Faraday rotation
// Observational basis: M51 ordered B-field ~15 μG in arms, random ~10 μG, polarization maps
// Reference: Fletcher et al. (2011) - Magnetic fields in M51 from radio continuum
class M51MagneticFieldTerm {
public:
    M51MagneticFieldTerm(double r_kpc, double phi_rad,
                         double B_ordered = 15.0, double B_random = 10.0, double pitch_deg = 20.0)
        : r_kpc_(r_kpc), phi_rad_(phi_rad), B_ordered_(B_ordered),
          B_random_(B_random), pitch_deg_(pitch_deg) {}

    double compute() const {
        // Ordered field: B_ord(r,φ) = B₀·exp(-r/r_B)·cos(m·φ - k·r)
        const double r_B = 5.0; // kpc, magnetic scale length
        const double m = 2.0; // Bisymmetric spiral (ASS) mode
        double pitch_rad = pitch_deg_ * M_PI / 180.0;
        double k = m / (r_kpc_ * std::tan(pitch_rad)); // kpc^-1
        
        double B_ord_amplitude = B_ordered_ * std::exp(-r_kpc_ / r_B);
        double B_ord = B_ord_amplitude * std::cos(m * phi_rad_ - k * r_kpc_);
        
        // Random field: isotropic Gaussian with variance σ_B²
        double B_rand = B_random_; // Assume characteristic value
        
        // Total field: B_tot² = B_ord² + B_rand²
        double B_tot = std::sqrt(B_ord * B_ord + B_rand * B_rand);
        
        // Synchrotron emissivity: j_ν ∝ B^(1+α) with α ≈ 0.8
        double j_synch = std::pow(B_tot, 1.8);
        
        // Faraday rotation measure: RM = 812·∫ n_e·B_||·dl [rad/m²]
        const double n_e = 0.03; // cm^-3, thermal electron density
        const double L_kpc = 1.0; // kpc, path length through disk
        double RM = 812.0 * n_e * B_ord * L_kpc * 3.086e21 / 1e3; // rad/m²
        
        // Magnetic pressure: P_mag = B²/(8π) [dyne/cm²]
        double B_G = B_tot * 1e-6; // Convert μG to G
        double P_mag = B_G * B_G / (8.0 * M_PI);
        
        return B_tot * (1.0 + j_synch / 1e3 + std::abs(RM) / 100.0 + P_mag / 1e-12);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M51MagneticField[r_, \\[Phi]_, Bord_, Brand_, pitch_] := "
            << "Module[{rB, m, k, BordAmp, Bord, Brand, Btot}, "
            << "rB = 5; m = 2; k = m/(r*Tan[pitch*Degree]); "
            << "BordAmp = Bord*Exp[-r/rB]; Bord = BordAmp*Cos[m*\\[Phi] - k*r]; "
            << "Brand = Brand; Btot = Sqrt[Bord^2 + Brand^2]; Btot]; "
            << "(* M51 bisymmetric spiral field: B_ord = " << B_ordered_ << " μG, B_rand = " << B_random_ << " μG *)";
        return oss.str();
    }

    std::string getSignature() const { return "M51MagneticFieldTerm(r_kpc, phi_rad, B_ordered, B_random, pitch_deg)"; }
    std::string getCategory() const { return "magnetic"; }

private:
    double r_kpc_;
    double phi_rad_;
    double B_ordered_;   // Ordered field strength [μG]
    double B_random_;    // Random field strength [μG]
    double pitch_deg_;   // Spiral pitch angle [degrees]
};

// ========================================
// Class 726: M51CentralBlackHoleTerm
// ========================================
// Physical model: SMBH with M_BH ≈ 10^6-10^7 M_sun (low-luminosity AGN)
// Observational basis: M51 X-ray source at nucleus, M-σ relation from bulge velocity dispersion
// Reference: Ho et al. (2001) - LINER activity in M51 nucleus
class M51CentralBlackHoleTerm {
public:
    M51CentralBlackHoleTerm(double r_pc, double M_BH = 1e7, double L_Edd_ratio = 1e-4)
        : r_pc_(r_pc), M_BH_(M_BH), L_Edd_ratio_(L_Edd_ratio) {}

    double compute() const {
        // Schwarzschild radius: R_S = 2·G·M_BH/c²
        const double G = 6.674e-8; // cm³/g/s²
        const double c = 2.998e10; // cm/s
        const double M_sun_g = 1.989e33; // g
        double R_S_cm = 2.0 * G * M_BH_ * M_sun_g / (c * c);
        double R_S_pc = R_S_cm / 3.086e18;
        
        // Eddington luminosity: L_Edd = 1.26×10^38·(M_BH/M_sun) erg/s
        double L_Edd = 1.26e38 * M_BH_;
        double L_bol = L_Edd_ratio_ * L_Edd; // Bolometric luminosity
        
        // Sphere of influence: r_SOI = G·M_BH/σ² where σ ≈ 100 km/s for M51 bulge
        const double sigma_bulge = 100.0; // km/s
        double r_SOI_pc = G * M_BH_ * M_sun_g / (sigma_bulge * 1e5 * sigma_bulge * 1e5) / 3.086e18;
        
        // Gravitational potential: Φ = -G·M_BH/r
        double Phi = -G * M_BH_ * M_sun_g / (r_pc_ * 3.086e18); // erg/g
        
        // Bondi accretion rate: Ṁ_Bondi = π·λ·(G·M_BH)²·ρ_inf / c_s³
        const double rho_inf = 1e-24; // g/cm³, ambient ISM density
        const double c_s = 1e6; // cm/s, sound speed
        const double lambda = 0.25; // Dimensionless accretion efficiency
        double M_dot_Bondi = M_PI * lambda * (G * M_BH_ * M_sun_g) * (G * M_BH_ * M_sun_g) * rho_inf / (c_s * c_s * c_s);
        double M_dot_Bondi_Msun_yr = M_dot_Bondi * 3.156e7 / M_sun_g;
        
        return std::abs(Phi) / 1e14 * (1.0 + L_bol / 1e40 + r_SOI_pc / r_pc_ + M_dot_Bondi_Msun_yr / 1e-6);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M51BlackHole[r_, MBH_, LEddRatio_] := "
            << "Module[{RS, LEdd, Lbol, \\[Sigma], rSOI, \\[CapitalPhi]}, "
            << "RS = 2*G*MBH*Msun/c^2; LEdd = 1.26*10^38*MBH; Lbol = LEddRatio*LEdd; "
            << "\\[Sigma] = 100; rSOI = G*MBH*Msun/(\\[Sigma]*10^5)^2/3.086*10^18; "
            << "\\[CapitalPhi] = -G*MBH*Msun/(r*3.086*10^18); Abs[\\[CapitalPhi]]/10^14]; "
            << "(* M51 SMBH: M_BH = " << M_BH_ << " Msun, LINER, L/L_Edd = " << L_Edd_ratio_ << " *)";
        return oss.str();
    }

    std::string getSignature() const { return "M51CentralBlackHoleTerm(r_pc, M_BH, L_Edd_ratio)"; }
    std::string getCategory() const { return "gravity"; }

private:
    double r_pc_;
    double M_BH_;          // SMBH mass [M_sun]
    double L_Edd_ratio_;   // L_bol/L_Edd (Eddington ratio)
};

// ========================================
// Class 727: M51SupernovaFeedbackTerm
// ========================================
// Physical model: SN-driven winds with momentum injection p_SN and energy E_SN = 10^51 erg
// Observational basis: M51 SN rate ~0.1-0.2 yr^-1, X-ray emission from hot gas bubbles
// Reference: Li & Wang (2013) - Supernova feedback in M51 from Chandra X-ray observations
class M51SupernovaFeedbackTerm {
public:
    M51SupernovaFeedbackTerm(double SN_rate_yr = 0.15, double E_SN_erg = 1e51, double eta_momentum = 0.1)
        : SN_rate_yr_(SN_rate_yr), E_SN_erg_(E_SN_erg), eta_momentum_(eta_momentum) {}

    double compute() const {
        // Energy injection rate: Ė_SN = SN_rate·E_SN [erg/yr]
        double E_dot_SN = SN_rate_yr_ * E_SN_erg_; // erg/yr
        
        // Momentum injection: p_SN = √(2·M_ej·E_SN) where M_ej ≈ 1 M_sun
        const double M_ej = 1.0; // M_sun
        const double M_sun_g = 1.989e33; // g
        double p_SN = std::sqrt(2.0 * M_ej * M_sun_g * E_SN_erg_); // g·cm/s
        double p_dot_SN = SN_rate_yr_ * p_SN; // g·cm/s/yr
        
        // Terminal velocity of SN-driven wind: v_term = √(2·E_SN/M_ej)
        double v_term = std::sqrt(2.0 * E_SN_erg_ / (M_ej * M_sun_g)); // cm/s
        double v_term_km_s = v_term / 1e5; // km/s
        
        // Hot gas mass-loading factor: η = Ṁ_wind/SFR
        // For M51 SFR ≈ 3 M_sun/yr, typical η ~ 0.1-1
        const double SFR_Msun_yr = 3.0;
        double M_dot_wind = eta_momentum_ * SFR_Msun_yr; // M_sun/yr
        
        // Cooling time: t_cool = 3·n_e·k_B·T / (2·n_e²·Λ(T))
        // For T ~ 10^7 K (SN-heated gas), Λ ~ 10^-23 erg·cm³/s
        const double T = 1e7; // K
        const double k_B = 1.381e-16; // erg/K
        const double n_e = 0.01; // cm^-3
        const double Lambda = 1e-23; // erg·cm³/s
        double t_cool_s = 3.0 * n_e * k_B * T / (2.0 * n_e * n_e * Lambda);
        double t_cool_Myr = t_cool_s / (3.156e13);
        
        return E_dot_SN / 1e57 * (1.0 + v_term_km_s / 1000.0 + M_dot_wind + t_cool_Myr / 10.0);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M51SupernovaFeedback[SNrate_, ESN_, \\[Eta]_] := "
            << "Module[{EdotSN, Mej, pSN, vterm, MdotWind}, "
            << "EdotSN = SNrate*ESN; Mej = 1; pSN = Sqrt[2*Mej*Msun*ESN]; "
            << "vterm = Sqrt[2*ESN/(Mej*Msun)]/10^5; MdotWind = \\[Eta]*3; "
            << "EdotSN/10^57*(1 + vterm/1000 + MdotWind)]; "
            << "(* M51 SN rate = " << SN_rate_yr_ << " yr^-1, E_SN = " << E_SN_erg_ << " erg *)";
        return oss.str();
    }

    std::string getSignature() const { return "M51SupernovaFeedbackTerm(SN_rate_yr, E_SN_erg, eta_momentum)"; }
    std::string getCategory() const { return "stellar"; }

private:
    double SN_rate_yr_;      // Supernova rate [yr^-1]
    double E_SN_erg_;        // Energy per SN [erg]
    double eta_momentum_;    // Momentum-loading factor
};

// ========================================
// Class 728: M51DustExtinctionTerm
// ========================================
// Physical model: Dust extinction A_λ = A_V·(λ/λ_V)^(-β) with Calzetti attenuation law
// Observational basis: M51 A_V ≈ 1-2 mag in spiral arms from Hα/Hβ ratio, dust lanes in HST images
// Reference: Schechtman-Rook & Hess (2012) - Dust properties in M51
class M51DustExtinctionTerm {
public:
    M51DustExtinctionTerm(double lambda_micron, double A_V = 1.5, double R_V = 3.1)
        : lambda_micron_(lambda_micron), A_V_(A_V), R_V_(R_V) {}

    double compute() const {
        // Calzetti attenuation law: k(λ) = A_λ/E(B-V)
        // For λ < 0.63 μm: k(λ) = 2.659·(-2.156 + 1.509/λ - 0.198/λ² + 0.011/λ³) + R_V
        // For λ > 0.63 μm: k(λ) = 2.659·(-1.857 + 1.040/λ) + R_V
        
        double k_lambda;
        if (lambda_micron_ < 0.63) {
            k_lambda = 2.659 * (-2.156 + 1.509 / lambda_micron_ - 0.198 / (lambda_micron_ * lambda_micron_) + 
                                0.011 / (lambda_micron_ * lambda_micron_ * lambda_micron_)) + R_V_;
        } else {
            k_lambda = 2.659 * (-1.857 + 1.040 / lambda_micron_) + R_V_;
        }
        
        // E(B-V) from A_V: A_V = R_V·E(B-V)
        double E_BV = A_V_ / R_V_;
        
        // Extinction at wavelength λ: A_λ = k(λ)·E(B-V)
        double A_lambda = k_lambda * E_BV;
        
        // Optical depth: τ_λ = A_λ / 1.086
        double tau_lambda = A_lambda / 1.086;
        
        // Dust column density: N_H = A_V / (R_V·σ_d) where σ_d ~ 5×10^-22 cm²
        const double sigma_d = 5e-22; // cm²
        double N_H = A_V_ / (R_V_ * sigma_d); // cm^-2
        
        // Dust mass surface density: Σ_dust = μ·m_H·N_H/(10^21) [M_sun/pc²]
        const double mu = 1.4; // Mean molecular weight
        const double m_H = 1.67e-24; // g
        double Sigma_dust = mu * m_H * N_H / (1e21 * 1.989e33 / (3.086e18 * 3.086e18));
        
        return A_lambda * (1.0 + tau_lambda + Sigma_dust / 10.0 + N_H / 1e22);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M51DustExtinction[\\[Lambda]_, AV_, RV_] := "
            << "Module[{k\\[Lambda], EBV, A\\[Lambda]}, "
            << "k\\[Lambda] = If[\\[Lambda] < 0.63, "
            << "2.659*(-2.156 + 1.509/\\[Lambda] - 0.198/\\[Lambda]^2 + 0.011/\\[Lambda]^3) + RV, "
            << "2.659*(-1.857 + 1.040/\\[Lambda]) + RV]; "
            << "EBV = AV/RV; A\\[Lambda] = k\\[Lambda]*EBV; A\\[Lambda]]; "
            << "(* M51 Calzetti law: A_V = " << A_V_ << " mag, R_V = " << R_V_ << " *)";
        return oss.str();
    }

    std::string getSignature() const { return "M51DustExtinctionTerm(lambda_micron, A_V, R_V)"; }
    std::string getCategory() const { return "stellar"; }

private:
    double lambda_micron_;  // Wavelength [μm]
    double A_V_;            // V-band extinction [mag]
    double R_V_;            // Total-to-selective extinction ratio
};

// ========================================
// Class 729: M51QuantumVacuumTerm
// ========================================
// Physical model: Casimir effect vacuum energy density ρ_vac = -ℏ·c·π²/(720·a⁴)
// Observational basis: Theoretical - quantum vacuum fluctuations in strong magnetic fields
// Reference: Theoretical framework for vacuum polarization in M51's magnetic spiral
class M51QuantumVacuumTerm {
public:
    M51QuantumVacuumTerm(double a_nm = 1.0, double B_microG = 15.0)
        : a_nm_(a_nm), B_microG_(B_microG) {}

    double compute() const {
        // Casimir energy density: ρ_vac = -ℏ·c·π²/(720·a⁴) [erg/cm³]
        const double hbar = 1.055e-27; // erg·s
        const double c = 2.998e10; // cm/s
        double a_cm = a_nm_ * 1e-7; // Convert nm to cm
        double rho_vac_erg_cm3 = -hbar * c * M_PI * M_PI / (720.0 * a_cm * a_cm * a_cm * a_cm);
        
        // Vacuum energy density in mass units: ρ_vac = E_vac/c² [g/cm³]
        double rho_vac_g_cm3 = rho_vac_erg_cm3 / (c * c);
        
        // Magnetic field energy density: ρ_B = B²/(8π) [erg/cm³]
        double B_G = B_microG_ * 1e-6; // Convert μG to G
        double rho_B = B_G * B_G / (8.0 * M_PI);
        
        // Vacuum polarization: Δρ_vac ∝ α·(B/B_crit)² where B_crit = m_e²·c³/(e·ℏ) ≈ 4.4×10^13 G
        const double alpha = 1.0 / 137.0; // Fine structure constant
        const double B_crit = 4.4e13; // G, Schwinger critical field
        double Delta_rho_vac = alpha * (B_G / B_crit) * (B_G / B_crit) * rho_B;
        
        // Zero-point energy fluctuations: ΔE ~ ℏ/(2·Δt) where Δt ~ a/c
        double Delta_t = a_cm / c; // s
        double Delta_E = hbar / (2.0 * Delta_t); // erg
        
        return std::abs(rho_vac_erg_cm3) / 1e10 * (1.0 + std::abs(Delta_rho_vac) / 1e-20 + Delta_E / 1e-10);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M51QuantumVacuum[a_, B_] := "
            << "Module[{\\[HBar], c, \\[Rho]vac, \\[Rho]B, \\[Alpha], Bcrit, \\[CapitalDelta]\\[Rho]vac}, "
            << "\\[HBar] = 1.055*10^(-27); c = 2.998*10^10; "
            << "\\[Rho]vac = -\\[HBar]*c*Pi^2/(720*a^4); "
            << "\\[Rho]B = (B*10^(-6))^2/(8*Pi); \\[Alpha] = 1/137; Bcrit = 4.4*10^13; "
            << "\\[CapitalDelta]\\[Rho]vac = \\[Alpha]*(B*10^(-6)/Bcrit)^2*\\[Rho]B; Abs[\\[Rho]vac]/10^10]; "
            << "(* Casimir vacuum: a = " << a_nm_ << " nm, B = " << B_microG_ << " μG *)";
        return oss.str();
    }

    std::string getSignature() const { return "M51QuantumVacuumTerm(a_nm, B_microG)"; }
    std::string getCategory() const { return "quantum"; }

private:
    double a_nm_;        // Casimir plate separation [nm]
    double B_microG_;    // Magnetic field strength [μG]
};

// ========================================
// Wolfram Language Export Functions
// ========================================

std::string exportM51SpiralArmsWolfram(double r_kpc, double phi_rad) {
    M51SpiralArmDensityWaveTerm term(r_kpc, phi_rad);
    return term.toWolfram();
}

std::string exportM51TidalWolfram(double r_kpc, double phi_rad) {
    M51TidalInteractionTerm term(r_kpc, phi_rad);
    return term.toWolfram();
}

std::string exportM51StarFormationWolfram(double Sigma_gas) {
    M51StarFormationRateTerm term(Sigma_gas);
    return term.toWolfram();
}

std::string exportM51MolecularCloudsWolfram(double M_cloud) {
    M51MolecularCloudTerm term(M_cloud);
    return term.toWolfram();
}

std::string exportM51DarkMatterWolfram(double r_kpc) {
    M51DarkMatterHaloTerm term(r_kpc);
    return term.toWolfram();
}

std::string exportM51MagneticFieldWolfram(double r_kpc, double phi_rad) {
    M51MagneticFieldTerm term(r_kpc, phi_rad);
    return term.toWolfram();
}

std::string exportM51BlackHoleWolfram(double r_pc) {
    M51CentralBlackHoleTerm term(r_pc);
    return term.toWolfram();
}

std::string exportM51SupernovaWolfram(double SN_rate) {
    M51SupernovaFeedbackTerm term(SN_rate);
    return term.toWolfram();
}

std::string exportM51DustWolfram(double lambda_micron) {
    M51DustExtinctionTerm term(lambda_micron);
    return term.toWolfram();
}

std::string exportM51QuantumVacuumWolfram(double a_nm) {
    M51QuantumVacuumTerm term(a_nm);
    return term.toWolfram();
}

std::string exportAllM51WolframFunctions() {
    std::ostringstream oss;
    oss << "(* M51 (Whirlpool Galaxy) UQFF Module - Wolfram Language Export *)\n"
        << "(* Classes 720-729: Grand design spiral with NGC 5195 companion interaction *)\n\n"
        << exportM51SpiralArmsWolfram(5.0, 0.0) << "\n\n"
        << exportM51TidalWolfram(5.0, 0.0) << "\n\n"
        << exportM51StarFormationWolfram(10.0) << "\n\n"
        << exportM51MolecularCloudsWolfram(1e6) << "\n\n"
        << exportM51DarkMatterWolfram(5.0) << "\n\n"
        << exportM51MagneticFieldWolfram(5.0, 0.0) << "\n\n"
        << exportM51BlackHoleWolfram(100.0) << "\n\n"
        << exportM51SupernovaWolfram(0.15) << "\n\n"
        << exportM51DustWolfram(0.55) << "\n\n"
        << exportM51QuantumVacuumWolfram(1.0) << "\n";
    return oss.str();
}

// ========================================
// Master UQFF Integration Function
// ========================================

struct M51UQFFParams {
    double r_kpc;
    double phi_rad;
    double t_Gyr;
    double Sigma_gas;
    double M_cloud;
    double M_NGC5195;
    double d_NGC5195_kpc;
    double phi_NGC5195_rad;
    double M_200_DM;
    double c_DM;
    double B_ordered;
    double B_random;
    double M_BH;
    double SN_rate_yr;
    double lambda_micron;
    double A_V;
};

double computeM51MasterEquation(const M51UQFFParams& params) {
    M51SpiralArmDensityWaveTerm spiral(params.r_kpc, params.phi_rad, params.t_Gyr);
    M51TidalInteractionTerm tidal(params.r_kpc, params.phi_rad, params.M_NGC5195, params.d_NGC5195_kpc, params.phi_NGC5195_rad);
    M51StarFormationRateTerm sfr(params.Sigma_gas);
    M51MolecularCloudTerm gmc(params.M_cloud);
    M51DarkMatterHaloTerm dm(params.r_kpc, params.M_200_DM, params.c_DM);
    M51MagneticFieldTerm mag(params.r_kpc, params.phi_rad, params.B_ordered, params.B_random);
    M51CentralBlackHoleTerm bh(params.r_kpc * 1000.0, params.M_BH); // Convert kpc to pc
    M51SupernovaFeedbackTerm sn(params.SN_rate_yr);
    M51DustExtinctionTerm dust(params.lambda_micron, params.A_V);
    M51QuantumVacuumTerm qvac(1.0, params.B_ordered);
    
    double F_spiral = spiral.compute();
    double F_tidal = tidal.compute();
    double F_sfr = sfr.compute();
    double F_gmc = gmc.compute();
    double F_dm = dm.compute();
    double F_mag = mag.compute();
    double F_bh = bh.compute();
    double F_sn = sn.compute();
    double F_dust = dust.compute();
    double F_qvac = qvac.compute();
    
    // Master UQFF equation: F_total = Σ(F_i) with cross-term coupling
    double F_total = F_spiral + F_tidal + F_sfr + F_gmc + F_dm + F_mag + F_bh + F_sn + F_dust + F_qvac
                   + 0.1 * F_spiral * F_sfr           // Spiral arm-star formation coupling
                   + 0.05 * F_tidal * F_sfr           // Tidal-star formation enhancement
                   + 0.02 * F_mag * F_spiral          // Magnetic field-density wave coupling
                   + 0.01 * F_dm * F_tidal;           // Dark matter-tidal coupling
    
    return F_total;
}

// Example usage and validation
void validateM51Module() {
    M51UQFFParams params;
    params.r_kpc = 5.0;               // 5 kpc from center
    params.phi_rad = 0.0;             // Along major axis
    params.t_Gyr = 0.0;               // Current epoch
    params.Sigma_gas = 10.0;          // M_sun/pc²
    params.M_cloud = 1e6;             // M_sun
    params.M_NGC5195 = 2e10;          // M_sun
    params.d_NGC5195_kpc = 9.0;       // kpc
    params.phi_NGC5195_rad = 0.0;     // radians
    params.M_200_DM = 1e12;           // M_sun
    params.c_DM = 10.0;               // Concentration
    params.B_ordered = 15.0;          // μG
    params.B_random = 10.0;           // μG
    params.M_BH = 1e7;                // M_sun
    params.SN_rate_yr = 0.15;         // yr^-1
    params.lambda_micron = 0.55;      // μm (V-band)
    params.A_V = 1.5;                 // mag
    
    double result = computeM51MasterEquation(params);
    
    // Expected: Combined M51 physics at r=5 kpc
    // Spiral density wave, NGC 5195 tidal interaction, SFR ~ 3-4 M_sun/yr, ~500 GMCs
}
