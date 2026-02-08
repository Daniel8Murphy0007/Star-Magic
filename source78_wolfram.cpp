// source78_wolfram.cpp - M81/M82 Galaxy Pair UQFF Module
// Classes 730-739: M81 grand design spiral + M82 starburst galaxy interaction
// Physical basis: Prototypical galaxy pair with tidal interaction, M82 starburst triggered by M81 encounter

#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ========================================
// Class 730: M81SpiralStructureTerm
// ========================================
// Physical model: Grand design spiral with m=2 mode, similar to M51 but less prominent
// Observational basis: M81 rotation curve v_rot ≈ 200 km/s, spiral pitch angle i ≈ 15°
// Reference: Kaufman et al. (1989) - HI observations of M81
class M81SpiralStructureTerm {
public:
    M81SpiralStructureTerm(double r_kpc, double phi_rad, double t_Gyr = 0.0,
                           double A = 0.2, double m = 2.0, double pitch_deg = 15.0,
                           double Omega_p = 30.0)
        : r_kpc_(r_kpc), phi_rad_(phi_rad), t_Gyr_(t_Gyr), A_(A), m_(m),
          pitch_deg_(pitch_deg), Omega_p_(Omega_p) {}

    double compute() const {
        // Logarithmic spiral: r = r₀·exp(φ·tan(i))
        double pitch_rad = pitch_deg_ * M_PI / 180.0;
        double k = m_ / (r_kpc_ * std::tan(pitch_rad)); // kpc^-1
        
        // Phase: ψ = m·φ - k·r - Ω_p·t
        double psi = m_ * phi_rad_ - k * r_kpc_ - Omega_p_ * t_Gyr_;
        
        // Surface density: Σ(r,φ,t) = Σ₀(r)·[1 + A·cos(ψ)]
        double Sigma_0 = 150.0 * std::exp(-r_kpc_ / 4.0); // M_sun/pc²
        double Sigma = Sigma_0 * (1.0 + A_ * std::cos(psi));
        
        // Corotation radius: r_CR where Ω(r_CR) = Ω_p
        double v_circ = 200.0; // km/s, approximately flat
        double r_CR = v_circ / Omega_p_; // kpc
        
        // Lindblad resonances: Ω ± κ/m = Ω_p
        // For flat rotation curve, κ ≈ √2·Ω
        double Omega_r = v_circ / r_kpc_; // km/s/kpc
        double kappa = std::sqrt(2.0) * Omega_r;
        double r_ILR = v_circ / (Omega_p_ + kappa / m_); // Inner Lindblad resonance
        double r_OLR = v_circ / (Omega_p_ - kappa / m_); // Outer Lindblad resonance
        
        return Sigma * (1.0 + A_ + std::abs(r_kpc_ - r_CR) / r_CR + std::abs(r_kpc_ - r_ILR) / 10.0);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M81SpiralDensity[r_, \\[Phi]_, t_, A_, m_, pitch_, \\[CapitalOmega]p_] := "
            << "Module[{k, \\[Psi], \\[CapitalSigma]0}, "
            << "k = m/(r*Tan[pitch*Degree]); \\[Psi] = m*\\[Phi] - k*r - \\[CapitalOmega]p*t; "
            << "\\[CapitalSigma]0 = 150*Exp[-r/4]; \\[CapitalSigma]0*(1 + A*Cos[\\[Psi]])]; "
            << "(* M81 m=2 spiral, pitch = " << pitch_deg_ << "°, v_rot ~ 200 km/s *)";
        return oss.str();
    }

    std::string getSignature() const { return "M81SpiralStructureTerm(r_kpc, phi_rad, t_Gyr, A, m, pitch_deg, Omega_p)"; }
    std::string getCategory() const { return "dynamics"; }

private:
    double r_kpc_;
    double phi_rad_;
    double t_Gyr_;
    double A_;           // Amplitude
    double m_;           // Mode number
    double pitch_deg_;   // Pitch angle [degrees]
    double Omega_p_;     // Pattern speed [km/s/kpc]
};

// ========================================
// Class 731: M82StarburstTerm
// ========================================
// Physical model: Nuclear starburst with SFR ≈ 10 M_sun/yr concentrated in central ~500 pc
// Observational basis: M82 infrared luminosity L_IR ≈ 4×10^10 L_sun, superwind outflow
// Reference: Förster Schreiber et al. (2003) - M82 starburst kinematics from SINFONI
class M82StarburstTerm {
public:
    M82StarburstTerm(double r_pc, double SFR_central = 10.0, double r_starburst_pc = 500.0)
        : r_pc_(r_pc), SFR_central_(SFR_central), r_starburst_pc_(r_starburst_pc) {}

    double compute() const {
        // Starburst surface density profile: Σ_SFR(r) = Σ_SFR,0·exp(-r²/r_sb²)
        double Sigma_SFR_0 = SFR_central_ / (M_PI * r_starburst_pc_ * r_starburst_pc_); // M_sun/yr/pc²
        double Sigma_SFR = Sigma_SFR_0 * std::exp(-r_pc_ * r_pc_ / (r_starburst_pc_ * r_starburst_pc_));
        
        // Star formation efficiency: ε = SFR/(M_gas/t_ff)
        const double M_gas_central = 5e8; // M_sun, molecular gas in central region
        const double t_ff_Myr = 10.0; // Myr, free-fall time for dense clouds
        double epsilon = SFR_central_ / (M_gas_central / (t_ff_Myr * 1e6));
        
        // Infrared luminosity: L_IR = η·SFR·10^10 L_sun where η ~ 4 for M82
        double L_IR = 4.0 * SFR_central_ * 1e10; // L_sun
        
        // Dust temperature: T_dust ~ 40-60 K in starburst core
        const double T_dust = 50.0; // K
        
        // Star formation rate density: ρ_SFR = Σ_SFR/h where h ~ 100 pc scale height
        const double h_gas = 100.0; // pc
        double rho_SFR = Sigma_SFR / h_gas; // M_sun/yr/pc³
        
        return Sigma_SFR * (1.0 + epsilon + L_IR / 1e11 + T_dust / 100.0 + rho_SFR * 1e3);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M82StarburstSFR[r_, SFRcentral_, rsb_] := "
            << "Module[{\\[CapitalSigma]SFR0, \\[CapitalSigma]SFR, LIR}, "
            << "\\[CapitalSigma]SFR0 = SFRcentral/(Pi*rsb^2); "
            << "\\[CapitalSigma]SFR = \\[CapitalSigma]SFR0*Exp[-r^2/rsb^2]; "
            << "LIR = 4*SFRcentral*10^10; {\\[CapitalSigma]SFR, LIR}]; "
            << "(* M82 nuclear starburst: SFR = " << SFR_central_ << " Msun/yr, r_sb = " << r_starburst_pc_ << " pc *)";
        return oss.str();
    }

    std::string getSignature() const { return "M82StarburstTerm(r_pc, SFR_central, r_starburst_pc)"; }
    std::string getCategory() const { return "stellar"; }

private:
    double r_pc_;
    double SFR_central_;      // Central SFR [M_sun/yr]
    double r_starburst_pc_;   // Starburst radius [pc]
};

// ========================================
// Class 732: M82SuperwindTerm
// ========================================
// Physical model: Bipolar superwind with v_wind ≈ 500-1000 km/s, Ṁ_wind ≈ 3·SFR
// Observational basis: M82 Hα filaments extend to ~10 kpc, X-ray hot gas T ≈ 10^7 K
// Reference: Strickland & Heckman (2009) - M82 superwind energetics
class M82SuperwindTerm {
public:
    M82SuperwindTerm(double z_kpc, double v_wind = 750.0, double M_dot_wind = 30.0, double T_wind = 1e7)
        : z_kpc_(z_kpc), v_wind_(v_wind), M_dot_wind_(M_dot_wind), T_wind_(T_wind) {}

    double compute() const {
        // Wind density profile: ρ_wind(z) = ρ_0·(r_0/r)² for momentum-driven wind
        const double r_0 = 0.5; // kpc, launch radius
        double r = std::sqrt(r_0 * r_0 + z_kpc_ * z_kpc_); // kpc, distance from center
        
        // Mass-loss rate to density: ρ_0 = Ṁ_wind/(4π·r_0²·v_wind)
        const double M_sun_g = 1.989e33; // g
        const double kpc_cm = 3.086e21; // cm
        const double yr_s = 3.156e7; // s
        double M_dot_wind_g_s = M_dot_wind_ * M_sun_g / yr_s; // g/s
        double rho_0 = M_dot_wind_g_s / (4.0 * M_PI * r_0 * r_0 * kpc_cm * kpc_cm * v_wind_ * 1e5); // g/cm³
        
        double rho_wind = rho_0 * (r_0 / r) * (r_0 / r); // g/cm³
        
        // Kinetic energy flux: Ė_kin = ½·Ṁ_wind·v_wind² [erg/s]
        double E_dot_kin = 0.5 * M_dot_wind_g_s * v_wind_ * 1e5 * v_wind_ * 1e5; // erg/s
        
        // Thermal energy flux: Ė_th = Ṁ_wind·c_p·T where c_p = 5k_B/(2μ·m_H)
        const double k_B = 1.381e-16; // erg/K
        const double mu = 0.6; // Mean molecular weight
        const double m_H = 1.67e-24; // g
        double c_p = 5.0 * k_B / (2.0 * mu * m_H); // erg/g/K
        double E_dot_th = M_dot_wind_g_s * c_p * T_wind_; // erg/s
        
        // Total wind power: Ė_wind = Ė_kin + Ė_th
        double E_dot_wind = E_dot_kin + E_dot_th; // erg/s
        
        // Sound speed in hot gas: c_s = √(γ·k_B·T/(μ·m_H)) where γ = 5/3
        double c_s = std::sqrt(5.0 / 3.0 * k_B * T_wind_ / (mu * m_H)); // cm/s
        double c_s_km_s = c_s / 1e5; // km/s
        
        return rho_wind * 1e24 * (1.0 + E_dot_wind / 1e42 + v_wind_ / 1000.0 + c_s_km_s / 100.0);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M82SuperwindDensity[z_, vwind_, Mdotwind_, Twind_] := "
            << "Module[{r0, r, \\[Rho]0, \\[Rho]wind, EdotKin, EdotTh}, "
            << "r0 = 0.5; r = Sqrt[r0^2 + z^2]; "
            << "\\[Rho]0 = Mdotwind*Msun/(4*Pi*r0^2*kpc^2*vwind*10^5*yr); "
            << "\\[Rho]wind = \\[Rho]0*(r0/r)^2; "
            << "EdotKin = 0.5*Mdotwind*Msun*vwind^2*10^10/yr; "
            << "{\\[Rho]wind, EdotKin}]; "
            << "(* M82 superwind: v = " << v_wind_ << " km/s, Ṁ = " << M_dot_wind_ << " Msun/yr *)";
        return oss.str();
    }

    std::string getSignature() const { return "M82SuperwindTerm(z_kpc, v_wind, M_dot_wind, T_wind)"; }
    std::string getCategory() const { return "dynamics"; }

private:
    double z_kpc_;
    double v_wind_;       // Wind velocity [km/s]
    double M_dot_wind_;   // Mass-loss rate [M_sun/yr]
    double T_wind_;       // Wind temperature [K]
};

// ========================================
// Class 733: M81M82TidalInteractionTerm
// ========================================
// Physical model: M81-M82 separation ~150 kpc, tidal bridge/tail, closest approach ~300 Myr ago
// Observational basis: HI tidal features extend to ~200 kpc, M82 starburst triggered by M81
// Reference: Yun et al. (1994) - HI shells and tidal features in M81/M82 group
class M81M82TidalInteractionTerm {
public:
    M81M82TidalInteractionTerm(double r_kpc, double phi_rad,
                               double M_M81 = 7e10, double M_M82 = 3e10,
                               double d_kpc = 150.0, double phi_M82_rad = 0.0)
        : r_kpc_(r_kpc), phi_rad_(phi_rad), M_M81_(M_M81), M_M82_(M_M82),
          d_kpc_(d_kpc), phi_M82_rad_(phi_M82_rad) {}

    double compute() const {
        // M82 position relative to M81 center
        double x_M82 = d_kpc_ * std::cos(phi_M82_rad_);
        double y_M82 = d_kpc_ * std::sin(phi_M82_rad_);
        
        // Test particle position
        double x = r_kpc_ * std::cos(phi_rad_);
        double y = r_kpc_ * std::sin(phi_rad_);
        
        // Separation from M82
        double dx = x - x_M82;
        double dy = y - y_M82;
        double r_sep = std::sqrt(dx * dx + dy * dy);
        
        // Tidal force from M82 on test particle in M81
        const double G = 4.3e-6; // kpc·(km/s)²/M_sun
        double F_tid_mag = G * M_M82_ / (r_sep * r_sep);
        
        // Tidal radius: r_tid = d·(M_M81/(3·M_M82))^(1/3)
        double r_tid = d_kpc_ * std::pow(M_M81_ / (3.0 * M_M82_), 1.0/3.0);
        
        // Relative velocity at closest approach: Δv ~ √(G·(M_M81 + M_M82)/d)
        double Delta_v = std::sqrt(G * (M_M81_ + M_M82_) / d_kpc_);
        
        // Dynamical time: t_dyn = d/Δv ~ 300 Myr (time since closest approach)
        double t_dyn_yr = d_kpc_ * 3.086e16 / (Delta_v * 1e5 * 3.156e7);
        double t_dyn_Myr = t_dyn_yr / 1e6;
        
        // Tidal heating parameter: Q = (r/r_tid)³·(M_M82/M_M81)
        double Q_tidal = std::pow(r_kpc_ / r_tid, 3.0) * (M_M82_ / M_M81_);
        
        return F_tid_mag * (1.0 + Q_tidal + Delta_v / 100.0 + t_dyn_Myr / 300.0);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M81M82TidalForce[r_, \\[Phi]_, MM81_, MM82_, d_, \\[Phi]M82_] := "
            << "Module[{xM82, yM82, x, y, dx, dy, rsep, Ftid, rtid}, "
            << "xM82 = d*Cos[\\[Phi]M82]; yM82 = d*Sin[\\[Phi]M82]; "
            << "x = r*Cos[\\[Phi]]; y = r*Sin[\\[Phi]]; "
            << "dx = x - xM82; dy = y - yM82; rsep = Sqrt[dx^2 + dy^2]; "
            << "Ftid = G*MM82/rsep^2; rtid = d*(MM81/(3*MM82))^(1/3); "
            << "{Ftid, rtid}]; "
            << "(* M81-M82 d = " << d_kpc_ << " kpc, closest approach ~300 Myr ago *)";
        return oss.str();
    }

    std::string getSignature() const { return "M81M82TidalInteractionTerm(r_kpc, phi_rad, M_M81, M_M82, d_kpc, phi_M82_rad)"; }
    std::string getCategory() const { return "dynamics"; }

private:
    double r_kpc_;
    double phi_rad_;
    double M_M81_;         // M81 mass [M_sun]
    double M_M82_;         // M82 mass [M_sun]
    double d_kpc_;         // M81-M82 separation [kpc]
    double phi_M82_rad_;   // M82 position angle [radians]
};

// ========================================
// Class 734: M81DarkMatterHaloTerm
// ========================================
// Physical model: NFW halo with M_200 ~ 10^12 M_sun, concentration c ~ 12
// Observational basis: M81 rotation curve requires dark matter beyond ~15 kpc
// Reference: de Blok et al. (2008) - HI kinematics of M81
class M81DarkMatterHaloTerm {
public:
    M81DarkMatterHaloTerm(double r_kpc, double M_200 = 1e12, double c = 12.0)
        : r_kpc_(r_kpc), M_200_(M_200), c_(c) {}

    double compute() const {
        // NFW parameters
        const double H_0 = 70.0; // km/s/Mpc
        const double rho_crit = 3.0 * H_0 * H_0 / (8.0 * M_PI * 4.3e-6 * 1e6); // M_sun/kpc³
        double r_200 = std::pow(3.0 * M_200_ / (4.0 * M_PI * 200.0 * rho_crit), 1.0/3.0); // kpc
        double r_s = r_200 / c_;
        
        double f_c = std::log(1.0 + c_) - c_ / (1.0 + c_);
        double rho_s = M_200_ / (4.0 * M_PI * r_s * r_s * r_s * f_c);
        
        // NFW density
        double x = r_kpc_ / r_s;
        double rho_DM = rho_s / (x * (1.0 + x) * (1.0 + x));
        
        // Enclosed mass
        double M_DM = 4.0 * M_PI * rho_s * r_s * r_s * r_s * (std::log(1.0 + x) - x / (1.0 + x));
        
        // Circular velocity
        const double G = 4.3e-6; // kpc·(km/s)²/M_sun
        double v_DM = std::sqrt(G * M_DM / r_kpc_);
        
        // Jeans equation: σ_r² = (1/ρ)·∫[ρ·G·M(<r)/r²]·dr
        // For NFW, σ_r ~ v_circ at most radii
        double sigma_r = v_DM / std::sqrt(3.0); // 1D velocity dispersion
        
        return rho_DM * (1.0 + v_DM / 200.0 + M_DM / 1e11 + sigma_r / 100.0);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M81NFWDensity[r_, M200_, c_] := "
            << "Module[{H0, \\[Rho]crit, r200, rs, fc, \\[Rho]s, x, \\[Rho]DM}, "
            << "H0 = 70; \\[Rho]crit = 3*H0^2/(8*Pi*G*10^6); "
            << "r200 = (3*M200/(4*Pi*200*\\[Rho]crit))^(1/3); rs = r200/c; "
            << "fc = Log[1 + c] - c/(1 + c); \\[Rho]s = M200/(4*Pi*rs^3*fc); "
            << "x = r/rs; \\[Rho]DM = \\[Rho]s/(x*(1 + x)^2); \\[Rho]DM]; "
            << "(* M81 NFW halo: M_200 = " << M_200_ << " Msun, c = " << c_ << " *)";
        return oss.str();
    }

    std::string getSignature() const { return "M81DarkMatterHaloTerm(r_kpc, M_200, c)"; }
    std::string getCategory() const { return "dark_matter"; }

private:
    double r_kpc_;
    double M_200_;
    double c_;
};

// ========================================
// Class 735: M82MolecularOutflowTerm
// ========================================
// Physical model: Molecular gas outflow Ṁ_mol ~ 10 M_sun/yr with v_out ~ 200 km/s
// Observational basis: M82 CO observations show outflowing molecular gas to z ~ 1 kpc
// Reference: Walter et al. (2002) - Molecular outflow in M82
class M82MolecularOutflowTerm {
public:
    M82MolecularOutflowTerm(double z_kpc, double M_dot_mol = 10.0, double v_out = 200.0)
        : z_kpc_(z_kpc), M_dot_mol_(M_dot_mol), v_out_(v_out) {}

    double compute() const {
        // Molecular gas surface density in outflow: Σ_mol(z) = Σ_0·exp(-z/z_0)
        const double z_0 = 1.0; // kpc, scale height
        const double Sigma_0 = 50.0; // M_sun/pc², central density
        double Sigma_mol = Sigma_0 * std::exp(-std::abs(z_kpc_) / z_0);
        
        // Mass-loss rate: Ṁ_mol = Σ_mol·v_out·A where A ~ π·r²
        const double r_out = 0.5; // kpc, outflow radius
        double A_kpc2 = M_PI * r_out * r_out; // kpc²
        double M_dot_expected = Sigma_mol * 1e6 * v_out * A_kpc2 * 3.086e16 / (3.156e7 * 1.989e33); // M_sun/yr
        
        // Depletion time: t_dep = M_mol/Ṁ_mol
        const double M_mol_total = 5e8; // M_sun, total molecular gas
        double t_dep_Myr = M_mol_total / M_dot_mol_ / 1e6;
        
        // Kinetic energy flux: Ė_kin = ½·Ṁ_mol·v_out²
        const double M_sun_g = 1.989e33; // g
        const double yr_s = 3.156e7; // s
        double E_dot_kin = 0.5 * M_dot_mol_ * M_sun_g / yr_s * v_out_ * 1e5 * v_out_ * 1e5; // erg/s
        
        // Momentum flux: ṗ = Ṁ_mol·v_out
        double p_dot = M_dot_mol_ * M_sun_g / yr_s * v_out_ * 1e5; // g·cm/s²
        
        return Sigma_mol * (1.0 + M_dot_expected / 10.0 + E_dot_kin / 1e41 + p_dot / 1e35);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M82MolecularOutflow[z_, Mdotmol_, vout_] := "
            << "Module[{z0, \\[CapitalSigma]0, \\[CapitalSigma]mol, EdotKin}, "
            << "z0 = 1; \\[CapitalSigma]0 = 50; \\[CapitalSigma]mol = \\[CapitalSigma]0*Exp[-Abs[z]/z0]; "
            << "EdotKin = 0.5*Mdotmol*Msun*vout^2*10^10/yr; "
            << "{\\[CapitalSigma]mol, EdotKin}]; "
            << "(* M82 molecular outflow: Ṁ = " << M_dot_mol_ << " Msun/yr, v = " << v_out_ << " km/s *)";
        return oss.str();
    }

    std::string getSignature() const { return "M82MolecularOutflowTerm(z_kpc, M_dot_mol, v_out)"; }
    std::string getCategory() const { return "stellar"; }

private:
    double z_kpc_;
    double M_dot_mol_;  // Molecular outflow rate [M_sun/yr]
    double v_out_;      // Outflow velocity [km/s]
};

// ========================================
// Class 736: M81AGNActivityTerm
// ========================================
// Physical model: LINER with M_BH ~ 10^7 M_sun, low-luminosity AGN L_bol/L_Edd ~ 10^-5
// Observational basis: M81 X-ray source, radio core at 1.3 cm
// Reference: Ho et al. (1996) - M81 nuclear emission
class M81AGNActivityTerm {
public:
    M81AGNActivityTerm(double r_pc, double M_BH = 7e7, double L_Edd_ratio = 1e-5)
        : r_pc_(r_pc), M_BH_(M_BH), L_Edd_ratio_(L_Edd_ratio) {}

    double compute() const {
        // Schwarzschild radius
        const double G = 6.674e-8; // cm³/g/s²
        const double c = 2.998e10; // cm/s
        const double M_sun_g = 1.989e33; // g
        double R_S = 2.0 * G * M_BH_ * M_sun_g / (c * c); // cm
        
        // Eddington luminosity
        double L_Edd = 1.26e38 * M_BH_; // erg/s
        double L_bol = L_Edd_ratio_ * L_Edd;
        
        // Bondi accretion radius: r_Bondi = G·M_BH/c_s²
        const double c_s = 1e6; // cm/s, sound speed in ISM
        double r_Bondi = G * M_BH_ * M_sun_g / (c_s * c_s); // cm
        double r_Bondi_pc = r_Bondi / 3.086e18;
        
        // LINER spectrum: [OI]/Hα > 1/3, [NII]/Hα > 0.6
        const double ratio_OI_Ha = 0.5;
        const double ratio_NII_Ha = 0.8;
        
        // Radio luminosity: L_radio ~ 10^38 erg/s at 1.4 GHz
        const double L_radio = 1e38; // erg/s
        
        // Jet power (if present): P_jet ~ 10^(-3)·L_Edd for LINER
        double P_jet = 1e-3 * L_Edd; // erg/s
        
        return L_bol / 1e38 * (1.0 + r_Bondi_pc / r_pc_ + ratio_OI_Ha + P_jet / 1e40);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M81AGNActivity[r_, MBH_, LEddRatio_] := "
            << "Module[{RS, LEdd, Lbol, rBondi}, "
            << "RS = 2*G*MBH*Msun/c^2; LEdd = 1.26*10^38*MBH; "
            << "Lbol = LEddRatio*LEdd; rBondi = G*MBH*Msun/(10^6)^2/3.086*10^18; "
            << "{Lbol, rBondi}]; "
            << "(* M81 LINER: M_BH = " << M_BH_ << " Msun, L/L_Edd = " << L_Edd_ratio_ << " *)";
        return oss.str();
    }

    std::string getSignature() const { return "M81AGNActivityTerm(r_pc, M_BH, L_Edd_ratio)"; }
    std::string getCategory() const { return "gravity"; }

private:
    double r_pc_;
    double M_BH_;
    double L_Edd_ratio_;
};

// ========================================
// Class 737: M82PAHEmissionTerm
// ========================================
// Physical model: Polycyclic Aromatic Hydrocarbon (PAH) emission at 3.3, 6.2, 7.7, 8.6, 11.3 μm
// Observational basis: M82 strong PAH features from UV-excited regions
// Reference: Engelbracht et al. (2006) - Spitzer PAH observations of M82
class M82PAHEmissionTerm {
public:
    M82PAHEmissionTerm(double lambda_micron, double I_UV = 100.0)
        : lambda_micron_(lambda_micron), I_UV_(I_UV) {}

    double compute() const {
        // PAH band wavelengths and relative strengths
        std::map<double, double> PAH_bands;
        PAH_bands[3.3] = 0.1;   // C-H stretch
        PAH_bands[6.2] = 0.3;   // C-C stretch
        PAH_bands[7.7] = 1.0;   // C-C stretch complex (strongest)
        PAH_bands[8.6] = 0.15;  // C-H in-plane bend
        PAH_bands[11.3] = 0.3;  // C-H out-of-plane bend
        
        // Find closest PAH band
        double I_PAH = 0.0;
        double min_delta = 1e10;
        for (const auto& band : PAH_bands) {
            double delta = std::abs(lambda_micron_ - band.first);
            if (delta < min_delta) {
                min_delta = delta;
                I_PAH = band.second;
            }
        }
        
        // PAH intensity scales with UV radiation field: I_PAH ∝ I_UV
        // Normalized to Galactic I_UV = 1 (Habing field)
        double I_PAH_scaled = I_PAH * I_UV_;
        
        // PAH ionization fraction: f_ion = I_UV/(I_UV + 10³) (Draine & Li 2007)
        double f_ion = I_UV_ / (I_UV_ + 1000.0);
        
        // Neutral vs. ionized PAH emission ratio
        double neutral_to_ion_ratio = (1.0 - f_ion) / f_ion;
        
        // PAH mass fraction: ~4.6% of dust mass in starburst galaxies
        const double f_PAH = 0.046;
        
        return I_PAH_scaled * (1.0 + f_ion + neutral_to_ion_ratio + f_PAH);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M82PAHEmission[\\[Lambda]_, IUV_] := "
            << "Module[{PAHbands, IPAH, fion}, "
            << "PAHbands = {3.3 -> 0.1, 6.2 -> 0.3, 7.7 -> 1.0, 8.6 -> 0.15, 11.3 -> 0.3}; "
            << "IPAH = Lookup[PAHbands, Nearest[Keys[PAHbands], \\[Lambda]][[1]]]; "
            << "fion = IUV/(IUV + 1000); IPAH*IUV*(1 + fion)]; "
            << "(* M82 PAH bands: 3.3, 6.2, 7.7 (strongest), 8.6, 11.3 μm *)";
        return oss.str();
    }

    std::string getSignature() const { return "M82PAHEmissionTerm(lambda_micron, I_UV)"; }
    std::string getCategory() const { return "stellar"; }

private:
    double lambda_micron_;
    double I_UV_;  // UV radiation field intensity (Habing units)
};

// ========================================
// Class 738: M81M82HIDiskTerm
// ========================================
// Physical model: Extended HI disk with M_HI ~ 10^9 M_sun, tidal streamers
// Observational basis: M81 group HI mass ~10^10 M_sun, filaments connect M81/M82
// Reference: Chynoweth et al. (2008) - HI in M81 group
class M81M82HIDiskTerm {
public:
    M81M82HIDiskTerm(double r_kpc, double M_HI = 1e9, double r_HI = 30.0)
        : r_kpc_(r_kpc), M_HI_(M_HI), r_HI_(r_HI) {}

    double compute() const {
        // HI surface density: Σ_HI(r) = Σ_HI,0·exp(-r/r_HI)
        double Sigma_HI_0 = M_HI_ / (2.0 * M_PI * r_HI_ * r_HI_ * 1e6); // M_sun/pc²
        double Sigma_HI = Sigma_HI_0 * std::exp(-r_kpc_ / r_HI_);
        
        // HI column density: N_HI = Σ_HI/(μ·m_H) where μ = 1.4
        const double mu = 1.4;
        const double m_H = 1.67e-24; // g
        const double M_sun_g = 1.989e33; // g
        const double pc_cm = 3.086e18; // cm
        double N_HI = Sigma_HI * M_sun_g / (pc_cm * pc_cm * mu * m_H); // cm^-2
        
        // HI 21-cm optical depth: τ_21 = 5.5×10^(-19)·T_s^(-1)·N_HI [dimensionless]
        const double T_s = 100.0; // K, spin temperature
        double tau_21 = 5.5e-19 * N_HI / T_s;
        
        // HI mass within radius r: M_HI(<r) = M_HI·[1 - exp(-r/r_HI)·(1 + r/r_HI)]
        double M_HI_enclosed = M_HI_ * (1.0 - std::exp(-r_kpc_ / r_HI_) * (1.0 + r_kpc_ / r_HI_));
        
        // HI velocity dispersion: σ_HI ~ 10 km/s (turbulent + thermal)
        const double sigma_HI = 10.0; // km/s
        
        return Sigma_HI * (1.0 + tau_21 + M_HI_enclosed / 1e9 + sigma_HI / 10.0);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M81M82HIDisk[r_, MHI_, rHI_] := "
            << "Module[{\\[CapitalSigma]HI0, \\[CapitalSigma]HI, NHI, \\[Tau]21}, "
            << "\\[CapitalSigma]HI0 = MHI/(2*Pi*rHI^2*10^6); "
            << "\\[CapitalSigma]HI = \\[CapitalSigma]HI0*Exp[-r/rHI]; "
            << "NHI = \\[CapitalSigma]HI*Msun/(pc^2*1.4*mH); "
            << "\\[Tau]21 = 5.5*10^(-19)*NHI/100; {\\[CapitalSigma]HI, NHI}]; "
            << "(* M81 group M_HI ~ " << M_HI_ << " Msun, r_HI = " << r_HI_ << " kpc *)";
        return oss.str();
    }

    std::string getSignature() const { return "M81M82HIDiskTerm(r_kpc, M_HI, r_HI)"; }
    std::string getCategory() const { return "stellar"; }

private:
    double r_kpc_;
    double M_HI_;
    double r_HI_;
};

// ========================================
// Class 739: M82CosmicRayTerm
// ========================================
// Physical model: Cosmic ray pressure P_CR = (γ_CR - 1)·u_CR where u_CR is energy density
// Observational basis: M82 gamma-ray emission suggests enhanced CR production
// Reference: VERITAS Collaboration (2009) - TeV gamma rays from M82
class M82CosmicRayTerm {
public:
    M82CosmicRayTerm(double SN_rate_yr = 0.1, double E_CR_per_SN = 1e49)
        : SN_rate_yr_(SN_rate_yr), E_CR_per_SN_(E_CR_per_SN) {}

    double compute() const {
        // CR injection rate: Ė_CR = SN_rate·E_CR_per_SN [erg/yr]
        double E_dot_CR = SN_rate_yr_ * E_CR_per_SN_; // erg/yr
        
        // CR energy density: u_CR = Ė_CR·t_loss/V
        // For M82, t_loss ~ 10 Myr (advection + diffusion)
        const double t_loss_yr = 1e7; // yr
        const double V_kpc3 = 1.0; // kpc³, starburst volume
        const double kpc_cm = 3.086e21; // cm
        double u_CR = E_dot_CR * t_loss_yr / (V_kpc3 * kpc_cm * kpc_cm * kpc_cm); // erg/cm³
        
        // CR pressure: P_CR = (γ_CR - 1)·u_CR where γ_CR = 4/3 for relativistic particles
        const double gamma_CR = 4.0 / 3.0;
        double P_CR = (gamma_CR - 1.0) * u_CR; // dyne/cm²
        
        // Equipartition magnetic field: B_eq = √(8π·P_CR) assuming P_B ~ P_CR
        double B_eq = std::sqrt(8.0 * M_PI * P_CR); // G
        double B_eq_microG = B_eq * 1e6; // μG
        
        // CR diffusion coefficient: D_CR ~ 10^28 cm²/s (Bohm limit)
        const double D_CR = 1e28; // cm²/s
        
        // CR streaming velocity: v_stream = D_CR/L where L ~ 500 pc
        const double L_pc = 500.0;
        double v_stream = D_CR / (L_pc * 3.086e18); // cm/s
        double v_stream_km_s = v_stream / 1e5; // km/s
        
        return u_CR / 1e-11 * (1.0 + P_CR / 1e-12 + B_eq_microG / 100.0 + v_stream_km_s / 100.0);
    }

    std::string toWolfram() const {
        std::ostringstream oss;
        oss << "M82CosmicRayPressure[SNrate_, ECRperSN_] := "
            << "Module[{EdotCR, tloss, V, uCR, \\[Gamma]CR, PCR, Beq}, "
            << "EdotCR = SNrate*ECRperSN; tloss = 10^7; V = 1*kpc^3; "
            << "uCR = EdotCR*tloss/V; \\[Gamma]CR = 4/3; PCR = (\\[Gamma]CR - 1)*uCR; "
            << "Beq = Sqrt[8*Pi*PCR]; {uCR, PCR, Beq}]; "
            << "(* M82 CR injection: SN rate = " << SN_rate_yr_ << " yr^-1, E_CR = " << E_CR_per_SN_ << " erg/SN *)";
        return oss.str();
    }

    std::string getSignature() const { return "M82CosmicRayTerm(SN_rate_yr, E_CR_per_SN)"; }
    std::string getCategory() const { return "stellar"; }

private:
    double SN_rate_yr_;
    double E_CR_per_SN_;  // CR energy per supernova [erg]
};

// ========================================
// Wolfram Language Export Functions
// ========================================

std::string exportM81SpiralWolfram(double r_kpc, double phi_rad) {
    M81SpiralStructureTerm term(r_kpc, phi_rad);
    return term.toWolfram();
}

std::string exportM82StarburstWolfram(double r_pc) {
    M82StarburstTerm term(r_pc);
    return term.toWolfram();
}

std::string exportM82SuperwindWolfram(double z_kpc) {
    M82SuperwindTerm term(z_kpc);
    return term.toWolfram();
}

std::string exportM81M82TidalWolfram(double r_kpc, double phi_rad) {
    M81M82TidalInteractionTerm term(r_kpc, phi_rad);
    return term.toWolfram();
}

std::string exportM81DarkMatterWolfram(double r_kpc) {
    M81DarkMatterHaloTerm term(r_kpc);
    return term.toWolfram();
}

std::string exportM82MolecularOutflowWolfram(double z_kpc) {
    M82MolecularOutflowTerm term(z_kpc);
    return term.toWolfram();
}

std::string exportM81AGNWolfram(double r_pc) {
    M81AGNActivityTerm term(r_pc);
    return term.toWolfram();
}

std::string exportM82PAHWolfram(double lambda_micron) {
    M82PAHEmissionTerm term(lambda_micron);
    return term.toWolfram();
}

std::string exportM81M82HIWolfram(double r_kpc) {
    M81M82HIDiskTerm term(r_kpc);
    return term.toWolfram();
}

std::string exportM82CosmicRayWolfram(double SN_rate) {
    M82CosmicRayTerm term(SN_rate);
    return term.toWolfram();
}

std::string exportAllM81M82WolframFunctions() {
    std::ostringstream oss;
    oss << "(* M81/M82 Galaxy Pair UQFF Module - Wolfram Language Export *)\n"
        << "(* Classes 730-739: M81 spiral + M82 starburst interaction *)\n\n"
        << exportM81SpiralWolfram(5.0, 0.0) << "\n\n"
        << exportM82StarburstWolfram(300.0) << "\n\n"
        << exportM82SuperwindWolfram(2.0) << "\n\n"
        << exportM81M82TidalWolfram(10.0, 0.0) << "\n\n"
        << exportM81DarkMatterWolfram(10.0) << "\n\n"
        << exportM82MolecularOutflowWolfram(1.0) << "\n\n"
        << exportM81AGNWolfram(100.0) << "\n\n"
        << exportM82PAHWolfram(7.7) << "\n\n"
        << exportM81M82HIWolfram(20.0) << "\n\n"
        << exportM82CosmicRayWolfram(0.1) << "\n";
    return oss.str();
}

// ========================================
// Master UQFF Integration Function
// ========================================

struct M81M82UQFFParams {
    double r_kpc;
    double phi_rad;
    double z_kpc;
    double t_Gyr;
    double SFR_M82;
    double M_M81;
    double M_M82;
    double d_M81M82_kpc;
    double v_wind;
    double M_BH_M81;
    double lambda_micron;
    double SN_rate_yr;
};

double computeM81M82MasterEquation(const M81M82UQFFParams& params) {
    M81SpiralStructureTerm m81_spiral(params.r_kpc, params.phi_rad, params.t_Gyr);
    M82StarburstTerm m82_sb(params.r_kpc * 1000.0, params.SFR_M82); // Convert kpc to pc
    M82SuperwindTerm m82_wind(params.z_kpc, params.v_wind);
    M81M82TidalInteractionTerm tidal(params.r_kpc, params.phi_rad, params.M_M81, params.M_M82, params.d_M81M82_kpc);
    M81DarkMatterHaloTerm m81_dm(params.r_kpc);
    M82MolecularOutflowTerm m82_mol_out(params.z_kpc);
    M81AGNActivityTerm m81_agn(params.r_kpc * 1000.0, params.M_BH_M81); // Convert kpc to pc
    M82PAHEmissionTerm m82_pah(params.lambda_micron);
    M81M82HIDiskTerm hi_disk(params.r_kpc);
    M82CosmicRayTerm m82_cr(params.SN_rate_yr);
    
    double F_m81_spiral = m81_spiral.compute();
    double F_m82_sb = m82_sb.compute();
    double F_m82_wind = m82_wind.compute();
    double F_tidal = tidal.compute();
    double F_m81_dm = m81_dm.compute();
    double F_m82_mol = m82_mol_out.compute();
    double F_m81_agn = m81_agn.compute();
    double F_m82_pah = m82_pah.compute();
    double F_hi = hi_disk.compute();
    double F_m82_cr = m82_cr.compute();
    
    // Master UQFF equation
    double F_total = F_m81_spiral + F_m82_sb + F_m82_wind + F_tidal + F_m81_dm + 
                     F_m82_mol + F_m81_agn + F_m82_pah + F_hi + F_m82_cr
                   + 0.1 * F_tidal * F_m82_sb         // Tidal-starburst coupling
                   + 0.05 * F_m82_sb * F_m82_wind     // Starburst-wind coupling
                   + 0.02 * F_m82_cr * F_m82_wind;    // CR-wind coupling
    
    return F_total;
}

// Example usage and validation
void validateM81M82Module() {
    M81M82UQFFParams params;
    params.r_kpc = 5.0;
    params.phi_rad = 0.0;
    params.z_kpc = 1.0;
    params.t_Gyr = 0.0;
    params.SFR_M82 = 10.0;           // M_sun/yr
    params.M_M81 = 7e10;             // M_sun
    params.M_M82 = 3e10;             // M_sun
    params.d_M81M82_kpc = 150.0;     // kpc
    params.v_wind = 750.0;           // km/s
    params.M_BH_M81 = 7e7;           // M_sun
    params.lambda_micron = 7.7;      // μm
    params.SN_rate_yr = 0.1;         // yr^-1
    
    double result = computeM81M82MasterEquation(params);
    
    // Expected: M81 spiral + M82 starburst + tidal interaction
}
