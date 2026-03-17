// HYDROGEN_ATOM_UQFF_MODULE.cpp â€” UQFF 2.0 Full Implementation
// Session 85: 27th C++ UQFF module â€” FIRST atomic-scale module; Hydrogen ground state (Bohr model)
// System: Proton M_p=1.6726e-27 kg, r_Bohr=5.2918e-11 m, electron v_orb=alpha*c=2.1877e6 m/s, z=0
// 9-term co-sum: base_expanded + Lorentz_EM + quantum_HUP + fluid_elec + osc_Lyman + GR_min + Lambda + Ug + DM
// g_H = (Î£ 9 terms) x (1 - B/B_crit) x (1 + f_TRZ)
//
// PAPER_299: Hydrogen Atom Electrogravitational Dominance Ratio
//            a_Lorentz = q*v_orb*B/m_e = 3.85e13 m/s^2 DOMINANT
//            g_base = G*M_p/r_Bohr^2 = 3.99e-17 m/s^2  (SMALLEST in all UQFF modules)
//            eta_EM = a_Lorentz/g_base = 9.65e29 (EM exceeds gravity by 30 orders at atomic scale)
//            FIRST UQFF atomic-scale module; FIRST eta_EM computation
// PAPER_300: Lyman-Alpha Cosmic Bridge â€” T/S = pi/13.8 = 0.2277 at atomic scale
//            omega_Lyman = 2pi*c/lambda = 1.549e16 rad/s (UV)
//            cosmic-age T/S ratio pi/T_U_gyr = pi/13.8 = 0.2277 (= PAPER_288 RSC value)
//            chi_bridge = omega_Lyman * t_H = 6.745e33 (Lyman-universe coupling factor)
//            FIRST UQFF demonstration that pi/T_U is universal at atomic scale
// PAPER_301: Proton GR Spectral Minimum â€” epsilon_GR = 3GM/(r*c^2) = 7.04e-44
//            r_S(proton) = 2.484e-54 m; r_Bohr/r_S = 2.13e43
//            UQFF GR spectral span: epsilon_GR(H) 7.04e-44 -> epsilon_GR(Universe PAPER_298) 5.056
//            Dynamic range: 7.18e43 (44 orders) â€” FIRST sub-Newtonian epsilon_GR module
//
// Copyright â€” Daniel T. Murphy, original Oct 09 2025. UQFF 2.0 upgrade Session 85.

#ifndef HYDROGEN_ATOM_UQFF_MODULE_H
#define HYDROGEN_ATOM_UQFF_MODULE_H

#include <cmath>
#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <complex>
#include <sstream>

// WOLFRAM_TERM markers (Wolfram KB / WSTP export)
#ifndef WOLFRAM_TERM
#define WOLFRAM_TERM(name, expr) /* name: expr */
#endif

WOLFRAM_TERM(HYDROGEN_BASE,
    "g_base = G*M_p/r_Bohr^2 = 3.99e-17 m/s^2; FIRST atomic UQFF; SMALLEST g_base (5 orders below M16 1.454e-12)")
WOLFRAM_TERM(HYDROGEN_LORENTZ,
    "a_Lorentz = q*v_orb*B/m_e = 3.85e13 m/s^2 [dominant]; eta_EM = a_Lorentz/g_base = 9.65e29 [PAPER_299]")
WOLFRAM_TERM(HYDROGEN_LYMAN,
    "a_osc = 2A*cos(omega_L*t) + (2pi/T_U_gyr)*A*cos(omega_L*t); T/S=pi/13.8=0.2277; chi_bridge=6.745e33 [PAPER_300]")
WOLFRAM_TERM(HYDROGEN_GR_MIN,
    "epsilon_GR = 3*G*M_p/(r_Bohr*c^2) = 7.04e-44; r_S=2.484e-54 m; GR span H->Universe = 7.18e43 [PAPER_301]")


// ============================================================================
class HydrogenAtomUQFFModule {
public:
    // â”€â”€ Static physical constants â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    static constexpr double C_LIGHT     = 2.99792e8;          // m/s
    static constexpr double G_NEWTON    = 6.6743e-11;         // m^3 kg^-1 s^-2
    static constexpr double HBAR        = 1.0546e-34;         // JÂ·s
    static constexpr double Q_ELEM      = 1.6022e-19;         // C  (electron charge magnitude)
    static constexpr double M_PROTON    = 1.6726e-27;         // kg (proton mass)
    static constexpr double M_ELEC      = 9.1094e-31;         // kg (electron mass)
    static constexpr double R_BOHR_DEF  = 5.2918e-11;         // m  (Bohr radius a_0)
    static constexpr double V_ORB_DEF   = 2.1877e6;           // m/s (alpha*c, alpha=1/137.036)
    static constexpr double LAMBDA_LY   = 1.2160e-7;          // m  (Lyman-alpha wavelength)
    static constexpr double PI          = 3.141592653589793;
    static constexpr double T_HUBBLE_S  = 13.8e9 * 3.15576e7; // s  (13.8 Gyr)
    static constexpr double T_COSMIC_GYR = 13.8;              // Gyr (for T/S ratio)
    static constexpr double ALPHA_FS    = 7.2974e-3;          // fine-structure constant 1/137.036

private:
    // â”€â”€ Named system parameters (UQFF 2.0 â€” NOT std::map) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    double r_Bohr;          // m    Bohr radius (orbital radius)
    double M_proton;        // kg   proton mass
    double m_elec;          // kg   electron mass
    double v_orb;           // m/s  electron orbital velocity (= ALPHA_FS * C_LIGHT)
    double B_atom;          // T    atomic internal magnetic field
    double B_crit;          // T    critical superconducting field
    double Delta_x;         // m    Heisenberg position uncertainty proxy (Compton scale)
    double rho_fluid;       // kg/m^3 electron cloud density estimate
    double delta_rho_frac;  // dimensionless cloud density contrast
    double A_osc;           // m/s^2 oscillation amplitude (wavefunction normalised)
    double omega_Lyman;     // rad/s Lyman-alpha angular frequency [PAPER_300]
    double k_Lyman;         // m^-1  Lyman-alpha wave vector       [PAPER_300]
    double H0_kms;          // km/s/Mpc Hubble constant
    double Omega_m;         // Omega_matter
    double Omega_L;         // Omega_Lambda
    double Lambda_cosm;     // m^-2  cosmological constant (negligible at atomic scale)
    double f_TRZ;           // TRZ correction factor
    double f_sc;            // superconducting scaling    (Ug4 = Ug1 * f_sc)

    // â”€â”€ Cached derived quantities â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    double g_base_cache;        // G*M_p/r_Bohr^2           [PAPER_299 reference]
    double a_Lorentz_cache;     // q*v_orb*B/m_e            [PAPER_299 dominant term]
    double eta_EM_cache;        // a_Lorentz / g_base       [PAPER_299: 9.65e29]
    double omega_L_cache;       // 2pi*c/lambda_Lyman        [PAPER_300]
    double k_L_cache;           // 2pi/lambda_Lyman          [PAPER_300]
    double chi_bridge_cache;    // omega_L * t_H             [PAPER_300: 6.745e33]
    double T_over_S_cache;      // pi / T_COSMIC_GYR         [PAPER_300: 0.2277]
    double epsilon_GR_cache;    // 3GM/(r*c^2)               [PAPER_301: 7.04e-44]
    double r_S_cache;           // 2GM/c^2                   [PAPER_301: 2.484e-54 m]
    double r_over_rS_cache;     // r_Bohr / r_S              [PAPER_301: 2.13e43]
    double Delta_p_cache;       // HBAR / Delta_x
    double V_orbital_cache;     // (4/3) pi r_Bohr^3

    // â”€â”€ UQFF 2.0 self-expanding fields â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

    // â”€â”€ Private compute helpers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
    void   updateCache();
    double computeHz()                    const;
    double computeBaseTerm(double t)      const;
    double computeLorentzEMTerm()         const;   // [PAPER_299] DOMINANT
    double computeHUPQuantumTerm()        const;
    double computeFluidElectronTerm()     const;
    double computeLymanResonantTerm(double t) const; // [PAPER_300]
    double computeGRMinTerm()             const;   // [PAPER_301]
    double computeLambdaCosmTerm()        const;
    double computeUgTerm()                const;   // Ug1 + Ug4 = 2*g_base
    double computeDMTerm()                const;   // = 0 (no DM in H atom)

public:
    // Constructor
    explicit HydrogenAtomUQFFModule();

    // UQFF 2.0 API
    void   setEnableLogging(bool enabled)              { logging_enabled = enabled; }
    bool   getLoggingEnabled() const                   { return logging_enabled; }
    void   setDynamicParameter(const std::string& key, double val) { dynamic_params[key] = val; }
    double getDynamicParameter(const std::string& key) const {
        auto it = dynamic_params.find(key);
        return (it != dynamic_params.end()) ? it->second
                                            : std::numeric_limits<double>::quiet_NaN();
    }

    // Core computation: g_H(t) â€” 9-term + SC x TRZ
    double computeG(double t);

    // State persistence
    void exportState(const std::string& filename) const;

    // Cross-validation template (header-only)
    template<typename OtherModule>
    double cross_validate(OtherModule& other, double t) {
        return std::abs(computeG(t) - other.computeG(t));
    }

    // Legacy variable interface (maps old string keys to named members)
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Descriptive output
    std::string getEquationText() const;
    void printVariables() const;
};


// â”€â”€ Constructor â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
HydrogenAtomUQFFModule::HydrogenAtomUQFFModule()
    : r_Bohr        (R_BOHR_DEF),
      M_proton      (M_PROTON),
      m_elec        (M_ELEC),
      v_orb         (V_ORB_DEF),                         // alpha * c = 2.1877e6 m/s
      B_atom        (1.0e-4),                             // T  internal atomic B estimate
      B_crit        (1.0e11),                             // T  Meissner critical field
      Delta_x       (1.0e-10),                            // m  Compton-scale uncertainty
      rho_fluid     (1.0e-25),                            // kg/m^3 electron cloud density
      delta_rho_frac(0.1),                                // 10% density contrast
      A_osc         (1.0e-10),                            // m/s^2 normalised amplitude
      omega_Lyman   (2.0 * PI * C_LIGHT / LAMBDA_LY),    // 1.549e16 rad/s [PAPER_300]
      k_Lyman       (2.0 * PI / LAMBDA_LY),              // 5.166e7 m^-1   [PAPER_300]
      H0_kms        (70.0),
      Omega_m       (0.3),
      Omega_L       (0.7),
      Lambda_cosm   (1.1e-52),
      f_TRZ         (0.1),
      f_sc          (1.0),
      g_base_cache  (0.0),
      a_Lorentz_cache(0.0),
      eta_EM_cache  (0.0),
      omega_L_cache (0.0),
      k_L_cache     (0.0),
      chi_bridge_cache(0.0),
      T_over_S_cache (0.0),
      epsilon_GR_cache(0.0),
      r_S_cache     (0.0),
      r_over_rS_cache(0.0),
      Delta_p_cache (0.0),
      V_orbital_cache(0.0),
      logging_enabled(false)
{
    updateCache();
}


// â”€â”€ updateCache â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
void HydrogenAtomUQFFModule::updateCache() {
    // Core Newtonian gravity at Bohr radius
    g_base_cache     = G_NEWTON * M_proton / (r_Bohr * r_Bohr);              // 3.99e-17 m/s^2

    // [PAPER_299] Lorentz acceleration of electron in atomic magnetic field
    a_Lorentz_cache  = Q_ELEM * v_orb * B_atom / m_elec;                     // 3.85e13  m/s^2
    eta_EM_cache     = (g_base_cache > 0.0)
                       ? a_Lorentz_cache / g_base_cache : 0.0;               // 9.65e29  [PAPER_299]

    // [PAPER_300] Lyman-alpha resonance quantities
    omega_L_cache    = omega_Lyman;                                           // 1.549e16 rad/s
    k_L_cache        = k_Lyman;                                               // 5.166e7  m^-1
    chi_bridge_cache = omega_L_cache * T_HUBBLE_S;                            // 6.745e33 [PAPER_300]
    T_over_S_cache   = PI / T_COSMIC_GYR;                                     // 0.2277   [PAPER_300]

    // [PAPER_301] GR curvature minimum at Bohr radius
    epsilon_GR_cache = 3.0 * G_NEWTON * M_proton
                       / (r_Bohr * C_LIGHT * C_LIGHT);                       // 7.04e-44 [PAPER_301]
    r_S_cache        = 2.0 * G_NEWTON * M_proton
                       / (C_LIGHT * C_LIGHT);                                 // 2.484e-54 m [PAPER_301]
    r_over_rS_cache  = (r_S_cache > 0.0) ? r_Bohr / r_S_cache : 0.0;         // 2.13e43  [PAPER_301]

    // Uncertainty momentum + orbital volume
    Delta_p_cache    = HBAR / Delta_x;
    V_orbital_cache  = (4.0 / 3.0) * PI * r_Bohr * r_Bohr * r_Bohr;
}


// â”€â”€ computeHz  â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
double HydrogenAtomUQFFModule::computeHz() const {
    double H0_si = H0_kms * 1.0e3 / 3.086e22;
    return H0_si * std::sqrt(Omega_m + Omega_L);    // z = 0 (atom is local)
}


// â”€â”€ Term 1: Hubble-expanded Newtonian base â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
double HydrogenAtomUQFFModule::computeBaseTerm(double t) const {
    double Hz = computeHz();
    return g_base_cache * (1.0 + Hz * t);            // ~4e-17 m/s^2 (negligible vs EM)
}


// â”€â”€ Term 2: Lorentz EM â€” DOMINANT [PAPER_299] â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
double HydrogenAtomUQFFModule::computeLorentzEMTerm() const {
    // a_Lorentz = q * v_orb * B / m_e = 3.85e13 m/s^2
    // This term dominates all others by 30 orders: eta_EM = 9.65e29
    // Represents the centripetal/Lorentz force on orbiting electron
    return a_Lorentz_cache;
}


// â”€â”€ Term 3: HUP Quantum (Heisenberg uncertainty) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
double HydrogenAtomUQFFModule::computeHUPQuantumTerm() const {
    // (hbar / sqrt(Delta_x * Delta_p)) * integral_psi * (2pi / t_H)
    double unc          = std::sqrt(Delta_x * Delta_p_cache);   // = sqrt(hbar) ~ 1.027e-17
    double integral_psi = 1.0;                                  // ground state normalised
    return (HBAR / unc) * integral_psi * (2.0 * PI / T_HUBBLE_S);
}


// â”€â”€ Term 4: Fluid electron cloud â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
double HydrogenAtomUQFFModule::computeFluidElectronTerm() const {
    // rho_fluid * V_orbital * g_base  (electron cloud density coupling)
    return rho_fluid * V_orbital_cache * g_base_cache;
}


// â”€â”€ Term 5: Lyman-alpha resonant oscillation [PAPER_300] â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
double HydrogenAtomUQFFModule::computeLymanResonantTerm(double t) const {
    // Standing wave:  2A * cos(k x)|_{x=0} * cos(omega_L * t) = 2A cos(omega_L t)
    double a_standing  = 2.0 * A_osc * std::cos(omega_L_cache * t);

    // Traveling wave with cosmic-age normalisation (T/S = pi/T_U_gyr = 0.2277)
    // Same T/S constant as PAPER_288 RSC module â€” universal across 27 frequency orders
    double a_traveling = T_over_S_cache * A_osc
                         * std::cos(k_L_cache * 0.0 - omega_L_cache * t);

    // [PAPER_300] chi_bridge = omega_L * t_H = 6.745e33 (Lyman-universe coupling)
    return a_standing + a_traveling;
}


// â”€â”€ Term 6: GR curvature minimum [PAPER_301] â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
double HydrogenAtomUQFFModule::computeGRMinTerm() const {
    // a_GR = g_base * epsilon_GR  (epsilon_GR = 7.04e-44 â€” UQFF GR spectral minimum)
    // Counterpart to PAPER_298 Universe epsilon_GR = 5.056 (44-order spectral range)
    return g_base_cache * epsilon_GR_cache;
}


// â”€â”€ Term 7: Cosmological Lambda (negligible at atomic scale) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
double HydrogenAtomUQFFModule::computeLambdaCosmTerm() const {
    return Lambda_cosm * C_LIGHT * C_LIGHT / 3.0;   // ~3.3e-36 m/s^2
}


// â”€â”€ Term 8: Ug1 + Ug4 co-sum (26-layer triadic) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
double HydrogenAtomUQFFModule::computeUgTerm() const {
    // Ug1 = G*M/r^2 = g_base;  Ug4 = Ug1 * f_sc;  Ug2 = Ug3 = 0 (charge/string ~0 for proton)
    double Ug1 = g_base_cache;
    double Ug4 = Ug1 * f_sc;
    return Ug1 + Ug4;    // = 2 * g_base  at f_sc=1
}


// â”€â”€ Term 9: DM perturbation (negligible â€” M_DM = 0 for hydrogen) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
double HydrogenAtomUQFFModule::computeDMTerm() const {
    // Standard: g_base * f_DM * delta; f_DM=0 for isolated hydrogen atom
    return 0.0;
}


// â”€â”€ computeG: 9-term co-sum + SC x TRZ â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
double HydrogenAtomUQFFModule::computeG(double t) {
    double a_base     = computeBaseTerm(t);           // ~4.0e-17  m/s^2
    double a_Lorentz  = computeLorentzEMTerm();       // ~3.85e13  m/s^2  [PAPER_299 DOMINANT]
    double a_quantum  = computeHUPQuantumTerm();      // ~1.5e-34  m/s^2
    double a_fluid    = computeFluidElectronTerm();   // ~2.5e-72  m/s^2
    double a_osc      = computeLymanResonantTerm(t);  // ~2.0e-10  m/s^2  [PAPER_300]
    double a_GR_min   = computeGRMinTerm();           // ~2.8e-60  m/s^2  [PAPER_301]
    double a_Lambda   = computeLambdaCosmTerm();      // ~3.3e-36  m/s^2
    double a_Ug       = computeUgTerm();              // ~8.0e-17  m/s^2
    double a_DM       = computeDMTerm();              // = 0

    double sigma = a_base + a_Lorentz + a_quantum + a_fluid
                 + a_osc + a_GR_min + a_Lambda + a_Ug + a_DM;

    double SCm   = 1.0 - B_atom / B_crit;
    double g_H   = sigma * SCm * (1.0 + f_TRZ);

    if (logging_enabled) {
        std::cout << "[LOG] HydrogenAtomUQFFModule computeG(t=" << std::scientific << t << ")\n"
                  << "  a_base     = " << std::setw(14) << a_base
                  << " m/s^2  (Newtonian; FIRST atomic UQFF g_base=3.99e-17)\n"
                  << "  a_Lorentz  = " << std::setw(14) << a_Lorentz
                  << " m/s^2  [PAPER_299 DOMINANT; eta_EM=" << eta_EM_cache << "]\n"
                  << "  a_quantum  = " << std::setw(14) << a_quantum
                  << " m/s^2  (HUP Heisenberg)\n"
                  << "  a_fluid    = " << std::setw(14) << a_fluid
                  << " m/s^2  (electron cloud)\n"
                  << "  a_osc_Ly   = " << std::setw(14) << a_osc
                  << " m/s^2  [PAPER_300 Lyman-alpha; T/S=" << T_over_S_cache << "]\n"
                  << "  a_GR_min   = " << std::setw(14) << a_GR_min
                  << " m/s^2  [PAPER_301 eps_GR=" << epsilon_GR_cache << "]\n"
                  << "  a_Lambda   = " << std::setw(14) << a_Lambda
                  << " m/s^2  (cosmological; negligible)\n"
                  << "  a_Ug       = " << std::setw(14) << a_Ug
                  << " m/s^2  (Ug1+Ug4 triadic co-sum)\n"
                  << "  a_DM       = " << std::setw(14) << a_DM
                  << " m/s^2  (M_DM=0; negligible)\n"
                  << "  sigma      = " << std::setw(14) << sigma  << " m/s^2\n"
                  << "  SCm        = " << SCm  << "\n"
                  << "  g_H        = " << g_H  << " m/s^2\n";
    }

    return g_H;
}


// â”€â”€ exportState â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
void HydrogenAtomUQFFModule::exportState(const std::string& filename) const {
    std::ofstream f(filename);
    if (!f) return;
    f << std::scientific << std::setprecision(6);
    f << "# HydrogenAtomUQFFModule State Export â€” Session 85 â€” UQFF 2.0\n";
    f << "# 27th C++ UQFF module; FIRST atomic-scale\n";
    // System parameters
    f << "r_Bohr_m             = " << r_Bohr          << "\n";
    f << "M_proton_kg          = " << M_proton         << "\n";
    f << "m_elec_kg            = " << m_elec           << "\n";
    f << "v_orb_m_s            = " << v_orb            << "\n";  // alpha*c
    f << "B_atom_T             = " << B_atom           << "\n";
    f << "B_crit_T             = " << B_crit           << "\n";
    f << "Delta_x_m            = " << Delta_x          << "\n";
    f << "Delta_p_kg_m_s       = " << Delta_p_cache    << "\n";
    f << "rho_fluid_kg_m3      = " << rho_fluid        << "\n";
    f << "delta_rho_frac       = " << delta_rho_frac   << "\n";
    f << "V_orbital_m3         = " << V_orbital_cache  << "\n";
    f << "A_osc_m_s2           = " << A_osc            << "\n";
    f << "Lambda_cosm_m-2      = " << Lambda_cosm      << "\n";
    f << "H0_kms_Mpc           = " << H0_kms           << "\n";
    f << "Omega_m              = " << Omega_m          << "\n";
    f << "Omega_L              = " << Omega_L          << "\n";
    f << "f_TRZ                = " << f_TRZ            << "\n";
    f << "f_sc                 = " << f_sc             << "\n";
    // PAPER_299 â€” Electrogravitational dominance
    f << "g_base_m_s2          = " << g_base_cache     << "\n";  // 3.99e-17
    f << "a_Lorentz_m_s2       = " << a_Lorentz_cache  << "\n";  // 3.85e13  DOMINANT
    f << "eta_EM               = " << eta_EM_cache     << "\n";  // 9.65e29  [PAPER_299]
    f << "eta_EM_log10         = " << std::log10(eta_EM_cache > 0.0 ? eta_EM_cache : 1.0) << "\n";
    // PAPER_300 â€” Lyman-Alpha Cosmic Bridge
    f << "omega_Lyman_rad_s    = " << omega_L_cache    << "\n";  // 1.549e16
    f << "k_Lyman_m-1          = " << k_L_cache        << "\n";  // 5.166e7
    f << "chi_bridge_omega_tH  = " << chi_bridge_cache << "\n";  // 6.745e33
    f << "T_over_S_ratio       = " << T_over_S_cache   << "\n";  // 0.2277 = pi/13.8
    f << "scale_TS_from_PAPER288 = " << T_over_S_cache / (PI / T_COSMIC_GYR) << "\n"; // = 1.0
    // PAPER_301 â€” GR Spectral Minimum
    f << "epsilon_GR           = " << epsilon_GR_cache << "\n";  // 7.04e-44
    f << "r_S_proton_m         = " << r_S_cache        << "\n";  // 2.484e-54
    f << "r_Bohr_over_r_S      = " << r_over_rS_cache  << "\n";  // 2.13e43
    f << "eps_GR_Universe      = " << 5.056            << "\n";  // PAPER_298 reference
    f << "eps_GR_spectral_span = " << 5.056 / epsilon_GR_cache << "\n";  // 7.18e43
    f << "eps_GR_log10_span    = " << std::log10(5.056 / epsilon_GR_cache) << "\n"; // ~43.9
    // Dynamic params
    for (const auto& kv : dynamic_params)
        f << "dynamic." << kv.first << " = " << kv.second << "\n";
}


// â”€â”€ Legacy variable interface â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
void HydrogenAtomUQFFModule::updateVariable(const std::string& name, double value) {
    bool need_cache = true;
    if      (name == "r"         || name == "r_Bohr")   { r_Bohr    = value; }
    else if (name == "M"         || name == "M_proton") { M_proton  = value; }
    else if (name == "m_elec")                          { m_elec    = value; }
    else if (name == "v_orbital" || name == "v_orb")    { v_orb     = value; }
    else if (name == "B")                               { B_atom    = value; need_cache = false; }
    else if (name == "B_crit")                          { B_crit    = value; need_cache = false; }
    else if (name == "Delta_x")                         { Delta_x   = value; }
    else if (name == "rho_fluid")                       { rho_fluid = value; }
    else if (name == "A" || name == "A_osc")            { A_osc     = value; need_cache = false; }
    else if (name == "omega")                           { omega_Lyman = value; }
    else if (name == "k")                               { k_Lyman   = value; need_cache = false; }
    else if (name == "H0")                              { H0_kms    = value; need_cache = false; }
    else if (name == "f_TRZ")                           { f_TRZ     = value; need_cache = false; }
    else if (name == "f_sc")                            { f_sc      = value; need_cache = false; }
    else if (name == "Lambda")                          { Lambda_cosm = value; need_cache = false; }
    else { setDynamicParameter(name, value); need_cache = false; }
    if (need_cache) updateCache();
}

void HydrogenAtomUQFFModule::addToVariable(const std::string& name, double delta) {
    if      (name == "r"    || name == "r_Bohr")    updateVariable("r_Bohr",    r_Bohr    + delta);
    else if (name == "M"    || name == "M_proton")  updateVariable("M_proton",  M_proton  + delta);
    else if (name == "v_orb"|| name == "v_orbital") updateVariable("v_orb",     v_orb     + delta);
    else if (name == "B")    { B_atom  += delta; }
    else if (name == "f_TRZ"){ f_TRZ   += delta; }
    else {
        double cur = getDynamicParameter(name);
        setDynamicParameter(name, (std::isnan(cur) ? 0.0 : cur) + delta);
    }
}

void HydrogenAtomUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}


// â”€â”€ getEquationText â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
std::string HydrogenAtomUQFFModule::getEquationText() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "g_H(r,t) = [a_base + a_Lorentz + a_quantum + a_fluid + a_osc + a_GR_min + a_Lambda + a_Ug + a_DM]\n"
        << "           x (1 - B/B_crit) x (1 + f_TRZ)\n"
        << "\nWOLFRAM_TERM HYDROGEN_BASE:\n"
        << "  g_base = G*M_p/r_Bohr^2 = " << g_base_cache << " m/s^2\n"
        << "  FIRST atomic UQFF module; SMALLEST g_base in all 27 modules\n"
        << "\nWOLFRAM_TERM HYDROGEN_LORENTZ [PAPER_299]:\n"
        << "  a_Lorentz = q*v_orb*B/m_e = " << a_Lorentz_cache << " m/s^2  DOMINANT\n"
        << "  eta_EM = a_Lorentz/g_base = " << eta_EM_cache << "  (EM/gravity = 9.65e29)\n"
        << "\nWOLFRAM_TERM HYDROGEN_LYMAN [PAPER_300]:\n"
        << "  a_osc = 2A cos(omega_L t) + (pi/13.8)*A cos(omega_L t)\n"
        << "  omega_Lyman = " << omega_L_cache << " rad/s\n"
        << "  T/S = pi/T_U_gyr = " << T_over_S_cache << "  (= PAPER_288 RSC: pi/13.8)\n"
        << "  chi_bridge = omega_L * t_H = " << chi_bridge_cache << "\n"
        << "\nWOLFRAM_TERM HYDROGEN_GR_MIN [PAPER_301]:\n"
        << "  epsilon_GR = 3GM/(r*c^2) = " << epsilon_GR_cache << "  (GR spectral minimum)\n"
        << "  r_S = " << r_S_cache << " m  (proton Schwarzschild radius)\n"
        << "  r_Bohr/r_S = " << r_over_rS_cache << "\n"
        << "  GR spectral span (H -> Universe) = " << 5.056/epsilon_GR_cache << "  (44 orders)\n";
    return oss.str();
}


// â”€â”€ printVariables â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
void HydrogenAtomUQFFModule::printVariables() const {
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "=== HydrogenAtomUQFFModule UQFF 2.0 â€” Session 85 (27th C++ module) ===\n";
    std::cout << "  System: Hydrogen ground state, Bohr model, z=0\n";
    std::cout << "  r_Bohr          = " << r_Bohr          << " m   (Bohr radius)\n";
    std::cout << "  M_proton        = " << M_proton         << " kg\n";
    std::cout << "  m_elec          = " << m_elec           << " kg\n";
    std::cout << "  v_orb           = " << v_orb            << " m/s  (= alpha*c)\n";
    std::cout << "  B_atom          = " << B_atom           << " T\n";
    std::cout << "  B_crit          = " << B_crit           << " T\n";
    std::cout << "  Delta_x         = " << Delta_x          << " m\n";
    std::cout << "  Delta_p         = " << Delta_p_cache    << " kg*m/s\n";
    std::cout << "  omega_Lyman     = " << omega_L_cache    << " rad/s  [PAPER_300]\n";
    std::cout << "  k_Lyman         = " << k_L_cache        << " m^-1   [PAPER_300]\n";
    std::cout << "  f_TRZ           = " << f_TRZ            << "\n";
    std::cout << "  f_sc            = " << f_sc             << "\n";
    std::cout << "  --- PAPER_299: Electrogravitational Dominance ---\n";
    std::cout << "  g_base          = " << g_base_cache     << " m/s^2  (3.99e-17; UQFF minimum)\n";
    std::cout << "  a_Lorentz       = " << a_Lorentz_cache  << " m/s^2  (3.85e13; DOMINANT)\n";
    std::cout << "  eta_EM          = " << eta_EM_cache     << "  (9.65e29; 30 orders)\n";
    std::cout << "  --- PAPER_300: Lyman-Alpha Cosmic Bridge ---\n";
    std::cout << "  chi_bridge      = " << chi_bridge_cache << "  (omega_L * t_H = 6.745e33)\n";
    std::cout << "  T/S ratio       = " << T_over_S_cache   << "  (pi/13.8 = 0.2277 = PAPER_288)\n";
    std::cout << "  --- PAPER_301: GR Spectral Minimum ---\n";
    std::cout << "  epsilon_GR      = " << epsilon_GR_cache << "  (7.04e-44; UQFF minimum)\n";
    std::cout << "  r_S_proton      = " << r_S_cache        << " m   (2.484e-54)\n";
    std::cout << "  r_Bohr / r_S    = " << r_over_rS_cache  << "  (2.13e43)\n";
    std::cout << "  GR span H->U    = " << 5.056/epsilon_GR_cache << "  (7.18e43; 44 orders)\n";
}

#endif // HYDROGEN_ATOM_UQFF_MODULE_H
