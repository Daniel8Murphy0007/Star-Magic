// UNIVERSE_DIAMETER_UQFF_MODULE.cpp
// UQFF 2.0 Observable Universe Diameter Module — Session 84 (Universe as System)
// 26th C++ UQFF module — FIRST dedicated Universe-as-system scale module.
// System: Observable Universe — r=4.4×10²⁶ m, M=1×10⁵⁴ kg, H₀=70 km/s/Mpc, z=0
// 9-term co-sum: base_expanded + Λ + quantum + EM + fluid + osc + DM_pert + GR_curv
// g_UVDIAM(r,t) = (all 9 terms) × (1-B/B_crit) × (1+f_TRZ)
// PAPER_296: a_Λ = Λc²/3 = 3.30×10⁻³⁶ m/s² — FIRST UQFF explicit dark-energy vacuum acceleration
// PAPER_297: η_exp = v_exp/c = H₀×r/c = 3.328 > 1 — FIRST UQFF superluminal Hubble expansion ratio
// PAPER_298: ε_GR = 3GM/(rc²) = 5.056 > 1 — FIRST UQFF GR curvature dominance (all 25 prior: ε_GR << 1)
// Watermark: Copyright - Daniel T. Murphy, upgraded UQFF 2.0 Session 84 (March 17, 2026)

#ifndef UNIVERSE_DIAMETER_UQFF_MODULE_H
#define UNIVERSE_DIAMETER_UQFF_MODULE_H

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>

// ---------------------------------------------------------------------------
// WOLFRAM_TERM macros — 4 unique physics anchors (Session 84, PAPER_296–298)
// ---------------------------------------------------------------------------
#define WOLFRAM_TERM_UVDIAM_BASE \
    "g_UVDIAM(r,t)=(g_base*(1+H(z)*t)+a_Lambda+a_q+a_EM+a_fluid+a_osc+a_DM_pert+a_GR)*(1-B/B_crit)*(1+f_TRZ);9-term universe-scale co-sum;g_base=GM/r2=3.447e-10;H0=70km/s/Mpc=2.269e-18s-1;M=1e54kg;r=4.4e26m"

#define WOLFRAM_TERM_UVDIAM_LAMBDA \
    "a_Lambda=Lambda*c2/3=1.1e-52*9e16/3=3.30e-36m/s2;FIRST UQFF explicit dark-energy term;Gamma_Lambda=a_Lambda/g_base=9.57e-27;d_Lambda=0.5*a_Lambda*t_H2=0.313m cosmic displacement;all 25 prior modules Lambda implicit in H(z) only [PAPER_296]"

#define WOLFRAM_TERM_UVDIAM_HUBBLE \
    "eta_exp=v_exp/c=H0*r/c=9.984e8/3e8=3.328>1;FIRST UQFF eta_exp>1;r_H=c/H0=1.322e26m Hubble sphere;r_obs=3.328*r_H;expansion_factor(t_H)=1+H*t=1.988(near-doubling);EM term a_EM prop eta_exp [PAPER_297]"

#define WOLFRAM_TERM_UVDIAM_GR_CURV \
    "epsilon_GR=3GM/(r*c2)=3*6.674e-11*1e54/(4.4e26*9e16)=5.056>1;FIRST UQFF epsilon_GR>1;a_GR=g_base*epsilon_GR=1.743e-9m/s2(5x Newtonian!);r_S/r_obs=2*epsilon_GR/3=3.371;r_obs=0.297*r_S;all 25 prior UQFF epsilon_GR<<1 [PAPER_298]"

// ===========================================================================
// Class declaration
// ===========================================================================

class UniverseDiameterUQFFModule {
private:
    // --- Universal physical constants ---
    static constexpr double C_LIGHT    = 3.0e8;           // m/s
    static constexpr double G_NEWTON   = 6.6743e-11;      // m³/(kg·s²)
    static constexpr double HBAR       = 1.0546e-34;      // J·s reduced Planck
    static constexpr double Q_ELEM     = 1.602e-19;       // C elementary charge
    static constexpr double M_PROTON   = 1.673e-27;       // kg proton mass
    static constexpr double PI         = 3.141592653589793;
    static constexpr double T_HUBBLE_S = 13.8e9 * 3.15576e7; // s = 4.355e17 s (canonical)
    static constexpr double T_COSMIC_GYR = 13.8;          // Gyr — PAPER_288 cosmic-age normaliser

    // --- Observable Universe parameters ---
    double r_obs;          // m = 4.4e26   half-diameter of observable universe
    double M_total;        // kg = 1.0e54  total matter+DM (from ρ_c × Ω_m × V_obs)
    double f_baryon;       // = 0.05       baryonic fraction of M_total
    double f_DM;           // = 0.27       dark matter fraction of M_total
    double rho_c;          // kg/m³ = 9.21e-27   critical density at z=0
    double delta_contrast; // = 0.1        δρ/ρ DM density perturbation contrast

    // --- Cosmological parameters ---
    double H0_kms;         // km/s/Mpc = 70.0
    double H0_si;          // s⁻¹  = 2.269e-18
    double Lambda;         // m⁻² = 1.1e-52 (cosmological constant Λ)
    double Omega_m;        // = 0.3
    double Omega_L;        // = 0.7
    double z_obs;          // = 0.0 (observable universe at z=0)

    // --- EM / quantum / oscillatory parameters ---
    double B_cosmic;       // T = 1e-15    cosmic magnetic field
    double B_crit;         // T = 1e11
    double Delta_x;        // m = 1e-10    cosmic quantum scale proxy
    double A_cmb;          // m = 1e-10    CMB / large-scale oscillation amplitude
    double k_cmb;          // m⁻¹ = 1e20  CMB wave vector
    double omega_cmb;      // rad/s = 1e11 CMB / large-scale structure frequency
    double scale_EM;       // = 1e-12      macro scale factor for cosmic EM term

    // --- UQFF correction factors ---
    double f_TRZ;          // = 0.1
    double f_sc;           // = 1.0

    // --- Cached computed quantities (set by updateCache) ---
    double g_base_cache;       // m/s²  = 3.447e-10  G×M/r²
    double Hz_cache;           // s⁻¹   = 2.269e-18  H(z=0)
    double v_exp_cache;        // m/s   = 9.984e8    H₀×r_obs (superluminal)
    double Delta_p_cache;      // kg⋅m/s = 1.0546e-24  ħ/Δx
    double a_Lambda_cache;     // m/s²  = 3.30e-36   Λ×c²/3  [PAPER_296]
    double eta_exp_cache;      // dimensionless = 3.328   v_exp/c  [PAPER_297]
    double r_H_cache;          // m     = 1.322e26   Hubble sphere c/H₀
    double epsilon_GR_cache;   // dimensionless = 5.056   3GM/(r×c²)  [PAPER_298]
    double a_GR_cache;         // m/s²  = 1.743e-9   g_base×ε_GR  [PAPER_298]

    // --- UQFF 2.0 runtime infrastructure ---
    bool logging_enabled;
    std::map<std::string, double> dynamic_params;

    // --- Private compute helpers ---
    double computeBaseTerm(double t) const;        // g_base × (1 + H(z) × t) Newtonian + Hubble coupling
    double computeLambdaTerm() const;              // [PAPER_296]  Λ×c²/3 dark energy vacuum acceleration
    double computeQuantumTerm() const;             // ħ/√(Δx⋅Δp) × 2π/t_H Heisenberg cosmic uncertainty
    double computeEMLorentzTerm() const;           // [PAPER_297]  q×v_exp×B/m_p × (1+η_exp) superluminal EM
    double computeFluidJeansTerm() const;          // ρ_c × G × δ × r_obs  Jeans density drive
    double computeResonantCMBTerm(double t) const; // PAPER_288 standing+traveling pair (CMB-scale)
    double computeDMPertTerm() const;              // g_base × f_DM × δ_contrast  DM density perturbation
    double computeGRCurvatureTerm() const;         // [PAPER_298]  g_base × ε_GR  GR curvature dominance

    double computeHz(double z_val) const;
    void   updateCache();

public:
    UniverseDiameterUQFFModule();

    // --- Core computation (legacy API preserved) ---
    double computeG(double t);

    // --- Legacy variable and add/subtract operations ---
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // --- Descriptive output ---
    std::string getEquationText() const;
    void printVariables() const;

    // --- UQFF 2.0 runtime configuration ---
    void   setEnableLogging(bool enable)                         { logging_enabled = enable; }
    bool   getLoggingEnabled()                             const { return logging_enabled; }
    void   setDynamicParameter(const std::string& nm, double v)  { dynamic_params[nm] = v; }
    double getDynamicParameter(const std::string& nm)      const {
        auto it = dynamic_params.find(nm);
        return (it != dynamic_params.end()) ? it->second
                                            : std::numeric_limits<double>::quiet_NaN();
    }

    // --- State persistence ---
    void exportState(const std::string& filename =
                     "UniverseDiameterUQFFModule_UQFF_state.txt") const;

    // --- Cross-validation template (UQFF 2.0) ---
    template <typename OtherModule>
    double cross_validate(OtherModule& other, double t) {
        double g_self  = computeG(t);
        double g_other = other.computeG(t);
        double ratio   = (g_other != 0.0) ? g_self / g_other : 0.0;
        if (logging_enabled) {
            std::cout << "[UVDIAM XVAL] g_self=" << std::scientific << g_self
                      << " g_other=" << g_other
                      << " ratio="   << ratio << "\n";
        }
        return ratio;
    }
};

// ===========================================================================
// Implementation
// ===========================================================================

double UniverseDiameterUQFFModule::computeHz(double z_val) const {
    // Friedmann equation: H(z) = H₀ × √(Ω_m×(1+z)³ + Ω_Λ)
    return H0_si * std::sqrt(Omega_m * std::pow(1.0 + z_val, 3.0) + Omega_L);
}

void UniverseDiameterUQFFModule::updateCache() {
    Hz_cache        = computeHz(z_obs);
    g_base_cache    = G_NEWTON * M_total / (r_obs * r_obs);
    v_exp_cache     = H0_si * r_obs;
    Delta_p_cache   = HBAR / Delta_x;
    // [PAPER_296] a_Λ = Λ × c² / 3
    a_Lambda_cache  = Lambda * C_LIGHT * C_LIGHT / 3.0;
    // [PAPER_297] η_exp = v_exp / c — Universe boundary recedes superluminally
    eta_exp_cache   = v_exp_cache / C_LIGHT;
    r_H_cache       = C_LIGHT / H0_si;
    // [PAPER_298] ε_GR = 3GM/(r×c²) — post-Newtonian curvature parameter
    epsilon_GR_cache = 3.0 * G_NEWTON * M_total / (r_obs * C_LIGHT * C_LIGHT);
    a_GR_cache      = g_base_cache * epsilon_GR_cache;
}

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

UniverseDiameterUQFFModule::UniverseDiameterUQFFModule()
    : r_obs(4.4e26), M_total(1.0e54),
      f_baryon(0.05), f_DM(0.27),
      rho_c(9.21e-27), delta_contrast(0.1),
      H0_kms(70.0), H0_si(70.0e3 / 3.086e22),
      Lambda(1.1e-52), Omega_m(0.3), Omega_L(0.7), z_obs(0.0),
      B_cosmic(1.0e-15), B_crit(1.0e11),
      Delta_x(1.0e-10),
      A_cmb(1.0e-10), k_cmb(1.0e20), omega_cmb(1.0e11),
      scale_EM(1.0e-12),
      f_TRZ(0.1), f_sc(1.0),
      g_base_cache(0.0), Hz_cache(0.0), v_exp_cache(0.0),
      Delta_p_cache(0.0), a_Lambda_cache(0.0),
      eta_exp_cache(0.0), r_H_cache(0.0),
      epsilon_GR_cache(0.0), a_GR_cache(0.0),
      logging_enabled(false)
{
    updateCache();
}

// ---------------------------------------------------------------------------
// Private compute helpers
// ---------------------------------------------------------------------------

double UniverseDiameterUQFFModule::computeBaseTerm(double t) const {
    // Newtonian base + Hubble expansion coupling: g_base × (1 + H(z)×t)
    // At t=t_H: expansion factor = 1.9882 (near-doubling of base gravity over universe age)
    const double expansion = 1.0 + Hz_cache * t;
    return g_base_cache * expansion;
}

double UniverseDiameterUQFFModule::computeLambdaTerm() const {
    // [PAPER_296] a_Λ = Λ × c² / 3 = 3.30e-36 m/s²
    // FIRST UQFF explicit dark-energy vacuum acceleration term.
    // All 25 prior modules: Λ only implicit inside H(z) Friedmann equation.
    // Γ_Λ = a_Λ / g_base = 9.57e-27 (dark energy 27 orders below gravity)
    // Cosmic displacement: d_Λ = ½ × a_Λ × t_H² = 0.313 m (macroscopic cumulative)
    return a_Lambda_cache;
}

double UniverseDiameterUQFFModule::computeQuantumTerm() const {
    // Heisenberg cosmic quantum uncertainty coupled to Hubble time
    // a_q = (ħ / √(Δx × Δp)) × 2π / t_Hubble
    const double uncertainty = std::sqrt(Delta_x * Delta_p_cache);
    return (HBAR / uncertainty) * (2.0 * PI / T_HUBBLE_S);
}

double UniverseDiameterUQFFModule::computeEMLorentzTerm() const {
    // [PAPER_297] a_EM = (q × v_exp × B) / m_p × (1 + η_exp) × scale_EM
    // v_exp = H₀ × r_obs = 9.984e8 m/s = 3.328c (superluminal Hubble flow)
    // η_exp = 3.328 > 1 — FIRST UQFF superluminal expansion parameter.
    // Factor (1 + η_exp) = 4.328 encodes the super-c Hubble scaling explicitly.
    // Physical: Lorentz force on CMB photon-equivalent charge at boundary recession velocity.
    return (Q_ELEM * v_exp_cache * B_cosmic / M_PROTON)
           * (1.0 + eta_exp_cache) * scale_EM;
}

double UniverseDiameterUQFFModule::computeFluidJeansTerm() const {
    // Jeans-mode density perturbation acceleration:
    // a_fluid = ρ_c × G × δ_contrast × r_obs
    // Encodes the large-scale structure gravitational instability drive.
    return rho_c * G_NEWTON * delta_contrast * r_obs;
}

double UniverseDiameterUQFFModule::computeResonantCMBTerm(double t) const {
    // PAPER_288-inherited standing + traveling wave pair (CMB-scale extension)
    // a_osc = 2A cos(kx) cos(ωt) + (2π/T_13.8) × A × Re[exp(i(kx - ωt))]
    // x = 0 (centre-of-universe reference frame)
    const double a_standing  = 2.0 * A_cmb
                               * std::cos(k_cmb * 0.0) * std::cos(omega_cmb * t);
    std::complex<double> phase(0.0, k_cmb * 0.0 - omega_cmb * t);
    const double a_traveling = (2.0 * PI / T_COSMIC_GYR) * A_cmb
                               * std::exp(phase).real();
    return a_standing + a_traveling;
}

double UniverseDiameterUQFFModule::computeDMPertTerm() const {
    // Dark matter density perturbation contribution:
    // a_DM_pert = g_base × f_DM × δ_contrast
    return g_base_cache * f_DM * delta_contrast;
}

double UniverseDiameterUQFFModule::computeGRCurvatureTerm() const {
    // [PAPER_298] Post-Newtonian GR curvature: a_GR = g_base × ε_GR
    // ε_GR = 3GM/(r × c²) = 5.056 > 1 — FIRST UQFF GR curvature dominance.
    // At Universe scale the GR correction (1.743e-9 m/s²) EXCEEDS the Newtonian base.
    // All 25 prior UQFF modules: ε_GR << 1 (Saturn: 1.4e-8; Andromeda: 2.8e-6).
    // r_S / r_obs = 2×ε_GR/3 = 3.371 → observable universe is at ~30% of its
    // own Schwarzschild radius — the only UQFF system in the GR-Dominant Regime.
    return a_GR_cache;
}

// ---------------------------------------------------------------------------
// Main computation — UQFF 2.0 9-term co-sum
// ---------------------------------------------------------------------------

double UniverseDiameterUQFFModule::computeG(double t) {
    // === 9-term co-sum ===
    const double a_base     = computeBaseTerm(t);
    const double a_Lambda   = computeLambdaTerm();            // [PAPER_296] Λc²/3
    const double a_quantum  = computeQuantumTerm();
    const double a_EM       = computeEMLorentzTerm();         // [PAPER_297] η_exp > 1
    const double a_fluid    = computeFluidJeansTerm();
    const double a_osc      = computeResonantCMBTerm(t);
    const double a_DM_pert  = computeDMPertTerm();
    const double a_GR       = computeGRCurvatureTerm();       // [PAPER_298] ε_GR > 1

    const double sum_all    = a_base + a_Lambda + a_quantum + a_EM
                            + a_fluid + a_osc + a_DM_pert + a_GR;

    // === SC correction + TRZ global multiplier ===
    const double SCm        = 1.0 - B_cosmic / B_crit;
    const double g_UVDIAM   = sum_all * SCm * (1.0 + f_TRZ);

    if (logging_enabled) {
        std::cout << std::scientific << std::setprecision(6)
                  << "[UVDIAM][BASE]      a_base="    << a_base
                  << "  (g_base=" << g_base_cache
                  << "  expansion=" << (1.0 + Hz_cache * t) << ")\n"
                  << "[UVDIAM][LAMBDA]    a_Lambda="  << a_Lambda
                  << "  Gamma_Lambda=" << (g_base_cache != 0.0 ? a_Lambda / g_base_cache : 0.0)
                  << "  [PAPER_296 FIRST UQFF dark-energy term]\n"
                  << "[UVDIAM][QUANTUM]   a_quantum=" << a_quantum << "\n"
                  << "[UVDIAM][EM]        a_EM="      << a_EM
                  << "  eta_exp=" << eta_exp_cache
                  << "  v_exp=" << v_exp_cache << " m/s"
                  << "  [PAPER_297 eta_exp>1]\n"
                  << "[UVDIAM][FLUID]     a_fluid="   << a_fluid << "\n"
                  << "[UVDIAM][OSC]       a_osc="     << a_osc << "\n"
                  << "[UVDIAM][DM_PERT]   a_DM_pert=" << a_DM_pert << "\n"
                  << "[UVDIAM][GR_CURV]   a_GR="      << a_GR
                  << "  epsilon_GR=" << epsilon_GR_cache
                  << "  [PAPER_298 FIRST UQFF GR>Newton]\n"
                  << "[UVDIAM][SC]        SCm="       << SCm
                  << "  B=" << B_cosmic << " B_crit=" << B_crit << "\n"
                  << "[UVDIAM][RESULT]    g_UVDIAM="  << g_UVDIAM << " m/s^2\n";
    }

    return g_UVDIAM;
}

// ---------------------------------------------------------------------------
// UQFF 2.0 state export (~38 parameters)
// ---------------------------------------------------------------------------

void UniverseDiameterUQFFModule::exportState(const std::string& filename) const {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "[UVDIAM] Cannot open state file: " << filename << "\n";
        return;
    }
    ofs << std::scientific << std::setprecision(10);
    ofs << "# UniverseDiameterUQFFModule — UQFF 2.0 State (Session 84)\n";
    ofs << "# Papers: PAPER_296 PAPER_297 PAPER_298\n";
    ofs << "# Module:  26th C++ UQFF module — FIRST Universe-as-system scale module\n";
    ofs << "C_LIGHT              = " << C_LIGHT            << "  # m/s\n";
    ofs << "G_NEWTON             = " << G_NEWTON           << "  # m^3/(kg*s^2)\n";
    ofs << "HBAR                 = " << HBAR               << "  # J*s\n";
    ofs << "T_HUBBLE_S           = " << T_HUBBLE_S         << "  # s canonical 13.8 Gyr\n";
    ofs << "r_obs                = " << r_obs              << "  # m half-diameter obs universe\n";
    ofs << "M_total              = " << M_total            << "  # kg from rho_c*Omega_m*V_obs\n";
    ofs << "f_baryon             = " << f_baryon           << "\n";
    ofs << "f_DM                 = " << f_DM               << "\n";
    ofs << "rho_c                = " << rho_c              << "  # kg/m^3 critical density\n";
    ofs << "delta_contrast       = " << delta_contrast     << "  # dm/rho density perturbation\n";
    ofs << "H0_kms               = " << H0_kms             << "  # km/s/Mpc\n";
    ofs << "H0_si                = " << H0_si              << "  # s^-1\n";
    ofs << "Lambda               = " << Lambda             << "  # m^-2 cosmological constant\n";
    ofs << "Omega_m              = " << Omega_m            << "\n";
    ofs << "Omega_L              = " << Omega_L            << "\n";
    ofs << "z_obs                = " << z_obs              << "\n";
    ofs << "B_cosmic             = " << B_cosmic           << "  # T cosmic magnetic field\n";
    ofs << "B_crit               = " << B_crit             << "  # T Meissner quench\n";
    ofs << "Delta_x              = " << Delta_x            << "  # m quantum scale\n";
    ofs << "A_cmb                = " << A_cmb              << "  # m CMB amplitude\n";
    ofs << "k_cmb                = " << k_cmb              << "  # m^-1 CMB wave vector\n";
    ofs << "omega_cmb            = " << omega_cmb          << "  # rad/s\n";
    ofs << "scale_EM             = " << scale_EM           << "  # macro EM scaling factor\n";
    ofs << "f_TRZ                = " << f_TRZ              << "\n";
    ofs << "f_sc                 = " << f_sc               << "\n";
    ofs << "# --- Cached / derived ---\n";
    ofs << "g_base_cache         = " << g_base_cache       << "  # m/s^2 = GM/r^2 = 3.447e-10\n";
    ofs << "Hz_cache             = " << Hz_cache           << "  # s^-1 H(z=0)\n";
    ofs << "v_exp_cache          = " << v_exp_cache        << "  # m/s H0*r_obs = 9.984e8 = 3.33c\n";
    ofs << "Delta_p_cache        = " << Delta_p_cache      << "  # kg*m/s hbar/Delta_x\n";
    ofs << "a_Lambda_cache       = " << a_Lambda_cache     << "  # m/s^2 PAPER_296 3.30e-36\n";
    ofs << "eta_exp_cache        = " << eta_exp_cache      << "  # PAPER_297 superluminal 3.328\n";
    ofs << "r_H_cache            = " << r_H_cache          << "  # m Hubble sphere 1.322e26\n";
    ofs << "epsilon_GR_cache     = " << epsilon_GR_cache   << "  # PAPER_298 GR curv param 5.056\n";
    ofs << "a_GR_cache           = " << a_GR_cache         << "  # m/s^2 PAPER_298 1.743e-9\n";
    ofs << "r_S_over_r_obs       = " << (2.0 * epsilon_GR_cache / 3.0)
        << "  # Schwarzschild/obs ratio 3.371\n";
    ofs << "Gamma_Lambda         = " << (g_base_cache != 0.0 ? a_Lambda_cache / g_base_cache : 0.0)
        << "  # PAPER_296 9.57e-27\n";
    ofs << "expansion_at_t_H     = " << (1.0 + Hz_cache * T_HUBBLE_S)
        << "  # PAPER_297 1.988 near-doubling\n";
    ofs << "notes_PAPER_296      = Lambda*c2/3=3.30e-36;FIRST UQFF explicit dark-energy;27 orders below g_base\n";
    ofs << "notes_PAPER_297      = eta_exp=3.328>1;FIRST UQFF superluminal;r_H=1.322e26m Hubble boundary\n";
    ofs << "notes_PAPER_298      = epsilon_GR=5.056>1;FIRST UQFF GR>Newton;a_GR=1.743e-9=5*g_base\n";
    for (const auto& p : dynamic_params) {
        ofs << "dyn:" << p.first << " = " << p.second << "\n";
    }
    ofs.close();
}

// ---------------------------------------------------------------------------
// Legacy variable helpers — map old string-key API to named members
// ---------------------------------------------------------------------------

void UniverseDiameterUQFFModule::updateVariable(const std::string& name, double value) {
    dynamic_params[name] = value;
    bool need_cache = false;
    if      (name == "r" || name == "r_obs")           { r_obs          = value; need_cache = true; }
    else if (name == "M" || name == "M_total")         { M_total        = value; need_cache = true; }
    else if (name == "H0" || name == "H0_si")          { H0_si          = value; need_cache = true; }
    else if (name == "Lambda")                         { Lambda         = value; need_cache = true; }
    else if (name == "Omega_m")                        { Omega_m        = value; need_cache = true; }
    else if (name == "Omega_Lambda" || name == "Omega_L") { Omega_L     = value; need_cache = true; }
    else if (name == "z" || name == "z_obs")           { z_obs          = value; need_cache = true; }
    else if (name == "rho_fluid" || name == "rho_c")   { rho_c          = value; }
    else if (name == "delta_rho" || name == "delta_contrast") { delta_contrast = value; }
    else if (name == "B" || name == "B_cosmic")        { B_cosmic       = value; }
    else if (name == "B_crit")                         { B_crit         = value; }
    else if (name == "Delta_x")                        { Delta_x        = value; need_cache = true; }
    else if (name == "A" || name == "A_cmb")           { A_cmb          = value; }
    else if (name == "k" || name == "k_cmb")           { k_cmb          = value; }
    else if (name == "omega" || name == "omega_cmb")   { omega_cmb      = value; }
    else if (name == "scale_macro" || name == "scale_EM") { scale_EM    = value; }
    else if (name == "f_TRZ")                          { f_TRZ          = value; }
    else if (name == "f_sc")                           { f_sc           = value; }
    else if (name == "f_DM")                           { f_DM           = value; }
    else if (name == "f_baryon")                       { f_baryon       = value; }
    if (need_cache) updateCache();
}

void UniverseDiameterUQFFModule::addToVariable(const std::string& name, double delta) {
    auto it = dynamic_params.find(name);
    double current = (it != dynamic_params.end()) ? it->second : 0.0;
    updateVariable(name, current + delta);
}

void UniverseDiameterUQFFModule::subtractFromVariable(const std::string& name,
                                                       double delta) {
    addToVariable(name, -delta);
}

// ---------------------------------------------------------------------------
// Descriptive output
// ---------------------------------------------------------------------------

std::string UniverseDiameterUQFFModule::getEquationText() const {
    return
        "=== WOLFRAM_TERM_UVDIAM_BASE ===\n"
        WOLFRAM_TERM_UVDIAM_BASE "\n\n"
        "=== WOLFRAM_TERM_UVDIAM_LAMBDA ===\n"
        WOLFRAM_TERM_UVDIAM_LAMBDA "\n\n"
        "=== WOLFRAM_TERM_UVDIAM_HUBBLE ===\n"
        WOLFRAM_TERM_UVDIAM_HUBBLE "\n\n"
        "=== WOLFRAM_TERM_UVDIAM_GR_CURV ===\n"
        WOLFRAM_TERM_UVDIAM_GR_CURV "\n\n"
        "--- 9-Term Architecture ---\n"
        "  1. a_base    = G*M/r^2 * (1+H(z)*t)   Newtonian + Hubble expansion coupling\n"
        "  2. a_Lambda  = Lambda*c^2/3             [PAPER_296] dark-energy vacuum accel\n"
        "  3. a_quantum = hbar/sqrt(Dx*Dp)*2pi/tH Heisenberg cosmic uncertainty\n"
        "  4. a_EM      = q*v_exp*B/m_p*(1+eta)   [PAPER_297] superluminal Lorentz\n"
        "  5. a_fluid   = rho_c*G*delta*r_obs      Jeans density perturbation drive\n"
        "  6. a_osc     = 2A cos(kx)cos(wt)+...   CMB standing+traveling (PAPER_288)\n"
        "  7. a_DM_pert = g_base*f_DM*delta        DM fractional density perturbation\n"
        "  8. a_GR      = g_base*epsilon_GR        [PAPER_298] GR curvature dominance\n"
        "  g_UVDIAM = (sum all 8) * (1-B/B_crit) * (1+f_TRZ)\n\n"
        "--- PAPER_296 Cosmological Constant Direct Vacuum Acceleration ---\n"
        "  a_Lambda = Lambda*c^2/3 = 3.30e-36 m/s^2  (FIRST UQFF explicit Lambda term)\n"
        "  Gamma_Lambda = a_Lambda/g_base = 9.57e-27\n"
        "  d_Lambda = 0.5*a_Lambda*t_H^2 = 0.313 m  (macroscopic cosmic displacement)\n\n"
        "--- PAPER_297 Superluminal Hubble Expansion Ratio ---\n"
        "  eta_exp = v_exp/c = H0*r_obs/c = 3.328 > 1  (FIRST UQFF eta_exp > 1)\n"
        "  r_H = c/H0 = 1.322e26 m  r_obs = 3.328 * r_H (3.328 Hubble lengths)\n"
        "  Hubble expansion factor at t_H = 1.988 (near-doubling of base gravity)\n\n"
        "--- PAPER_298 GR Curvature Dominance ---\n"
        "  epsilon_GR = 3GM/(r*c^2) = 5.056 > 1  (FIRST UQFF epsilon_GR > 1)\n"
        "  a_GR = g_base * epsilon_GR = 1.743e-9 m/s^2 (5x Newtonian base!)\n"
        "  r_S/r_obs = 2*epsilon_GR/3 = 3.371  (universe at ~30% of Schwarzschild)";
}

void UniverseDiameterUQFFModule::printVariables() const {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "UniverseDiameterUQFFModule — UQFF 2.0 (Session 84, 26th C++ module)\n";
    std::cout << "  [Universe]  r_obs="     << r_obs
              << "  M_total="    << M_total
              << "  f_DM="       << f_DM     << "\n";
    std::cout << "  [Cosmo]     H0_si="     << H0_si
              << "  Lambda="     << Lambda
              << "  z_obs="      << z_obs    << "\n";
    std::cout << "  [EM/Quant]  B_cosmic="  << B_cosmic
              << "  Delta_x="    << Delta_x
              << "  A_cmb="      << A_cmb    << "\n";
    std::cout << "  [UQFF]      f_TRZ="     << f_TRZ
              << "  f_sc="       << f_sc
              << "  B_crit="     << B_crit   << "\n";
    std::cout << "  [Cache]     g_base="    << g_base_cache
              << "  Hz="         << Hz_cache
              << "  v_exp="      << v_exp_cache << " m/s\n";
    std::cout << "  [PAPER_296] a_Lambda="  << a_Lambda_cache
              << "  Gamma_Lambda=" << (g_base_cache != 0.0 ? a_Lambda_cache / g_base_cache : 0.0)
              << "  d_Lambda=" << (0.5 * a_Lambda_cache * T_HUBBLE_S * T_HUBBLE_S) << " m\n";
    std::cout << "  [PAPER_297] eta_exp="   << eta_exp_cache
              << "  r_H="        << r_H_cache
              << "  r_obs/r_H="  << (r_obs / r_H_cache) << "\n";
    std::cout << "  [PAPER_298] epsilon_GR=" << epsilon_GR_cache
              << "  a_GR="       << a_GR_cache
              << "  r_S/r_obs="  << (2.0 * epsilon_GR_cache / 3.0) << "\n";
    if (!dynamic_params.empty()) {
        for (const auto& p : dynamic_params) {
            std::cout << "  [dyn]       " << p.first << "=" << p.second << "\n";
        }
    }
}

#endif // UNIVERSE_DIAMETER_UQFF_MODULE_H
