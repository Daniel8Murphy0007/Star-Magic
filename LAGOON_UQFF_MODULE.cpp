#ifndef LAGOON_UQFF_MODULE_H
#define LAGOON_UQFF_MODULE_H
// =============================================================================
// LAGOON_UQFF_MODULE.cpp â€” Lagoon Nebula (M8 / NGC 6523) UQFF 2.0
// Session 87 | 29th C++ UQFF module | FIRST H II Region UQFF module
// Papers: PAPER_305 (SFR Mass Runaway Amplifier)
//         PAPER_306 (Herschel 36 Radiation Erosion â€” eta_rad = 1.53e18)
//         PAPER_307 (Dual Radiation-EM Barrier â€” a_EM/a_rad = 12.77)
// 9-term pipeline:
//   g_base(t) Ã— (1+HzÂ·t) Ã— (1-B/B_crit) Ã— (1+f_TRZ)
//   + Ug_sum + Î›/3 + â„/(m_HÂ·Î”xÂ²) + a_EM_turb + g_fluid
//   + A_oscÂ·cos(kÂ·r+Ï‰Â·t) + GÂ·M_DM/rÂ² âˆ’ P_rad
// =============================================================================

#include <cmath>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

class LagoonUQFFModule {
public:
    // ---- Physical constants ------------------------------------------------
    static constexpr double G_NEWTON   = 6.6743e-11;   // mÂ³/(kgÂ·sÂ²)
    static constexpr double C_LIGHT    = 2.998e8;       // m/s
    static constexpr double HBAR       = 1.0546e-34;    // JÂ·s
    static constexpr double LAMBDA_CC  = 1.089e-52;     // mâ»Â²  (cosmological constant Î›)
    static constexpr double Q_PROTON   = 1.602e-19;     // C
    static constexpr double M_H        = 1.6726e-27;    // kg  (hydrogen atom mass)
    static constexpr double M_SUN      = 1.989e30;      // kg
    static constexpr double PI         = 3.14159265358979323846;
    static constexpr double YR_TO_S    = 3.15576e7;     // s/yr
    static constexpr double T_HUBBLE_S = 4.355e17;      // s  (13.8 Gyr)
    static constexpr double E_VAC_DEF  = 7.09e-36;      // J/mÂ³  (vacuum energy density)
    static constexpr double MPC_TO_M   = 3.0857e22;     // m/Mpc

    // ---- Constructor -------------------------------------------------------
    LagoonUQFFModule()
        : M0(1.0e4 * M_SUN),              // 1e4 solar masses â€” molecular cloud mass
          SFR_kg_s(0.1 * M_SUN / YR_TO_S),// 0.1 M_sun/yr â†’ 6.303e21 kg/s
          r(5.2e17),                       // m  (~55 ly half-span)
          z(0.0013),                       // redshift  (~1.25 kpc)
          H0(2.268e-18),                   // sâ»Â¹  (70.0 km/s/Mpc)
          Omega_m(0.3),
          Omega_L(0.7),
          rho_fluid(1.0e-20),              // kg/mÂ³  (nebula gas density)
          B(1.0e-5),                       // T  (nebula B-field)
          B_crit(4.4e13),                  // T  (magnetar reference)
          L_H36(7.65e31),                  // W  (Herschel 36 O7V luminosity)
          v_gas(1.0e5),                    // m/s  (turbulent velocity)
          f_TRZ(0.01),                     // TRZ resonance correction
          f_sc(1.0),                       // SC amplitude factor
          E_vac(E_VAC_DEF),               // J/mÂ³
          A_osc(1.0e-13),                 // resonant oscillation amplitude (m/sÂ²)
          k_osc(1.0e-17),                 // wave number (mâ»Â¹)
          omega_osc(2.0e-15),             // angular frequency (rad/s)
          delta_rho(1.0e-22),             // DM density perturbation (kg/mÂ³)
          M_DM(5.0e4 * M_SUN),           // dark matter halo mass (kg)
          Delta_x(1.0e15),               // HUP position uncertainty (m)
          V_fluid(1.0e3),                // fluid velocity Navier-Stokes proxy (m/s)
          logging_enabled(false)
    {
        updateCache();
    }

    // ---- Public setters ----------------------------------------------------
    void setEnableLogging(bool enable) { logging_enabled = enable; }

    void setDynamicParameter(const std::string& key, double value) {
        dynamic_params[key] = value;
        updateCache();
    }

    double getDynamicParameter(const std::string& key) const {
        auto it = dynamic_params.find(key);
        if (it != dynamic_params.end()) return it->second;
        throw std::runtime_error("LagoonUQFFModule: parameter not found: " + key);
    }

    void updateVariable(const std::string& key, double value) {
        if      (key == "M0")        { M0        = value; }
        else if (key == "SFR_kg_s")  { SFR_kg_s  = value; }
        else if (key == "r")         { r          = value; }
        else if (key == "z")         { z          = value; }
        else if (key == "B")         { B          = value; }
        else if (key == "L_H36")     { L_H36      = value; }
        else if (key == "v_gas")     { v_gas      = value; }
        else if (key == "rho_fluid") { rho_fluid  = value; }
        else if (key == "f_TRZ")     { f_TRZ      = value; }
        else                         { dynamic_params[key] = value; }
        updateCache();
    }

    void addToVariable(const std::string& key, double delta) {
        updateVariable(key, getVariableValue(key) + delta);
    }

    void subtractFromVariable(const std::string& key, double delta) {
        updateVariable(key, getVariableValue(key) - delta);
    }

    // ---- Main compute ------------------------------------------------------
    // g_Lagoon(t) â€” full 9-term UQFF pipeline
    double computeG(double t) const {
        // M(t) = M0 + SFR_kg_s * t  [PAPER_305 mass growth]
        double M_t       = M0 + SFR_kg_s * t;
        double g_base    = G_NEWTON * M_t / (r * r);

        // Hubble expansion correction: (1 + HzÂ·t)
        double Hz        = Hz_cache;
        double hz_corr   = 1.0 + Hz * t;

        // Superconducting B correction: (1 âˆ’ B/B_crit)
        double sc_corr   = 1.0 - B / B_crit;

        // TRZ resonance enhancement: (1 + f_TRZ)
        double trz_corr  = 1.0 + f_TRZ;

        // Stage 1 â€” base pipeline
        double g_pipeline = g_base * hz_corr * sc_corr * trz_corr;

        // Stage 2 â€” additive terms
        double ug        = computeUgSum();
        double lambda    = (LAMBDA_CC * C_LIGHT * C_LIGHT) / 3.0;  // Î›cÂ²/3
        double quantum   = computeQuantumTerm();
        double em_turb   = a_EM_cache;                              // [PAPER_307]
        double fluid     = computeFluidTerm(g_base);
        double resonant  = computeResonantTerm(t);
        double dm        = computeDMTerm();

        // Stage 3 â€” P_rad subtracted [PAPER_306]
        double p_rad     = a_rad_cache;

        double g_total   = g_pipeline + ug + lambda + quantum + em_turb
                          + fluid + resonant + dm - p_rad;

        if (logging_enabled) {
            std::cout << "[LagoonUQFFModule] t=" << t
                      << "  M_t=" << M_t
                      << "  g_pipeline=" << g_pipeline
                      << "  em_turb=" << em_turb
                      << "  p_rad=" << p_rad
                      << "  g_total=" << g_total << "\n";
        }
        return g_total;
    }

    // ---- Unique-physics cached accessors -----------------------------------
    // PAPER_305 â€” SFR Mass Runaway Amplifier
    double getSFRDeltaMoverM0_1Myr()   const { return msf_1Myr_cache; }      // 10.0
    double getMFactor_1Myr()           const { return m_factor_1Myr_cache; } // 11.0
    double getTConsume_yr()            const { return t_consume_cache; }      // 1e5 yr
    double getSFRoverM0_yr()           const { return SFR_over_M0_cache; }    // 1e-5 yrâ»Â¹
    // PAPER_306 â€” Herschel 36 Radiation Erosion
    double getARad_m_s2()              const { return a_rad_cache; }          // 7.51e6
    double getEtaRad()                 const { return eta_rad_cache; }        // 1.53e18
    // PAPER_307 â€” Dual Radiation-EM Barrier
    double getAEM_m_s2()               const { return a_EM_cache; }           // 9.59e7
    double getEtaEM()                  const { return eta_EM_cache; }         // 1.96e19
    double getAEMoverARad()            const { return a_EM_over_a_rad_cache; }// 12.77

    // ---- Text output -------------------------------------------------------
    std::string getEquationText() const {
        std::ostringstream oss;
        oss << "=== Lagoon Nebula (M8/NGC 6523) UQFF 2.0 â€” Session 87 ===\n"
            << "System: H II region at z=0.0013 (~1.25 kpc), ionised by Herschel 36 (O7V)\n"
            << "Master equation:\n"
            << "  g_Lagoon(t) = [GÂ·M(t)/rÂ²]Â·(1+HzÂ·t)Â·(1-B/B_crit)Â·(1+f_TRZ)\n"
            << "              + Ug_sum + Î›cÂ²/3 + â„/(m_HÂ·Î”xÂ²) + a_EM_turb\n"
            << "              + g_fluid + A_oscÂ·cos(kÂ·r+Ï‰Â·t) + GÂ·M_DM/rÂ² âˆ’ a_rad\n\n"
            << "PAPER_305 â€” SFR Mass Runaway Amplifier:\n"
            << "  M(t) = M0 + SFR_kg_sÂ·t\n"
            << "  Î”M/M0 at 1 Myr = " << msf_1Myr_cache << "  (= 10.0)\n"
            << "  m_factor(1 Myr) = " << m_factor_1Myr_cache << "  (= 11.0)\n"
            << "  t_consume       = " << t_consume_cache << " yr  (= 100 kyr)\n"
            << "  SFR/M0          = " << SFR_over_M0_cache << " yrâ»Â¹  (= 1e-5 yrâ»Â¹)\n\n"
            << "PAPER_306 â€” Herschel 36 Radiation Erosion:\n"
            << "  flux = L_H36 / (4Ï€Â·rÂ²Â·c) = " << L_H36/(4.0*PI*r*r*C_LIGHT) << " Pa\n"
            << "  a_rad = flux / Ï_fluid   = " << a_rad_cache << " m/sÂ²\n"
            << "  Î·_rad = a_rad / g_base   = " << eta_rad_cache << "\n\n"
            << "PAPER_307 â€” Dual Radiation-EM Barrier:\n"
            << "  a_EM = qÂ·v_gasÂ·B / m_H  = " << a_EM_cache << " m/sÂ²\n"
            << "  Î·_EM = a_EM / g_base    = " << eta_EM_cache << "\n"
            << "  a_EM / a_rad            = " << a_EM_over_a_rad_cache << "  (EM leads radiation by 12.77Ã—)\n";
        return oss.str();
    }

    void printParameters() const {
        std::cout << getEquationText();
        std::cout << "Parameters:\n"
                  << "  M0        = " << M0        << " kg  ("  << M0/M_SUN  << " M_sun)\n"
                  << "  SFR_kg_s  = " << SFR_kg_s  << " kg/s\n"
                  << "  r         = " << r          << " m\n"
                  << "  z         = " << z          << "\n"
                  << "  B         = " << B          << " T\n"
                  << "  L_H36     = " << L_H36      << " W\n"
                  << "  v_gas     = " << v_gas      << " m/s\n"
                  << "  rho_fluid = " << rho_fluid  << " kg/mÂ³\n"
                  << "  g_base(0) = " << g_base_cache << " m/sÂ²\n";
    }

    void exportState(const std::string& filename) const {
        std::ofstream f(filename);
        if (!f) {
            std::cerr << "[LagoonUQFFModule] Cannot open file: " << filename << "\n";
            return;
        }
        f << "# LagoonUQFFModule state export â€” Session 87 â€” UQFF 2.0\n"
          << "M0="           << M0           << "\n"
          << "SFR_kg_s="     << SFR_kg_s     << "\n"
          << "r="            << r            << "\n"
          << "z="            << z            << "\n"
          << "H0="           << H0           << "\n"
          << "Omega_m="      << Omega_m      << "\n"
          << "Omega_L="      << Omega_L      << "\n"
          << "rho_fluid="    << rho_fluid    << "\n"
          << "B="            << B            << "\n"
          << "B_crit="       << B_crit       << "\n"
          << "L_H36="        << L_H36        << "\n"
          << "v_gas="        << v_gas        << "\n"
          << "f_TRZ="        << f_TRZ        << "\n"
          << "E_vac="        << E_vac        << "\n"
          << "g_base_t0="    << g_base_cache << "\n"
          << "Hz="           << Hz_cache     << "\n"
          // PAPER_305
          << "msf_1Myr="        << msf_1Myr_cache        << "\n"
          << "m_factor_1Myr="   << m_factor_1Myr_cache   << "\n"
          << "t_consume_yr="    << t_consume_cache        << "\n"
          << "SFR_over_M0_yr="  << SFR_over_M0_cache     << "\n"
          // PAPER_306
          << "a_rad="           << a_rad_cache            << "\n"
          << "eta_rad="         << eta_rad_cache          << "\n"
          // PAPER_307
          << "a_EM="            << a_EM_cache             << "\n"
          << "eta_EM="          << eta_EM_cache           << "\n"
          << "a_EM_over_a_rad=" << a_EM_over_a_rad_cache  << "\n"
          // WOLFRAM_TERM annotations
          << "# WOLFRAM_TERM: LAGOON_BASE  g_base=" << g_base_cache
          << "  M=" << M0 << "  r=" << r << "\n"
          << "# WOLFRAM_TERM: LAGOON_SFR [P305]  msf_1Myr=" << msf_1Myr_cache
          << "  m_factor=" << m_factor_1Myr_cache
          << "  t_consume=" << t_consume_cache << " yr\n"
          << "# WOLFRAM_TERM: LAGOON_RAD [P306]  a_rad=" << a_rad_cache
          << "  eta_rad=" << eta_rad_cache << "\n"
          << "# WOLFRAM_TERM: LAGOON_EM_TURB [P307]  a_EM=" << a_EM_cache
          << "  eta_EM=" << eta_EM_cache
          << "  a_EM_over_a_rad=" << a_EM_over_a_rad_cache << "\n";
        f.close();
        if (logging_enabled)
            std::cout << "[LagoonUQFFModule] State exported to " << filename << "\n";
    }

    template<typename OtherModule>
    double cross_validate(OtherModule& other, double t) {
        double g_this  = this->computeG(t);
        double g_other = other.computeG(t);
        return (g_this - g_other) / g_this;  // fractional residual
    }

private:
    // ---- Named typed members -----------------------------------------------
    double M0, SFR_kg_s, r, z, H0, Omega_m, Omega_L;
    double rho_fluid, B, B_crit, L_H36, v_gas;
    double f_TRZ, f_sc, E_vac;
    double A_osc, k_osc, omega_osc;
    double delta_rho, M_DM, Delta_x, V_fluid;
    bool   logging_enabled;
    std::map<std::string, double> dynamic_params;

    // ---- Cache members -----------------------------------------------------
    double g_base_cache          = 0.0;
    double Hz_cache              = 0.0;
    // PAPER_305
    double msf_1Myr_cache        = 0.0;
    double m_factor_1Myr_cache   = 0.0;
    double t_consume_cache       = 0.0;
    double SFR_over_M0_cache     = 0.0;
    // PAPER_306
    double a_rad_cache           = 0.0;
    double eta_rad_cache         = 0.0;
    // PAPER_307
    double a_EM_cache            = 0.0;
    double eta_EM_cache          = 0.0;
    double a_EM_over_a_rad_cache = 0.0;

    // ---- Cache update ------------------------------------------------------
    void updateCache() {
        // Base gravity at t=0
        g_base_cache = G_NEWTON * M0 / (r * r);

        // Hubble parameter at z
        Hz_cache = computeHz();

        // PAPER_305 â€” SFR Mass Runaway Amplifier
        double M0_sun       = M0 / M_SUN;                       // solar masses
        double SFR_sun_yr   = SFR_kg_s * YR_TO_S / M_SUN;      // M_sun/yr
        msf_1Myr_cache      = SFR_sun_yr * 1.0e6 / M0_sun;     // Î”M/M0 at 1 Myr = 10.0
        m_factor_1Myr_cache = 1.0 + msf_1Myr_cache;            // 11.0
        t_consume_cache     = M0_sun / SFR_sun_yr;              // 1e5 yr (100 kyr)
        SFR_over_M0_cache   = SFR_sun_yr / M0_sun;              // 1e-5 yrâ»Â¹

        // PAPER_306 â€” Herschel 36 Radiation Erosion
        double flux  = L_H36 / (4.0 * PI * r * r * C_LIGHT);   // Pa
        a_rad_cache  = flux / rho_fluid;                        // 7.51e6 m/sÂ²
        eta_rad_cache = a_rad_cache / g_base_cache;             // 1.53e18

        // PAPER_307 â€” Dual Radiation-EM Barrier
        a_EM_cache           = Q_PROTON * v_gas * B / M_H;      // 9.59e7 m/sÂ²
        eta_EM_cache         = a_EM_cache / g_base_cache;       // 1.96e19
        a_EM_over_a_rad_cache = a_EM_cache / a_rad_cache;       // 12.77
    }

    // ---- Private helpers ---------------------------------------------------
    double computeHz() const {
        double factor = Omega_m * std::pow(1.0 + z, 3.0) + Omega_L;
        return H0 * std::sqrt(factor);   // H(z) = H0*sqrt(Î©_mÂ·(1+z)Â³ + Î©_Î›)
    }

    double computeUgSum() const {
        // Ug1: magnetic pressure / fluid density
        double mu0 = 4.0e-7 * PI;
        double Ug1 = (B * B) / (2.0 * mu0 * rho_fluid);
        // Ug2: charge-velocity / (mass Â· radius)
        double Ug2 = Q_PROTON * v_gas / (M_H * r);
        // Ug3: vacuum energy / (mass Â· rÂ²)
        double Ug3 = E_vac / (M_H * r * r);
        // Ug4: vacuum energy Â· r / mass
        double Ug4 = E_vac * r / M_H;
        return Ug1 + Ug2 + Ug3 + Ug4;
    }

    double computeQuantumTerm() const {
        // HUP: a_HUP = â„ / (m_H Â· Î”xÂ²)
        return HBAR / (M_H * Delta_x * Delta_x);
    }

    double computeFluidTerm(double g_base) const {
        // Navier-Stokes proxy: Î½Â·(v_gas/r) / (Ï_cloud Â· rÂ²)
        double nu         = 1.0e10;                             // mÂ²/s turbulent viscosity
        double grad_v     = v_gas / r;
        double rho_cloud  = M0 / ((4.0 / 3.0) * PI * r * r * r);
        return nu * grad_v / (rho_cloud * r * r);
    }

    double computeResonantTerm(double t) const {
        return A_osc * std::cos(k_osc * r + omega_osc * t);
    }

    double computeDMTerm() const {
        return G_NEWTON * M_DM / (r * r);
    }

    double getVariableValue(const std::string& key) const {
        if (key == "M0")        return M0;
        if (key == "SFR_kg_s")  return SFR_kg_s;
        if (key == "r")         return r;
        if (key == "z")         return z;
        if (key == "B")         return B;
        if (key == "L_H36")     return L_H36;
        if (key == "v_gas")     return v_gas;
        if (key == "rho_fluid") return rho_fluid;
        if (key == "f_TRZ")     return f_TRZ;
        auto it = dynamic_params.find(key);
        if (it != dynamic_params.end()) return it->second;
        throw std::runtime_error("LagoonUQFFModule: unknown variable: " + key);
    }
};

// =============================================================================
// WOLFRAM_TERM: LAGOON_BASE
//   g_base(t=0) = G * M0 / r^2
//              = 6.6743e-11 * 1.989e34 / (5.2e17)^2 = 4.91e-12 m/s^2
//   M0 = 1.989e34 kg (1e4 M_sun);  r = 5.2e17 m;  z = 0.0013
// =============================================================================

// =============================================================================
// WOLFRAM_TERM: LAGOON_SFR  [PAPER_305 â€” SFR Mass Runaway Amplifier]
//   msf(1 Myr) = SFR_sun * 1e6 yr / M0_sun = 0.1 * 1e6 / 1e4 = 10.0
//   m_factor(1 Myr) = 11.0  â†’  g(1 Myr) = 11 Ã— g(0)
//   t_consume = M0/SFR = 1e5 yr = 100 kyr
//   SFR/M0 = 1e-5 yr^-1  [FIRST UQFF SFR runaway: Î”M > M0 at 1 Myr]
// =============================================================================

// =============================================================================
// WOLFRAM_TERM: LAGOON_RAD  [PAPER_306 â€” Herschel 36 Radiation Erosion]
//   flux = L_H36 / (4*pi*r^2*c) = 7.65e31/(4*pi*(5.2e17)^2*3e8) = 7.511e-14 Pa
//   a_rad = flux / rho_fluid = 7.511e-14 / 1e-20 = 7.51e6 m/s^2
//   eta_rad = a_rad / g_base = 7.51e6 / 4.91e-12 = 1.53e18
//   [FIRST UQFF single-star (Herschel 36 O7V) radiation pressure parameter]
// =============================================================================

// =============================================================================
// WOLFRAM_TERM: LAGOON_EM_TURB  [PAPER_307 â€” Dual Radiation-EM Barrier]
//   a_EM = q * v_gas * B / m_H = 1.602e-19 * 1e5 * 1e-5 / 1.6726e-27 = 9.59e7 m/s^2
//   eta_EM = a_EM / g_base = 9.59e7 / 4.91e-12 = 1.96e19
//   a_EM / a_rad = 9.59e7 / 7.51e6 = 12.77
//   [FIRST UQFF dual-barrier: a_EM AND a_rad both >> g_base simultaneously]
//   [EM leads radiation by 12.77x â€” explains Lagoon extended HII morphology]
// =============================================================================

#endif // LAGOON_UQFF_MODULE_H
