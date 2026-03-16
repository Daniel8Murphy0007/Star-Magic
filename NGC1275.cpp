/**
 * ================================================================================================
 * Header: NGC1275.h
 * 
 * Description: C++ Module for NGC 1275 (Magnetic Monster Perseus A) Class
 *              This is the sixteenth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on active galactic nucleus evolution
 *              and gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025,
 *              upgraded to full UQFF 2.0 + 3-tier buoyancy pipeline (10→13-term MUGE) on March 16, 2026).
 * 
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for NGC 1275 evolution.
 *          Includes ALL terms: base gravity (static M), cosmic expansion (H(z)), magnetic field B(t),
 *          filament support F(t), black hole influence, UQFF Ug components with f_TRZ, Lambda,
 *          quantum uncertainty, scaled EM with [UA], fluid dynamics, oscillatory waves,
 *          DM/density perturbations, and cooling flow term. Supports dynamic variable updates.
 * 
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: NGC1275 ngc1275;
 *              Compute: double g = ngc1275.compute_g_NGC1275(t);
 * 
 * Key Features:
 *   - Default values from UQFF document: M = 1e11 Msun, r = 1.893e21 m (200k ly), M_BH = 8e8 Msun,
 *     z = 0.0176, Hz ≈ 2.20e-18 s^-1, B0 = 5e-9 T, tau_B = 100 Myr, F0 = 0.1, tau_fil = 100 Myr,
 *     rho_cool = 1e-20 kg/m^3, v_cool = 3e3 m/s.
 *   - Units handled: Msun to kg, ly to m; cooling term as (rho * v_cool^2) / rho_fluid for acceleration.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_NGC1275(r, t) with every term explicitly included.
 * 
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: March 16, 2026
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef NGC_1275_H
#define NGC_1275_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <limits>
#include <map>
#include <fstream>

class NGC1275 {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M;               // Total galaxy mass (kg)
    double r;               // Radius (m)
    double Hz;              // Hubble parameter at z (s^-1)
    double B0;              // Initial magnetic field (T)
    double tau_B;           // B decay timescale (s)
    double B_crit;          // Critical B field (T)
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double gas_v;           // Gas velocity for EM (m/s)
    double f_TRZ;           // Time-reversal factor
    double M_BH;            // Black hole mass (kg)
    double r_BH;            // Black hole radius (m)
    double F0;              // Initial filament factor
    double tau_fil;         // Filament timescale (s)
    double rho_cool;        // Cooling flow density (kg/m^3)
    double v_cool;          // Cooling flow velocity (m/s)
    double rho_fluid;       // Fluid density (kg/m^3)
    double rho_vac_UA;      // UA vacuum density (J/m^3)
    double rho_vac_SCm;     // SCm vacuum density (J/m^3)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass for EM acceleration
    double z_gal;           // Galaxy redshift

    // Additional parameters for full inclusion of terms
    double hbar;            // Reduced Planck's constant
    double t_Hubble;        // Hubble time (s)
    double delta_x;         // Position uncertainty (m)
    double delta_p;         // Momentum uncertainty (kg m/s)
    double integral_psi;    // Wavefunction integral approximation
    double A_osc;           // Oscillatory amplitude (m/s^2)
    double k_osc;           // Wave number (1/m)
    double omega_osc;       // Angular frequency (rad/s)
    double x_pos;           // Position for oscillation (m)
    double t_Hubble_gyr;    // Hubble time in Gyr
    double M_DM_factor;     // Dark matter mass fraction
    double delta_rho_over_rho; // Density perturbation fraction

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 = G*M/r^2
    double g_BH;            // Cached BH acceleration

    // UQFF 2.0: 3-tier buoyancy parameters
    double beta_i;          // Buoyancy coupling coefficient (=0.61)
    double omega_g;         // Gravitational angular frequency (=7.3e-16)
    double U_UA;            // UA vacuum coupling constant (=1e-11)
    double M_ext_vc;        // Virgo Cluster mass kg [outer frame ~1.2e15 M_sun]
    double r_ext_vc;        // Distance NGC 1275 (Perseus) to Virgo Cluster m [~77 Mpc]
    bool logging_enabled;   // UQFF 2.0 diagnostic logging flag
    std::map<std::string, double> dynamic_params; // UQFF 2.0 runtime parameters

public:
    // Constructor with default UQFF values
    NGC1275() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~NGC1275() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M = 1e11 * M_sun;
        double ly_to_m = 9.461e15;
        r = 200000.0 * ly_to_m;
        z_gal = 0.0176;
        double Hz_kms = 70 * sqrt(0.3 * pow(1 + z_gal, 3) + 0.7);  // km/s/Mpc
        Hz = (Hz_kms * 1000 / 3.086e19);  // s^-1
        B0 = 5e-9;
        tau_B = 100e6 * 3.156e7;
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        gas_v = 1e5;
        f_TRZ = 0.1;
        M_BH = 8e8 * M_sun;
        r_BH = 1e18;  // Approximate influence radius
        F0 = 0.1;
        tau_fil = 100e6 * 3.156e7;
        rho_cool = 1e-20;
        v_cool = 3e3;
        rho_fluid = 1e-20;
        rho_vac_UA = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        scale_EM = 1e-12;
        proton_mass = 1.673e-27;

        // Full terms defaults
        hbar = 1.0546e-34;
        t_Hubble = 13.8e9 * 3.156e7;
        t_Hubble_gyr = 13.8;
        delta_x = 1e-10;
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
        A_osc = 1e-10;
        k_osc = 1.0 / r;
        omega_osc = 2 * M_PI / (r / c_light);
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;

        // UQFF 2.0: 3-tier buoyancy + logging init
        beta_i = 0.61;
        omega_g = 7.3e-16;
        U_UA = 1e-11;
        M_ext_vc = 2.387e45;    // Virgo Cluster ~1.2e15 M_sun
        r_ext_vc = 2.38e24;     // ~77 Mpc NGC 1275 (Perseus cluster) to Virgo Cluster
        logging_enabled = false;

        updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        ug1_base = (G * M) / (r * r);
        g_BH = (G * M_BH) / (r_BH * r_BH);
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M") { M = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "Hz") { Hz = newValue; }
        else if (varName == "B0") { B0 = newValue; }
        else if (varName == "tau_B") { tau_B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "gas_v") { gas_v = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "M_BH") { M_BH = newValue; }
        else if (varName == "r_BH") { r_BH = newValue; }
        else if (varName == "F0") { F0 = newValue; }
        else if (varName == "tau_fil") { tau_fil = newValue; }
        else if (varName == "rho_cool") { rho_cool = newValue; }
        else if (varName == "v_cool") { v_cool = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
        else if (varName == "proton_mass") { proton_mass = newValue; }
        else if (varName == "z_gal") { z_gal = newValue; }
        // Full terms
        else if (varName == "hbar") { hbar = newValue; }
        else if (varName == "t_Hubble") { t_Hubble = newValue; }
        else if (varName == "t_Hubble_gyr") { t_Hubble_gyr = newValue; }
        else if (varName == "delta_x") { delta_x = newValue; }
        else if (varName == "delta_p") { delta_p = newValue; }
        else if (varName == "integral_psi") { integral_psi = newValue; }
        else if (varName == "A_osc") { A_osc = newValue; }
        else if (varName == "k_osc") { k_osc = newValue; }
        else if (varName == "omega_osc") { omega_osc = newValue; }
        else if (varName == "x_pos") { x_pos = newValue; }
        else if (varName == "M_DM_factor") { M_DM_factor = newValue; }
        else if (varName == "delta_rho_over_rho") { delta_rho_over_rho = newValue; }
        // UQFF 2.0 buoyancy parameters
        else if (varName == "beta_i") { beta_i = newValue; }
        else if (varName == "omega_g") { omega_g = newValue; }
        else if (varName == "U_UA") { U_UA = newValue; }
        else if (varName == "M_ext_vc") { M_ext_vc = newValue; }
        else if (varName == "r_ext_vc") { r_ext_vc = newValue; }
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return false;
        }
        updateCache();
        return true;
    }

    // Addition method for variables
    bool addToVariable(const std::string& varName, double delta) {
        return setVariable(varName, getVariable(varName) + delta);
    }

    // Subtraction method for variables
    bool subtractFromVariable(const std::string& varName, double delta) {
        return addToVariable(varName, -delta);
    }

    // Getter for any variable (helper for add/subtract)
    double getVariable(const std::string& varName) const {
        if (varName == "G") return G;
        else if (varName == "M") return M;
        else if (varName == "r") return r;
        else if (varName == "Hz") return Hz;
        else if (varName == "B0") return B0;
        else if (varName == "tau_B") return tau_B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "gas_v") return gas_v;
        else if (varName == "f_TRZ") return f_TRZ;
        else if (varName == "M_BH") return M_BH;
        else if (varName == "r_BH") return r_BH;
        else if (varName == "F0") return F0;
        else if (varName == "tau_fil") return tau_fil;
        else if (varName == "rho_cool") return rho_cool;
        else if (varName == "v_cool") return v_cool;
        else if (varName == "rho_fluid") return rho_fluid;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "scale_EM") return scale_EM;
        else if (varName == "proton_mass") return proton_mass;
        else if (varName == "z_gal") return z_gal;
        // Full terms
        else if (varName == "hbar") return hbar;
        else if (varName == "t_Hubble") return t_Hubble;
        else if (varName == "t_Hubble_gyr") return t_Hubble_gyr;
        else if (varName == "delta_x") return delta_x;
        else if (varName == "delta_p") return delta_p;
        else if (varName == "integral_psi") return integral_psi;
        else if (varName == "A_osc") return A_osc;
        else if (varName == "k_osc") return k_osc;
        else if (varName == "omega_osc") return omega_osc;
        else if (varName == "x_pos") return x_pos;
        else if (varName == "M_DM_factor") return M_DM_factor;
        else if (varName == "delta_rho_over_rho") return delta_rho_over_rho;
        // UQFF 2.0 buoyancy parameters
        else if (varName == "beta_i") return beta_i;
        else if (varName == "omega_g") return omega_g;
        else if (varName == "U_UA") return U_UA;
        else if (varName == "M_ext_vc") return M_ext_vc;
        else if (varName == "r_ext_vc") return r_ext_vc;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return 0.0;
        }
    }

    // B(t) computation
    double B_t(double t) const {
        return B0 * exp(-t / tau_B);
    }

    // F(t) computation
    double F_t(double t) const {
        return F0 * exp(-t / tau_fil);
    }

    // Ug terms computation
    double compute_Ug(double Bt, double Ft) const {
        double Ug1 = ug1_base;
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - Bt / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ) * (1 + Ft);
    }

    // Volume computation for fluid
    double compute_V() const {
        return (4.0 / 3.0) * M_PI * r * r * r;
    }

    // Main MUGE computation (includes ALL terms)
    double compute_g_NGC1275(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Bt = B_t(t);
        double Ft = F_t(t);

        // Term 1: Base + Hz + B + F corrections
        double corr_H = 1 + Hz * t;
        double corr_B = 1 - Bt / B_crit;
        double corr_F = 1 + Ft;
        double term1 = ug1_base * corr_H * corr_B * corr_F;

        // BH term
        double term_BH = g_BH;

        // Term 2: UQFF Ug with f_TRZ, B, F
        double term2 = compute_Ug(Bt, Ft);

        // Term 3: Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // Term 4: Scaled EM with UA
        double cross_vB = gas_v * Bt;  // Magnitude, assuming perpendicular
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
        double term4 = (em_base * corr_UA) * scale_EM;

        // Quantum uncertainty term
        double sqrt_unc = sqrt(delta_x * delta_p);
        double term_q = (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);

        // Fluid term (effective acceleration)
        double V = compute_V();
        double term_fluid = (rho_fluid * V * ug1_base) / M;

        // Oscillatory terms (real parts)
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        double term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term (converted to acceleration)
        double M_dm = M * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M / (r * r * r);
        double term_dm_force_like = (M + M_dm) * (pert1 + pert2);
        double term_DM = term_dm_force_like / M;

        // Cooling flow term (pressure / density for acceleration)
        double cool_pressure = rho_cool * v_cool * v_cool;
        double term_cool = cool_pressure / rho_fluid;

        // Term 11 (Tier-1): CP3/PAPER_198 static buoyancy half-gravity (static M variant)
        double term_Ubi = 0.5 * ug1_base;

        // Term 12 (Tier-2): PAPER_198 dynamic compact buoyancy modulation
        double term_F_UBii = -beta_i * ug1_base * omega_g * (M / r) * U_UA * cos(M_PI * t);

        // Term 13 (Tier-3): CP1 outer-frame buoyancy via Virgo Cluster external body (~77 Mpc)
        double term_Ub_i = -beta_i * ug1_base * omega_g * (M_ext_vc / r_ext_vc) * U_UA * cos(M_PI * t);

        // FU diagnostic (monitoring only — not summed into result)
        double FU_diag = -(term2 + term_Ubi) * ug1_base;

        // Total g_NGC1275 (13 terms: 10 original MUGE + 3 buoyancy tiers)
        double result_g = term1 + term_BH + term2 + term3 + term4 + term_q
                        + term_fluid + term_osc + term_DM + term_cool
                        + term_Ubi + term_F_UBii + term_Ub_i;

        if (logging_enabled) {
            std::cout << "[LOG] NGC1275::compute_g_NGC1275(t=" << t << ")\n"
                      << "  term1(base+Hz+B+F)=" << term1
                      << "  term_BH=" << term_BH
                      << "  term2(Ug+fTRZ)=" << term2
                      << "  term3(Lambda)=" << term3
                      << "  term4(EM+UA)=" << term4
                      << "  term_q(quantum)=" << term_q
                      << "  term_fluid=" << term_fluid
                      << "  term_osc=" << term_osc
                      << "  term_DM=" << term_DM
                      << "  term_cool(cooling_flow)=" << term_cool
                      << "  term_Ubi(T1)=" << term_Ubi
                      << "  term_F_UBii(T2)=" << term_F_UBii
                      << "  term_Ub_i(T3_Virgo)=" << term_Ub_i
                      << "  FU_diag=" << FU_diag
                      << "  result_g=" << result_g << "\n";
        }

        return result_g;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "NGC 1275 Parameters:" << std::endl;
        os << "G: " << G << ", M: " << M << ", r: " << r << std::endl;
        os << "Hz: " << Hz << ", B0: " << B0 << ", tau_B: " << tau_B << std::endl;
        os << "f_TRZ: " << f_TRZ << ", M_BH: " << M_BH << ", r_BH: " << r_BH << std::endl;
        os << "F0: " << F0 << ", tau_fil: " << tau_fil << std::endl;
        os << "rho_fluid: " << rho_fluid << ", rho_cool: " << rho_cool << ", v_cool: " << v_cool << std::endl;
        os << "gas_v: " << gas_v << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        os << "ug1_base: " << ug1_base << ", g_BH: " << g_BH << std::endl;
        os << "beta_i: " << beta_i << ", omega_g: " << omega_g << ", U_UA: " << U_UA << std::endl;
        os << "M_ext_vc: " << M_ext_vc << " [Virgo Cluster], r_ext_vc: " << r_ext_vc << " [~77 Mpc Perseus->Virgo]" << std::endl;
    }

    // Example computation at t=50 Myr (for testing)
    double exampleAt50Myr() const {
        double t_example = 50e6 * 3.156e7;
        return compute_g_NGC1275(t_example);
    }

    // UQFF 2.0: Diagnostic logging control
    void setEnableLogging(bool enable) { logging_enabled = enable; }
    bool getLoggingEnabled() const { return logging_enabled; }

    // UQFF 2.0: Dynamic parameter access
    void setDynamicParameter(const std::string& key, double value) {
        dynamic_params[key] = value;
    }
    double getDynamicParameter(const std::string& key) const {
        auto it = dynamic_params.find(key);
        if (it != dynamic_params.end()) return it->second;
        std::cerr << "[WARN] NGC1275: dynamic param '" << key << "' not found.\n";
        return std::numeric_limits<double>::quiet_NaN();
    }

    // UQFF 2.0: State export (44 params + dynamic_params)
    void exportState(const std::string& filename = "NGC1275_state.txt") const {
        std::ofstream ofs(filename);
        if (!ofs) { std::cerr << "[ERR] NGC1275: cannot open " << filename << "\n"; return; }
        ofs << std::scientific << std::setprecision(6);
        ofs << "# NGC1275 UQFF 2.0 State Export — March 16, 2026\n";
        ofs << "G=" << G << "\nM=" << M << "\nr=" << r << "\nHz=" << Hz << "\n";
        ofs << "B0=" << B0 << "\ntau_B=" << tau_B << "\nB_crit=" << B_crit << "\n";
        ofs << "Lambda=" << Lambda << "\nc_light=" << c_light << "\n";
        ofs << "q_charge=" << q_charge << "\ngas_v=" << gas_v << "\nf_TRZ=" << f_TRZ << "\n";
        ofs << "M_BH=" << M_BH << "\nr_BH=" << r_BH << "\n";
        ofs << "F0=" << F0 << "\ntau_fil=" << tau_fil << "\n";
        ofs << "rho_cool=" << rho_cool << "\nv_cool=" << v_cool << "\nrho_fluid=" << rho_fluid << "\n";
        ofs << "rho_vac_UA=" << rho_vac_UA << "\nrho_vac_SCm=" << rho_vac_SCm << "\n";
        ofs << "scale_EM=" << scale_EM << "\nproton_mass=" << proton_mass << "\nz_gal=" << z_gal << "\n";
        ofs << "hbar=" << hbar << "\nt_Hubble=" << t_Hubble << "\nt_Hubble_gyr=" << t_Hubble_gyr << "\n";
        ofs << "delta_x=" << delta_x << "\ndelta_p=" << delta_p << "\nintegral_psi=" << integral_psi << "\n";
        ofs << "A_osc=" << A_osc << "\nk_osc=" << k_osc << "\nomega_osc=" << omega_osc << "\nx_pos=" << x_pos << "\n";
        ofs << "M_DM_factor=" << M_DM_factor << "\ndelta_rho_over_rho=" << delta_rho_over_rho << "\n";
        ofs << "beta_i=" << beta_i << "\nomega_g=" << omega_g << "\nU_UA=" << U_UA << "\n";
        ofs << "M_ext_vc=" << M_ext_vc << "\nr_ext_vc=" << r_ext_vc << "\n";
        ofs << "logging_enabled=" << logging_enabled << "\n";
        for (const auto& kv : dynamic_params)
            ofs << "dyn:" << kv.first << "=" << kv.second << "\n";
        ofs << "ug1_base=" << ug1_base << "\ng_BH=" << g_BH << "\n";
    }

    // UQFF 2.0: Cross-validation template
    template <typename OtherModuleT>
    double cross_validate(const OtherModuleT& other, double t) const {
        return compute_g_NGC1275(t) - other.compute_g_NGC1275(t);
    }
};

#endif // NGC_1275_H