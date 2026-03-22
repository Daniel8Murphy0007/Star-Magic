// MagnetarDualUQFFModule.cpp
// Doc 39b: SGR 1745-2900 dual-mode UQFF.
// Expected results: compressed ~1.782e39 m/s², frequency ~1.773e-9 m/s²
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#include "MagnetarDualUQFFModule.h"
#include "UQFFConstants.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <sstream>

MagnetarDualUQFFModule::MagnetarDualUQFFModule(const std::string& mode)
    : current_mode(mode)
{
    // Universal constants
    variables["G"]          = UQFF::G;
    variables["c"]          = UQFF::c;
    variables["hbar"]       = UQFF::hbar;
    variables["Lambda"]     = UQFF::Lambda;
    variables["q"]          = UQFF::q_e;
    variables["pi"]         = UQFF::pi;
    variables["t_Hubble"]   = UQFF::t_Hubble;
    variables["year_to_s"]  = UQFF::year_to_s;
    variables["H0"]         = UQFF::H0;
    variables["Mpc_to_m"]   = UQFF::Mpc_to_m;
    variables["Omega_m"]    = UQFF::Omega_m;
    variables["Omega_Lambda"]= UQFF::Omega_Lambda;
    variables["m_p"]        = UQFF::m_proton;

    // Magnetar SGR 1745-2900 params
    variables["M"]          = 2.8 * UQFF::M_sun;     // 5.57e30 kg
    variables["M0"]         = variables["M"];
    variables["M_visible"]  = variables["M"];
    variables["M_DM"]       = 0.0;
    variables["r"]          = 1.0e4;                  // m (10 km neutron star radius)
    variables["z"]          = 0.026;

    // Black hole (Sgr A*) influence
    variables["M_BH"]       = 4.0e6 * UQFF::M_sun;
    variables["r_BH"]       = 8.0e9;                  // m (0.026 pc separation)

    // Time
    variables["t"]          = 1.0e3 * UQFF::year_to_s;   // 1 kyr default

    // Fluid (neutron star density)
    variables["rho_fluid"]  = 1.0e17;                 // kg/m^3
    variables["V"]          = 1.0 / variables["rho_fluid"];
    variables["delta_rho"]  = 1.0e-5 * variables["rho_fluid"];
    variables["rho"]        = variables["rho_fluid"];

    // EM
    variables["B"]          = 1.0e11;                 // T  (extreme magnetar)
    variables["B_crit"]     = 1.0e11;                 // T
    variables["tau_decay"]  = 1.0e3 * UQFF::year_to_s;

    // Quantum
    variables["Delta_x"]    = 1.0e-15;               // m (nuclear scale)
    variables["Delta_p"]    = UQFF::hbar / variables["Delta_x"];
    variables["integral_psi_total"] = 1.0;
    variables["f_sc"]       = 10.0;
    variables["f_TRZ"]      = 0.1;

    // Cosmological
    variables["Ug1"] = 0.0;
    variables["Ug4"] = 0.0;

    // --- Frequency mode parameters ---
    variables["F_DPM"]      = 1.702e56;               // A m^2
    variables["V_sys"]      = 5.913e53;               // m^3
    variables["v_exp"]      = 1.0e5;                  // m/s
    variables["f_DPM"]      = 1.0e-8;                 // Hz
    variables["f_THz"]      = 1.0e12;                 // Hz
    variables["rho_vac_UA"] = UQFF::rho_vac_UA;
    variables["rho_vac_SCm"]= UQFF::rho_vac_SCm;
    variables["f_super"]    = 1.0e15;                 // Hz
    variables["f_aether"]   = 1.576e-35;              // Hz
    variables["f_quantum"]  = 1.445e-17;              // Hz
    variables["f_fluid"]    = 8.457e-11;              // Hz
    variables["f_osc"]      = 4.57e14;                // Hz
    variables["f_exp"]      = 2.279e-5;               // Hz
    variables["A"]          = 1.0e-10;
    variables["k"]          = 1.0e20;
    variables["omega"]      = 1.0e15;
    variables["x"]          = 0.0;
}

void MagnetarDualUQFFModule::setMode(const std::string& mode) {
    current_mode = mode;
}

void MagnetarDualUQFFModule::onVariableUpdated(const std::string& name, double value) {
    if (name == "Delta_x" && value > 0.0)
        variables["Delta_p"] = UQFF::hbar / value;
    else if (name == "M") {
        variables["M0"]       = value;
        variables["M_visible"]= value;
    } else if (name == "rho_fluid" && value > 0.0) {
        variables["V"]         = 1.0 / value;
        variables["delta_rho"] = 1.0e-5 * value;
        variables["rho"]       = value;
    }
}

// ===== Compressed mode =====

double MagnetarDualUQFFModule::computeHtz(double z) const {
    return variables.at("H0") * std::sqrt(
        variables.at("Omega_m") * std::pow(1.0 + z, 3) +
        variables.at("Omega_Lambda")) * 1.0e3 / variables.at("Mpc_to_m");
}

// F_env(t) = 1 + M_mag/(Mc²) + D(t) + G M_BH/r_BH²
double MagnetarDualUQFFModule::computeF_env(double t) const {
    double M_mag = 1.0e40;   // J (magnetic energy estimate)
    double D_t   = std::exp(-t / variables.at("tau_decay"));
    double BH_term = variables.at("G") * variables.at("M_BH") /
                     (variables.at("r_BH") * variables.at("r_BH"));
    return 1.0 + (M_mag / (variables.at("M") * variables.at("c") * variables.at("c"))) +
           D_t + BH_term;
}

double MagnetarDualUQFFModule::computeQuantumTerm() const {
    double unc = std::sqrt(variables.at("Delta_x") * variables.at("Delta_p"));
    return (variables.at("hbar") / unc) *
           variables.at("integral_psi_total") *
           (2.0 * variables.at("pi") / variables.at("t_Hubble"));
}

double MagnetarDualUQFFModule::computeFluidTerm(double g_base) const {
    return variables.at("rho_fluid") * variables.at("V") * g_base;
}

double MagnetarDualUQFFModule::computeUgSum() const {
    double G   = variables.at("G");
    double M   = variables.at("M");
    double r   = variables.at("r");
    double Ug1 = G * M / (r * r);
    double Ug3_prime = G * variables.at("M_BH") / (variables.at("r_BH") * variables.at("r_BH"));
    double Ug4 = Ug1 * variables.at("f_sc");
    return Ug1 + Ug3_prime + Ug4;
}

double MagnetarDualUQFFModule::computeDMPertTerm() const {
    double r = variables.at("r");
    double pert = variables.at("delta_rho") / variables.at("rho") +
                  3.0 * variables.at("G") * variables.at("M") / std::pow(r, 3);
    return (variables.at("M_visible") + variables.at("M_DM")) * pert;
}

double MagnetarDualUQFFModule::computeCompressed(double t) {
    double z         = variables.at("z");
    double Hz        = computeHtz(z);
    double expansion = 1.0 + Hz * t;
    double sc_corr   = 1.0 - variables["B"] / variables["B_crit"];
    double f_env     = computeF_env(t);
    double r         = variables["r"];

    double g_base  = (variables["G"] * variables["M"] / (r * r))
                     * expansion * sc_corr * f_env;

    double ug_sum  = computeUgSum();
    double lam     = variables["Lambda"] * variables["c"] * variables["c"] / 3.0;
    double quantum = computeQuantumTerm();
    double fluid   = computeFluidTerm(g_base);
    double dm_pert = computeDMPertTerm();

    return g_base + ug_sum + lam + quantum + fluid + dm_pert;
}

// ===== Frequency mode =====

double MagnetarDualUQFFModule::computeADPM() const {
    double F = variables.at("F_DPM"), f = variables.at("f_DPM");
    double V = variables.at("V_sys"), ve = variables.at("v_exp");
    double q = variables.at("q"),   r  = variables.at("r"),   c = variables.at("c");
    return F * f * V * ve / (q * r * r * c);
}

double MagnetarDualUQFFModule::computeATHz() const {
    return variables.at("q") * variables.at("f_THz") *
           variables.at("v_exp") / variables.at("r");
}

double MagnetarDualUQFFModule::computeAvac_diff() const {
    double vac_r = variables.at("rho_vac_UA") / variables.at("rho_vac_SCm");
    return variables.at("q") * variables.at("B") *
           variables.at("v_exp") / variables.at("m_p") * vac_r;
}

double MagnetarDualUQFFModule::computeAsuper_freq() const {
    double corr = 1.0 - variables.at("B") / variables.at("B_crit");
    return variables.at("f_super") * variables.at("hbar") / variables.at("m_p") * corr;
}

double MagnetarDualUQFFModule::computeAaether_res(double t) const {
    double k = variables.at("k"), x = variables.at("x"), om = variables.at("omega");
    return variables.at("f_aether") * variables.at("c") * std::cos(k * x - om * t);
}

double MagnetarDualUQFFModule::computeUg4i(double t) const {
    double Ug1 = variables.at("G") * variables.at("M") /
                 (variables.at("r") * variables.at("r"));
    double Ug4 = Ug1 * variables.at("f_sc");
    return Ug4 * std::sin(variables.at("omega") * t);
}

double MagnetarDualUQFFModule::computeAquantum_freq() const {
    return variables.at("f_quantum") * variables.at("hbar") / variables.at("Delta_x");
}

double MagnetarDualUQFFModule::computeAAether_freq() const {
    return variables.at("f_aether") * variables.at("c") / variables.at("r");
}

double MagnetarDualUQFFModule::computeAfluid_freq() const {
    return variables.at("f_fluid") * variables.at("rho_fluid") *
           variables.at("V") * variables.at("G");
}

double MagnetarDualUQFFModule::computeOsc_term(double t) const {
    double A  = variables.at("A"), k = variables.at("k");
    double om = variables.at("omega"), x = variables.at("x");
    return 2.0 * A * std::cos(k * x) * std::cos(om * t);
}

double MagnetarDualUQFFModule::computeAexp_freq(double t) const {
    double A  = variables.at("A"), k = variables.at("k");
    double om = variables.at("omega"), x = variables.at("x");
    std::complex<double> e = std::exp(std::complex<double>(0.0, k * x - om * t));
    return (2.0 * variables.at("pi") / 13.8) * A * e.real();
}

double MagnetarDualUQFFModule::computeFTRZ() const {
    return variables.at("f_TRZ");
}

double MagnetarDualUQFFModule::computeFrequency(double t) {
    return computeADPM()
         + computeATHz()
         + computeAvac_diff()
         + computeAsuper_freq()
         + computeAaether_res(t)
         + computeUg4i(t)
         + computeAquantum_freq()
         + computeAAether_freq()
         + computeAfluid_freq()
         + computeOsc_term(t)
         + computeAexp_freq(t)
         + computeFTRZ();
}

double MagnetarDualUQFFModule::computeG(double t) {
    variables["t"] = t;
    if (current_mode == "compressed")  return computeCompressed(t);
    if (current_mode == "frequency")   return computeFrequency(t);
    std::cerr << "[MagnetarDualUQFFModule] Unknown mode: " << current_mode << "\n";
    return 0.0;
}

std::string MagnetarDualUQFFModule::getEquationText() const {
    if (current_mode == "compressed") {
        return
            "COMPRESSED MODE — SGR 1745-2900:\n"
            "g(r,t) = [G M/r²] * (1+H(t,z)) * (1-B/B_crit) * F_env(t)\n"
            "       + (Ug1 + Ug3'_BH + Ug4) + (Λc²/3) + quantum_psi + fluid + DM_pert\n"
            "F_env(t) = 1 + M_mag/(Mc²) + exp(-t/τ_decay) + G M_BH/r_BH²\n"
            "Ug3' = G M_BH/r_BH²  (Sgr A* influence at r_BH=8e9 m)\n"
            "Expected: ~1.782e39 m/s²  (SM-like dominant via expansion+SC+BH)";
    }
    return
        "FREQUENCY MODE — SGR 1745-2900:\n"
        "g = a_DPM + a_THz + a_vac_diff + a_super_freq + a_aether_res\n"
        "  + Ug4i + a_quantum_freq + a_Aether_freq + a_fluid_freq + Osc_term\n"
        "  + a_exp_freq + f_TRZ\n"
        "a_DPM = F_DPM f_DPM V_sys v_exp / (q r² c)    [pseudo-monopole]\n"
        "a_THz = q f_THz v_exp / r                      [THz coupling]\n"
        "a_vac_diff = q B v_exp/m_p * (ρ_vac_UA/ρ_vac_SCm)  [vacuum ratio]\n"
        "a_super_freq = f_super ℏ/m_p * (1-B/B_crit)   [superconductivity]\n"
        "Expected: ~1.773e-9 m/s²  (pure resonance, no SM illusion)";
}
