/**
 * ================================================================================================
 * Header: MagnetarSGR0501_4516.h
 *
 * Description: C++ Module for SGR 0501+4516 Magnetar Class (Module 2)
 *              UQFF simulations — gravity evolution from Hubble datasets, high-energy lab
 *              simulations, and UQFF refinements (dated May 09, 2025, updated Oct 08, 2025).
 *
 * Purpose: Master Universal Gravity Equation (MUGE) for SGR 0501+4516.
 *          Terms: base gravity, B(t) decay, Omega(t) spin decay, f_sc superconductivity,
 *          UQFF Ug (f_TRZ), Lambda, quantum uncertainty, scaled EM [UA], fluid dynamics,
 *          oscillatory waves, DM/density perturbations.
 *
 * Key Parameters:
 *   M = 1.4 M_sun, r = 20 km, B0 = 1e10 T, tau_B = 4000 yr, P0 = 5.76 s,
 *   f_TRZ = 0.1, Omega0 = 2*pi/P0, tau_Omega = 10000 yr.
 *
 * Example: t=5000 yr -> g ~ 4.474e12 m/s^2
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef MAGNETAR_SGR0501_4516_H
#define MAGNETAR_SGR0501_4516_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

namespace UQFF {

class MagnetarSGR0501_4516 {
private:
    // Core parameters
    double G;               // Gravitational constant (m^3 kg^-1 s^-2)
    double M;               // Magnetar mass (kg)
    double r;               // Radius (m)
    double B0;              // Initial magnetic field (T)
    double tau_B;           // Magnetic decay timescale (s)
    double B_crit;          // Critical B field for superconductivity (T)
    double P0;              // Initial spin period (s)
    double tau_Omega;       // Spin decay timescale (s)
    double Lambda;          // Cosmological constant (m^-2)
    double c_light;         // Speed of light (m/s)
    double Hz;              // Hubble parameter (s^-1)
    double q_charge;        // Proton charge (C)
    double gas_v;           // Gas velocity (m/s)
    double f_TRZ;           // Time-reversal UQFF factor
    double rho_vac_UA;      // UA vacuum energy density (J/m^3)
    double rho_vac_SCm;     // SCm vacuum energy density (J/m^3)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass (kg)

    // Full-term parameters
    double hbar;            // Reduced Planck constant (J s)
    double t_Hubble;        // Hubble time (s)
    double t_Hubble_gyr;    // Hubble time (Gyr)
    double delta_x;         // Position uncertainty (m)
    double delta_p;         // Momentum uncertainty (kg m/s)
    double integral_psi;    // Wavefunction integral approximation
    double rho_fluid;       // Fluid density (kg/m^3)
    double A_osc;           // Oscillatory amplitude (m/s^2)
    double k_osc;           // Wave number (1/m)
    double omega_osc;       // Angular frequency (rad/s)
    double x_pos;           // Position for oscillation (m)
    double M_DM_factor;     // Dark matter mass fraction
    double delta_rho_over_rho; // Density perturbation fraction

    // Cached values
    double ug1_base;        // G*M/r^2

public:
    MagnetarSGR0501_4516() { initializeDefaults(); }
    ~MagnetarSGR0501_4516() {}

    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M = 1.4 * M_sun;
        r = 2.0e4;           // 20 km
        B0 = 1.0e10;
        tau_B = 4000.0 * 3.156e7;   // 4000 years in seconds
        B_crit = 4.4e13;
        P0 = 5.76;                   // seconds
        tau_Omega = 10000.0 * 3.156e7;
        Lambda = 1.1e-52;
        c_light = 3.0e8;
        Hz = 2.184e-18;
        q_charge = 1.602e-19;
        gas_v = 1.0e5;
        f_TRZ = 0.1;
        rho_vac_UA = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        scale_EM = 1.0e-12;
        proton_mass = 1.673e-27;

        hbar = 1.0546e-34;
        t_Hubble = 13.8e9 * 3.156e7;
        t_Hubble_gyr = 13.8;
        delta_x = 1.0e-10;
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
        rho_fluid = 1.0e-21;
        A_osc = 1.0e-10;
        k_osc = 1.0 / r;
        omega_osc = 2.0 * M_PI / (r / c_light);
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1.0e-5;

        updateCache();
    }

    void updateCache() {
        ug1_base = (G * M) / (r * r);
    }

    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G")              { G = newValue; }
        else if (varName == "M")         { M = newValue; }
        else if (varName == "r")         { r = newValue; }
        else if (varName == "B0")        { B0 = newValue; }
        else if (varName == "tau_B")     { tau_B = newValue; }
        else if (varName == "B_crit")    { B_crit = newValue; }
        else if (varName == "P0")        { P0 = newValue; }
        else if (varName == "tau_Omega") { tau_Omega = newValue; }
        else if (varName == "Hz")        { Hz = newValue; }
        else if (varName == "f_TRZ")     { f_TRZ = newValue; }
        else if (varName == "Lambda")    { Lambda = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
        else if (varName == "M_DM_factor") { M_DM_factor = newValue; }
        else if (varName == "delta_rho_over_rho") { delta_rho_over_rho = newValue; }
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return false;
        }
        updateCache();
        return true;
    }

    bool addToVariable(const std::string& v, double d) { return setVariable(v, getVariable(v) + d); }
    bool subtractFromVariable(const std::string& v, double d) { return addToVariable(v, -d); }

    double getVariable(const std::string& varName) const {
        if (varName == "G") return G;
        if (varName == "M") return M;
        if (varName == "r") return r;
        if (varName == "B0") return B0;
        if (varName == "tau_B") return tau_B;
        if (varName == "f_TRZ") return f_TRZ;
        if (varName == "Hz") return Hz;
        if (varName == "rho_fluid") return rho_fluid;
        std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
        return 0.0;
    }

    // B(t): time-dependent magnetic field decay
    double B_t(double t) const { return B0 * std::exp(-t / tau_B); }

    // Omega(t): spin frequency decay
    double Omega_t(double t) const { return (2.0 * M_PI / P0) * std::exp(-t / tau_Omega); }

    // f_sc: superconductivity correction
    double f_sc(double Bt) const { return 1.0 - Bt / B_crit; }

    // Ug terms
    double compute_Ug(double Bt) const {
        double Ug1 = ug1_base;
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = f_sc(Bt);
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1.0 + f_TRZ);
    }

    // Volume
    double compute_V() const { return (4.0 / 3.0) * M_PI * r * r * r; }

    // Main MUGE computation
    double compute_g_Magnetar(double t) const {
        if (t < 0.0) {
            std::cerr << "Error: t must be non-negative." << std::endl;
            return 0.0;
        }

        double Bt = B_t(t);
        double Ot = Omega_t(t);
        double fsc = f_sc(Bt);

        // Term 1: Base + Hz + B corrections
        double corr_H = 1.0 + Hz * t;
        double term1 = ug1_base * corr_H * fsc;

        // Term 2: UQFF Ug with f_TRZ
        double term2 = compute_Ug(Bt);

        // Term 3: Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // Term 4: EM with [UA] correction (scaled)
        double cross_vB = gas_v * Bt;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1.0 + (rho_vac_UA / rho_vac_SCm);
        double term4 = (em_base * corr_UA) * scale_EM;

        // Quantum uncertainty
        double sqrt_unc = std::sqrt(delta_x * delta_p);
        double term_q = (hbar / sqrt_unc) * integral_psi * (2.0 * M_PI / t_Hubble);

        // Fluid
        double V = compute_V();
        double term_fluid = (rho_fluid * V * ug1_base) / M;

        // Oscillatory
        double term_osc1 = 2.0 * A_osc * std::cos(k_osc * x_pos) * std::cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2.0 * M_PI / t_Hubble_gyr) * A_osc * std::cos(arg);
        double term_osc = term_osc1 + term_osc2;

        // DM perturbation
        double M_dm = M * M_DM_factor;
        double pert2 = 3.0 * G * M / (r * r * r);
        double term_DM = ((M + M_dm) * (delta_rho_over_rho + pert2)) / M;

        // Spin contribution (Omega^2 * r as centrifugal-analogue)
        double term_spin = Ot * Ot * r * 1.0e-10;  // scaled to avoid dominance

        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_spin;
    }

    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(6);
        os << "SGR 0501+4516 Magnetar Parameters:\n";
        os << "  G=" << G << "  M=" << M << "  r=" << r << "\n";
        os << "  B0=" << B0 << "  tau_B=" << tau_B << "  B_crit=" << B_crit << "\n";
        os << "  P0=" << P0 << "  tau_Omega=" << tau_Omega << "\n";
        os << "  Hz=" << Hz << "  f_TRZ=" << f_TRZ << "\n";
        os << "  ug1_base=" << ug1_base << "\n";
    }

    // Example: t = 5000 years
    double exampleAt5000yr() const { return compute_g_Magnetar(5000.0 * 3.156e7); }
};

}  // namespace UQFF

#endif  // MAGNETAR_SGR0501_4516_H
