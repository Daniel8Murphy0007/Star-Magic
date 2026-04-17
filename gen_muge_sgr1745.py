"""
gen_muge_sgr1745.py — Generator for MagnetarSGR1745_2900.h
Module 4 in MUGE series: SGR 1745-2900 magnetar near Sgr A* gravity evolution.
Unique physics: g_BH (BH proximity), M_mag (magnetic energy), cum_D (cumulative decay energy),
                f_sc (superconductivity), H(z) redshift-dependent Hubble parameter.

Run:  python gen_muge_sgr1745.py
Output: MagnetarSGR1745_2900.h
"""
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "MagnetarSGR1745_2900.h")


def get_content():
    return r"""/**
 * ================================================================================================
 * Header: MagnetarSGR1745_2900.h
 *
 * Description: C++ Module for SGR 1745-2900 Magnetar (near Sgr A*) Class (Module 4)
 *              UQFF simulations — gravity evolution including BH proximity, magnetic energy,
 *              cumulative luminosity decay energy, and superconductivity.
 *
 * Unique Terms vs SGR 0501:
 *   - g_BH = G*M_BH/r_BH^2 (black hole proximity at 2.83e16 m)
 *   - M_mag = (B^2/2*mu_0)*V (magnetic energy as effective acceleration)
 *   - cum_D(t) = L0*tau_decay*(1-exp(-t/tau_decay)) cumulative decay energy
 *   - f_sc = 1 - B/B_crit (superconductivity factor)
 *   - H(z) = H0*sqrt(0.3*(1+z)^3 + 0.7) in s^-1
 *
 * Key Parameters:
 *   M=1.4 M_sun, r=12 km, B=2e10 T, P=3.76 s, L0=5e28 W, tau_decay=3.5 yr,
 *   M_BH=4e6 M_sun, r_BH=2.83e16 m, z=0.01, Hz=2.269e-18 s^-1
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef MAGNETAR_SGR1745_2900_H
#define MAGNETAR_SGR1745_2900_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

namespace UQFF {

class MagnetarSGR1745_2900 {
private:
    double G, M, r;
    double B, B_crit;
    double P;               // Spin period (s)
    double L0;              // Initial luminosity (W)
    double tau_decay;       // Luminosity decay timescale (s)
    double M_BH;            // Sgr A* BH mass (kg)
    double r_BH;            // Distance to BH (m)
    double Lambda, c_light;
    double Hz;              // H(z) at z=0.01 (s^-1)
    double z_gal;           // Redshift
    double q_charge, gas_v;
    double f_TRZ;
    double rho_vac_UA, rho_vac_SCm, scale_EM, proton_mass;
    double mu_0;            // Permeability of free space
    // Full-term parameters
    double hbar, t_Hubble, t_Hubble_gyr;
    double delta_x, delta_p, integral_psi;
    double rho_fluid, A_osc, k_osc, omega_osc, x_pos;
    double M_DM_factor, delta_rho_over_rho;
    // Cached
    double ug1_base, g_BH;

public:
    MagnetarSGR1745_2900() { initializeDefaults(); }
    ~MagnetarSGR1745_2900() {}

    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M = 1.4 * M_sun;
        r = 1.2e4;           // 12 km
        B = 2.0e10;
        B_crit = 4.4e13;
        P = 3.76;
        L0 = 5.0e28;
        tau_decay = 3.5 * 3.156e7;  // 3.5 years
        M_BH = 4.0e6 * M_sun;
        r_BH = 2.83e16;
        Lambda = 1.1e-52;
        c_light = 3.0e8;
        z_gal = 0.01;
        double H0 = 70.0;  // km/s/Mpc
        double Hz_kms = H0 * std::sqrt(0.3 * std::pow(1.0 + z_gal, 3) + 0.7);
        Hz = Hz_kms * 1000.0 / 3.086e22;  // Convert to s^-1
        q_charge = 1.602e-19;
        gas_v = 1.0e5;
        f_TRZ = 0.1;
        rho_vac_UA = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        scale_EM = 1.0e-12;
        proton_mass = 1.673e-27;
        mu_0 = 4.0 * M_PI * 1.0e-7;

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
        ug1_base = B_field * r * G * M;
        g_BH = (G * M_BH) / (r_BH * r_BH);
    }

    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G")         { G = newValue; }
        else if (varName == "M")    { M = newValue; }
        else if (varName == "r")    { r = newValue; }
        else if (varName == "B")    { B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Hz")   { Hz = newValue; }
        else if (varName == "L0")   { L0 = newValue; }
        else if (varName == "tau_decay") { tau_decay = newValue; }
        else if (varName == "M_BH") { M_BH = newValue; }
        else if (varName == "r_BH") { r_BH = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'.\n";
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
        if (varName == "B") return B;
        if (varName == "Hz") return Hz;
        if (varName == "f_TRZ") return f_TRZ;
        std::cerr << "Unknown variable '" << varName << "'.\n";
        return 0.0;
    }

    // f_sc: superconductivity
    double f_sc() const { return 1.0 - B / B_crit; }

    // M_mag: magnetic energy as effective acceleration (J/m^3 -> m/s^2 via / rho_fluid)
    double compute_M_mag() const {
        double V = (4.0 / 3.0) * M_PI * r * r * r;
        double mag_energy = (B * B / (2.0 * mu_0)) * V;
        return mag_energy / (M * r);  // Effective acceleration
    }

    // cum_D(t): cumulative luminosity decay (J -> scaled accel)
    double cum_D(double t) const {
        return L0 * tau_decay * (1.0 - std::exp(-t / tau_decay)) / (M * c_light * c_light);
    }

    // Ug terms
    double compute_Ug() const {
        double fsc = f_sc();
        double Ug1 = ug1_base;
        double Ug4 = Ug1 * fsc;
        return (Ug1 + Ug4) * (1.0 + f_TRZ);
    }

    double compute_V() const { return (4.0 / 3.0) * M_PI * r * r * r; }

    double compute_g_SGR1745(double t) const {
        if (t < 0.0) { std::cerr << "Error: t must be non-negative.\n"; return 0.0; }

        double fsc = f_sc();
        double corr_H = 1.0 + Hz * t;
        double term1 = ug1_base * corr_H * fsc;

        // BH proximity
        double term_BH = g_BH;

        // UQFF Ug
        double term2 = compute_Ug();

        // Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // EM
        double cross_vB = gas_v * B;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1.0 + (rho_vac_UA / rho_vac_SCm);
        double term4 = em_base * corr_UA * scale_EM;

        // Quantum
        double sqrt_unc = std::sqrt(delta_x * delta_p);
        double term_q = (hbar / sqrt_unc) * integral_psi * (2.0 * M_PI / t_Hubble);

        // Fluid
        double V = compute_V();
        double term_fluid = (rho_fluid * V * ug1_base) / M;

        // Oscillatory
        double term_osc = 2.0 * A_osc * std::cos(k_osc * x_pos) * std::cos(omega_osc * t)
                        + (2.0 * M_PI / t_Hubble_gyr) * A_osc * std::cos(k_osc * x_pos - omega_osc * t);

        // DM
        double M_dm = M * M_DM_factor;
        double term_DM = ((M + M_dm) * (delta_rho_over_rho + 3.0 * G * M / (r * r * r))) / M;

        // Magnetic energy term (unique to SGR1745)
        double term_Mmag = compute_M_mag();

        // Cumulative decay energy (unique to SGR1745)
        double term_cumD = cum_D(t);

        return term1 + term_BH + term2 + term3 + term4 + term_q + term_fluid + term_osc
             + term_DM + term_Mmag + term_cumD;
    }

    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(6);
        os << "SGR 1745-2900 Parameters:\n";
        os << "  M=" << M << "  r=" << r << "  B=" << B << "\n";
        os << "  Hz=" << Hz << "  M_BH=" << M_BH << "  r_BH=" << r_BH << "\n";
        os << "  L0=" << L0 << "  tau_decay=" << tau_decay << "\n";
        os << "  f_sc=" << f_sc() << "  g_BH=" << g_BH << "\n";
    }

    double exampleAt10yr() const { return compute_g_SGR1745(10.0 * 3.156e7); }
};

}  // namespace UQFF

#endif  // MAGNETAR_SGR1745_2900_H
"""


def main():
    print("gen_muge_sgr1745.py — Generating MagnetarSGR1745_2900.h ...")
    content = get_content()
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  Written: {OUTPUT_FILE}  ({len(content.splitlines())} lines)")


if __name__ == "__main__":
    main()
