"""
gen_muge_ngc2525.py — Generator for GalaxyNGC2525.h
Module 10: NGC 2525 barred spiral galaxy with Type Ia supernova (SN 2018gv).
Unique physics: SN mass loss term -(G*M_SN(t))/r^2, BH influence term.

Run:  python gen_muge_ngc2525.py
Output: GalaxyNGC2525.h
"""
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "GalaxyNGC2525.h")


def get_content():
    return r"""/**
 * ================================================================================================
 * Header: GalaxyNGC2525.h
 *
 * Description: C++ Module for NGC 2525 Barred Spiral Galaxy + SN 2018gv Class (Module 10)
 *              UQFF simulations — galaxy gravity with Ia supernova mass loss and BH influence.
 *
 * Unique Terms:
 *   - SN mass loss: term_SN = -(G * M_SN(t)) / r^2  (negative: mass leaves system)
 *   - M_SN(t) = M_SN0 * exp(-t / tau_SN)
 *   - BH term: g_BH = G * M_BH / r_BH^2
 *
 * Key Parameters:
 *   M=1e10 M_sun, M_BH=2.25e7 M_sun, r_BH=1.496e11 m (~1 AU),
 *   M_SN0=1.4 M_sun, tau_SN=1 yr, z=0.016
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef GALAXY_NGC2525_H
#define GALAXY_NGC2525_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

namespace UQFF {

class GalaxyNGC2525 {
private:
    double G, M, r;
    double B, B_crit;
    double Lambda, c_light;
    double Hz;              // H(z) at z=0.016
    double z_gal;
    double q_charge, gas_v, f_TRZ;
    double M_BH;            // Black hole mass (kg)
    double r_BH;            // BH influence radius (m)
    double M_SN0;           // Initial SN ejecta mass (kg)
    double tau_SN;          // SN mass loss timescale (s)
    double rho_vac_UA, rho_vac_SCm, scale_EM, proton_mass;
    double rho_fluid;
    double hbar, t_Hubble, t_Hubble_gyr;
    double delta_x, delta_p, integral_psi;
    double A_osc, k_osc, omega_osc, x_pos;
    double M_DM_factor, delta_rho_over_rho;
    // Cached
    double ug1_base, g_BH;

public:
    GalaxyNGC2525() { initializeDefaults(); }
    ~GalaxyNGC2525() {}

    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M = 1.0e10 * M_sun;
        double ly = 9.461e15;
        r = 50000.0 * ly;   // ~50 kly galaxy radius
        B = 1.0e-8;
        B_crit = 1.0e11;
        Lambda = 1.1e-52;
        c_light = 3.0e8;
        z_gal = 0.016;
        double H0 = 70.0;
        double Hz_kms = H0 * std::sqrt(0.3 * std::pow(1.0 + z_gal, 3) + 0.7);
        Hz = Hz_kms * 1000.0 / 3.086e22;
        q_charge = 1.602e-19;
        gas_v = 1.0e5;
        f_TRZ = 0.1;
        M_BH = 2.25e7 * M_sun;
        r_BH = 1.496e11;    // ~1 AU
        M_SN0 = 1.4 * M_sun;
        tau_SN = 1.0 * 3.156e7;  // 1 year
        rho_vac_UA = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        scale_EM = 1.0e-12;
        proton_mass = 1.673e-27;
        rho_fluid = 1.0e-26;

        hbar = 1.0546e-34;
        t_Hubble = 13.8e9 * 3.156e7;
        t_Hubble_gyr = 13.8;
        delta_x = 1.0e-10;
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
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
        g_BH = (G * M_BH) / (r_BH * r_BH);
    }

    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G")        { G = newValue; }
        else if (varName == "M")   { M = newValue; }
        else if (varName == "r")   { r = newValue; }
        else if (varName == "Hz")  { Hz = newValue; }
        else if (varName == "M_BH")  { M_BH = newValue; }
        else if (varName == "r_BH")  { r_BH = newValue; }
        else if (varName == "M_SN0") { M_SN0 = newValue; }
        else if (varName == "tau_SN") { tau_SN = newValue; }
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
        if (varName == "M_BH") return M_BH;
        if (varName == "M_SN0") return M_SN0;
        if (varName == "tau_SN") return tau_SN;
        std::cerr << "Unknown variable '" << varName << "'.\n";
        return 0.0;
    }

    // M_SN(t): SN ejecta mass decay
    double M_SN_t(double t) const { return M_SN0 * std::exp(-t / tau_SN); }

    double compute_Ug() const {
        double corr_B = 1.0 - B / B_crit;
        return (ug1_base + ug1_base * corr_B) * (1.0 + f_TRZ);
    }

    double compute_V() const { return (4.0 / 3.0) * M_PI * r * r * r; }

    double compute_g_NGC2525(double t) const {
        if (t < 0.0) { std::cerr << "Error: t must be non-negative.\n"; return 0.0; }

        double MSNt = M_SN_t(t);

        // Term 1: Base + Hz + B
        double corr_H = 1.0 + Hz * t;
        double corr_B = 1.0 - B / B_crit;
        double term1 = ug1_base * corr_H * corr_B;

        // BH term (unique)
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

        // SN mass loss (negative — ejecta reduces effective gravity, unique to NGC2525)
        double term_SN = -(G * MSNt) / (r * r);

        return term1 + term_BH + term2 + term3 + term4 + term_q + term_fluid
             + term_osc + term_DM + term_SN;
    }

    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(6);
        os << "NGC 2525 Parameters:\n";
        os << "  M=" << M << "  r=" << r << "  z_gal=" << z_gal << "\n";
        os << "  M_BH=" << M_BH << "  r_BH=" << r_BH << "\n";
        os << "  M_SN0=" << M_SN0 << "  tau_SN=" << tau_SN << "\n";
        os << "  g_BH=" << g_BH << "\n";
    }

    double exampleAt7yr() const { return compute_g_NGC2525(7.0 * 3.156e7); }
};

}  // namespace UQFF

#endif  // GALAXY_NGC2525_H
"""


def main():
    print("gen_muge_ngc2525.py — Generating GalaxyNGC2525.h ...")
    content = get_content()
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  Written: {OUTPUT_FILE}  ({len(content.splitlines())} lines)")


if __name__ == "__main__":
    main()
