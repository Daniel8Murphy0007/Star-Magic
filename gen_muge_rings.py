"""
gen_muge_rings.py — Generator for RingsOfRelativity.h
Module 9: GAL-CLUS-022058s Einstein ring / gravitational lensing system.
Unique physics: L_t = (G*M)/(c^2*r) * L_factor as lensing amplification in base correction,
                redshift-dependent H(z), z=0.5.

Run:  python gen_muge_rings.py
Output: RingsOfRelativity.h
"""
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "RingsOfRelativity.h")


def get_content():
    return r"""/**
 * ================================================================================================
 * Header: RingsOfRelativity.h
 *
 * Description: C++ Module for GAL-CLUS-022058s Einstein Ring / Gravitational Lensing (Module 9)
 *              UQFF simulations — lensing cluster gravity with relativistic amplification.
 *
 * Unique Terms:
 *   - Lensing amplification: L_t = (G*M) / (c^2*r) * L_factor
 *   - Applied as additional multiplicative correction to base gravity term
 *   - H(z) at z=0.5: Hz = H0*sqrt(0.3*(1+z)^3 + 0.7)*1e3/3.086e22
 *
 * Key Parameters:
 *   M=1e14 M_sun, r=10 kpc=3.086e20 m, z=0.5, L_factor=0.67
 *   Hz=2.42e-18 s^-1
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef RINGS_OF_RELATIVITY_H
#define RINGS_OF_RELATIVITY_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

namespace UQFF {

class RingsOfRelativity {
private:
    double G, M, r;
    double B, B_crit;
    double Lambda, c_light;
    double Hz;              // H(z) at z=0.5
    double z_gal;
    double L_factor;        // Lensing geometry factor (0.67 for GAL-CLUS-022058s)
    double q_charge, gas_v, f_TRZ;
    double rho_vac_UA, rho_vac_SCm, scale_EM, proton_mass;
    double rho_fluid;
    double hbar, t_Hubble, t_Hubble_gyr;
    double delta_x, delta_p, integral_psi;
    double A_osc, k_osc, omega_osc, x_pos;
    double M_DM_factor, delta_rho_over_rho;
    double ug1_base;

public:
    RingsOfRelativity() { initializeDefaults(); }
    ~RingsOfRelativity() {}

    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M = 1.0e14 * M_sun;
        r = 3.086e20;        // 10 kpc
        B = 1.0e-7;
        B_crit = 1.0e11;
        Lambda = 1.1e-52;
        c_light = 3.0e8;
        z_gal = 0.5;
        double H0 = 70.0;
        double Hz_kms = H0 * std::sqrt(0.3 * std::pow(1.0 + z_gal, 3) + 0.7);
        Hz = Hz_kms * 1000.0 / 3.086e22;
        L_factor = 0.67;
        q_charge = 1.602e-19;
        gas_v = 1.0e5;
        f_TRZ = 0.1;
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

    void updateCache() { ug1_base = B_field * r * G * M; }

    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G")           { G = newValue; }
        else if (varName == "M")      { M = newValue; }
        else if (varName == "r")      { r = newValue; }
        else if (varName == "Hz")     { Hz = newValue; }
        else if (varName == "z_gal")  { z_gal = newValue; }
        else if (varName == "L_factor") { L_factor = newValue; }
        else if (varName == "f_TRZ")  { f_TRZ = newValue; }
        else if (varName == "B")      { B = newValue; }
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
        if (varName == "Hz") return Hz;
        if (varName == "L_factor") return L_factor;
        std::cerr << "Unknown variable '" << varName << "'.\n";
        return 0.0;
    }

    // L_t: lensing amplification factor (unique to RingsOfRelativity)
    double L_t() const {
        return (G * M) / (c_light * c_light * r) * L_factor;
    }

    double compute_Ug() const {
        double corr_B = 1.0 - B / B_crit;
        return (ug1_base + ug1_base * corr_B) * (1.0 + f_TRZ);
    }

    double compute_V() const { return (4.0 / 3.0) * M_PI * r * r * r; }

    double compute_g_Rings(double t) const {
        if (t < 0.0) { std::cerr << "Error: t must be non-negative.\n"; return 0.0; }

        double Lt = L_t();

        // Term 1: Base + Hz + B, with lensing amplification (1 + L_t)
        double corr_H = 1.0 + Hz * t;
        double corr_B = 1.0 - B / B_crit;
        double term1 = ug1_base * corr_H * corr_B * (1.0 + Lt);

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

        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM;
    }

    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(6);
        os << "Rings of Relativity (GAL-CLUS-022058s) Parameters:\n";
        os << "  M=" << M << "  r=" << r << "  z_gal=" << z_gal << "\n";
        os << "  Hz=" << Hz << "  L_factor=" << L_factor << "  L_t=" << L_t() << "\n";
    }

    double exampleAt100Myr() const { return compute_g_Rings(100.0e6 * 3.156e7); }
};

}  // namespace UQFF

#endif  // RINGS_OF_RELATIVITY_H
"""


def main():
    print("gen_muge_rings.py — Generating RingsOfRelativity.h ...")
    content = get_content()
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  Written: {OUTPUT_FILE}  ({len(content.splitlines())} lines)")


if __name__ == "__main__":
    main()
