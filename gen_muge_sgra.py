"""
gen_muge_sgra.py — Generator for SMBHSgrAStar.h
Module 5 in MUGE series: Sagittarius A* SMBH gravity evolution.
Unique physics: M(t) mass accretion, B in Gauss->Tesla conversion, spin_factor=0.3,
                sin(30 deg) precession on DM, H(z) at z=0.

Run:  python gen_muge_sgra.py
Output: SMBHSgrAStar.h
"""
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "SMBHSgrAStar.h")


def get_content():
    return r"""/**
 * ================================================================================================
 * Header: SMBHSgrAStar.h
 *
 * Description: C++ Module for Sagittarius A* SMBH Class (Module 5)
 *              UQFF simulations — SMBH gravity evolution with mass accretion,
 *              spin effects, and Galactic Center environment.
 *
 * Unique Terms:
 *   - M(t) = M_initial*(1 + M_dot_0*exp(-t/tau_acc))  mass accretion growth
 *   - B in Gauss converted to Tesla: B_T = B_G*1e-4
 *   - spin_factor = 0.3 applied to Omega calculation
 *   - sin(30 deg) precession applied to DM perturbation term
 *
 * Key Parameters:
 *   M=4.3e6 M_sun, r=1.27e10 m, B0=1e4 G, spin_factor=0.3,
 *   M_dot_0=0.01, tau_acc=9 Gyr, z~0
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef SMBH_SGR_A_STAR_H
#define SMBH_SGR_A_STAR_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

namespace UQFF {

class SMBHSgrAStar {
private:
    double G, M_initial, r;
    double B0_G;            // Magnetic field in Gauss
    double B_crit;
    double Lambda, c_light;
    double Hz;              // Hubble param (s^-1)
    double q_charge, gas_v;
    double f_TRZ;
    double spin_factor;     // 0.3 for Sgr A* spin
    double M_dot_0;         // Initial accretion rate (dimensionless)
    double tau_acc;         // Accretion timescale (s)
    double rho_vac_UA, rho_vac_SCm, scale_EM, proton_mass;

    // Full-term parameters
    double hbar, t_Hubble, t_Hubble_gyr;
    double delta_x, delta_p, integral_psi;
    double rho_fluid, A_osc, k_osc, omega_osc, x_pos;
    double M_DM_factor, delta_rho_over_rho;
    double precession_angle; // radians (30 deg)

    // Cached
    double ug1_base_initial;

public:
    SMBHSgrAStar() { initializeDefaults(); }
    ~SMBHSgrAStar() {}

    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M_initial = 4.3e6 * M_sun;
        r = 1.27e10;         // ~0.08 AU (Schwarzschild radius vicinity)
        B0_G = 1.0e4;        // Gauss at center
        B_crit = 4.4e13;
        Lambda = 1.1e-52;
        c_light = 3.0e8;
        Hz = 2.184e-18;      // ~H0, z~0
        q_charge = 1.602e-19;
        gas_v = 1.0e5;
        f_TRZ = 0.1;
        spin_factor = 0.3;
        M_dot_0 = 0.01;
        tau_acc = 9.0e9 * 3.156e7;  // 9 Gyr
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
        precession_angle = 30.0 * M_PI / 180.0;  // 30 degrees

        updateCache();
    }

    void updateCache() {
        ug1_base_initial = B_field * r * G * M_initial;
    }

    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G")            { G = newValue; }
        else if (varName == "M_initial") { M_initial = newValue; }
        else if (varName == "r")       { r = newValue; }
        else if (varName == "B0_G")    { B0_G = newValue; }
        else if (varName == "Hz")      { Hz = newValue; }
        else if (varName == "f_TRZ")   { f_TRZ = newValue; }
        else if (varName == "spin_factor") { spin_factor = newValue; }
        else if (varName == "M_dot_0") { M_dot_0 = newValue; }
        else if (varName == "tau_acc") { tau_acc = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
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
        if (varName == "M_initial") return M_initial;
        if (varName == "r") return r;
        if (varName == "B0_G") return B0_G;
        if (varName == "Hz") return Hz;
        if (varName == "f_TRZ") return f_TRZ;
        std::cerr << "Unknown variable '" << varName << "'.\n";
        return 0.0;
    }

    // M(t): mass with accretion growth
    double M_t(double t) const {
        return M_initial * (1.0 + M_dot_0 * std::exp(-t / tau_acc));
    }

    // B(t): in Tesla (convert from Gauss, assume static here)
    double B_T() const { return B0_G * 1.0e-4; }

    // Ug terms using M(t)
    double compute_Ug(double Mt) const {
        // DPM-seeded: mu_s x grad(M_s/r) (not Newtonian GM/r^2)\ndouble ug1 = B_field * r * G * Mt;
        double corr_B = 1.0 - B_T() / B_crit;
        double Ug4 = ug1 * corr_B;
        return (ug1 + Ug4) * (1.0 + f_TRZ);
    }

    double compute_V() const { return (4.0 / 3.0) * M_PI * r * r * r; }

    double compute_g_SgrAStar(double t) const {
        if (t < 0.0) { std::cerr << "Error: t must be non-negative.\n"; return 0.0; }

        double Mt = M_t(t);
        double Bt = B_T();
        double ug1_t = B_field * r * G * Mt;

        // Spin-modified Omega (spin_factor reduces effective Omega)
        double Omega_t = (2.0 * M_PI / 0.1) * spin_factor;  // ~0.1 s period, spin_factor=0.3

        // Term 1: Base + Hz + B
        double corr_H = 1.0 + Hz * t;
        double corr_B = 1.0 - Bt / B_crit;
        double term1 = ug1_t * corr_H * corr_B;

        // UQFF Ug
        double term2 = compute_Ug(Mt);

        // Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // EM
        double cross_vB = gas_v * Bt;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1.0 + (rho_vac_UA / rho_vac_SCm);
        double term4 = em_base * corr_UA * scale_EM;

        // Quantum
        double sqrt_unc = std::sqrt(delta_x * delta_p);
        double term_q = (hbar / sqrt_unc) * integral_psi * (2.0 * M_PI / t_Hubble);

        // Fluid
        double V = compute_V();
        double term_fluid = (rho_fluid * V * ug1_t) / Mt;

        // Oscillatory
        double term_osc = 2.0 * A_osc * std::cos(k_osc * x_pos) * std::cos(omega_osc * t)
                        + (2.0 * M_PI / t_Hubble_gyr) * A_osc * std::cos(k_osc * x_pos - omega_osc * t);

        // DM with sin(30 deg) precession (unique to SgrA*)
        double M_dm = Mt * M_DM_factor;
        double pert2 = 3.0 * G * Mt / (r * r * r);
        double precession_proj = std::sin(precession_angle);
        double term_DM = ((Mt + M_dm) * (delta_rho_over_rho + pert2) * precession_proj) / Mt;

        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM;
    }

    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(6);
        os << "Sgr A* SMBH Parameters:\n";
        os << "  M_initial=" << M_initial << "  r=" << r << "\n";
        os << "  B0_G=" << B0_G << "  B_T=" << B_T() << "\n";
        os << "  spin_factor=" << spin_factor << "  M_dot_0=" << M_dot_0 << "\n";
        os << "  tau_acc=" << tau_acc << "  Hz=" << Hz << "\n";
    }

    double exampleAt1Gyr() const { return compute_g_SgrAStar(1.0e9 * 3.156e7); }
};

}  // namespace UQFF

#endif  // SMBH_SGR_A_STAR_H
"""


def main():
    print("gen_muge_sgra.py — Generating SMBHSgrAStar.h ...")
    content = get_content()
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  Written: {OUTPUT_FILE}  ({len(content.splitlines())} lines)")


if __name__ == "__main__":
    main()
