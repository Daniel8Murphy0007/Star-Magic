"""
gen_muge_ngc3603.py — Generator for NGC3603.h
Module 11: NGC 3603 extreme star cluster.
Unique physics: cavity pressure P(t)=P0*exp(-t/tau_exp)/rho_fluid,
                stellar wind term, M(t) with SFR growth.

Run:  python gen_muge_ngc3603.py
Output: NGC3603.h
"""
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "NGC3603.h")


def get_content():
    return r"""/**
 * ================================================================================================
 * Header: NGC3603.h
 *
 * Description: C++ Module for NGC 3603 Extreme Star Cluster Class (Module 11)
 *              UQFF simulations — young massive star cluster with cavity expansion.
 *
 * Unique Terms:
 *   - Cavity pressure: term_pressure = P(t) / rho_fluid
 *   - P(t) = P0 * exp(-t / tau_exp)  (expanding cavity dynamics)
 *   - Stellar wind: term_wind = rho_wind * v_wind^2 / rho_fluid
 *   - Mass growth: M(t) = M0 * (1 + M_dot_factor * exp(-t / tau_SF))
 *
 * Key Parameters:
 *   M0=400000 M_sun, r=9.5 ly, tau_SF=1 Myr, P0=4e-8 Pa,
 *   tau_exp=1 Myr, rho_wind=1e-20, v_wind=2e6 m/s
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef NGC_3603_H
#define NGC_3603_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

namespace UQFF {

class NGC3603 {
private:
    double G, M0, r;
    double B, B_crit;
    double Lambda, c_light, H0;
    double q_charge, gas_v, f_TRZ;
    double M_dot_factor, tau_SF;
    double rho_wind, v_wind;
    double rho_fluid;
    double P0;              // Initial cavity pressure (Pa)
    double tau_exp;         // Cavity expansion timescale (s)
    double rho_vac_UA, rho_vac_SCm, scale_EM, proton_mass;
    double hbar, t_Hubble, t_Hubble_gyr;
    double delta_x, delta_p, integral_psi;
    double A_osc, k_osc, omega_osc, x_pos;
    double M_DM_factor, delta_rho_over_rho;
    double ug1_base;

public:
    NGC3603() { initializeDefaults(); }
    ~NGC3603() {}

    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M0 = 400000.0 * M_sun;
        double ly = 9.461e15;
        r = 9.5 * ly;
        B = 1.0e-5;
        B_crit = 1.0e11;
        Lambda = 1.1e-52;
        c_light = 3.0e8;
        H0 = 2.184e-18;
        q_charge = 1.602e-19;
        gas_v = 1.0e5;
        f_TRZ = 0.1;
        M_dot_factor = 1.0;       // SFR growth factor
        tau_SF = 1.0e6 * 3.156e7; // 1 Myr
        rho_wind = 1.0e-20;
        v_wind = 2.0e6;
        rho_fluid = 1.0e-20;
        P0 = 4.0e-8;              // Pa
        tau_exp = 1.0e6 * 3.156e7;
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
        A_osc = 1.0e-10;
        k_osc = 1.0 / r;
        omega_osc = 2.0 * M_PI / (r / c_light);
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1.0e-5;

        updateCache();
    }

    void updateCache() { ug1_base = B_field * r * G * M0; }

    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G")           { G = newValue; }
        else if (varName == "M0")     { M0 = newValue; }
        else if (varName == "r")      { r = newValue; }
        else if (varName == "B")      { B = newValue; }
        else if (varName == "P0")     { P0 = newValue; }
        else if (varName == "tau_exp") { tau_exp = newValue; }
        else if (varName == "rho_wind") { rho_wind = newValue; }
        else if (varName == "v_wind") { v_wind = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
        else if (varName == "M_dot_factor") { M_dot_factor = newValue; }
        else if (varName == "tau_SF") { tau_SF = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "H0")    { H0 = newValue; }
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
        if (varName == "M0") return M0;
        if (varName == "r") return r;
        if (varName == "P0") return P0;
        if (varName == "tau_exp") return tau_exp;
        std::cerr << "Unknown variable '" << varName << "'.\n";
        return 0.0;
    }

    // M(t): with SFR growth
    double M_t(double t) const {
        return M0 * (1.0 + M_dot_factor * std::exp(-t / tau_SF));
    }

    // P(t): cavity pressure decay
    double P_t(double t) const { return P0 * std::exp(-t / tau_exp); }

    double compute_Ug(double Mt) const {
        // DPM-seeded: mu_s x grad(M_s/r) (not Newtonian GM/r^2)\ndouble ug1 = B_field * r * G * Mt;
        double corr_B = 1.0 - B / B_crit;
        return (ug1 + ug1 * corr_B) * (1.0 + f_TRZ);
    }

    double compute_V() const { return (4.0 / 3.0) * M_PI * r * r * r; }

    double compute_g_NGC3603(double t) const {
        if (t < 0.0) { std::cerr << "Error: t must be non-negative.\n"; return 0.0; }

        double Mt = M_t(t);
        double Pt = P_t(t);
        double ug1_t = B_field * r * G * Mt;

        // Term 1
        double corr_H = 1.0 + H0 * t;
        double corr_B = 1.0 - B / B_crit;
        double term1 = ug1_t * corr_H * corr_B;

        // UQFF Ug
        double term2 = compute_Ug(Mt);

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
        double term_fluid = (rho_fluid * V * ug1_t) / Mt;

        // Oscillatory
        double term_osc = 2.0 * A_osc * std::cos(k_osc * x_pos) * std::cos(omega_osc * t)
                        + (2.0 * M_PI / t_Hubble_gyr) * A_osc * std::cos(k_osc * x_pos - omega_osc * t);

        // DM
        double M_dm = Mt * M_DM_factor;
        double term_DM = ((Mt + M_dm) * (delta_rho_over_rho + 3.0 * G * Mt / (r * r * r))) / Mt;

        // Wind feedback
        double term_wind = (rho_wind * v_wind * v_wind) / rho_fluid;

        // Cavity pressure (unique to NGC3603)
        double term_pressure = Pt / rho_fluid;

        return term1 + term2 + term3 + term4 + term_q + term_fluid
             + term_osc + term_DM + term_wind + term_pressure;
    }

    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(6);
        os << "NGC 3603 Parameters:\n";
        os << "  M0=" << M0 << "  r=" << r << "\n";
        os << "  P0=" << P0 << "  tau_exp=" << tau_exp << "\n";
        os << "  rho_wind=" << rho_wind << "  v_wind=" << v_wind << "\n";
    }

    double exampleAt500kyr() const { return compute_g_NGC3603(5.0e5 * 3.156e7); }
};

}  // namespace UQFF

#endif  // NGC_3603_H
"""


def main():
    print("gen_muge_ngc3603.py — Generating NGC3603.h ...")
    content = get_content()
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  Written: {OUTPUT_FILE}  ({len(content.splitlines())} lines)")


if __name__ == "__main__":
    main()
