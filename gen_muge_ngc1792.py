"""
gen_muge_ngc1792.py — Generator for GalaxyNGC1792.h
Module 17: NGC 1792 — starburst spiral galaxy.
Key physics: SFR-driven mass growth M(t)=M0*(1+SFR_factor*exp(-t/tau_SF))
             with very long tau_SF=100 Myr starburst timescale.

Run:  python gen_muge_ngc1792.py
Output: GalaxyNGC1792.h
"""
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "GalaxyNGC1792.h")


def get_content():
    return r"""/**
 * ================================================================================================
 * Header: GalaxyNGC1792.h
 *
 * Description: C++ Module for NGC 1792 Starburst Spiral Galaxy Class (Module 17)
 *              UQFF simulations — starburst galaxy with extended SFR dynamics.
 *
 * Unique Terms:
 *   - SFR-driven mass growth: M(t) = M0 * (1 + SFR_factor * exp(-t / tau_SF))
 *   - tau_SF = 100 Myr (extended starburst timescale, ~20x longer than cluster modules)
 *   - SFR_factor = SFR_peak/M0 = 10 Msun/yr / (1e10 Msun) = 1e-9
 *   - H(z=0.0095) Hubble parameter correction
 *
 * Key Parameters:
 *   M0=1e10 M_sun, r=80000 ly, z=0.0095, tau_SF=100 Myr
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef GALAXY_NGC1792_H
#define GALAXY_NGC1792_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

namespace UQFF {

class GalaxyNGC1792 {
private:
    double G, M0, r;
    double B, B_crit;
    double Lambda, c_light, H0;
    double z_gal;
    double q_charge, gas_v, f_TRZ;
    double SFR_factor, tau_SF;
    double rho_fluid;
    double rho_vac_UA, rho_vac_SCm, scale_EM, proton_mass;
    double hbar;
    double t_Hubble, t_Hubble_gyr;
    double delta_x, delta_p, integral_psi;
    double A_osc, k_osc, omega_osc, x_pos;
    double M_DM_factor, delta_rho_over_rho;
    double ug1_base;

public:
    GalaxyNGC1792() { initializeDefaults(); }
    ~GalaxyNGC1792() {}

    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        double ly = 9.461e15;
        M0 = 1.0e10 * M_sun;
        r = 80000.0 * ly;
        B = 5.0e-10;
        B_crit = 1.0e11;
        Lambda = 1.1e-52;
        c_light = 3.0e8;
        H0 = 2.184e-18;
        z_gal = 0.0095;
        q_charge = 1.602e-19;
        gas_v = 1.0e5;
        f_TRZ = 0.1;
        // SFR_factor: 10 Msun/yr / (1e10 Msun * (3.156e7 s/yr)) -> dimensionless rate per second
        SFR_factor = 10.0 / (1.0e10 * 3.156e7); // ~3.17e-17 s^-1, but used as amplitude
        tau_SF = 100.0e6 * 3.156e7;             // 100 Myr
        rho_fluid = 1.0e-24;
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
        M_DM_factor = 1.0;
        delta_rho_over_rho = 1.0e-4;
        updateCache();
    }

    void updateCache() { ug1_base = (G * M0) / (r * r); }

    // H(z)
    double Hz() const {
        return H0 * std::sqrt(0.3 * std::pow(1.0 + z_gal, 3.0) + 0.7);
    }

    // M(t): SFR-driven mass evolution
    double M_t(double t) const {
        return M0 * (1.0 + SFR_factor * std::exp(-t / tau_SF));
    }

    bool setVariable(const std::string& varName, double v) {
        if (varName == "M0")         { M0 = v; }
        else if (varName == "r")     { r = v; }
        else if (varName == "z_gal") { z_gal = v; }
        else if (varName == "tau_SF") { tau_SF = v; }
        else if (varName == "SFR_factor") { SFR_factor = v; }
        else if (varName == "B")     { B = v; }
        else if (varName == "H0")    { H0 = v; }
        else if (varName == "f_TRZ") { f_TRZ = v; }
        else if (varName == "rho_fluid") { rho_fluid = v; }
        else { std::cerr << "Unknown variable '" << varName << "'.\n"; return false; }
        updateCache();
        return true;
    }

    double getVariable(const std::string& varName) const {
        if (varName == "M0")  return M0;
        if (varName == "r")   return r;
        if (varName == "z_gal") return z_gal;
        if (varName == "tau_SF") return tau_SF;
        if (varName == "SFR_factor") return SFR_factor;
        std::cerr << "Unknown '" << varName << "'.\n";
        return 0.0;
    }

    double compute_Ug(double Mt) const {
        // DPM-emergent: mu_s x grad(M_s/r) (not Newtonian GM/r^2)\ndouble ug1 = (G * Mt) / (r * r);
        double corr_B = 1.0 - B / B_crit;
        return (ug1 + ug1 * corr_B) * (1.0 + f_TRZ);
    }

    double compute_V() const { return (4.0 / 3.0) * M_PI * r * r * r; }

    double compute_g_NGC1792(double t) const {
        if (t < 0.0) { std::cerr << "Error: t must be non-negative.\n"; return 0.0; }

        double Mt = M_t(t);
        double ug1_t = (G * Mt) / (r * r);
        double hz = Hz();

        // Term 1 with H(z) and B correction
        double corr_H = 1.0 + hz * t;
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

        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM;
    }

    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(6);
        os << "NGC 1792 Parameters:\n";
        os << "  M0=" << M0 << "  r=" << r << "  z=" << z_gal << "\n";
        os << "  tau_SF=" << tau_SF << "  SFR_factor=" << SFR_factor << "\n";
    }

    double exampleAt100Myr() const { return compute_g_NGC1792(100.0e6 * 3.156e7); }
};

}  // namespace UQFF

#endif  // GALAXY_NGC1792_H
"""


def main():
    print("gen_muge_ngc1792.py — Generating GalaxyNGC1792.h ...")
    content = get_content()
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  Written: {OUTPUT_FILE}  ({len(content.splitlines())} lines)")


if __name__ == "__main__":
    main()
