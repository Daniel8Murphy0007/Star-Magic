/**
 * ================================================================================================
 * Header: HUDFGalaxies.h
 *
 * Description: C++ Module for Hubble Ultra Deep Field (HUDF) Galaxy Ensemble Class (Module 16)
 *              UQFF simulations â€” cosmic-scale deep field galaxy population.
 *
 * Unique Terms:
 *   - H(z) at z_avg=3.5: full cosmological Hubble expansion correction
 *   - H(z) = H0 * sqrt(0.3*(1+z)^3 + 0.7) * 1000/3.086e22
 *   - Ensemble mass: M0 = 1e12 M_sun (average across ~10,000 galaxies)
 *   - Cosmic comoving distance r at z~3.5 (~30 Gly)
 *
 * Key Parameters:
 *   M0=1e12 M_sun, r=1.3e11 ly (comoving), z_avg=3.5, Hz~2.5e-18 s^-1
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * PAPER: PAPER_444 (Session 119, grok_share_68eb34022.txt, Daniel T. Murphy)
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef HUDF_GALAXIES_H
#define HUDF_GALAXIES_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

namespace UQFF {

class HUDFGalaxies {
private:
    double G, M0, r;
    double B, B_crit;
    double Lambda, c_light, H0;
    double z_avg;
    double q_charge, gas_v, f_TRZ;
    double rho_fluid;
    double rho_vac_UA, rho_vac_SCm, scale_EM, proton_mass;
    double hbar;
    double t_Hubble, t_Hubble_gyr;
    double delta_x, delta_p, integral_psi;
    double A_osc, k_osc, omega_osc, x_pos;
    double M_DM_factor, delta_rho_over_rho;
    double ug1_base;

public:
    HUDFGalaxies() { initializeDefaults(); }
    ~HUDFGalaxies() {}

    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        double ly = 9.461e15;
        M0 = 1.0e12 * M_sun;                // Ensemble mass
        r = 1.3e11 * ly;                    // ~30 Gly comoving
        B = 1.0e-10;
        B_crit = 1.0e11;
        Lambda = 1.1e-52;
        c_light = 3.0e8;
        H0 = 2.184e-18;
        z_avg = 3.5;
        q_charge = 1.602e-19;
        gas_v = 1.0e5;
        f_TRZ = 0.1;
        rho_fluid = 1.0e-30;                // Cosmic mean density
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
        M_DM_factor = 5.0;                  // High DM fraction at z~3.5
        delta_rho_over_rho = 1.0e-3;
        updateCache();
    }

    void updateCache() { ug1_base = B * r * G * M0  /* DPM: mu_s * grad(M_s/r) */; }

    // H(z) at z_avg â€” full cosmological Hubble parameter (s^-1)
    double Hz() const {
        return H0 * std::sqrt(0.3 * std::pow(1.0 + z_avg, 3.0) + 0.7);
    }

    bool setVariable(const std::string& varName, double v) {
        if (varName == "M0")        { M0 = v; }
        else if (varName == "r")    { r = v; }
        else if (varName == "z_avg") { z_avg = v; }
        else if (varName == "B")    { B = v; }
        else if (varName == "H0")   { H0 = v; }
        else if (varName == "f_TRZ") { f_TRZ = v; }
        else if (varName == "rho_fluid") { rho_fluid = v; }
        else if (varName == "M_DM_factor") { M_DM_factor = v; }
        else { std::cerr << "Unknown variable '" << varName << "'.\n"; return false; }
        updateCache();
        return true;
    }

    double getVariable(const std::string& varName) const {
        if (varName == "M0")  return M0;
        if (varName == "r")   return r;
        if (varName == "z_avg") return z_avg;
        if (varName == "Hz")  return Hz();
        std::cerr << "Unknown '" << varName << "'.\n";
        return 0.0;
    }

    double compute_Ug() const {
        // DPM-seeded: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        // DPM-seeded: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        double ug1 = B * r * G * M0  /* DPM: mu_s * grad(M_s/r) */;
        double corr_B = 1.0 - B / B_crit;
        return (ug1 + ug1 * corr_B) * (1.0 + f_TRZ);
    }

    double compute_V() const { return (4.0 / 3.0) * M_PI * r * r * r; }

    double compute_g_HUDF(double t) const {
        if (t < 0.0) { std::cerr << "Error: t must be non-negative.\n"; return 0.0; }

        double hz = Hz();
        double ug1_t = B * r * G * M0  /* DPM: mu_s * grad(M_s/r) */;

        // Term 1: full H(z) Hubble expansion at z_avg=3.5
        double corr_H = 1.0 + hz * t;
        double corr_B = 1.0 - B / B_crit;
        double term1 = ug1_t * corr_H * corr_B;

        // UQFF Ug
        double term2 = compute_Ug();

        // Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // EM (weak at cosmic scales)
        double cross_vB = gas_v * B;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1.0 + (rho_vac_UA / rho_vac_SCm);
        double term4 = em_base * corr_UA * scale_EM;

        // Quantum
        double sqrt_unc = std::sqrt(delta_x * delta_p);
        double term_q = (hbar / sqrt_unc) * integral_psi * (2.0 * M_PI / t_Hubble);

        // Fluid
        double V = compute_V();
        double term_fluid = (rho_fluid * V * ug1_t) / M0;

        // Oscillatory
        double term_osc = 2.0 * A_osc * std::cos(k_osc * x_pos) * std::cos(omega_osc * t)
                        + (2.0 * M_PI / t_Hubble_gyr) * A_osc * std::cos(k_osc * x_pos - omega_osc * t);

        // DM (large at z~3.5)
        double M_dm = M0 * M_DM_factor;
        double term_DM = ((M0 + M_dm) * (delta_rho_over_rho + 3.0 * B * G * M0  /* DPM tidal */)) / M0;

        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM;
    }

    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(6);
        os << "HUDF Galaxies Parameters:\n";
        os << "  M0=" << M0 << "  r=" << r << "\n";
        os << "  z_avg=" << z_avg << "  Hz=" << Hz() << " s^-1\n";
    }

    double exampleAt1Gyr() const { return compute_g_HUDF(1.0e9 * 3.156e7); }
};

}  // namespace UQFF

#endif  // HUDF_GALAXIES_H
