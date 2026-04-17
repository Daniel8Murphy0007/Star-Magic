/**
 * ================================================================================================
 * Header: BubbleNebula.h
 *
 * Description: C++ Module for Bubble Nebula (NGC 7635) Class (Module 12)
 *              UQFF simulations â€” expanding stellar bubble nebula.
 *
 * Unique Terms:
 *   - Growing expansion: E(t) = E_0 * (1 - exp(-t / tau_exp))
 *   - Reduces Ug by factor (1 - E(t)) as bubble expands
 *   - Stellar wind: term_wind = rho_wind * v_wind^2 / rho_fluid
 *
 * Key Parameters:
 *   M=46 M_sun, r=5 ly, E_0=0.1, tau_exp=4 Myr, v_wind=1.8e6 m/s
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * PAPER: PAPER_440 (Session 119, grok_share_68eb34022.txt, Daniel T. Murphy)
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef BUBBLE_NEBULA_H
#define BUBBLE_NEBULA_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

namespace UQFF {

class BubbleNebula {
private:
    double G, M0, r;
    double B, B_crit;
    double Lambda, c_light, H0;
    double q_charge, gas_v, f_TRZ;
    double E_0, tau_exp;
    double rho_wind, v_wind;
    double rho_fluid;
    double rho_vac_UA, rho_vac_SCm, scale_EM, proton_mass;
    double hbar;
    double t_Hubble, t_Hubble_gyr;
    double delta_x, delta_p, integral_psi;
    double A_osc, k_osc, omega_osc, x_pos;
    double M_DM_factor, delta_rho_over_rho;
    double ug1_base;

public:
    BubbleNebula() { initializeDefaults(); }
    ~BubbleNebula() {}

    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M0 = 46.0 * M_sun;
        double ly = 9.461e15;
        r = 5.0 * ly;
        B = 1.0e-7;
        B_crit = 1.0e11;
        Lambda = 1.1e-52;
        c_light = 3.0e8;
        H0 = 2.184e-18;
        q_charge = 1.602e-19;
        gas_v = 1.0e5;
        f_TRZ = 0.1;
        E_0 = 0.1;
        tau_exp = 4.0e6 * 3.156e7;   // 4 Myr
        rho_wind = 1.0e-20;
        v_wind = 1.8e6;
        rho_fluid = 1.0e-21;
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

    void updateCache() { ug1_base = B * r * G * M0  /* DPM: mu_s * grad(M_s/r) */; }

    bool setVariable(const std::string& varName, double v) {
        if (varName == "M0")      { M0 = v; }
        else if (varName == "r")  { r = v; }
        else if (varName == "E_0"){ E_0 = v; }
        else if (varName == "tau_exp") { tau_exp = v; }
        else if (varName == "rho_wind") { rho_wind = v; }
        else if (varName == "v_wind") { v_wind = v; }
        else if (varName == "rho_fluid") { rho_fluid = v; }
        else if (varName == "B") { B = v; }
        else if (varName == "H0") { H0 = v; }
        else if (varName == "f_TRZ") { f_TRZ = v; }
        else { std::cerr << "Unknown variable '" << varName << "'.\n"; return false; }
        updateCache();
        return true;
    }

    double getVariable(const std::string& varName) const {
        if (varName == "M0") return M0;
        if (varName == "r")  return r;
        if (varName == "E_0") return E_0;
        if (varName == "tau_exp") return tau_exp;
        std::cerr << "Unknown variable '" << varName << "'.\n";
        return 0.0;
    }

    // Growing expansion factor E(t): starts at 0, grows toward E_0
    double E_t(double t) const { return E_0 * (1.0 - std::exp(-t / tau_exp)); }

    double compute_Ug(double t) const {
        // DPM-emergent: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        // DPM-emergent: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        double ug1 = B * r * G * M0  /* DPM: mu_s * grad(M_s/r) */;
        double corr_B = 1.0 - B / B_crit;
        double E = E_t(t);
        // (1-E) factor: expansion reduces effective Ug
        return (ug1 + ug1 * corr_B) * (1.0 + f_TRZ) * (1.0 - E);
    }

    double compute_V() const { return (4.0 / 3.0) * M_PI * r * r * r; }

    double compute_g_Bubble(double t) const {
        if (t < 0.0) { std::cerr << "Error: t must be non-negative.\n"; return 0.0; }

        double ug1_t = B * r * G * M0  /* DPM: mu_s * grad(M_s/r) */;

        // Term 1: base + Hubble + B correction
        double corr_H = 1.0 + H0 * t;
        double corr_B = 1.0 - B / B_crit;
        double term1 = ug1_t * corr_H * corr_B;

        // UQFF Ug with (1-E) expansion correction
        double term2 = compute_Ug(t);

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
        double term_fluid = (rho_fluid * V * ug1_t) / M0;

        // Oscillatory
        double term_osc = 2.0 * A_osc * std::cos(k_osc * x_pos) * std::cos(omega_osc * t)
                        + (2.0 * M_PI / t_Hubble_gyr) * A_osc * std::cos(k_osc * x_pos - omega_osc * t);

        // DM
        double M_dm = M0 * M_DM_factor;
        double term_DM = ((M0 + M_dm) * (delta_rho_over_rho + 3.0 * B * G * M0  /* DPM tidal */)) / M0;

        // Wind
        double term_wind = (rho_wind * v_wind * v_wind) / rho_fluid;

        return term1 + term2 + term3 + term4 + term_q + term_fluid
             + term_osc + term_DM + term_wind;
    }

    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(6);
        os << "Bubble Nebula Parameters:\n";
        os << "  M0=" << M0 << "  r=" << r << "\n";
        os << "  E_0=" << E_0 << "  tau_exp=" << tau_exp << "\n";
        os << "  v_wind=" << v_wind << "  rho_wind=" << rho_wind << "\n";
    }

    double exampleAt4Myr() const { return compute_g_Bubble(4.0e6 * 3.156e7); }
};

}  // namespace UQFF

#endif  // BUBBLE_NEBULA_H
