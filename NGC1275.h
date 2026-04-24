/**
 * ================================================================================================
 * Header: NGC1275.h
 *
 * Description: C++ Module for NGC 1275 (Perseus A) Galaxy Cluster BCG Class (Module 15)
 *              UQFF simulations — BCG with cooling flows, filaments, and AGN feedback.
 *
 * Unique Terms:
 *   - B(t) = B0 * exp(-t / tau_B)            (decaying AGN magnetic field)
 *   - F(t) = F0 * exp(-t / tau_fil)           (decaying filament factor)
 *     * Applied as (1+F(t)) multiplier on BOTH term1 AND Ug
 *   - Cooling flow: term_cool = rho_cool * v_cool^2 / rho_fluid
 *   - Central BH: g_BH = G * M_BH / r_BH^2
 *   - H(z) = H0 * sqrt(0.3*(1+z)^3 + 0.7) * 1000/3.086e22
 *
 * Key Parameters:
 *   M=1e11 M_sun, r=200000 ly, M_BH=8e8 M_sun, z=0.0176
 *   B0=5e-9 T, tau_B=100 Myr, F0=0.1, tau_fil=100 Myr
 *   rho_cool=1e-20, v_cool=3e3 m/s
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * PAPER: PAPER_443 (Session 119, grok_share_68eb34022.txt, Daniel T. Murphy)
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef NGC_1275_H
#define NGC_1275_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

namespace UQFF {

class NGC1275 {
private:
    double G, M0, r;
    double M_BH, r_BH;
    double B0, B_crit, tau_B;
    double F0, tau_fil;
    double Lambda, c_light, H0;
    double z_gal;
    double q_charge, gas_v, f_TRZ;
    double rho_cool, v_cool, rho_fluid;
    double rho_vac_UA, rho_vac_SCm, scale_EM, proton_mass;
    double hbar;
    double t_Hubble, t_Hubble_gyr;
    double delta_x, delta_p, integral_psi;
    double A_osc, k_osc, omega_osc, x_pos;
    double M_DM_factor, delta_rho_over_rho;
    double ug1_base;

public:
    NGC1275() { initializeDefaults(); }
    ~NGC1275() {}

    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        double ly = 9.461e15;
        double kpc = 3.086e19;
        M0 = 1.0e11 * M_sun;
        r = 200000.0 * ly;
        M_BH = 8.0e8 * M_sun;
        r_BH = 0.1 * kpc;
        B0 = 5.0e-9;
        B_crit = 1.0e11;
        tau_B = 100.0e6 * 3.156e7;   // 100 Myr
        F0 = 0.1;
        tau_fil = 100.0e6 * 3.156e7; // 100 Myr
        Lambda = 1.1e-52;
        c_light = 3.0e8;
        H0 = 2.184e-18;
        z_gal = 0.0176;
        q_charge = 1.602e-19;
        gas_v = 1.0e5;
        f_TRZ = 0.1;
        rho_cool = 1.0e-20;
        v_cool = 3.0e3;
        rho_fluid = 1.0e-23;
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
        M_DM_factor = 5.0;
        delta_rho_over_rho = 1.0e-3;
        updateCache();
    }

    void updateCache() { ug1_base = (G * M0) / (r * r); }

    // H(z): full Hubble parameter
    double Hz() const {
        return H0 * std::sqrt(0.3 * std::pow(1.0 + z_gal, 3.0) + 0.7);
    }

    // B(t): decaying AGN magnetic field
    double B_t(double t) const { return B0 * std::exp(-t / tau_B); }

    // F(t): decaying filament contribution factor
    double F_t(double t) const { return F0 * std::exp(-t / tau_fil); }

    bool setVariable(const std::string& varName, double v) {
        if (varName == "M0")         { M0 = v; }
        else if (varName == "r")     { r = v; }
        else if (varName == "M_BH")  { M_BH = v; }
        else if (varName == "r_BH")  { r_BH = v; }
        else if (varName == "B0")    { B0 = v; }
        else if (varName == "tau_B") { tau_B = v; }
        else if (varName == "F0")    { F0 = v; }
        else if (varName == "tau_fil") { tau_fil = v; }
        else if (varName == "z_gal") { z_gal = v; }
        else if (varName == "rho_cool") { rho_cool = v; }
        else if (varName == "v_cool")   { v_cool = v; }
        else if (varName == "rho_fluid") { rho_fluid = v; }
        else if (varName == "H0")    { H0 = v; }
        else if (varName == "f_TRZ") { f_TRZ = v; }
        else { std::cerr << "Unknown variable '" << varName << "'.\n"; return false; }
        updateCache();
        return true;
    }

    double getVariable(const std::string& varName) const {
        if (varName == "M0")  return M0;
        if (varName == "r")   return r;
        if (varName == "M_BH") return M_BH;
        if (varName == "B0")  return B0;
        if (varName == "F0")  return F0;
        if (varName == "z_gal") return z_gal;
        std::cerr << "Unknown '" << varName << "'.\n";
        return 0.0;
    }

    double compute_Ug(double Bt, double Ft) const {
        // DPM-seeded: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        // DPM-seeded: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        double ug1 = (G * M0) / (r * r);
        double corr_B = 1.0 - Bt / B_crit;
        // (1+F) filament enhancement
        return (ug1 + ug1 * corr_B) * (1.0 + f_TRZ) * (1.0 + Ft);
    }

    double compute_V() const { return (4.0 / 3.0) * M_PI * r * r * r; }

    double compute_g_NGC1275(double t) const {
        if (t < 0.0) { std::cerr << "Error: t must be non-negative.\n"; return 0.0; }

        double Bt = B_t(t);
        double Ft = F_t(t);
        double ug1_t = (G * M0) / (r * r);
        double hz = Hz();

        // Term 1: base gravity with H(z), B(t) correction, and (1+F) filament factor
        double corr_H = 1.0 + hz * t;
        double corr_B = 1.0 - Bt / B_crit;
        double term1 = ug1_t * corr_H * corr_B * (1.0 + Ft);

        // UQFF Ug with filament enhancement
        double term2 = compute_Ug(Bt, Ft);

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
        double term_fluid = (rho_fluid * V * ug1_t) / M0;

        // Oscillatory
        double term_osc = 2.0 * A_osc * std::cos(k_osc * x_pos) * std::cos(omega_osc * t)
                        + (2.0 * M_PI / t_Hubble_gyr) * A_osc * std::cos(k_osc * x_pos - omega_osc * t);

        // DM
        double M_dm = M0 * M_DM_factor;
        double term_DM = ((M0 + M_dm) * (delta_rho_over_rho + 3.0 * G * M0 / (r * r * r))) / M0;

        // Cooling flow (unique to NGC1275)
        double term_cool = (rho_cool * v_cool * v_cool) / rho_fluid;

        // Central black hole gravity
        double g_BH = (G * M_BH) / (r_BH * r_BH);

        return term1 + term2 + term3 + term4 + term_q + term_fluid
             + term_osc + term_DM + term_cool + g_BH;
    }

    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(6);
        os << "NGC 1275 (Perseus A) Parameters:\n";
        os << "  M0=" << M0 << "  r=" << r << "  z=" << z_gal << "\n";
        os << "  M_BH=" << M_BH << "  r_BH=" << r_BH << "\n";
        os << "  B0=" << B0 << "  tau_B=" << tau_B << "\n";
        os << "  F0=" << F0 << "  tau_fil=" << tau_fil << "\n";
        os << "  rho_cool=" << rho_cool << "  v_cool=" << v_cool << "\n";
    }

    double exampleAt50Myr() const { return compute_g_NGC1275(50.0e6 * 3.156e7); }
};

}  // namespace UQFF

#endif  // NGC_1275_H
