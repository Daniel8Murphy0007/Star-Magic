/**
 * ================================================================================================
 * Header: AntennaeGalaxies.h
 *
 * Description: C++ Module for Antennae Galaxies (NGC 4038/4039) Class (Module 13)
 *              UQFF simulations — galaxy merger interaction dynamics.
 *
 * Unique Terms:
 *   - Merger interaction: I(t) = I_0 * exp(-t / tau_merger)  (decaying interaction)
 *   - Applied as (1+I(t)) to BOTH term1 AND Ug  (boost during peak merger)
 *   - Star formation mass growth: M(t) = M0 * (1 + SFR_factor * exp(-t / tau_SF))
 *   - Hubble correction: H(z) with z=0.0105
 *
 * Key Parameters:
 *   M0=2e11 M_sun, r=30000 ly, z=0.0105, I_0=0.1, tau_merger=400 Myr
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * PAPER: PAPER_441 (Session 119, grok_share_68eb34022.txt, Daniel T. Murphy)
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef ANTENNAE_GALAXIES_H
#define ANTENNAE_GALAXIES_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

namespace UQFF {

class AntennaeGalaxies {
private:
    double G, M0, r;
    double B, B_crit;
    double Lambda, c_light, H0;
    double q_charge, gas_v, f_TRZ;
    double I_0, tau_merger;
    double SFR_factor, tau_SF;
    double z_gal;
    double rho_fluid;
    double rho_vac_UA, rho_vac_SCm, scale_EM, proton_mass;
    double hbar;
    double t_Hubble, t_Hubble_gyr;
    double delta_x, delta_p, integral_psi;
    double A_osc, k_osc, omega_osc, x_pos;
    double M_DM_factor, delta_rho_over_rho;
    double ug1_base;

public:
    AntennaeGalaxies() { initializeDefaults(); }
    ~AntennaeGalaxies() {}

    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M0 = 2.0e11 * M_sun;
        double ly = 9.461e15;
        r = 30000.0 * ly;
        B = 1.0e-10;
        B_crit = 1.0e11;
        Lambda = 1.1e-52;
        c_light = 3.0e8;
        H0 = 2.184e-18;
        q_charge = 1.602e-19;
        gas_v = 1.0e5;
        f_TRZ = 0.1;
        I_0 = 0.1;
        tau_merger = 400.0e6 * 3.156e7;   // 400 Myr
        SFR_factor = 20.0 / (2.0e11);     // 20 Msun/yr / M0
        tau_SF = 100.0e6 * 3.156e7;
        z_gal = 0.0105;
        rho_fluid = 1.0e-25;
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
        M_DM_factor = 0.5;
        delta_rho_over_rho = 1.0e-4;
        updateCache();
    }

    void updateCache() { ug1_base = (G * M0) / (r * r); }

    // Hz(z): approximate H(z) in s^-1
    double Hz() const {
        return H0 * std::sqrt(0.3 * std::pow(1.0 + z_gal, 3.0) + 0.7);
    }

    // Merger interaction I(t): decays from I_0
    double I_t(double t) const { return I_0 * std::exp(-t / tau_merger); }

    // M(t): SFR growth
    double M_t(double t) const { return M0 * (1.0 + SFR_factor * std::exp(-t / tau_SF)); }

    bool setVariable(const std::string& varName, double v) {
        if (varName == "M0")        { M0 = v; }
        else if (varName == "r")    { r = v; }
        else if (varName == "I_0")  { I_0 = v; }
        else if (varName == "tau_merger") { tau_merger = v; }
        else if (varName == "z_gal") { z_gal = v; }
        else if (varName == "SFR_factor") { SFR_factor = v; }
        else if (varName == "tau_SF") { tau_SF = v; }
        else if (varName == "B")    { B = v; }
        else if (varName == "H0")   { H0 = v; }
        else if (varName == "f_TRZ"){ f_TRZ = v; }
        else { std::cerr << "Unknown variable '" << varName << "'.\n"; return false; }
        updateCache();
        return true;
    }

    double getVariable(const std::string& varName) const {
        if (varName == "M0") return M0;
        if (varName == "r")  return r;
        if (varName == "I_0") return I_0;
        if (varName == "tau_merger") return tau_merger;
        std::cerr << "Unknown variable '" << varName << "'.\n";
        return 0.0;
    }

    double compute_Ug(double Mt, double It) const {
        // DPM-emergent: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        // DPM-emergent: gravity from magnetic moment x mass gradient (not Newtonian GM/r^2)
        double ug1 = (G * Mt) / (r * r);
        double corr_B = 1.0 - B / B_crit;
        // (1+I) boosts Ug during merger
        return (ug1 + ug1 * corr_B) * (1.0 + f_TRZ) * (1.0 + It);
    }

    double compute_V() const { return (4.0 / 3.0) * M_PI * r * r * r; }

    double compute_g_Antennae(double t) const {
        if (t < 0.0) { std::cerr << "Error: t must be non-negative.\n"; return 0.0; }

        double Mt = M_t(t);
        double It = I_t(t);
        double ug1_t = (G * Mt) / (r * r);
        double hz = Hz();

        // Term 1 with H(z) and B correction, boosted by (1+I)
        double corr_H = 1.0 + hz * t;
        double corr_B = 1.0 - B / B_crit;
        double term1 = ug1_t * corr_H * corr_B * (1.0 + It);

        // UQFF Ug with (1+I) merger boost
        double term2 = compute_Ug(Mt, It);

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
        os << "Antennae Galaxies Parameters:\n";
        os << "  M0=" << M0 << "  r=" << r << "  z=" << z_gal << "\n";
        os << "  I_0=" << I_0 << "  tau_merger=" << tau_merger << "\n";
    }

    double exampleAt400Myr() const { return compute_g_Antennae(400.0e6 * 3.156e7); }
};

}  // namespace UQFF

#endif  // ANTENNAE_GALAXIES_H
