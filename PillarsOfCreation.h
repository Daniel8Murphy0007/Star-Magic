/**
 * ================================================================================================
 * Header: PillarsOfCreation.h
 *
 * Description: C++ Module for Eagle Nebula Pillars of Creation Class (Module 8)
 *              UQFF simulations — pillar erosion by photoionization and radiation pressure.
 *
 * Unique Terms:
 *   - E(t) = E_0 * exp(-t / tau_erosion)  [decaying erosion function]
 *   - Correction: (1 - E(t)) applied to term1 AND Ug (reduces gravity as pillars erode)
 *   - This is DECAYING erosion (E decreases over time as pillars are destroyed)
 *
 * Key Parameters:
 *   M=10100 M_sun, r=5 ly, B=1e-6 T, E_0=0.1, tau_erosion=1 Myr
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef PILLARS_OF_CREATION_H
#define PILLARS_OF_CREATION_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

namespace UQFF {

class PillarsOfCreation {
private:
    double G, M, r;
    double B, B_crit;
    double Lambda, c_light, Hz;
    double q_charge, gas_v, f_TRZ;
    double E_0;             // Initial erosion factor (dimensionless)
    double tau_erosion;     // Erosion timescale (s)
    double rho_fluid;
    double rho_vac_UA, rho_vac_SCm, scale_EM, proton_mass;
    double hbar, t_Hubble, t_Hubble_gyr;
    double delta_x, delta_p, integral_psi;
    double A_osc, k_osc, omega_osc, x_pos;
    double M_DM_factor, delta_rho_over_rho;
    double ug1_base;

public:
    PillarsOfCreation() { initializeDefaults(); }
    ~PillarsOfCreation() {}

    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M = 10100.0 * M_sun;
        double ly = 9.461e15;
        r = 5.0 * ly;
        B = 1.0e-6;
        B_crit = 1.0e11;
        Lambda = 1.1e-52;
        c_light = 3.0e8;
        Hz = 2.184e-18;
        q_charge = 1.602e-19;
        gas_v = 1.0e5;
        f_TRZ = 0.1;
        E_0 = 0.1;
        tau_erosion = 1.0e6 * 3.156e7;  // 1 Myr
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

    void updateCache() { ug1_base = (G * M) / (r * r); }

    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G")            { G = newValue; }
        else if (varName == "M")       { M = newValue; }
        else if (varName == "r")       { r = newValue; }
        else if (varName == "B")       { B = newValue; }
        else if (varName == "E_0")     { E_0 = newValue; }
        else if (varName == "tau_erosion") { tau_erosion = newValue; }
        else if (varName == "f_TRZ")   { f_TRZ = newValue; }
        else if (varName == "Hz")      { Hz = newValue; }
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
        if (varName == "M") return M;
        if (varName == "r") return r;
        if (varName == "E_0") return E_0;
        if (varName == "tau_erosion") return tau_erosion;
        std::cerr << "Unknown variable '" << varName << "'.\n";
        return 0.0;
    }

    // E(t): decaying erosion function (UNIQUE TO PILLARS / HORSEHEAD)
    double E_t(double t) const { return E_0 * std::exp(-t / tau_erosion); }

    double compute_Ug(double Et) const {
        double Ug1 = ug1_base;
        double corr_B = 1.0 - B / B_crit;
        double Ug4 = Ug1 * corr_B;
        // Apply (1-E(t)) erosion correction to Ug
        return (Ug1 + Ug4) * (1.0 + f_TRZ) * (1.0 - Et);
    }

    double compute_V() const { return (4.0 / 3.0) * M_PI * r * r * r; }

    double compute_g_Pillars(double t) const {
        if (t < 0.0) { std::cerr << "Error: t must be non-negative.\n"; return 0.0; }

        double Et = E_t(t);

        // Term 1: Base + Hz + B, with erosion (1-E(t)) correction
        double corr_H = 1.0 + Hz * t;
        double corr_B = 1.0 - B / B_crit;
        double term1 = ug1_base * corr_H * corr_B * (1.0 - Et);

        // UQFF Ug with erosion
        double term2 = compute_Ug(Et);

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
        os << "Pillars of Creation Parameters:\n";
        os << "  M=" << M << "  r=" << r << "  B=" << B << "\n";
        os << "  E_0=" << E_0 << "  tau_erosion=" << tau_erosion << "\n";
    }

    double exampleAt1Myr() const { return compute_g_Pillars(1.0e6 * 3.156e7); }
};

}  // namespace UQFF

#endif  // PILLARS_OF_CREATION_H
