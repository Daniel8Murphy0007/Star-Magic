/**
 * ================================================================================================
 * Header: StarbirthTapestry.h
 *
 * Description: C++ Module for Tapestry of Blazing Starbirth (NGC 2014/2020, LMC) Class (Module 6)
 *              UQFF simulations — starbirth region gravity with wind feedback and SFR evolution.
 *
 * Unique Terms:
 *   - Wind feedback: term_wind = rho_wind*v_wind^2 / rho_fluid
 *   - M_dot_factor = gas_mass / M_initial (dimensionless SFR ratio)
 *   - M(t) = M_initial*(1 + M_dot_factor*exp(-t/tau_SF))
 *
 * Key Parameters:
 *   M=240 M_sun, r=10 ly, B=1e-6 T, rho_wind=1e-21 kg/m^3,
 *   v_wind=2e6 m/s, tau_SF=5 Myr, gas_mass=2400 M_sun
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy
 * ================================================================================================
 */

#ifndef STARBIRTH_TAPESTRY_H
#define STARBIRTH_TAPESTRY_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

namespace UQFF {

class StarbirthTapestry {
private:
    double G, M_initial, r;
    double B, B_crit;
    double Lambda, c_light, Hz;
    double q_charge, gas_v, f_TRZ;
    double rho_wind;        // Wind density (kg/m^3)
    double v_wind;          // Wind velocity (m/s)
    double rho_fluid;       // Ambient fluid density
    double M_dot_factor;    // Dimensionless SFR ratio (gas_mass/M_initial)
    double tau_SF;          // Star formation timescale (s)
    double rho_vac_UA, rho_vac_SCm, scale_EM, proton_mass;
    // Full-term
    double hbar, t_Hubble, t_Hubble_gyr;
    double delta_x, delta_p, integral_psi;
    double A_osc, k_osc, omega_osc, x_pos;
    double M_DM_factor, delta_rho_over_rho;
    // Cached
    double ug1_base;

public:
    StarbirthTapestry() { initializeDefaults(); }
    ~StarbirthTapestry() {}

    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M_initial = 240.0 * M_sun;
        double ly = 9.461e15;
        r = 10.0 * ly;
        B = 1.0e-6;
        B_crit = 1.0e11;
        Lambda = 1.1e-52;
        c_light = 3.0e8;
        Hz = 2.184e-18;
        q_charge = 1.602e-19;
        gas_v = 1.0e5;
        f_TRZ = 0.1;
        rho_wind = 1.0e-21;
        v_wind = 2.0e6;
        rho_fluid = 1.0e-21;
        double gas_mass = 2400.0 * M_sun;
        M_dot_factor = gas_mass / M_initial;  // = 10.0
        tau_SF = 5.0e6 * 3.156e7;
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

    void updateCache() { ug1_base = (G * M_initial) / (r * r); }

    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G")           { G = newValue; }
        else if (varName == "M_initial") { M_initial = newValue; }
        else if (varName == "r")      { r = newValue; }
        else if (varName == "B")      { B = newValue; }
        else if (varName == "rho_wind") { rho_wind = newValue; }
        else if (varName == "v_wind") { v_wind = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
        else if (varName == "M_dot_factor") { M_dot_factor = newValue; }
        else if (varName == "tau_SF") { tau_SF = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "Hz")    { Hz = newValue; }
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
        if (varName == "rho_wind") return rho_wind;
        if (varName == "v_wind") return v_wind;
        if (varName == "tau_SF") return tau_SF;
        std::cerr << "Unknown variable '" << varName << "'.\n";
        return 0.0;
    }

    // M(t): mass with SFR growth
    double M_t(double t) const {
        return M_initial * (1.0 + M_dot_factor * std::exp(-t / tau_SF));
    }

    double compute_Ug(double Mt) const {
        double ug1 = (G * Mt) / (r * r);
        double corr_B = 1.0 - B / B_crit;
        return (ug1 + ug1 * corr_B) * (1.0 + f_TRZ);
    }

    double compute_V() const { return (4.0 / 3.0) * M_PI * r * r * r; }

    double compute_g_Tapestry(double t) const {
        if (t < 0.0) { std::cerr << "Error: t must be non-negative.\n"; return 0.0; }

        double Mt = M_t(t);
        double ug1_t = (G * Mt) / (r * r);

        // Term 1: Base + Hz + B
        double term1 = ug1_t * (1.0 + Hz * t) * (1.0 - B / B_crit);

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

        // Wind feedback term (unique to Tapestry/WD2)
        double term_wind = (rho_wind * v_wind * v_wind) / rho_fluid;

        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind;
    }

    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(6);
        os << "Starbirth Tapestry (NGC 2014/2020) Parameters:\n";
        os << "  M_initial=" << M_initial << "  r=" << r << "\n";
        os << "  B=" << B << "  rho_wind=" << rho_wind << "  v_wind=" << v_wind << "\n";
        os << "  M_dot_factor=" << M_dot_factor << "  tau_SF=" << tau_SF << "\n";
    }

    double exampleAt5Myr() const { return compute_g_Tapestry(5.0e6 * 3.156e7); }
};

}  // namespace UQFF

#endif  // STARBIRTH_TAPESTRY_H
