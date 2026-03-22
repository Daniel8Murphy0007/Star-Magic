// YoungStarsOutflowsUQFFModule.cpp
// Doc 35: NGC 346-analogue — P_outflow = ρ v_out² (1 + t/t_evolve).
//
// Key new term:
//   P_outflow(t) = rho * v_out^2 * (1 + t/t_evolve)
//   This replaces the simple wind term with a time-growing outflow pressure.
//
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#include "YoungStarsOutflowsUQFFModule.h"
#include "UQFFConstants.h"
#include <cmath>

YoungStarsOutflowsUQFFModule::YoungStarsOutflowsUQFFModule() {
    variables["G"]          = UQFF::G;
    variables["c"]          = UQFF::c;
    variables["hbar"]       = UQFF::hbar;
    variables["Lambda"]     = UQFF::Lambda;
    variables["pi"]         = UQFF::pi;
    variables["t_Hubble"]   = UQFF::t_Hubble;
    variables["year_to_s"]  = UQFF::year_to_s;
    variables["H0"]         = UQFF::H0;
    variables["Mpc_to_m"]   = UQFF::Mpc_to_m;
    variables["Omega_m"]    = UQFF::Omega_m;
    variables["Omega_Lambda"]= UQFF::Omega_Lambda;
    variables["B_crit"]     = UQFF::B_crit_magnetar;

    // NGC 346 analogue system parameters (Doc 35)
    variables["M"]          = 1.0e3 * UQFF::M_sun;   // 1000 M☉
    variables["M0"]         = variables["M"];
    variables["r"]          = 2.365e17;               // m  (~7.66 pc)
    variables["z"]          = 0.05;
    variables["t_evolve"]   = 5.0e6 * UQFF::year_to_s;  // 5 Myr evolution timescale
    variables["SFR"]        = 0.05 * UQFF::M_sun;    // 0.05 M☉/yr

    // Outflow parameters
    variables["rho_fluid"]  = 5.0e-21;               // kg/m^3
    variables["v_out"]      = 1.0e5;                 // m/s (100 km/s)
    variables["V"]          = 1.0 / variables["rho_fluid"];

    // EM
    variables["B"]          = 5.0e-5;
    variables["f_sc"]       = 10.0;

    // Quantum
    variables["Delta_x"]    = 1.0e-10;
    variables["Delta_p"]    = UQFF::hbar / variables["Delta_x"];
    variables["integral_psi"] = 1.0;
    variables["delta_rho_over_rho"] = 1.0e-5;

    variables["M_visible"]  = variables["M"];
    variables["M_DM"]       = 0.0;
    variables["M_ext"]      = 0.0;
    variables["r_ext"]      = 1.0;
}

void YoungStarsOutflowsUQFFModule::onVariableUpdated(const std::string& name, double value) {
    if (name == "Delta_x" && value > 0.0)
        variables["Delta_p"] = UQFF::hbar / value;
    else if (name == "M") {
        variables["M0"]       = value;
        variables["M_visible"]= value;
    } else if (name == "rho_fluid" && value > 0.0)
        variables["V"] = 1.0 / value;
}

double YoungStarsOutflowsUQFFModule::computeHtz() const {
    double z = variables.at("z");
    return variables.at("H0") * std::sqrt(
        variables.at("Omega_m") * std::pow(1.0 + z, 3) +
        variables.at("Omega_Lambda")) * 1.0e3 / variables.at("Mpc_to_m");
}

// P_outflow(t) = rho * v_out^2 * (1 + t/t_evolve)
double YoungStarsOutflowsUQFFModule::computePoutflow(double t) const {
    double rho  = variables.at("rho_fluid");
    double vout = variables.at("v_out");
    double tev  = variables.at("t_evolve");
    return rho * vout * vout * (1.0 + t / tev);
}

double YoungStarsOutflowsUQFFModule::computeMsfFactor(double t) const {
    double SFR = variables.at("SFR");
    double M0  = variables.at("M0");
    if (SFR == 0.0 || M0 == 0.0) return 0.0;
    return SFR * (t / variables.at("year_to_s")) / M0;
}

double YoungStarsOutflowsUQFFModule::computeUgSum() const {
    double G   = variables.at("G");
    double M   = variables.at("M");
    double r   = variables.at("r");
    double Ug1 = G * M / (r * r);
    double Ug3 = (variables.at("M_ext") > 0.0) ?
        G * variables.at("M_ext") / (variables.at("r_ext") * variables.at("r_ext")) : 0.0;
    double Ug4 = Ug1 * variables.at("f_sc");
    return Ug1 + Ug3 + Ug4;
}

double YoungStarsOutflowsUQFFModule::computeQuantumPsi() const {
    double unc = std::sqrt(variables.at("Delta_x") * variables.at("Delta_p"));
    return (variables.at("hbar") / unc) *
           variables.at("integral_psi") *
           (2.0 * variables.at("pi") / variables.at("t_Hubble"));
}

double YoungStarsOutflowsUQFFModule::computeFluidTerm(double g_base) const {
    return variables.at("rho_fluid") * variables.at("V") * g_base;
}

double YoungStarsOutflowsUQFFModule::computeDMPert() const {
    double r   = variables.at("r");
    double drr = variables.at("delta_rho_over_rho");
    double pert = drr + 3.0 * variables.at("G") * variables.at("M") / std::pow(r, 3);
    return (variables.at("M_visible") + variables.at("M_DM")) * pert;
}

double YoungStarsOutflowsUQFFModule::computeG(double t) {
    variables["t"] = t;

    double Hz         = computeHtz();
    double expansion  = 1.0 + Hz * t;
    double sc_corr    = 1.0 - variables["B"] / variables["B_crit"];
    double msf        = computeMsfFactor(t);
    double m_factor   = 1.0 + msf;
    double P_out      = computePoutflow(t);

    // F_env: outflow pressure as environmental modifier
    double env_factor = 1.0 + P_out;

    double r       = variables["r"];
    double g_base  = (variables["G"] * variables["M"] * m_factor / (r * r))
                     * expansion * sc_corr * env_factor;

    double ug_sum  = computeUgSum();
    double lam     = variables["Lambda"] * variables["c"] * variables["c"] / 3.0;
    double quantum = computeQuantumPsi();
    double fluid   = computeFluidTerm(g_base);
    double dm_pert = computeDMPert();

    return g_base + ug_sum + lam + quantum + fluid + dm_pert;
}

std::string YoungStarsOutflowsUQFFModule::getEquationText() const {
    return
        "g_YoungStars(r,t) = [G M(t)/r²] * (1+H(t,z)) * (1-B/B_crit) * (1+P_outflow(t))\n"
        "                  + (Ug1 + Ug3' + Ug4) + (Λc²/3)\n"
        "                  + quantum_psi + ρ_fluid V g_base + DM_pert\n"
        "\n"
        "KEY TERM — P_outflow(t) = ρ v_out² * (1 + t/t_evolve)\n"
        "  where t_evolve = 5 Myr, v_out = 100 km/s (NGC346 analogue)\n"
        "  Outflow pressure grows linearly with time, modifying the effective gravity.\n"
        "\n"
        "M(t) = M*(1+SFR*t_yr/M0)  [star formation]\n"
        "H(t,z) = H0 sqrt(Ω_m(1+z)³+Ω_Λ)  (Friedmann)";
}
