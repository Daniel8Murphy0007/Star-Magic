// V838MonUQFFModule.cpp
// Implementation of MUGE for V838 Monocerotis Light Echo Evolution.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "V838MonUQFFModule.h"
#include <complex>

V838MonUQFFModule::V838MonUQFFModule() {
    variables["c"] = 3e8;
    variables["G"] = 6.6743e-11;
    variables["hbar"] = 1.0546e-34;
    variables["pi"] = 3.141592653589793;
    double M_sun_val = 1.989e30;
    double L_sun_val = 3.826e26;

    // V838 Mon parameters
    variables["M_s"] = 8 * M_sun_val;
    variables["L_outburst"] = 600000 * L_sun_val;  // ~2.3e38 W
    variables["rho_0"] = 1e-22;                    // kg/m^3 dust
    variables["sigma_scatter"] = 1e-12;            // m^2 dust grain
    variables["k1"] = 1.0;
    variables["mu_s"] = 1.0;
    variables["alpha"] = 0.0005;
    variables["beta"] = 1.0;
    variables["t_n"] = 0.0;
    variables["delta_def"] = 0.01 * std::sin(0.001 * 1e7);
    variables["f_TRZ"] = 0.1;
    variables["rho_vac_UA"] = 7.09e-36;
    variables["rho_vac_SCm"] = 7.09e-37;
    variables["t"] = 3 * 3.156e7;                  // 3 years s

    variables["scale_macro"] = 1e-12;
}

void V838MonUQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding." << std::endl;
        variables[name] = value;
    }
    if (name == "t") {
        variables["delta_def"] = 0.01 * std::sin(0.001 * value);
    }
}

void V838MonUQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) variables[name] += delta;
    else variables[name] = delta;
}
void V838MonUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

double V838MonUQFFModule::computeUg1(double t, double r) {
    double grad_term = variables["M_s"] / (r * r * r);
    double exp_decay = std::exp(-variables["alpha"] * t);
    double cos_phase = std::cos(variables["pi"] * variables["t_n"]);
    double delta = variables["delta_def"];
    return variables["k1"] * variables["mu_s"] * grad_term * exp_decay * cos_phase * (1 + delta);
}

double V838MonUQFFModule::computeRhodust(double r, double t) {
    double ug1 = computeUg1(t, r);
    return variables["rho_0"] * std::exp(-variables["beta"] * ug1);
}

double V838MonUQFFModule::computeIechoBase(double r) {
    return variables["L_outburst"] / (4 * variables["pi"] * r * r);
}

double V838MonUQFFModule::computeTRZCorrection() {
    return 1.0 + variables["f_TRZ"];
}

double V838MonUQFFModule::computeUAscCorrection() {
    return 1.0 + (variables["rho_vac_UA"] / variables["rho_vac_SCm"]);
}

double V838MonUQFFModule::computeIecho(double t, double r) {
    variables["t"] = t;
    double rho_dust = computeRhodust(r, t);
    double i_base = computeIechoBase(r);
    double trz = computeTRZCorrection();
    double ua_sc = computeUAscCorrection();
    return i_base * variables["sigma_scatter"] * rho_dust * trz * ua_sc;
}

std::string V838MonUQFFModule::getEquationText() {
    return "I_echo(r, t) = [L_outburst / (4 pi (c t)^2)] * sigma_scatter * rho_0 * "
           "exp(-beta [k1 mu_s(t, rho_vac,SCm) grad(M_s / (c t)) e^{-alpha t} cos(pi t_n) (1 + delta_def)]) "
           "* (1 + f_TRZ) * (1 + rho_vac,UA / rho_vac,SCm)\n"
           "Where: r_echo(t) = c t; delta_def = 0.01 sin(0.001 t); grad(M_s / r) ~ M_s / r^3;\n"
           "L_outburst ~2.3e38 W; rho_0=1e-22 kg/m^3; f_TRZ=0.1.\n"
           "Insights: Attractive Ug1 modulates dust density; repulsive [UA] corrects propagation.\n"
           "Adaptations: Hubble ACS 2004. I_echo ~1e-20 W/m^2 at t=3 yr, r=9e15 m.";
}

void V838MonUQFFModule::printVariables() {
    std::cout << "V838 Mon Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}
