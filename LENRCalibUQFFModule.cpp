// LENRCalibUQFFModule.cpp
// K_n neutron calibration with non-local exponential and quantum state coupling.
// NOVEL: eta(t,n) = k_eta * exp(-[SS_q]^n * 2^6 * exp(-pi - t/yr)) * Um / rho_vac_UA
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#include "LENRCalibUQFFModule.h"

LENRCalibUQFFModule::LENRCalibUQFFModule() : m_scenario("hydride") {
    variables["pi"] = 3.141592653589793;
    variables["q"] = 1.602e-19;
    variables["c"] = 3e8;
    variables["year_to_s"] = 3.156e7;
    variables["hbar"] = 1.0546e-34;

    // [SS_q] calibrated quantum state constant
    variables["SS_q"] = 0.57;

    // Vacuum densities
    variables["rho_vac_UA"] = 7.09e-36;
    variables["rho_vac_SCm"] = 7.09e-37;

    // LENR k_eta scenarios
    variables["k_eta"] = 1e13;      // hydride default
    variables["B"] = 0.1;
    variables["q_charge"] = variables["q"];
    variables["v_drift"] = 1e6;
    variables["t"] = 0.0;
}

void LENRCalibUQFFModule::setScenario(const std::string& scenario) {
    m_scenario = scenario;
    if (scenario == "hydride") {
        variables["k_eta"] = 1e13;
        variables["B"] = 0.1;
    } else if (scenario == "wires") {
        variables["k_eta"] = 1e8;
        variables["B"] = 1.0;
    } else if (scenario == "corona") {
        variables["k_eta"] = 7e-3;
        variables["B"] = 1e-4;
    }
}

void LENRCalibUQFFModule::updateVariable(const std::string& name, double value) { variables[name] = value; }
void LENRCalibUQFFModule::addToVariable(const std::string& name, double delta) { variables[name] += delta; }
void LENRCalibUQFFModule::subtractFromVariable(const std::string& name, double delta) { variables[name] -= delta; }

// NOVEL: non-local exponential = exp(-[SS_q]^n * 2^6 * exp(-pi - t/yr))
double LENRCalibUQFFModule::computeNonLocalExp(int n, double t) {
    double t_yr = t / variables["year_to_s"];
    double SS_q_n = std::pow(variables["SS_q"], n);
    double inner = std::exp(-variables["pi"] - t_yr);
    return std::exp(-SS_q_n * 64.0 * inner);   // 2^6 = 64
}

// delta_n = (2*pi)^(n/6)
double LENRCalibUQFFModule::computeDeltaN(int n) {
    return std::pow(2.0 * variables["pi"], n / 6.0);
}

// rho_vac_state(n,t) = 1e-23 * (0.1)^n * nonLocalExp(n,t)
double LENRCalibUQFFModule::computeRhoVacState(int n, double t) {
    return 1e-23 * std::pow(0.1, n) * computeNonLocalExp(n, t);
}

double LENRCalibUQFFModule::computeUm() {
    return variables["q_charge"] * variables["v_drift"] * variables["B"];
}

// eta(t,n) = k_eta * nonLocalExp(n,t) * Um / rho_vac_UA
double LENRCalibUQFFModule::computeEta(double t, int n) {
    variables["t"] = t;
    double non_local = computeNonLocalExp(n, t);
    double Um = computeUm();
    return variables["k_eta"] * non_local * std::abs(Um) / variables["rho_vac_UA"];
}

// rho_vac_UA':SCm(n,t) — quantum state-coupled vacuum density
double LENRCalibUQFFModule::computeRhoVacUASCm(int n, double t) {
    return computeRhoVacState(n, t);
}

std::string LENRCalibUQFFModule::getEquationText() {
    return "eta(t,n) = k_eta * exp(-[SS_q]^n * 2^6 * exp(-pi - t/yr)) * Um / rho_vac_UA\n"
           "NOVEL non-local exponential: exp(-[SS_q]^n * 2^6 * exp(-pi - t/yr))\n"
           "delta_n = (2*pi)^(n/6); [SS_q]=0.57; k_eta: hydride=1e13, wires=1e8, corona=7e-3\n"
           "rho_vac_state(n,t) = 1e-23*(0.1)^n * nonLocalExp(n,t); 3 LENR scenarios";
}
void LENRCalibUQFFModule::printVariables() {
    std::cout << "LENR Calibration Variables [scenario: " << m_scenario << "]:\n";
    for (const auto& pair : variables) std::cout << pair.first << " = " << std::scientific << pair.second << "\n";
}
