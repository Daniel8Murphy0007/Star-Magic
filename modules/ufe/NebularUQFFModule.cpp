// NebularUQFFModule.cpp
// UQFF for Nebular Cloud Analysis (Drawing 32), LENR, Higgs, NGC 346.
// Equations from grok_share_e70525fa.txt Doc 43.b.
// Copyright - Daniel T. Murphy.

#include "NebularUQFFModule.h"

NebularUQFFModule::NebularUQFFModule(NebularSystemType sys) : current_system(sys) {
    // ---- Universal UQFF vacuum / quantum constants ----
    variables["rho_vac_SCm"]   = 2.39e-22;   // J/m^3 level-13 vacuum energy [SCm]
    variables["rho_vac_UA"]    = 2.39e-21;   // J/m^3 level-13 + 1 decade [UA]
    variables["SSq"]           = 0.57;        // UQFF calibrated [SSq]
    variables["n26"]           = 26.0;        // 26-D UQFF layers
    variables["k_eta"]         = 1e-113;      // neutron-rate coupling k_η
    variables["k_trans"]       = 1.0;         // transmutation prefactor
    variables["k_Higgs"]       = 1.0;         // Higgs coupling
    variables["k3"]            = 1.0;         // Ug3 prefactor
    variables["kappa_V"]       = 1.0;         // volume coupling
    variables["kappa_F"]       = 1.0;         // frequency coupling
    variables["pi"]            = M_PI;
    variables["c"]             = 2.998e8;     // m/s
    variables["hbar"]          = 1.055e-34;   // J·s
    variables["m_e"]           = 9.109e-31;   // kg
    variables["e_charge"]      = 1.602e-19;   // C
    variables["mu"]            = 1.0;         // reduced mass / coupling (dimensionless)
    variables["omega_c"]       = 1.0;         // DNA oscillation frequency
    variables["omega_i"]       = 1.0;         // inertial frequency
    variables["M_stars"]       = 1.0e30;      // stellar mass [kg]

    setSystem(sys);
}

void NebularUQFFModule::setSystem(NebularSystemType sys) {
    current_system = sys;
    switch (sys) {
    case NebularSystemType::NEBULA_CLOUD:
        variables["n_e"]      = 1e10;         // electron/m^3
        variables["sigma"]    = 1e-28;        // cross-section m^2
        variables["v"]        = 1e5;          // relative velocity m/s
        variables["Omega"]    = 1e9;          // frequency Hz
        variables["V_vol"]    = 1e6;          // volume m^3
        variables["non_local_t"] = 0.0;
        break;
    case NebularSystemType::NGC346:
        variables["n_e"]      = 5e10;
        variables["sigma"]    = 2e-28;
        variables["v"]        = 2e5;
        variables["Omega"]    = 5e9;
        variables["SFR"]      = 0.06;         // solar masses / yr
        variables["M_stars"]  = 5e4 * 1.989e30;
        variables["non_local_t"] = 0.0;
        variables["delta_lambda_over_lambda"] = -3.33e-5; // blueshift
        break;
    case NebularSystemType::LENR_CELL:
        variables["n_e"]      = 1e14;
        variables["sigma"]    = 1e-24;
        variables["v"]        = 1e3;
        variables["Omega"]    = 1e12;
        variables["V_vol"]    = 1e-6;
        variables["non_local_t"] = 0.0;
        break;
    case NebularSystemType::HIGGS_PHYSICS:
        variables["n_e"]      = 1e15;
        variables["sigma"]    = 1e-24;
        variables["v"]        = 3e7;
        variables["Omega"]    = 1.25e34;
        variables["kappa_F"]  = 1.0;
        variables["mu"]       = 1.0;
        variables["non_local_t"] = 0.0;
        break;
    case NebularSystemType::GENERIC:
    default:
        variables["n_e"]      = 1e10;
        variables["sigma"]    = 1e-28;
        variables["v"]        = 1e5;
        variables["Omega"]    = 1e9;
        variables["V_vol"]    = 1e6;
        variables["non_local_t"] = 0.0;
        break;
    }
}

void NebularUQFFModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
}
void NebularUQFFModule::addToVariable(const std::string& name, double delta) {
    variables[name] += delta;
}
void NebularUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    variables[name] -= delta;
}

// Non-local correction: [SSq]^{n26} * exp(-(pi + t))
double NebularUQFFModule::computeNonLocalTerm(double t, int n26) {
    double SSq = variables["SSq"];
    double pi  = variables["pi"];
    return std::pow(SSq, n26) * std::exp(-(pi + t));
}

// Ug3 star formation temperature (internal)
// Ug3 = k3 * M_stars * 3.38e20 / r^3 * cos(theta) * 1e46 * (1 + non_local)^n
double NebularUQFFModule::computeUg3internal(double t, double r, double theta, int n) {
    auto nl = computeNonLocalTerm(t, (int)variables["n26"]);
    double k3 = variables["k3"];
    double M  = variables["M_stars"];
    if (r <= 0.0) r = 1.0;
    return k3 * M * 3.38e20 / (r*r*r) * std::cos(theta) * 1.0e46 * std::pow(1.0 + nl, n);
}

// --- Public methods ---

// Electric field: k_eta * e * Omega / m_e * sqrt(n_e * sigma * v) * kappa_V
double NebularUQFFModule::computeElectricField() {
    double ke    = variables["k_eta"];
    double e     = variables["e_charge"];
    double Omega = variables["Omega"];
    double me    = variables["m_e"];
    double ne    = variables["n_e"];
    double sig   = variables["sigma"];
    double v     = variables["v"];
    double kappV = variables["kappa_V"];
    return ke * e * Omega / me * std::sqrt(ne * sig * v) * kappV;
}

// Neutron rate: k_eta * n_e * sigma * v
double NebularUQFFModule::computeNeutronRate() {
    return variables["k_eta"] * variables["n_e"] * variables["sigma"] * variables["v"];
}

// Transmutation energy: k_trans * rho_vac_Ug4 * [SSq]^26 * exp(-(pi+t))
// At t=0: exp(-(pi)) ~ 0.0432
double NebularUQFFModule::computeTransmutationEnergy() {
    double kt   = variables["k_trans"];
    double rho  = variables["rho_vac_SCm"];
    double SSq  = variables["SSq"];
    double pi   = variables["pi"];
    return kt * rho * std::pow(SSq, 26.0) * std::exp(-pi);
}

// Higgs pole mass coupling: k_Higgs * 125 * mu * kappa_F
double NebularUQFFModule::computeHiggsMass() {
    return variables["k_Higgs"] * 125.0 * variables["mu"] * variables["kappa_F"];
}

// Star formation / Ug3 temperature
double NebularUQFFModule::computeStarFormationTemp(double t, double r) {
    return computeUg3internal(t, r, 0.0, 1);
}

// Blueshift radial velocity: v_r = c * (delta_lambda / lambda)
double NebularUQFFModule::computeRadialVelocity(double dlambda_over_lambda) {
    return variables["c"] * dlambda_over_lambda;
}

// Neutrino proto energy: rho_vac_UA_prime_SCm * exp(-non_local) * Um / rho_vac_UA
double NebularUQFFModule::computeNeutrinoEnergy(double t) {
    double nl  = computeNonLocalTerm(t, (int)variables["n26"]);
    double rho_SCm = variables["rho_vac_SCm"];
    double rho_UA  = variables["rho_vac_UA"];
    double Um  = rho_SCm * variables["sigma"] * variables["v"];
    return rho_SCm * std::exp(-nl) * Um / rho_UA;
}

// Universal decay: (rho_SCm / rho_UA) * exp(-non_local) * 0.963
double NebularUQFFModule::computeUniversalDecay(double t) {
    double nl  = computeNonLocalTerm(t, (int)variables["n26"]);
    return (variables["rho_vac_SCm"] / variables["rho_vac_UA"]) * std::exp(-nl) * 0.963;
}

// DNA energy: Um * cos(omega_c * t)
double NebularUQFFModule::computeDNAEnergy(double t) {
    double Um  = variables["rho_vac_SCm"] * variables["sigma"] * variables["v"];
    return Um * std::cos(variables["omega_c"] * t);
}

// Buoyancy ratio: (rho_UA / rho_SCm) * V_little / V_big
double NebularUQFFModule::computeBuoyancyRatio(double V_little, double V_big) {
    return (variables["rho_vac_UA"] / variables["rho_vac_SCm"]) * V_little / V_big;
}

// Geometric angle: average atan2 over star pairs
double NebularUQFFModule::computeGeometricAngle(const std::vector<std::pair<double,double>>& positions) {
    if (positions.size() < 2) return 0.0;
    double sum = 0.0;
    int count  = 0;
    for (size_t i = 0; i + 1 < positions.size(); ++i) {
        double dy = positions[i+1].second - positions[i].second;
        double dx = positions[i+1].first  - positions[i].first;
        sum += std::atan2(dy, dx);
        ++count;
    }
    return count > 0 ? sum / count : 0.0;
}

// Scenario accuracy table (demonstration)
double NebularUQFFModule::computeAccuracy(const std::string& scenario) {
    if (scenario == "LENR_neutron_rate")  return 0.95;
    if (scenario == "blueshift_NGC346")   return 0.97;
    if (scenario == "Higgs_125GeV")       return 0.99;
    if (scenario == "DNA_oscillation")    return 0.90;
    if (scenario == "buoyancy")           return 0.93;
    return 0.85; // generic
}

// Weighted UQFF sum: combines key terms at time t
double NebularUQFFModule::computeUQFF(double t) {
    double nl  = computeNonLocalTerm(t, (int)variables["n26"]);
    double E   = computeElectricField();
    double nu  = computeNeutrinoEnergy(t);
    double dec = computeUniversalDecay(t);
    double dna = computeDNAEnergy(t);
    double Ug3 = computeUg3internal(t, 1.0e18, 0.0, 1);
    // Weighted superposition
    return (E + nu + dec + dna + Ug3) * (1.0 + nl);
}

std::string NebularUQFFModule::getEquationText() {
    return
        "NebularUQFFModule — Drawing 32 Equations:\n"
        "  E_field  = k_eta * e * Omega / m_e * sqrt(n_e * sigma * v) * kappa_V\n"
        "  n_rate   = k_eta * n_e * sigma * v\n"
        "  E_trans  = k_trans * rho_vac_SCm * [SSq]^26 * exp(-(pi+t))\n"
        "  M_Higgs  = k_Higgs * 125 * mu * kappa_F\n"
        "  Ug3      = k3 * M_stars * 3.38e20 / r^3 * cos(theta) * 1e46 * (1+NL)^n\n"
        "  v_r      = c * (delta_lambda / lambda)\n"
        "  NuProto  = rho_SCm * exp(-NL) * Um / rho_UA\n"
        "  Decay    = (rho_SCm/rho_UA) * exp(-NL) * 0.963\n"
        "  E_DNA    = Um * cos(omega_c * t)\n"
        "  F_buoy   = (rho_UA/rho_SCm) * V_little/V_big\n"
        "  NL       = [SSq]^n26 * exp(-(pi+t))\n";
}

std::string NebularUQFFModule::getSolutions(double t) {
    double nl  = computeNonLocalTerm(t, (int)variables["n26"]);
    double E   = computeElectricField();
    double nr  = computeNeutronRate();
    double Et  = computeTransmutationEnergy();
    double Mh  = computeHiggsMass();
    double vr  = computeRadialVelocity(variables.count("delta_lambda_over_lambda") ?
                     variables["delta_lambda_over_lambda"] : -3.33e-5);
    double nu  = computeNeutrinoEnergy(t);
    double dec = computeUniversalDecay(t);
    double dna = computeDNAEnergy(t);
    double fb  = computeBuoyancyRatio(1.0, 33.0);
    double uqff = computeUQFF(t);

    std::ostringstream oss;
    oss << std::scientific << std::setprecision(4);
    oss << "=== NebularUQFFModule Solutions (t=" << t << " s) ===\n";
    oss << "Non-local term NL : " << nl    << "\n";
    oss << "Electric field    : " << E     << " V/m\n";
    oss << "Neutron rate      : " << nr    << " s^-1\n";
    oss << "Transmutation E   : " << Et    << " J\n";
    oss << "Higgs mass (arb)  : " << Mh    << " GeV\n";
    oss << "Radial velocity   : " << vr    << " m/s\n";
    oss << "Neutrino energy   : " << nu    << " J\n";
    oss << "Universal decay   : " << dec   << " (ratio)\n";
    oss << "DNA energy        : " << dna   << " J\n";
    oss << "Buoyancy ratio    : " << fb    << "\n";
    oss << "UQFF combined     : " << uqff  << "\n";
    return oss.str();
}

void NebularUQFFModule::printVariables() {
    std::cout << "=== NebularUQFFModule Variables ===\n";
    for (auto& [k,v] : variables)
        std::cout << "  " << k << " = " << v << "\n";
}
