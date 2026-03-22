// UFEOrbModule.cpp
// UFE for Red Dwarf Reactor Plasma Orb Experiment (Doc 43).
// t^- = -t_n * exp(pi - t_n) time transformation for plasmoid dynamics.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#include "UFEOrbModule.h"

UFEOrbModule::UFEOrbModule(BatchType batch) : current_batch(batch) {
    variables["G"]      = 6.6743e-11;
    variables["c"]      = 3e8;
    variables["hbar"]   = 1.0546e-34;
    variables["pi"]     = 3.141592653589793;
    variables["gamma"]  = 0.001;         // Decay rate
    variables["fps"]    = 33.3;
    variables["cylinder_r"] = 0.0445;   // m (1.75" radius)
    variables["cylinder_h"] = 0.254;    // m (10")

    // SCm & UA (Red Dwarf specific)
    variables["SCm"]       = 1e15;       // kg/m^3
    variables["SCm_prime"] = 1e15;       // m^{-3}
    variables["UA"]        = 1e-11;      // C

    // Vacuum energies (J/m^3)
    variables["rho_vac_SCm_atomic"] = 1.60e19;
    variables["rho_vac_UA_atomic"]  = 1.60e20;
    variables["E_vac_neb"]          = 7.09e-36;
    variables["E_vac_ISM"]          = 7.09e-37;
    variables["rho_vac_Ug"]         = 5e-89;   // Cosmic
    variables["rho_vac_Um"]         = 1.42e-36;
    variables["rho_vac_Ub"]         = 2.13e-36;
    variables["rho_vac_Ui"]         = 2.84e-36;

    // Ug/Um coefficients
    variables["k1"]       = 1.0;
    variables["beta1"]    = 0.1;
    variables["Omega_g"]  = 1.0;
    variables["M_bh"]     = 1e6 * 1.989e30;  // SMBH kg
    variables["E_react"]  = 1e-20;            // J
    variables["mu1"]      = 1.0;
    variables["phi1"]     = 1.0;
    variables["eta"]      = 1.0;
    variables["lambda1"]  = 0.1;

    // Experiment params
    variables["B_s"]    = 1e-3;    // T
    variables["t_n"]    = 1.0;     // Normalized time
    variables["omega_s"]= 1e3;     // rad/s
    variables["T_s"]    = 300.0;   // K
    variables["RM"]     = 1.0;
    variables["SM"]     = 1.0;
    variables["r"]      = 0.0445;  // m
    variables["t"]      = 9.03;    // s

    variables["plasmoid_count"]    = 40.0;
    variables["energy_per_frame"]  = 0.019;  // J

    setBatch(batch);
}

void UFEOrbModule::setBatch(BatchType batch) {
    current_batch = batch;
    switch (batch) {
        case BatchType::BATCH_31:
            variables["t"]             = 9.03;
            variables["frame_start"]   = 301;
            variables["plasmoid_count"]= 45.0;
            break;
        case BatchType::BATCH_39:
            variables["t"]             = 13.53;
            variables["frame_start"]   = 451;
            variables["plasmoid_count"]= 50.0;
            break;
        case BatchType::EARLY_SEQUENCE:
            variables["t"]             = 0.24;
            variables["plasmoid_count"]= 30.0;
            break;
        case BatchType::MID_SEQUENCE:
            variables["t"]             = 8.73;
            variables["plasmoid_count"]= 40.0;
            break;
        case BatchType::LATE_SEQUENCE:
            variables["t"]             = 13.68;
            variables["plasmoid_count"]= 50.0;
            break;
        default:
            break;
    }
    // Normalized time: t_n = t * fps / total_frames
    variables["t_n"] = variables["t"] * variables["fps"] / 496.0;
}

void UFEOrbModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
    if (name == "SCm") {
        variables["rho_vac_SCm_atomic"] = value * 1e4;
    }
}

void UFEOrbModule::addToVariable(const std::string& name, double delta) {
    if (variables.count(name)) variables[name] += delta;
}

void UFEOrbModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// t^- = -t_n * exp(pi - t_n)  [core plasmoid time transformation]
double UFEOrbModule::computeTminus(double t_n) {
    return -t_n * std::exp(variables["pi"] - t_n);
}

// Σ_i k_i Ug_i (i=1 gravity mode)
double UFEOrbModule::computeUgSum(double t, double r) {
    double t_minus = computeTminus(variables["t_n"]);
    double Ug1 = variables["k1"] *
                 (variables["G"] * variables["M_bh"] / (r * r)) *
                 std::exp(-variables["gamma"] * t_minus) *
                 std::cos(variables["pi"] * variables["t_n"]);
    double beta_t = variables["beta1"] * Ug1 * variables["Omega_g"] *
                    variables["E_react"] / variables["M_bh"];
    return Ug1 - beta_t;
}

// Σ_j μ_j/r_j * (1 - e^{-γt^-} cos(πt_n)) ϕ^j Um_j
double UFEOrbModule::computeUmSum(double t, double r) {
    double t_minus  = computeTminus(variables["t_n"]);
    double exp_cos  = 1.0 - std::exp(-variables["gamma"] * t_minus) *
                             std::cos(variables["pi"] * variables["t_n"]);
    double Um1      = (variables["mu1"] / r) * exp_cos *
                      std::pow(variables["phi1"], 1.0) * variables["rho_vac_Um"];
    return Um1;
}

// Metric + stress-energy term (g_μν + η T_s μν)
double UFEOrbModule::computeMetricTerm() {
    return variables["eta"] * variables["T_s"] * variables["rho_vac_Ug"];
}

// Ub(t^-) buoyancy vacuum term
double UFEOrbModule::computeUbTerm(double t_minus) {
    return variables["rho_vac_Ub"] * std::exp(t_minus);
}

// FU extension: -Σ λ_i Ui E_react
double UFEOrbModule::computeFUExtension() {
    return -variables["lambda1"] * variables["rho_vac_Ui"] * variables["E_react"];
}

// Vacuum energy by type
double UFEOrbModule::computeVacEnergy(const std::string& type) {
    if (type == "SCm") return variables["rho_vac_SCm_atomic"];
    if (type == "UA")  return variables["rho_vac_UA_atomic"];
    return variables["E_vac_neb"];
}

// Plasmoid count estimate: linear 20-50 across 496 frames
double UFEOrbModule::computePlasmoidCount(double timestamp) {
    return 20.0 + 2.0 * (timestamp / 149.88) * 30.0;
}

// Full UP(t)
double UFEOrbModule::computeUP(double t) {
    variables["t"] = t;
    double r        = variables["r"];
    double ug_sum   = computeUgSum(t, r);
    double um_sum   = computeUmSum(t, r);
    double metric   = computeMetricTerm();
    double t_minus  = computeTminus(variables["t_n"]);
    double ub       = computeUbTerm(t_minus);
    double vac_sc   = computeVacEnergy("SCm");
    double vac_ua   = computeVacEnergy("UA");
    double spin_fac = std::cos(variables["omega_s"] * t) * variables["T_s"] * variables["B_s"];
    double sc_fac   = variables["SCm"] * variables["SCm_prime"] * variables["UA"];
    return ug_sum + um_sum + metric + ub + spin_fac * (vac_sc + vac_ua) * sc_fac;
}

double UFEOrbModule::computeFU(double t) {
    return computeUP(t) + computeFUExtension();
}

std::string UFEOrbModule::getEquationText() {
    return
        "UFE Orb Module (Doc 43 — Red Dwarf Reactor):\n"
        "UP(t) = Σ_i [k_i Ug_i(r,t^-,ω_s,T_s,B_s,SCm,UA)]\n"
        "      + Σ_j [μ_j/r_j (1 - e^{-γ t^-} cos(π t_n)) ϕ^j Um_j]\n"
        "      + (g_μν + η T_s μν) + Ub(t^-) + spin_fac*(vac_SCm+vac_UA)*SC_fac\n\n"
        "FU = UP(t) - Σ λ_i Ui E_react\n\n"
        "t^- = -t_n * exp(π - t_n)   [plasmoid time transformation]\n"
        "t_n = t * fps / 496         [normalized frame time]\n\n"
        "Vacuum: ρ_vac[SCm]=1.60e19 J/m³(atomic), ρ_vac[UA]=1.60e20 J/m³\n"
        "Cylinder: r=0.0445m, h=0.254m; SCm=1e15 kg/m³; UA=1e-11 C; fps=33.3";
}

std::string UFEOrbModule::getSolutions(double t) {
    double r        = variables["r"];
    double t_n      = variables["t_n"];
    double t_minus  = computeTminus(t_n);
    double ug       = computeUgSum(t, r);
    double um       = computeUmSum(t, r);
    double metric   = computeMetricTerm();
    double ub       = computeUbTerm(t_minus);
    double fu_ext   = computeFUExtension();
    double up_total = computeUP(t);
    double fu_total = computeFU(t);
    double plasm    = computePlasmoidCount(t);

    std::ostringstream ss;
    ss << std::scientific << std::setprecision(4);
    ss << "UFEOrbModule Solutions t=" << t << " s (Batch=" << static_cast<int>(current_batch) << "):\n";
    ss << "  t_n    = " << t_n    << "\n";
    ss << "  t^-    = " << t_minus << "\n";
    ss << "  Ug_sum = " << ug      << " J/m^3\n";
    ss << "  Um_sum = " << um      << " J/m^3\n";
    ss << "  Metric = " << metric  << " J/m^3\n";
    ss << "  Ub     = " << ub      << " J/m^3\n";
    ss << "  UP(t)  = " << up_total<< " J/m^3\n";
    ss << "  FU(t)  = " << fu_total<< " J/m^3\n";
    ss << "  Plasmoids ~ " << plasm << "/frame\n";
    ss << "  E/frame  = " << variables["energy_per_frame"] << " J\n";
    return ss.str();
}

void UFEOrbModule::printVariables() {
    std::cout << "UFEOrbModule variables (Batch=" << static_cast<int>(current_batch) << "):\n";
    for (const auto& p : variables) {
        std::cout << "  " << p.first << " = " << p.second << "\n";
    }
}
