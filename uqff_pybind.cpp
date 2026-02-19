/**
 * @file uqff_pybind.cpp
 * @brief PyBind11 bindings for UQFF Star-Magic C++ core
 * 
 * Exposes key PhysicsTerm classes and computation functions to Python.
 * 
 * Build:
 *   pip install pybind11
 *   cmake -DBUILD_PYTHON_BINDINGS=ON ..
 * 
 * Usage (Python):
 *   import uqff_core
 *   
 *   # Create system parameters
 *   sys = uqff_core.SystemParams("Sagittarius A*", M=8.155e36, r=4.4e19)
 *   
 *   # Compute UQFF
 *   result = uqff_core.compute_F_U_Bi_i(sys)
 *   print(f"F_U_Bi_i = {result.F_U_Bi_i} N")
 *   
 *   # Access PhysicsTermRegistry
 *   registry = uqff_core.PhysicsTermRegistry.get_instance()
 *   print(f"Registered terms: {registry.list_terms()}")
 * 
 * Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
 * Framework: UQFF Star-Magic v3.0
 * Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
 */

// Enable M_PI on MSVC
#define _USE_MATH_DEFINES
#include <cmath>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <cmath>

namespace py = pybind11;

// Forward declarations from MAIN_1_CoAnQi.cpp
// These would normally be in a header file

namespace UQFF {

// ═══════════════════════════════════════════════════════════════════════════════
// SYSTEM PARAMETERS (matches SystemParams in MAIN_1_CoAnQi.cpp)
// ═══════════════════════════════════════════════════════════════════════════════

struct SystemParams {
    std::string name;
    double M = 0.0;           // Mass (kg)
    double r = 0.0;           // Radius/distance (m)
    double v = 0.0;           // Velocity (m/s)
    double B0 = 1e-4;         // Magnetic field (T)
    double t = 0.0;           // Time (s)
    double t_n = 0.1;         // Normalized time
    double theta = 0.0;       // Angle (rad)
    double SFR = 0.0;         // Star formation rate (M_sun/yr)
    double z = 0.0;           // Redshift
    double L = 0.0;           // Luminosity (W)
    double rho_vac_UA = 7.09e-36;   // Aether density
    double rho_vac_SCm = 7.09e-37;  // SCm density
    
    SystemParams() = default;
    
    SystemParams(const std::string& name, double M, double r)
        : name(name), M(M), r(r) {}
};


// ═══════════════════════════════════════════════════════════════════════════════
// COMPUTATION RESULTS
// ═══════════════════════════════════════════════════════════════════════════════

struct ComputationResult {
    std::string system_name;
    double F_U_Bi_i = 0.0;
    double g_compressed = 0.0;
    double F_U = 0.0;
    double Ug1 = 0.0;
    double Ug2 = 0.0;
    double Ug3 = 0.0;
    double Ug4 = 0.0;
    double Ug_sum = 0.0;
    std::string status = "success";
    std::string error_message;
};


// ═══════════════════════════════════════════════════════════════════════════════
// PHYSICS CONSTANTS (from shared_constants.h)
// ═══════════════════════════════════════════════════════════════════════════════

namespace Constants {
    constexpr double G = 6.67430e-11;
    constexpr double c = 2.99792458e8;
    constexpr double hbar = 1.054571817e-34;
    constexpr double M_sun = 1.98892e30;
    constexpr double mu_0 = 1.25663706212e-6;
    constexpr double epsilon_0 = 8.8541878128e-12;
    constexpr double rho_vac_UA = 7.09e-36;
    constexpr double rho_vac_SCm = 7.09e-37;
    constexpr double kappa = 0.0005;
    constexpr double SSq = 0.57;
    constexpr double beta_i = 0.603;
    constexpr double k_LENR = 1e-10;
    constexpr double omega_LENR = 1.2e12;
}


// ═══════════════════════════════════════════════════════════════════════════════
// CORE UQFF COMPUTATIONS
// ═══════════════════════════════════════════════════════════════════════════════

/**
 * Compute Ug1: Magnetic dipole gravity component
 */
inline double compute_Ug1(const SystemParams& sys) {
    double M = sys.M;
    double r = sys.r;
    double B = sys.B0;
    
    if (r <= 0 || M <= 0) return 0.0;
    
    // Ug1 = μ₀ * μ² / (4π r⁵) where μ = B * r³
    double mu = B * r * r * r;  // Magnetic moment
    double Ug1 = (Constants::mu_0 * mu * mu) / (4 * M_PI * std::pow(r, 5));
    
    return Ug1;
}

/**
 * Compute Ug2: Charge/reactivity gravity component
 */
inline double compute_Ug2(const SystemParams& sys) {
    double M = sys.M;
    double r = sys.r;
    
    if (r <= 0 || M <= 0) return 0.0;
    
    // Ug2 ~ G M q² / (r² * 4πε₀ r)
    // Simplified for astrophysical scales
    double q_ratio = 1e-10;  // Typical charge-to-mass ratio
    double Ug2 = Constants::G * M * q_ratio / (r * r);
    
    return Ug2;
}

/**
 * Compute Ug3: String rotation gravity component
 */
inline double compute_Ug3(const SystemParams& sys) {
    double M = sys.M;
    double r = sys.r;
    double t = sys.t;
    
    if (r <= 0 || M <= 0) return 0.0;
    
    // Ug3 includes rotational modulation
    double omega = 2 * M_PI / (24 * 3600);  // Reference daily rotation
    double Ug3 = Constants::G * M / (r * r) * std::cos(omega * t) * 1e-6;
    
    return Ug3;
}

/**
 * Compute Ug4: Vacuum concentration gravity component
 */
inline double compute_Ug4(const SystemParams& sys) {
    double M = sys.M;
    double r = sys.r;
    
    if (r <= 0) return 0.0;
    
    // Ug4 = (rho_vac_UA / rho_vac_SCm - 1) * c² / r
    double rho_ratio = Constants::rho_vac_UA / Constants::rho_vac_SCm;
    double Ug4 = (rho_ratio - 1) * Constants::c * Constants::c / r * 1e-20;
    
    return Ug4;
}

/**
 * Compute F_U_Bi_i: Master buoyancy-gravity force (Outside-In)
 * 
 * F_U_Bi_i = Σᵢ (Ug_i × β_i × correction_terms)
 */
inline double compute_F_U_Bi_i(const SystemParams& sys) {
    double Ug1 = compute_Ug1(sys);
    double Ug2 = compute_Ug2(sys);
    double Ug3 = compute_Ug3(sys);
    double Ug4 = compute_Ug4(sys);
    
    double Ug_sum = Ug1 + Ug2 + Ug3 + Ug4;
    
    // Apply buoyancy correction
    double beta_term = Constants::beta_i * Ug_sum;
    
    // LENR coupling
    double lenr_factor = 1.0 + Constants::k_LENR * std::sin(Constants::omega_LENR * sys.t);
    
    // SCm vacuum modulation
    double scm_factor = Constants::SSq * (1.0 - std::exp(-Constants::kappa * sys.t));
    
    // Combine
    double F_U_Bi_i = sys.M * beta_term * lenr_factor * (1.0 + scm_factor);
    
    return F_U_Bi_i;
}

/**
 * Compute compressed gravity g (26-layer triadic)
 */
inline double compute_compressed_g(const SystemParams& sys) {
    double M = sys.M;
    double r = sys.r;
    
    if (r <= 0 || M <= 0) return 0.0;
    
    // Base Newtonian
    double g_N = Constants::G * M / (r * r);
    
    // UQFF corrections
    double Ug1 = compute_Ug1(sys);
    double Ug2 = compute_Ug2(sys);
    double Ug3 = compute_Ug3(sys);
    double Ug4 = compute_Ug4(sys);
    
    // 26-layer sum (simplified - full version in MAIN_1_CoAnQi.cpp)
    double layer_sum = 0.0;
    for (int i = 1; i <= 26; ++i) {
        double Q_i = 1.0;  // Quantum state factor
        double UA_i = 1.0 / (1.0 + i * 0.01);  // Aether correction
        double SCm_i = 0.99;  // SCm modulation
        
        layer_sum += (Ug1 + Ug2 + Ug3 + Ug4) * Q_i * UA_i * SCm_i / 26.0;
    }
    
    return g_N + layer_sum;
}

/**
 * Full UQFF computation returning all components
 */
inline ComputationResult compute_full(const SystemParams& sys) {
    ComputationResult result;
    result.system_name = sys.name;
    
    try {
        result.Ug1 = compute_Ug1(sys);
        result.Ug2 = compute_Ug2(sys);
        result.Ug3 = compute_Ug3(sys);
        result.Ug4 = compute_Ug4(sys);
        result.Ug_sum = result.Ug1 + result.Ug2 + result.Ug3 + result.Ug4;
        result.F_U_Bi_i = compute_F_U_Bi_i(sys);
        result.g_compressed = compute_compressed_g(sys);
        result.F_U = result.F_U_Bi_i;  // Simplified
        result.status = "success";
    } catch (const std::exception& e) {
        result.status = "error";
        result.error_message = e.what();
    }
    
    return result;
}


// ═══════════════════════════════════════════════════════════════════════════════
// PHYSICS TERM BASE CLASS (matches PhysicsTerm in MAIN_1_CoAnQi.cpp)
// ═══════════════════════════════════════════════════════════════════════════════

class PhysicsTerm {
public:
    virtual ~PhysicsTerm() = default;
    
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
    virtual double evaluate(const SystemParams& sys) const = 0;
    virtual bool validate() const { return true; }
};


// ═══════════════════════════════════════════════════════════════════════════════
// PHYSICS TERM REGISTRY (singleton pattern)
// ═══════════════════════════════════════════════════════════════════════════════

class PhysicsTermRegistry {
public:
    static PhysicsTermRegistry& getInstance() {
        static PhysicsTermRegistry instance;
        return instance;
    }
    
    void registerTerm(std::shared_ptr<PhysicsTerm> term) {
        terms_[term->getName()] = term;
    }
    
    std::shared_ptr<PhysicsTerm> getTerm(const std::string& name) const {
        auto it = terms_.find(name);
        if (it != terms_.end()) {
            return it->second;
        }
        return nullptr;
    }
    
    std::vector<std::string> listTerms() const {
        std::vector<std::string> names;
        for (const auto& pair : terms_) {
            names.push_back(pair.first);
        }
        return names;
    }
    
    size_t size() const { return terms_.size(); }
    
private:
    PhysicsTermRegistry() = default;
    std::map<std::string, std::shared_ptr<PhysicsTerm>> terms_;
};


// ═══════════════════════════════════════════════════════════════════════════════
// EXAMPLE PHYSICS TERMS
// ═══════════════════════════════════════════════════════════════════════════════

class DynamicVacuumTerm : public PhysicsTerm {
public:
    DynamicVacuumTerm(double base_energy = 1e-10, double coupling = 1e-15)
        : base_energy_(base_energy), coupling_(coupling) {}
    
    std::string getName() const override { return "DynamicVacuumTerm"; }
    
    std::string getDescription() const override {
        return "Dynamic vacuum energy contribution from quantum fluctuations";
    }
    
    double evaluate(const SystemParams& sys) const override {
        return base_energy_ * std::exp(-coupling_ * sys.r);
    }
    
private:
    double base_energy_;
    double coupling_;
};


class MagneticDipoleTerm : public PhysicsTerm {
public:
    std::string getName() const override { return "MagneticDipoleTerm"; }
    
    std::string getDescription() const override {
        return "Ug1: Magnetic dipole gravity contribution";
    }
    
    double evaluate(const SystemParams& sys) const override {
        return compute_Ug1(sys);
    }
};


} // namespace UQFF


// ═══════════════════════════════════════════════════════════════════════════════
// PYBIND11 MODULE DEFINITION
// ═══════════════════════════════════════════════════════════════════════════════

PYBIND11_MODULE(uqff_core, m) {
    m.doc() = "UQFF Star-Magic C++ Core - Python bindings";
    
    // ─────────────────────────────────────────────────────────────────────────
    // Constants submodule
    // ─────────────────────────────────────────────────────────────────────────
    py::module constants = m.def_submodule("constants", "Physical constants");
    constants.attr("G") = UQFF::Constants::G;
    constants.attr("c") = UQFF::Constants::c;
    constants.attr("hbar") = UQFF::Constants::hbar;
    constants.attr("M_sun") = UQFF::Constants::M_sun;
    constants.attr("mu_0") = UQFF::Constants::mu_0;
    constants.attr("epsilon_0") = UQFF::Constants::epsilon_0;
    constants.attr("rho_vac_UA") = UQFF::Constants::rho_vac_UA;
    constants.attr("rho_vac_SCm") = UQFF::Constants::rho_vac_SCm;
    constants.attr("kappa") = UQFF::Constants::kappa;
    constants.attr("SSq") = UQFF::Constants::SSq;
    constants.attr("beta_i") = UQFF::Constants::beta_i;
    
    // ─────────────────────────────────────────────────────────────────────────
    // SystemParams
    // ─────────────────────────────────────────────────────────────────────────
    py::class_<UQFF::SystemParams>(m, "SystemParams")
        .def(py::init<>())
        .def(py::init<const std::string&, double, double>(),
             py::arg("name"), py::arg("M"), py::arg("r"))
        .def_readwrite("name", &UQFF::SystemParams::name)
        .def_readwrite("M", &UQFF::SystemParams::M)
        .def_readwrite("r", &UQFF::SystemParams::r)
        .def_readwrite("v", &UQFF::SystemParams::v)
        .def_readwrite("B0", &UQFF::SystemParams::B0)
        .def_readwrite("t", &UQFF::SystemParams::t)
        .def_readwrite("t_n", &UQFF::SystemParams::t_n)
        .def_readwrite("theta", &UQFF::SystemParams::theta)
        .def_readwrite("SFR", &UQFF::SystemParams::SFR)
        .def_readwrite("z", &UQFF::SystemParams::z)
        .def_readwrite("L", &UQFF::SystemParams::L)
        .def_readwrite("rho_vac_UA", &UQFF::SystemParams::rho_vac_UA)
        .def_readwrite("rho_vac_SCm", &UQFF::SystemParams::rho_vac_SCm)
        .def("__repr__", [](const UQFF::SystemParams& s) {
            return "<SystemParams '" + s.name + "' M=" + std::to_string(s.M) + " r=" + std::to_string(s.r) + ">";
        });
    
    // ─────────────────────────────────────────────────────────────────────────
    // ComputationResult
    // ─────────────────────────────────────────────────────────────────────────
    py::class_<UQFF::ComputationResult>(m, "ComputationResult")
        .def(py::init<>())
        .def_readwrite("system_name", &UQFF::ComputationResult::system_name)
        .def_readwrite("F_U_Bi_i", &UQFF::ComputationResult::F_U_Bi_i)
        .def_readwrite("g_compressed", &UQFF::ComputationResult::g_compressed)
        .def_readwrite("F_U", &UQFF::ComputationResult::F_U)
        .def_readwrite("Ug1", &UQFF::ComputationResult::Ug1)
        .def_readwrite("Ug2", &UQFF::ComputationResult::Ug2)
        .def_readwrite("Ug3", &UQFF::ComputationResult::Ug3)
        .def_readwrite("Ug4", &UQFF::ComputationResult::Ug4)
        .def_readwrite("Ug_sum", &UQFF::ComputationResult::Ug_sum)
        .def_readwrite("status", &UQFF::ComputationResult::status)
        .def_readwrite("error_message", &UQFF::ComputationResult::error_message)
        .def("__repr__", [](const UQFF::ComputationResult& r) {
            return "<ComputationResult '" + r.system_name + "' F_U_Bi_i=" + std::to_string(r.F_U_Bi_i) + ">";
        });
    
    // ─────────────────────────────────────────────────────────────────────────
    // Computation functions
    // ─────────────────────────────────────────────────────────────────────────
    m.def("compute_Ug1", &UQFF::compute_Ug1, "Compute Ug1 magnetic dipole component",
          py::arg("sys"));
    m.def("compute_Ug2", &UQFF::compute_Ug2, "Compute Ug2 charge/reactivity component",
          py::arg("sys"));
    m.def("compute_Ug3", &UQFF::compute_Ug3, "Compute Ug3 string rotation component",
          py::arg("sys"));
    m.def("compute_Ug4", &UQFF::compute_Ug4, "Compute Ug4 vacuum concentration component",
          py::arg("sys"));
    m.def("compute_F_U_Bi_i", &UQFF::compute_F_U_Bi_i, 
          "Compute F_U_Bi_i master buoyancy-gravity force",
          py::arg("sys"));
    m.def("compute_compressed_g", &UQFF::compute_compressed_g,
          "Compute compressed gravity (26-layer triadic)",
          py::arg("sys"));
    m.def("compute_full", &UQFF::compute_full,
          "Full UQFF computation returning all components",
          py::arg("sys"));
    
    // ─────────────────────────────────────────────────────────────────────────
    // PhysicsTerm base class (trampoline for Python subclassing)
    // ─────────────────────────────────────────────────────────────────────────
    py::class_<UQFF::PhysicsTerm, std::shared_ptr<UQFF::PhysicsTerm>>(m, "PhysicsTerm")
        .def("getName", &UQFF::PhysicsTerm::getName)
        .def("getDescription", &UQFF::PhysicsTerm::getDescription)
        .def("evaluate", &UQFF::PhysicsTerm::evaluate)
        .def("validate", &UQFF::PhysicsTerm::validate);
    
    // ─────────────────────────────────────────────────────────────────────────
    // Concrete PhysicsTerm classes
    // ─────────────────────────────────────────────────────────────────────────
    py::class_<UQFF::DynamicVacuumTerm, UQFF::PhysicsTerm, 
               std::shared_ptr<UQFF::DynamicVacuumTerm>>(m, "DynamicVacuumTerm")
        .def(py::init<double, double>(),
             py::arg("base_energy") = 1e-10,
             py::arg("coupling") = 1e-15);
    
    py::class_<UQFF::MagneticDipoleTerm, UQFF::PhysicsTerm,
               std::shared_ptr<UQFF::MagneticDipoleTerm>>(m, "MagneticDipoleTerm")
        .def(py::init<>());
    
    // ─────────────────────────────────────────────────────────────────────────
    // PhysicsTermRegistry
    // ─────────────────────────────────────────────────────────────────────────
    py::class_<UQFF::PhysicsTermRegistry>(m, "PhysicsTermRegistry")
        .def_static("get_instance", &UQFF::PhysicsTermRegistry::getInstance,
                    py::return_value_policy::reference)
        .def("register_term", &UQFF::PhysicsTermRegistry::registerTerm)
        .def("get_term", &UQFF::PhysicsTermRegistry::getTerm)
        .def("list_terms", &UQFF::PhysicsTermRegistry::listTerms)
        .def("size", &UQFF::PhysicsTermRegistry::size);
    
    // Module version
    m.attr("__version__") = "3.0.0";
}
