/**
 * @file physics_service.h
 * @brief Physics Service interface — C++ facade over all UQFF compute backends
 *
 * PhysicsService provides a single, synchronous call-site for all physics
 * computations in the UQFF pipeline:
 *
 *   source2.cpp (Qt6 GUI) ──► PhysicsService::compute()
 *       │
 *       ├─► MAIN_1_CoAnQi.exe (in-process via ModuleRegistry) — fastest
 *       ├─► physics_backend.cpp (CPU-bound, headless server)
 *       ├─► Python bridge (QCalc / CP2 / CP1 via uqff_ipc.py)
 *       └─► uqff_server.js REST (port 3141, fallback)
 *
 * All concrete backends implement IPhysicsBackend.  PhysicsService
 * selects the backend based on query category and availability.
 *
 * Calibrated constants (GrokThread7b0e namespace in shared_constants.h):
 *   κ  = 0.0005/day,  [SSq] = 0.57,  H_SCm ≈ 0.99,
 *   U_UA = 0.0001,    F_rel = 4.31e33 N  (2024 LEP reanalysis)
 *
 * Author: Daniel T. Murphy
 * Date: March 2026 (v3.1 — GrokThread7b0e961f integration)
 */

#ifndef PHYSICS_SERVICE_H
#define PHYSICS_SERVICE_H

#include "uqff_ipc.h"
#include "python_bridge.h"
#include "shared_constants.h"

#include <string>
#include <memory>
#include <unordered_map>
#include <vector>
#include <functional>

namespace UQFF {

// ─────────────────────────────────────────────────────────────
// PHYSICS QUERY
// ─────────────────────────────────────────────────────────────

/// Caller fills PhysicsQuery and passes it to PhysicsService::compute().
struct PhysicsQuery {
    // ── System parameters ───────────────────────────────────
    std::string system_name;    ///< e.g. "SgrA*", "NGC1365", custom name
    double      mass_kg         = 0.0;
    double      radius_m        = 0.0;
    double      distance_m      = 0.0;
    double      redshift        = 0.0;
    double      t_n             = -2512.0; ///< normalised time

    // ── Calibrated UQFF constants (defaults from GrokThread7b0e) ─
    double      U_UA            = Constants::GrokThread7b0e::F_REL_2024 * 0.0;  // placeholder
    double      SSq             = Constants::SSq;
    double      kappa           = Constants::kappa;
    double      F_rel           = Constants::GrokThread7b0e::F_REL_2024;
    double      H_SCm           = Constants::H_SCm;

    // ── Routing hints ────────────────────────────────────────
    std::string trigger;        ///< Free-text routing hint (see IPC::route_trigger)
    bool        force_python    = false; ///< Skip in-process backend
    bool        force_cp2       = false; ///< Route directly to CP2 (skip QCalc)
    uint32_t    session_id      = 0;

    // ── Extra key-value parameters ───────────────────────────
    std::unordered_map<std::string, double> extra;
};

// ─────────────────────────────────────────────────────────────
// PHYSICS RESULT
// ─────────────────────────────────────────────────────────────

struct PhysicsResult {
    bool        success         = false;
    std::string system_name;
    double      F_U_Bi_i        = 0.0;   ///< Buoyancy integral (N)
    double      compressed_g    = 0.0;   ///< MUGE compressed gravity (m/s²)
    double      resonance_g     = 0.0;   ///< MUGE resonance gravity  (m/s²)
    double      Q_wave          = 0.0;   ///< Q_wave energy density   (J/m³)
    double      rho_vac         = 0.0;   ///< Vacuum density          (kg/m³)
    double      chi_squared     = 0.0;   ///< Goodness of fit

    std::string backend_used;           ///< "in-process", "python:QCalc", "python:CP2", "rest"
    double      elapsed_ms      = 0.0;
    std::string equation_text;          ///< Long-form equation if available
    std::string json_raw;               ///< Full JSON from backend (for inspection)
    std::string error_message;

    explicit operator bool() const { return success; }
};

// ─────────────────────────────────────────────────────────────
// BACKEND INTERFACE
// ─────────────────────────────────────────────────────────────

class IPhysicsBackend {
public:
    virtual ~IPhysicsBackend() = default;

    /// Returns true if this backend can handle the given query.
    virtual bool can_handle(const PhysicsQuery& q) const = 0;

    /// Synchronous compute. Thread-safe implementations preferred.
    virtual PhysicsResult compute(const PhysicsQuery& q) = 0;

    /// Human-readable name for diagnostics.
    virtual const char* name() const = 0;
};

// ─────────────────────────────────────────────────────────────
// CONCRETE BACKEND: Python Bridge (QCalc / CP2 / CP1)
// ─────────────────────────────────────────────────────────────

class PythonBackend final : public IPhysicsBackend {
public:
    explicit PythonBackend(std::string python_exe)
        : python_exe_(std::move(python_exe)) {}

    bool can_handle(const PhysicsQuery& q) const override { return true; } // fallback

    PhysicsResult compute(const PhysicsQuery& q) override;

    const char* name() const override { return "python"; }

private:
    std::string python_exe_;
};

// ─────────────────────────────────────────────────────────────
// PHYSICS SERVICE  (Façade / Router)
// ─────────────────────────────────────────────────────────────

class PhysicsService {
public:
    /// Singleton accessor (source2.cpp calls this once on startup).
    static PhysicsService& instance();

    // ── Backend registration ─────────────────────────────────
    void register_backend(std::shared_ptr<IPhysicsBackend> backend);
    void set_python_exe(const std::string& exe);

    // ── Primary entry-point ──────────────────────────────────
    PhysicsResult compute(const PhysicsQuery& query);

    /// Batch compute over a vector of queries; returns results in same order.
    std::vector<PhysicsResult> compute_batch(const std::vector<PhysicsQuery>& queries);

    // ── Dynamic parameter control ────────────────────────────
    void set_global_parameter(const std::string& key, double value);
    double get_global_parameter(const std::string& key, double default_val = 0.0) const;

    // ── Health / diagnostics ─────────────────────────────────
    bool ping_python(uint32_t timeout_ms = 3000);
    std::string diagnostics() const;

    // ── Calibrated constant accessors ───────────────────────
    static constexpr double KAPPA()  { return Constants::kappa; }
    static constexpr double SSQ()    { return Constants::SSq; }
    static constexpr double H_SCM()  { return Constants::H_SCm; }
    static constexpr double U_UA()   { return Constants::U_UA; }
    static constexpr double F_REL()  { return Constants::GrokThread7b0e::F_REL_2024; }

private:
    PhysicsService() = default;

    std::vector<std::shared_ptr<IPhysicsBackend>> backends_;
    std::unordered_map<std::string, double>       global_params_;
    std::string                                   python_exe_ = "python";
};

// ─────────────────────────────────────────────────────────────
// INLINE IMPLEMENTATION — PhysicsService
// ─────────────────────────────────────────────────────────────

inline PhysicsService& PhysicsService::instance() {
    static PhysicsService s;
    return s;
}

inline void PhysicsService::register_backend(std::shared_ptr<IPhysicsBackend> backend) {
    backends_.push_back(std::move(backend));
}

inline void PhysicsService::set_python_exe(const std::string& exe) {
    python_exe_ = exe;
}

inline void PhysicsService::set_global_parameter(const std::string& key, double value) {
    global_params_[key] = value;
}

inline double PhysicsService::get_global_parameter(const std::string& key, double default_val) const {
    auto it = global_params_.find(key);
    return (it != global_params_.end()) ? it->second : default_val;
}

inline PhysicsResult PhysicsService::compute(const PhysicsQuery& query) {
    for (auto& backend : backends_) {
        if (backend->can_handle(query)) {
            PhysicsResult r = backend->compute(query);
            if (r.success) return r;
            // Fall through to next backend on failure
        }
    }
    PhysicsResult err;
    err.error_message = "No backend could handle query: " + query.system_name;
    return err;
}

inline std::vector<PhysicsResult> PhysicsService::compute_batch(
        const std::vector<PhysicsQuery>& queries)
{
    std::vector<PhysicsResult> results;
    results.reserve(queries.size());
    for (const auto& q : queries)
        results.push_back(compute(q));
    return results;
}

inline bool PhysicsService::ping_python(uint32_t timeout_ms) {
#ifdef _WIN32
    IPC::ComputePayload ping_payload{};
    std::memset(&ping_payload, 0, sizeof(ping_payload));
    auto res = Bridge::call_via_pipe(
        IPC::MessageType::PING, 0,
        &ping_payload, sizeof(ping_payload),
        timeout_ms);
    return res.success;
#else
    (void)timeout_ms;
    return false; // Named pipe not available on non-Windows
#endif
}

inline std::string PhysicsService::diagnostics() const {
    std::string out = "PhysicsService v3.1\n";
    out += "  Backends: " + std::to_string(backends_.size()) + "\n";
    out += "  Python exe: " + python_exe_ + "\n";
    out += "  F_rel:  4.31e33 N (GrokThread7b0e)\n";
    out += "  kappa:  0.0005/day\n";
    out += "  [SSq]:  0.57\n";
    out += "  H_SCm:  0.99\n";
    out += "  U_UA:   0.0001\n";
    return out;
}

} // namespace UQFF

#endif // PHYSICS_SERVICE_H
