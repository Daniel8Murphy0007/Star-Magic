/**
 * @file physics_service.h
 * @brief Physics Backend Service for UQFF Star-Magic
 * 
 * Provides headless service mode for source2(HEAD PROGRAM).cpp, enabling:
 * - IPC server accepting physics calculation requests
 * - VR runtime field data streaming via shared memory
 * - gRPC endpoints for structured commands (Phase 3 complete)
 * - Integration with MAIN_1_CoAnQi physics engine
 * - Self-Expand: Dynamic term registration from VR/external (v3.1)
 * - Self-Update: Runtime parameter tuning for κ, [SSq], etc. (v3.1)
 * - Self-Simulate: Time evolution with VR streaming (v3.1)
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.1
 * Phase: 3.1 - Self-Expanding Physics Backend
 */

#ifndef PHYSICS_SERVICE_H
#define PHYSICS_SERVICE_H

#include "uqff_ipc.h"
#include <string>
#include <vector>
#include <memory>
#include <thread>
#include <atomic>
#include <functional>
#include <unordered_map>
#include <map>
#include <shared_mutex>

// Forward declarations for Qt types (optional GUI support)
class QCoreApplication;

namespace UQFF {

// Forward declarations
struct SystemParams;
struct FieldResult;

/**
 * @struct ServiceConfig
 * @brief Configuration for physics service startup
 */
struct ServiceConfig {
    // IPC Settings
    std::string ipc_channel_name = "uqff_physics";
    size_t shm_buffer_size = 16 * 1024 * 1024;  // 16 MB default
    
    // gRPC Settings (Phase 3 - Full Implementation)
    std::string grpc_address = "localhost:50051";
    bool enable_grpc = true;  // Enabled in Phase 3
    
    // Threading
    int worker_threads = 4;
    bool use_thread_pool = true;
    
    // Logging
    bool verbose = false;
    std::string log_file = "physics_service.log";
    
    // Integration
    bool enable_main1_integration = true;  // Use MAIN_1_CoAnQi calculators
    bool enable_python_integration = true; // Use Python calculators
    
    // Parse from command line
    static ServiceConfig from_args(int argc, char* argv[]);
};

/**
 * @struct FieldRequest
 * @brief Request for field calculation
 */
struct FieldRequest {
    std::string system_name;
    double r;       // Radial distance [m]
    double t;       // Time [s]
    double theta;   // Angle [rad]
    uint32_t flags; // Calculation options
    
    // Optional: system parameters override
    double mass = 0;
    double radius = 0;
    double magnetic_field = 0;
    double spin_period = 0;
};

/**
 * @struct FieldResponse
 * @brief Response containing calculated field values
 */
struct FieldResponse {
    bool success = false;
    std::string error_message;
    
    // UQFF Field Components
    double F_U = 0;         // Total unified field
    double Ug1 = 0;         // Magnetic dipole field
    double Ug2 = 0;         // Charge-reactivity field
    double Ug3 = 0;         // String rotation field
    double Ug4 = 0;         // Vacuum concentration field
    double Um = 0;          // Magnetism field
    double Ubi = 0;         // Buoyancy opposition
    
    // MUGE Compressed gravity
    double g_compressed = 0;
    
    // Validation metrics
    double residual = 0;    // |calculated - observed|
    double confidence = 0;  // 0-1 confidence score
    
    // Timing
    double compute_time_ms = 0;
};

// ============================================================================
// SELF-EXPAND: Dynamic Physics Terms (v3.1)
// ============================================================================

/**
 * @enum DynamicTermType
 * @brief Types of dynamically registered physics terms
 */
enum class DynamicTermType : uint32_t {
    GRAVITY_MODIFIER    = 0x0001,  // Modifies Ug1-Ug4
    VACUUM_ENERGY       = 0x0002,  // Vacuum fluctuation term
    QUANTUM_COUPLING    = 0x0003,  // Quantum field coupling
    DARK_MATTER         = 0x0004,  // Dark matter halo contribution
    TORSION_FIELD       = 0x0005,  // Spacetime torsion
    CUSTOM              = 0xFFFF   // User-defined
};

/**
 * @struct DynamicTerm
 * @brief Runtime-registered physics term from VR or external source
 */
struct DynamicTerm {
    std::string name;
    DynamicTermType type = DynamicTermType::CUSTOM;
    std::string description;
    double coefficient = 1.0;           // Multiplier
    double base_value = 0.0;            // Base contribution
    double r_dependence = -2.0;         // Power law: term ∝ r^exponent
    double t_dependence = 0.0;          // Time dependence exponent
    bool enabled = true;
    uint64_t registered_at = 0;         // Timestamp
    std::string source;                 // "VR", "external", "script"
    
    // Calculate term contribution at (r, t)
    double evaluate(double r, double t) const {
        if (!enabled) return 0.0;
        double r_factor = (r_dependence != 0.0) ? std::pow(r, r_dependence) : 1.0;
        double t_factor = (t_dependence != 0.0) ? std::pow(t, t_dependence) : 1.0;
        return coefficient * base_value * r_factor * t_factor;
    }
};

// ============================================================================
// SELF-UPDATE: Calibrated Parameters (v3.1)
// ============================================================================

/**
 * @struct CalibratedParameters
 * @brief UQFF calibrated constants (runtime-tunable)
 * 
 * Default values from UQFF validation (κ=0.0005/day, [SSq]=0.57, etc.)
 */
struct CalibratedParameters {
    // Core UQFF parameters
    double kappa = 0.0005;              // κ - decay rate [1/day]
    double SSq = 0.57;                  // [SSq] - string squeeze factor
    double H_SCm = 0.99;                // H_SCm - superconducting efficiency
    double U_UA = 0.0001;               // U_UA - vacuum coupling
    double k_eta = 1e-113;              // k_η - Planck area ratio
    double beta_i = 0.603;              // β_i - buoyancy coefficient
    
    // MUGE parameters
    double alpha_DPM = 1.0;             // Dark matter amplitude
    double Omega_Lambda = 0.69;         // Cosmological constant fraction
    double H_0 = 67.4;                  // Hubble constant [km/s/Mpc]
    
    // 26-layer compression
    int active_layers = 26;             // Number of active compression layers
    double layer_coupling[26] = {1.0};  // Per-layer coupling factors
    
    // Learning / optimization
    double learning_rate = 0.001;       // For self-optimization
    bool auto_optimize = false;
    
    // Metadata
    uint64_t last_updated = 0;
    std::string update_source;
};

// ============================================================================
// SELF-SIMULATE: Time Evolution (v3.1)
// ============================================================================

/**
 * @struct SimulationConfig
 * @brief Configuration for time evolution simulation
 */
struct SimulationConfig {
    std::string system_name;
    double r_start = 1e6;               // Starting radius [m]
    double r_end = 1e12;                // Ending radius [m]
    double t_start = 0.0;               // Start time [s]
    double t_end = 3.156e7;             // End time [s] (1 year default)
    double dt = 3600.0;                 // Time step [s] (1 hour default)
    double dr = 0.0;                    // Radial increment (0 = fixed r)
    int frames = 1000;                  // Number of output frames
    bool stream_to_vr = true;           // Stream results to VR
    bool save_to_file = false;          // Save to CSV/JSON
    std::string output_path;
    
    // Derived
    int total_steps() const { 
        return static_cast<int>((t_end - t_start) / dt); 
    }
    int steps_per_frame() const { 
        return std::max(1, total_steps() / frames); 
    }
};

/**
 * @struct SimulationFrame
 * @brief Single frame of simulation output (streamed to VR)
 */
struct SimulationFrame {
    uint64_t frame_number = 0;
    double t = 0.0;                     // Current time [s]
    double r = 0.0;                     // Current radius [m]
    double theta = 0.0;                 // Current angle [rad]
    
    // Field values at this point
    double F_U = 0.0;
    double Ug1 = 0.0, Ug2 = 0.0, Ug3 = 0.0, Ug4 = 0.0;
    double Um = 0.0, Ubi = 0.0;
    double g_compressed = 0.0;
    
    // Derived quantities
    double dF_dt = 0.0;                 // Field time derivative
    double orbital_velocity = 0.0;      // v = sqrt(F_U * r)
    double escape_velocity = 0.0;       // v_esc = sqrt(2 * F_U * r)
    
    // Progress
    double progress = 0.0;              // 0-1 completion
    bool is_final = false;
};

/**
 * @class PhysicsService
 * @brief Headless physics computation service
 * 
 * Runs as a backend service, accepting IPC requests from:
 * - vr_runtime (VR frame updates, gesture input)
 * - External tools (testing, validation)
 * - Distributed computing nodes
 */
class PhysicsService {
public:
    /**
     * @brief Construct a new Physics Service
     * @param config Service configuration
     */
    explicit PhysicsService(const ServiceConfig& config);
    ~PhysicsService();
    
    // Non-copyable
    PhysicsService(const PhysicsService&) = delete;
    PhysicsService& operator=(const PhysicsService&) = delete;
    
    /**
     * @brief Start the physics service
     * @return true if started successfully
     */
    bool start();
    
    /**
     * @brief Stop the physics service
     */
    void stop();
    
    /**
     * @brief Check if service is running
     * @return true if running
     */
    bool is_running() const { return running_.load(); }
    
    /**
     * @brief Run service until shutdown signal
     * Blocks until stop() is called or shutdown message received
     */
    void run();
    
    /**
     * @brief Calculate field synchronously
     * @param request Field calculation request
     * @return Field calculation response
     */
    FieldResponse calculate_field(const FieldRequest& request);
    
    /**
     * @brief Get service statistics
     */
    struct Stats {
        uint64_t requests_processed = 0;
        uint64_t requests_failed = 0;
        double total_compute_time_ms = 0;
        double avg_compute_time_ms = 0;
    };
    Stats get_stats() const;
    
    // Handler registration for custom physics
    using FieldHandler = std::function<FieldResponse(const FieldRequest&)>;
    void register_system_handler(const std::string& system_name, FieldHandler handler);
    
    // ========================================================================
    // SELF-EXPAND: Dynamic Term Registration (v3.1)
    // ========================================================================
    
    /**
     * @brief Register a new physics term at runtime (from VR or external)
     * @param term The dynamic term to register
     * @return true if registered successfully
     */
    bool onRegisterTerm(const DynamicTerm& term);
    
    /**
     * @brief Unregister a dynamic term by name
     * @param name Term name to remove
     * @return true if found and removed
     */
    bool unregisterTerm(const std::string& name);
    
    /**
     * @brief Get all registered dynamic terms
     * @return Map of term name -> DynamicTerm
     */
    std::map<std::string, DynamicTerm> getDynamicTerms() const;
    
    /**
     * @brief Calculate total contribution from dynamic terms at (r, t)
     */
    double evaluateDynamicTerms(double r, double t) const;
    
    // Callback for term registration events
    using TermCallback = std::function<void(const DynamicTerm&, bool registered)>;
    void setTermCallback(TermCallback callback) { term_callback_ = std::move(callback); }
    
    // ========================================================================
    // SELF-UPDATE: Runtime Parameter Tuning (v3.1)
    // ========================================================================
    
    /**
     * @brief Update a calibrated parameter at runtime
     * @param param_name Parameter name (kappa, SSq, beta_i, etc.)
     * @param value New value
     * @param source Update source ("VR", "optimization", "user")
     * @return true if parameter exists and was updated
     */
    bool onUpdateParameter(const std::string& param_name, double value, 
                          const std::string& source = "external");
    
    /**
     * @brief Get current value of a calibrated parameter
     */
    double getParameter(const std::string& param_name) const;
    
    /**
     * @brief Get all calibrated parameters
     */
    CalibratedParameters getParameters() const;
    
    /**
     * @brief Reset parameters to validated defaults
     */
    void resetParameters();
    
    // Callback for parameter updates
    using ParameterCallback = std::function<void(const std::string& name, double old_val, double new_val)>;
    void setParameterCallback(ParameterCallback callback) { param_callback_ = std::move(callback); }
    
    // ========================================================================
    // SELF-SIMULATE: Time Evolution (v3.1)
    // ========================================================================
    
    /**
     * @brief Start a time evolution simulation
     * @param config Simulation configuration
     * @return Simulation ID (for tracking/cancellation)
     */
    uint64_t startSimulation(const SimulationConfig& config);
    
    /**
     * @brief Stop an active simulation
     * @param sim_id Simulation ID from startSimulation
     */
    void stopSimulation(uint64_t sim_id);
    
    /**
     * @brief Check if simulation is running
     */
    bool isSimulationRunning(uint64_t sim_id) const;
    
    /**
     * @brief Get simulation progress (0-1)
     */
    double getSimulationProgress(uint64_t sim_id) const;
    
    // Callback for simulation frame output (streamed to VR)
    using FrameCallback = std::function<void(const SimulationFrame&)>;
    void setFrameCallback(FrameCallback callback) { frame_callback_ = std::move(callback); }
    
private:
    ServiceConfig config_;
    std::atomic<bool> running_{false};
    std::atomic<bool> shutdown_requested_{false};
    
    // IPC channels
    std::unique_ptr<IPC::SharedMemoryChannel> shm_channel_;
    std::unique_ptr<IPC::GrpcChannel> grpc_channel_;
    IPC::MessageDispatcher dispatcher_;
    
    // Worker threads
    std::vector<std::thread> workers_;
    
    // Statistics
    mutable std::mutex stats_mutex_;
    Stats stats_;
    
    // Custom handlers
    std::unordered_map<std::string, FieldHandler> system_handlers_;
    
    // ========== Self-Expand: Dynamic Terms ==========
    mutable std::shared_mutex terms_mutex_;
    std::map<std::string, DynamicTerm> dynamic_terms_;
    TermCallback term_callback_;
    
    // ========== Self-Update: Calibrated Parameters ==========
    mutable std::shared_mutex params_mutex_;
    CalibratedParameters calibrated_params_;
    ParameterCallback param_callback_;
    
    // ========== Self-Simulate: Time Evolution ==========
    struct SimulationState {
        SimulationConfig config;
        std::atomic<bool> running{false};
        std::atomic<double> progress{0.0};
        std::thread thread;
        uint64_t start_time = 0;
    };
    mutable std::mutex sim_mutex_;
    std::map<uint64_t, std::unique_ptr<SimulationState>> simulations_;
    std::atomic<uint64_t> next_sim_id_{1};
    FrameCallback frame_callback_;
    
    // Simulation worker
    void simulation_loop(uint64_t sim_id);
    
    // Internal methods
    void worker_loop();
    void handle_calculate_field(const IPC::MessageHeader& header, 
                               const std::vector<uint8_t>& payload);
    void handle_register_term(const IPC::MessageHeader& header,
                             const std::vector<uint8_t>& payload);
    void handle_update_parameter(const IPC::MessageHeader& header,
                                const std::vector<uint8_t>& payload);
    void handle_vr_frame_update(const IPC::MessageHeader& header,
                               const std::vector<uint8_t>& payload);
    void handle_start_simulation(const IPC::MessageHeader& header,
                                const std::vector<uint8_t>& payload);
    
    // Physics calculations (delegates to MAIN_1_CoAnQi or Python)
    FieldResponse calculate_uqff(const FieldRequest& request);
    FieldResponse calculate_muge(const FieldRequest& request);
};

/**
 * @brief Run physics service from command line
 * 
 * Called from main() when --service flag is passed.
 * Creates QCoreApplication (no GUI) for event loop.
 * 
 * @param argc Argument count
 * @param argv Argument values
 * @return Exit code (0 = normal shutdown)
 */
int run_physics_service(int argc, char* argv[]);

/**
 * @brief Check if running in service mode
 * @param argc Argument count
 * @param argv Argument values
 * @return true if --service flag present
 */
bool is_service_mode(int argc, char* argv[]);

} // namespace UQFF

#endif // PHYSICS_SERVICE_H
