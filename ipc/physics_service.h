/**
 * @file physics_service.h
 * @brief Physics Backend Service for UQFF Star-Magic
 * 
 * Provides headless service mode for source2.cpp, enabling:
 * - IPC server accepting physics calculation requests
 * - VR runtime field data streaming via shared memory
 * - gRPC endpoints for structured commands
 * - Integration with MAIN_1_CoAnQi physics engine
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 2 - Physics Backend Service Mode
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
    
    // gRPC Settings (Phase 2 stub, full in Phase 3)
    std::string grpc_address = "localhost:50051";
    bool enable_grpc = false;  // Disabled until Phase 3
    
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
