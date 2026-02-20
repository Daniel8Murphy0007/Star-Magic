/**
 * @file physics_service.cpp
 * @brief Physics Backend Service Implementation
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 2 - Physics Backend Service Mode
 */

#include "physics_service.h"
#include <iostream>
#include <chrono>
#include <csignal>
#include <cstring>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace UQFF {

// ============================================================================
// SERVICE CONFIG
// ============================================================================

ServiceConfig ServiceConfig::from_args(int argc, char* argv[]) {
    ServiceConfig config;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--verbose" || arg == "-v") {
            config.verbose = true;
        }
        else if (arg == "--channel" && i + 1 < argc) {
            config.ipc_channel_name = argv[++i];
        }
        else if (arg == "--buffer-size" && i + 1 < argc) {
            config.shm_buffer_size = std::stoull(argv[++i]);
        }
        else if (arg == "--threads" && i + 1 < argc) {
            config.worker_threads = std::stoi(argv[++i]);
        }
        else if (arg == "--grpc" && i + 1 < argc) {
            config.grpc_address = argv[++i];
            config.enable_grpc = true;
        }
        else if (arg == "--log" && i + 1 < argc) {
            config.log_file = argv[++i];
        }
        else if (arg == "--no-main1") {
            config.enable_main1_integration = false;
        }
        else if (arg == "--no-python") {
            config.enable_python_integration = false;
        }
    }
    
    return config;
}

// ============================================================================
// PHYSICS SERVICE IMPLEMENTATION
// ============================================================================

PhysicsService::PhysicsService(const ServiceConfig& config)
    : config_(config)
{
    // Register message handlers
    dispatcher_.register_handler(IPC::MessageType::CALCULATE_FIELD,
        [this](const IPC::MessageHeader& h, const std::vector<uint8_t>& p) {
            handle_calculate_field(h, p);
        });
    
    dispatcher_.register_handler(IPC::MessageType::REGISTER_TERM,
        [this](const IPC::MessageHeader& h, const std::vector<uint8_t>& p) {
            handle_register_term(h, p);
        });
    
    dispatcher_.register_handler(IPC::MessageType::UPDATE_PARAMETER,
        [this](const IPC::MessageHeader& h, const std::vector<uint8_t>& p) {
            handle_update_parameter(h, p);
        });
    
    dispatcher_.register_handler(IPC::MessageType::VR_FRAME_UPDATE,
        [this](const IPC::MessageHeader& h, const std::vector<uint8_t>& p) {
            handle_vr_frame_update(h, p);
        });
    
    if (config_.verbose) {
        std::cout << "[PhysicsService] Initialized with config:\n"
                  << "  Channel: " << config_.ipc_channel_name << "\n"
                  << "  Buffer: " << config_.shm_buffer_size / (1024*1024) << " MB\n"
                  << "  Threads: " << config_.worker_threads << "\n"
                  << "  gRPC: " << (config_.enable_grpc ? config_.grpc_address : "disabled") << "\n"
                  << std::endl;
    }
}

PhysicsService::~PhysicsService() {
    stop();
}

bool PhysicsService::start() {
    if (running_.load()) {
        std::cerr << "[PhysicsService] Already running" << std::endl;
        return false;
    }
    
    // Create shared memory channel (server mode = create)
    shm_channel_ = std::make_unique<IPC::SharedMemoryChannel>(
        config_.ipc_channel_name,
        config_.shm_buffer_size,
        true  // Create new
    );
    
    if (!shm_channel_->is_connected()) {
        std::cerr << "[PhysicsService] Failed to create shared memory channel" << std::endl;
        return false;
    }
    
    // Create gRPC channel if enabled (stub for Phase 2)
    if (config_.enable_grpc) {
        grpc_channel_ = std::make_unique<IPC::GrpcChannel>(config_.grpc_address);
        std::cout << "[PhysicsService] gRPC stub created (full impl in Phase 3)" << std::endl;
    }
    
    running_ = true;
    shutdown_requested_ = false;
    
    // Start worker threads
    if (config_.use_thread_pool) {
        for (int i = 0; i < config_.worker_threads; ++i) {
            workers_.emplace_back(&PhysicsService::worker_loop, this);
        }
    }
    
    std::cout << "[PhysicsService] Started with " << config_.worker_threads 
              << " worker threads" << std::endl;
    
    return true;
}

void PhysicsService::stop() {
    if (!running_.load()) return;
    
    std::cout << "[PhysicsService] Stopping..." << std::endl;
    
    shutdown_requested_ = true;
    running_ = false;
    
    // Wait for workers to finish
    for (auto& worker : workers_) {
        if (worker.joinable()) {
            worker.join();
        }
    }
    workers_.clear();
    
    // Close channels
    if (shm_channel_) {
        shm_channel_->close();
        shm_channel_.reset();
    }
    
    if (grpc_channel_) {
        grpc_channel_->close();
        grpc_channel_.reset();
    }
    
    // Print final stats
    auto stats = get_stats();
    std::cout << "[PhysicsService] Shutdown complete. Stats:\n"
              << "  Requests processed: " << stats.requests_processed << "\n"
              << "  Requests failed: " << stats.requests_failed << "\n"
              << "  Avg compute time: " << stats.avg_compute_time_ms << " ms\n"
              << std::endl;
}

void PhysicsService::run() {
    if (!running_.load()) {
        if (!start()) {
            return;
        }
    }
    
    std::cout << "[PhysicsService] Running (Ctrl+C to stop)..." << std::endl;
    
    // Main loop - process messages from shared memory
    IPC::MessageHeader header;
    std::vector<uint8_t> payload;
    
    while (running_.load() && !shutdown_requested_.load()) {
        if (shm_channel_ && shm_channel_->receive(header, payload, 100)) {
            // Handle shutdown message
            if (header.type == IPC::MessageType::SHUTDOWN) {
                std::cout << "[PhysicsService] Received shutdown command" << std::endl;
                break;
            }
            
            // Handle ping
            if (header.type == IPC::MessageType::PING) {
                IPC::MessageHeader pong(IPC::MessageType::PONG);
                shm_channel_->send(pong);
                continue;
            }
            
            // Dispatch to handler
            dispatcher_.dispatch(header, payload);
        }
    }
    
    stop();
}

void PhysicsService::worker_loop() {
    if (config_.verbose) {
        std::cout << "[PhysicsService] Worker thread started" << std::endl;
    }
    
    // Workers currently share the main message loop
    // In future: separate request queue per worker
    while (running_.load() && !shutdown_requested_.load()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    
    if (config_.verbose) {
        std::cout << "[PhysicsService] Worker thread stopped" << std::endl;
    }
}

// ============================================================================
// MESSAGE HANDLERS
// ============================================================================

void PhysicsService::handle_calculate_field(const IPC::MessageHeader& header,
                                            const std::vector<uint8_t>& payload) {
    auto start = std::chrono::high_resolution_clock::now();
    
    // Deserialize request
    FieldRequest request;
    if (payload.size() >= sizeof(IPC::CalculateFieldRequest)) {
        const auto* ipc_req = reinterpret_cast<const IPC::CalculateFieldRequest*>(payload.data());
        request.system_name = std::string(ipc_req->system_name, 
            strnlen(ipc_req->system_name, sizeof(ipc_req->system_name)));
        request.r = ipc_req->r;
        request.t = ipc_req->t;
        request.theta = ipc_req->theta;
        request.flags = ipc_req->flags;
    }
    
    // Check for custom handler
    FieldResponse response;
    auto handler_it = system_handlers_.find(request.system_name);
    if (handler_it != system_handlers_.end()) {
        response = handler_it->second(request);
    } else {
        // Use default UQFF calculation
        response = calculate_uqff(request);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    response.compute_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    // Update stats
    {
        std::lock_guard<std::mutex> lock(stats_mutex_);
        if (response.success) {
            stats_.requests_processed++;
        } else {
            stats_.requests_failed++;
        }
        stats_.total_compute_time_ms += response.compute_time_ms;
        stats_.avg_compute_time_ms = stats_.total_compute_time_ms / 
            (stats_.requests_processed + stats_.requests_failed);
    }
    
    // Send response
    IPC::CalculateFieldResponse ipc_resp;
    ipc_resp.status = response.success ? 0 : 1;  // 0 = success
    ipc_resp.F_U = response.F_U;
    ipc_resp.Ug1 = response.Ug1;
    ipc_resp.Ug2 = response.Ug2;
    ipc_resp.Ug3 = response.Ug3;
    ipc_resp.Ug4 = response.Ug4;
    ipc_resp.Um = response.Um;
    ipc_resp.Ub = response.Ubi;                  // Ubi -> Ub
    ipc_resp.compressed_g = response.g_compressed;
    ipc_resp.reserved = 0;
    
    IPC::MessageHeader resp_header(IPC::MessageType::CALCULATE_FIELD);
    resp_header.payload_size = sizeof(ipc_resp);
    resp_header.flags = 0x01;  // Response flag
    
    if (shm_channel_) {
        shm_channel_->send(resp_header, &ipc_resp);
    }
    
    if (config_.verbose) {
        std::cout << "[PhysicsService] CALCULATE_FIELD: " << request.system_name
                  << " r=" << request.r << " -> F_U=" << response.F_U
                  << " (" << response.compute_time_ms << " ms)" << std::endl;
    }
}

void PhysicsService::handle_register_term(const IPC::MessageHeader& header,
                                          const std::vector<uint8_t>& payload) {
    // TODO: Implement dynamic term registration
    // This will integrate with MAIN_1_CoAnQi's PhysicsTermRegistry
    if (config_.verbose) {
        std::cout << "[PhysicsService] REGISTER_TERM received (not yet implemented)" << std::endl;
    }
}

void PhysicsService::handle_update_parameter(const IPC::MessageHeader& header,
                                             const std::vector<uint8_t>& payload) {
    // TODO: Implement runtime parameter updates
    if (config_.verbose) {
        std::cout << "[PhysicsService] UPDATE_PARAMETER received (not yet implemented)" << std::endl;
    }
}

void PhysicsService::handle_vr_frame_update(const IPC::MessageHeader& header,
                                            const std::vector<uint8_t>& payload) {
    // VR frame update - high-priority field recalculation for current viewpoint
    // Will be fully implemented in Phase 3
    if (config_.verbose) {
        std::cout << "[PhysicsService] VR_FRAME_UPDATE received" << std::endl;
    }
}

// ============================================================================
// PHYSICS CALCULATIONS
// ============================================================================

FieldResponse PhysicsService::calculate_uqff(const FieldRequest& request) {
    FieldResponse response;
    response.success = true;
    
    // TODO: Integrate with MAIN_1_CoAnQi physics engine
    // For now, use simplified UQFF approximation
    
    const double G = 6.67430e-11;
    const double c = 2.99792458e8;
    const double M_sun = 1.989e30;
    
    // Default to solar mass if not specified
    double M = request.mass > 0 ? request.mass : M_sun;
    double r = request.r > 0 ? request.r : 1e6;  // Default 1000 km
    
    // Simplified Ug components (placeholder - real calc from MAIN_1_CoAnQi)
    response.Ug1 = (G * M) / (r * r);                    // Newtonian base
    response.Ug2 = response.Ug1 * 1e-6;                  // Charge-reactivity (small)
    response.Ug3 = response.Ug1 * std::sin(request.theta) * 1e-3;  // Angular dependence
    response.Ug4 = response.Ug1 * 1e-12;                 // Vacuum (very small)
    
    // Buoyancy opposition (β_i ≈ 0.603 calibrated)
    const double beta_i = 0.603;
    response.Ubi = -beta_i * (response.Ug1 + response.Ug2 + response.Ug3 + response.Ug4);
    
    // Magnetism (placeholder)
    response.Um = 0;
    
    // Total unified field
    response.F_U = response.Ug1 + response.Ug2 + response.Ug3 + response.Ug4 + 
                   response.Um + response.Ubi;
    
    // MUGE compressed gravity
    response.g_compressed = response.Ug1;  // Simplified
    
    // Validation
    response.residual = 0.0;  // Unknown without observation
    response.confidence = 0.85;  // Moderate confidence for placeholder
    
    return response;
}

FieldResponse PhysicsService::calculate_muge(const FieldRequest& request) {
    FieldResponse response;
    response.success = true;
    
    // TODO: Full MUGE implementation from SOURCE4
    // For now, delegate to calculate_uqff
    return calculate_uqff(request);
}

// ============================================================================
// STATISTICS
// ============================================================================

PhysicsService::Stats PhysicsService::get_stats() const {
    std::lock_guard<std::mutex> lock(stats_mutex_);
    return stats_;
}

void PhysicsService::register_system_handler(const std::string& system_name, 
                                             FieldHandler handler) {
    system_handlers_[system_name] = std::move(handler);
    if (config_.verbose) {
        std::cout << "[PhysicsService] Registered handler for: " << system_name << std::endl;
    }
}

// ============================================================================
// SERVICE MODE ENTRY POINTS
// ============================================================================

// Global signal handler for graceful shutdown
static std::atomic<bool> g_shutdown_requested{false};

#ifdef _WIN32
BOOL WINAPI console_handler(DWORD signal) {
    if (signal == CTRL_C_EVENT || signal == CTRL_BREAK_EVENT) {
        g_shutdown_requested = true;
        return TRUE;
    }
    return FALSE;
}
#else
void signal_handler(int signal) {
    if (signal == SIGINT || signal == SIGTERM) {
        g_shutdown_requested = true;
    }
}
#endif

bool is_service_mode(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--service" || arg == "-s") {
            return true;
        }
    }
    return false;
}

int run_physics_service(int argc, char* argv[]) {
    std::cout << "======================================\n"
              << "  UQFF Physics Backend Service\n"
              << "  Star-Magic v3.0 - Phase 2\n"
              << "======================================\n" << std::endl;
    
    // Install signal handlers
#ifdef _WIN32
    SetConsoleCtrlHandler(console_handler, TRUE);
#else
    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);
#endif
    
    // Parse configuration
    ServiceConfig config = ServiceConfig::from_args(argc, argv);
    
    // Create and run service
    PhysicsService service(config);
    
    if (!service.start()) {
        std::cerr << "[ERROR] Failed to start physics service" << std::endl;
        return 1;
    }
    
    // Run until shutdown
    while (!g_shutdown_requested.load() && service.is_running()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        
        // Check for messages in main thread as well
        // (workers handle computation, main handles control)
    }
    
    service.stop();
    
    return 0;
}

} // namespace UQFF
