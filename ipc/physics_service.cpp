/**
 * @file physics_service.cpp
 * @brief Physics Backend Service Implementation
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.1
 * Phase: 3.1 - Self-Expanding Physics Backend
 * 
 * v3.1 Features:
 * - Self-Expand: onRegisterTerm() for dynamic physics term registration
 * - Self-Update: onUpdateParameter() for runtime κ, [SSq], β_i tuning
 * - Self-Simulate: Time evolution with VR streaming
 */

#include "physics_service.h"
#include "python_bridge.h"
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
    
    // v3.1: Simulation control handler
    dispatcher_.register_handler(IPC::MessageType::SIM_START,
        [this](const IPC::MessageHeader& h, const std::vector<uint8_t>& p) {
            handle_start_simulation(h, p);
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
    
    // Create gRPC channel if enabled (Phase 3 - full implementation)
    if (config_.enable_grpc) {
        grpc_channel_ = std::make_unique<IPC::GrpcChannel>(config_.grpc_address);
        if (grpc_channel_->is_connected()) {
            std::cout << "[PhysicsService] gRPC channel ready: " << config_.grpc_address << std::endl;
        } else {
            std::cout << "[PhysicsService] gRPC channel created (pending connection): " 
                      << config_.grpc_address << std::endl;
        }
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

void PhysicsService::handle_vr_frame_update(const IPC::MessageHeader& header,
                                            const std::vector<uint8_t>& payload) {
    // VR frame update - high-priority field recalculation for current viewpoint
    if (payload.size() >= sizeof(IPC::VRFrameUpdate)) {
        const auto* update = reinterpret_cast<const IPC::VRFrameUpdate*>(payload.data());
        
        // Calculate field at probe position
        FieldRequest request;
        request.r = std::sqrt(update->field_probe_position[0] * update->field_probe_position[0] +
                              update->field_probe_position[1] * update->field_probe_position[1] +
                              update->field_probe_position[2] * update->field_probe_position[2]);
        request.t = static_cast<double>(update->frame_number) / 90.0;  // Assume 90 FPS
        request.theta = std::atan2(update->field_probe_position[1], 
                                    update->field_probe_position[0]);
        
        FieldResponse response = calculate_uqff(request);
        
        // Add dynamic term contributions (v3.1)
        response.F_U += evaluateDynamicTerms(request.r, request.t);
        
        // Send response back to VR
        if (shm_channel_) {
            IPC::CalculateFieldResponse ipc_resp;
            ipc_resp.status = 0;
            ipc_resp.F_U = response.F_U;
            ipc_resp.Ug1 = response.Ug1;
            ipc_resp.Ug2 = response.Ug2;
            ipc_resp.Ug3 = response.Ug3;
            ipc_resp.Ug4 = response.Ug4;
            ipc_resp.Um = response.Um;
            ipc_resp.Ub = response.Ubi;
            ipc_resp.compressed_g = response.g_compressed;
            
            IPC::MessageHeader resp_header(IPC::MessageType::VR_FRAME_UPDATE);
            resp_header.payload_size = sizeof(ipc_resp);
            resp_header.flags = 0x01;
            shm_channel_->try_send(resp_header, &ipc_resp);
        }
    }
    
    if (config_.verbose) {
        std::cout << "[PhysicsService] VR_FRAME_UPDATE processed" << std::endl;
    }
}

// ============================================================================
// PHYSICS CALCULATIONS (integrated with Python calculators)
// ============================================================================

FieldResponse PhysicsService::calculate_uqff(const FieldRequest& request) {
    FieldResponse response;
    
    // Try Python calculators first (CondensedPhysics.py → QCalc.py fallback)
    if (config_.enable_python_integration) {
        PythonBridge& python = PythonBridge::instance();
        
        if (!python.is_ready()) {
            python.initialize();
        }
        
        if (python.is_ready()) {
            PythonQueryParams params;
            params.system_name = request.system_name;
            params.r = request.r;
            params.t = request.t;
            params.theta = request.theta;
            params.mass = request.mass;
            params.radius = request.radius;
            params.magnetic_field = request.magnetic_field;
            params.spin_period = request.spin_period;
            
            // Try CondensedPhysics.py first (81K lines, full equations)
            PythonFieldResult py_result = python.calculate_condensed_physics(params);
            
            if (!py_result.success) {
                // Fallback to QCalc.py (9K lines, pure solver)
                py_result = python.calculate_qcalc(params);
            }
            
            if (py_result.success) {
                response.success = true;
                response.F_U = py_result.F_U;
                response.Ug1 = py_result.Ug1;
                response.Ug2 = py_result.Ug2;
                response.Ug3 = py_result.Ug3;
                response.Ug4 = py_result.Ug4;
                response.Um = py_result.Um;
                response.Ubi = py_result.Ubi;
                response.g_compressed = py_result.g_compressed;
                response.compute_time_ms = py_result.compute_time_ms;
                response.confidence = 0.999;  // 99.9% solvability
                
                if (config_.verbose) {
                    std::cout << "[PhysicsService] Python calculation: F_U=" << response.F_U
                              << " (" << py_result.compute_time_ms << " ms)" << std::endl;
                }
                
                return response;
            } else {
                std::cerr << "[PhysicsService] Python error: " << py_result.error_message << std::endl;
            }
        }
    }
    
    // Fallback: C++ simplified UQFF calculation (DPM-based, matches source4.cpp canonical)
    // IMPORTANT: Newtonian gravity is EMERGENT from DPM substrate, not foundational.
    // Order: DPM → Ug2 (quantum shell traps magnetics) → mass emerges → Ug1 (emergent)
    response.success = true;
    
    const double PI = 3.14159265358979323846;
    const double rho_A = 7.09e-37;   // [SCm] vacuum density [J/m³]
    const double rho_UA = 7.09e-36;  // [UA] vacuum density [J/m³]
    
    double R = request.radius > 0 ? request.radius : 6.96e8;  // Body radius [m]
    double r = request.r > 0 ? request.r : 1e6;
    double t = request.t;
    double B = request.magnetic_field > 0 ? request.magnetic_field : 1e-4;  // [T]
    
    // Get calibrated parameters (v3.1 - runtime tunable)
    CalibratedParameters params = getParameters();
    
    // ===== DPM-based Ug components (canonical source4.cpp form) =====
    
    // Ug1: Magnetic dipole rotation (emergent gravity from magnetic moment)
    // μ_s = ρ_A × V_body (magnetic moment from [SCm] vacuum density)
    // Ug1 = k1 × μ_s × ∇(M_s/r) × decay × oscillation × deformation
    double V_body = (4.0 / 3.0) * PI * R * R * R;
    double mu_s = rho_A * V_body;               // Magnetic moment [J/T]
    double M = request.mass > 0 ? request.mass : 1.989e30;
    double grad_M = M / (r * r);                // Gradient of M_s/r [kg/m²]
    double alpha = params.kappa > 0 ? params.kappa : 0.0005;
    double exp_decay = std::exp(-alpha * t);
    double tn = std::fmod(t, 1.0);              // Normalized cycle time
    double cos_ptn = std::cos(PI * tn);
    response.Ug1 = mu_s * grad_M * exp_decay * cos_ptn * 1.1;  // 1.1 = deformation factor
    
    // Ug2: Charge-reactivity coupling (quantum shell trapping magnetics)
    // Dual charges Q_SCm and Q_UA from vacuum densities
    double Q_SCm = rho_A * V_body;
    double Q_UA = rho_UA * V_body;
    double v_sw = 4e5;                           // Solar wind velocity [m/s]
    double E_react = rho_A * v_sw * v_sw / rho_UA * std::exp(-alpha * t);
    double R_b = R * 100.0;                      // Bubble radius
    double S_rb = (r > R_b) ? 1.0 : 0.0;        // Heliosphere step function
    double sw_factor = 1.0 + 0.01 * v_sw;       // Solar wind enhancement
    response.Ug2 = (Q_SCm + Q_UA) * M / (r * r) * S_rb * sw_factor * params.H_SCm * E_react;
    
    // Ug3: Magnetic string rotation (B_disk × cos(ω_s t π) × P_core × E_react)
    double omega_s = 2.0 * PI / (25.0 * 86400.0);  // Solar rotation rate
    response.Ug3 = B * std::cos(omega_s * t * PI) * E_react * params.SSq;
    
    // Ug4: Vacuum concentration (ρ_vac × C_concentration × decay × oscillation)
    double rho_v = 7.09e-37;                     // [SCm] vacuum density
    double C_concentration = 1e30;               // Concentration factor
    response.Ug4 = rho_v * C_concentration * exp_decay * cos_ptn * params.k_eta;
    
    // Buoyancy opposition (β_i × Σ Ug × galactic factors)
    double sum_Ug = response.Ug1 + response.Ug2 + response.Ug3 + response.Ug4;
    response.Ubi = -params.beta_i * sum_Ug * rho_A * cos_ptn;
    
    // Magnetism: μ_j/r³ dipole field with H_SCm superconducting factor
    double mu_j = B * R * R * R;
    response.Um = mu_j / (r * r * r) * params.H_SCm;
    
    // Total unified field
    response.F_U = sum_Ug + response.Ubi + response.Um;
    
    // Apply time decay factor κ
    if (request.t > 0 && params.kappa > 0) {
        double decay = std::exp(-params.kappa * request.t / 86400.0);  // κ is per day
        response.F_U *= decay;
    }
    
    // MUGE compressed gravity: DPM foundation (a_DPM = F_DPM × f_DPM × E_vac,neb / (c × V_sys))
    double c_light = 2.998e8;
    double V_sys = V_body * 1e6;                 // System volume estimate
    double F_DPM = B * V_body * omega_s;         // DPM force proxy
    double f_DPM = 1.0;                          // DPM frequency modulation
    double E_vac_neb = rho_A * V_sys;            // Vacuum energy in nebula
    response.g_compressed = F_DPM * f_DPM * E_vac_neb / (c_light * V_sys);
    
    // Validation
    response.residual = 0.0;  // Unknown without observation
    response.confidence = 0.999;  // 99.9% solvability
    
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
// SELF-EXPAND: Dynamic Term Registration (v3.1)
// ============================================================================

bool PhysicsService::onRegisterTerm(const DynamicTerm& term) {
    std::unique_lock<std::shared_mutex> lock(terms_mutex_);
    
    // Validate term
    if (term.name.empty()) {
        std::cerr << "[PhysicsService] Cannot register term: empty name" << std::endl;
        return false;
    }
    
    // Check for duplicate
    bool is_update = (dynamic_terms_.find(term.name) != dynamic_terms_.end());
    
    // Register/update term
    DynamicTerm stored_term = term;
    stored_term.registered_at = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    
    dynamic_terms_[term.name] = stored_term;
    
    if (config_.verbose) {
        std::cout << "[PhysicsService] " << (is_update ? "Updated" : "Registered") 
                  << " dynamic term: " << term.name 
                  << " (type=" << static_cast<uint32_t>(term.type) 
                  << ", coeff=" << term.coefficient << ")" << std::endl;
    }
    
    // Notify callback
    if (term_callback_) {
        lock.unlock();  // Release lock before callback
        term_callback_(stored_term, true);
    }
    
    return true;
}

bool PhysicsService::unregisterTerm(const std::string& name) {
    std::unique_lock<std::shared_mutex> lock(terms_mutex_);
    
    auto it = dynamic_terms_.find(name);
    if (it == dynamic_terms_.end()) {
        return false;
    }
    
    DynamicTerm removed_term = it->second;
    dynamic_terms_.erase(it);
    
    if (config_.verbose) {
        std::cout << "[PhysicsService] Unregistered dynamic term: " << name << std::endl;
    }
    
    // Notify callback
    if (term_callback_) {
        lock.unlock();
        term_callback_(removed_term, false);
    }
    
    return true;
}

std::map<std::string, DynamicTerm> PhysicsService::getDynamicTerms() const {
    std::shared_lock<std::shared_mutex> lock(terms_mutex_);
    return dynamic_terms_;
}

double PhysicsService::evaluateDynamicTerms(double r, double t) const {
    std::shared_lock<std::shared_mutex> lock(terms_mutex_);
    
    double total = 0.0;
    for (const auto& [name, term] : dynamic_terms_) {
        total += term.evaluate(r, t);
    }
    return total;
}

// ============================================================================
// SELF-UPDATE: Runtime Parameter Tuning (v3.1)
// ============================================================================

bool PhysicsService::onUpdateParameter(const std::string& param_name, double value,
                                        const std::string& source) {
    std::unique_lock<std::shared_mutex> lock(params_mutex_);
    
    double old_value = 0.0;
    bool found = true;
    
    // Map parameter name to struct field
    if (param_name == "kappa" || param_name == "κ") {
        old_value = calibrated_params_.kappa;
        calibrated_params_.kappa = value;
    } else if (param_name == "SSq" || param_name == "[SSq]") {
        old_value = calibrated_params_.SSq;
        calibrated_params_.SSq = value;
    } else if (param_name == "H_SCm") {
        old_value = calibrated_params_.H_SCm;
        calibrated_params_.H_SCm = value;
    } else if (param_name == "U_UA") {
        old_value = calibrated_params_.U_UA;
        calibrated_params_.U_UA = value;
    } else if (param_name == "k_eta" || param_name == "k_η") {
        old_value = calibrated_params_.k_eta;
        calibrated_params_.k_eta = value;
    } else if (param_name == "beta_i" || param_name == "β_i") {
        old_value = calibrated_params_.beta_i;
        calibrated_params_.beta_i = value;
    } else if (param_name == "alpha_DPM") {
        old_value = calibrated_params_.alpha_DPM;
        calibrated_params_.alpha_DPM = value;
    } else if (param_name == "Omega_Lambda" || param_name == "Ω_Λ") {
        old_value = calibrated_params_.Omega_Lambda;
        calibrated_params_.Omega_Lambda = value;
    } else if (param_name == "H_0") {
        old_value = calibrated_params_.H_0;
        calibrated_params_.H_0 = value;
    } else if (param_name == "learning_rate") {
        old_value = calibrated_params_.learning_rate;
        calibrated_params_.learning_rate = value;
    } else if (param_name == "active_layers") {
        old_value = static_cast<double>(calibrated_params_.active_layers);
        calibrated_params_.active_layers = static_cast<int>(value);
    } else {
        found = false;
    }
    
    if (!found) {
        std::cerr << "[PhysicsService] Unknown parameter: " << param_name << std::endl;
        return false;
    }
    
    calibrated_params_.last_updated = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    calibrated_params_.update_source = source;
    
    if (config_.verbose) {
        std::cout << "[PhysicsService] Updated parameter: " << param_name 
                  << " = " << value << " (was " << old_value << ")"
                  << " [source: " << source << "]" << std::endl;
    }
    
    // Notify callback
    if (param_callback_) {
        lock.unlock();
        param_callback_(param_name, old_value, value);
    }
    
    return true;
}

double PhysicsService::getParameter(const std::string& param_name) const {
    std::shared_lock<std::shared_mutex> lock(params_mutex_);
    
    if (param_name == "kappa" || param_name == "κ") return calibrated_params_.kappa;
    if (param_name == "SSq" || param_name == "[SSq]") return calibrated_params_.SSq;
    if (param_name == "H_SCm") return calibrated_params_.H_SCm;
    if (param_name == "U_UA") return calibrated_params_.U_UA;
    if (param_name == "k_eta" || param_name == "k_η") return calibrated_params_.k_eta;
    if (param_name == "beta_i" || param_name == "β_i") return calibrated_params_.beta_i;
    if (param_name == "alpha_DPM") return calibrated_params_.alpha_DPM;
    if (param_name == "Omega_Lambda" || param_name == "Ω_Λ") return calibrated_params_.Omega_Lambda;
    if (param_name == "H_0") return calibrated_params_.H_0;
    if (param_name == "learning_rate") return calibrated_params_.learning_rate;
    if (param_name == "active_layers") return static_cast<double>(calibrated_params_.active_layers);
    
    return std::nan("");  // Unknown parameter
}

CalibratedParameters PhysicsService::getParameters() const {
    std::shared_lock<std::shared_mutex> lock(params_mutex_);
    return calibrated_params_;
}

void PhysicsService::resetParameters() {
    std::unique_lock<std::shared_mutex> lock(params_mutex_);
    calibrated_params_ = CalibratedParameters();  // Reset to defaults
    
    if (config_.verbose) {
        std::cout << "[PhysicsService] Parameters reset to validated defaults" << std::endl;
    }
}

// ============================================================================
// SELF-SIMULATE: Time Evolution (v3.1)
// ============================================================================

uint64_t PhysicsService::startSimulation(const SimulationConfig& config) {
    std::lock_guard<std::mutex> lock(sim_mutex_);
    
    uint64_t sim_id = next_sim_id_++;
    
    auto state = std::make_unique<SimulationState>();
    state->config = config;
    state->running = true;
    state->progress = 0.0;
    state->start_time = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    
    // Start simulation thread
    SimulationState* state_ptr = state.get();
    state->thread = std::thread([this, sim_id]() {
        simulation_loop(sim_id);
    });
    
    simulations_[sim_id] = std::move(state);
    
    if (config_.verbose) {
        std::cout << "[PhysicsService] Started simulation #" << sim_id 
                  << " for " << config.system_name
                  << " (" << config.frames << " frames)" << std::endl;
    }
    
    return sim_id;
}

void PhysicsService::stopSimulation(uint64_t sim_id) {
    std::unique_lock<std::mutex> lock(sim_mutex_);
    
    auto it = simulations_.find(sim_id);
    if (it == simulations_.end()) return;
    
    it->second->running = false;
    
    // Join the thread
    lock.unlock();
    if (it->second->thread.joinable()) {
        it->second->thread.join();
    }
    
    lock.lock();
    simulations_.erase(it);
    
    if (config_.verbose) {
        std::cout << "[PhysicsService] Stopped simulation #" << sim_id << std::endl;
    }
}

bool PhysicsService::isSimulationRunning(uint64_t sim_id) const {
    std::lock_guard<std::mutex> lock(sim_mutex_);
    
    auto it = simulations_.find(sim_id);
    if (it == simulations_.end()) return false;
    return it->second->running.load();
}

double PhysicsService::getSimulationProgress(uint64_t sim_id) const {
    std::lock_guard<std::mutex> lock(sim_mutex_);
    
    auto it = simulations_.find(sim_id);
    if (it == simulations_.end()) return 1.0;  // Completed
    return it->second->progress.load();
}

void PhysicsService::simulation_loop(uint64_t sim_id) {
    SimulationState* state = nullptr;
    {
        std::lock_guard<std::mutex> lock(sim_mutex_);
        auto it = simulations_.find(sim_id);
        if (it == simulations_.end()) return;
        state = it->second.get();
    }
    
    const SimulationConfig& config = state->config;
    
    double t = config.t_start;
    double r = config.r_start;
    double theta = 0.0;
    uint64_t frame_number = 0;
    
    double prev_F_U = 0.0;
    
    int total_steps = config.total_steps();
    int steps_per_frame = config.steps_per_frame();
    int step = 0;
    
    while (state->running.load() && t <= config.t_end) {
        // Calculate field at current (r, t)
        FieldRequest request;
        request.system_name = config.system_name;
        request.r = r;
        request.t = t;
        request.theta = theta;
        
        FieldResponse response = calculate_uqff(request);
        
        // Build frame
        SimulationFrame frame;
        frame.frame_number = frame_number;
        frame.t = t;
        frame.r = r;
        frame.theta = theta;
        frame.F_U = response.F_U;
        frame.Ug1 = response.Ug1;
        frame.Ug2 = response.Ug2;
        frame.Ug3 = response.Ug3;
        frame.Ug4 = response.Ug4;
        frame.Um = response.Um;
        frame.Ubi = response.Ubi;
        frame.g_compressed = response.g_compressed;
        
        // Add dynamic term contributions
        frame.F_U += evaluateDynamicTerms(r, t);
        
        // Derived quantities
        frame.dF_dt = (step > 0) ? (frame.F_U - prev_F_U) / config.dt : 0.0;
        frame.orbital_velocity = (frame.F_U > 0 && r > 0) ? std::sqrt(frame.F_U * r) : 0.0;
        frame.escape_velocity = frame.orbital_velocity * std::sqrt(2.0);
        
        // Progress
        frame.progress = static_cast<double>(step) / total_steps;
        frame.is_final = (step >= total_steps - 1);
        
        // Update state
        state->progress = frame.progress;
        prev_F_U = frame.F_U;
        
        // Output frame (every steps_per_frame)
        if (step % steps_per_frame == 0 || frame.is_final) {
            // Stream to VR via callback
            if (frame_callback_) {
                frame_callback_(frame);
            }
            
            // Stream to VR via IPC
            if (config.stream_to_vr && shm_channel_ && shm_channel_->is_connected()) {
                IPC::MessageHeader header(IPC::MessageType::SIM_FRAME);
                header.payload_size = sizeof(SimulationFrame);
                shm_channel_->try_send(header, &frame);
            }
            
            frame_number++;
        }
        
        // Advance time
        t += config.dt;
        if (config.dr != 0.0) {
            r += config.dr;
        }
        theta += 0.01;  // Slow rotation
        step++;
    }
    
    // Mark as complete
    state->running = false;
    state->progress = 1.0;
    
    // Send completion message
    if (shm_channel_ && shm_channel_->is_connected()) {
        IPC::MessageHeader header(IPC::MessageType::SIM_COMPLETE);
        shm_channel_->send(header);
    }
    
    if (config_.verbose) {
        std::cout << "[PhysicsService] Simulation #" << sim_id 
                  << " completed (" << frame_number << " frames)" << std::endl;
    }
}

// ============================================================================
// IPC MESSAGE HANDLERS (v3.1 updates)
// ============================================================================

void PhysicsService::handle_register_term(const IPC::MessageHeader& header,
                                          const std::vector<uint8_t>& payload) {
    if (payload.size() < sizeof(IPC::RegisterTermRequest)) {
        std::cerr << "[PhysicsService] REGISTER_TERM: payload too small" << std::endl;
        return;
    }
    
    const auto* req = reinterpret_cast<const IPC::RegisterTermRequest*>(payload.data());
    
    DynamicTerm term;
    term.name = std::string(req->term_name, strnlen(req->term_name, sizeof(req->term_name)));
    term.description = std::string(req->description, strnlen(req->description, sizeof(req->description)));
    term.base_value = req->initial_value;
    term.type = static_cast<DynamicTermType>(req->term_type);
    term.source = "IPC";
    
    bool success = onRegisterTerm(term);
    
    // Send response
    IPC::MessageHeader resp_header(success ? IPC::MessageType::RESPONSE_SUCCESS 
                                           : IPC::MessageType::RESPONSE_ERROR);
    if (shm_channel_) {
        shm_channel_->send(resp_header);
    }
}

void PhysicsService::handle_update_parameter(const IPC::MessageHeader& header,
                                             const std::vector<uint8_t>& payload) {
    if (payload.size() < sizeof(IPC::UpdateParameterRequest)) {
        std::cerr << "[PhysicsService] UPDATE_PARAMETER: payload too small" << std::endl;
        return;
    }
    
    const auto* req = reinterpret_cast<const IPC::UpdateParameterRequest*>(payload.data());
    
    std::string param_name(req->param_name, strnlen(req->param_name, sizeof(req->param_name)));
    
    bool success = onUpdateParameter(param_name, req->value, "IPC");
    
    // Send response
    IPC::MessageHeader resp_header(success ? IPC::MessageType::RESPONSE_SUCCESS 
                                           : IPC::MessageType::RESPONSE_ERROR);
    if (shm_channel_) {
        shm_channel_->send(resp_header);
    }
}

void PhysicsService::handle_start_simulation(const IPC::MessageHeader& header,
                                            const std::vector<uint8_t>& payload) {
    // Parse simulation config from payload
    SimulationConfig config;
    
    // For now, use basic struct serialization
    // TODO: Use proper protocol buffer or JSON serialization
    if (payload.size() >= 64) {
        const char* data = reinterpret_cast<const char*>(payload.data());
        config.system_name = std::string(data, strnlen(data, 64));
        
        if (payload.size() >= 64 + 7 * sizeof(double) + 2 * sizeof(int)) {
            const double* doubles = reinterpret_cast<const double*>(payload.data() + 64);
            config.r_start = doubles[0];
            config.r_end = doubles[1];
            config.t_start = doubles[2];
            config.t_end = doubles[3];
            config.dt = doubles[4];
            config.dr = doubles[5];
            
            const int* ints = reinterpret_cast<const int*>(doubles + 6);
            config.frames = ints[0];
            config.stream_to_vr = (ints[1] != 0);
        }
    } else {
        // Use defaults with system name from payload
        config.system_name = "SGR1745";
    }
    
    uint64_t sim_id = startSimulation(config);
    
    // Send response with simulation ID
    IPC::MessageHeader resp_header(IPC::MessageType::RESPONSE_DATA);
    resp_header.payload_size = sizeof(sim_id);
    if (shm_channel_) {
        shm_channel_->send(resp_header, &sim_id);
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
