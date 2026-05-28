/**
 * @file vr_runtime.cpp
 * @brief VR Runtime Layer Implementation for UQFF Star-Magic
 * 
 * Main VR compositor evolved from source2(HEAD PROGRAM).cpp
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 5 - Full VR Experience
 */

#include "vr_runtime.h"
#include "openxr_session.h"
#include "vulkan_compositor.h"
#include "task_bot.h"
#include "astro_graphics.h"
#include "../ipc/uqff_ipc.h"
#include "CoAnQi_bot.h"  // MAIN_1 exclusive specialist (rec #5 wiring)
// CoAnQi_bot hook (VR_ARCHITECTURE + PIPELINE): instantiate for MAIN_1 PhysicsTerm mgmt / self-expand
// Example (guarded, link with vr/CoAnQi_bot.cpp if building unified VR):
//   UQFF::VR::CoAnQiBot coanqiBot(".");
//   auto terms = coanqiBot.ListPhysicsTerms();  // 6,688+ including new S233 VDS/DVP/BH26
// Library extension: coanqiBot also owns whitepapers/PDF query surface (MAIN_1 menu Option 23).
// Per-frame AstroGraphics (FLOW_DIAGRAM v5.x): probe QCALCGEOM::vds_branches for VDS overlay on VR render

#include <iostream>
#include <sstream>
#include <chrono>
#include <cmath>

namespace VR {

// ============================================================================
// RuntimeConfig
// ============================================================================

RuntimeConfig RuntimeConfig::from_args(int argc, char* argv[]) {
    RuntimeConfig config;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--verbose" || arg == "-v") {
            config.verbose = true;
        }
        else if (arg == "--physics-channel" && i + 1 < argc) {
            config.physics_ipc_channel = argv[++i];
        }
        else if (arg == "--physics-port" && i + 1 < argc) {
            config.physics_port = static_cast<uint16_t>(std::stoi(argv[++i]));
        }
        else if (arg == "--width" && i + 1 < argc) {
            config.render_width = static_cast<uint32_t>(std::stoi(argv[++i]));
        }
        else if (arg == "--height" && i + 1 < argc) {
            config.render_height = static_cast<uint32_t>(std::stoi(argv[++i]));
        }
        else if (arg == "--refresh" && i + 1 < argc) {
            config.refresh_rate = std::stof(argv[++i]);
        }
        else if (arg == "--no-voice") {
            config.enable_voice_commands = false;
        }
        else if (arg == "--no-gesture") {
            config.enable_gesture_recognition = false;
        }
        else if (arg == "--astro-program" && i + 1 < argc) {
            config.astro_program_path = argv[++i];
        }
        else if (arg == "--perf-overlay") {
            config.enable_performance_overlay = true;
        }
        else if (arg == "--no-shm") {
            config.use_shared_memory = false;
        }
    }
    
    return config;
}

// ============================================================================
// VRRuntime Singleton
// ============================================================================

VRRuntime& VRRuntime::instance() {
    static VRRuntime instance;
    return instance;
}

VRRuntime::~VRRuntime() {
    shutdown();
}

// ============================================================================
// Lifecycle
// ============================================================================

bool VRRuntime::initialize(const RuntimeConfig& config) {
    if (state_ != RuntimeState::Uninitialized) {
        std::cerr << "VRRuntime: Already initialized" << std::endl;
        return false;
    }
    
    config_ = config;
    state_ = RuntimeState::Initializing;
    
    if (config_.verbose) {
        std::cout << "VRRuntime: Initializing..." << std::endl;
        std::cout << "  Application: " << config_.application_name << std::endl;
        std::cout << "  Resolution: " << config_.render_width << "x" << config_.render_height << std::endl;
        std::cout << "  Refresh: " << config_.refresh_rate << " Hz" << std::endl;
    }
    
    // Step 1: Initialize OpenXR
    if (!initOpenXR()) {
        std::cerr << "VRRuntime: Failed to initialize OpenXR" << std::endl;
        state_ = RuntimeState::Error;
        return false;
    }
    
    // Step 2: Initialize Vulkan compositor
    if (!initVulkan()) {
        std::cerr << "VRRuntime: Failed to initialize Vulkan compositor" << std::endl;
        state_ = RuntimeState::Error;
        return false;
    }
    
    // Step 3: Initialize Task Bot
    if (!initTaskBot()) {
        std::cerr << "VRRuntime: Failed to initialize Task Bot (non-fatal)" << std::endl;
        // Continue - Task Bot is optional
    }
    
    // Step 4: Connect to Physics Backend
    if (!initPhysicsConnection()) {
        std::cerr << "VRRuntime: Failed to connect to Physics Backend (non-fatal)" << std::endl;
        // Continue - can work without physics
    }
    
    // Step 5: Load Astro Graphics (if configured)
    if (!config_.astro_program_path.empty()) {
        if (!loadAstroProgram(config_.astro_program_path)) {
            std::cerr << "VRRuntime: Failed to load Astro Graphics (non-fatal)" << std::endl;
        }
    }
    
    state_ = RuntimeState::Ready;
    if (config_.verbose) {
        std::cout << "VRRuntime: Initialization complete" << std::endl;
    }
    
    return true;
}

void VRRuntime::shutdown() {
    if (state_ == RuntimeState::Uninitialized) {
        return;
    }
    
    state_ = RuntimeState::ShuttingDown;
    shutdown_requested_ = true;
    
    if (config_.verbose) {
        std::cout << "VRRuntime: Shutting down..." << std::endl;
    }
    
    // Wait for render thread
    if (render_thread_.joinable()) {
        render_thread_.join();
    }
    
    // Wait for physics thread
    if (physics_thread_.joinable()) {
        physics_thread_.join();
    }
    
    // Cleanup in reverse order
    astro_graphics_.reset();
    task_bot_.reset();
    vulkan_.reset();
    openxr_.reset();
    physics_channel_.reset();
    
    state_ = RuntimeState::Uninitialized;
    
    if (config_.verbose) {
        std::cout << "VRRuntime: Shutdown complete" << std::endl;
    }
}

// ============================================================================
// Initialization Steps
// ============================================================================

bool VRRuntime::initOpenXR() {
    openxr_ = std::make_unique<OpenXRSession>();
    
    if (!openxr_->initialize(config_.application_name, config_.engine_version)) {
        std::cerr << "VRRuntime: OpenXR initialization failed" << std::endl;
        std::cerr << "  Hint: Is an OpenXR runtime (SteamVR, Oculus, WMR) running?" << std::endl;
        return false;
    }
    
    if (config_.verbose) {
        std::cout << "  OpenXR: Initialized" << std::endl;
    }
    
    return true;
}

bool VRRuntime::initVulkan() {
    vulkan_ = std::make_unique<VulkanCompositor>();
    
    if (!vulkan_->initialize(openxr_.get())) {
        std::cerr << "VRRuntime: Vulkan compositor initialization failed" << std::endl;
        return false;
    }
    
    if (config_.verbose) {
        std::cout << "  Vulkan: Initialized" << std::endl;
    }
    
    return true;
}

bool VRRuntime::initTaskBot() {
    if (!config_.enable_voice_commands && !config_.enable_gesture_recognition) {
        if (config_.verbose) {
            std::cout << "  Task Bot: Disabled (no voice/gesture)" << std::endl;
        }
        return true;  // Not an error
    }
    
    task_bot_ = std::make_unique<TaskBot>();
    
    if (!task_bot_->initialize(config_.voice_model_path, config_.enable_gesture_recognition)) {
        std::cerr << "VRRuntime: Task Bot initialization failed" << std::endl;
        task_bot_.reset();
        return false;
    }
    
    // Register built-in commands
    task_bot_->registerBuiltinCommands(this);
    
    if (config_.verbose) {
        std::cout << "  Task Bot: Initialized" << std::endl;
        if (config_.enable_voice_commands) {
            std::cout << "    Voice: Enabled" << std::endl;
        }
        if (config_.enable_gesture_recognition) {
            std::cout << "    Gesture: Enabled" << std::endl;
        }
    }
    
    return true;
}

bool VRRuntime::initPhysicsConnection() {
    // Create IPC channel to physics backend
    try {
        physics_channel_ = std::make_unique<UQFF::IPC::SharedMemoryChannel>(
            config_.physics_ipc_channel, 1024 * 1024, false);  // 1MB, open existing
        
        if (!physics_channel_->is_connected()) {
            std::cerr << "VRRuntime: Could not connect to physics backend" << std::endl;
            std::cerr << "  Hint: Start physics backend with: Source2.exe --service" << std::endl;
            physics_channel_.reset();
            return false;
        }
    } catch (const std::exception& e) {
        std::cerr << "VRRuntime: IPC initialization failed: " << e.what() << std::endl;
        std::cerr << "  Hint: Start physics backend with: Source2.exe --service" << std::endl;
        physics_channel_.reset();
        return false;
    }
    
    if (config_.verbose) {
        std::cout << "  Physics: Connected to " << config_.physics_ipc_channel << std::endl;
    }
    
    return true;
}

// ============================================================================
// Main Loop
// ============================================================================

bool VRRuntime::run() {
    if (state_ != RuntimeState::Ready) {
        std::cerr << "VRRuntime: Not in Ready state" << std::endl;
        return false;
    }
    
    // Begin OpenXR session
    if (!openxr_->beginSession()) {
        std::cerr << "VRRuntime: Failed to begin OpenXR session" << std::endl;
        return false;
    }
    
    state_ = RuntimeState::Running;
    
    // Start voice recognition if enabled
    if (task_bot_ && config_.enable_voice_commands) {
        task_bot_->startListening();
    }
    
    // Main render loop
    renderLoop();
    
    // Stop voice recognition
    if (task_bot_) {
        task_bot_->stopListening();
    }
    
    // End OpenXR session
    openxr_->endSession();
    
    return true;
}

void VRRuntime::requestShutdown() {
    shutdown_requested_ = true;
}

void VRRuntime::renderLoop() {
    auto last_time = std::chrono::steady_clock::now();
    uint64_t frame_count = 0;
    
    while (!shutdown_requested_) {
        // Process OpenXR/system events
        processEvents();
        
        if (state_ == RuntimeState::Paused) {
            // Headset off - wait
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            continue;
        }
        
        // Begin frame
        FrameState frame = beginFrame();
        
        if (frame.frame_number == 0) {
            // Frame not ready - skip
            continue;
        }
        
        // Update physics for this frame
        updatePhysics(frame);
        
        // Render for each eye
        int eye_count = openxr_->getViewCount();
        for (int eye = 0; eye < eye_count; ++eye) {
            vulkan_->beginFrame(eye);
            
            // Get view/projection matrices from OpenXR
            ViewInfo view = openxr_->getView(eye);
            float view_matrix[16], proj_matrix[16];
            // TODO: Convert view to matrices
            
            renderScene(frame);
            
            vulkan_->endFrame(eye);
        }
        
        // Submit frame
        endFrame(frame);
        
        // Update metrics
        auto now = std::chrono::steady_clock::now();
        double dt = std::chrono::duration<double, std::milli>(now - last_time).count();
        last_time = now;
        
        metrics_.frame_time_ms = dt;
        metrics_.frames_rendered = ++frame_count;
        metrics_.avg_fps = 1000.0 / dt;
    }
}

void VRRuntime::processEvents() {
    // OpenXR events are processed in waitFrame()
    
    // Process Task Bot commands
    if (task_bot_ && task_bot_->getPendingTaskCount() > 0) {
        task_bot_->executeNextTask();
    }
}

void VRRuntime::updatePhysics(FrameState& frame) {
    if (!physics_channel_ || !physics_channel_->is_connected()) {
        frame.physics_data.valid = false;
        return;
    }
    
    // Request field at head position
    requestFieldUpdate(frame.head_pose.position, 1);
    
    // Get result (non-blocking)
    getFieldResult(frame.physics_data);
}

void VRRuntime::renderScene(const FrameState& frame) {
    // Clear background (space black)
    float bg[4] = {0.0f, 0.0f, 0.02f, 1.0f};
    vulkan_->clearBackground(bg);
    
    // Render field visualization
    vulkan_->renderFieldVisualization();
    
    // Render astronomical objects
    vulkan_->renderAstronomicalObjects();
    
    // Render point cloud if available
    vulkan_->renderPointCloud();
    
    // Render astro graphics program if loaded
    if (astro_graphics_) {
        // Get view/proj from current eye
        float view[16], proj[16];
        // TODO: Get actual matrices
        astro_graphics_->render(0, view, proj);
    }
    
    // Performance overlay
    if (config_.enable_performance_overlay) {
        vulkan_->beginUIOverlay();
        
        std::ostringstream oss;
        oss << "FPS: " << static_cast<int>(metrics_.avg_fps);
        vulkan_->renderText(oss.str(), 10, 10, 16);
        
        oss.str("");
        oss << "Frame: " << metrics_.frame_time_ms << " ms";
        vulkan_->renderText(oss.str(), 10, 30, 16);
        
        if (frame.physics_data.valid) {
            oss.str("");
            oss << "F_U: " << frame.physics_data.F_U;
            vulkan_->renderText(oss.str(), 10, 50, 16);
        }
        
        vulkan_->endUIOverlay();
    }
}

// ============================================================================
// Frame Control
// ============================================================================

FrameState VRRuntime::beginFrame() {
    FrameState frame;
    
    // Wait for OpenXR frame timing
    auto timing = openxr_->waitFrame();
    
    if (!timing.should_render) {
        frame.frame_number = 0;  // Signal skip
        return frame;
    }
    
    openxr_->beginFrame();
    
    frame.frame_number = timing.frame_index;
    frame.timestamp = std::chrono::steady_clock::now();
    frame.predicted_display_time = timing.predicted_display_time / 1e9;  // ns to s
    
    // Get head pose
    openxr_->getHeadPose(frame.head_pose.position, frame.head_pose.orientation);
    
    // Get controller poses
    auto left_pose = openxr_->getControllerPose(true);
    if (left_pose.pose_valid) {
        frame.left_hand.valid = true;
        std::copy(left_pose.position, left_pose.position + 3, frame.left_hand.position);
        std::copy(left_pose.orientation, left_pose.orientation + 4, frame.left_hand.orientation);
    }
    
    auto right_pose = openxr_->getControllerPose(false);
    if (right_pose.pose_valid) {
        frame.right_hand.valid = true;
        std::copy(right_pose.position, right_pose.position + 3, frame.right_hand.position);
        std::copy(right_pose.orientation, right_pose.orientation + 4, frame.right_hand.orientation);
    }
    
    // Process hand tracking for gesture recognition
    if (task_bot_ && openxr_->isHandTrackingSupported()) {
        auto left_joints = openxr_->getHandJoints(true);
        auto right_joints = openxr_->getHandJoints(false);
        task_bot_->processHandTracking(left_joints, right_joints);
    }
    
    return frame;
}

void VRRuntime::endFrame(const FrameState& frame) {
    // Submit frame to OpenXR
    openxr_->endFrame(frame.frame_number > 0);
}

// ============================================================================
// Physics Integration
// ============================================================================

bool VRRuntime::connectPhysicsBackend(const std::string& channel) {
    try {
        // Open existing shared memory channel created by physics backend
        physics_channel_ = std::make_unique<UQFF::IPC::SharedMemoryChannel>(
            channel, 1024 * 1024, false  // 1MB buffer, open existing
        );
        return physics_channel_->is_connected();
    } catch (const std::exception& e) {
        std::cerr << "VRRuntime: Failed to connect to physics backend: " << e.what() << std::endl;
        return false;
    }
}

void VRRuntime::requestFieldUpdate(const float* position, int count) {
    if (!physics_channel_ || !physics_channel_->is_connected()) {
        return;
    }
    
    // Create field request message header
    UQFF::IPC::MessageHeader header(UQFF::IPC::MessageType::CALCULATE_FIELD, 
                                     static_cast<uint32_t>(count * sizeof(float) * 3));
    
    // Send header + position payload
    physics_channel_->send(header, position);
}

bool VRRuntime::getFieldResult(FrameState::PhysicsData& result) {
    if (!physics_channel_ || !physics_channel_->is_connected()) {
        result.valid = false;
        return false;
    }
    
    UQFF::IPC::MessageHeader header;
    if (!physics_channel_->receive(header, physics_payload_buffer_, 0)) {  // Non-blocking
        result.valid = false;
        return false;
    }
    
    if (header.type != UQFF::IPC::MessageType::RESPONSE_DATA) {
        result.valid = false;
        return false;
    }
    
    // Unpack field data from payload
    if (physics_payload_buffer_.size() >= sizeof(double) * 4) {
        const double* data = reinterpret_cast<const double*>(physics_payload_buffer_.data());
        result.F_U = data[0];
        result.field_gradient[0] = data[1];
        result.field_gradient[1] = data[2];
        result.field_gradient[2] = data[3];
        result.valid = true;
        result.timestamp = std::chrono::duration<double>(
            std::chrono::steady_clock::now().time_since_epoch()).count();
        
        // Update metrics
        metrics_.physics_latency_ms = 0;  // TODO: Calculate actual latency
        
        return true;
    }
    
    result.valid = false;
    return false;
}

// ============================================================================
// Task Bot
// ============================================================================

void VRRuntime::registerCommand(const std::string& trigger, std::function<void()> callback) {
    if (task_bot_) {
        task_bot_->registerCommand(trigger, [callback](const VoiceCommand&) {
            callback();
            return TaskResult{true, "OK", "", 0};
        });
    }
}

// ============================================================================
// Astro Graphics
// ============================================================================

bool VRRuntime::loadAstroProgram(const std::string& path) {
    astro_graphics_ = std::make_unique<AstroGraphics>();
    
    if (!astro_graphics_->initialize(vulkan_.get())) {
        std::cerr << "VRRuntime: Failed to initialize Astro Graphics" << std::endl;
        astro_graphics_.reset();
        return false;
    }
    
    // Phase 4: Connect IPC channel to AstroGraphics for physics backend communication
    if (physics_channel_ && physics_channel_->is_connected()) {
        astro_graphics_->setPhysicsChannel(physics_channel_.get());
        if (config_.verbose) {
            std::cout << "  Astro Graphics: IPC channel connected" << std::endl;
        }
    }
    
    if (!astro_graphics_->loadExternalProgram(path)) {
        std::cerr << "VRRuntime: Failed to load program: " << path << std::endl;
        // Continue - can work without external program
    }
    
    if (config_.verbose) {
        std::cout << "  Astro Graphics: Loaded" << std::endl;
    }
    
    return true;
}

} // namespace VR

// ============================================================================
// Main Entry Point (Standalone VR Runtime)
// ============================================================================

#ifdef VR_RUNTIME_STANDALONE

int main(int argc, char* argv[]) {
    std::cout << "Star-Magic UQFF VR Runtime v3.0" << std::endl;
    std::cout << "Phase 5 - Full VR Experience" << std::endl;
    std::cout << "=============================" << std::endl;
    
    // Parse command line
    VR::RuntimeConfig config = VR::RuntimeConfig::from_args(argc, argv);
    
    // Initialize runtime
    VR::VRRuntime& runtime = VR::VRRuntime::instance();
    
    if (!runtime.initialize(config)) {
        std::cerr << "Failed to initialize VR runtime" << std::endl;
        return 1;
    }
    
    // Run main loop
    if (!runtime.run()) {
        std::cerr << "VR runtime exited with error" << std::endl;
        return 1;
    }
    
    // Cleanup
    runtime.shutdown();
    
    std::cout << "VR runtime exited normally" << std::endl;
    return 0;
}

#endif // VR_RUNTIME_STANDALONE
