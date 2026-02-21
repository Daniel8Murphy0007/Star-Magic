/**
 * @file vr_runtime.h
 * @brief VR Runtime Layer for UQFF Star-Magic (Phase 4)
 * 
 * VR Compositor evolved from source2(HEAD PROGRAM).cpp
 * Provides:
 * - OpenXR session management for VR headset interaction
 * - Vulkan/DirectX12 GPU compositing
 * - Task Bot agent for voice/gesture commands
 * - Astro Graphics Engine integration with IPC (Phase 4)
 * - IPC connection to Physics Backend Service
 * 
 * Architecture:
 *   User (VR Headset) → VR Runtime → Physics Backend (source2.cpp --service)
 *                                  → MAIN_1_CoAnQi (physics engine)
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 5 - Full VR Experience
 */

#ifndef VR_RUNTIME_H
#define VR_RUNTIME_H

// Standard Library
#include <string>
#include <vector>
#include <memory>
#include <thread>
#include <atomic>
#include <functional>
#include <chrono>

// IPC support
#include "../ipc/uqff_ipc.h"

namespace VR {

// Forward declarations
class OpenXRSession;
class VulkanCompositor;
class TaskBot;
class AstroGraphics;

/**
 * @enum RuntimeState
 * @brief VR runtime lifecycle states
 */
enum class RuntimeState {
    Uninitialized,      // Not yet initialized
    Initializing,       // Loading dependencies
    Ready,              // Ready to start VR session
    Running,            // Active VR session
    Paused,             // Session paused (headset removed)
    ShuttingDown,       // Cleaning up resources
    Error               // Error state
};

/**
 * @struct RuntimeConfig
 * @brief Configuration for VR runtime startup
 */
struct RuntimeConfig {
    // OpenXR Settings
    std::string application_name = "Star-Magic UQFF VR";
    uint32_t engine_version = 0x00030000;  // v3.0
    bool require_head_tracking = true;
    bool require_hand_tracking = false;    // Optional
    
    // Rendering Settings
    uint32_t render_width = 2048;          // Per-eye width
    uint32_t render_height = 2048;         // Per-eye height
    float refresh_rate = 90.0f;            // Target refresh rate (Hz)
    bool enable_msaa = true;               // Multi-sample anti-aliasing
    int msaa_samples = 4;
    
    // Physics Backend
    std::string physics_ipc_channel = "uqff_physics";
    std::string physics_host = "localhost";
    uint16_t physics_port = 50051;
    bool use_shared_memory = true;         // Prefer SharedMem over gRPC
    
    // Task Bot Settings
    bool enable_voice_commands = true;
    bool enable_gesture_recognition = true;
    std::string voice_model_path = "models/pocketsphinx/";
    
    // Astro Graphics
    std::string astro_program_path = "";   // Path to astronomical graphics program
    bool enable_field_visualization = true;
    double field_sample_rate = 60.0;       // Hz - field updates per second
    
    // Debug
    bool verbose = false;
    bool enable_performance_overlay = false;
    
    // Parse from command line
    static RuntimeConfig from_args(int argc, char* argv[]);
};

/**
 * @struct FrameState
 * @brief Per-frame rendering state
 */
struct FrameState {
    uint64_t frame_number = 0;
    std::chrono::steady_clock::time_point timestamp;
    double predicted_display_time = 0.0;   // In seconds
    
    // Head pose (6DOF)
    struct {
        float position[3] = {0, 0, 0};     // XYZ in meters
        float orientation[4] = {0, 0, 0, 1}; // Quaternion XYZW
    } head_pose;
    
    // Controller poses (optional)
    struct {
        bool valid = false;
        float position[3] = {0, 0, 0};
        float orientation[4] = {0, 0, 0, 1};
        float trigger = 0.0f;
        float grip = 0.0f;
        bool buttons[4] = {false};
    } left_hand, right_hand;
    
    // Physics field data from backend
    struct PhysicsData {
        bool valid = false;
        double F_U = 0;                    // Total unified field at head position
        double field_gradient[3] = {0};    // Field gradient vector
        double timestamp = 0;              // Physics calculation timestamp
    } physics_data;
};

/**
 * @struct PerformanceMetrics
 * @brief VR runtime performance tracking
 */
struct PerformanceMetrics {
    double frame_time_ms = 0;              // Time to render frame
    double physics_latency_ms = 0;         // IPC round-trip time
    double gpu_util_percent = 0;           // GPU utilization
    uint64_t frames_rendered = 0;
    uint64_t frames_dropped = 0;
    double avg_fps = 0;
};

/**
 * @class VRRuntime
 * @brief Main VR runtime controller
 * 
 * Manages the VR experience lifecycle:
 * 1. Initialize OpenXR session
 * 2. Create Vulkan compositor
 * 3. Connect to Physics Backend
 * 4. Start Task Bot agent
 * 5. Load Astro Graphics program
 * 6. Run render loop
 */
class VRRuntime {
public:
    /**
     * @brief Get singleton instance
     */
    static VRRuntime& instance();
    
    // Lifecycle
    bool initialize(const RuntimeConfig& config);
    void shutdown();
    RuntimeState getState() const { return state_; }
    
    // Main loop
    bool run();  // Blocking - runs until shutdown
    void requestShutdown();
    
    // Frame control
    FrameState beginFrame();
    void endFrame(const FrameState& frame);
    
    // Physics integration
    bool connectPhysicsBackend(const std::string& channel);
    void requestFieldUpdate(const float* position, int count);
    bool getFieldResult(FrameState::PhysicsData& result);
    
    // Task Bot
    TaskBot* getTaskBot() { return task_bot_.get(); }
    void registerCommand(const std::string& trigger, std::function<void()> callback);
    
    // Astro Graphics
    bool loadAstroProgram(const std::string& path);
    AstroGraphics* getAstroGraphics() { return astro_graphics_.get(); }
    
    // Metrics
    PerformanceMetrics getMetrics() const { return metrics_; }
    
    // Accessors
    OpenXRSession* getOpenXRSession() { return openxr_.get(); }
    VulkanCompositor* getVulkanCompositor() { return vulkan_.get(); }
    const RuntimeConfig& getConfig() const { return config_; }

private:
    VRRuntime() = default;
    ~VRRuntime();
    VRRuntime(const VRRuntime&) = delete;
    VRRuntime& operator=(const VRRuntime&) = delete;
    
    // Initialization steps
    bool initOpenXR();
    bool initVulkan();
    bool initTaskBot();
    bool initPhysicsConnection();
    
    // Render loop internals
    void renderLoop();
    void processEvents();
    void updatePhysics(FrameState& frame);
    void renderScene(const FrameState& frame);
    void submitFrame(const FrameState& frame);
    
    // State
    RuntimeState state_ = RuntimeState::Uninitialized;
    RuntimeConfig config_;
    std::atomic<bool> shutdown_requested_{false};
    
    // Components
    std::unique_ptr<OpenXRSession> openxr_;
    std::unique_ptr<VulkanCompositor> vulkan_;
    std::unique_ptr<TaskBot> task_bot_;
    std::unique_ptr<AstroGraphics> astro_graphics_;
    
    // Physics IPC (uses IChannel interface from uqff_ipc.h)
    std::unique_ptr<UQFF::IPC::IChannel> physics_channel_;
    std::vector<uint8_t> physics_payload_buffer_;  // Reusable buffer
    
    // Metrics
    PerformanceMetrics metrics_;
    
    // Threading
    std::thread render_thread_;
    std::thread physics_thread_;
};

} // namespace VR

#endif // VR_RUNTIME_H
