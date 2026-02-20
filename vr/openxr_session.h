/**
 * @file openxr_session.h
 * @brief OpenXR Session Management for UQFF Star-Magic VR Runtime
 * 
 * Handles OpenXR runtime initialization, session lifecycle, and input:
 * - XR instance creation with extensions
 * - System enumeration and selection
 * - Session creation and state management
 * - Reference space setup (Local, Stage, Unbounded)
 * - Input action binding (controllers, hand tracking)
 * - Frame timing and prediction
 * 
 * OpenXR Extensions Used:
 * - XR_KHR_vulkan_enable2 (Vulkan graphics)
 * - XR_EXT_hand_tracking (optional)
 * - XR_FB_passthrough (optional, Quest)
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 3 - VR Runtime Scaffold
 */

#ifndef OPENXR_SESSION_H
#define OPENXR_SESSION_H

// Conditional compilation for OpenXR availability
#ifdef HAVE_OPENXR
#include <openxr/openxr.h>
#include <openxr/openxr_platform.h>
#endif

#include <string>
#include <vector>
#include <memory>
#include <functional>

namespace VR {

/**
 * @enum SessionState
 * @brief OpenXR session lifecycle states
 */
enum class SessionState {
    Unknown,
    Idle,           // Session created but not running
    Ready,          // Ready to begin session
    Synchronized,   // Frames synced with compositor
    Visible,        // App visible but not focused
    Focused,        // App has full control
    Stopping,       // Session ending
    LossPending,    // Runtime is shutting down
    Exiting         // Session exited
};

/**
 * @struct SwapchainInfo
 * @brief Swapchain configuration for rendering
 */
struct SwapchainInfo {
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t sample_count = 1;
    int64_t format = 0;         // Vulkan format
    uint32_t image_count = 0;
    
#ifdef HAVE_OPENXR
    XrSwapchain handle = XR_NULL_HANDLE;
    std::vector<XrSwapchainImageVulkanKHR> images;
#else
    void* handle = nullptr;
    std::vector<void*> images;
#endif
};

/**
 * @struct ViewInfo
 * @brief Per-eye view configuration
 */
struct ViewInfo {
    // Pose
    float position[3] = {0, 0, 0};
    float orientation[4] = {0, 0, 0, 1};
    
    // Field of View
    float fov_left = -0.8f;     // Radians
    float fov_right = 0.8f;
    float fov_up = 0.8f;
    float fov_down = -0.8f;
    
    // Projection
    float near_z = 0.05f;       // 5 cm
    float far_z = 1000.0f;      // 1 km
};

/**
 * @struct ActionState
 * @brief State of an input action
 */
struct ActionState {
    bool active = false;
    
    // For boolean actions
    bool current_state = false;
    bool changed_since_last_sync = false;
    
    // For float actions
    float current_float = 0.0f;
    
    // For pose actions
    bool pose_valid = false;
    float position[3] = {0, 0, 0};
    float orientation[4] = {0, 0, 0, 1};
    float velocity[3] = {0, 0, 0};
    float angular_velocity[3] = {0, 0, 0};
};

/**
 * @struct HandJoints
 * @brief Hand tracking joint data (26 joints per hand)
 */
struct HandJoints {
    static constexpr int JOINT_COUNT = 26;  // XR_HAND_JOINT_COUNT_EXT
    
    bool valid = false;
    struct Joint {
        float position[3] = {0, 0, 0};
        float orientation[4] = {0, 0, 0, 1};
        float radius = 0.01f;  // Joint radius in meters
    } joints[JOINT_COUNT];
    
    // Gesture recognition results
    bool pinching = false;     // Thumb + Index close
    bool pointing = false;     // Index extended
    bool fist = false;         // All fingers closed
    float pinch_strength = 0;  // 0.0 = open, 1.0 = closed
};

/**
 * @class OpenXRSession
 * @brief Manages OpenXR instance, session, and input
 */
class OpenXRSession {
public:
    OpenXRSession();
    ~OpenXRSession();
    
    // Initialization
    bool initialize(const std::string& app_name, uint32_t engine_version);
    void shutdown();
    
    // Session lifecycle
    bool beginSession();
    void endSession();
    SessionState getState() const { return session_state_; }
    
    // Frame timing
    struct FrameTiming {
        bool should_render = false;
        double predicted_display_time = 0;  // Nanoseconds
        double predicted_display_period = 0;
        uint64_t frame_index = 0;
    };
    
    FrameTiming waitFrame();
    void beginFrame();
    void endFrame(bool did_render);
    
    // Views (per-eye rendering)
    int getViewCount() const { return static_cast<int>(views_.size()); }
    ViewInfo getView(int index) const;
    SwapchainInfo& getSwapchain(int eye_index);
    int acquireSwapchainImage(int eye_index);
    void releaseSwapchainImage(int eye_index);
    
    // Input
    void syncInput();
    ActionState getControllerPose(bool left_hand);
    ActionState getTrigger(bool left_hand);
    ActionState getGrip(bool left_hand);
    ActionState getButton(bool left_hand, int button_index);
    
    // Hand tracking (optional)
    bool isHandTrackingSupported() const { return hand_tracking_supported_; }
    HandJoints getHandJoints(bool left_hand);
    
    // Reference spaces
    bool setReferenceSpace(const std::string& type);  // "local", "stage", "unbounded"
    void getHeadPose(float* position, float* orientation);
    
    // Extensions
    bool hasExtension(const std::string& name) const;
    void* getGraphicsBinding() const { return graphics_binding_; }
    
    // Vulkan integration
    void* getVulkanInstance() const;
    void* getVulkanPhysicalDevice() const;
    void setVulkanDevice(void* device, uint32_t queue_family);
    
    // Callbacks
    using SessionStateCallback = std::function<void(SessionState old_state, SessionState new_state)>;
    void setSessionStateCallback(SessionStateCallback callback) { state_callback_ = callback; }

private:
    // OpenXR handles
#ifdef HAVE_OPENXR
    XrInstance instance_ = XR_NULL_HANDLE;
    XrSystemId system_id_ = XR_NULL_SYSTEM_ID;
    XrSession session_ = XR_NULL_HANDLE;
    XrSpace reference_space_ = XR_NULL_HANDLE;
    XrSpace view_space_ = XR_NULL_HANDLE;
    
    // Action system
    XrActionSet action_set_ = XR_NULL_HANDLE;
    XrAction pose_action_ = XR_NULL_HANDLE;
    XrAction trigger_action_ = XR_NULL_HANDLE;
    XrAction grip_action_ = XR_NULL_HANDLE;
    XrAction button_action_ = XR_NULL_HANDLE;
    XrSpace hand_spaces_[2] = {XR_NULL_HANDLE, XR_NULL_HANDLE};
    
    // Frame state
    XrFrameState frame_state_ = {};
    
    // Hand tracking
    XrHandTrackerEXT hand_trackers_[2] = {XR_NULL_HANDLE, XR_NULL_HANDLE};
#endif
    
    // State
    SessionState session_state_ = SessionState::Unknown;
    bool hand_tracking_supported_ = false;
    
    // Graphics binding (Vulkan)
    void* graphics_binding_ = nullptr;
    void* vulkan_instance_ = nullptr;
    void* vulkan_physical_device_ = nullptr;
    void* vulkan_device_ = nullptr;
    uint32_t vulkan_queue_family_ = 0;
    
    // Views and swapchains
    std::vector<ViewInfo> views_;
    std::vector<SwapchainInfo> swapchains_;
    
    // Available extensions
    std::vector<std::string> enabled_extensions_;
    
    // Callback
    SessionStateCallback state_callback_;
    
    // Helper methods
    bool createInstance(const std::string& app_name, uint32_t version);
    bool getSystem();
    bool createSession();
    bool createSwapchains();
    bool createActionSet();
    bool createReferenceSpaces();
    void processEvents();
    void onSessionStateChanged(SessionState new_state);
};

} // namespace VR

#endif // OPENXR_SESSION_H
