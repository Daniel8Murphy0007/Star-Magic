/**
 * @file task_bot.h
 * @brief Task Bot Agent for UQFF Star-Magic VR Runtime
 * 
 * Voice and gesture command processing:
 * - PocketSphinx speech recognition
 * - OpenCV gesture recognition
 * - Natural language command parsing
 * - Physics calculation dispatch
 * - Automation task queuing
 * 
 * Command Categories:
 * - Navigation: "go to", "orbit", "zoom"
 * - Physics: "calculate field", "show gradient", "trace field line"
 * - Visualization: "change color", "toggle wireframe", "show vectors"
 * - System: "take screenshot", "record", "save state"
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 3 - VR Runtime Scaffold
 */

#ifndef TASK_BOT_H
#define TASK_BOT_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>
#include <queue>
#include <mutex>
#include <thread>
#include <atomic>

namespace VR {

// Forward declarations
class VRRuntime;

/**
 * @enum GestureType
 * @brief Recognized hand gestures
 */
enum class GestureType {
    None,
    PointAt,        // Index finger extended
    Pinch,          // Thumb + index together
    Grab,           // Full hand grip
    Open,           // Open palm
    ThumbsUp,       // Thumbs up
    Wave,           // Side-to-side wave
    Swipe,          // Quick directional movement
    CircleGesture   // Circular motion
};

/**
 * @struct GestureEvent
 * @brief Detected gesture event
 */
struct GestureEvent {
    GestureType type = GestureType::None;
    bool is_left_hand = false;
    float confidence = 0.0f;
    
    // For directional gestures
    float direction[3] = {0, 0, 0};     // Unit vector
    float velocity = 0.0f;               // m/s
    
    // For pointing
    float target_point[3] = {0, 0, 0};  // World space hit point
    std::string target_object;           // Object name if hit
};

/**
 * @struct VoiceCommand
 * @brief Parsed voice command
 */
struct VoiceCommand {
    std::string raw_text;               // Original transcription
    std::string intent;                 // Parsed intent (e.g., "navigate", "calculate")
    std::map<std::string, std::string> parameters;  // Extracted parameters
    float confidence = 0.0f;
    int64_t timestamp = 0;              // Unix timestamp
};

/**
 * @struct TaskResult
 * @brief Result of task execution
 */
struct TaskResult {
    bool success = false;
    std::string message;
    std::string data;                   // JSON or other result data
    double execution_time_ms = 0;
};

/**
 * @class TaskBot
 * @brief Voice/gesture command processor and task executor
 */
class TaskBot {
public:
    TaskBot();
    ~TaskBot();
    
    // Initialization
    bool initialize(const std::string& voice_model_path, bool enable_gestures);
    void shutdown();
    
    // Voice recognition
    void startListening();
    void stopListening();
    bool isListening() const { return listening_.load(); }
    void setWakeWord(const std::string& word) { wake_word_ = word; }
    
    // Gesture recognition
    void processHandTracking(const struct HandJoints& left, const struct HandJoints& right);
    GestureEvent getLastGesture() const { return last_gesture_; }
    
    // Command registration
    using CommandHandler = std::function<TaskResult(const VoiceCommand&)>;
    void registerCommand(const std::string& intent, CommandHandler handler);
    void registerGestureAction(GestureType gesture, std::function<void(const GestureEvent&)> action);
    
    // Command parsing
    VoiceCommand parseCommand(const std::string& text);
    
    // Task queue
    void queueTask(const std::string& task_json);
    void executeNextTask();
    int getPendingTaskCount() const { return static_cast<int>(task_queue_.size()); }
    
    // Built-in commands
    void registerBuiltinCommands(VRRuntime* runtime);
    
    // Feedback
    void speak(const std::string& text);     // Text-to-speech
    void playSound(const std::string& name); // Audio feedback
    
    // Status
    struct Status {
        bool voice_active = false;
        bool gesture_active = false;
        std::string last_command;
        std::string last_result;
        int commands_processed = 0;
    };
    Status getStatus() const { return status_; }
    
    // Callbacks
    using CommandCallback = std::function<void(const VoiceCommand&, const TaskResult&)>;
    void setCommandCallback(CommandCallback cb) { command_callback_ = cb; }

private:
    // Voice recognition thread
    std::thread voice_thread_;
    std::atomic<bool> listening_{false};
    std::atomic<bool> shutdown_requested_{false};
    
    // PocketSphinx (optional)
#ifdef HAVE_POCKETSPHINX
    void* ps_decoder_ = nullptr;
#endif
    std::string voice_model_path_;
    std::string wake_word_ = "computer";
    
    // Gesture state
    GestureEvent last_gesture_;
    bool gesture_enabled_ = false;
    
    // Command handlers
    std::map<std::string, CommandHandler> command_handlers_;
    std::map<GestureType, std::function<void(const GestureEvent&)>> gesture_handlers_;
    
    // Task queue (thread-safe)
    std::queue<std::string> task_queue_;
    mutable std::mutex queue_mutex_;
    
    // Status
    Status status_;
    
    // Callback
    CommandCallback command_callback_;
    
    // Internal methods
    void voiceRecognitionLoop();
    GestureType classifyGesture(const struct HandJoints& hand);
    std::vector<std::string> tokenize(const std::string& text);
    std::string extractParameter(const std::vector<std::string>& tokens, const std::string& key);
};

// ============================================================================
// Built-in command intents
// ============================================================================

namespace Intent {
    // Navigation
    constexpr const char* NAVIGATE = "navigate";      // "go to [object]"
    constexpr const char* ORBIT = "orbit";            // "orbit [object]"
    constexpr const char* ZOOM = "zoom";              // "zoom in/out"
    constexpr const char* RESET_VIEW = "reset_view";  // "reset view"
    
    // Physics
    constexpr const char* CALCULATE = "calculate";    // "calculate field at [position]"
    constexpr const char* TRACE = "trace";            // "trace field line"
    constexpr const char* SHOW_GRADIENT = "gradient"; // "show gradient"
    constexpr const char* COMPARE = "compare";        // "compare [system1] and [system2]"
    
    // Visualization
    constexpr const char* TOGGLE_MODE = "toggle";     // "toggle wireframe"
    constexpr const char* SET_COLOR = "color";        // "change color to [color]"
    constexpr const char* SCALE = "scale";            // "scale up/down"
    constexpr const char* HIGHLIGHT = "highlight";    // "highlight [object]"
    
    // System
    constexpr const char* SCREENSHOT = "screenshot";  // "take screenshot"
    constexpr const char* RECORD = "record";          // "start/stop recording"
    constexpr const char* SAVE = "save";              // "save state"
    constexpr const char* LOAD = "load";              // "load state"
    constexpr const char* HELP = "help";              // "help"/"what can I do"
    constexpr const char* EXIT = "exit";              // "exit VR"
}

} // namespace VR

#endif // TASK_BOT_H
