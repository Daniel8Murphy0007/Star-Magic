/**
 * @file task_bot.cpp
 * @brief Task Bot Implementation for UQFF Star-Magic VR Runtime
 * 
 * Voice and gesture command processing
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 3 - VR Runtime Scaffold
 */

#include "task_bot.h"
#include "openxr_session.h"
#include "vulkan_compositor.h"
#include "astro_graphics.h"
#include "vr_runtime.h"  // Include last to get full type definitions above

#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cmath>

namespace VR {

// ============================================================================
// Constructor / Destructor
// ============================================================================

TaskBot::TaskBot() = default;

TaskBot::~TaskBot() {
    shutdown();
}

// ============================================================================
// Initialization
// ============================================================================

bool TaskBot::initialize(const std::string& voice_model_path, bool enable_gestures) {
    voice_model_path_ = voice_model_path;
    gesture_enabled_ = enable_gestures;
    
#ifdef HAVE_POCKETSPHINX
    // Initialize PocketSphinx
    ps_config_t* config = ps_config_init(nullptr);
    ps_config_set_str(config, "hmm", (voice_model_path + "/en-us").c_str());
    ps_config_set_str(config, "lm", (voice_model_path + "/en-us.lm.bin").c_str());
    ps_config_set_str(config, "dict", (voice_model_path + "/cmudict-en-us.dict").c_str());
    
    ps_decoder_ = ps_init(config);
    if (!ps_decoder_) {
        std::cerr << "TaskBot: Failed to initialize PocketSphinx" << std::endl;
        return false;
    }
    
    std::cout << "TaskBot: Voice recognition initialized" << std::endl;
#else
    std::cout << "TaskBot: PocketSphinx not available - voice commands disabled" << std::endl;
    std::cout << "  To enable: vcpkg install pocketsphinx:x64-windows" << std::endl;
#endif
    
    if (gesture_enabled_) {
        std::cout << "TaskBot: Gesture recognition enabled" << std::endl;
    }
    
    return true;
}

void TaskBot::shutdown() {
    shutdown_requested_ = true;
    
    if (voice_thread_.joinable()) {
        voice_thread_.join();
    }
    
#ifdef HAVE_POCKETSPHINX
    if (ps_decoder_) {
        ps_free(static_cast<ps_decoder_t*>(ps_decoder_));
        ps_decoder_ = nullptr;
    }
#endif
}

// ============================================================================
// Voice Recognition
// ============================================================================

void TaskBot::startListening() {
    if (listening_.load()) {
        return;  // Already listening
    }
    
#ifdef HAVE_POCKETSPHINX
    listening_ = true;
    shutdown_requested_ = false;
    
    voice_thread_ = std::thread(&TaskBot::voiceRecognitionLoop, this);
    status_.voice_active = true;
#else
    std::cout << "TaskBot: Voice recognition not available" << std::endl;
#endif
}

void TaskBot::stopListening() {
    listening_ = false;
    status_.voice_active = false;
}

void TaskBot::voiceRecognitionLoop() {
#ifdef HAVE_POCKETSPHINX
    // Audio capture and speech recognition loop
    while (!shutdown_requested_.load() && listening_.load()) {
        // Capture audio
        // Process with PocketSphinx
        // Parse commands
        // Execute handlers
        
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
#endif
}

// ============================================================================
// Gesture Recognition
// ============================================================================

void TaskBot::processHandTracking(const HandJoints& left, const HandJoints& right) {
    if (!gesture_enabled_) {
        return;
    }
    
    // Process left hand
    if (left.valid) {
        GestureType gesture = classifyGesture(left);
        if (gesture != GestureType::None) {
            GestureEvent event;
            event.type = gesture;
            event.is_left_hand = true;
            event.confidence = 0.9f;  // Placeholder
            
            last_gesture_ = event;
            
            // Execute gesture handler
            auto it = gesture_handlers_.find(gesture);
            if (it != gesture_handlers_.end()) {
                it->second(event);
            }
        }
    }
    
    // Process right hand
    if (right.valid) {
        GestureType gesture = classifyGesture(right);
        if (gesture != GestureType::None) {
            GestureEvent event;
            event.type = gesture;
            event.is_left_hand = false;
            event.confidence = 0.9f;
            
            last_gesture_ = event;
            
            auto it = gesture_handlers_.find(gesture);
            if (it != gesture_handlers_.end()) {
                it->second(event);
            }
        }
    }
    
    status_.gesture_active = true;
}

GestureType TaskBot::classifyGesture(const HandJoints& hand) {
    // Simple gesture classification based on joint positions
    // In production, would use ML model or more sophisticated tracking
    
    if (!hand.valid) {
        return GestureType::None;
    }
    
    // Check for pinch (thumb tip close to index tip)
    // Joint indices: 4 = thumb tip, 8 = index tip
    const auto& thumb_tip = hand.joints[4];
    const auto& index_tip = hand.joints[8];
    
    float dx = thumb_tip.position[0] - index_tip.position[0];
    float dy = thumb_tip.position[1] - index_tip.position[1];
    float dz = thumb_tip.position[2] - index_tip.position[2];
    float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    if (dist < 0.02f) {  // 2 cm threshold
        return GestureType::Pinch;
    }
    
    // Check for pointing (index extended, others curled)
    // This is simplified - production would check all fingers
    const auto& index_mcp = hand.joints[5];
    float index_extension = std::sqrt(
        std::pow(index_tip.position[0] - index_mcp.position[0], 2) +
        std::pow(index_tip.position[1] - index_mcp.position[1], 2) +
        std::pow(index_tip.position[2] - index_mcp.position[2], 2)
    );
    
    if (index_extension > 0.08f) {  // Index extended
        return GestureType::PointAt;
    }
    
    // Check for open palm
    // (would check all fingers extended)
    
    // Check for fist
    // (would check all fingers curled)
    
    return GestureType::None;
}

// ============================================================================
// Command Registration
// ============================================================================

void TaskBot::registerCommand(const std::string& intent, CommandHandler handler) {
    command_handlers_[intent] = std::move(handler);
}

void TaskBot::registerGestureAction(GestureType gesture, std::function<void(const GestureEvent&)> action) {
    gesture_handlers_[gesture] = std::move(action);
}

// ============================================================================
// Command Parsing
// ============================================================================

VoiceCommand TaskBot::parseCommand(const std::string& text) {
    VoiceCommand cmd;
    cmd.raw_text = text;
    cmd.timestamp = std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    
    // Tokenize
    auto tokens = tokenize(text);
    if (tokens.empty()) {
        return cmd;
    }
    
    // Simple intent classification
    std::string first = tokens[0];
    std::transform(first.begin(), first.end(), first.begin(), ::tolower);
    
    // Navigation intents
    if (first == "go" || first == "navigate" || first == "fly") {
        cmd.intent = Intent::NAVIGATE;
        cmd.parameters["target"] = extractParameter(tokens, "to");
    }
    else if (first == "orbit") {
        cmd.intent = Intent::ORBIT;
        if (tokens.size() > 1) {
            cmd.parameters["target"] = tokens[1];
        }
    }
    else if (first == "zoom") {
        cmd.intent = Intent::ZOOM;
        if (tokens.size() > 1) {
            cmd.parameters["direction"] = tokens[1];  // "in" or "out"
        }
    }
    
    // Physics intents
    else if (first == "calculate") {
        cmd.intent = Intent::CALCULATE;
        cmd.parameters["field"] = extractParameter(tokens, "field");
        cmd.parameters["at"] = extractParameter(tokens, "at");
    }
    else if (first == "trace") {
        cmd.intent = Intent::TRACE;
    }
    else if (first == "show" || first == "display") {
        if (tokens.size() > 1 && tokens[1] == "gradient") {
            cmd.intent = Intent::SHOW_GRADIENT;
        } else {
            cmd.intent = Intent::HIGHLIGHT;
            cmd.parameters["target"] = extractParameter(tokens, "");
        }
    }
    
    // Visualization intents
    else if (first == "toggle") {
        cmd.intent = Intent::TOGGLE_MODE;
        if (tokens.size() > 1) {
            cmd.parameters["mode"] = tokens[1];
        }
    }
    else if (first == "change" && tokens.size() > 1 && tokens[1] == "color") {
        cmd.intent = Intent::SET_COLOR;
        cmd.parameters["color"] = extractParameter(tokens, "to");
    }
    else if (first == "scale") {
        cmd.intent = Intent::SCALE;
        if (tokens.size() > 1) {
            cmd.parameters["direction"] = tokens[1];
        }
    }
    
    // System intents
    else if (first == "screenshot" || (first == "take" && tokens.size() > 1 && tokens[1] == "screenshot")) {
        cmd.intent = Intent::SCREENSHOT;
    }
    else if (first == "record") {
        cmd.intent = Intent::RECORD;
    }
    else if (first == "save") {
        cmd.intent = Intent::SAVE;
    }
    else if (first == "load") {
        cmd.intent = Intent::LOAD;
        if (tokens.size() > 1) {
            cmd.parameters["file"] = tokens[1];
        }
    }
    else if (first == "help" || first == "what") {
        cmd.intent = Intent::HELP;
    }
    else if (first == "exit" || first == "quit") {
        cmd.intent = Intent::EXIT;
    }
    
    cmd.confidence = cmd.intent.empty() ? 0.0f : 0.8f;
    status_.last_command = text;
    status_.commands_processed++;
    
    return cmd;
}

std::vector<std::string> TaskBot::tokenize(const std::string& text) {
    std::vector<std::string> tokens;
    std::istringstream iss(text);
    std::string token;
    
    while (iss >> token) {
        // Remove punctuation
        token.erase(std::remove_if(token.begin(), token.end(), ::ispunct), token.end());
        if (!token.empty()) {
            tokens.push_back(token);
        }
    }
    
    return tokens;
}

std::string TaskBot::extractParameter(const std::vector<std::string>& tokens, const std::string& key) {
    for (size_t i = 0; i < tokens.size(); ++i) {
        if (tokens[i] == key && i + 1 < tokens.size()) {
            return tokens[i + 1];
        }
    }
    
    // Return last token if no key match
    if (!key.empty() && tokens.size() > 1) {
        return tokens.back();
    }
    
    return "";
}

// ============================================================================
// Task Queue
// ============================================================================

void TaskBot::queueTask(const std::string& task_json) {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    task_queue_.push(task_json);
}

void TaskBot::executeNextTask() {
    std::string task_json;
    {
        std::lock_guard<std::mutex> lock(queue_mutex_);
        if (task_queue_.empty()) {
            return;
        }
        task_json = task_queue_.front();
        task_queue_.pop();
    }
    
    // Parse JSON task and execute
    // For now, just parse as voice command
    VoiceCommand cmd = parseCommand(task_json);
    
    auto it = command_handlers_.find(cmd.intent);
    if (it != command_handlers_.end()) {
        TaskResult result = it->second(cmd);
        status_.last_result = result.message;
        
        if (command_callback_) {
            command_callback_(cmd, result);
        }
    }
}

// ============================================================================
// Built-in Commands
// ============================================================================

void TaskBot::registerBuiltinCommands(VRRuntime* runtime) {
    if (!runtime) {
        return;
    }
    
    // Exit command
    registerCommand(Intent::EXIT, [runtime](const VoiceCommand&) {
        runtime->requestShutdown();
        return TaskResult{true, "Exiting VR", "", 0};
    });
    
    // Help command
    registerCommand(Intent::HELP, [](const VoiceCommand&) {
        std::string help = R"(
Available commands:
- Navigation: "go to [object]", "orbit [object]", "zoom in/out"
- Physics: "calculate field", "trace field line", "show gradient"
- Visualization: "toggle wireframe", "change color", "scale up/down"
- System: "screenshot", "record", "save state", "exit"
)";
        return TaskResult{true, help, "", 0};
    });
    
    // Navigate command
    registerCommand(Intent::NAVIGATE, [runtime](const VoiceCommand& cmd) {
        auto it = cmd.parameters.find("target");
        if (it != cmd.parameters.end()) {
            // TODO: Implement navigation when AstroGraphics is fully integrated
            // runtime->getAstroGraphics()->flyTo(it->second, 2.0);
            return TaskResult{true, "Flying to " + it->second, "", 0};
        }
        return TaskResult{false, "No target specified", "", 0};
    });
    
    // Toggle wireframe
    registerCommand(Intent::TOGGLE_MODE, [runtime](const VoiceCommand& cmd) {
        auto it = cmd.parameters.find("mode");
        (void)runtime;  // Suppress unused warning - will use when VulkanCompositor is integrated
        if (it != cmd.parameters.end()) {
            if (it->second == "wireframe") {
                // TODO: Implement render mode toggle when VulkanCompositor is fully integrated
                // VulkanCompositor* compositor = runtime->getVulkanCompositor();
                // VR::RenderMode current = compositor->getRenderMode();
                // VR::RenderMode newMode = (current == VR::RenderMode::Wireframe) ? 
                //                VR::RenderMode::Standard : VR::RenderMode::Wireframe;
                // compositor->setRenderMode(newMode);
                return TaskResult{true, "Toggled wireframe", "", 0};
            }
        }
        return TaskResult{false, "Unknown mode", "", 0};
    });
    
    // Screenshot (placeholder)
    registerCommand(Intent::SCREENSHOT, [](const VoiceCommand&) {
        return TaskResult{true, "Screenshot saved", "screenshot.png", 0};
    });
    
    // Register gesture handlers
    registerGestureAction(GestureType::Pinch, [](const GestureEvent& event) {
        std::cout << "Pinch detected on " << (event.is_left_hand ? "left" : "right") << " hand" << std::endl;
    });
    
    registerGestureAction(GestureType::PointAt, [runtime](const GestureEvent& event) {
        // Ray cast from hand to find pointed object
        (void)runtime;  // To be used when AstroGraphics is fully integrated
        (void)event;
        // TODO: Implement ray picking when VR is fully integrated
        std::cout << "Pointing detected" << std::endl;
    });
}

// ============================================================================
// Feedback
// ============================================================================

void TaskBot::speak(const std::string& text) {
    // TODO: Implement text-to-speech
    std::cout << "TaskBot says: " << text << std::endl;
}

void TaskBot::playSound(const std::string& name) {
    // TODO: Implement audio playback
    std::cout << "TaskBot plays sound: " << name << std::endl;
}

} // namespace VR
