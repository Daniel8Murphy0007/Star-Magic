/**
 * @file uqff_ipc.h
 * @brief Inter-Process Communication Foundation for UQFF Star-Magic
 * 
 * Provides message types and channel abstractions for communication between:
 * - source2.cpp (GUI Orchestrator)
 * - source2(HEAD PROGRAM).cpp (VR/VM Backend)
 * - MAIN_1_CoAnQi.exe (Physics Engine)
 * - Python calculators (CondensedPhysics.py, QCalc.py)
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 3 - Full gRPC Implementation
 */

#ifndef UQFF_IPC_H
#define UQFF_IPC_H

#include <cstdint>
#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <mutex>
#include <atomic>
#include <chrono>
#include <thread>
#include <queue>
#include <condition_variable>
#include <unordered_map>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX  // Prevent min/max macro conflicts with std::min/std::max
#endif
#include <windows.h>
#else
#include <sys/mman.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <fcntl.h>
#include <unistd.h>
#include <poll.h>
#endif

// gRPC support (Phase 3)
#ifdef USE_GRPC
#include <grpcpp/grpcpp.h>
#include "uqff_service.grpc.pb.h"
#endif

namespace UQFF {
namespace IPC {

// ============================================================================
// MESSAGE TYPES
// ============================================================================

/**
 * @enum MessageType
 * @brief IPC message types for UQFF physics operations
 */
enum class MessageType : uint32_t {
    // Field calculations
    CALCULATE_FIELD         = 0x0001,  // Request F_U calculation
    CALCULATE_GRAVITY       = 0x0002,  // Request Ug1-Ug4 calculation
    CALCULATE_BUOYANCY      = 0x0003,  // Request Ub calculation
    CALCULATE_MAGNETISM     = 0x0004,  // Request Um calculation
    
    // Dynamic term management
    REGISTER_TERM           = 0x0010,  // Register new PhysicsTerm
    UNREGISTER_TERM         = 0x0011,  // Remove a PhysicsTerm
    UPDATE_PARAMETER        = 0x0012,  // Update runtime parameter
    GET_PARAMETER           = 0x0013,  // Query parameter value
    
    // System management
    GET_SYSTEM_LIST         = 0x0020,  // List available systems
    SELECT_SYSTEM           = 0x0021,  // Set active system
    ADD_CUSTOM_SYSTEM       = 0x0022,  // Register custom system
    
    // State management
    EXPORT_STATE            = 0x0030,  // Export module state
    IMPORT_STATE            = 0x0031,  // Import module state
    SYNC_STATE              = 0x0032,  // Synchronize state across processes
    
    // VR-specific
    VR_FRAME_UPDATE         = 0x0100,  // Per-frame VR data update
    VR_GESTURE_INPUT        = 0x0101,  // Gesture command from VR
    VR_RENDER_REQUEST       = 0x0102,  // Request visualization data
    
    // Simulation control (v3.1 - Self-Simulate)
    SIM_START               = 0x0200,  // Start time evolution simulation
    SIM_STOP                = 0x0201,  // Stop active simulation
    SIM_PAUSE               = 0x0202,  // Pause simulation
    SIM_RESUME              = 0x0203,  // Resume paused simulation
    SIM_FRAME               = 0x0210,  // Simulation frame output (streamed)
    SIM_PROGRESS            = 0x0211,  // Simulation progress update
    SIM_COMPLETE            = 0x0212,  // Simulation completed
    
    // Pipeline operations (Feb 24, 2026 - CondensedPhysics integration)
    PIPELINE_PROCESS        = 0x0300,  // Process object through UQFFPipeline
    PIPELINE_RESULT         = 0x0301,  // Pipeline computation result
    PIPELINE_STORE          = 0x0302,  // Store result to OutputDataStore
    PIPELINE_EXPORT         = 0x0303,  // Export result as JSON for IPC
    PIPELINE_CALLBACK       = 0x0310,  // Real-time callback event
    
    // Epoch Framework (March 4, 2026 - Grok Thread 4e0ecf23)
    EPOCH_GET_CURRENT       = 0x0400,  // Get current cosmic epoch (1-5)
    EPOCH_SET               = 0x0401,  // Set epoch for calculations
    EPOCH_CALCULATE_F_U     = 0x0402,  // Calculate F_U at specific epoch
    EPOCH_GET_UG_ACTIVE     = 0x0403,  // Query which Ug ranges active at epoch
    EPOCH_VALIDATION_DATA   = 0x0410,  // Request epoch validation dataset
    
    // Responses
    RESPONSE_SUCCESS        = 0x1000,  // Operation completed
    RESPONSE_ERROR          = 0x1001,  // Operation failed
    RESPONSE_DATA           = 0x1002,  // Data payload follows
    
    // Control
    PING                    = 0xFF00,  // Health check
    PONG                    = 0xFF01,  // Health response
    SHUTDOWN                = 0xFFFF   // Graceful shutdown
};

/**
 * @struct MessageHeader
 * @brief Fixed-size header for all IPC messages (32 bytes)
 */
struct alignas(8) MessageHeader {
    uint32_t magic;          // 0x55514646 = "UQFF"
    uint32_t version;        // Protocol version (1)
    MessageType type;        // Message type
    uint32_t payload_size;   // Size of payload in bytes
    uint64_t timestamp;      // Microseconds since epoch
    uint32_t sequence;       // Sequence number for ordering
    uint32_t flags;          // Reserved flags
    
    static constexpr uint32_t MAGIC = 0x55514646;
    static constexpr uint32_t VERSION = 1;
    
    MessageHeader() 
        : magic(MAGIC), version(VERSION), type(MessageType::PING),
          payload_size(0), timestamp(0), sequence(0), flags(0) {}
    
    MessageHeader(MessageType t, uint32_t size = 0)
        : magic(MAGIC), version(VERSION), type(t), payload_size(size),
          timestamp(std::chrono::duration_cast<std::chrono::microseconds>(
              std::chrono::system_clock::now().time_since_epoch()).count()),
          sequence(0), flags(0) {}
    
    bool is_valid() const {
        return magic == MAGIC && version == VERSION;
    }
};

static_assert(sizeof(MessageHeader) == 32, "MessageHeader must be 32 bytes");

// ============================================================================
// PAYLOAD STRUCTURES
// ============================================================================

/**
 * @struct CalculateFieldRequest
 * @brief Payload for CALCULATE_FIELD message
 */
struct CalculateFieldRequest {
    char system_name[64];    // System identifier
    double r;                // Radial distance (m)
    double t;                // Time (s)
    double tn;               // Negative time factor
    double theta;            // Angular position (rad)
    uint32_t flags;          // Calculation flags
    uint32_t reserved;
};

/**
 * @struct CalculateFieldResponse
 * @brief Response payload for field calculations
 */
struct CalculateFieldResponse {
    double F_U;              // Unified field strength
    double Ug1, Ug2, Ug3, Ug4;   // Gravity components
    double Um;               // Magnetism
    double Ub;               // Buoyancy
    double compressed_g;     // MUGE compressed gravity
    uint32_t status;         // 0 = success
    uint32_t reserved;
};

/**
 * @struct RegisterTermRequest
 * @brief Payload for REGISTER_TERM message
 */
struct RegisterTermRequest {
    char term_name[64];      // Term identifier
    char description[256];   // Human-readable description
    double initial_value;    // Initial contribution
    uint32_t term_type;      // Type enum
    uint32_t flags;
};

/**
 * @struct UpdateParameterRequest
 * @brief Payload for UPDATE_PARAMETER message
 */
struct UpdateParameterRequest {
    char param_name[64];     // Parameter name
    double value;            // New value
    uint32_t module_id;      // Target module (0 = global)
    uint32_t flags;
};

/**
 * @struct VRFrameUpdate
 * @brief Per-frame VR state for real-time visualization
 */
struct VRFrameUpdate {
    uint64_t frame_number;
    double head_position[3];     // x, y, z
    double head_orientation[4];  // quaternion
    double left_hand[7];         // pos(3) + quat(4)
    double right_hand[7];        // pos(3) + quat(4)
    double field_probe_position[3];  // Where to sample field
    uint32_t gesture_flags;
    uint32_t reserved;
};

/**
 * @struct PipelineProcessRequest
 * @brief Payload for PIPELINE_PROCESS message (Feb 24, 2026)
 */
struct PipelineProcessRequest {
    char object_name[128];       // Object identifier (e.g., "SGR 1935+2154")
    uint32_t flags;              // Processing flags
    uint32_t timeout_ms;         // Timeout in milliseconds (0 = default)
    char callback_id[64];        // Optional callback identifier
    uint32_t reserved[4];
};

/**
 * @struct PipelineResultPayload
 * @brief Response payload for PIPELINE_RESULT message (Feb 24, 2026)
 */
struct PipelineResultPayload {
    char query_id[64];           // Pipeline query identifier
    char object_type[32];        // Classification (magnetar, pulsar, etc.)
    uint32_t calculators_run;    // Number of calculators executed
    uint32_t calculation_success; // Number successful
    double compute_time_ms;      // Total computation time
    
    // Core UQFF field results
    double F_U;                  // Total unified field
    double Ug1, Ug2, Ug3, Ug4;  // Gravity components
    double Um;                   // Magnetism
    double Ubi;                  // Buoyancy opposition
    double g_compressed;         // MUGE Compressed gravity
    double g_resonant;           // MUGE Resonant gravity
    
    uint32_t status;             // 0 = success, non-zero = error code
    uint32_t json_payload_follows; // 1 if full JSON follows this header
    
    // Variable-length JSON payload follows if json_payload_follows == 1
};

/**
 * @struct PipelineCallbackPayload
 * @brief Event payload for PIPELINE_CALLBACK message (Feb 24, 2026)
 */
struct PipelineCallbackPayload {
    char event_type[32];         // fetch_start, calc_complete, process_complete
    char query_id[64];           // Associated query
    char object_name[128];       // Object being processed
    char calculator_name[64];    // Calculator involved (for calc_complete)
    double elapsed_ms;           // Time since process start
    uint32_t status;             // 0 = success
    uint32_t reserved[2];
};

// ============================================================================
// EPOCH FRAMEWORK PAYLOADS (March 4, 2026 - Grok Thread 4e0ecf23)
// ============================================================================

/**
 * @struct EpochGetCurrentRequest
 * @brief Query for current cosmic epoch (MESSAGE_TYPE: EPOCH_GET_CURRENT)
 */
struct EpochGetCurrentRequest {
    char system_name[64];        // System to query (use "default" for global)
    double cosmic_time;          // Cosmic time scale (1.0-5.9)
    uint32_t flags;              // Reserved
    uint32_t reserved;
};

/**
 * @struct EpochGetCurrentResponse
 * @brief Response with current epoch information
 */
struct EpochGetCurrentResponse {
    int epoch_number;            // Current epoch (1-5)
    char epoch_name[64];         // "Fisile Nuclei", "Star/Planetary", etc.
    char scm_state[16];          // SCm, SCm'', SCm''', SCm'''', SCm'''''
    bool ug1_active;             // Internal Dipole active?
    bool ug2_active;             // Heliosphere active?
    bool ug3_active;             // Magnetic Strings active?
    bool ug4_active;             // Star-Black Hole active?
    char cosmic_structure[64];   // "Periodic Table", "Stars and Planets", etc.
    uint32_t status;             // 0 = success
};

/**
 * @struct EpochSetRequest
 * @brief Set epoch for calculations (MESSAGE_TYPE: EPOCH_SET)
 */
struct EpochSetRequest {
    int epoch_number;            // Epoch to set (1-5)
    char module_name[64];        // Module to apply to ("all" for global)
    uint32_t flags;              // Reserved
    uint32_t reserved;
};

/**
 * @struct EpochCalculateFURequest
 * @brief Calculate F_U at specific epoch (MESSAGE_TYPE: EPOCH_CALCULATE_F_U)
 */
struct EpochCalculateFURequest {
    int epoch_number;            // Epoch (1-5)
    double rho_vac_UA;           // Universal Aether vacuum density (J/m³)
    double omega_LENR;           // LENR resonance frequency (Hz)
    double sigma_n;              // Neutron cross-section (m²)
    uint32_t flags;              // Calculation flags
    uint32_t reserved;
};

/**
 * @struct EpochCalculateFUResponse
 * @brief Response with epoch-specific F_U calculation
 */
struct EpochCalculateFUResponse {
    int epoch_number;            // Epoch calculated
    double F_U_total_N;          // Total unified field (N)
    double F_core_N;             // Core contribution (N)
    double Ui_sum_N;             // Internal energy sum (N)
    double Fp_sum_N;             // Pressure sum (N)
    char ug_ranges_active[64];   // "Ug1,Ug2,Ug3" or "All", etc.
    uint32_t status;             // 0 = success
    uint32_t reserved;
};

/**
 * @struct EpochGetUgActiveRequest
 * @brief Query which Ug ranges active at epoch (MESSAGE_TYPE: EPOCH_GET_UG_ACTIVE)
 */
struct EpochGetUgActiveRequest {
    int epoch_number;            // Epoch to query (1-5)
    uint32_t flags;              // Reserved
    uint32_t reserved[2];
};

/**
 * @struct EpochGetUgActiveResponse
 * @brief Response with Ug range activation status
 */
struct EpochGetUgActiveResponse {
    int epoch_number;            // Epoch queried
    bool ug1_active;             // Ug1 (Internal Dipole)
    bool ug2_active;             // Ug2 (Heliosphere)
    bool ug3_active;             // Ug3 (Magnetic Strings)
    bool ug4_active;             // Ug4 (Star-Black Hole)
    char activation_context[128]; // Explanation of activation state
    uint32_t status;             // 0 = success
};

/**
 * @struct EpochValidationDataRequest
 * @brief Request epoch validation dataset (MESSAGE_TYPE: EPOCH_VALIDATION_DATA)
 */
struct EpochValidationDataRequest {
    int epoch_number;            // Epoch to validate (1-5, or 0 for all)
    char validation_type[32];    // "gaia", "fermi", "cmb", "all"
    uint32_t flags;              // Reserved
    uint32_t reserved;
};

/**
 * @struct EpochValidationDataResponse
 * @brief Response with validation references
 * 
 * Variable-length JSON payload follows this header with full validation data:
 * {
 *   "epoch": 4,
 *   "validation_targets": [
 *     {"source": "Gaia DR4", "observable": "Sgr A* orbits", "prediction": "Ug4 dominance"},
 *     {"source": "Fermi", "observable": "Solar flares", "prediction": "\u03b1 = 0.001 day⁻¹"}
 *   ]
 * }
 */
struct EpochValidationDataResponse {
    int epoch_number;            // Epoch validated
    uint32_t num_targets;        // Number of validation targets
    uint32_t json_payload_size;  // Size of following JSON (bytes)
    uint32_t status;             // 0 = success
    uint32_t reserved;
    // Variable-length JSON follows
};

// ============================================================================
// CHANNEL INTERFACE
// ============================================================================

/**
 * @class IChannel
 * @brief Abstract interface for IPC channels
 */
class IChannel {
public:
    virtual ~IChannel() = default;
    
    /**
     * Send a message with optional payload
     */
    virtual bool send(const MessageHeader& header, const void* payload = nullptr) = 0;
    
    /**
     * Receive a message (blocking with timeout)
     * @param timeout_ms Timeout in milliseconds (-1 = infinite)
     * @return true if message received
     */
    virtual bool receive(MessageHeader& header, std::vector<uint8_t>& payload, 
                        int timeout_ms = -1) = 0;
    
    /**
     * Check if channel is connected/ready
     */
    virtual bool is_connected() const = 0;
    
    /**
     * Close the channel
     */
    virtual void close() = 0;
    
    /**
     * Get channel name for debugging
     */
    virtual std::string name() const = 0;
};

// ============================================================================
// SHARED MEMORY CHANNEL
// ============================================================================

/**
 * @class SharedMemoryChannel
 * @brief Low-latency IPC via shared memory (for VR real-time data)
 * 
 * Uses a ring buffer in shared memory for lock-free communication.
 * Designed for single-producer single-consumer (SPSC) pattern.
 */
class SharedMemoryChannel : public IChannel {
public:
    /**
     * Create or open a shared memory channel
     * @param name Unique channel name
     * @param buffer_size Size of ring buffer (power of 2)
     * @param create True to create, false to open existing
     */
    SharedMemoryChannel(const std::string& name, 
                        size_t buffer_size = 1024 * 1024,
                        bool create = false);
    
    ~SharedMemoryChannel() override;
    
    bool send(const MessageHeader& header, const void* payload = nullptr) override;
    bool receive(MessageHeader& header, std::vector<uint8_t>& payload, 
                int timeout_ms = -1) override;
    bool is_connected() const override;
    void close() override;
    std::string name() const override { return channel_name_; }
    
    /**
     * Non-blocking send for VR frame updates
     */
    bool try_send(const MessageHeader& header, const void* payload = nullptr);
    
    /**
     * Non-blocking receive
     */
    bool try_receive(MessageHeader& header, std::vector<uint8_t>& payload);
    
private:
    struct RingBuffer {
        std::atomic<uint64_t> write_pos;
        std::atomic<uint64_t> read_pos;
        size_t buffer_size;
        uint8_t data[];  // Flexible array member
    };
    
    std::string channel_name_;
    size_t buffer_size_;
    bool is_creator_;
    
#ifdef _WIN32
    HANDLE mapping_handle_;
    HANDLE mutex_handle_;
#else
    int shm_fd_;
#endif
    
    RingBuffer* buffer_;
    std::atomic<bool> connected_;
};

// ============================================================================
// GRPC CHANNEL (Phase 3 - Full Implementation)
// ============================================================================

/**
 * @class GrpcChannel
 * @brief Structured command channel via gRPC (for complex operations)
 * 
 * Phase 3: Full gRPC implementation for structured physics commands.
 * Uses uqff_service.proto for message definitions.
 */
class GrpcChannel : public IChannel {
public:
    GrpcChannel(const std::string& endpoint);
    ~GrpcChannel() override;
    
    bool send(const MessageHeader& header, const void* payload = nullptr) override;
    bool receive(MessageHeader& header, std::vector<uint8_t>& payload, 
                int timeout_ms = -1) override;
    bool is_connected() const override;
    void close() override;
    std::string name() const override { return endpoint_; }
    
    // gRPC-specific methods (Phase 3)
    bool connect(int timeout_ms = 5000);
    
    // Field calculation via gRPC
    struct FieldResult {
        bool success = false;
        std::string error;
        double F_U = 0, Ug1 = 0, Ug2 = 0, Ug3 = 0, Ug4 = 0;
        double Um = 0, Ubi = 0, g_compressed = 0;
        double compute_time_ms = 0;
    };
    FieldResult calculateField(const std::string& system_name, double r, double t, double theta);
    
    // Service status
    struct ServiceStatus {
        bool healthy = false;
        std::string version;
        uint64_t requests_processed = 0;
        uint64_t uptime_seconds = 0;
    };
    ServiceStatus getStatus();
    
private:
    std::string endpoint_;
    std::atomic<bool> connected_;
    std::mutex channel_mutex_;
    
    // Incoming message queue (for async receive)
    std::queue<std::pair<MessageHeader, std::vector<uint8_t>>> incoming_queue_;
    std::mutex queue_mutex_;
    std::condition_variable queue_cv_;
    
#ifdef USE_GRPC
    std::shared_ptr<grpc::Channel> grpc_channel_;
    std::unique_ptr<uqff::PhysicsService::Stub> stub_;
#endif
};

// ============================================================================
// NAMED PIPE CHANNEL (Cross-Platform)
// ============================================================================

/**
 * @class NamedPipeChannel
 * @brief Cross-platform named pipe IPC channel
 * 
 * Windows: Uses native named pipes (CreateNamedPipe/CreateFile)
 * Linux/macOS: Uses Unix domain sockets (equivalent semantics)
 * 
 * Provides reliable, connection-oriented communication between processes.
 * Suitable for command/response patterns (vs SharedMemory for streaming).
 */
class NamedPipeChannel : public IChannel {
public:
    enum class Mode { SERVER, CLIENT };
    
    /**
     * Create a named pipe channel
     * @param name Pipe name (automatically prefixed with platform path)
     * @param mode SERVER to create/listen, CLIENT to connect
     */
    NamedPipeChannel(const std::string& name, Mode mode);
    ~NamedPipeChannel() override;
    
    bool send(const MessageHeader& header, const void* payload = nullptr) override;
    bool receive(MessageHeader& header, std::vector<uint8_t>& payload, 
                int timeout_ms = -1) override;
    bool is_connected() const override;
    void close() override;
    std::string name() const override { return pipe_name_; }
    
    /**
     * Accept a client connection (SERVER mode only)
     * @param timeout_ms Wait timeout (-1 = infinite)
     * @return true if client connected
     */
    bool accept_connection(int timeout_ms = -1);
    
    /**
     * Connect to server (CLIENT mode only)
     * @param timeout_ms Connection timeout
     * @return true if connected
     */
    bool connect(int timeout_ms = 5000);
    
private:
    std::string pipe_name_;
    Mode mode_;
    std::atomic<bool> connected_;
    
#ifdef _WIN32
    HANDLE pipe_handle_;
    HANDLE connect_event_;  // For overlapped I/O
#else
    int listen_fd_;         // Server: listening socket
    int conn_fd_;           // Connected socket (server: accepted, client: connected)
    std::string socket_path_;
#endif
    
    // Platform-specific initialization
    bool init_server();
    bool init_client();
};

// ============================================================================
// MESSAGE DISPATCHER
// ============================================================================

using MessageHandler = std::function<void(const MessageHeader&, const std::vector<uint8_t>&)>;

/**
 * @class MessageDispatcher
 * @brief Routes incoming messages to registered handlers
 */
class MessageDispatcher {
public:
    /**
     * Register a handler for a message type
     */
    void register_handler(MessageType type, MessageHandler handler);
    
    /**
     * Dispatch a message to its handler
     */
    void dispatch(const MessageHeader& header, const std::vector<uint8_t>& payload);
    
    /**
     * Process messages from a channel until shutdown
     */
    void run(IChannel& channel);
    
    /**
     * Request shutdown
     */
    void stop() { running_ = false; }
    
private:
    std::unordered_map<MessageType, MessageHandler> handlers_;
    std::mutex handlers_mutex_;
    std::atomic<bool> running_{false};
};

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

/**
 * Create a ping message
 */
inline MessageHeader make_ping() {
    return MessageHeader(MessageType::PING);
}

/**
 * Create a shutdown message
 */
inline MessageHeader make_shutdown() {
    return MessageHeader(MessageType::SHUTDOWN);
}

/**
 * Serialize a payload structure to bytes
 */
template<typename T>
std::vector<uint8_t> serialize(const T& payload) {
    std::vector<uint8_t> bytes(sizeof(T));
    std::memcpy(bytes.data(), &payload, sizeof(T));
    return bytes;
}

/**
 * Deserialize bytes to a payload structure
 */
template<typename T>
bool deserialize(const std::vector<uint8_t>& bytes, T& payload) {
    if (bytes.size() < sizeof(T)) return false;
    std::memcpy(&payload, bytes.data(), sizeof(T));
    return true;
}

// ============================================================================
// PIPELINE IPC HELPERS (Feb 24, 2026)
// ============================================================================

/**
 * Create a pipeline process request
 */
inline std::pair<MessageHeader, std::vector<uint8_t>> 
make_pipeline_process(const std::string& object_name, uint32_t timeout_ms = 30000) {
    PipelineProcessRequest req{};
    strncpy(req.object_name, object_name.c_str(), sizeof(req.object_name) - 1);
    req.timeout_ms = timeout_ms;
    req.flags = 0;
    
    MessageHeader hdr(MessageType::PIPELINE_PROCESS, sizeof(req));
    return {hdr, serialize(req)};
}

/**
 * Create a pipeline store request
 */
inline std::pair<MessageHeader, std::vector<uint8_t>> 
make_pipeline_store(const std::string& query_id) {
    PipelineProcessRequest req{};  // Reuse structure for query_id transport
    strncpy(req.callback_id, query_id.c_str(), sizeof(req.callback_id) - 1);
    
    MessageHeader hdr(MessageType::PIPELINE_STORE, sizeof(req));
    return {hdr, serialize(req)};
}

/**
 * Create a pipeline export request
 */
inline std::pair<MessageHeader, std::vector<uint8_t>> 
make_pipeline_export(const std::string& query_id) {
    PipelineProcessRequest req{};
    strncpy(req.callback_id, query_id.c_str(), sizeof(req.callback_id) - 1);
    
    MessageHeader hdr(MessageType::PIPELINE_EXPORT, sizeof(req));
    return {hdr, serialize(req)};
}

/**
 * Parse a pipeline result response
 */
inline bool parse_pipeline_result(const std::vector<uint8_t>& bytes, 
                                   PipelineResultPayload& result,
                                   std::string& json_payload) {
    if (!deserialize(bytes, result)) return false;
    
    // Check if JSON payload follows
    if (result.json_payload_follows && bytes.size() > sizeof(result)) {
        size_t json_start = sizeof(result);
        size_t json_len = bytes.size() - json_start;
        json_payload.assign(reinterpret_cast<const char*>(bytes.data() + json_start), json_len);
    }
    return true;
}

} // namespace IPC
} // namespace UQFF

#endif // UQFF_IPC_H
