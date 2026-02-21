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

} // namespace IPC
} // namespace UQFF

#endif // UQFF_IPC_H
