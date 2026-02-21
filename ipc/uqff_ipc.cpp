/**
 * @file uqff_ipc.cpp
 * @brief Implementation of IPC channels for UQFF Star-Magic
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 1 - IPC Foundation
 */

#include "uqff_ipc.h"
#include <iostream>
#include <thread>
#include <cstring>
#include <unordered_map>

namespace UQFF {
namespace IPC {

// ============================================================================
// SHARED MEMORY CHANNEL IMPLEMENTATION
// ============================================================================

SharedMemoryChannel::SharedMemoryChannel(const std::string& name, 
                                         size_t buffer_size,
                                         bool create)
    : channel_name_(name)
    , buffer_size_(buffer_size)
    , is_creator_(create)
    , buffer_(nullptr)
    , connected_(false)
#ifdef _WIN32
    , mapping_handle_(nullptr)
    , mutex_handle_(nullptr)
#else
    , shm_fd_(-1)
#endif
{
    const size_t total_size = sizeof(RingBuffer) + buffer_size;
    
#ifdef _WIN32
    std::wstring wide_name(name.begin(), name.end());
    std::wstring mapping_name = L"Local\\UQFF_SHM_" + wide_name;
    std::wstring mutex_name = L"Local\\UQFF_MTX_" + wide_name;
    
    if (create) {
        // Create new shared memory
        mapping_handle_ = CreateFileMappingW(
            INVALID_HANDLE_VALUE,
            nullptr,
            PAGE_READWRITE,
            0,
            static_cast<DWORD>(total_size),
            mapping_name.c_str()
        );
        
        if (mapping_handle_ && GetLastError() != ERROR_ALREADY_EXISTS) {
            buffer_ = static_cast<RingBuffer*>(
                MapViewOfFile(mapping_handle_, FILE_MAP_ALL_ACCESS, 0, 0, total_size)
            );
            
            if (buffer_) {
                // Initialize ring buffer
                buffer_->write_pos.store(0);
                buffer_->read_pos.store(0);
                buffer_->buffer_size = buffer_size;
            }
        }
    } else {
        // Open existing shared memory
        mapping_handle_ = OpenFileMappingW(
            FILE_MAP_ALL_ACCESS,
            FALSE,
            mapping_name.c_str()
        );
        
        if (mapping_handle_) {
            buffer_ = static_cast<RingBuffer*>(
                MapViewOfFile(mapping_handle_, FILE_MAP_ALL_ACCESS, 0, 0, total_size)
            );
        }
    }
    
    // Create/open mutex for synchronization
    mutex_handle_ = CreateMutexW(nullptr, FALSE, mutex_name.c_str());
    
    connected_ = (buffer_ != nullptr && mutex_handle_ != nullptr);
    
#else
    // POSIX shared memory
    std::string shm_name = "/uqff_shm_" + name;
    
    if (create) {
        shm_fd_ = shm_open(shm_name.c_str(), O_CREAT | O_RDWR, 0666);
        if (shm_fd_ >= 0) {
            ftruncate(shm_fd_, total_size);
        }
    } else {
        shm_fd_ = shm_open(shm_name.c_str(), O_RDWR, 0666);
    }
    
    if (shm_fd_ >= 0) {
        buffer_ = static_cast<RingBuffer*>(
            mmap(nullptr, total_size, PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd_, 0)
        );
        
        if (buffer_ == MAP_FAILED) {
            buffer_ = nullptr;
        } else if (create) {
            buffer_->write_pos.store(0);
            buffer_->read_pos.store(0);
            buffer_->buffer_size = buffer_size;
        }
    }
    
    connected_ = (buffer_ != nullptr);
#endif
    
    if (connected_) {
        std::cout << "[IPC] SharedMemoryChannel '" << name << "' " 
                  << (create ? "created" : "opened") << " successfully" << std::endl;
    } else {
        std::cerr << "[IPC] Failed to " << (create ? "create" : "open") 
                  << " SharedMemoryChannel '" << name << "'" << std::endl;
    }
}

SharedMemoryChannel::~SharedMemoryChannel() {
    close();
}

void SharedMemoryChannel::close() {
    if (!connected_) return;
    
    connected_ = false;
    
#ifdef _WIN32
    if (buffer_) {
        UnmapViewOfFile(buffer_);
        buffer_ = nullptr;
    }
    if (mapping_handle_) {
        CloseHandle(mapping_handle_);
        mapping_handle_ = nullptr;
    }
    if (mutex_handle_) {
        CloseHandle(mutex_handle_);
        mutex_handle_ = nullptr;
    }
#else
    if (buffer_) {
        munmap(buffer_, sizeof(RingBuffer) + buffer_size_);
        buffer_ = nullptr;
    }
    if (shm_fd_ >= 0) {
        ::close(shm_fd_);
        if (is_creator_) {
            std::string shm_name = "/uqff_shm_" + channel_name_;
            shm_unlink(shm_name.c_str());
        }
        shm_fd_ = -1;
    }
#endif
    
    std::cout << "[IPC] SharedMemoryChannel '" << channel_name_ << "' closed" << std::endl;
}

bool SharedMemoryChannel::is_connected() const {
    return connected_;
}

bool SharedMemoryChannel::send(const MessageHeader& header, const void* payload) {
    if (!connected_ || !buffer_) return false;
    
#ifdef _WIN32
    WaitForSingleObject(mutex_handle_, INFINITE);
#endif
    
    const size_t msg_size = sizeof(MessageHeader) + header.payload_size;
    const uint64_t write_pos = buffer_->write_pos.load(std::memory_order_relaxed);
    const uint64_t read_pos = buffer_->read_pos.load(std::memory_order_acquire);
    
    // Check available space (ring buffer)
    const size_t used = (write_pos >= read_pos) ? 
        (write_pos - read_pos) : 
        (buffer_->buffer_size - read_pos + write_pos);
    
    if (used + msg_size + sizeof(size_t) > buffer_->buffer_size) {
#ifdef _WIN32
        ReleaseMutex(mutex_handle_);
#endif
        return false; // Buffer full
    }
    
    // Write message size
    const size_t offset = write_pos % buffer_->buffer_size;
    std::memcpy(buffer_->data + offset, &msg_size, sizeof(size_t));
    
    // Write header
    size_t pos = (offset + sizeof(size_t)) % buffer_->buffer_size;
    std::memcpy(buffer_->data + pos, &header, sizeof(MessageHeader));
    
    // Write payload
    if (payload && header.payload_size > 0) {
        pos = (pos + sizeof(MessageHeader)) % buffer_->buffer_size;
        std::memcpy(buffer_->data + pos, payload, header.payload_size);
    }
    
    buffer_->write_pos.store(write_pos + sizeof(size_t) + msg_size, 
                             std::memory_order_release);
    
#ifdef _WIN32
    ReleaseMutex(mutex_handle_);
#endif
    
    return true;
}

bool SharedMemoryChannel::try_send(const MessageHeader& header, const void* payload) {
#ifdef _WIN32
    if (WaitForSingleObject(mutex_handle_, 0) != WAIT_OBJECT_0) {
        return false;
    }
    bool result = send(header, payload);
    // Mutex already released in send()
    return result;
#else
    return send(header, payload);
#endif
}

bool SharedMemoryChannel::receive(MessageHeader& header, std::vector<uint8_t>& payload, 
                                  int timeout_ms) {
    if (!connected_ || !buffer_) return false;
    
    auto start = std::chrono::steady_clock::now();
    
    while (true) {
        const uint64_t write_pos = buffer_->write_pos.load(std::memory_order_acquire);
        const uint64_t read_pos = buffer_->read_pos.load(std::memory_order_relaxed);
        
        if (write_pos > read_pos) {
            // Data available
#ifdef _WIN32
            WaitForSingleObject(mutex_handle_, INFINITE);
#endif
            
            // Re-check after acquiring lock
            const size_t offset = read_pos % buffer_->buffer_size;
            
            // Read message size
            size_t msg_size;
            std::memcpy(&msg_size, buffer_->data + offset, sizeof(size_t));
            
            // Read header
            size_t pos = (offset + sizeof(size_t)) % buffer_->buffer_size;
            std::memcpy(&header, buffer_->data + pos, sizeof(MessageHeader));
            
            // Validate
            if (!header.is_valid()) {
#ifdef _WIN32
                ReleaseMutex(mutex_handle_);
#endif
                return false;
            }
            
            // Read payload
            if (header.payload_size > 0) {
                payload.resize(header.payload_size);
                pos = (pos + sizeof(MessageHeader)) % buffer_->buffer_size;
                std::memcpy(payload.data(), buffer_->data + pos, header.payload_size);
            } else {
                payload.clear();
            }
            
            buffer_->read_pos.store(read_pos + sizeof(size_t) + msg_size, 
                                   std::memory_order_release);
            
#ifdef _WIN32
            ReleaseMutex(mutex_handle_);
#endif
            return true;
        }
        
        // Check timeout
        if (timeout_ms >= 0) {
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - start).count();
            if (elapsed >= timeout_ms) {
                return false;
            }
        }
        
        // Sleep briefly before retry
        std::this_thread::sleep_for(std::chrono::microseconds(100));
    }
}

bool SharedMemoryChannel::try_receive(MessageHeader& header, std::vector<uint8_t>& payload) {
    return receive(header, payload, 0);
}

// ============================================================================
// GRPC CHANNEL IMPLEMENTATION (Phase 3 - Full gRPC)
// ============================================================================

GrpcChannel::GrpcChannel(const std::string& endpoint)
    : endpoint_(endpoint)
    , connected_(false)
{
#ifdef USE_GRPC
    // Create gRPC channel with default credentials
    grpc_channel_ = grpc::CreateChannel(endpoint, grpc::InsecureChannelCredentials());
    stub_ = uqff::PhysicsService::NewStub(grpc_channel_);
    
    // Check connection immediately
    auto state = grpc_channel_->GetState(true);
    if (state == GRPC_CHANNEL_READY || state == GRPC_CHANNEL_IDLE) {
        connected_ = true;
        std::cout << "[IPC/gRPC] Connected to: " << endpoint << std::endl;
    } else {
        std::cout << "[IPC/gRPC] Created channel (pending connection): " << endpoint << std::endl;
    }
#else
    std::cout << "[IPC/gRPC] gRPC support disabled (USE_GRPC not defined)" << std::endl;
    std::cout << "[IPC/gRPC] Endpoint: " << endpoint << " (stub mode)" << std::endl;
    connected_ = true;  // Stub mode
#endif
}

GrpcChannel::~GrpcChannel() {
    close();
}

bool GrpcChannel::connect(int timeout_ms) {
#ifdef USE_GRPC
    if (!grpc_channel_) return false;
    
    auto deadline = std::chrono::system_clock::now() + std::chrono::milliseconds(timeout_ms);
    bool success = grpc_channel_->WaitForConnected(deadline);
    connected_ = success;
    
    if (success) {
        std::cout << "[IPC/gRPC] Connection established: " << endpoint_ << std::endl;
    } else {
        std::cerr << "[IPC/gRPC] Connection timeout: " << endpoint_ << std::endl;
    }
    return success;
#else
    connected_ = true;
    return true;
#endif
}

bool GrpcChannel::send(const MessageHeader& header, const void* payload) {
    if (!connected_) return false;
    
#ifdef USE_GRPC
    // Convert IPC message to gRPC call based on message type
    if (header.type == MessageType::CALCULATE_FIELD && payload && header.payload_size >= sizeof(CalculateFieldRequest)) {
        const auto* req = static_cast<const CalculateFieldRequest*>(payload);
        
        grpc::ClientContext context;
        uqff::FieldRequest grpc_req;
        uqff::FieldResponse grpc_resp;
        
        grpc_req.set_system_name(req->system_name);
        grpc_req.set_r(req->r);
        grpc_req.set_t(req->t);
        grpc_req.set_theta(req->theta);
        grpc_req.set_flags(req->flags);
        
        grpc::Status status = stub_->CalculateField(&context, grpc_req, &grpc_resp);
        
        if (status.ok() && grpc_resp.success()) {
            // Queue the response for receive()
            MessageHeader resp_header(MessageType::RESPONSE_DATA);
            CalculateFieldResponse resp_payload;
            resp_payload.status = 0;
            resp_payload.F_U = grpc_resp.f_u();
            resp_payload.Ug1 = grpc_resp.ug1();
            resp_payload.Ug2 = grpc_resp.ug2();
            resp_payload.Ug3 = grpc_resp.ug3();
            resp_payload.Ug4 = grpc_resp.ug4();
            resp_payload.Um = grpc_resp.um();
            resp_payload.Ub = grpc_resp.ubi();
            resp_payload.compressed_g = grpc_resp.g_compressed();
            
            std::vector<uint8_t> resp_bytes(sizeof(resp_payload));
            std::memcpy(resp_bytes.data(), &resp_payload, sizeof(resp_payload));
            
            {
                std::lock_guard<std::mutex> lock(queue_mutex_);
                incoming_queue_.push({resp_header, std::move(resp_bytes)});
            }
            queue_cv_.notify_one();
            return true;
        } else {
            std::cerr << "[IPC/gRPC] CalculateField failed: " << status.error_message() << std::endl;
            return false;
        }
    }
    
    // Log unhandled message types
    std::cout << "[IPC/gRPC] send: type=" << static_cast<uint32_t>(header.type) 
              << " size=" << header.payload_size << std::endl;
    return true;
#else
    // Stub mode: just log
    std::cout << "[IPC/gRPC] STUB send: type=" << static_cast<uint32_t>(header.type) 
              << " size=" << header.payload_size << std::endl;
    return true;
#endif
}

bool GrpcChannel::receive(MessageHeader& header, std::vector<uint8_t>& payload, 
                         int timeout_ms) {
    if (!connected_) return false;
    
    std::unique_lock<std::mutex> lock(queue_mutex_);
    
    if (timeout_ms < 0) {
        // Block indefinitely
        queue_cv_.wait(lock, [this] { return !incoming_queue_.empty(); });
    } else if (timeout_ms > 0) {
        // Wait with timeout
        bool got_msg = queue_cv_.wait_for(lock, std::chrono::milliseconds(timeout_ms),
                                          [this] { return !incoming_queue_.empty(); });
        if (!got_msg) return false;
    } else {
        // Non-blocking
        if (incoming_queue_.empty()) return false;
    }
    
    if (!incoming_queue_.empty()) {
        auto& front = incoming_queue_.front();
        header = front.first;
        payload = std::move(front.second);
        incoming_queue_.pop();
        return true;
    }
    
    return false;
}

bool GrpcChannel::is_connected() const {
#ifdef USE_GRPC
    if (!grpc_channel_) return false;
    auto state = grpc_channel_->GetState(false);
    return state == GRPC_CHANNEL_READY || state == GRPC_CHANNEL_IDLE;
#else
    return connected_;
#endif
}

void GrpcChannel::close() {
    if (connected_) {
        std::cout << "[IPC/gRPC] Channel closed: " << endpoint_ << std::endl;
        connected_ = false;
#ifdef USE_GRPC
        stub_.reset();
        grpc_channel_.reset();
#endif
    }
}

GrpcChannel::FieldResult GrpcChannel::calculateField(const std::string& system_name, 
                                                     double r, double t, double theta) {
    FieldResult result;
    
#ifdef USE_GRPC
    if (!stub_) {
        result.error = "gRPC stub not initialized";
        return result;
    }
    
    grpc::ClientContext context;
    uqff::FieldRequest request;
    uqff::FieldResponse response;
    
    request.set_system_name(system_name);
    request.set_r(r);
    request.set_t(t);
    request.set_theta(theta);
    
    auto start = std::chrono::high_resolution_clock::now();
    grpc::Status status = stub_->CalculateField(&context, request, &response);
    auto end = std::chrono::high_resolution_clock::now();
    
    result.compute_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    if (status.ok()) {
        result.success = response.success();
        result.error = response.error_message();
        result.F_U = response.f_u();
        result.Ug1 = response.ug1();
        result.Ug2 = response.ug2();
        result.Ug3 = response.ug3();
        result.Ug4 = response.ug4();
        result.Um = response.um();
        result.Ubi = response.ubi();
        result.g_compressed = response.g_compressed();
    } else {
        result.error = status.error_message();
    }
#else
    // Stub: Return placeholder
    result.success = true;
    result.F_U = 1.83e71;  // Placeholder base force
    result.compute_time_ms = 1.0;
#endif
    
    return result;
}

GrpcChannel::ServiceStatus GrpcChannel::getStatus() {
    ServiceStatus status;
    
#ifdef USE_GRPC
    if (!stub_) return status;
    
    grpc::ClientContext context;
    uqff::StatusRequest request;
    uqff::StatusResponse response;
    
    request.set_include_stats(true);
    
    grpc::Status grpc_status = stub_->GetStatus(&context, request, &response);
    
    if (grpc_status.ok()) {
        status.healthy = response.healthy();
        status.version = response.version();
        status.requests_processed = response.requests_processed();
        status.uptime_seconds = response.uptime_seconds();
    }
#else
    status.healthy = true;
    status.version = "3.0-stub";
#endif
    
    return status;
}

// ============================================================================
// MESSAGE DISPATCHER IMPLEMENTATION
// ============================================================================

void MessageDispatcher::register_handler(MessageType type, MessageHandler handler) {
    std::lock_guard<std::mutex> lock(handlers_mutex_);
    handlers_[type] = std::move(handler);
}

void MessageDispatcher::dispatch(const MessageHeader& header, 
                                 const std::vector<uint8_t>& payload) {
    std::lock_guard<std::mutex> lock(handlers_mutex_);
    
    auto it = handlers_.find(header.type);
    if (it != handlers_.end()) {
        it->second(header, payload);
    } else {
        std::cerr << "[IPC] No handler for message type: " 
                  << static_cast<uint32_t>(header.type) << std::endl;
    }
}

void MessageDispatcher::run(IChannel& channel) {
    running_ = true;
    std::cout << "[IPC] MessageDispatcher started on channel: " << channel.name() << std::endl;
    
    MessageHeader header;
    std::vector<uint8_t> payload;
    
    while (running_) {
        if (channel.receive(header, payload, 100)) {
            // Handle shutdown
            if (header.type == MessageType::SHUTDOWN) {
                std::cout << "[IPC] Received SHUTDOWN message" << std::endl;
                running_ = false;
                break;
            }
            
            // Handle ping
            if (header.type == MessageType::PING) {
                MessageHeader pong(MessageType::PONG);
                channel.send(pong);
                continue;
            }
            
            // Dispatch to handler
            dispatch(header, payload);
        }
    }
    
    std::cout << "[IPC] MessageDispatcher stopped" << std::endl;
}

} // namespace IPC
} // namespace UQFF
