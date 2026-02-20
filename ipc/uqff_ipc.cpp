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
// GRPC CHANNEL STUB IMPLEMENTATION
// ============================================================================

GrpcChannel::GrpcChannel(const std::string& endpoint)
    : endpoint_(endpoint)
    , connected_(false)
{
    // Stub implementation - full gRPC integration in Phase 2
    std::cout << "[IPC] GrpcChannel stub created for endpoint: " << endpoint << std::endl;
    std::cout << "[IPC] NOTE: Full gRPC implementation pending Phase 2" << std::endl;
    connected_ = true;  // Stub always "connected"
}

GrpcChannel::~GrpcChannel() {
    close();
}

bool GrpcChannel::send(const MessageHeader& header, const void* payload) {
    if (!connected_) return false;
    
    // Stub: Log the message
    std::cout << "[IPC/gRPC] STUB send: type=" << static_cast<uint32_t>(header.type) 
              << " size=" << header.payload_size << std::endl;
    
    // In Phase 2, this will serialize and send via gRPC
    return true;
}

bool GrpcChannel::receive(MessageHeader& header, std::vector<uint8_t>& payload, 
                         int timeout_ms) {
    if (!connected_) return false;
    
    // Stub: No messages available
    // In Phase 2, this will receive from gRPC stream
    std::this_thread::sleep_for(std::chrono::milliseconds(timeout_ms > 0 ? timeout_ms : 100));
    return false;
}

bool GrpcChannel::is_connected() const {
    return connected_;
}

void GrpcChannel::close() {
    if (connected_) {
        std::cout << "[IPC] GrpcChannel closed: " << endpoint_ << std::endl;
        connected_ = false;
    }
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
