/**
 * @file uqff_ipc.h
 * @brief UQFF Inter-Process Communication Protocol (v3.1)
 *
 * Defines the complete message protocol for all UQFF process boundaries:
 *   source2.cpp  ←→  MAIN_1_CoAnQi.exe
 *   source2.cpp  ←→  QCalc.py / CondensedPhysics.py / CondensedPhysics2.py
 *   source2.cpp  ←→  uqff_server.js (REST, port 3141)
 *   physics_backend.cpp  ←→  source2(HEAD PROGRAM).cpp
 *
 * Transport options (selected at build time):
 *   UQFF_IPC_NAMED_PIPE  — \\.\pipe\StarMagic_UQFF  (Windows default)
 *   UQFF_IPC_SHARED_MEM  — mapped file fallback
 *   UQFF_IPC_SOCKET      — TCP loopback (cross-platform)
 *
 * 45 message types across 6 functional groups:
 *   0x01xx — Query / Compute
 *   0x02xx — Data transfer
 *   0x03xx — Control / Session
 *   0x04xx — Result / Response
 *   0x05xx — Error / Status
 *   0x0Axx — Extended / Plugin
 *
 * Author: Daniel T. Murphy
 * Date: March 2026 (v3.1 — GrokThread7b0e961f integration)
 * Constants: shared_constants.h (GrokThread7b0e namespace)
 */

#ifndef UQFF_IPC_H
#define UQFF_IPC_H

#include <cstdint>
#include <string>
#include <vector>
#include <functional>

#ifdef _WIN32
#  include <windows.h>
#else
#  include <sys/mman.h>
#  include <sys/stat.h>
#  include <fcntl.h>
#  include <unistd.h>
#endif

namespace UQFF {
namespace IPC {

// ─────────────────────────────────────────────────────────────
// NAMED PIPE / TRANSPORT CONSTANTS
// ─────────────────────────────────────────────────────────────
#ifdef _WIN32
constexpr const char*    PIPE_NAME        = "\\\\.\\pipe\\StarMagic_UQFF";
constexpr DWORD          PIPE_BUFFER_SIZE = 65536;   // 64 KB per direction
constexpr DWORD          PIPE_TIMEOUT_MS  = 5000;
#else
constexpr const char*    SOCKET_HOST      = "127.0.0.1";
constexpr uint16_t       SOCKET_PORT      = 41000;   // loopback IPC port
#endif

constexpr const char*    SHM_NAME         = "/StarMagic_UQFF_Mem";
constexpr size_t         SHM_SIZE         = 512 * 1024; // 512 KB shared region

// ─────────────────────────────────────────────────────────────
// MESSAGE TYPE ENUMERATION  (45 types, 6 groups)
// ─────────────────────────────────────────────────────────────
enum class MessageType : uint16_t {

    // ── Group 0x01xx : Query / Compute ──────────────────────
    COMPUTE_SINGLE        = 0x0101, ///< Compute single astronomical system
    COMPUTE_BATCH         = 0x0102, ///< Compute all registered systems (parallel)
    QUERY_PHYSICS_TERM    = 0x0103, ///< Look up a named PhysicsTerm result
    QUERY_MODULE          = 0x0104, ///< Run a named SOURCE module
    COMPUTE_F_U_BI_I      = 0x0105, ///< Buoyancy integral F_U_Bi_i
    COMPUTE_COMPRESSED_G  = 0x0106, ///< Compressed MUGE gravity
    COMPUTE_RESONANCE_G   = 0x0107, ///< Resonance MUGE gravity
    COMPUTE_SOURCE4       = 0x0108, ///< Full SOURCE4 Unified Field validation
    COMPUTE_26D           = 0x0109, ///< 26-Layer Cosmic Quantum Egg simulation

    // ── Group 0x02xx : Data Transfer ────────────────────────
    SEND_BODY_PARAMS      = 0x0201, ///< Push SystemParams / bodies CSV row
    SEND_DATASET          = 0x0202, ///< Push arbitrary JSON dataset
    SEND_BINARY_BLOB      = 0x0203, ///< Push raw binary payload (time-series, etc.)
    REQUEST_BODY_PARAMS   = 0x0204, ///< Pull SystemParams by system name
    REQUEST_CSV_ROW       = 0x0205, ///< Pull raw CSV row from bodies_*.csv
    STREAM_TIME_SERIES    = 0x0206, ///< Chunked time-series stream (multi-frame)
    SYNC_SHARED_MEM       = 0x0207, ///< Invalidate / flush shared memory region

    // ── Group 0x03xx : Control / Session ────────────────────
    SESSION_START         = 0x0301, ///< Open a compute session (returns session_id)
    SESSION_END           = 0x0302, ///< Close session, flush caches
    PING                  = 0x0303, ///< Heartbeat / liveness check
    CANCEL                = 0x0304, ///< Abort in-flight computation
    SET_PARAMETER         = 0x0305, ///< Set a named dynamic parameter
    GET_PARAMETER         = 0x0306, ///< Get a named dynamic parameter
    REGISTER_SYSTEM       = 0x0307, ///< Register a new astrophysical system at runtime
    LOAD_MODULE           = 0x0308, ///< Dynamically load a UQFF module
    EXPORT_STATE          = 0x0309, ///< Trigger module exportState() to file
    IMPORT_STATE          = 0x030A, ///< Load module state from file

    // ── Group 0x04xx : Result / Response ────────────────────
    RESPONSE_DATA         = 0x0401, ///< Scalar / vector result payload
    RESPONSE_JSON         = 0x0402, ///< Full JSON result object
    RESPONSE_EQUATION     = 0x0403, ///< Long-form equation string (CondensedPhysics)
    RESPONSE_EQUATION_SET = 0x0404, ///< All solvable equations for query
    RESPONSE_BINARY       = 0x0405, ///< Binary result blob
    RESPONSE_PONG         = 0x0406, ///< Heartbeat reply (to PING)
    RESPONSE_PARTIAL      = 0x0407, ///< Chunked response (indicates more to follow)
    RESPONSE_FINAL        = 0x0408, ///< Last chunk of a multi-chunk response

    // ── Group 0x05xx : Error / Status ───────────────────────
    ERROR_GENERAL         = 0x0501, ///< Generic error (see error_code field)
    ERROR_UNKNOWN_SYSTEM  = 0x0502, ///< System name not found in registry
    ERROR_UNKNOWN_MODULE  = 0x0503, ///< Module / SOURCE not found
    ERROR_TIMEOUT         = 0x0504, ///< Computation timeout
    ERROR_INVALID_PARAMS  = 0x0505, ///< Parameter validation failure
    ERROR_IPC_TRANSPORT   = 0x0506, ///< Pipe / socket / shm transport failure
    STATUS_BUSY           = 0x0507, ///< Server busy, retry after retry_ms
    STATUS_READY          = 0x0508, ///< Server ready (startup acknowledgement)

    // ── Group 0x0Axx : Extended / Plugin ────────────────────
    CP2_TRIGGER           = 0x0A00, ///< Route query to CondensedPhysics2 calculator
    CP2_RESPONSE          = 0x0A01, ///< CondensedPhysics2 result payload
    WOLFRAM_EVAL          = 0x0A02, ///< Forward expression to Wolfram WSTP kernel
    WOLFRAM_RESULT        = 0x0A03, ///< Wolfram evaluation result
    GROK_QUERY            = 0x0A04, ///< Forward query to Grok xAI endpoint
    GROK_RESPONSE         = 0x0A05, ///< Grok xAI response
    PLUGIN_CUSTOM         = 0x0A06, ///< User-defined plugin message (opaque payload)
};

// ─────────────────────────────────────────────────────────────
// MESSAGE HEADER
// ─────────────────────────────────────────────────────────────
#pragma pack(push, 1)
struct MessageHeader {
    uint32_t    magic;          ///< 0x55514646 == 'UQFF'
    uint16_t    version;        ///< Protocol version: 0x0301
    MessageType type;           ///< Message type (see enum)
    uint32_t    session_id;     ///< Session cookie (SESSION_START → SESSION_END)
    uint32_t    sequence;       ///< Monotone sequence number (per direction)
    uint32_t    payload_length; ///< Bytes following this header
    uint32_t    checksum;       ///< CRC32 of payload bytes
};
#pragma pack(pop)

static_assert(sizeof(MessageHeader) == 22, "MessageHeader must be 22 bytes");

constexpr uint32_t HEADER_MAGIC   = 0x55514646u; // 'UQFF'
constexpr uint16_t PROTOCOL_VER   = 0x0301u;

// ─────────────────────────────────────────────────────────────
// COMMON PAYLOAD STRUCTURES
// ─────────────────────────────────────────────────────────────

/// Payload for COMPUTE_SINGLE / SEND_BODY_PARAMS
struct ComputePayload {
    char    system_name[64];    ///< Null-terminated system identifier
    double  mass_kg;            ///< kg
    double  radius_m;           ///< m
    double  distance_m;         ///< m
    double  redshift;           ///< dimensionless
    double  t_n;                ///< normalised time
    double  theta;              ///< angle parameter (rad)
    double  U_UA;               ///< buoyancy-absorption factor (calibrated: 0.0001)
    double  SSq;                ///< quantum coupling ([SSq]=0.57)
    double  kappa;              ///< calibrated decay constant (0.0005/day)
    double  F_rel;              ///< relativistic force (4.31e33 N, 2024 LEP)
    double  extra[8];           ///< spare fields for future use
};

/// Payload for RESPONSE_DATA
struct ResponsePayload {
    char    system_name[64];
    double  F_U_Bi_i;           ///< Buoyancy integral result (N)
    double  compressed_g;       ///< MUGE compressed gravity (m/s²)
    double  resonance_g;        ///< MUGE resonance gravity (m/s²)
    double  Q_wave;             ///< Q_wave energy density (J/m³)
    double  rho_vac;            ///< Vacuum density (kg/m³)
    double  chi_squared;        ///< Goodness-of-fit statistic
    double  elapsed_ms;         ///< Wall-clock compute time
    uint8_t status;             ///< 0 = OK, non-zero = error code
    char    message[128];       ///< Human-readable status / equation string
};

/// Payload for SET_PARAMETER / GET_PARAMETER
struct ParameterPayload {
    char    key[64];
    double  value;
    char    module_name[64];    ///< Target module (empty = global)
};

/// Error detail payload (follows MessageHeader when type is ERROR_*)
struct ErrorPayload {
    uint32_t error_code;
    char     description[192];
    double   retry_after_ms;    ///< > 0 when STATUS_BUSY
};

// ─────────────────────────────────────────────────────────────
// FRAMING HELPERS
// ─────────────────────────────────────────────────────────────
inline uint32_t compute_crc32(const uint8_t* data, size_t length) {
    uint32_t crc = 0xFFFFFFFFu;
    for (size_t i = 0; i < length; ++i) {
        crc ^= data[i];
        for (int j = 0; j < 8; ++j)
            crc = (crc >> 1) ^ (0xEDB88320u & -(crc & 1u));
    }
    return crc ^ 0xFFFFFFFFu;
}

/// Build a complete framed message (header + payload) into out_buf.
/// Returns total bytes written.
inline size_t build_message(
        MessageType          type,
        uint32_t             session_id,
        uint32_t             sequence,
        const void*          payload,
        uint32_t             payload_len,
        std::vector<uint8_t>& out_buf)
{
    out_buf.resize(sizeof(MessageHeader) + payload_len);
    auto* hdr = reinterpret_cast<MessageHeader*>(out_buf.data());
    hdr->magic          = HEADER_MAGIC;
    hdr->version        = PROTOCOL_VER;
    hdr->type           = type;
    hdr->session_id     = session_id;
    hdr->sequence       = sequence;
    hdr->payload_length = payload_len;
    hdr->checksum       = payload_len
        ? compute_crc32(static_cast<const uint8_t*>(payload), payload_len)
        : 0u;
    if (payload_len && payload)
        std::memcpy(out_buf.data() + sizeof(MessageHeader), payload, payload_len);
    return out_buf.size();
}

/// Validate an inbound buffer (at least header must be present).
/// Returns false if magic, version, or checksum does not match.
inline bool validate_message(const uint8_t* buf, size_t buf_len) {
    if (buf_len < sizeof(MessageHeader)) return false;
    const auto* hdr = reinterpret_cast<const MessageHeader*>(buf);
    if (hdr->magic   != HEADER_MAGIC) return false;
    if (hdr->version != PROTOCOL_VER) return false;
    if (buf_len < sizeof(MessageHeader) + hdr->payload_length) return false;
    if (hdr->payload_length > 0) {
        uint32_t actual = compute_crc32(buf + sizeof(MessageHeader), hdr->payload_length);
        if (actual != hdr->checksum) return false;
    }
    return true;
}

// ─────────────────────────────────────────────────────────────
// NAMED PIPE TRANSPORT (Windows)
// ─────────────────────────────────────────────────────────────
#ifdef _WIN32
/// Create the server-side named pipe instance.
/// Call once at startup in MAIN_1_CoAnQi.exe / source2.cpp server thread.
inline HANDLE create_server_pipe() {
    return CreateNamedPipeA(
        PIPE_NAME,
        PIPE_ACCESS_DUPLEX | FILE_FLAG_OVERLAPPED,
        PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT,
        PIPE_UNLIMITED_INSTANCES,
        PIPE_BUFFER_SIZE, PIPE_BUFFER_SIZE,
        PIPE_TIMEOUT_MS,
        nullptr);
}

/// Connect client side (source2.cpp → MAIN_1).
inline HANDLE connect_client_pipe() {
    if (!WaitNamedPipeA(PIPE_NAME, PIPE_TIMEOUT_MS)) return INVALID_HANDLE_VALUE;
    return CreateFileA(
        PIPE_NAME,
        GENERIC_READ | GENERIC_WRITE,
        0, nullptr,
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL | FILE_FLAG_OVERLAPPED,
        nullptr);
}
#endif // _WIN32

} // namespace IPC
} // namespace UQFF

#endif // UQFF_IPC_H
