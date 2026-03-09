/**
 * @file python_bridge.h
 * @brief Python ↔ C++ bridge for UQFF calculators
 *
 * Provides a thin C++ interface for launching and communicating with
 * the Python UQFF calculators (CondensedPhysics.py, CondensedPhysics2.py,
 * QCalc.py) over the UQFF IPC protocol (uqff_ipc.h).
 *
 * Two usage modes:
 *   1. Subprocess mode (default) — spawns Python via QProcess / CreateProcess,
 *      exchanges JSON over stdin/stdout. No Python embedding overhead.
 *   2. Named-pipe mode — communicates with a pre-started uqff_ipc.py server
 *      using the binary UQFF IPC framing defined in uqff_ipc.h.
 *
 * Calculator routing:
 *   "standard"      → QCalc.py::UnifiedFieldSolver           (~920 ms)
 *   "cp2:<trigger>" → CondensedPhysics2.py calculator class  (~2500 ms)
 *   "cp1:<trigger>" → CondensedPhysics.py calculator class   (~30 s – avoid interactive)
 *
 * Author: Daniel T. Murphy
 * Date: March 2026 (v1.0 — initial creation, GrokThread7b0e961f)
 */

#ifndef PYTHON_BRIDGE_H
#define PYTHON_BRIDGE_H

#include "uqff_ipc.h"

#include <string>
#include <vector>
#include <functional>
#include <memory>
#include <cstdint>

#ifdef _WIN32
#  include <windows.h>
#else
#  include <sys/wait.h>
#  include <unistd.h>
#endif

namespace UQFF {
namespace Bridge {

// ─────────────────────────────────────────────────────────────
// CALCULATOR TARGET IDENTIFIERS
// ─────────────────────────────────────────────────────────────
enum class CalculatorTarget {
    QCalc,              ///< QCalc.py — fast unified field solver
    CP2,                ///< CondensedPhysics2.py — extended UQFF (548 classes)
    CP1,                ///< CondensedPhysics.py  — primary calculator (176 classes)
    UQFF_SERVER,        ///< uqff_server.js REST API (port 3141)
};

/// Map an IPC CP2_TRIGGER description string to the right target.
/// Mirrors the routing logic in ipc_pipeline_handler.h.
inline CalculatorTarget route_trigger(const std::string& trigger) {
    static const std::vector<std::string> CP2_KEYWORDS = {
        "Orb10", "Orb11", "Orb12", "Orb13", "Orb14", "Orb15", "Orb16",
        "Red Mercury", "Plasmoid", "UFEQFET", "Monte Carlo", "Relativistic",
        "ResonanceGravity", "AsymCap", "FractalTime", "VacuumProb",
        "26Layer", "CompressedGravity", "BuoyancyProof",
        "ShapiroWilk", "qwaveNormal", "rotorCS", "H2OH2PES",
        "DPMfreqMUGE", "BECalpha", "complexUi", "SuperCondUI",
    };
    for (const auto& kw : CP2_KEYWORDS)
        if (trigger.find(kw) != std::string::npos) return CalculatorTarget::CP2;
    return CalculatorTarget::QCalc;
}

// ─────────────────────────────────────────────────────────────
// RESULT TYPE
// ─────────────────────────────────────────────────────────────
struct BridgeResult {
    bool        success       = false;
    std::string json_payload;           ///< Raw JSON string from Python
    std::string equation_text;          ///< Long-form equation (if provided)
    double      elapsed_ms    = 0.0;
    std::string error_message;

    /// Convenience: check for success without inspecting json_payload.
    explicit operator bool() const { return success; }
};

// ─────────────────────────────────────────────────────────────
// SUBPROCESS-MODE BRIDGE
// ─────────────────────────────────────────────────────────────

/// Build the command line for the chosen calculator.
/// python_exe should be the full path to the project venv interpreter.
inline std::string make_subprocess_cmd(
        const std::string&  python_exe,
        CalculatorTarget    target,
        const std::string&  input_json_b64)  ///< base-64 encoded JSON dataset
{
    const char* script = nullptr;
    switch (target) {
        case CalculatorTarget::QCalc:       script = "QCalc.py";             break;
        case CalculatorTarget::CP2:         script = "CondensedPhysics2.py"; break;
        case CalculatorTarget::CP1:         script = "CondensedPhysics.py";  break;
        case CalculatorTarget::UQFF_SERVER: script = "uqff_server.js";       break;
    }
    if (!script) return {};
    return python_exe + " " + script + " --ipc-json-b64=" + input_json_b64;
}

// ─────────────────────────────────────────────────────────────
// NAMED-PIPE MODE BRIDGE  (Windows)
// ─────────────────────────────────────────────────────────────
#ifdef _WIN32

/// Synchronous call: send one IPC message over the StarMagic_UQFF named pipe
/// and wait for a matching RESPONSE_* message.  Returns BridgeResult.
///
/// Requires uqff_ipc.py server to be running (started by production_pipeline.py
/// or manually via:  python uqff_ipc.py --serve).
inline BridgeResult call_via_pipe(
        IPC::MessageType          request_type,
        uint32_t                  session_id,
        const void*               payload,
        uint32_t                  payload_len,
        DWORD                     timeout_ms = 8000)
{
    BridgeResult result;

    HANDLE pipe = IPC::connect_client_pipe();
    if (pipe == INVALID_HANDLE_VALUE) {
        result.error_message = "Failed to connect to StarMagic_UQFF named pipe";
        return result;
    }

    // Build and send framed request
    std::vector<uint8_t> frame;
    static uint32_t seq = 0;
    IPC::build_message(request_type, session_id, ++seq, payload, payload_len, frame);

    DWORD written = 0;
    if (!WriteFile(pipe, frame.data(), static_cast<DWORD>(frame.size()), &written, nullptr)
        || written != static_cast<DWORD>(frame.size()))
    {
        result.error_message = "WriteFile to named pipe failed";
        CloseHandle(pipe);
        return result;
    }

    // Read response header
    std::vector<uint8_t> resp_buf(sizeof(IPC::MessageHeader) + 64 * 1024);
    DWORD read_bytes = 0;
    if (!ReadFile(pipe, resp_buf.data(), static_cast<DWORD>(resp_buf.size()), &read_bytes, nullptr)) {
        result.error_message = "ReadFile from named pipe failed";
        CloseHandle(pipe);
        return result;
    }
    CloseHandle(pipe);

    if (!IPC::validate_message(resp_buf.data(), read_bytes)) {
        result.error_message = "Invalid UQFF IPC framing in response";
        return result;
    }

    const auto* hdr  = reinterpret_cast<const IPC::MessageHeader*>(resp_buf.data());
    const char* body = reinterpret_cast<const char*>(resp_buf.data() + sizeof(IPC::MessageHeader));

    if (hdr->type == IPC::MessageType::RESPONSE_JSON  ||
        hdr->type == IPC::MessageType::CP2_RESPONSE)
    {
        result.json_payload  = std::string(body, hdr->payload_length);
        result.success       = true;
    } else if (hdr->type == IPC::MessageType::RESPONSE_EQUATION) {
        result.equation_text = std::string(body, hdr->payload_length);
        result.success       = true;
    } else if (static_cast<uint16_t>(hdr->type) >= 0x0500 &&
               static_cast<uint16_t>(hdr->type) <  0x0600)
    {
        // Error group
        if (hdr->payload_length >= sizeof(IPC::ErrorPayload)) {
            const auto* err = reinterpret_cast<const IPC::ErrorPayload*>(body);
            result.error_message = err->description;
        }
    }

    return result;
}

#endif // _WIN32

// ─────────────────────────────────────────────────────────────
// CALLBACK / ASYNC  (platform-agnostic)
// ─────────────────────────────────────────────────────────────
using BridgeCallback = std::function<void(BridgeResult)>;

/// Launch a Python calculator asynchronously in a detached thread.
/// result_cb is invoked on completion (from the worker thread — caller
/// must synchronise if updating UI).
void call_async(
        const std::string&   python_exe,
        CalculatorTarget     target,
        const std::string&   json_dataset,
        BridgeCallback       result_cb);

} // namespace Bridge
} // namespace UQFF

#endif // PYTHON_BRIDGE_H
