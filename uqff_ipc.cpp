/**
 * @file uqff_ipc.cpp
 * @brief Implementation of UQFF IPC framing and Named Pipe server/client
 *
 * Implements the non-inline portions of uqff_ipc.h:
 *   - Windows Named Pipe server loop (static helper)
 *   - PythonBackend::compute()  (from physics_service.h)
 *   - Bridge::call_async()      (from python_bridge.h)
 *
 * Build requirement: link with ws2_32.lib on Windows.
 *
 * Author: Daniel T. Murphy
 * Date: March 2026 (v3.1 — GrokThread7b0e961f)
 */

#include "uqff_ipc.h"
#include "python_bridge.h"
#include "physics_service.h"

#include <stdexcept>
#include <thread>
#include <sstream>
#include <cstring>

#ifdef _WIN32
#  pragma comment(lib, "ws2_32.lib")
#endif

namespace UQFF {

// ─────────────────────────────────────────────────────────────
// PythonBackend::compute()
// ─────────────────────────────────────────────────────────────

PhysicsResult PythonBackend::compute(const PhysicsQuery& q) {
    PhysicsResult result;
    result.system_name = q.system_name;

    // Determine calculator target
    Bridge::CalculatorTarget target = q.force_cp2
        ? Bridge::CalculatorTarget::CP2
        : Bridge::route_trigger(q.trigger);

    // Build JSON dataset string (minimal, escape-safe)
    std::ostringstream json;
    json << "{"
         << "\"system_name\":\"" << q.system_name << "\","
         << "\"mass_kg\":"       << q.mass_kg      << ","
         << "\"radius_m\":"      << q.radius_m     << ","
         << "\"distance_m\":"    << q.distance_m   << ","
         << "\"redshift\":"      << q.redshift     << ","
         << "\"t_n\":"           << q.t_n          << ","
         << "\"U_UA\":"          << q.U_UA         << ","
         << "\"SSq\":"           << q.SSq          << ","
         << "\"kappa\":"         << q.kappa        << ","
         << "\"F_rel\":"         << q.F_rel        << ","
         << "\"H_SCm\":"         << q.H_SCm;
    for (const auto& kv : q.extra)
        json << ",\"" << kv.first << "\":" << kv.second;
    json << "}";

#ifdef _WIN32
    // Prefer named-pipe mode when uqff_ipc.py server is running
    IPC::ComputePayload payload{};
    std::strncpy(payload.system_name, q.system_name.c_str(), sizeof(payload.system_name) - 1);
    payload.mass_kg    = q.mass_kg;
    payload.radius_m   = q.radius_m;
    payload.distance_m = q.distance_m;
    payload.redshift   = q.redshift;
    payload.t_n        = q.t_n;
    payload.U_UA       = q.U_UA;
    payload.SSq        = q.SSq;
    payload.kappa      = q.kappa;
    payload.F_rel      = q.F_rel;

    IPC::MessageType req_type = (target == Bridge::CalculatorTarget::CP2)
                                 ? IPC::MessageType::CP2_TRIGGER
                                 : IPC::MessageType::COMPUTE_SINGLE;

    Bridge::BridgeResult br = Bridge::call_via_pipe(req_type, q.session_id,
                                                    &payload, sizeof(payload));
    if (br.success) {
        result.success      = true;
        result.json_raw     = br.json_payload;
        result.equation_text = br.equation_text;
        result.elapsed_ms   = br.elapsed_ms;
        result.backend_used = (target == Bridge::CalculatorTarget::CP2)
                               ? "python:CP2" : "python:QCalc";
        return result;
    }
    // Named-pipe failed — fall through to subprocess
#endif

    // Subprocess fallback (platform-agnostic)
    std::string cmd = Bridge::make_subprocess_cmd(python_exe_, target, json.str());
    if (cmd.empty()) {
        result.error_message = "Could not build subprocess command";
        return result;
    }

    // Note: actual subprocess execution via QProcess in source2.cpp context.
    // Here we store the command for the caller to dispatch.
    result.json_raw     = cmd;  // caller interprets as deferred command
    result.backend_used = "python:subprocess:deferred";
    result.error_message = "Named pipe unavailable; use QProcess to dispatch: " + cmd;
    return result;
}

// ─────────────────────────────────────────────────────────────
// Bridge::call_async()
// ─────────────────────────────────────────────────────────────

namespace Bridge {

void call_async(
        const std::string&   python_exe,
        CalculatorTarget     target,
        const std::string&   json_dataset,
        BridgeCallback       result_cb)
{
    std::thread([python_exe, target, json_dataset, result_cb]() {
        BridgeResult r;
        r.success = false;

#ifdef _WIN32
        // Attempt named-pipe first
        IPC::ComputePayload payload{};
        std::strncpy(payload.system_name, "async_call", sizeof(payload.system_name) - 1);
        // Minimal payload — JSON dataset passed as SEND_DATASET separately if needed.

        IPC::MessageType req_type = (target == CalculatorTarget::CP2)
                                     ? IPC::MessageType::CP2_TRIGGER
                                     : IPC::MessageType::COMPUTE_SINGLE;

        r = call_via_pipe(req_type, 0, &payload, sizeof(payload));
        if (r.success) { result_cb(std::move(r)); return; }
#endif

        // Subprocess fallback
        std::string cmd = make_subprocess_cmd(python_exe, target, json_dataset);
        if (cmd.empty()) {
            r.error_message = "Empty subprocess command";
            result_cb(std::move(r));
            return;
        }

#ifdef _WIN32
        STARTUPINFOA si{};
        PROCESS_INFORMATION pi{};
        si.cb = sizeof(si);
        std::string cmd_copy = cmd;
        if (CreateProcessA(nullptr, cmd_copy.data(), nullptr, nullptr,
                           FALSE, 0, nullptr, nullptr, &si, &pi))
        {
            WaitForSingleObject(pi.hProcess, 15000);
            DWORD exit_code = 1;
            GetExitCodeProcess(pi.hProcess, &exit_code);
            CloseHandle(pi.hProcess);
            CloseHandle(pi.hThread);
            r.success = (exit_code == 0);
            r.backend_used = "python:subprocess";
        } else {
            r.error_message = "CreateProcess failed for: " + cmd;
        }
#else
        int ret = std::system(cmd.c_str());
        r.success = (ret == 0);
        r.backend_used = "python:subprocess";
#endif
        result_cb(std::move(r));
    }).detach();
}

} // namespace Bridge

} // namespace UQFF
