# Phase 0: IPC Unification - Integration Guide
## Wiring Instructions for source2.cpp ↔ source2(HEAD PROGRAM).cpp ↔ CondensedPhysics.py

**Date:** March 2, 2026  
**Status:** Implementation Ready  
**Estimated Time:** 2-3 hours for full integration

---

## OVERVIEW

This guide provides step-by-step instructions to wire the 3 critical IPC connections:

1. **Backend IPC Server** (source2(HEAD PROGRAM).cpp) - Listens for requests
2. **Frontend IPC Client** (source2.cpp) - Sends calculation requests  
3. **Python Bridge** (condensed_physics_subprocess.py) - Executes physics calculations

---

## FILES CREATED

| File | Purpose | Status |
|------|---------|--------|
| `qcalc_subprocess.py` | **PRODUCTION** Python subprocess wrapper (fast: ~920ms) | ✓ Tested |
| `ipc_pipeline_handler.h` | C++ IPC handler header (Qt-based) | ✓ Created |
| `test_ipc_pipeline.py` | Test suite to verify wiring | ⚠ Use manual test instead |
| `PHASE0_UNIFICATION_INTEGRATION_GUIDE.md` | This file | ✓ Created |
| `PHASE0_PROGRESS_REPORT.md` | Progress summary and decision log | ✓ Created |

**Note:** We're using `qcalc_subprocess.py` (not `condensed_physics_subprocess.py`) because:
- **QCalc.py:** 9,149 lines, imports in 1.09s, subprocess completes in ~920ms
- **CondensedPhysics.py:** 168,494 lines, imports in 30s+, causes timeouts
- Both provide same 8 UQFF Master Equations via UnifiedFieldSolver
- QCalc provides 30x speed improvement for Phase 0

---

## STEP 1: Test the Python Subprocess (2 minutes)

Before integrating, verify the Python subprocess works:

```powershell
# Activate virtual environment
& .\.venv_py314_backup\Scripts\Activate.ps1

# Test QCalc subprocess (should complete in ~1 second)
echo '{"object_name":"SGR 1745+29","M":2.0,"r":1000000.0,"B":1e13}' | python qcalc_subprocess.py 2>$null | python -m json.tool
```

**Expected Output (abbreviated):**
```json
{
    "query_id": "SGR 1745+29_20260302...",
    "timestamp": "2026-03-02T...",
    "long_form_equations": [
        {"name": "Ug1", "result": 2.00229e-22, "unit": "m/s²"},
        {"name": "Ug2", "result": 1.5858e-22, "unit": "m/s²"},
        {"name": "Ug3", "result": 2.4027e-22, "unit": "m/s²"},
        {"name": "Ug4", "result": 1.3349e-22, "unit": "m/s²"},
        {"name": "Ug", "result": 7.3257e-22, "unit": "m/s²"}
    ],
    "solutions": {
        "Ug": 7.32571e-22,
        "Ub_i": 1.87347e-12
    },
    "total_time_ms": 919.5
}
```

**✓ Success Criteria:**
- JSON output is valid (python -m json.tool formats it)
- `total_time_ms` is ~900-1200ms (fast!)
- Equations computed: Ug1, Ug2, Ug3, Ug4, Ug, Ub_i
- No errors or exceptions

If test passes, proceed to Step 2. If it fails:
- Check virtual environment is activated
- Verify QCalc.py exists in current directory
- Run: `pip install -r requirements.txt` to install dependencies

---

## STEP 2: Add IPC Server to source2(HEAD PROGRAM).cpp (30-45 min)

### 2.1 Add Include at Top

After line 100 (after existing includes), add:

```cpp
// IPC Pipeline Handler (Phase 0 - Unification)
#include "ipc_pipeline_handler.h"
```

### 2.2 Add Global IPC Handler Instance

After line 200 (in global variables section), add:

```cpp
// ============================================================================
// IPC PIPELINE (Phase 0 - Unification)
// ============================================================================
std::unique_ptr<UQFF::IPCPipelineHandler> g_ipc_handler = nullptr;
std::unique_ptr<UQFF::NamedPipeServer> g_pipe_server = nullptr;
std::thread g_ipc_thread;
bool g_ipc_running = false;
```

### 2.3 Add IPC Initialization Function

After line 400 (after helper functions), add:

```cpp
// ============================================================================
// IPC PIPELINE INITIALIZATION (Phase 0)
// ============================================================================

/**
 * @brief Initialize IPC server for physics calculation requests
 */
void InitializeIPCServer() {
    qDebug() << "[IPC Server] Initializing...";
    
    // Create handler
    g_ipc_handler = std::make_unique<UQFF::IPCPipelineHandler>();
    
    // Optional: Set custom Python path (if not in PATH)
    // g_ipc_handler->setPythonPath("C:/Python312/python.exe");
    
    // Set script path (default: qcalc_subprocess.py in current directory)
    // g_ipc_handler->setScriptPath("qcalc_subprocess.py");
    
    // Create Named Pipe server
    g_pipe_server = std::make_unique<UQFF::NamedPipeServer>("StarMagic_UQFF");
    
    // Set message handler
    g_pipe_server->setMessage Handler([](const QJsonObject& request) -> QJsonObject {
        if (!g_ipc_handler) {
            QJsonObject error;
            error["success"] = false;
            error["error"] = "IPC handler not initialized";
            return error;
        }
        
        // Extract message type
        QString msg_type = request["type"].toString();
        
        if (msg_type == "PIPELINE_PROCESS") {
            // Convert JSON to PipelineProcessRequest
            UQFF::PipelineProcessRequest req = {};
            
            // Copy object name
            QString object_name = request["object_name"].toString();
            strncpy(req.object_name, object_name.toUtf8().constData(), 127);
            req.object_name[127] = '\0';
            
            // Optional parameters
            req.M = request["M"].toDouble(0);
            req.r = request["r"].toDouble(0);
            req.z = request["z"].toDouble(0);
            req.B = request["B"].toDouble(0);
            req.T = request["T"].toDouble(0);
            req.SFR = request["SFR"].toDouble(0);
            req.timeout_ms = request["timeout_ms"].toInt(30000);
            
            QString callback_id = request["callback_id"].toString();
            strncpy(req.callback_id, callback_id.toUtf8().constData(), 63);
            req.callback_id[63] = '\0';
            
            // Process request
            return g_ipc_handler->processPipelineRequest(req);
        } else {
            QJsonObject error;
            error["success"] = false;
            error["error"] = "Unknown message type: " + msg_type;
            return error;
        }
    });
    
    // Start server in background thread
    g_ipc_running = true;
    g_ipc_thread = std::thread([]() {
        qDebug() << "[IPC Server] Starting server thread...";
        g_pipe_server->start();  // Blocks until stop()
        qDebug() << "[IPC Server] Server thread stopped";
    });
    
    qDebug() << "[IPC Server] Initialization complete";
}

/**
 * @brief Shutdown IPC server
 */
void ShutdownIPCServer() {
    qDebug() << "[IPC Server] Shutting down...";
    
    if (g_pipe_server) {
        g_pipe_server->stop();
    }
    
    if (g_ipc_thread.joinable()) {
        g_ipc_thread.join();
    }
    
    g_ipc_handler.reset();
    g_pipe_server.reset();
    g_ipc_running = false;
    
    qDebug() << "[IPC Server] Shutdown complete";
}
```

### 2.4 Call Initialization in main()

Find the `main()` function (around line 3400), and add **AFTER** `QApplication app(argc, argv);` but **BEFORE** `app.exec()`:

```cpp
int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    
    // ... existing setup code ...
    
    // Initialize IPC Server (Phase 0 - Unification)
    InitializeIPCServer();
    
    // ... existing code ...
    
    int result = app.exec();
    
    // Shutdown IPC Server before exit
    ShutdownIPCServer();
    
    return result;
}
```

---

## STEP 3: Add IPC Client to source2.cpp (30-45 min)

### 3.1 Add IPC Client Helper Class

After line 500 (after helper functions), add:

```cpp
// ============================================================================
// IPC CLIENT (Phase 0 - Unification)
// ============================================================================

/**
 * @class IPCClient
 * @brief Simple IPC client to send requests to backend
 */
class IPCClient {
public:
    IPCClient(const QString& pipeName = "StarMagic_UQFF") : pipe_name_(pipeName) {}
    
    /**
     * @brief Send PIPELINE_PROCESS request to backend
     */
    QJsonObject sendPipelineRequest(const QString& objectName, const QJsonObject& params) {
        QJsonObject request;
        request["type"] = "PIPELINE_PROCESS";
        request["object_name"] = objectName;
        request["callback_id"] = QUuid::createUuid().toString();
        request["timeout_ms"] = 30000;
        
        // Add optional parameters
        if (params.contains("M")) request["M"] = params["M"];
        if (params.contains("r")) request["r"] = params["r"];
        if (params.contains("z")) request["z"] = params["z"];
        if (params.contains("B")) request["B"] = params["B"];
        if (params.contains("T")) request["T"] = params["T"];
        if (params.contains("SFR")) request["SFR"] = params["SFR"];
        
        qDebug() << "[IPC Client] Sending request for:" << objectName;
        
#ifdef _WIN32
        // Windows Named Pipe client
        QString fullPipeName = QString("\\\\.\\pipe\\") + pipe_name_;
        
        HANDLE hPipe = CreateFileW(
            (LPCWSTR)fullPipeName.utf16(),
            GENERIC_READ | GENERIC_WRITE,
            0,
            NULL,
            OPEN_EXISTING,
            0,
            NULL
        );
        
        if (hPipe == INVALID_HANDLE_VALUE) {
            qWarning() << "[IPC Client] Failed to connect to pipe:" << GetLastError();
            QJsonObject error;
            error["success"] = false;
            error["error"] = "Failed to connect to backend (pipe not available)";
            return error;
        }
        
        // Send request
        QJsonDocument requestDoc(request);
        QByteArray requestData = requestDoc.toJson(QJsonDocument::Compact);
        DWORD bytesWritten = 0;
        
        BOOL success = WriteFile(hPipe, requestData.data(), requestData.size(), &bytesWritten, NULL);
        
        if (!success) {
            qWarning() << "[IPC Client] WriteFile failed:" << GetLastError();
            CloseHandle(hPipe);
            QJsonObject error;
            error["success"] = false;
            error["error"] = "Failed to write to pipe";
            return error;
        }
        
        qDebug() << "[IPC Client] Sent" << bytesWritten << "bytes";
        
        // Read response
        const size_t BUFFER_SIZE = 65536;  // 64KB buffer
        char buffer[BUFFER_SIZE];
        DWORD bytesRead = 0;
        
        success = ReadFile(hPipe, buffer, BUFFER_SIZE - 1, &bytesRead, NULL);
        
        CloseHandle(hPipe);
        
        if (!success) {
            qWarning() << "[IPC Client] ReadFile failed:" << GetLastError();
            QJsonObject error;
            error["success"] = false;
            error["error"] = "Failed to read from pipe";
            return error;
        }
        
        buffer[bytesRead] = '\0';
        
        qDebug() << "[IPC Client] Received" << bytesRead << "bytes";
        
        // Parse response
        QJsonDocument responseDoc = QJsonDocument::fromJson(QByteArray(buffer, bytesRead));
        
        if (responseDoc.isNull() || !responseDoc.isObject()) {
            qWarning() << "[IPC Client] Invalid JSON response";
            QJsonObject error;
            error["success"] = false;
            error["error"] = "Invalid JSON response from backend";
            return error;
        }
        
        return responseDoc.object();
#else
        // POSIX implementation (TODO)
        QJsonObject error;
        error["success"] = false;
        error["error"] = "IPC not yet implemented for POSIX";
        return error;
#endif
    }
    
private:
    QString pipe_name_;
};
```

### 3.2 Add IPC Call in PerformSearch()

Find the `PerformSearch()` function (around line 2500). After **APIFetch.py is called and dataset is retrieved**, add:

```cpp
// Inside PerformSearch() function, after APIFetch completes:

// PHASE 0 UNIFICATION: Send to backend for physics calculation
IPCClient ipc_client;

QJsonObject params;
// TODO: Extract parameters from APIFetch result (bodies_*.csv)
// For now, use query string as object name
params["object_name"] = QString::fromStdString(query);

// Send request to backend
QJsonObject result = ipc_client.sendPipelineRequest(
    QString::fromStdString(query), 
    params
);

// Check result
if (result["success"].toBool()) {
    qDebug() << "✓ Backend calculation successful";
    qDebug() << "  compute_time_ms:" << result["compute_time_ms"].toDouble();
    
    // TODO: Display results in Tab 8 (VTK) and Tab 9 (Session Logger)
    // For now, just log
    if (result.contains("solutions")) {
        qDebug() << "  Solutions computed";
    }
} else {
    qWarning() << "✗ Backend calculation failed:" << result["error"].toString();
}
```

---

## STEP 4: Test End-to-End (15 minutes)

### 4.1 Build source2(HEAD PROGRAM).cpp

```powershell
# Using CMake
cmake --build build_msvc --config Release --target source2

# Or use Visual Studio solution
```

### 4.2 Run Backend First

```powershell
# Terminal 1: Start backend server
.\build_msvc\Release\source2(HEAD PROGRAM).exe
```

**Expected Output:**
```
[IPC Server] Initializing...
[IPC Server] Starting server thread...
[Named Pipe] Starting server: \\.\pipe\StarMagic_UQFF
[Named Pipe] Waiting for client connection...
```

### 4.3 Run Frontend

```powershell
# Terminal 2: Start frontend GUI
.\build_msvc\Release\source2.exe
```

### 4.4 Test Query

1. In source2.exe GUI, go to **Tab 1** (Query field)
2. Enter: `SGR 1745+29`
3. Press **Enter**
4. Watch **Terminal 1** (backend) for IPC messages
5. Watch **Terminal 2** (frontend) for results

**Expected Backend Output:**
```
[Named Pipe] Client connected
[Named Pipe] Received 245 bytes
[IPC Pipeline] Processing request for: SGR 1745+29
[IPC Pipeline] Starting subprocess: python condensed_physics_subprocess.py
[IPC Pipeline] Sending input: 87 bytes
[IPC Pipeline] Subprocess finished. Exit code: 0
[IPC Pipeline] stdout size: 1523 bytes
[IPC Pipeline] Request completed in 347ms
[Named Pipe] Sent response: 1523 bytes
```

**Expected Frontend Output:**
```
[IPC Client] Sending request for: SGR 1745+29
[IPC Client] Sent 245 bytes
[IPC Client] Received 1523 bytes
✓ Backend calculation successful
  compute_time_ms: 347
  Solutions computed
```

---

## STEP 5: Verify Data Flow (5 minutes)

Check the complete pipeline:

```
USER QUERY "SGR 1745+29"
    ↓
source2.cpp Tab 1 → APIFetch.py → bodies_20260302_*.csv
    ↓
IPCClient.sendPipelineRequest() → Named Pipe \\.\pipe\StarMagic_UQFF
    ↓
source2(HEAD PROGRAM).cpp receives PIPELINE_PROCESS
    ↓
IPCPipelineHandler.processPipelineRequest()
    ↓
QProcess → python condensed_physics_subprocess.py
    ↓
CondensedPhysics.solve({"object_name": "SGR 1745+29", ...})
    ↓
JSON result → stdout → IPCPipelineHandler
    ↓
Named Pipe response → source2.cpp
    ↓
Display in Tab 8 (VTK) + Tab 9 (Session Logger)
```

---

## TROUBLESHOOTING

### Issue: "Failed to connect to backend (pipe not available)"

**Cause:** Backend not running or Named Pipe not created  
**Solution:** 
1. Verify backend (source2(HEAD PROGRAM).exe) is running
2. Check backend console shows "Waiting for client connection..."
3. On Windows, verify pipe exists: `Get-ChildItem \\.\pipe\`

### Issue: "Failed to start Python subprocess"

**Cause:** Python not in PATH or script not found  
**Solution:**
1. Verify Python path: `python --version`
2. Check script exists: `Test-Path condensed_physics_subprocess.py`
3. In IPC handler, set explicit path:
   ```cpp
   g_ipc_handler->setPythonPath("C:/Python312/python.exe");
   ```

### Issue: "Subprocess timeout after 30000ms"

**Cause:** CondensedPhysics taking too long or hanging  
**Solution:**
1. Test subprocess manually: `python test_ipc_pipeline.py`
2. Increase timeout in `sendPipelineRequest()`:
   ```cpp
   request["timeout_ms"] = 60000;  // 60 seconds
   ```

### Issue: "CondensedPhysics.py not available"

**Cause:** Python imports failing  
**Solution:**
1. Activate venv: `& .\.venv_py314_backup\Scripts\Activate.ps1`
2. Install dependencies: `pip install -r requirements.txt`
3. Verify imports: `python -c "from CondensedPhysics import UnifiedFieldSolver"`

---

## NEXT STEPS (Post-Unification)

After Phase 0 is complete and working:

1. **Phase 1: Code Quality** (3-4 hours)
   - Deduplicate QCalc vs CondensedPhysics (remove duplicate UnifiedFieldSolver)
   - Implement Session Logger (Tab 9) recall functionality

2. **Phase 2: Extensions** (1-2 hours)
   - Wire CondensedPhysics2.py to parent pipeline
   - Optional: Expose QCalc_API.py Flask REST as alternative to subprocess

3. **Phase 3: Polish** (2-3 hours)
   - Add error handling, retry logic
   - Implement progress bars for long calculations
   - Add caching for repeated queries

---

## COMPLETION CHECKLIST

- [ ] Test subprocess: `python test_ipc_pipeline.py` → 3/3 tests pass
- [ ] Backend includes `ipc_pipeline_handler.h`
- [ ] Backend calls `InitializeIPCServer()` in main()
- [ ] Backend shows "Waiting for client connection..." on startup
- [ ] Frontend includes `IPCClient` class
- [ ] Frontend calls `ipc_client.sendPipelineRequest()` in PerformSearch()
- [ ] End-to-end test: Query in GUI → Backend processes → Results returned
- [ ] Backend console shows IPC traffic
- [ ] Frontend console shows calculation results

When all items checked: **Phase 0 Unification Complete** ✓

---

**Estimated Total Time:** 2-3 hours for full integration + testing  
**Next Phase:** Code Quality (Deduplication + Session Logger)
