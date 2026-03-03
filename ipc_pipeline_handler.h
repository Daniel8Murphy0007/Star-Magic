/**
 * @file ipc_pipeline_handler.h
 * @brief IPC Pipeline Handler for source2(HEAD PROGRAM).cpp
 * 
 * Implements PIPELINE_PROCESS message handling:
 * 1. Receives PIPELINE_PROCESS message from source2.cpp
 * 2. Calls QCalc.solve() via Python subprocess (qcalc_subprocess.py)
 * 3. Returns RESPONSE_DATA with results
 * 
 * NOTE: We use QCalc (not CondensedPhysics) for Phase 0 because:
 * - QCalc: 9,149 lines, imports in 1.09s, subprocess ~920ms
 * - CondensedPhysics: 168,494 lines, imports in 30s+, causes timeouts
 * - Both provide same 8 UQFF Master Equations
 * - QCalc gives 30x speed improvement
 * 
 * Author: Daniel T. Murphy
 * Date: March 2, 2026
 * Phase: 0 - Unification (IPC Wiring)
 */

#ifndef IPC_PIPELINE_HANDLER_H
#define IPC_PIPELINE_HANDLER_H

#include <QProcess>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QString>
#include <QDebug>
#include <QElapsedTimer>
#include <memory>
#include <string>
#include <vector>

// Platform-specific Named Pipe includes
#ifdef _WIN32
#include <windows.h>
#include <string>
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#endif

namespace UQFF {

/**
 * @struct PipelineProcessRequest
 * @brief Request structure for PIPELINE_PROCESS message
 */
struct PipelineProcessRequest {
    char object_name[128];      // "SGR 1935+2154"
    uint32_t timeout_ms;        // Computation timeout (default 30000)
    char callback_id[64];       // For async callbacks
    
    // Optional parameters (0 = not set)
    double M;                   // Mass [solar masses]
    double r;                   // Radius [meters]
    double z;                   // Redshift
    double B;                   // Magnetic field [Tesla]
    double T;                   // Temperature [K]
    double SFR;                 // Star formation rate
};

/**
 * @struct PipelineResultPayload
 * @brief Response structure for PIPELINE_RESULT message
 */
struct PipelineResultPayload {
    double F_U;                 // Unified field [N/kg]
    double Ug1, Ug2, Ug3, Ug4;  // Gravity components [N/kg]
    double Um;                  // Magnetism [N/kg]
    double Ubi;                 // Buoyancy [N/kg]
    uint32_t calculators_run;   // Number of calculators executed
    uint32_t calculation_success; // Number successful
    double compute_time_ms;     // Total computation time
    uint32_t json_payload_follows; // 1 if full JSON attached
};

/**
 * @class IPCPipelineHandler
 * @brief Handles PIPELINE_PROCESS messages and Python subprocess communication
 */
class IPCPipelineHandler {
public:
    IPCPipelineHandler() : 
        python_path_("python"),
        script_path_("qcalc_subprocess.py"),  // Fast: ~920ms (vs CondensedPhysics 30s+)
        default_timeout_ms_(5000)  // Reduced from 30s to 5s (QCalc is fast!)
    {
        qDebug() << "[IPC Pipeline] Handler initialized (using QCalc)";
    }
    
    ~IPCPipelineHandler() {
        qDebug() << "[IPC Pipeline] Handler destroyed";
    }
    
    /**
     * @brief Process PIPELINE_PROCESS request
     * @param request Deserialized request payload
     * @return JSON result from QCalc (8 UQFF Master Equations)
     */
    QJsonObject processPipelineRequest(const PipelineProcessRequest& request) {
        QElapsedTimer timer;
        timer.start();
        
        qDebug() << "[IPC Pipeline] Processing request for:" << request.object_name;
        
        // Build input JSON for Python subprocess
        QJsonObject input;
        input["object_name"] = QString::fromUtf8(request.object_name);
        
        // Add optional parameters if set
        if (request.M > 0) input["M"] = request.M;
        if (request.r > 0) input["r"] = request.r;
        if (request.z != 0) input["z"] = request.z;
        if (request.B > 0) input["B"] = request.B;
        if (request.T > 0) input["T"] = request.T;
        if (request.SFR > 0) input["SFR"] = request.SFR;
        
        // Determine timeout
        uint32_t timeout = (request.timeout_ms > 0) ? request.timeout_ms : default_timeout_ms_;
        
        // Call Python subprocess
        QJsonObject result = callCondensedPhysics(input, timeout);
        
        // Add metadata
        result["backend_time_ms"] = timer.elapsed();
        result["callback_id"] = QString::fromUtf8(request.callback_id);
        
        qDebug() << "[IPC Pipeline] Request completed in" << timer.elapsed() << "ms";
        
        return result;
    }
    
    /**
     * @brief Set Python interpreter path (optional, default "python")
     */
    void setPythonPath(const QString& path) {
        python_path_ = path;
        qDebug() << "[IPC Pipeline] Python path set to:" << python_path_;
    }
    
    /**
     * @brief Set script path (optional, default "qcalc_subprocess.py")
     */
    void setScriptPath(const QString& path) {
        script_path_ = path;
        qDebug() << "[IPC Pipeline] Script path set to:" << script_path_;
    }
    
private:
    QString python_path_;
    QString script_path_;
    uint32_t default_timeout_ms_;
    
    /**
     * @brief Call Python subprocess with JSON input, return JSON output
     */
    QJsonObject callCondensedPhysics(const QJsonObject& input, uint32_t timeout_ms) {
        QProcess process;
        
        // Set up process
        process.setProgram(python_path_);
        process.setArguments(QStringList() << script_path_);
        
        qDebug() << "[IPC Pipeline] Starting subprocess:" << python_path_ << script_path_;
        
        // Start process
        process.start();
        
        if (!process.waitForStarted(5000)) {
            qWarning() << "[IPC Pipeline] Failed to start Python subprocess:" << process.errorString();
            return createErrorResponse("Failed to start Python subprocess: " + process.errorString());
        }
        
        // Write JSON input to stdin
        QJsonDocument inputDoc(input);
        QByteArray inputData = inputDoc.toJson(QJsonDocument::Compact) + "\n";
        
        qDebug() << "[IPC Pipeline] Sending input:" << inputData.size() << "bytes";
        
        process.write(inputData);
        process.closeWriteChannel(); // Signal EOF to Python
        
        // Wait for completion
        if (!process.waitForFinished(timeout_ms)) {
            qWarning() << "[IPC Pipeline] Subprocess timeout after" << timeout_ms << "ms";
            process.kill();
            process.waitForFinished(1000);
            return createErrorResponse("Subprocess timeout after " + QString::number(timeout_ms) + "ms");
        }
        
        // Read output
        QByteArray stdoutData = process.readAllStandardOutput();
        QByteArray stderrData = process.readAllStandardError();
        
        qDebug() << "[IPC Pipeline] Subprocess finished. Exit code:" << process.exitCode();
        qDebug() << "[IPC Pipeline] stdout size:" << stdoutData.size() << "bytes";
        
        if (!stderrData.isEmpty()) {
            qDebug() << "[IPC Pipeline] stderr:" << QString::fromUtf8(stderrData);
        }
        
        // Parse JSON output
        QJsonDocument outputDoc = QJsonDocument::fromJson(stdoutData);
        
        if (outputDoc.isNull()) {
            qWarning() << "[IPC Pipeline] Failed to parse JSON output";
            qWarning() << "[IPC Pipeline] Raw output:" << QString::fromUtf8(stdoutData);
            return createErrorResponse("Failed to parse subprocess JSON output");
        }
        
        if (!outputDoc.isObject()) {
            qWarning() << "[IPC Pipeline] Output is not a JSON object";
            return createErrorResponse("Subprocess returned invalid JSON (not an object)");
        }
        
        return outputDoc.object();
    }
    
    /**
     * @brief Create error response JSON
     */
    QJsonObject createErrorResponse(const QString& error) {
        QJsonObject obj;
        obj["success"] = false;
        obj["error"] = error;
        obj["compute_time_ms"] = 0;
        return obj;
    }
};

/**
 * @class NamedPipeServer
 * @brief Simple Named Pipe server for Windows IPC
 */
class NamedPipeServer {
public:
    NamedPipeServer(const std::string& pipe_name) : 
        pipe_name_(pipe_name),
        running_(false)
    {
#ifdef _WIN32
        full_pipe_name_ = "\\\\.\\pipe\\" + pipe_name_;
#else
        full_pipe_name_ = "/tmp/" + pipe_name_;
#endif
    }
    
    ~NamedPipeServer() {
        stop();
    }
    
    /**
     * @brief Start listening for connections (blocks)
     */
    bool start() {
        running_ = true;
        
        qDebug() << "[Named Pipe] Starting server:" << QString::fromStdString(full_pipe_name_);
        
#ifdef _WIN32
        // Windows Named Pipe implementation
        while (running_) {
            // Create named pipe
            HANDLE hPipe = CreateNamedPipeA(
                full_pipe_name_.c_str(),
                PIPE_ACCESS_DUPLEX,
                PIPE_TYPE_MESSAGE | PIPE_READMODE_MESSAGE | PIPE_WAIT,
                PIPE_UNLIMITED_INSTANCES,
                4096,  // Output buffer size
                4096,  // Input buffer size
                0,     // Default timeout
                NULL   // Default security
            );
            
            if (hPipe == INVALID_HANDLE_VALUE) {
                qWarning() << "[Named Pipe] CreateNamedPipe failed:" << GetLastError();
                return false
;
            }
            
            qDebug() << "[Named Pipe] Waiting for client connection...";
            
            // Wait for client
BOOL connected = ConnectNamedPipe(hPipe, NULL) ? TRUE : (GetLastError() == ERROR_PIPE_CONNECTED);
            
            if (!connected) {
                qWarning() << "[Named Pipe] ConnectNamedPipe failed:" << GetLastError();
                CloseHandle(hPipe);
                continue;
            }
            
            qDebug() << "[Named Pipe] Client connected";
            
            // Handle connection (read message, process, respond)
            handleConnection(hPipe);
            
            // Disconnect and close
            DisconnectNamedPipe(hPipe);
            CloseHandle(hPipe);
        }
#else
        // POSIX Named Pipe (FIFO) implementation
        // TODO: Implement for Linux/macOS
        qWarning() << "[Named Pipe] POSIX implementation not yet available";
        return false;
#endif
        
        return true;
    }
    
    /**
     * @brief Stop server
     */
    void stop() {
        running_ = false;
        qDebug() << "[Named Pipe] Server stopped";
    }
    
    /**
     * @brief Set message handler callback
     */
    void setMessageHandler(std::function<QJsonObject(const QJsonObject&)> handler) {
        message_handler_ = handler;
    }
    
private:
    std::string pipe_name_;
    std::string full_pipe_name_;
    bool running_;
    std::function<QJsonObject(const QJsonObject&)> message_handler_;
    
#ifdef _WIN32
    void handleConnection(HANDLE hPipe) {
        const size_t BUFFER_SIZE = 4096;
        char buffer[BUFFER_SIZE];
        DWORD bytesRead = 0;
        
        // Read message
        BOOL success = ReadFile(hPipe, buffer, BUFFER_SIZE - 1, &bytesRead, NULL);
        
        if (!success || bytesRead == 0) {
            qWarning() << "[Named Pipe] ReadFile failed or no data";
            return;
        }
        
        buffer[bytesRead] = '\0';
        
        qDebug() << "[Named Pipe] Received" << bytesRead << "bytes";
        
        // Parse JSON request
        QJsonDocument requestDoc = QJsonDocument::fromJson(QByteArray(buffer, bytesRead));
        
        if (requestDoc.isNull() || !requestDoc.isObject()) {
            qWarning() << "[Named Pipe] Invalid JSON request";
            return;
        }
        
        QJsonObject request = requestDoc.object();
        
        // Call handler
        QJsonObject response;
        if (message_handler_) {
            response = message_handler_(request);
        } else {
            response["error"] = "No message handler configured";
        }
        
        // Send response
        QJsonDocument responseDoc(response);
        QByteArray responseData = responseDoc.toJson(QJsonDocument::Compact);
        DWORD bytesWritten = 0;
        
        success = WriteFile(hPipe, responseData.data(), responseData.size(), &bytesWritten, NULL);
        
        if (!success) {
            qWarning() << "[Named Pipe] WriteFile failed";
        } else {
            qDebug() << "[Named Pipe] Sent response:" << bytesWritten << "bytes";
        }
    }
#endif
};

} // namespace UQFF

#endif // IPC_PIPELINE_HANDLER_H
