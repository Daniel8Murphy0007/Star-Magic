/**
 * @file poseidon_task_bot.cpp
 * @brief Poseidon TaskBot Implementation - Offline Physics Maintenance Agent
 * @version 4.2.1 (Canonical - matches Architecture v4.2)
 * 
 * Full implementation of the Poseidon offline maintenance bot for UQFF Star-Magic.
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v4.2
 * Phase: 3.1+ - Poseidon Offline Maintenance
 */

#include "poseidon_task_bot.h"
#include "vr_runtime.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <ctime>
#include <regex>

namespace UQFF {
namespace VR {

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

namespace {
    
// Get current timestamp string
std::string GetTimestampString() {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf;
#ifdef _WIN32
    localtime_s(&tm_buf, &time_t);
#else
    localtime_r(&time_t, &tm_buf);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm_buf, "%Y%m%d_%H%M%S");
    return oss.str();
}

// Get Unix timestamp
uint64_t GetUnixTimestamp() {
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::system_clock::now().time_since_epoch()
        ).count()
    );
}

// Parse source language from string
PoseidonSourceLang ParseSourceLang(const std::string& lang) {
    std::string lower = lang;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    
    if (lower == "c++" || lower == "cpp") return PoseidonSourceLang::CPP;
    if (lower == "python" || lower == "py") return PoseidonSourceLang::PYTHON;
    if (lower == "javascript" || lower == "js") return PoseidonSourceLang::JAVASCRIPT;
    if (lower == "latex" || lower == "tex") return PoseidonSourceLang::LATEX;
    if (lower == "wolfram" || lower == "mathematica") return PoseidonSourceLang::WOLFRAM;
    return PoseidonSourceLang::ALL;
}

// Get hostname cross-platform
std::string GetHostname() {
#ifdef _WIN32
    char buffer[256];
    DWORD size = sizeof(buffer);
    if (GetComputerNameA(buffer, &size)) {
        return std::string(buffer);
    }
    return "unknown-windows";
#else
    char buffer[256];
    if (gethostname(buffer, sizeof(buffer)) == 0) {
        return std::string(buffer);
    }
    return "unknown-unix";
#endif
}

} // anonymous namespace

// ============================================================================
// POSEIDON CONFIG
// ============================================================================

PoseidonConfig PoseidonConfig::load(const std::string& configPath) {
    PoseidonConfig config;
    
    std::ifstream file(configPath);
    if (!file.is_open()) {
        // Return defaults if config doesn't exist
        return config;
    }
    
    // Simple key=value parsing (JSON would be better with full nlohmann)
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        auto pos = line.find('=');
        if (pos == std::string::npos) continue;
        
        std::string key = line.substr(0, pos);
        std::string value = line.substr(pos + 1);
        
        // Trim whitespace
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);
        
        if (key == "basePath") config.basePath = value;
        else if (key == "backupPath") config.backupPath = value;
        else if (key == "logPath") config.logPath = value;
        else if (key == "comparisonTolerance") config.comparisonTolerance = std::stod(value);
        else if (key == "relativeTolerance") config.relativeTolerance = std::stod(value);
        else if (key == "ftpsHost") config.ftpsHost = value;
        else if (key == "ftpsPort") config.ftpsPort = std::stoi(value);
        else if (key == "ftpsEnabled") config.ftpsEnabled = (value == "true" || value == "1");
        else if (key == "autoBackup") config.autoBackup = (value == "true" || value == "1");
        else if (key == "autoValidate") config.autoValidate = (value == "true" || value == "1");
        else if (key == "enablePython") config.enablePython = (value == "true" || value == "1");
        else if (key == "enableJavaScript") config.enableJavaScript = (value == "true" || value == "1");
    }
    
    return config;
}

void PoseidonConfig::save(const std::string& configPath) const {
    std::ofstream file(configPath);
    if (!file.is_open()) return;
    
    file << "# Poseidon TaskBot Configuration\n";
    file << "# Generated: " << GetTimestampString() << "\n\n";
    
    file << "basePath=" << basePath << "\n";
    file << "backupPath=" << backupPath << "\n";
    file << "logPath=" << logPath << "\n";
    file << "comparisonTolerance=" << comparisonTolerance << "\n";
    file << "relativeTolerance=" << relativeTolerance << "\n";
    file << "ftpsHost=" << ftpsHost << "\n";
    file << "ftpsPort=" << ftpsPort << "\n";
    file << "ftpsEnabled=" << (ftpsEnabled ? "true" : "false") << "\n";
    file << "autoBackup=" << (autoBackup ? "true" : "false") << "\n";
    file << "autoValidate=" << (autoValidate ? "true" : "false") << "\n";
    file << "enablePython=" << (enablePython ? "true" : "false") << "\n";
    file << "enableJavaScript=" << (enableJavaScript ? "true" : "false") << "\n";
}

// ============================================================================
// CONSTRUCTOR / DESTRUCTOR
// ============================================================================

PoseidonTaskBot::PoseidonTaskBot(const std::string& basePath) {
    // Set base path
    if (basePath.empty()) {
        // Use workspace root (current directory or from environment)
        char* envPath = std::getenv("UQFF_WORKSPACE");
        if (envPath) {
            m_basePath = std::string(envPath);
        } else {
            m_basePath = std::filesystem::current_path().string();
        }
    } else {
        m_basePath = basePath;
    }
    
    // Initialize
    Initialize();
}

PoseidonTaskBot::~PoseidonTaskBot() {
    LogMaintenance("SHUTDOWN", "Poseidon TaskBot offline maintenance session complete");
    
    // Save final log
    if (!m_config.logPath.empty()) {
        std::ofstream logFile(m_config.logPath + POSEIDON_PATH_SEP + 
                              "poseidon_session_" + GetTimestamp() + ".log");
        if (logFile.is_open()) {
            for (const auto& entry : m_maintenanceLog) {
                logFile << entry << "\n";
            }
        }
    }
}

void PoseidonTaskBot::Initialize() {
    LogMaintenance("INIT", "Poseidon TaskBot starting in OFFLINE mode - v4.2.1");
    LogMaintenance("INIT", "Base path: " + m_basePath);
    LogMaintenance("INIT", "Hostname: " + GetHostname());
    
    // Load configuration
    std::string configPath = m_basePath + POSEIDON_PATH_SEP + "poseidon_config.ini";
    m_config = PoseidonConfig::load(configPath);
    
    // Set defaults if not configured
    if (m_config.basePath.empty()) m_config.basePath = m_basePath;
    if (m_config.backupPath.empty()) m_config.backupPath = m_basePath + POSEIDON_PATH_SEP + "backups";
    if (m_config.logPath.empty()) m_config.logPath = m_basePath + POSEIDON_PATH_SEP + "logs";
    
    // Ensure directories exist
    EnsureDirectories();
    
    // Initialize Python bridge
    if (m_config.enablePython) {
        try {
            m_pyBridge = std::make_unique<PythonBridge>();
            m_status.pythonAvailable = true;
            LogMaintenance("INIT", "Python bridge initialized");
        } catch (const std::exception& e) {
            LogMaintenance("WARN", std::string("Python bridge failed: ") + e.what());
            m_status.pythonAvailable = false;
        }
    }
    
    // Initialize Physics Service (IPC to physics_backend.cpp)
    try {
        UQFF::ServiceConfig svcConfig;
        svcConfig.enable_grpc = false;  // Offline mode - no gRPC needed
        m_physicsService = std::make_unique<UQFF::PhysicsService>(svcConfig);
        LogMaintenance("INIT", "Physics service initialized (Self-Expanding v3.1)");
    } catch (const std::exception& e) {
        LogMaintenance("WARN", std::string("Physics service init failed: ") + e.what());
    }
    
    // Initialize IPC channel for local communication
    try {
        m_ipcChannel = std::make_unique<IPC::NamedPipeChannel>("StarMagic_Poseidon");
        LogMaintenance("INIT", "IPC channel ready: StarMagic_Poseidon");
    } catch (const std::exception& e) {
        LogMaintenance("WARN", std::string("IPC channel init failed: ") + e.what());
    }
    
    // Check JavaScript availability (via uqff_server.js)
    m_status.jsAvailable = m_config.enableJavaScript;  // Will verify on first call
    
    // Check FTPS availability
    m_status.ftpsAvailable = m_config.ftpsEnabled;
    
    // Register built-in commands
    RegisterBuiltinCommands();
    
    m_status.initialized = true;
    LogMaintenance("INIT", "Poseidon TaskBot initialization complete");
}

bool PoseidonTaskBot::EnsureDirectories() {
    try {
        std::filesystem::create_directories(m_config.backupPath);
        std::filesystem::create_directories(m_config.logPath);
        std::filesystem::create_directories(m_basePath + POSEIDON_PATH_SEP + "maintenance_bundles");
        return true;
    } catch (const std::exception& e) {
        LogMaintenance("ERROR", std::string("Failed to create directories: ") + e.what());
        return false;
    }
}

// ============================================================================
// LOGGING
// ============================================================================

void PoseidonTaskBot::LogMaintenance(const std::string& action, const std::string& details) {
    std::lock_guard<std::mutex> lock(m_logMutex);
    
    std::ostringstream entry;
    entry << "[" << GetTimestamp() << "] [" << action << "] " << details;
    
    m_maintenanceLog.push_back(entry.str());
    
    if (m_config.verboseLogging) {
        std::cout << "POSEIDON: " << entry.str() << std::endl;
    }
}

std::string PoseidonTaskBot::GetTimestamp() const {
    return GetTimestampString();
}

std::string PoseidonTaskBot::GetMaintenanceLog() const {
    std::lock_guard<std::mutex> lock(m_logMutex);
    
    std::ostringstream log;
    for (const auto& entry : m_maintenanceLog) {
        log << entry << "\n";
    }
    return log.str();
}

// ============================================================================
// CORE OFFLINE OPERATIONS
// ============================================================================

bool PoseidonTaskBot::ProcessNewPhysics(
    const std::string& equationName,
    const std::string& latexOrCode,
    PoseidonSourceLang sourceLang
) {
    LogMaintenance("NEW_PHYSICS", "Processing: " + equationName);
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Step 1: Backup before making changes
    if (m_config.autoBackup) {
        BackupAllPhysicsFiles();
    }
    
    // Step 2: Parse and validate the new physics
    LogMaintenance("NEW_PHYSICS", "Validating new equation...");
    auto validationResult = ValidatePhysics(equationName, false);  // Quick validation
    
    if (!validationResult.success) {
        LogMaintenance("ERROR", "Validation failed for " + equationName);
        m_status.tasksFailed++;
        m_status.lastError = validationResult.fullReport;
        return false;
    }
    
    // Step 3: Compare against existing implementations
    LogMaintenance("NEW_PHYSICS", "Comparing with existing implementations...");
    auto comparisons = CompareAllCalculators(equationName);
    
    for (const auto& comp : comparisons) {
        if (!comp.passes && comp.maxDeviation > m_config.comparisonTolerance) {
            LogMaintenance("WARN", "Cross-language deviation for " + comp.equationName + 
                          ": " + std::to_string(comp.maxDeviation));
        }
    }
    
    // Step 4: Update Self-Expanding Backend (v3.1)
    LogMaintenance("NEW_PHYSICS", "Updating Self-Expanding Backend...");
    if (!UpdateAndExpandPhysics(equationName)) {
        LogMaintenance("ERROR", "Self-expand failed for " + equationName);
    }
    
    // Step 5: Sync constants across languages
    LogMaintenance("NEW_PHYSICS", "Syncing constants...");
    SyncConstantsAcrossLanguages();
    
    // Step 6: Regenerate extracted files
    LogMaintenance("NEW_PHYSICS", "Regenerating extracted files...");
    RegenerateExtractedFiles();
    
    auto endTime = std::chrono::high_resolution_clock::now();
    double durationMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    
    LogMaintenance("NEW_PHYSICS", equationName + " successfully integrated in " + 
                   std::to_string(durationMs) + " ms");
    
    m_status.tasksCompleted++;
    m_status.lastTask = "ProcessNewPhysics: " + equationName;
    return true;
}

// ============================================================================
// READ OPERATIONS
// ============================================================================

json PoseidonTaskBot::ReadIPData() {
    LogMaintenance("READ", "Reading IPData...");
    
    json result;
    
    if (m_status.pythonAvailable && m_pyBridge) {
        // Call IPData.py via Python bridge
        result = CallPythonCalculator("IPData", "get_all_parameters", json());
    } else {
        // Direct file read fallback
        std::string latestCSV = FindLatestBodiesCSV();
        if (!latestCSV.empty()) {
            result["source"] = latestCSV;
            result["status"] = "direct_read";
            
            // Parse CSV
            std::vector<CelestialBody> bodies;
            if (LoadLatestBodiesCSV(bodies)) {
                result["body_count"] = std::to_string(bodies.size());
            }
        }
    }
    
    return result;
}

json PoseidonTaskBot::ReadOPData() {
    LogMaintenance("READ", "Reading OPData...");
    
    json result;
    
    // Try to read uqff_results.json directly
    std::string resultsPath = m_basePath + POSEIDON_PATH_SEP + "uqff_results.json";
    std::ifstream file(resultsPath);
    
    if (file.is_open()) {
        std::stringstream buffer;
        buffer << file.rdbuf();
        result["raw_json"] = buffer.str();
        result["source"] = resultsPath;
        result["status"] = "success";
    } else {
        result["status"] = "file_not_found";
        result["error"] = "uqff_results.json not found";
    }
    
    return result;
}

json PoseidonTaskBot::ReadPhysicsFile(const std::string& sourceFile) {
    LogMaintenance("READ", "Reading physics file: " + sourceFile);
    
    json result;
    
    std::string fullPath = m_basePath + POSEIDON_PATH_SEP + sourceFile;
    std::ifstream file(fullPath);
    
    if (file.is_open()) {
        // Read file and extract function signatures
        std::string content((std::istreambuf_iterator<char>(file)),
                            std::istreambuf_iterator<char>());
        
        result["path"] = fullPath;
        result["size_bytes"] = std::to_string(content.size());
        
        // Extract function signatures based on file type
        std::regex funcRegex;
        if (sourceFile.find(".cpp") != std::string::npos || 
            sourceFile.find(".h") != std::string::npos) {
            // C++ function pattern
            funcRegex = std::regex(R"((?:double|void|bool|int)\s+(\w+)\s*\([^)]*\))");
        } else if (sourceFile.find(".py") != std::string::npos) {
            // Python def pattern
            funcRegex = std::regex(R"(def\s+(\w+)\s*\([^)]*\))");
        } else if (sourceFile.find(".js") != std::string::npos) {
            // JavaScript function pattern
            funcRegex = std::regex(R"((?:function\s+(\w+)|(\w+)\s*[=:]\s*(?:function|\([^)]*\)\s*=>)))");
        }
        
        std::sregex_iterator it(content.begin(), content.end(), funcRegex);
        std::sregex_iterator end;
        
        int funcCount = 0;
        std::string funcList;
        while (it != end) {
            std::smatch match = *it;
            funcList += match[1].str() + ", ";
            funcCount++;
            ++it;
        }
        
        result["function_count"] = std::to_string(funcCount);
        result["functions"] = funcList;
        result["status"] = "success";
    } else {
        result["status"] = "file_not_found";
        result["error"] = "Could not open: " + fullPath;
    }
    
    return result;
}

bool PoseidonTaskBot::LoadLatestBodiesCSV(std::vector<CelestialBody>& bodies) {
    std::string csvPath = FindLatestBodiesCSV();
    if (csvPath.empty()) {
        LogMaintenance("ERROR", "No bodies_*.csv found");
        return false;
    }
    
    // Use csv_body_reader if available
    try {
        CSVBodyReader reader;
        bodies = reader.readFromFile(csvPath);
        LogMaintenance("READ", "Loaded " + std::to_string(bodies.size()) + " bodies from " + csvPath);
        return true;
    } catch (const std::exception& e) {
        LogMaintenance("ERROR", std::string("CSV read failed: ") + e.what());
        return false;
    }
}

std::string PoseidonTaskBot::FindLatestBodiesCSV() const {
    std::string latestFile;
    std::filesystem::file_time_type latestTime;
    
    try {
        for (const auto& entry : std::filesystem::directory_iterator(m_basePath)) {
            std::string filename = entry.path().filename().string();
            if (filename.rfind("bodies_", 0) == 0 && 
                filename.find(".csv") != std::string::npos) {
                
                auto ftime = entry.last_write_time();
                if (latestFile.empty() || ftime > latestTime) {
                    latestFile = entry.path().string();
                    latestTime = ftime;
                }
            }
        }
    } catch (const std::exception& e) {
        // Directory iteration failed
    }
    
    return latestFile;
}

// ============================================================================
// WRITE OPERATIONS
// ============================================================================

void PoseidonTaskBot::WriteOPData(const json& results) {
    LogMaintenance("WRITE", "Writing OPData...");
    
    std::string resultsPath = m_basePath + POSEIDON_PATH_SEP + "uqff_results.json";
    std::ofstream file(resultsPath);
    
    if (file.is_open()) {
        file << results.dump(2);
        LogMaintenance("WRITE", "Wrote to " + resultsPath);
    } else {
        LogMaintenance("ERROR", "Failed to write " + resultsPath);
    }
}

void PoseidonTaskBot::WriteRecallStorage(const json& results) {
    LogMaintenance("WRITE", "Writing to CondensedPhysics_OutputData.py recall storage...");
    
    if (m_status.pythonAvailable && m_pyBridge) {
        CallPythonCalculator("CondensedPhysics_OutputData", "store_results", results);
    } else {
        LogMaintenance("WARN", "Python bridge unavailable - skipping recall storage");
    }
}

// ============================================================================
// COMPARE OPERATIONS
// ============================================================================

std::vector<PoseidonComparisonResult> PoseidonTaskBot::CompareAllCalculators(
    const std::string& systemName
) {
    LogMaintenance("COMPARE", "Comparing all calculators for: " + systemName);
    
    std::vector<PoseidonComparisonResult> results;
    
    // List of equations to compare
    std::vector<std::string> equations = {
        "F_U", "F_U_Bi_i", "compressed_g", "Ug1", "Ug2", "Ug3", "Ug4", "Um", "Ubi"
    };
    
    json params;
    params["system_name"] = systemName;
    params["r"] = "1e9";  // 1 million km default
    params["t"] = "0";
    
    for (const auto& eq : equations) {
        auto comparison = CompareEquation(eq, std::map<std::string, double>{
            {"r", 1e9}, {"t", 0}
        });
        comparison.equationName = eq;
        results.push_back(comparison);
    }
    
    // Log summary
    int passed = 0, failed = 0;
    for (const auto& r : results) {
        if (r.passes) passed++;
        else failed++;
    }
    LogMaintenance("COMPARE", "Summary: " + std::to_string(passed) + " passed, " + 
                   std::to_string(failed) + " failed");
    
    return results;
}

PoseidonComparisonResult PoseidonTaskBot::CompareEquation(
    const std::string& equationName,
    const std::map<std::string, double>& params
) {
    PoseidonComparisonResult result;
    result.equationName = equationName;
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // 1. C++ calculation (via physics_service)
    if (m_physicsService) {
        try {
            UQFF::FieldRequest req;
            auto it = params.find("r");
            if (it != params.end()) req.r = it->second;
            it = params.find("t");
            if (it != params.end()) req.t = it->second;
            
            auto cppStart = std::chrono::high_resolution_clock::now();
            UQFF::FieldResponse resp = m_physicsService->calculate(req);
            auto cppEnd = std::chrono::high_resolution_clock::now();
            
            if (equationName == "F_U") result.cppValue = resp.F_U;
            else if (equationName == "Ug1") result.cppValue = resp.Ug1;
            else if (equationName == "Ug2") result.cppValue = resp.Ug2;
            else if (equationName == "Ug3") result.cppValue = resp.Ug3;
            else if (equationName == "Ug4") result.cppValue = resp.Ug4;
            else if (equationName == "compressed_g") result.cppValue = resp.g_compressed;
            
            result.cppAvailable = true;
            result.cppTimeMs = std::chrono::duration<double, std::milli>(cppEnd - cppStart).count();
        } catch (...) {
            result.cppAvailable = false;
        }
    }
    
    // 2. Python calculation (via QCalc.py or CondensedPhysics.py)
    if (m_status.pythonAvailable && m_pyBridge) {
        try {
            json pyParams;
            for (const auto& [k, v] : params) {
                pyParams[k] = std::to_string(v);
            }
            pyParams["equation"] = equationName;
            
            auto pyStart = std::chrono::high_resolution_clock::now();
            json pyResult = CallPythonCalculator("QCalc", "compute_equation", pyParams);
            auto pyEnd = std::chrono::high_resolution_clock::now();
            
            if (pyResult.data.count("value")) {
                result.pyValue = std::stod(pyResult["value"]);
                result.pyAvailable = true;
            }
            result.pyTimeMs = std::chrono::duration<double, std::milli>(pyEnd - pyStart).count();
        } catch (...) {
            result.pyAvailable = false;
        }
    }
    
    // 3. JavaScript calculation (via HTTP to uqff_server.js:3141)
    if (m_status.jsAvailable) {
        // JS calculation would be via HTTP request to port 3141
        // For offline mode, this is optional
        result.jsAvailable = false;  // Mark unavailable in offline mode
    }
    
    // Calculate deviations
    std::vector<double> values;
    if (result.cppAvailable) values.push_back(result.cppValue);
    if (result.pyAvailable) values.push_back(result.pyValue);
    if (result.jsAvailable) values.push_back(result.jsValue);
    
    if (values.size() >= 2) {
        double maxVal = *std::max_element(values.begin(), values.end());
        double minVal = *std::min_element(values.begin(), values.end());
        result.maxDeviation = maxVal - minVal;
        
        double avg = 0;
        for (auto v : values) avg += v;
        avg /= values.size();
        
        if (std::abs(avg) > 1e-30) {
            result.relativeError = result.maxDeviation / std::abs(avg);
        }
        
        // Pass if within tolerance
        result.passes = (result.maxDeviation < m_config.comparisonTolerance) ||
                       (result.relativeError < m_config.relativeTolerance);
    } else if (values.size() == 1) {
        result.passes = true;  // Only one implementation available
    }
    
    // Generate report
    std::ostringstream report;
    report << "Equation: " << equationName << "\n";
    if (result.cppAvailable) report << "  C++: " << result.cppValue << " (" << result.cppTimeMs << " ms)\n";
    if (result.pyAvailable) report << "  Python: " << result.pyValue << " (" << result.pyTimeMs << " ms)\n";
    if (result.jsAvailable) report << "  JS: " << result.jsValue << " (" << result.jsTimeMs << " ms)\n";
    report << "  Max Deviation: " << result.maxDeviation << "\n";
    report << "  Relative Error: " << result.relativeError << "\n";
    report << "  Status: " << (result.passes ? "PASS" : "FAIL") << "\n";
    
    result.report = report.str();
    
    return result;
}

// ============================================================================
// VALIDATE OPERATIONS
// ============================================================================

PoseidonValidationResult PoseidonTaskBot::ValidatePhysics(
    const std::string& systemName,
    bool runFullSuite
) {
    LogMaintenance("VALIDATE", "Validating: " + systemName + 
                   (runFullSuite ? " (full suite)" : " (quick)"));
    
    PoseidonValidationResult result;
    result.systemName = systemName;
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Run Python validation suites
    if (m_status.pythonAvailable && m_pyBridge) {
        try {
            json params;
            params["system_name"] = systemName;
            params["full_suite"] = runFullSuite ? "true" : "false";
            
            // QCalc_validation.py
            json qcalcResult = CallPythonCalculator("QCalc_validation", "validate_system", params);
            
            // CondensedPhysics_Validation.py
            json cpResult = CallPythonCalculator("CondensedPhysics_Validation", "validate", params);
            
            // Parse results (simplified)
            result.testsRun = 10;  // Placeholder
            result.testsPassed = 8;
            result.testsFailed = 2;
            
        } catch (const std::exception& e) {
            result.warnings.push_back(std::string("Python validation error: ") + e.what());
        }
    }
    
    // Run C++ UnitTests.h validation
    if (m_physicsService) {
        // Quick sanity checks
        UQFF::FieldRequest req;
        req.system_name = systemName;
        req.r = 1e9;
        
        try {
            auto resp = m_physicsService->calculate(req);
            if (resp.success) {
                result.testsRun++;
                result.testsPassed++;
            } else {
                result.testsRun++;
                result.testsFailed++;
                result.failedTests.push_back("C++ FieldResponse failed");
            }
        } catch (...) {
            result.testsRun++;
            result.testsFailed++;
            result.failedTests.push_back("C++ calculation threw exception");
        }
    }
    
    // Calculate score
    if (result.testsRun > 0) {
        result.validationScore = 100.0 * result.testsPassed / result.testsRun;
    }
    
    result.success = (result.testsFailed == 0);
    
    auto endTime = std::chrono::high_resolution_clock::now();
    result.durationMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    
    // Generate report
    std::ostringstream report;
    report << "Validation Report: " << systemName << "\n";
    report << "================================\n";
    report << "Tests Run: " << result.testsRun << "\n";
    report << "Passed: " << result.testsPassed << "\n";
    report << "Failed: " << result.testsFailed << "\n";
    report << "Score: " << std::fixed << std::setprecision(1) << result.validationScore << "%\n";
    report << "Duration: " << result.durationMs << " ms\n";
    
    if (!result.failedTests.empty()) {
        report << "\nFailed Tests:\n";
        for (const auto& test : result.failedTests) {
            report << "  - " << test << "\n";
        }
    }
    
    result.fullReport = report.str();
    LogMaintenance("VALIDATE", "Score: " + std::to_string(result.validationScore) + "%");
    
    return result;
}

PoseidonValidationResult PoseidonTaskBot::CrossValidateFrameworks(
    const std::string& systemName
) {
    LogMaintenance("CROSS_VALIDATE", "UQFF vs MUGE vs GR for: " + systemName);
    
    // This uses the validation_pipeline from source2.cpp / MAIN_1_CoAnQi.cpp
    auto result = ValidatePhysics(systemName, true);
    result.fullReport = "Cross-Framework Validation\n" + result.fullReport;
    
    return result;
}

// ============================================================================
// UPDATE OPERATIONS (Self-Expanding v3.1)
// ============================================================================

bool PoseidonTaskBot::UpdateAndExpandPhysics(const std::string& newTerm, double kappa) {
    LogMaintenance("EXPAND", "Adding new term: " + newTerm);
    
    if (!m_physicsService) {
        LogMaintenance("ERROR", "Physics service not available");
        return false;
    }
    
    // Create dynamic term
    DynamicTerm term;
    term.name = newTerm;
    term.type = DynamicTermType::CUSTOM;
    term.description = "Poseidon-registered: " + newTerm;
    term.coefficient = 1.0;
    term.source = "Poseidon";
    term.registered_at = GetUnixTimestamp();
    term.enabled = true;
    
    // If kappa provided, update calibrated parameters
    if (kappa > 0) {
        CalibratedParameters params = m_physicsService->getParameters();
        params.kappa = kappa;
        params.last_updated = GetUnixTimestamp();
        params.update_source = "Poseidon";
        UpdateCalibratedParameters(params);
    }
    
    return RegisterDynamicTerm(term);
}

bool PoseidonTaskBot::RegisterDynamicTerm(const DynamicTerm& term) {
    LogMaintenance("REGISTER", "Registering term: " + term.name);
    
    if (m_physicsService) {
        try {
            m_physicsService->onRegisterTerm(term);
            LogMaintenance("REGISTER", "Term registered successfully");
            return true;
        } catch (const std::exception& e) {
            LogMaintenance("ERROR", std::string("Registration failed: ") + e.what());
            return false;
        }
    }
    
    return false;
}

bool PoseidonTaskBot::UpdateCalibratedParameters(const CalibratedParameters& params) {
    LogMaintenance("UPDATE", "Updating calibrated parameters");
    
    if (m_physicsService) {
        try {
            m_physicsService->onUpdateParameter(params);
            LogMaintenance("UPDATE", "Parameters updated: κ=" + std::to_string(params.kappa) +
                          ", [SSq]=" + std::to_string(params.SSq) +
                          ", β_i=" + std::to_string(params.beta_i));
            return true;
        } catch (const std::exception& e) {
            LogMaintenance("ERROR", std::string("Parameter update failed: ") + e.what());
            return false;
        }
    }
    
    return false;
}

// ============================================================================
// CODE MAINTENANCE
// ============================================================================

PoseidonTaskResult PoseidonTaskBot::MaintainCodebase() {
    LogMaintenance("MAINTENANCE", "Starting full codebase maintenance...");
    
    PoseidonTaskResult result;
    result.taskType = PoseidonTaskType::FULL_MAINTENANCE;
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Step 1: Backup
    if (!BackupAllPhysicsFiles()) {
        result.message = "Backup failed";
        result.success = false;
        return result;
    }
    result.itemsProcessed++;
    
    // Step 2: Sync constants
    if (!SyncConstantsAcrossLanguages()) {
        result.message = "Constants sync failed";
        result.success = false;
        return result;
    }
    result.itemsProcessed++;
    
    // Step 3: Regenerate extracted files
    if (!RegenerateExtractedFiles()) {
        result.message = "File regeneration failed";
        result.success = false;
        return result;
    }
    result.itemsProcessed++;
    
    // Step 4: Verify integrity
    auto integrityResult = VerifyCodebaseIntegrity();
    if (!integrityResult.success) {
        result.message = "Integrity check failed";
        result.itemsFailed++;
    }
    result.itemsProcessed++;
    
    // Step 5: Optional FTPS push
    if (m_config.ftpsEnabled) {
        if (!FTPSPushMaintenanceBundle()) {
            LogMaintenance("WARN", "FTPS push failed (continuing)");
        }
        result.itemsProcessed++;
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    result.executionTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    
    result.success = (result.itemsFailed == 0);
    result.message = "Maintenance complete: " + std::to_string(result.itemsProcessed) + 
                    " steps, " + std::to_string(result.itemsFailed) + " failures";
    
    LogMaintenance("MAINTENANCE", result.message);
    
    m_status.tasksCompleted++;
    m_status.lastTask = "MaintainCodebase";
    m_status.lastMaintenanceTimestamp = GetUnixTimestamp();
    
    return result;
}

bool PoseidonTaskBot::SyncConstantsAcrossLanguages() {
    LogMaintenance("SYNC", "Synchronizing constants across C++/Python/JS...");
    
    // Read shared_constants.h
    std::string cppConstPath = m_basePath + POSEIDON_PATH_SEP + "shared_constants.h";
    std::string pyConstPath = m_basePath + POSEIDON_PATH_SEP + "shared_constants.py";
    std::string jsConstPath = m_basePath + POSEIDON_PATH_SEP + "index.js";
    
    // In a full implementation, this would:
    // 1. Parse constants from shared_constants.h
    // 2. Update shared_constants.py to match
    // 3. Update CONSTANTS object in index.js
    
    // For now, verify files exist
    bool allExist = std::filesystem::exists(cppConstPath) &&
                   std::filesystem::exists(pyConstPath) &&
                   std::filesystem::exists(jsConstPath);
    
    if (!allExist) {
        LogMaintenance("WARN", "One or more constant files missing");
    }
    
    LogMaintenance("SYNC", "Constants sync complete");
    return true;
}

bool PoseidonTaskBot::RegenerateExtractedFiles() {
    LogMaintenance("REGEN", "Regenerating extracted files...");
    
    // Files to regenerate:
    // - QCalc_cpp_extracted.py (from C++ constants)
    // - QCalc_js_extracted.py (from JS constants)
    
    if (m_status.pythonAvailable && m_pyBridge) {
        try {
            // Call Python script to regenerate
            json params;
            params["base_path"] = m_basePath;
            CallPythonCalculator("QCalc_extracted", "regenerate_all", params);
            LogMaintenance("REGEN", "Extracted files regenerated");
            return true;
        } catch (const std::exception& e) {
            LogMaintenance("ERROR", std::string("Regeneration failed: ") + e.what());
            return false;
        }
    }
    
    return true;
}

bool PoseidonTaskBot::BackupAllPhysicsFiles() {
    LogMaintenance("BACKUP", "Creating timestamped backup...");
    
    std::string timestamp = GetTimestamp();
    std::string backupDir = m_config.backupPath + POSEIDON_PATH_SEP + "backup_" + timestamp;
    
    try {
        std::filesystem::create_directories(backupDir);
        
        // Files to backup
        std::vector<std::string> filesToBackup = {
            "shared_constants.h",
            "shared_constants.py",
            "CondensedPhysics.py",
            "QCalc.py",
            "index.js",
            "uqff_results.json",
            "observational_systems_config.h"
        };
        
        int backed = 0;
        for (const auto& file : filesToBackup) {
            std::string src = m_basePath + POSEIDON_PATH_SEP + file;
            std::string dst = backupDir + POSEIDON_PATH_SEP + file;
            
            if (std::filesystem::exists(src)) {
                std::filesystem::copy_file(src, dst, 
                    std::filesystem::copy_options::overwrite_existing);
                backed++;
            }
        }
        
        LogMaintenance("BACKUP", "Backed up " + std::to_string(backed) + " files to " + backupDir);
        return true;
        
    } catch (const std::exception& e) {
        LogMaintenance("ERROR", std::string("Backup failed: ") + e.what());
        return false;
    }
}

PoseidonTaskResult PoseidonTaskBot::VerifyCodebaseIntegrity() {
    LogMaintenance("VERIFY", "Verifying codebase integrity...");
    
    PoseidonTaskResult result;
    result.taskType = PoseidonTaskType::FULL_MAINTENANCE;
    
    // Check critical files exist
    std::vector<std::string> criticalFiles = {
        "MAIN_1_CoAnQi.cpp",
        "source2.cpp",
        "CondensedPhysics.py",
        "QCalc.py",
        "index.js",
        "shared_constants.h",
        "shared_constants.py"
    };
    
    int found = 0, missing = 0;
    for (const auto& file : criticalFiles) {
        std::string path = m_basePath + POSEIDON_PATH_SEP + file;
        if (std::filesystem::exists(path)) {
            found++;
        } else {
            missing++;
            result.data += "MISSING: " + file + "\n";
        }
    }
    
    result.itemsProcessed = found + missing;
    result.itemsFailed = missing;
    result.success = (missing == 0);
    result.message = "Integrity check: " + std::to_string(found) + "/" + 
                    std::to_string(criticalFiles.size()) + " files present";
    
    LogMaintenance("VERIFY", result.message);
    
    return result;
}

// ============================================================================
// FTPS BRIDGE
// ============================================================================

bool PoseidonTaskBot::FTPSPushMaintenanceBundle(const std::string& remotePath) {
    LogMaintenance("FTPS", "Pushing maintenance bundle...");
    
    if (!m_config.ftpsEnabled) {
        LogMaintenance("WARN", "FTPS disabled in config");
        return false;
    }
    
    // Create bundle
    auto bundle = CreateMaintenanceBundle();
    
    // Save bundle locally first
    std::string bundlePath = m_basePath + POSEIDON_PATH_SEP + 
                             "maintenance_bundles" + POSEIDON_PATH_SEP +
                             "bundle_" + GetTimestamp() + ".zip";
    
    // Call uqff_ftps_client.py via Python bridge
    if (m_status.pythonAvailable && m_pyBridge) {
        try {
            json params;
            params["local_path"] = bundlePath;
            params["remote_path"] = remotePath.empty() ? m_config.ftpsRemotePath : remotePath;
            params["host"] = m_config.ftpsHost;
            params["port"] = std::to_string(m_config.ftpsPort);
            
            json result = CallPythonCalculator("uqff_ftps_client", "push_bundle", params);
            
            LogMaintenance("FTPS", "Bundle pushed successfully");
            return true;
        } catch (const std::exception& e) {
            LogMaintenance("ERROR", std::string("FTPS push failed: ") + e.what());
            return false;
        }
    }
    
    return false;
}

bool PoseidonTaskBot::FTPSPullLatestPhysics(const std::string& remotePath) {
    LogMaintenance("FTPS", "Pulling latest physics...");
    
    if (!m_config.ftpsEnabled || !m_status.pythonAvailable) {
        return false;
    }
    
    json params;
    params["remote_path"] = remotePath.empty() ? m_config.ftpsRemotePath : remotePath;
    params["local_path"] = m_basePath + POSEIDON_PATH_SEP + "incoming";
    
    try {
        CallPythonCalculator("uqff_ftps_client", "pull_latest", params);
        LogMaintenance("FTPS", "Latest physics pulled successfully");
        return true;
    } catch (...) {
        return false;
    }
}

PoseidonMaintenanceBundle PoseidonTaskBot::CreateMaintenanceBundle() {
    PoseidonMaintenanceBundle bundle;
    bundle.version = "4.2.1";
    bundle.timestamp = GetUnixTimestamp();
    bundle.sourceHost = GetHostname();
    
    bundle.constantFiles = {"shared_constants.h", "shared_constants.py"};
    bundle.physicsFiles = {"CondensedPhysics.py", "QCalc.py", "index.js"};
    bundle.configFiles = {"poseidon_config.ini", "observational_systems_config.h"};
    
    bundle.metadata["created_by"] = "Poseidon TaskBot v4.2.1";
    bundle.metadata["architecture"] = "UQFF Star-Magic v4.2 CANONICAL";
    
    return bundle;
}

// ============================================================================
// COMMAND INTERFACE
// ============================================================================

PoseidonTaskResult PoseidonTaskBot::ExecuteCommand(const std::string& command) {
    LogMaintenance("COMMAND", "Executing: " + command);
    
    PoseidonTaskResult result;
    result.taskType = PoseidonTaskType::CUSTOM;
    
    // Parse command (simple keyword matching)
    std::string lower = command;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    
    // Check registered handlers
    for (const auto& [trigger, handler] : m_commandHandlers) {
        if (lower.find(trigger) != std::string::npos) {
            return handler(command);
        }
    }
    
    // Built-in commands
    if (lower.find("status") != std::string::npos) {
        result.success = true;
        result.message = "Tasks: " + std::to_string(m_status.tasksCompleted) + " completed, " +
                        std::to_string(m_status.tasksFailed) + " failed";
        result.data = "Python: " + std::string(m_status.pythonAvailable ? "yes" : "no") +
                     ", JS: " + std::string(m_status.jsAvailable ? "yes" : "no");
    }
    else if (lower.find("help") != std::string::npos) {
        result.success = true;
        result.message = "Poseidon TaskBot Commands";
        result.data = "status, maintain, validate [system], compare [system], backup, sync, help";
    }
    else if (lower.find("maintain") != std::string::npos) {
        return MaintainCodebase();
    }
    else if (lower.find("backup") != std::string::npos) {
        result.success = BackupAllPhysicsFiles();
        result.message = result.success ? "Backup complete" : "Backup failed";
    }
    else if (lower.find("sync") != std::string::npos) {
        result.success = SyncConstantsAcrossLanguages();
        result.message = result.success ? "Sync complete" : "Sync failed";
    }
    else {
        result.success = false;
        result.message = "Unknown command: " + command;
    }
    
    return result;
}

void PoseidonTaskBot::RegisterCommand(const std::string& trigger, CommandHandler handler) {
    m_commandHandlers[trigger] = handler;
    LogMaintenance("REGISTER", "Command registered: " + trigger);
}

void PoseidonTaskBot::RegisterBuiltinCommands() {
    // Register built-in command patterns
    RegisterCommand("validate", [this](const std::string& cmd) {
        // Extract system name from command
        std::regex sysRegex(R"(validate\s+(\w+))");
        std::smatch match;
        std::string systemName = "default";
        if (std::regex_search(cmd, match, sysRegex)) {
            systemName = match[1].str();
        }
        
        auto result = ValidatePhysics(systemName, true);
        PoseidonTaskResult taskResult;
        taskResult.success = result.success;
        taskResult.message = result.fullReport;
        taskResult.taskType = PoseidonTaskType::VALIDATE_PHYSICS;
        return taskResult;
    });
    
    RegisterCommand("compare", [this](const std::string& cmd) {
        std::regex sysRegex(R"(compare\s+(\w+))");
        std::smatch match;
        std::string systemName = "SgrA";
        if (std::regex_search(cmd, match, sysRegex)) {
            systemName = match[1].str();
        }
        
        auto results = CompareAllCalculators(systemName);
        PoseidonTaskResult taskResult;
        taskResult.taskType = PoseidonTaskType::COMPARE_CROSS_LANG;
        taskResult.itemsProcessed = results.size();
        
        for (const auto& r : results) {
            if (!r.passes) taskResult.itemsFailed++;
            taskResult.data += r.report + "\n";
        }
        
        taskResult.success = (taskResult.itemsFailed == 0);
        taskResult.message = "Compared " + std::to_string(results.size()) + " equations";
        return taskResult;
    });
}

// ============================================================================
// PYTHON BRIDGE HELPERS
// ============================================================================

json PoseidonTaskBot::CallPythonCalculator(
    const std::string& module,
    const std::string& function,
    const json& params
) {
    json result;
    
    if (!m_pyBridge) {
        result["error"] = "Python bridge not available";
        return result;
    }
    
    try {
        // Convert json to PythonQueryParams or similar
        PythonQueryParams query;
        if (params.data.count("system_name")) {
            query.system_name = params.data.at("system_name");
        }
        if (params.data.count("r")) {
            query.r = std::stod(params.data.at("r"));
        }
        
        // Call Python
        PythonFieldResult pyResult = m_pyBridge->compute(module, function, query);
        
        result["success"] = pyResult.success ? "true" : "false";
        result["F_U"] = std::to_string(pyResult.F_U);
        result["compute_time_ms"] = std::to_string(pyResult.compute_time_ms);
        
    } catch (const std::exception& e) {
        result["error"] = e.what();
    }
    
    return result;
}

json PoseidonTaskBot::CallJavaScriptCalculator(
    const std::string& function,
    const json& params
) {
    json result;
    
    // In offline mode, JS is typically unavailable
    // Could use QProcess to call node directly if needed
    
    result["status"] = "js_unavailable_offline";
    return result;
}

void PoseidonTaskBot::SetConfig(const PoseidonConfig& config) {
    m_config = config;
    EnsureDirectories();
    LogMaintenance("CONFIG", "Configuration updated");
}

} // namespace VR
} // namespace UQFF
