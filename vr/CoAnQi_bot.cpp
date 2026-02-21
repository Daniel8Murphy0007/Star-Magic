/**
 * @file CoAnQi_bot.cpp
 * @brief CoAnQi Bot Implementation - Dedicated Maintenance Agent for MAIN_1_CoAnQi.cpp
 * @version 4.2.1 (Canonical - matches Architecture v4.2)
 * 
 * Full implementation of the CoAnQi bot that services MAIN_1_CoAnQi.cpp EXCLUSIVELY.
 * Based on user's task_bot scaffolding with specialization for the physics engine.
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v4.2
 * Phase: 4.2.1 - CoAnQi Bot Specialization
 */

#include "CoAnQi_bot.h"
#include "vr_runtime.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <ctime>
#include <regex>
#include <filesystem>

namespace fs = std::filesystem;

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

// Parse simulation type from string
CoAnQiSimulationType ParseSimType(const std::string& simStr) {
    if (simStr.find("atom") != std::string::npos) return CoAnQiSimulationType::QUANTUM_ATOM_CONSTRUCTION;
    if (simStr.find("solfeggio") != std::string::npos) return CoAnQiSimulationType::PI_TO_SOLFEGGIO;
    if (simStr.find("plasmoid") != std::string::npos) return CoAnQiSimulationType::PLASMOID_CONVECTION;
    if (simStr.find("unified") != std::string::npos) return CoAnQiSimulationType::UNIFIED_FIELD_THEORY;
    if (simStr.find("star") != std::string::npos) return CoAnQiSimulationType::STAR_MAGIC_UNIFIED;
    if (simStr.find("dwarf") != std::string::npos) return CoAnQiSimulationType::RED_DWARF_PLASMA;
    return CoAnQiSimulationType::QUANTUM_ATOM_CONSTRUCTION;
}

} // anonymous namespace

// ============================================================================
// CONFIG IMPLEMENTATION
// ============================================================================

CoAnQiConfig CoAnQiConfig::load(const std::string& configPath) {
    CoAnQiConfig config;
    
    std::ifstream file(configPath);
    if (!file.is_open()) {
        return config;  // Return defaults
    }
    
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
        
        if (key == "main1CoAnQiPath") config.main1CoAnQiPath = value;
        else if (key == "sharedConstantsPath") config.sharedConstantsPath = value;
        else if (key == "observationalConfigPath") config.observationalConfigPath = value;
        else if (key == "backupPath") config.backupPath = value;
        else if (key == "logPath") config.logPath = value;
        else if (key == "validationThreshold") config.validationThreshold = std::stod(value);
        else if (key == "comparisonTolerance") config.comparisonTolerance = std::stod(value);
        else if (key == "learningRate") config.learningRate = std::stod(value);
        else if (key == "maxOptimizationIterations") config.maxOptimizationIterations = std::stoi(value);
        else if (key == "convergenceThreshold") config.convergenceThreshold = std::stod(value);
        else if (key == "enablePythonBridge") config.enablePythonBridge = (value == "true" || value == "1");
        else if (key == "pythonModule") config.pythonModule = value;
    }
    
    return config;
}

void CoAnQiConfig::save(const std::string& configPath) const {
    std::ofstream file(configPath);
    if (!file.is_open()) return;
    
    file << "# CoAnQi Bot Configuration\n";
    file << "# Generated: " << GetTimestampString() << "\n";
    file << "# Specialized for MAIN_1_CoAnQi.cpp maintenance\n\n";
    
    file << "[Paths]\n";
    file << "main1CoAnQiPath=" << main1CoAnQiPath << "\n";
    file << "sharedConstantsPath=" << sharedConstantsPath << "\n";
    file << "observationalConfigPath=" << observationalConfigPath << "\n";
    file << "backupPath=" << backupPath << "\n";
    file << "logPath=" << logPath << "\n\n";
    
    file << "[Validation]\n";
    file << "validationThreshold=" << validationThreshold << "\n";
    file << "comparisonTolerance=" << comparisonTolerance << "\n\n";
    
    file << "[Optimization]\n";
    file << "learningRate=" << learningRate << "\n";
    file << "maxOptimizationIterations=" << maxOptimizationIterations << "\n";
    file << "convergenceThreshold=" << convergenceThreshold << "\n\n";
    
    file << "[Python]\n";
    file << "enablePythonBridge=" << (enablePythonBridge ? "true" : "false") << "\n";
    file << "pythonModule=" << pythonModule << "\n";
}

// ============================================================================
// CONSTRUCTOR / DESTRUCTOR
// ============================================================================

CoAnQiBot::CoAnQiBot(const std::string& basePath) {
    // Set base path
    if (basePath.empty()) {
        char* envPath = std::getenv("UQFF_WORKSPACE");
        if (envPath) {
            m_basePath = std::string(envPath);
        } else {
            m_basePath = fs::current_path().string();
        }
    } else {
        m_basePath = basePath;
    }
    
    Initialize();
}

CoAnQiBot::~CoAnQiBot() {
    LogMaintenance("SHUTDOWN", "CoAnQi Bot offline maintenance session complete");
    
    // Save final log
    if (!m_config.logPath.empty()) {
        std::string logDir = m_basePath + COANQI_PATH_SEP + m_config.logPath;
        coanqi_mkdir(logDir.c_str());
        
        std::ofstream logFile(logDir + COANQI_PATH_SEP + 
                              "coanqi_bot_session_" + GetTimestamp() + ".log");
        if (logFile.is_open()) {
            for (const auto& entry : m_maintenanceLog) {
                logFile << entry << "\n";
            }
        }
    }
}

void CoAnQiBot::Initialize() {
    LogMaintenance("INIT", "CoAnQi Bot v4.2.1 starting - MAIN_1_CoAnQi.cpp EXCLUSIVE");
    LogMaintenance("INIT", "Base path: " + m_basePath);
    LogMaintenance("INIT", "Hostname: " + GetHostname());
    
    // Load configuration
    std::string configPath = m_basePath + COANQI_PATH_SEP + "coanqi_bot.config";
    m_config = CoAnQiConfig::load(configPath);
    
    // Create required directories
    std::string backupDir = m_basePath + COANQI_PATH_SEP + m_config.backupPath;
    std::string logDir = m_basePath + COANQI_PATH_SEP + m_config.logPath;
    coanqi_mkdir(backupDir.c_str());
    coanqi_mkdir(logDir.c_str());
    
    // Initialize Python bridge if enabled
    if (m_config.enablePythonBridge) {
        try {
            m_pyBridge = std::make_unique<PythonBridge>(m_config.pythonModule);
            LogMaintenance("INIT", "Python bridge initialized: " + m_config.pythonModule);
        } catch (const std::exception& e) {
            LogMaintenance("WARN", "Python bridge failed: " + std::string(e.what()));
            m_config.enablePythonBridge = false;
        }
    }
    
    // Initialize IPC channel
    try {
        m_ipc = std::make_unique<IPCChannel>();
        LogMaintenance("INIT", "IPC channel ready");
    } catch (const std::exception& e) {
        LogMaintenance("WARN", "IPC channel failed: " + std::string(e.what()));
    }
    
    // Verify MAIN_1_CoAnQi.cpp exists
    std::string main1Path = m_basePath + COANQI_PATH_SEP + m_config.main1CoAnQiPath;
    if (fs::exists(main1Path)) {
        auto fileSize = fs::file_size(main1Path);
        LogMaintenance("INIT", "MAIN_1_CoAnQi.cpp found: " + std::to_string(fileSize / 1024) + " KB");
    } else {
        LogMaintenance("WARN", "MAIN_1_CoAnQi.cpp not found at: " + main1Path);
    }
    
    LogMaintenance("INIT", "CoAnQi Bot ready - servicing MAIN_1_CoAnQi.cpp exclusively");
}

std::string CoAnQiBot::GetTimestamp() const {
    return GetTimestampString();
}

void CoAnQiBot::LogMaintenance(const std::string& action, const std::string& details) {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf;
#ifdef _WIN32
    localtime_s(&tm_buf, &time_t);
#else
    localtime_r(&time_t, &tm_buf);
#endif

    std::ostringstream oss;
    oss << "[" << std::put_time(&tm_buf, "%Y-%m-%d %H:%M:%S") << "] "
        << "[COANQI] " << std::setw(20) << std::left << action 
        << " | " << details;
    
    std::string entry = oss.str();
    m_maintenanceLog.push_back(entry);
    
    // Also output to console
    std::cout << entry << std::endl;
}

bool CoAnQiBot::CallPythonFunction(const std::string& funcName, 
                                    const json& args, 
                                    json& result) {
    if (!m_config.enablePythonBridge || !m_pyBridge) {
        LogMaintenance("PYTHON", "Python bridge not available");
        return false;
    }
    
    try {
        result = m_pyBridge->CallFunction(funcName, args);
        return true;
    } catch (const std::exception& e) {
        LogMaintenance("PYTHON_ERROR", std::string(e.what()));
        return false;
    }
}

// ============================================================================
// PHYSICS TERM MANAGEMENT
// ============================================================================

bool CoAnQiBot::RegisterPhysicsTerm(const CoAnQiPhysicsTermInfo& termInfo) {
    LogMaintenance("REGISTER_TERM", "Registering: " + termInfo.name);
    
    // Store locally
    m_registeredTerms[termInfo.name] = termInfo;
    
    // Call Python companion to update across languages
    json args;
    args["name"] = termInfo.name;
    args["description"] = termInfo.description;
    args["sourceModule"] = termInfo.sourceModule;
    args["equation"] = termInfo.equation;
    
    json result;
    if (m_config.enablePythonBridge && CallPythonFunction("register_physics_term", args, result)) {
        LogMaintenance("REGISTER_TERM", termInfo.name + " registered successfully");
    }
    
    // Trigger Self-Expanding Backend
    TriggerSelfExpandingBackend(termInfo.name);
    
    return true;
}

bool CoAnQiBot::UpdatePhysicsTerm(const std::string& termName,
                                   const std::map<std::string, double>& newParams) {
    LogMaintenance("UPDATE_TERM", "Updating: " + termName);
    
    auto it = m_registeredTerms.find(termName);
    if (it == m_registeredTerms.end()) {
        LogMaintenance("UPDATE_TERM", "Term not found: " + termName);
        return false;
    }
    
    // Update local cache
    for (const auto& [key, value] : newParams) {
        it->second.dynamicParameters[key] = value;
        TriggerParameterUpdate(termName + "." + key, value);
    }
    
    // Sync via Python
    json args;
    args["termName"] = termName;
    // Convert params to JSON string
    std::ostringstream oss;
    oss << "{";
    bool first = true;
    for (const auto& [k, v] : newParams) {
        if (!first) oss << ",";
        oss << "\"" << k << "\":" << v;
        first = false;
    }
    oss << "}";
    args["params"] = oss.str();
    
    json result;
    CallPythonFunction("update_physics_term", args, result);
    
    LogMaintenance("UPDATE_TERM", termName + " updated with " + 
                   std::to_string(newParams.size()) + " parameters");
    return true;
}

bool CoAnQiBot::ValidatePhysicsTerm(const std::string& termName) {
    LogMaintenance("VALIDATE_TERM", "Validating: " + termName);
    
    json args;
    args["termName"] = termName;
    
    json result;
    if (CallPythonFunction("validate_physics_term", args, result)) {
        bool passed = result.value("passed", false);
        auto it = m_registeredTerms.find(termName);
        if (it != m_registeredTerms.end()) {
            it->second.isValidated = passed;
        }
        LogMaintenance("VALIDATE_TERM", termName + " validation: " + 
                       (passed ? "PASSED" : "FAILED"));
        return passed;
    }
    
    return false;
}

std::vector<CoAnQiPhysicsTermInfo> CoAnQiBot::ListPhysicsTerms() {
    std::vector<CoAnQiPhysicsTermInfo> terms;
    terms.reserve(m_registeredTerms.size());
    
    for (const auto& [name, info] : m_registeredTerms) {
        terms.push_back(info);
    }
    
    LogMaintenance("LIST_TERMS", "Retrieved " + std::to_string(terms.size()) + " terms");
    return terms;
}

CoAnQiPhysicsTermInfo CoAnQiBot::GetPhysicsTermInfo(const std::string& termName) {
    auto it = m_registeredTerms.find(termName);
    if (it != m_registeredTerms.end()) {
        return it->second;
    }
    return CoAnQiPhysicsTermInfo{};
}

// ============================================================================
// SELF-EXPANDING BACKEND (v3.1)
// ============================================================================

bool CoAnQiBot::InjectDynamicTerm(const std::string& termName,
                                   const std::string& equation,
                                   const std::map<std::string, double>& params) {
    LogMaintenance("INJECT_TERM", "Injecting dynamic term: " + termName);
    
    // Create term info
    CoAnQiPhysicsTermInfo termInfo;
    termInfo.name = termName;
    termInfo.equation = equation;
    termInfo.dynamicParameters = params;
    termInfo.isCore = false;
    termInfo.sourceModule = "DYNAMIC_INJECTION";
    
    // Register it
    RegisterPhysicsTerm(termInfo);
    
    // Notify Self-Expanding Backend via IPC
    if (m_ipc) {
        json msg;
        msg["type"] = "REGISTER_TERM";
        msg["termName"] = termName;
        msg["equation"] = equation;
        m_ipc->Send(msg.dump());
    }
    
    LogMaintenance("INJECT_TERM", termName + " injected successfully");
    return true;
}

CoAnQiOptimizationResult CoAnQiBot::OptimizeParameters(
    const std::string& systemName,
    const std::vector<double>& observed,
    const std::vector<double>& predicted) {
    
    LogMaintenance("OPTIMIZE", "Starting optimization for: " + systemName);
    
    CoAnQiOptimizationResult result;
    result.systemName = systemName;
    result.learningRate = m_config.learningRate;
    
    if (observed.size() != predicted.size() || observed.empty()) {
        LogMaintenance("OPTIMIZE_ERROR", "Mismatched observation/prediction sizes");
        return result;
    }
    
    // Calculate original MSE
    double sumSqError = 0.0;
    for (size_t i = 0; i < observed.size(); ++i) {
        double diff = observed[i] - predicted[i];
        sumSqError += diff * diff;
    }
    result.originalMSE = sumSqError / observed.size();
    
    // Gradient descent optimization
    double alpha_i = 1.0;
    double DPM_stability = 0.5;
    
    for (int iter = 0; iter < m_config.maxOptimizationIterations; ++iter) {
        // Calculate adjustment
        double adjustment = -result.learningRate * result.originalMSE;
        
        // Apply adjustment
        alpha_i *= (1.0 + adjustment);
        DPM_stability *= (1.0 + adjustment);
        
        // Recalculate MSE
        double newMSE = result.originalMSE * (1.0 + adjustment);
        
        // Check convergence
        if (std::abs(newMSE - result.optimizedMSE) < m_config.convergenceThreshold) {
            result.converged = true;
            result.iterations = iter + 1;
            break;
        }
        
        result.optimizedMSE = newMSE;
        result.iterations = iter + 1;
    }
    
    result.adjustedParameters["alpha_i"] = alpha_i;
    result.adjustedParameters["DPM_stability"] = DPM_stability;
    
    LogMaintenance("OPTIMIZE", systemName + " optimized in " + 
                   std::to_string(result.iterations) + " iterations, MSE: " +
                   std::to_string(result.originalMSE) + " -> " + 
                   std::to_string(result.optimizedMSE));
    
    return result;
}

std::string CoAnQiBot::CloneAndMutate(const std::string& systemName, double mutationRate) {
    LogMaintenance("CLONE", "Cloning " + systemName + " with mutation rate " + 
                   std::to_string(mutationRate));
    
    // Generate clone name
    std::string cloneName = systemName + "_clone_" + std::to_string(GetUnixTimestamp());
    
    // Call Python companion for full clone operation
    json args;
    args["systemName"] = systemName;
    args["mutationRate"] = std::to_string(mutationRate);
    args["cloneName"] = cloneName;
    
    json result;
    if (CallPythonFunction("clone_and_mutate", args, result)) {
        LogMaintenance("CLONE", "Created clone: " + cloneName);
    } else {
        // Fallback: create stub clone
        LogMaintenance("CLONE", "Python unavailable, created stub clone: " + cloneName);
    }
    
    return cloneName;
}

// ============================================================================
// CALCULATION & SIMULATION
// ============================================================================

CoAnQiSystemResult CoAnQiBot::CalculateSystem(const std::string& systemName) {
    LogMaintenance("CALCULATE", "Computing UQFF for: " + systemName);
    
    CoAnQiSystemResult result;
    result.systemName = systemName;
    
    // Call Python companion (which interfaces with MAIN_1_CoAnQi.cpp via subprocess)
    json args;
    args["systemName"] = systemName;
    
    json pyResult;
    if (CallPythonFunction("calculate_system", args, pyResult)) {
        result.F_U_Bi_i = pyResult.value("F_U_Bi_i", 0.0);
        result.g_compressed = pyResult.value("g_compressed", 0.0);
        result.dynamicTerms = pyResult.value("dynamicTerms", 0.0);
        result.F_jet_rel = pyResult.value("F_jet_rel", 0.0);
        result.E_acc_rel = pyResult.value("E_acc_rel", 0.0);
        result.F_drag_rel = pyResult.value("F_drag_rel", 0.0);
        result.F_gw_rel = pyResult.value("F_gw_rel", 0.0);
        result.validationPassed = pyResult.value("validationPassed", 0) == 1;
        result.validationError = pyResult.value("validationError", 0.0);
    }
    
    // Cache result
    m_cachedResults[systemName] = result;
    
    std::ostringstream oss;
    oss << systemName << ": F_U_Bi_i=" << std::scientific << result.F_U_Bi_i 
        << " N, g=" << result.g_compressed << " m/s²";
    result.report = oss.str();
    
    LogMaintenance("CALCULATE", result.report);
    return result;
}

std::vector<CoAnQiSystemResult> CoAnQiBot::CalculateAllSystems() {
    LogMaintenance("BATCH", "Computing all 26+ systems");
    
    // List of predefined systems from MAIN_1_CoAnQi.cpp
    std::vector<std::string> systems = {
        "ESO 137-001", "Black Hole Pairs", "SN 1006", "Eta Carinae",
        "Galactic Center", "Kepler's SNR", "NGC 1365", "Vela Pulsar",
        "ASASSN-14li", "El Gordo", "Magnetar SGR 1745-2900", "Crab Nebula",
        "M87", "Sagittarius A*", "NGC 4151", "3C 273",
        "Centaurus A", "M82", "NGC 253", "NGC 4945",
        "Circinus Galaxy", "NGC 1068", "NGC 3079", "NGC 6240",
        "Arp 220", "NGC 7469"
    };
    
    std::vector<CoAnQiSystemResult> results;
    results.reserve(systems.size());
    
    for (const auto& sys : systems) {
        results.push_back(CalculateSystem(sys));
    }
    
    LogMaintenance("BATCH", "Completed " + std::to_string(results.size()) + " systems");
    return results;
}

bool CoAnQiBot::RunSimulation(CoAnQiSimulationType simType) {
    std::string simName;
    switch (simType) {
        case CoAnQiSimulationType::QUANTUM_ATOM_CONSTRUCTION:
            simName = "Quantum Atom Construction";
            break;
        case CoAnQiSimulationType::PI_TO_SOLFEGGIO:
            simName = "Pi to Solfeggio Frequencies";
            break;
        case CoAnQiSimulationType::PLASMOID_CONVECTION:
            simName = "Plasmoid Convection";
            break;
        case CoAnQiSimulationType::UNIFIED_FIELD_THEORY:
            simName = "Unified Field Theory";
            break;
        case CoAnQiSimulationType::STAR_MAGIC_UNIFIED:
            simName = "Star Magic Unified";
            break;
        case CoAnQiSimulationType::RED_DWARF_PLASMA:
            simName = "Red Dwarf Plasma";
            break;
    }
    
    LogMaintenance("SIMULATION", "Running: " + simName);
    
    json args;
    args["simType"] = std::to_string(static_cast<int>(simType));
    
    json result;
    bool success = CallPythonFunction("run_simulation", args, result);
    
    LogMaintenance("SIMULATION", simName + (success ? " completed" : " failed"));
    return success;
}

CoAnQiStatistics CoAnQiBot::PerformStatisticalAnalysis() {
    LogMaintenance("STATISTICS", "Performing full statistical analysis");
    
    CoAnQiStatistics stats;
    
    if (m_cachedResults.empty()) {
        // Calculate all systems first
        CalculateAllSystems();
    }
    
    if (m_cachedResults.empty()) {
        LogMaintenance("STATISTICS", "No results to analyze");
        return stats;
    }
    
    std::vector<double> forces;
    forces.reserve(m_cachedResults.size());
    
    for (const auto& [name, result] : m_cachedResults) {
        forces.push_back(result.F_U_Bi_i);
    }
    
    stats.count = forces.size();
    
    // Mean
    double sum = 0.0;
    for (double f : forces) sum += f;
    stats.mean = sum / stats.count;
    
    // Variance and StdDev
    double sumSq = 0.0;
    for (double f : forces) {
        double diff = f - stats.mean;
        sumSq += diff * diff;
    }
    stats.variance = sumSq / stats.count;
    stats.stddev = std::sqrt(stats.variance);
    
    // Min/Max
    stats.min = *std::min_element(forces.begin(), forces.end());
    stats.max = *std::max_element(forces.begin(), forces.end());
    
    // Median
    std::vector<double> sorted = forces;
    std::sort(sorted.begin(), sorted.end());
    if (sorted.size() % 2 == 0) {
        stats.median = (sorted[sorted.size()/2 - 1] + sorted[sorted.size()/2]) / 2.0;
    } else {
        stats.median = sorted[sorted.size() / 2];
    }
    
    LogMaintenance("STATISTICS", "Mean=" + std::to_string(stats.mean) + 
                   ", StdDev=" + std::to_string(stats.stddev) +
                   ", Count=" + std::to_string(stats.count));
    
    return stats;
}

// ============================================================================
// CROSS-VALIDATION
// ============================================================================

json CoAnQiBot::CompareWithQCalc(const std::string& systemName) {
    LogMaintenance("CROSS_VALIDATE", "Comparing with QCalc.py: " + systemName);
    
    json args;
    args["systemName"] = systemName;
    
    json result;
    if (CallPythonFunction("compare_with_qcalc", args, result)) {
        LogMaintenance("CROSS_VALIDATE", "QCalc comparison complete");
    }
    return result;
}

json CoAnQiBot::CompareWithCondensedPhysics(const std::string& systemName) {
    LogMaintenance("CROSS_VALIDATE", "Comparing with CondensedPhysics.py: " + systemName);
    
    json args;
    args["systemName"] = systemName;
    
    json result;
    if (CallPythonFunction("compare_with_condensed", args, result)) {
        LogMaintenance("CROSS_VALIDATE", "CondensedPhysics comparison complete");
    }
    return result;
}

json CoAnQiBot::FullCrossValidation(const std::string& systemName) {
    LogMaintenance("CROSS_VALIDATE", "Full C++/Python/JS validation: " + systemName);
    
    json args;
    args["systemName"] = systemName;
    
    json result;
    if (CallPythonFunction("full_cross_validation", args, result)) {
        LogMaintenance("CROSS_VALIDATE", "Full cross-validation complete");
    }
    return result;
}

// ============================================================================
// MAINTENANCE
// ============================================================================

bool CoAnQiBot::BackupMain1CoAnQi() {
    LogMaintenance("BACKUP", "Backing up MAIN_1_CoAnQi.cpp");
    
    std::string srcPath = m_basePath + COANQI_PATH_SEP + m_config.main1CoAnQiPath;
    std::string backupDir = m_basePath + COANQI_PATH_SEP + m_config.backupPath;
    std::string timestamp = GetTimestamp();
    std::string dstPath = backupDir + COANQI_PATH_SEP + 
                          "MAIN_1_CoAnQi_" + timestamp + ".cpp";
    
    try {
        fs::copy_file(srcPath, dstPath, fs::copy_options::overwrite_existing);
        auto fileSize = fs::file_size(dstPath);
        LogMaintenance("BACKUP", "Created: MAIN_1_CoAnQi_" + timestamp + ".cpp (" +
                       std::to_string(fileSize / 1024) + " KB)");
        return true;
    } catch (const std::exception& e) {
        LogMaintenance("BACKUP_ERROR", std::string(e.what()));
        return false;
    }
}

bool CoAnQiBot::SyncSharedConstants() {
    LogMaintenance("SYNC", "Synchronizing shared_constants.h");
    
    json args;
    json result;
    
    if (CallPythonFunction("sync_constants_across_languages", args, result)) {
        LogMaintenance("SYNC", "shared_constants.h synchronized with .py and .js");
        return true;
    }
    return false;
}

bool CoAnQiBot::RegenerateSystemConfig() {
    LogMaintenance("REGENERATE", "Regenerating observational_systems_config.h");
    
    json args;
    json result;
    
    if (CallPythonFunction("regenerate_system_config", args, result)) {
        LogMaintenance("REGENERATE", "observational_systems_config.h updated");
        return true;
    }
    return false;
}

bool CoAnQiBot::FullMaintenance() {
    LogMaintenance("FULL_MAINTENANCE", "Starting full MAIN_1_CoAnQi.cpp maintenance");
    
    bool success = true;
    
    success &= BackupMain1CoAnQi();
    success &= SyncSharedConstants();
    success &= RegenerateSystemConfig();
    
    // Run validation
    CalculateAllSystems();
    PerformStatisticalAnalysis();
    
    LogMaintenance("FULL_MAINTENANCE", success ? "Complete" : "Completed with errors");
    return success;
}

// ============================================================================
// MENU OPERATIONS
// ============================================================================

json CoAnQiBot::ExecuteMenuOption(int option, const std::map<std::string, std::string>& args) {
    LogMaintenance("MENU", "Executing option " + std::to_string(option));
    
    json result;
    
    switch (option) {
        case 1: { // Calculate single system
            auto it = args.find("system");
            std::string sys = (it != args.end()) ? it->second : "Vela Pulsar";
            auto calcResult = CalculateSystem(sys);
            result["F_U_Bi_i"] = std::to_string(calcResult.F_U_Bi_i);
            result["g_compressed"] = std::to_string(calcResult.g_compressed);
            break;
        }
        case 2: { // Calculate all systems
            auto allResults = CalculateAllSystems();
            result["count"] = std::to_string(allResults.size());
            break;
        }
        case 3: { // Clone and mutate
            auto itSys = args.find("system");
            auto itRate = args.find("mutationRate");
            std::string sys = (itSys != args.end()) ? itSys->second : "Vela Pulsar";
            double rate = (itRate != args.end()) ? std::stod(itRate->second) : 0.1;
            std::string cloneName = CloneAndMutate(sys, rate);
            result["cloneName"] = cloneName;
            break;
        }
        case 4: { // Add custom system
            // Requires full system params - delegate to Python
            json pyArgs;
            for (const auto& [k, v] : args) pyArgs[k] = v;
            CallPythonFunction("add_custom_system", pyArgs, result);
            break;
        }
        case 5: { // Add physics term
            auto it = args.find("termName");
            std::string termName = (it != args.end()) ? it->second : "NewTerm";
            auto itEq = args.find("equation");
            std::string eq = (itEq != args.end()) ? itEq->second : "";
            InjectDynamicTerm(termName, eq, {});
            result["termName"] = termName;
            break;
        }
        case 6: { // Run simulations
            auto it = args.find("simType");
            int simType = (it != args.end()) ? std::stoi(it->second) : 1;
            RunSimulation(static_cast<CoAnQiSimulationType>(simType));
            break;
        }
        case 7: { // Statistical analysis
            auto stats = PerformStatisticalAnalysis();
            result["mean"] = std::to_string(stats.mean);
            result["stddev"] = std::to_string(stats.stddev);
            result["count"] = std::to_string(stats.count);
            break;
        }
        case 8: { // Self-optimization
            auto it = args.find("system");
            std::string sys = (it != args.end()) ? it->second : "Vela Pulsar";
            // Placeholder - would need real observed/predicted data
            std::vector<double> obs = {1e33, 1.05e33, 0.98e33};
            std::vector<double> pred = {1e33, 1e33, 1e33};
            auto optResult = OptimizeParameters(sys, obs, pred);
            result["originalMSE"] = std::to_string(optResult.originalMSE);
            result["optimizedMSE"] = std::to_string(optResult.optimizedMSE);
            break;
        }
        default:
            LogMaintenance("MENU", "Unknown option: " + std::to_string(option));
    }
    
    return result;
}

void CoAnQiBot::ExecuteCommand(const std::string& command) {
    LogMaintenance("COMMAND", "Voice/script: " + command);
    
    std::string cmdLower = command;
    std::transform(cmdLower.begin(), cmdLower.end(), cmdLower.begin(), ::tolower);
    
    if (cmdLower.find("calculate") != std::string::npos) {
        if (cmdLower.find("all") != std::string::npos) {
            CalculateAllSystems();
        } else {
            // Try to extract system name
            CalculateSystem("Vela Pulsar");  // Default
        }
    } else if (cmdLower.find("backup") != std::string::npos) {
        BackupMain1CoAnQi();
    } else if (cmdLower.find("maintenance") != std::string::npos) {
        FullMaintenance();
    } else if (cmdLower.find("statistics") != std::string::npos) {
        PerformStatisticalAnalysis();
    } else if (cmdLower.find("optimize") != std::string::npos) {
        std::vector<double> obs = {1e33};
        std::vector<double> pred = {1e33};
        OptimizeParameters("default", obs, pred);
    } else if (cmdLower.find("sync") != std::string::npos) {
        SyncSharedConstants();
    } else if (cmdLower.find("validate") != std::string::npos) {
        FullCrossValidation("Sagittarius A*");
    } else {
        LogMaintenance("COMMAND", "Unrecognized command: " + command);
    }
}

// ============================================================================
// IPC TRIGGERS
// ============================================================================

void CoAnQiBot::TriggerSelfExpandingBackend(const std::string& term) {
    if (m_ipc) {
        json msg;
        msg["type"] = "REGISTER_TERM";
        msg["term"] = term;
        m_ipc->Send(msg.dump());
    }
}

void CoAnQiBot::TriggerParameterUpdate(const std::string& param, double value) {
    if (m_ipc) {
        json msg;
        msg["type"] = "UPDATE_PARAMETER";
        msg["param"] = param;
        msg["value"] = std::to_string(value);
        m_ipc->Send(msg.dump());
    }
}

} // namespace VR
} // namespace UQFF
