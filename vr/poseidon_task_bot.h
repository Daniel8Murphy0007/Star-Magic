/**
 * @file poseidon_task_bot.h
 * @brief Poseidon TaskBot - Offline Physics Maintenance Agent for UQFF Star-Magic
 * @version 4.2.1 (Canonical - matches Architecture v4.2)
 * 
 * PURPOSE: Offline-first physics task automation bot that lives inside the VR/VM Backend
 * (source2(HEAD PROGRAM).cpp). It reads/writes/compares/validates/updates physics
 * calculations across all languages (C++ / Python / JS) and maintains the entire
 * codebase cross-platform while respecting the Recirculation Loop and Self-Expanding
 * Backend (v3.1).
 * 
 * CAPABILITIES:
 * - Process new physics equations (read/parse/integrate)
 * - Write validated physics to appropriate calculators
 * - Compare results across C++/Python/JS implementations
 * - Validate physics using cross-validation pipeline
 * - Update Self-Expanding Backend (v3.1) parameters
 * - Cross-platform code maintenance (constants sync, file regeneration)
 * - Offline FTPS bundle push/pull for distributed teams
 * 
 * INTEGRATION:
 * - Uses physics_service.h (Self-Expand/Update/Simulate v3.1)
 * - Uses python_bridge.h (pybind11 to Python ecosystem)
 * - Uses uqff_cross_platform.h (harmonized implementations)
 * - Calls uqff_ftps_client.py for secure offline transfer
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v4.2
 * Phase: 3.1+ - Poseidon Offline Maintenance
 */

#ifndef POSEIDON_TASK_BOT_H
#define POSEIDON_TASK_BOT_H

// Cross-platform compatibility
#ifdef _WIN32
    #include <windows.h>
    #include <direct.h>
    #define POSEIDON_PATH_SEP "\\"
    #define poseidon_mkdir(path) _mkdir(path)
#else
    #include <unistd.h>
    #include <sys/stat.h>
    #define POSEIDON_PATH_SEP "/"
    #define poseidon_mkdir(path) mkdir(path, 0755)
#endif

#include "../ipc/physics_service.h"
#include "../ipc/python_bridge.h"
#include "../ipc/uqff_ipc.h"
#include "../uqff_cross_platform.h"
#include "../shared_constants.h"
#include "../csv_body_reader.h"
#include "../uqff_tracing.h"

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <mutex>
#include <atomic>

// JSON support (nlohmann/json or fallback)
#ifdef HAVE_NLOHMANN_JSON
    #include <nlohmann/json.hpp>
    using json = nlohmann::json;
#else
    // Minimal JSON implementation fallback
    #include <sstream>
    namespace poseidon_json {
        struct json {
            std::map<std::string, std::string> data;
            std::string& operator[](const std::string& key) { return data[key]; }
            std::string dump(int = 0) const {
                std::ostringstream oss;
                oss << "{";
                bool first = true;
                for (const auto& [k, v] : data) {
                    if (!first) oss << ",";
                    oss << "\"" << k << "\":\"" << v << "\"";
                    first = false;
                }
                oss << "}";
                return oss.str();
            }
        };
    }
    using json = poseidon_json::json;
#endif

namespace UQFF {
namespace VR {

// ============================================================================
// FORWARD DECLARATIONS
// ============================================================================
class VRRuntime;

// ============================================================================
// POSEIDON ENUMS
// ============================================================================

/**
 * @enum PoseidonTaskType
 * @brief Types of maintenance tasks Poseidon can execute
 */
enum class PoseidonTaskType {
    READ_PHYSICS,           // Read from any source (CSV, JSON, code)
    WRITE_PHYSICS,          // Write to calculators/configs
    COMPARE_CROSS_LANG,     // Compare C++/Python/JS results
    VALIDATE_PHYSICS,       // Run validation pipeline
    UPDATE_BACKEND,         // Update Self-Expanding parameters
    SYNC_CONSTANTS,         // Synchronize shared_constants across languages
    REGENERATE_FILES,       // Regenerate extracted files
    BACKUP_PHYSICS,         // Create timestamped backup
    FTPS_PUSH,              // Push maintenance bundle via FTPS
    FTPS_PULL,              // Pull latest physics via FTPS
    FULL_MAINTENANCE,       // Run complete maintenance cycle
    CUSTOM                  // User-defined task
};

/**
 * @enum PoseidonSourceLang
 * @brief Source languages for physics implementations
 */
enum class PoseidonSourceLang {
    CPP,        // C++ (MAIN_1_CoAnQi.cpp, source*.cpp)
    PYTHON,     // Python (CondensedPhysics.py, QCalc*.py)
    JAVASCRIPT, // JavaScript (index.js)
    LATEX,      // LaTeX (input equations)
    WOLFRAM,    // Wolfram Mathematica
    ALL         // Cross-language operation
};

// ============================================================================
// POSEIDON STRUCTURES
// ============================================================================

/**
 * @struct PoseidonComparisonResult
 * @brief Result of cross-language physics comparison
 */
struct PoseidonComparisonResult {
    std::string equationName;
    double cppValue = 0.0;
    double pyValue = 0.0;
    double jsValue = 0.0;
    double maxDeviation = 0.0;
    double relativeError = 0.0;      // |max-min| / |average|
    bool passes = false;             // Within tolerance
    std::string report;
    
    // Per-language success flags
    bool cppAvailable = false;
    bool pyAvailable = false;
    bool jsAvailable = false;
    
    // Timing
    double cppTimeMs = 0.0;
    double pyTimeMs = 0.0;
    double jsTimeMs = 0.0;
};

/**
 * @struct PoseidonValidationResult
 * @brief Result of physics validation
 */
struct PoseidonValidationResult {
    bool success = false;
    std::string systemName;
    int testsRun = 0;
    int testsPassed = 0;
    int testsFailed = 0;
    double validationScore = 0.0;    // 0-100
    std::vector<std::string> failedTests;
    std::vector<std::string> warnings;
    std::string fullReport;
    double durationMs = 0.0;
};

/**
 * @struct PoseidonTaskResult
 * @brief Generic result of any Poseidon task
 */
struct PoseidonTaskResult {
    bool success = false;
    PoseidonTaskType taskType = PoseidonTaskType::CUSTOM;
    std::string message;
    std::string data;                // JSON or other result data
    double executionTimeMs = 0.0;
    int itemsProcessed = 0;
    int itemsFailed = 0;
};

/**
 * @struct PoseidonMaintenanceBundle
 * @brief Bundle for FTPS transfer containing maintenance data
 */
struct PoseidonMaintenanceBundle {
    std::string version = "4.2.1";
    uint64_t timestamp = 0;
    std::string sourceHost;
    
    // Files to include
    std::vector<std::string> constantFiles;     // shared_constants.h, .py
    std::vector<std::string> physicsFiles;      // Updated physics files
    std::vector<std::string> configFiles;       // Configuration updates
    std::vector<std::string> logFiles;          // Maintenance logs
    
    // Metadata
    std::map<std::string, std::string> metadata;
};

/**
 * @struct PoseidonConfig
 * @brief Configuration for Poseidon TaskBot
 */
struct PoseidonConfig {
    // Base paths (cross-platform)
    std::string basePath;
    std::string backupPath;
    std::string logPath;
    
    // Tolerance thresholds
    double comparisonTolerance = 1e-10;  // Absolute tolerance
    double relativeTolerance = 1e-6;      // Relative tolerance (1 ppm)
    
    // FTPS settings
    std::string ftpsHost = "localhost";
    int ftpsPort = 990;
    std::string ftpsRemotePath = "/uqff/maintenance/";
    bool ftpsEnabled = false;
    
    // Maintenance settings
    bool autoBackup = true;
    bool autoValidate = true;
    bool verboseLogging = true;
    int backupRetentionDays = 30;
    
    // Integration flags
    bool enablePython = true;
    bool enableJavaScript = true;
    bool enableSelfExpanding = true;
    
    // Parse from JSON config file
    static PoseidonConfig load(const std::string& configPath);
    void save(const std::string& configPath) const;
};

// ============================================================================
// POSEIDON TASK BOT CLASS
// ============================================================================

/**
 * @class PoseidonTaskBot
 * @brief Offline physics maintenance bot for UQFF Star-Magic
 * 
 * Lives inside VR/VM Backend (source2(HEAD PROGRAM).cpp) and handles:
 * - Reading new physics equations (LaTeX, code, observational data)
 * - Writing validated physics to all calculators
 * - Comparing results across C++/Python/JS
 * - Validating using cross-validation pipeline
 * - Updating Self-Expanding Backend (v3.1) parameters
 * - Cross-platform code maintenance
 */
class PoseidonTaskBot {
public:
    /**
     * Constructor
     * @param basePath Base path for UQFF workspace (defaults to CrossPlatform::GetDataDir())
     */
    explicit PoseidonTaskBot(const std::string& basePath = "");
    ~PoseidonTaskBot();
    
    // Disable copy (unique resources)
    PoseidonTaskBot(const PoseidonTaskBot&) = delete;
    PoseidonTaskBot& operator=(const PoseidonTaskBot&) = delete;
    
    // ========================================================================
    // CORE OFFLINE OPERATIONS
    // ========================================================================
    
    /**
     * Process new physics equation for integration
     * @param equationName Name of the equation (e.g., "F_U_Bi_i", "compressed_g")
     * @param latexOrCode Equation content (LaTeX, code snippet)
     * @param sourceLang Source language ("C++", "Python", "JS", "LaTeX")
     * @return true if successfully integrated
     */
    bool ProcessNewPhysics(const std::string& equationName,
                           const std::string& latexOrCode,
                           PoseidonSourceLang sourceLang);
    
    // ========================================================================
    // READ OPERATIONS
    // ========================================================================
    
    /**
     * Read latest input parameters from bodies_*.csv via IPData.py
     */
    json ReadIPData();
    
    /**
     * Read computation results from uqff_results.json / OPData.py
     */
    json ReadOPData();
    
    /**
     * Read physics equations from a specific source file
     * @param sourceFile File path (e.g., "source115.cpp", "QCalc.py")
     * @return JSON with extracted equation signatures
     */
    json ReadPhysicsFile(const std::string& sourceFile);
    
    /**
     * Load latest bodies CSV into vector
     */
    bool LoadLatestBodiesCSV(std::vector<CelestialBody>& bodies);
    
    // ========================================================================
    // WRITE OPERATIONS
    // ========================================================================
    
    /**
     * Write computation results to OPData / uqff_results.json
     */
    void WriteOPData(const json& results);
    
    /**
     * Write to CondensedPhysics_OutputData.py for user recall
     */
    void WriteRecallStorage(const json& results);
    
    // ========================================================================
    // COMPARE OPERATIONS (Cross-Language)
    // ========================================================================
    
    /**
     * Compare all calculator implementations for a given system
     * @param systemName Astronomical system (e.g., "Sagittarius A*", "M87")
     * @return Vector of comparison results for each equation
     */
    std::vector<PoseidonComparisonResult> CompareAllCalculators(
        const std::string& systemName);
    
    /**
     * Compare a specific equation across languages
     */
    PoseidonComparisonResult CompareEquation(
        const std::string& equationName,
        const std::map<std::string, double>& params);
    
    // ========================================================================
    // VALIDATE OPERATIONS
    // ========================================================================
    
    /**
     * Run validation pipeline on a system
     * @param systemName Target system
     * @param runFullSuite Run all tests (vs quick validation)
     */
    PoseidonValidationResult ValidatePhysics(
        const std::string& systemName, 
        bool runFullSuite = true);
    
    /**
     * Cross-validate UQFF vs MUGE vs GR
     */
    PoseidonValidationResult CrossValidateFrameworks(
        const std::string& systemName);
    
    // ========================================================================
    // UPDATE OPERATIONS (Self-Expanding v3.1)
    // ========================================================================
    
    /**
     * Update Self-Expanding Backend with new physics term
     * @param newTerm Name of new physics term
     * @param kappa Optional kappa parameter override
     */
    bool UpdateAndExpandPhysics(const std::string& newTerm, double kappa = 0.0);
    
    /**
     * Register dynamic physics term with backend
     */
    bool RegisterDynamicTerm(const DynamicTerm& term);
    
    /**
     * Update calibrated parameters (κ, [SSq], β_i, etc.)
     */
    bool UpdateCalibratedParameters(const CalibratedParameters& params);
    
    // ========================================================================
    // CODE MAINTENANCE (Cross-Platform)
    // ========================================================================
    
    /**
     * Run full codebase maintenance cycle
     */
    PoseidonTaskResult MaintainCodebase();
    
    /**
     * Synchronize constants across all languages
     * shared_constants.h <-> shared_constants.py <-> index.js CONSTANTS
     */
    bool SyncConstantsAcrossLanguages();
    
    /**
     * Regenerate extracted files
     * QCalc_cpp_extracted.py, QCalc_js_extracted.py, etc.
     */
    bool RegenerateExtractedFiles();
    
    /**
     * Create timestamped backup of all physics files
     */
    bool BackupAllPhysicsFiles();
    
    /**
     * Verify file integrity across the codebase
     */
    PoseidonTaskResult VerifyCodebaseIntegrity();
    
    // ========================================================================
    // FTPS BRIDGE (Offline Secure Transfer)
    // ========================================================================
    
    /**
     * Push maintenance bundle via FTPS
     * @param remotePath Target path on FTPS server
     */
    bool FTPSPushMaintenanceBundle(const std::string& remotePath = "");
    
    /**
     * Pull latest physics via FTPS
     */
    bool FTPSPullLatestPhysics(const std::string& remotePath = "");
    
    /**
     * Create maintenance bundle (without transfer)
     */
    PoseidonMaintenanceBundle CreateMaintenanceBundle();
    
    // ========================================================================
    // COMMAND INTERFACE (Voice/Script)
    // ========================================================================
    
    /**
     * Execute command from voice or script
     * @param command Natural language or structured command
     */
    PoseidonTaskResult ExecuteCommand(const std::string& command);
    
    /**
     * Register custom command handler
     */
    using CommandHandler = std::function<PoseidonTaskResult(const std::string&)>;
    void RegisterCommand(const std::string& trigger, CommandHandler handler);
    
    // ========================================================================
    // STATUS & DIAGNOSTICS
    // ========================================================================
    
    struct Status {
        bool initialized = false;
        bool pythonAvailable = false;
        bool jsAvailable = false;
        bool ftpsAvailable = false;
        int tasksCompleted = 0;
        int tasksFailed = 0;
        std::string lastTask;
        std::string lastError;
        uint64_t lastMaintenanceTimestamp = 0;
    };
    
    Status GetStatus() const { return m_status; }
    
    /**
     * Get maintenance log
     */
    std::string GetMaintenanceLog() const;
    
    /**
     * Get configuration
     */
    PoseidonConfig GetConfig() const { return m_config; }
    
    /**
     * Update configuration at runtime
     */
    void SetConfig(const PoseidonConfig& config);

private:
    // Configuration
    PoseidonConfig m_config;
    std::string m_basePath;
    Status m_status;
    
    // Integration handles
    std::unique_ptr<PythonBridge> m_pyBridge;
    std::unique_ptr<UQFF::PhysicsService> m_physicsService;
    std::unique_ptr<IPC::NamedPipeChannel> m_ipcChannel;
    
    // Maintenance log
    std::vector<std::string> m_maintenanceLog;
    mutable std::mutex m_logMutex;
    
    // Command handlers
    std::map<std::string, CommandHandler> m_commandHandlers;
    
    // Internal helpers
    void Initialize();
    void LogMaintenance(const std::string& action, const std::string& details);
    void TriggerSelfExpandingBackend(const std::string& term);
    std::string GetTimestamp() const;
    std::string FindLatestBodiesCSV() const;
    bool EnsureDirectories();
    
    // Python bridge helpers
    json CallPythonCalculator(const std::string& module, 
                              const std::string& function,
                              const json& params);
    
    // JavaScript bridge helper (via uqff_server.js)
    json CallJavaScriptCalculator(const std::string& function,
                                   const json& params);
    
    // Built-in command registration
    void RegisterBuiltinCommands();
};

// ============================================================================
// POSEIDON COMMAND INTENTS (for voice/script parsing)
// ============================================================================

namespace PoseidonIntent {
    // Physics Operations
    constexpr const char* READ_PHYSICS = "read_physics";
    constexpr const char* WRITE_PHYSICS = "write_physics";
    constexpr const char* COMPARE = "compare";
    constexpr const char* VALIDATE = "validate";
    constexpr const char* UPDATE = "update";
    
    // Maintenance Operations
    constexpr const char* SYNC_CONSTANTS = "sync_constants";
    constexpr const char* REGENERATE = "regenerate";
    constexpr const char* BACKUP = "backup";
    constexpr const char* FULL_MAINTENANCE = "full_maintenance";
    
    // Transfer Operations
    constexpr const char* FTPS_PUSH = "ftps_push";
    constexpr const char* FTPS_PULL = "ftps_pull";
    
    // Status
    constexpr const char* STATUS = "status";
    constexpr const char* HELP = "help";
}

} // namespace VR
} // namespace UQFF

#endif // POSEIDON_TASK_BOT_H
