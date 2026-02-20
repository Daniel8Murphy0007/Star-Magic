/**
 * @file python_bridge.h
 * @brief Python Integration Bridge for UQFF Physics Service
 * 
 * Provides pybind11-based integration with:
 * - CondensedPhysics.py (81K lines, 8 UQFF Master Equations)
 * - CondensedPhysics_InputData.py (3249 lines, observational parameters)
 * - CondensedPhysics_OutputData.py (5378 lines, query result storage)
 * - CondensedPhysics_Validation.py (5926 lines, validation references)
 * - QCalc.py (9K lines, Pure Physics Solver)
 * - IPData.py / OPData.py (Data layer)
 * - APIFetch.py (SIMBAD/NASA parameter fetching)
 * - Phase5_Consolidated.py (838 lines, SOURCE52-65: 57 funcs, 41 systems)
 * - Phase6_Consolidated.py (487 lines, SOURCE70-80: 31 funcs, 3 systems)
 * - Phase7_Consolidated.py (3581 lines, SOURCE81-95: 110 funcs, 14 systems)
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 2 Fix - Python Calculator Integration
 */

#ifndef PYTHON_BRIDGE_H
#define PYTHON_BRIDGE_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <mutex>

// Forward declare pybind11 types to avoid header dependency
namespace pybind11 {
    class scoped_interpreter;
    class object;
    class dict;
    class module_;
}

namespace UQFF {

/**
 * @struct PythonFieldResult
 * @brief Result from Python calculator invocation
 */
struct PythonFieldResult {
    bool success = false;
    std::string error_message;
    
    // UQFF Field Components (from QCalc/CondensedPhysics)
    double F_U = 0;             // Total unified field
    double Ug1 = 0, Ug2 = 0, Ug3 = 0, Ug4 = 0;  // Gravity components
    double Um = 0;              // Magnetism
    double Ubi = 0;             // Buoyancy opposition
    double g_compressed = 0;    // MUGE Compressed gravity
    double g_resonant = 0;      // MUGE Resonant gravity
    
    // Extended results
    std::map<std::string, double> additional_fields;
    
    // Long-form equation output
    std::vector<std::string> long_form_equations;
    
    // Available equations for this system
    std::vector<std::string> available_equations;
    
    // Timing
    double compute_time_ms = 0;
};

/**
 * @struct PythonQueryParams
 * @brief Parameters for Python calculator query
 */
struct PythonQueryParams {
    std::string system_name;    // e.g., "Sagittarius A*", "M87", "Betelgeuse"
    double r = 0;               // Radial distance [m]
    double t = 0;               // Time [s]
    double theta = 0;           // Angle [rad]
    
    // Optional manual overrides (0 = use API fetch)
    double mass = 0;
    double radius = 0;
    double distance = 0;
    double magnetic_field = 0;
    double spin_period = 0;
    double redshift = 0;
    
    // Calculation options
    bool fetch_from_api = true;         // Use APIFetch.py
    bool compute_all_equations = false; // Compute all 8 master equations
    bool include_long_form = true;      // Include equation strings
};

/**
 * @class PythonBridge
 * @brief Singleton bridge to Python calculator ecosystem
 * 
 * Thread-safe interface to:
 * - CondensedPhysics.py (master calculator)
 * - QCalc.py (pure physics solver)
 * - APIFetch.py (parameter fetching)
 * - IPData.py / OPData.py (data layer)
 */
class PythonBridge {
public:
    /**
     * @brief Get singleton instance
     */
    static PythonBridge& instance();
    
    /**
     * @brief Initialize Python interpreter and load modules
     * @return true if initialization successful
     */
    bool initialize();
    
    /**
     * @brief Check if Python bridge is ready
     */
    bool is_ready() const { return initialized_; }
    
    /**
     * @brief Shutdown Python interpreter
     */
    void shutdown();
    
    /**
     * @brief Calculate field using CondensedPhysics.py
     * @param params Query parameters
     * @return Field calculation result
     */
    PythonFieldResult calculate_condensed_physics(const PythonQueryParams& params);
    
    /**
     * @brief Calculate field using QCalc.py
     * @param params Query parameters
     * @return Field calculation result
     */
    PythonFieldResult calculate_qcalc(const PythonQueryParams& params);
    
    /**
     * @brief Fetch system parameters from APIs (SIMBAD, NASA)
     * @param system_name System identifier
     * @return Map of parameter name -> value
     */
    std::map<std::string, double> fetch_parameters(const std::string& system_name);
    
    /**
     * @brief Get list of available systems from bodies.csv
     * @return Vector of system names
     */
    std::vector<std::string> get_available_systems();
    
    /**
     * @brief Execute arbitrary Python code (for testing/debugging)
     * @param code Python code string
     * @return Result as string
     */
    std::string execute_python(const std::string& code);
    
    /**
     * @brief Get Python module versions
     */
    struct ModuleVersions {
        std::string condensed_physics;
        std::string qcalc;
        std::string numpy;
        std::string scipy;
        std::string sympy;
    };
    ModuleVersions get_module_versions();
    
private:
    PythonBridge() = default;
    ~PythonBridge();
    
    // Non-copyable
    PythonBridge(const PythonBridge&) = delete;
    PythonBridge& operator=(const PythonBridge&) = delete;
    
    bool initialized_ = false;
    std::mutex python_mutex_;  // GIL wrapper for thread safety
    
    // Cached module references (pimpl pattern to hide pybind11)
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

/**
 * @brief Check if Python dependencies are available
 * @return true if CondensedPhysics.py, QCalc.py, etc. are importable
 */
bool check_python_dependencies();

/**
 * @brief Get path to Python virtual environment
 * @return Path string (e.g., ".venv/Scripts/python.exe")
 */
std::string get_python_executable();

} // namespace UQFF

#endif // PYTHON_BRIDGE_H
