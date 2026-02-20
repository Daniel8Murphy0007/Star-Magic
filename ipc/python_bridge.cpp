/**
 * @file python_bridge.cpp
 * @brief Python Integration Bridge Implementation
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 2 Fix - Python Calculator Integration
 */

#include "python_bridge.h"

// PYTHON_BRIDGE_ENABLED is set when pybind11 is found
// This is separate from NO_PYTHON which controls source2.cpp's Qt-conflicting code
#ifdef PYTHON_BRIDGE_ENABLED
#include <pybind11/embed.h>
#include <pybind11/stl.h>
namespace py = pybind11;
#endif

#include <iostream>
#include <chrono>
#include <filesystem>
#include <fstream>

namespace UQFF {

// ============================================================================
// PYTHON BRIDGE IMPLEMENTATION
// ============================================================================

#ifdef PYTHON_BRIDGE_ENABLED

struct PythonBridge::Impl {
    std::unique_ptr<py::scoped_interpreter> interpreter;
    py::object condensed_physics_module;
    py::object condensed_input_module;   // CondensedPhysics_InputData.py (3249 lines)
    py::object condensed_output_module;  // CondensedPhysics_OutputData.py (5378 lines)
    py::object condensed_validation_module; // CondensedPhysics_Validation.py (5926 lines)
    py::object qcalc_module;
    // QCalc ecosystem (12 modules, ~157K lines)
    py::object qcalc_cpp_equations;      // 49984 lines - ALL equations from C++
    py::object qcalc_extracted;          // 43145 lines - Auto-generated extractors
    py::object qcalc_cpp_extracted;      // 30903 lines - C++ extracted physics
    py::object qcalc_js_extracted;       // 13367 lines - JS extracted physics
    py::object qcalc_wolfram_ext;        // 5561 lines - Wolfram functions
    py::object qcalc_core_uqff;          // 1223 lines - Core UQFF equations
    py::object qcalc_validation;         // 1109 lines - Validation data fetcher
    py::object qcalc_advanced;           // 789 lines - Advanced physics
    py::object qcalc_stat;               // 719 lines - Statistical analysis
    py::object qcalc_phase1_val;         // 527 lines - Phase 1 validation
    py::object qcalc_api;                // 492 lines - REST API
    py::object qcalc_performance;        // 480 lines - Performance optimization
    py::object apifetch_module;
    py::object ipdata_module;
    py::object opdata_module;
    py::object phase5_module;  // SOURCE52-65: 57 funcs, 41 systems
    py::object phase6_module;  // SOURCE70-80: 31 funcs, 3 systems
    py::object phase7_module;  // SOURCE81-95: 110 funcs, 14 systems
};

PythonBridge& PythonBridge::instance() {
    static PythonBridge instance;
    return instance;
}

PythonBridge::~PythonBridge() {
    shutdown();
}

bool PythonBridge::initialize() {
    std::lock_guard<std::mutex> lock(python_mutex_);
    
    if (initialized_) {
        return true;
    }
    
    try {
        impl_ = std::make_unique<Impl>();
        
        // Start Python interpreter
        impl_->interpreter = std::make_unique<py::scoped_interpreter>();
        
        // Add current directory to Python path
        py::module_ sys = py::module_::import("sys");
        py::list path = sys.attr("path").cast<py::list>();
        path.insert(0, ".");
        
        std::cout << "[PythonBridge] Interpreter started, loading modules..." << std::endl;
        
        // Import core modules
        try {
            impl_->condensed_physics_module = py::module_::import("CondensedPhysics");
            std::cout << "[PythonBridge] CondensedPhysics.py loaded (81K lines)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: CondensedPhysics.py failed: " << e.what() << std::endl;
        }
        
        // CondensedPhysics support modules
        try {
            impl_->condensed_input_module = py::module_::import("CondensedPhysics_InputData");
            std::cout << "[PythonBridge] CondensedPhysics_InputData.py loaded (3249 lines)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: CondensedPhysics_InputData.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->condensed_output_module = py::module_::import("CondensedPhysics_OutputData");
            std::cout << "[PythonBridge] CondensedPhysics_OutputData.py loaded (5378 lines)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: CondensedPhysics_OutputData.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->condensed_validation_module = py::module_::import("CondensedPhysics_Validation");
            std::cout << "[PythonBridge] CondensedPhysics_Validation.py loaded (5926 lines)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: CondensedPhysics_Validation.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->qcalc_module = py::module_::import("QCalc");
            std::cout << "[PythonBridge] QCalc.py loaded (9K lines)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: QCalc.py failed: " << e.what() << std::endl;
        }
        
        // QCalc ecosystem modules (~157K lines total)
        try {
            impl_->qcalc_cpp_equations = py::module_::import("QCalc_cpp_equations");
            std::cout << "[PythonBridge] QCalc_cpp_equations.py loaded (50K lines)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: QCalc_cpp_equations.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->qcalc_extracted = py::module_::import("QCalc_extracted");
            std::cout << "[PythonBridge] QCalc_extracted.py loaded (43K lines)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: QCalc_extracted.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->qcalc_cpp_extracted = py::module_::import("QCalc_cpp_extracted");
            std::cout << "[PythonBridge] QCalc_cpp_extracted.py loaded (31K lines)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: QCalc_cpp_extracted.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->qcalc_js_extracted = py::module_::import("QCalc_js_extracted");
            std::cout << "[PythonBridge] QCalc_js_extracted.py loaded (13K lines)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: QCalc_js_extracted.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->qcalc_wolfram_ext = py::module_::import("QCalc_Wolfram_Extensions");
            std::cout << "[PythonBridge] QCalc_Wolfram_Extensions.py loaded (6K lines)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: QCalc_Wolfram_Extensions.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->qcalc_core_uqff = py::module_::import("QCalc_core_uqff");
            std::cout << "[PythonBridge] QCalc_core_uqff.py loaded (1K lines)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: QCalc_core_uqff.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->qcalc_validation = py::module_::import("QCalc_validation");
            std::cout << "[PythonBridge] QCalc_validation.py loaded" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: QCalc_validation.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->qcalc_advanced = py::module_::import("QCalc_Advanced");
            std::cout << "[PythonBridge] QCalc_Advanced.py loaded" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: QCalc_Advanced.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->qcalc_stat = py::module_::import("QCalc_stat");
            std::cout << "[PythonBridge] QCalc_stat.py loaded" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: QCalc_stat.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->qcalc_phase1_val = py::module_::import("QCalc_Phase1_Validation");
            std::cout << "[PythonBridge] QCalc_Phase1_Validation.py loaded" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: QCalc_Phase1_Validation.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->qcalc_api = py::module_::import("QCalc_API");
            std::cout << "[PythonBridge] QCalc_API.py loaded" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: QCalc_API.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->qcalc_performance = py::module_::import("QCalc_Performance");
            std::cout << "[PythonBridge] QCalc_Performance.py loaded" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: QCalc_Performance.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->apifetch_module = py::module_::import("APIFetch");
            std::cout << "[PythonBridge] APIFetch.py loaded" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: APIFetch.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->ipdata_module = py::module_::import("IPData");
            std::cout << "[PythonBridge] IPData.py loaded" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: IPData.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->opdata_module = py::module_::import("OPData");
            std::cout << "[PythonBridge] OPData.py loaded" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: OPData.py failed: " << e.what() << std::endl;
        }
        
        // Phase 5-7 Consolidated modules (extracted C++ SOURCE files)
        try {
            impl_->phase5_module = py::module_::import("Phase5_Consolidated");
            std::cout << "[PythonBridge] Phase5_Consolidated.py loaded (SOURCE52-65: 57 funcs)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: Phase5_Consolidated.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->phase6_module = py::module_::import("Phase6_Consolidated");
            std::cout << "[PythonBridge] Phase6_Consolidated.py loaded (SOURCE70-80: 31 funcs)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: Phase6_Consolidated.py failed: " << e.what() << std::endl;
        }
        
        try {
            impl_->phase7_module = py::module_::import("Phase7_Consolidated");
            std::cout << "[PythonBridge] Phase7_Consolidated.py loaded (SOURCE81-95: 110 funcs)" << std::endl;
        } catch (const py::error_already_set& e) {
            std::cerr << "[PythonBridge] WARNING: Phase7_Consolidated.py failed: " << e.what() << std::endl;
        }
        
        initialized_ = true;
        std::cout << "[PythonBridge] Initialization complete" << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "[PythonBridge] FATAL: " << e.what() << std::endl;
        return false;
    }
}

void PythonBridge::shutdown() {
    std::lock_guard<std::mutex> lock(python_mutex_);
    
    if (!initialized_) return;
    
    impl_.reset();
    initialized_ = false;
    
    std::cout << "[PythonBridge] Shutdown complete" << std::endl;
}

PythonFieldResult PythonBridge::calculate_condensed_physics(const PythonQueryParams& params) {
    std::lock_guard<std::mutex> lock(python_mutex_);
    
    PythonFieldResult result;
    
    if (!initialized_ || !impl_ || impl_->condensed_physics_module.is_none()) {
        result.error_message = "CondensedPhysics module not loaded";
        return result;
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    try {
        // Build parameter dict
        py::dict py_params;
        py_params["system_name"] = params.system_name;
        py_params["r"] = params.r;
        py_params["t"] = params.t;
        py_params["theta"] = params.theta;
        
        if (params.mass > 0) py_params["mass"] = params.mass;
        if (params.radius > 0) py_params["radius"] = params.radius;
        if (params.distance > 0) py_params["distance"] = params.distance;
        if (params.magnetic_field > 0) py_params["magnetic_field"] = params.magnetic_field;
        if (params.spin_period > 0) py_params["spin_period"] = params.spin_period;
        
        // Call CondensedPhysics calculator
        // Look for UQFFCalculator class or compute_uqff function
        py::object calculator;
        
        if (py::hasattr(impl_->condensed_physics_module, "UQFFCalculator")) {
            py::object calc_class = impl_->condensed_physics_module.attr("UQFFCalculator");
            calculator = calc_class();
        }
        
        py::object py_result;
        
        if (py::hasattr(calculator, "compute")) {
            py_result = calculator.attr("compute")(py_params);
        } else if (py::hasattr(impl_->condensed_physics_module, "compute_uqff")) {
            py_result = impl_->condensed_physics_module.attr("compute_uqff")(py_params);
        } else if (py::hasattr(impl_->condensed_physics_module, "solve")) {
            py_result = impl_->condensed_physics_module.attr("solve")(params.system_name, py_params);
        } else {
            result.error_message = "No compute method found in CondensedPhysics";
            return result;
        }
        
        // Extract results
        if (py::isinstance<py::dict>(py_result)) {
            py::dict res = py_result.cast<py::dict>();
            
            if (res.contains("F_U")) result.F_U = res["F_U"].cast<double>();
            if (res.contains("Ug1")) result.Ug1 = res["Ug1"].cast<double>();
            if (res.contains("Ug2")) result.Ug2 = res["Ug2"].cast<double>();
            if (res.contains("Ug3")) result.Ug3 = res["Ug3"].cast<double>();
            if (res.contains("Ug4")) result.Ug4 = res["Ug4"].cast<double>();
            if (res.contains("Um")) result.Um = res["Um"].cast<double>();
            if (res.contains("Ubi")) result.Ubi = res["Ubi"].cast<double>();
            if (res.contains("g_compressed")) result.g_compressed = res["g_compressed"].cast<double>();
            if (res.contains("g_resonant")) result.g_resonant = res["g_resonant"].cast<double>();
            
            // Long form equations
            if (res.contains("long_form_equations")) {
                py::list eqs = res["long_form_equations"].cast<py::list>();
                for (auto& eq : eqs) {
                    result.long_form_equations.push_back(eq.cast<std::string>());
                }
            }
            
            // Available equations
            if (res.contains("available_equations")) {
                py::list avail = res["available_equations"].cast<py::list>();
                for (auto& eq : avail) {
                    result.available_equations.push_back(eq.cast<std::string>());
                }
            }
        }
        
        result.success = true;
        
    } catch (const py::error_already_set& e) {
        result.error_message = std::string("Python error: ") + e.what();
    } catch (const std::exception& e) {
        result.error_message = std::string("C++ error: ") + e.what();
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    result.compute_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    return result;
}

PythonFieldResult PythonBridge::calculate_qcalc(const PythonQueryParams& params) {
    std::lock_guard<std::mutex> lock(python_mutex_);
    
    PythonFieldResult result;
    
    if (!initialized_ || !impl_ || impl_->qcalc_module.is_none()) {
        result.error_message = "QCalc module not loaded";
        return result;
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    try {
        // Build parameter dict
        py::dict py_params;
        py_params["system_name"] = params.system_name;
        py_params["r"] = params.r;
        py_params["t"] = params.t;
        py_params["theta"] = params.theta;
        
        if (params.mass > 0) py_params["mass"] = params.mass;
        if (params.radius > 0) py_params["radius"] = params.radius;
        
        // Call QCalc.solve()
        py::object py_result;
        
        if (py::hasattr(impl_->qcalc_module, "solve")) {
            py_result = impl_->qcalc_module.attr("solve")(params.system_name, py_params);
        } else if (py::hasattr(impl_->qcalc_module, "compute")) {
            py_result = impl_->qcalc_module.attr("compute")(py_params);
        } else {
            result.error_message = "No solve/compute method in QCalc";
            return result;
        }
        
        // Extract results (similar to condensed physics)
        if (py::isinstance<py::dict>(py_result)) {
            py::dict res = py_result.cast<py::dict>();
            
            if (res.contains("solutions")) {
                py::dict sols = res["solutions"].cast<py::dict>();
                if (sols.contains("F_U")) result.F_U = sols["F_U"].cast<double>();
                if (sols.contains("Ug1")) result.Ug1 = sols["Ug1"].cast<double>();
                if (sols.contains("Ug2")) result.Ug2 = sols["Ug2"].cast<double>();
                if (sols.contains("Ug3")) result.Ug3 = sols["Ug3"].cast<double>();
                if (sols.contains("Ug4")) result.Ug4 = sols["Ug4"].cast<double>();
            }
        }
        
        result.success = true;
        
    } catch (const py::error_already_set& e) {
        result.error_message = std::string("Python error: ") + e.what();
    } catch (const std::exception& e) {
        result.error_message = std::string("C++ error: ") + e.what();
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    result.compute_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    return result;
}

std::map<std::string, double> PythonBridge::fetch_parameters(const std::string& system_name) {
    std::lock_guard<std::mutex> lock(python_mutex_);
    
    std::map<std::string, double> params;
    
    if (!initialized_ || !impl_) {
        return params;
    }
    
    try {
        if (!impl_->apifetch_module.is_none() && py::hasattr(impl_->apifetch_module, "fetch")) {
            py::dict py_result = impl_->apifetch_module.attr("fetch")(system_name).cast<py::dict>();
            
            for (auto& item : py_result) {
                std::string key = item.first.cast<std::string>();
                try {
                    double value = item.second.cast<double>();
                    params[key] = value;
                } catch (...) {
                    // Skip non-numeric values
                }
            }
        }
    } catch (const py::error_already_set& e) {
        std::cerr << "[PythonBridge] fetch_parameters error: " << e.what() << std::endl;
    }
    
    return params;
}

std::vector<std::string> PythonBridge::get_available_systems() {
    std::lock_guard<std::mutex> lock(python_mutex_);
    
    std::vector<std::string> systems;
    
    // Read from bodies.csv
    std::ifstream file("bodies.csv");
    if (file.is_open()) {
        std::string line;
        std::getline(file, line);  // Skip header
        while (std::getline(file, line)) {
            size_t comma = line.find(',');
            if (comma != std::string::npos) {
                systems.push_back(line.substr(0, comma));
            }
        }
    }
    
    return systems;
}

std::string PythonBridge::execute_python(const std::string& code) {
    std::lock_guard<std::mutex> lock(python_mutex_);
    
    if (!initialized_) {
        return "ERROR: Python not initialized";
    }
    
    try {
        py::exec(code);
        return "OK";
    } catch (const py::error_already_set& e) {
        return std::string("ERROR: ") + e.what();
    }
}

PythonBridge::ModuleVersions PythonBridge::get_module_versions() {
    std::lock_guard<std::mutex> lock(python_mutex_);
    
    ModuleVersions versions;
    
    if (!initialized_) {
        return versions;
    }
    
    try {
        // Try to get __version__ from each module
        if (impl_ && !impl_->condensed_physics_module.is_none()) {
            if (py::hasattr(impl_->condensed_physics_module, "__version__")) {
                versions.condensed_physics = impl_->condensed_physics_module.attr("__version__").cast<std::string>();
            } else {
                versions.condensed_physics = "available";
            }
        }
        
        if (impl_ && !impl_->qcalc_module.is_none()) {
            if (py::hasattr(impl_->qcalc_module, "__version__")) {
                versions.qcalc = impl_->qcalc_module.attr("__version__").cast<std::string>();
            } else {
                versions.qcalc = "available";
            }
        }
        
        // NumPy
        try {
            py::module_ np = py::module_::import("numpy");
            versions.numpy = np.attr("__version__").cast<std::string>();
        } catch (...) {}
        
        // SciPy
        try {
            py::module_ sp = py::module_::import("scipy");
            versions.scipy = sp.attr("__version__").cast<std::string>();
        } catch (...) {}
        
        // SymPy
        try {
            py::module_ sym = py::module_::import("sympy");
            versions.sympy = sym.attr("__version__").cast<std::string>();
        } catch (...) {}
        
    } catch (...) {}
    
    return versions;
}

#else // PYTHON_BRIDGE_ENABLED not defined

// Stub implementation when Python bridge is disabled

struct PythonBridge::Impl {};

PythonBridge& PythonBridge::instance() {
    static PythonBridge instance;
    return instance;
}

PythonBridge::~PythonBridge() {}

bool PythonBridge::initialize() {
    std::cout << "[PythonBridge] Python bridge disabled (pybind11 not found)" << std::endl;
    return false;
}

void PythonBridge::shutdown() {}

PythonFieldResult PythonBridge::calculate_condensed_physics(const PythonQueryParams&) {
    PythonFieldResult result;
    result.error_message = "Python bridge disabled";
    return result;
}

PythonFieldResult PythonBridge::calculate_qcalc(const PythonQueryParams&) {
    PythonFieldResult result;
    result.error_message = "Python bridge disabled";
    return result;
}

std::map<std::string, double> PythonBridge::fetch_parameters(const std::string&) {
    return {};
}

std::vector<std::string> PythonBridge::get_available_systems() {
    return {};
}

std::string PythonBridge::execute_python(const std::string&) {
    return "ERROR: Python bridge disabled";
}

PythonBridge::ModuleVersions PythonBridge::get_module_versions() {
    return {};
}

#endif // PYTHON_BRIDGE_ENABLED

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

bool check_python_dependencies() {
    // Check if required Python files exist
    return std::filesystem::exists("CondensedPhysics.py") &&
           std::filesystem::exists("QCalc.py") &&
           std::filesystem::exists("IPData.py") &&
           std::filesystem::exists("OPData.py");
}

std::string get_python_executable() {
#ifdef _WIN32
    if (std::filesystem::exists(".venv/Scripts/python.exe")) {
        return ".venv/Scripts/python.exe";
    }
    return "python.exe";
#else
    if (std::filesystem::exists(".venv/bin/python")) {
        return ".venv/bin/python";
    }
    return "python3";
#endif
}

} // namespace UQFF
