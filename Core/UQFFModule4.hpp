/**
 * ================================================================================================
 * UQFFModule4.hpp - Self-Expanding UQFF Physics Module
 * ================================================================================================
 * 
 * Extracted from: source4.cpp (lines 336-792)
 * Purpose: Adaptive physics engine with runtime term registration and auto-calibration
 * Phase: 1, Week 2 - Module Extraction
 * 
 * Contains:
 * - Variable management with history tracking
 * - Dynamic physics term registration system
 * - Auto-calibration via gradient descent
 * - Adaptive parameter tuning
 * - State persistence (export/import)
 * - Observational data scaling
 * 
 * Features:
 * - Register new physics terms at runtime (no recompilation)
 * - Track variable changes over time (1000-step history)
 * - Auto-calibrate to observational targets
 * - Export/import module state for reproducibility
 * - Learning rate adjustment for optimization
 * 
 * Dependencies:
 * - PhysicsTerms.hpp (for dynamic term framework)
 * - C++17 standard library
 * - <map>, <vector>, <string> for data structures
 * - <fstream> for state persistence
 * 
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef UQFFMODULE4_HPP
#define UQFFMODULE4_HPP

#include "PhysicsTerms.hpp"
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

namespace UQFFCore {

/**
 * UQFFModule4: Self-expanding physics computation engine
 * 
 * Enhanced framework (v2.0) with:
 * - Runtime physics term registration
 * - Variable history tracking (last 1000 updates)
 * - Auto-calibration to observational data
 * - State export/import for reproducibility
 * - Adaptive learning system
 */
class UQFFModule4
{
private:
    // Core variables storage with history tracking
    std::map<std::string, double> variables;
    std::map<std::string, std::vector<double>> variable_history;
    std::map<std::string, std::string> variable_dependencies;

    // Dynamic term system
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, double> dynamicParameters;

    // Metadata tracking
    std::map<std::string, std::string> metadata;

    // Configuration flags
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;

    // Tunable parameters for auto-calibration
    std::vector<std::string> tunableParams;

    // Update tracking
    int updateCounter;

public:
    /**
     * Constructor
     * Initializes module with default configuration and metadata
     */
    UQFFModule4();

    // ========================================================================
    // DYNAMIC TERM MANAGEMENT
    // ========================================================================

    /**
     * Register a new physics term for runtime computation
     * @param term Unique pointer to PhysicsTerm implementation
     */
    void registerDynamicTerm(std::unique_ptr<PhysicsTerm> term);

    /**
     * Compute sum of all registered dynamic terms
     * @param t Time (s)
     * @return Total contribution from all dynamic terms
     */
    double computeDynamicTerms(double t) const;

    // ========================================================================
    // VARIABLE MANAGEMENT WITH HISTORY
    // ========================================================================

    /**
     * Update variable value and record in history
     * @param name Variable name
     * @param value New value
     */
    void updateVariable(const std::string& name, double value);

    /**
     * Get current variable value
     * @param name Variable name
     * @return Current value (0.0 if not found)
     */
    double getVariable(const std::string& name) const;

    /**
     * Add custom variable with optional dependency tracking
     * @param name Variable name
     * @param value Initial value
     * @param dependency Optional dependency description
     */
    void addCustomVariable(const std::string& name, double value,
                           const std::string& dependency = "");

    /**
     * Get variable history for analysis
     * @param name Variable name
     * @param steps Number of steps to retrieve (-1 for all)
     * @return Vector of historical values
     */
    std::vector<double> getVariableHistory(const std::string& name, int steps = -1) const;

    // ========================================================================
    // DYNAMIC PARAMETER MANAGEMENT
    // ========================================================================

    /**
     * Set dynamic parameter for physics terms
     * @param name Parameter name
     * @param value Parameter value
     */
    void setDynamicParameter(const std::string& name, double value);

    /**
     * Get dynamic parameter value
     * @param name Parameter name
     * @return Parameter value (0.0 if not found)
     */
    double getDynamicParameter(const std::string& name) const;

    // ========================================================================
    // AUTO-CALIBRATION
    // ========================================================================

    /**
     * Add parameter to auto-calibration tunable set
     * @param name Parameter name
     */
    void addTunableParameter(const std::string& name);

    /**
     * Auto-calibrate parameters to match observational target
     * Uses gradient descent optimization
     * @param observable Observable variable name
     * @param targetValue Target value for observable
     * @param tolerance Relative tolerance for convergence
     * @param maxIterations Maximum optimization iterations
     * @return true if converged, false otherwise
     */
    bool autoCalibrate(const std::string& observable, double targetValue,
                       double tolerance = 0.01, int maxIterations = 100);

    /**
     * Compute gradient of observable with respect to parameter
     * @param param Parameter name
     * @param observable Observable name
     * @return Numerical gradient
     */
    double computeGradient(const std::string& param, const std::string& observable);

    // ========================================================================
    // ADAPTIVE UPDATES AND SELF-LEARNING
    // ========================================================================

    /**
     * Set learning rate for gradient descent
     * @param rate Learning rate (typically 0.001 - 0.1)
     */
    void setLearningRate(double rate);

    /**
     * Get current learning rate
     * @return Learning rate
     */
    double getLearningRate() const;

    /**
     * Perform adaptive update based on recent history
     * Adjusts learning rate based on convergence behavior
     */
    void adaptiveUpdate();

    /**
     * Scale module variables to observational data
     * @param obsData Map of observable names to measured values
     */
    void scaleToObservations(const std::map<std::string, double>& obsData);

    // ========================================================================
    // STATE PERSISTENCE
    // ========================================================================

    /**
     * Export module state to file
     * Format: INI-style with sections [Metadata], [Variables], [DynamicParameters], [Configuration]
     * @param filename Output filename
     */
    void exportState(const std::string& filename) const;

    /**
     * Import module state from file
     * @param filename Input filename
     */
    void importState(const std::string& filename);

    // ========================================================================
    // CONFIGURATION
    // ========================================================================

    /**
     * Enable/disable logging
     * @param enable true to enable, false to disable
     */
    void setEnableLogging(bool enable);

    /**
     * Enable/disable dynamic term computation
     * @param enable true to enable, false to disable
     */
    void setEnableDynamicTerms(bool enable);

    /**
     * Get module metadata
     * @return Map of metadata key-value pairs
     */
    std::map<std::string, std::string> getMetadata() const;

    /**
     * Get update counter (total number of variable updates)
     * @return Update count
     */
    int getUpdateCounter() const;
};

} // namespace UQFFCore

#endif // UQFFMODULE4_HPP
