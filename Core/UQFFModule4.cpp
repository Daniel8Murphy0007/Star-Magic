/**
 * ================================================================================================
 * UQFFModule4.cpp - Implementation of Self-Expanding UQFF Physics Module
 * ================================================================================================
 * 
 * Extracted from: source4.cpp (lines 336-792)
 * Purpose: Implementation of adaptive physics engine methods
 * Phase: 1, Week 2 - Module Extraction
 * 
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#include "UQFFModule4.hpp"
#include <array> // MSVC requirement

namespace UQFFCore {

UQFFModule4::UQFFModule4()
    : enableDynamicTerms(false),
      enableLogging(false),
      learningRate(0.001),
      updateCounter(0)
{
    // Initialize metadata
    metadata["version"] = "2.0-Enhanced";
    metadata["module"] = "UQFFModule4";
    metadata["created"] = "November 08, 2025";
    metadata["framework"] = "Self-Expanding UQFF";

    // Initialize default variables
    variables["mass"] = 1e30;
    variables["radius"] = 1e6;
    variables["temperature"] = 1e6;
    variables["magnetic_field"] = 1e-5;
}

// ========================================================================
// DYNAMIC TERM MANAGEMENT
// ========================================================================

void UQFFModule4::registerDynamicTerm(std::unique_ptr<PhysicsTerm> term)
{
    if (enableLogging)
    {
        std::cout << "[UQFFModule4] Registering dynamic term: "
                  << term->getName() << std::endl;
    }
    dynamicTerms.push_back(std::move(term));
}

double UQFFModule4::computeDynamicTerms(double t) const
{
    if (!enableDynamicTerms)
        return 0.0;

    double total = 0.0;
    std::map<std::string, double> params;
    for (const auto &kv : variables)
    {
        params[kv.first] = kv.second;
    }
    for (const auto &kv : dynamicParameters)
    {
        params[kv.first] = kv.second;
    }

    for (const auto &term : dynamicTerms)
    {
        if (term->validate(params))
        {
            total += term->compute(t, params);
        }
    }
    return total;
}

// ========================================================================
// VARIABLE MANAGEMENT WITH HISTORY
// ========================================================================

void UQFFModule4::updateVariable(const std::string &name, double value)
{
    variables[name] = value;
    variable_history[name].push_back(value);

    // Keep history limited to last 1000 entries
    if (variable_history[name].size() > 1000)
    {
        variable_history[name].erase(variable_history[name].begin());
    }

    updateCounter++;

    if (enableLogging)
    {
        std::cout << "[UQFFModule4] Updated " << name << " = " << value << std::endl;
    }
}

double UQFFModule4::getVariable(const std::string &name) const
{
    auto it = variables.find(name);
    return (it != variables.end()) ? it->second : 0.0;
}

void UQFFModule4::addCustomVariable(const std::string &name, double value,
                       const std::string &dependency)
{
    variables[name] = value;
    if (!dependency.empty())
    {
        variable_dependencies[name] = dependency;
    }
    if (enableLogging)
    {
        std::cout << "[UQFFModule4] Added custom variable: " << name
                  << " = " << value << std::endl;
    }
}

std::vector<double> UQFFModule4::getVariableHistory(const std::string &name, int steps) const
{
    auto it = variable_history.find(name);
    if (it == variable_history.end())
        return {};

    if (steps < 0 || static_cast<size_t>(steps) >= it->second.size())
    {
        return it->second;
    }

    return std::vector<double>(
        it->second.end() - steps,
        it->second.end());
}

// ========================================================================
// DYNAMIC PARAMETER MANAGEMENT
// ========================================================================

void UQFFModule4::setDynamicParameter(const std::string &name, double value)
{
    dynamicParameters[name] = value;
    if (enableLogging)
    {
        std::cout << "[UQFFModule4] Set dynamic parameter: " << name
                  << " = " << value << std::endl;
    }
}

double UQFFModule4::getDynamicParameter(const std::string &name) const
{
    auto it = dynamicParameters.find(name);
    return (it != dynamicParameters.end()) ? it->second : 0.0;
}

// ========================================================================
// AUTO-CALIBRATION
// ========================================================================

void UQFFModule4::addTunableParameter(const std::string &name)
{
    tunableParams.push_back(name);
}

bool UQFFModule4::autoCalibrate(const std::string &observable, double targetValue,
                   double tolerance, int maxIterations)
{
    if (enableLogging)
    {
        std::cout << "[UQFFModule4] Auto-calibrating " << observable
                  << " to target: " << targetValue << std::endl;
    }

    for (int iter = 0; iter < maxIterations; ++iter)
    {
        double currentValue = getVariable(observable);
        double error = targetValue - currentValue;

        if (std::abs(error / targetValue) < tolerance)
        {
            if (enableLogging)
            {
                std::cout << "[UQFFModule4] Calibration converged in "
                          << iter << " iterations" << std::endl;
            }
            return true;
        }

        // Adjust tunable parameters using gradient descent
        for (const auto &param : tunableParams)
        {
            double currentParam = getVariable(param);
            double gradient = computeGradient(param, observable);
            double adjustment = learningRate * error / (gradient + 1e-10);
            updateVariable(param, currentParam + adjustment);
        }
    }

    if (enableLogging)
    {
        std::cout << "[UQFFModule4] Calibration did not converge" << std::endl;
    }
    return false;
}

double UQFFModule4::computeGradient(const std::string &param, const std::string &observable)
{
    double epsilon = 1e-6;
    double originalValue = getVariable(param);
    double originalObservable = getVariable(observable);

    updateVariable(param, originalValue + epsilon);
    double perturbedObservable = getVariable(observable);
    updateVariable(param, originalValue);

    return (perturbedObservable - originalObservable) / epsilon;
}

// ========================================================================
// ADAPTIVE UPDATES AND SELF-LEARNING
// ========================================================================

void UQFFModule4::setLearningRate(double rate)
{
    learningRate = rate;
    if (enableLogging)
    {
        std::cout << "[UQFFModule4] Learning rate set to " << rate << std::endl;
    }
}

double UQFFModule4::getLearningRate() const
{
    return learningRate;
}

void UQFFModule4::adaptiveUpdate()
{
    if (!enableDynamicTerms)
        return;

    // Evolution timescale (example: 8e14 seconds)
    double evolution_timescale = 8e14;
    double dt = 1.0; // Default timestep
    double evolution_factor = std::exp(-dt / evolution_timescale);

    // Update key variables with adaptive evolution
    for (auto &kv : variables)
    {
        const std::string &varName = kv.first;
        double &varValue = kv.second;

        // Apply evolution factor
        varValue *= evolution_factor;

        // Record in history
        variable_history[varName].push_back(varValue);
    }

    updateCounter++;

    if (enableLogging && updateCounter % 5 == 0)
    {
        std::cout << "[UQFFModule4] Adaptive update #" << updateCounter << std::endl;
    }
}

void UQFFModule4::scaleToObservations(const std::map<std::string, double> &obsData)
{
    if (enableLogging)
    {
        std::cout << "[UQFFModule4] Scaling to observational data..." << std::endl;
    }

    for (const auto &kv : obsData)
    {
        const std::string &obsName = kv.first;
        double obsValue = kv.second;

        auto it = variables.find(obsName);
        if (it != variables.end())
        {
            double currentValue = it->second;
            double scaleFactor = obsValue / (currentValue + 1e-30);

            // Apply scaling to variable and related dependencies
            updateVariable(obsName, obsValue);

            if (enableLogging)
            {
                std::cout << "  Scaled " << obsName << " by factor "
                          << scaleFactor << std::endl;
            }
        }
    }
}

// ========================================================================
// STATE PERSISTENCE
// ========================================================================

void UQFFModule4::exportState(const std::string &filename) const
{
    std::ofstream out(filename);
    if (!out.is_open())
    {
        std::cerr << "[UQFFModule4] Failed to open " << filename << std::endl;
        return;
    }

    out << "# UQFFModule4 State Export\n";
    out << "# Generated: November 08, 2025\n\n";

    out << "[Metadata]\n";
    for (const auto &kv : metadata)
    {
        out << kv.first << "=" << kv.second << "\n";
    }

    out << "\n[Variables]\n";
    for (const auto &kv : variables)
    {
        out << kv.first << "=" << kv.second << "\n";
    }

    out << "\n[DynamicParameters]\n";
    for (const auto &kv : dynamicParameters)
    {
        out << kv.first << "=" << kv.second << "\n";
    }

    out << "\n[Configuration]\n";
    out << "enableDynamicTerms=" << enableDynamicTerms << "\n";
    out << "enableLogging=" << enableLogging << "\n";
    out << "learningRate=" << learningRate << "\n";
    out << "updateCounter=" << updateCounter << "\n";

    out.close();

    if (enableLogging)
    {
        std::cout << "[UQFFModule4] State exported to " << filename << std::endl;
    }
}

void UQFFModule4::importState(const std::string &filename)
{
    std::ifstream in(filename);
    if (!in.is_open())
    {
        std::cerr << "[UQFFModule4] Failed to open " << filename << std::endl;
        return;
    }

    std::string line, section;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#')
            continue;

        if (line[0] == '[')
        {
            section = line;
            continue;
        }

        size_t pos = line.find('=');
        if (pos == std::string::npos)
            continue;

        std::string key = line.substr(0, pos);
        std::string value = line.substr(pos + 1);

        if (section == "[Variables]")
        {
            variables[key] = std::stod(value);
        }
        else if (section == "[DynamicParameters]")
        {
            dynamicParameters[key] = std::stod(value);
        }
        else if (section == "[Metadata]")
        {
            metadata[key] = value;
        }
        else if (section == "[Configuration]")
        {
            if (key == "enableDynamicTerms")
                enableDynamicTerms = (value == "1");
            else if (key == "enableLogging")
                enableLogging = (value == "1");
            else if (key == "learningRate")
                learningRate = std::stod(value);
            else if (key == "updateCounter")
                updateCounter = std::stoi(value);
        }
    }

    in.close();

    if (enableLogging)
    {
        std::cout << "[UQFFModule4] State imported from " << filename << std::endl;
    }
}

// ========================================================================
// CONFIGURATION
// ========================================================================

void UQFFModule4::setEnableLogging(bool enable)
{
    enableLogging = enable;
}

void UQFFModule4::setEnableDynamicTerms(bool enable)
{
    enableDynamicTerms = enable;
}

std::map<std::string, std::string> UQFFModule4::getMetadata() const
{
    return metadata;
}

int UQFFModule4::getUpdateCounter() const
{
    return updateCounter;
}

} // namespace UQFFCore
