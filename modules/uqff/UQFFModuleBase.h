// UQFFModuleBase.h
// Abstract base class for all UQFF modules. Provides the map-based variable
// interface and pure virtual compute/text API.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef UQFF_MODULE_BASE_H
#define UQFF_MODULE_BASE_H

#include <map>
#include <string>
#include <iostream>

class UQFFModuleBase {
protected:
    std::map<std::string, double> variables;

    void _warnMissing(const std::string& name) const {
        std::cerr << "[UQFF] Variable '" << name << "' not found.\n";
    }

public:
    virtual ~UQFFModuleBase() = default;

    // --- Dynamic variable interface ---

    void updateVariable(const std::string& name, double value) {
        variables[name] = value;
        onVariableUpdated(name, value);
    }

    void addToVariable(const std::string& name, double delta) {
        auto it = variables.find(name);
        if (it != variables.end()) {
            it->second += delta;
        } else {
            _warnMissing(name);
            variables[name] = delta;
        }
    }

    void subtractFromVariable(const std::string& name, double delta) {
        addToVariable(name, -delta);
    }

    double getVariable(const std::string& name) const {
        auto it = variables.find(name);
        if (it != variables.end()) return it->second;
        _warnMissing(name);
        return 0.0;
    }

    bool hasVariable(const std::string& name) const {
        return variables.count(name) > 0;
    }

    // --- Core API (implement in each module) ---

    // Compute effective gravity / field value at time t
    virtual double computeG(double t) = 0;

    // Return human-readable equation text
    virtual std::string getEquationText() const = 0;

    // Print all variables to stdout
    virtual void printVariables() const {
        std::cout << "Variables:\n";
        for (const auto& kv : variables) {
            std::cout << "  " << kv.first << " = " << kv.second << "\n";
        }
    }

protected:
    // Hook called after updateVariable — subclasses override for derived deps
    virtual void onVariableUpdated(const std::string& /*name*/, double /*value*/) {}
};

#endif // UQFF_MODULE_BASE_H
