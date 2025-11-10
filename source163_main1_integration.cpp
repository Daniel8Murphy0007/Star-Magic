// source163_main1_integration.cpp
// Integration layer connecting Source163 (AstroSystemsUQFFModule) with MAIN_1.cpp
// Allows MAIN_1.cpp to use multi-system UQFF calculations from Source163
// Copyright - Daniel T. Murphy, Nov 10, 2025
// NOTE: Compile Source163.cpp separately and link, or use wrapper functions

#include <iostream>
#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <complex>

using cdouble = std::complex<double>;

// NOTE: This is a demonstration - in practice, create Source163.o first:
// g++ -c -std=c++17 Source163.cpp -o Source163.o
// Then: g++ source163_main1_integration.cpp Source163.o -o integration
// For this demo, we'll just document the interface

// Integration wrapper class
class Source163Integration
{
private:
    AstroSystemsUQFFModule module;

public:
    Source163Integration()
    {
        module.setEnableLogging(false); // Quiet mode for integration
    }

    // Compute F_U_Bi_i for a specific system using Source163
    double computeForce(const std::string &system, double t, double M, double r)
    {
        // Override system params if provided
        if (M > 0)
            module.updateVariable("M", cdouble(M, 0.0));
        if (r > 0)
            module.updateVariable("r", cdouble(r, 0.0));

        cdouble result = module.computeMasterEquations(system, t);
        return result.real(); // Return real part for MAIN_1.cpp
    }

    // Compute compressed gravity g(r,t)
    double computeGravity(const std::string &system, double t)
    {
        return module.computeCompressedG(system, t);
    }

    // Get all supported systems
    std::vector<std::string> getSystems()
    {
        return module.getSupportedSystems();
    }

    // Batch compute all systems
    std::map<std::string, double> computeAllSystems(double t)
    {
        auto results = module.computeAllSystems(t);
        std::map<std::string, double> real_results;
        for (const auto &pair : results)
        {
            real_results[pair.first] = pair.second.real();
        }
        return real_results;
    }

    // Time evolution for a system
    std::vector<double> evolveSystem(const std::string &system, double t_start, double t_end, double dt)
    {
        auto complex_results = module.simulateTimeEvolution(system, t_start, t_end, dt);
        std::vector<double> real_results;
        for (const auto &val : complex_results)
        {
            real_results.push_back(val.real());
        }
        return real_results;
    }

    // DPM evolution
    double computeDPMEvolution(const std::string &system, double t, double dt)
    {
        return module.computeDPMEvolution(system, t, dt).real();
    }

    // SMBH accretion
    double computeSMBHAccretion(const std::string &system, double t)
    {
        return module.computeSMBHAccretion(system, t).real();
    }

    // TDE dynamics
    double computeTDEDynamics(const std::string &system, double t)
    {
        return module.computeTDEDynamics(system, t).real();
    }

    // Get equation text
    std::string getEquation(const std::string &system)
    {
        return module.getEquationText(system);
    }

    // Export/import state
    void exportState(const std::string &filename)
    {
        module.exportState(filename);
    }

    void importState(const std::string &filename)
    {
        module.importState(filename);
    }

    // Update any variable
    void updateParameter(const std::string &name, double value)
    {
        module.updateVariable(name, cdouble(value, 0.0));
    }
};

// ============================================================================
// DEMONSTRATION: Using Source163 with MAIN_1.cpp-style interface
// ============================================================================

int main()
{
    Source163Integration integration;

    std::cout << "=== Source163 Integration with MAIN_1.cpp ===" << std::endl;
    std::cout << "Multi-system UQFF calculations\n"
              << std::endl;

    // Get supported systems
    auto systems = integration.getSystems();
    std::cout << "Supported systems: ";
    for (const auto &sys : systems)
    {
        std::cout << sys << " ";
    }
    std::cout << "\n"
              << std::endl;

    // Compute for NGC685 (similar to MAIN_1.cpp workflow)
    std::string system = "NGC685";
    double t = 1e12;
    double M = 1e41; // kg
    double r = 1e21; // m

    double F = integration.computeForce(system, t, M, r);
    double g = integration.computeGravity(system, t);

    std::cout << "System: " << system << std::endl;
    std::cout << "M = " << M << " kg" << std::endl;
    std::cout << "r = " << r << " m" << std::endl;
    std::cout << "t = " << t << " s" << std::endl;
    std::cout << "F_U_Bi_i = " << std::scientific << F << " N" << std::endl;
    std::cout << "g(r,t) = " << g << " m/s^2" << std::endl;
    std::cout << std::endl;

    // Batch computation (all systems at once)
    std::cout << "=== Batch Computation (all systems) ===" << std::endl;
    auto all_results = integration.computeAllSystems(t);
    for (const auto &pair : all_results)
    {
        std::cout << pair.first << ": F = " << std::scientific << pair.second << " N" << std::endl;
    }
    std::cout << std::endl;

    // Time evolution (like MAIN_1.cpp simulations)
    std::cout << "=== Time Evolution (NGC685) ===" << std::endl;
    auto evolution = integration.evolveSystem("NGC685", 0, 1e13, 2e12);
    std::cout << "Time steps: " << evolution.size() << std::endl;
    std::cout << "F(t=0): " << evolution[0] << " N" << std::endl;
    std::cout << "F(t=end): " << evolution.back() << " N" << std::endl;
    std::cout << std::endl;

    // DPM evolution
    std::cout << "=== DPM Evolution ===" << std::endl;
    double dF_dt = integration.computeDPMEvolution("NGC685", t, 1e11);
    std::cout << "dF/dt = " << dF_dt << " N/s" << std::endl;
    std::cout << std::endl;

    // SMBH accretion
    std::cout << "=== SMBH Accretion ===" << std::endl;
    double L_acc = integration.computeSMBHAccretion("NGC685", t);
    std::cout << "L_accretion = " << L_acc << " W" << std::endl;
    std::cout << std::endl;

    // TDE dynamics (AT2024tvd only)
    std::cout << "=== Tidal Disruption Event ===" << std::endl;
    double L_tde = integration.computeTDEDynamics("AT2024tvd", 1.5e6);
    std::cout << "L_TDE = " << L_tde << " W" << std::endl;
    std::cout << std::endl;

    // Get equation
    std::cout << "=== Equation ===" << std::endl;
    std::cout << integration.getEquation("NGC685") << std::endl;
    std::cout << std::endl;

    // Export state
    integration.exportState("integration_state.txt");
    std::cout << "State exported to integration_state.txt" << std::endl;

    std::cout << "\n=== Integration Complete ===" << std::endl;
    std::cout << "Source163 can now be used by MAIN_1.cpp" << std::endl;

    return 0;
}
