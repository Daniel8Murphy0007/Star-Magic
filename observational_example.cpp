// observational_example.cpp
// Example usage of observational systems config with SOURCE10 buoyancy physics
// Demonstrates how to use existing MAIN_1_CoAnQi.cpp physics terms with different systems
// Copyright - Daniel T. Murphy, 2025

#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>
#include "observational_systems_config.h"
#include <array> // MSVC requirement

// Example: Simulate using SOURCE10 BuoyancyUQFF term
// (In real usage, you'd link against MAIN_1_CoAnQi.cpp and use actual PhysicsTerm classes)

class BuoyancyUQFF_Example
{
public:
    double compute(double t, const std::map<std::string, double> &params)
    {
        // Simplified buoyancy calculation (real one is in SOURCE10)
        double beta_i = 0.6;
        double V_infl = 1e-6;
        double rho_vac = 1e-30;
        double a_universal = 1e12;

        // System-specific enhancements
        double M = params.count("M") ? params.at("M") : 1e30;
        double r = params.count("r") ? params.at("r") : 1e10;
        double L_X = params.count("L_X") ? params.at("L_X") : 1e30;
        double B0 = params.count("B0") ? params.at("B0") : 1e-9;
        double G = params.count("G") ? params.at("G") : 6.6743e-11;

        // Base buoyancy
        double base_buoyancy = beta_i * V_infl * rho_vac * a_universal;

        // Gravitational enhancement
        double grav_enhancement = G * M / (r * r);

        // Observational scaling
        double obs_scaling = std::sqrt(L_X / 1e30) * std::sqrt(B0 / 1e-9);

        return base_buoyancy * grav_enhancement * obs_scaling;
    }
};

void printSystemInfo(const std::string &name)
{
    const ObservationalSystem *sys = getSystem(name);
    if (!sys)
    {
        std::cout << "System '" << name << "' not found!" << std::endl;
        return;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "SYSTEM: " << sys->name << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Description: " << sys->description << std::endl;
    std::cout << "Category:    " << sys->category << std::endl;
    std::cout << "Telescope:   " << sys->telescope << std::endl;
    std::cout << "\nPhysical Parameters:" << std::endl;
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  Mass (M):              " << sys->M << " kg" << std::endl;
    std::cout << "  Radius (r):            " << sys->r << " m" << std::endl;
    std::cout << "  X-ray Luminosity:      " << sys->L_X << " W" << std::endl;
    std::cout << "  Magnetic Field (B0):   " << sys->B0 << " T" << std::endl;
    std::cout << "  Gas Density:           " << sys->rho_gas << " kg/m³" << std::endl;
    std::cout << "  Gas Temperature:       " << sys->T_gas << " K" << std::endl;
    std::cout << "  Angular Frequency:     " << sys->omega0 << " rad/s" << std::endl;
    std::cout << "  Age/Timescale:         " << sys->t_age << " s" << std::endl;
}

void computeBuoyancyForSystem(const std::string &system_name, double t)
{
    // Get system parameters
    auto params = systemToParams(system_name);
    if (params.empty())
    {
        std::cout << "System '" << system_name << "' not found!" << std::endl;
        return;
    }

    // Compute buoyancy using example term
    BuoyancyUQFF_Example buoyancy;
    double result = buoyancy.compute(t, params);

    std::cout << "\nBuoyancy Force Calculation:" << std::endl;
    std::cout << "  System: " << system_name << std::endl;
    std::cout << "  Time:   " << std::scientific << t << " s" << std::endl;
    std::cout << "  F_buoy: " << std::scientific << result << " N" << std::endl;
}

void compareSystemsByCategory(const std::string &category)
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "SYSTEMS IN CATEGORY: " << category << std::endl;
    std::cout << "========================================" << std::endl;

    auto systems = getSystemsByCategory(category);
    if (systems.empty())
    {
        std::cout << "No systems found in category '" << category << "'" << std::endl;
        return;
    }

    BuoyancyUQFF_Example buoyancy;
    double t = 1e12; // Example time

    std::cout << std::left << std::setw(20) << "System"
              << std::setw(15) << "Mass (kg)"
              << std::setw(15) << "L_X (W)"
              << std::setw(15) << "Buoyancy (N)" << std::endl;
    std::cout << std::string(65, '-') << std::endl;

    for (const auto &sys_name : systems)
    {
        auto params = systemToParams(sys_name);
        double F = buoyancy.compute(t, params);

        std::cout << std::left << std::setw(20) << sys_name
                  << std::scientific << std::setprecision(2)
                  << std::setw(15) << params["M"]
                  << std::setw(15) << params["L_X"]
                  << std::setw(15) << F << std::endl;
    }
}

int main()
{
    std::cout << "========================================" << std::endl;
    std::cout << "OBSERVATIONAL SYSTEMS CONFIGURATION" << std::endl;
    std::cout << "Using SOURCE10 Buoyancy Physics" << std::endl;
    std::cout << "========================================" << std::endl;

    // Example 1: List all available systems
    std::cout << "\nAvailable Systems (" << OBSERVATIONAL_SYSTEMS.size() << " total):" << std::endl;
    auto all_systems = listSystems();
    int count = 0;
    for (const auto &name : all_systems)
    {
        std::cout << "  " << std::setw(18) << std::left << name;
        if (++count % 3 == 0)
            std::cout << std::endl;
    }
    std::cout << std::endl;

    // Example 2: Show detailed info for specific systems
    printSystemInfo("Vela");
    printSystemInfo("ESO137");
    printSystemInfo("ASASSN14li");

    // Example 3: Compute buoyancy for different systems
    std::cout << "\n========================================" << std::endl;
    std::cout << "BUOYANCY CALCULATIONS (t = 1e12 s)" << std::endl;
    std::cout << "========================================" << std::endl;

    computeBuoyancyForSystem("Vela", 1e12);
    computeBuoyancyForSystem("NGC1365", 1e12);
    computeBuoyancyForSystem("ElGordo", 1e12);

    // Example 4: Compare systems by category
    compareSystemsByCategory("pulsar");
    compareSystemsByCategory("galaxy_cluster");
    compareSystemsByCategory("agn");

    // Example 5: Time evolution for a single system
    std::cout << "\n========================================" << std::endl;
    std::cout << "TIME EVOLUTION: Vela Pulsar" << std::endl;
    std::cout << "========================================" << std::endl;

    BuoyancyUQFF_Example buoyancy;
    auto vela_params = systemToParams("Vela");

    std::cout << std::scientific << std::setprecision(3);
    std::cout << std::left << std::setw(15) << "Time (s)"
              << std::setw(20) << "Buoyancy (N)" << std::endl;
    std::cout << std::string(35, '-') << std::endl;

    for (int i = 0; i <= 10; i++)
    {
        double t = 1e11 * std::pow(10, i / 10.0);
        double F = buoyancy.compute(t, vela_params);
        std::cout << std::setw(15) << t << std::setw(20) << F << std::endl;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "USAGE NOTE:" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "To use with actual MAIN_1_CoAnQi.cpp physics terms:" << std::endl;
    std::cout << "1. Include observational_systems_config.h" << std::endl;
    std::cout << "2. Get params: auto params = systemToParams(\"Vela\");" << std::endl;
    std::cout << "3. Use with SOURCE10 terms:" << std::endl;
    std::cout << "   BuoyancyUQFF buoy;" << std::endl;
    std::cout << "   double F = buoy.compute(t, params);" << std::endl;
    std::cout << "\nAll 35+ observational systems use the same" << std::endl;
    std::cout << "core physics (SOURCE10), just different parameters!" << std::endl;

    return 0;
}
