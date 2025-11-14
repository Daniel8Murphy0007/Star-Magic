/**
 * ================================================================================================
 * FluidSolver.hpp - Navier-Stokes Fluid Simulation for UQFF Jet Dynamics
 * ================================================================================================
 * 
 * Extracted from: source4.cpp (lines 792-975)
 * Purpose: 2D incompressible fluid solver for quasar jet and UQFF dynamics
 * Phase: 1, Week 2 - Module Extraction
 * 
 * Contains:
 * - FluidSolver class: Jos Stam's "Stable Fluids" method
 * - Navier-Stokes solver: diffuse, advect, project, boundary conditions
 * - UQFF integration: step() with gravity-like force
 * - Jet simulation: add_jet_force() for SCm expulsion
 * - Visualization: print_velocity_field()
 * 
 * Grid Configuration:
 * - N = 32 (grid size, configurable)
 * - dt = 0.1 (time step)
 * - visc = 0.0001 (viscosity)
 * 
 * Dependencies:
 * - C++17 standard library
 * - <vector> for grid storage
 * - <cmath> for sqrt
 * - <iostream> for visualization
 * 
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef FLUIDSOLVER_HPP
#define FLUIDSOLVER_HPP

#include <vector>
#include <cmath>
#include <iostream>

namespace UQFFCore {

// Fluid solver configuration constants
constexpr int GRID_SIZE = 32;              // Grid size (N+2) x (N+2)
constexpr double TIME_STEP = 0.1;          // Time step for integration
constexpr double VISCOSITY = 0.0001;       // Fluid viscosity
constexpr double JET_FORCE = 10.0;         // Default jet force magnitude

// Macro for 2D grid indexing: IX(i, j) = i + (N+2)*j
#define IX(i, j) ((i) + (GRID_SIZE + 2) * (j))

/**
 * FluidSolver: 2D incompressible Navier-Stokes solver
 * 
 * Implements Jos Stam's "Stable Fluids" method with:
 * - Diffusion (viscosity)
 * - Advection (velocity transport)
 * - Projection (incompressibility)
 * - UQFF gravity integration
 * 
 * Grid: (N+2) x (N+2) with boundary cells
 * Velocities: u (horizontal), v (vertical)
 * Density: dens (optional scalar field)
 */
class FluidSolver
{
public:
    // Velocity fields (current and previous timestep)
    std::vector<double> u, v, u_prev, v_prev;
    
    // Density fields (optional, for visualization)
    std::vector<double> dens, dens_prev;

    /**
     * Constructor
     * Initializes all grids to zero
     */
    FluidSolver();

    /**
     * Add source term to field
     * @param x Target field
     * @param s Source field
     */
    void add_source(std::vector<double>& x, std::vector<double>& s);

    /**
     * Diffusion step (viscosity)
     * @param b Boundary condition type (1=horizontal, 2=vertical, 0=none)
     * @param x Output field
     * @param x0 Input field
     * @param diff Diffusion coefficient (e.g., viscosity)
     */
    void diffuse(int b, std::vector<double>& x, std::vector<double>& x0, double diff);

    /**
     * Advection step (velocity transport)
     * @param b Boundary condition type
     * @param d Output field
     * @param d0 Input field
     */
    void advect(int b, std::vector<double>& d, std::vector<double>& d0);

    /**
     * Projection step (enforce incompressibility)
     * @param u Horizontal velocity
     * @param v Vertical velocity
     * @param p Pressure field (temporary)
     * @param div Divergence field (temporary)
     */
    void project(std::vector<double>& u, std::vector<double>& v, 
                 std::vector<double>& p, std::vector<double>& div);

    /**
     * Set boundary conditions
     * @param b Boundary type (1=horizontal flip, 2=vertical flip, 0=copy)
     * @param x Field to apply boundaries
     */
    void set_bnd(int b, std::vector<double>& x);

    /**
     * Advance simulation one timestep
     * Integrates UQFF gravity-like force
     * @param uqff_g UQFF gravity acceleration (m/s²)
     */
    void step(double uqff_g = 0.0);

    /**
     * Add jet force in center (simulates SCm expulsion)
     * @param force Jet force magnitude
     */
    void add_jet_force(double force);

    /**
     * Print velocity field magnitude as ASCII art
     * Symbols: # (>1.0), + (>0.5), . (>0.1), space (<0.1)
     */
    void print_velocity_field();

private:
    const int N = GRID_SIZE;
    const double dt_ns = TIME_STEP;
    const double visc = VISCOSITY;
};

} // namespace UQFFCore

#endif // FLUIDSOLVER_HPP
