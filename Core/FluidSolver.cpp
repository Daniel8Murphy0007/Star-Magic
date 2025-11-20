/**
 * ================================================================================================
 * FluidSolver.cpp - Implementation of Navier-Stokes Fluid Solver
 * ================================================================================================
 * 
 * Extracted from: source4.cpp (lines 802-975)
 * Purpose: Implementation of 2D incompressible fluid dynamics
 * Phase: 1, Week 2 - Module Extraction
 * 
 * Based on Jos Stam's "Stable Fluids" method (SIGGRAPH 1999)
 * 
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#include "FluidSolver.hpp"
#include <array> // MSVC requirement

namespace UQFFCore {

FluidSolver::FluidSolver()
{
    int size = (N + 2) * (N + 2);
    u.resize(size, 0.0);
    v.resize(size, 0.0);
    u_prev.resize(size, 0.0);
    v_prev.resize(size, 0.0);
    dens.resize(size, 0.0);
    dens_prev.resize(size, 0.0);
}

void FluidSolver::add_source(std::vector<double>& x, std::vector<double>& s)
{
    for (size_t i = 0; i < x.size(); ++i)
    {
        x[i] += dt_ns * s[i];
    }
}

void FluidSolver::diffuse(int b, std::vector<double>& x, std::vector<double>& x0, double diff)
{
    double a = dt_ns * diff * N * N;
    for (int k = 0; k < 20; ++k)
    {
        for (int i = 1; i <= N; ++i)
        {
            for (int j = 1; j <= N; ++j)
            {
                x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                                   x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
                              (1 + 4 * a);
            }
        }
        set_bnd(b, x);
    }
}

void FluidSolver::advect(int b, std::vector<double>& d, std::vector<double>& d0)
{
    int i0, j0, i1, j1;
    double x, y, s0, t0, s1, t1;
    for (int i = 1; i <= N; ++i)
    {
        for (int j = 1; j <= N; ++j)
        {
            x = i - dt_ns * N * u[IX(i, j)];
            y = j - dt_ns * N * v[IX(i, j)];
            if (x < 0.5)
                x = 0.5;
            if (x > N + 0.5)
                x = N + 0.5;
            if (y < 0.5)
                y = 0.5;
            if (y > N + 0.5)
                y = N + 0.5;
            i0 = (int)x;
            i1 = i0 + 1;
            j0 = (int)y;
            j1 = j0 + 1;
            s1 = x - i0;
            s0 = 1 - s1;
            t1 = y - j0;
            t0 = 1 - t1;
            d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                          s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }
    set_bnd(b, d);
}

void FluidSolver::project(std::vector<double>& u, std::vector<double>& v, 
                          std::vector<double>& p, std::vector<double>& div)
{
    double h = 1.0 / N;
    for (int i = 1; i <= N; ++i)
    {
        for (int j = 1; j <= N; ++j)
        {
            div[IX(i, j)] = -0.5 * h * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + 
                                         v[IX(i, j + 1)] - v[IX(i, j - 1)]);
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);
    
    for (int k = 0; k < 20; ++k)
    {
        for (int i = 1; i <= N; ++i)
        {
            for (int j = 1; j <= N; ++j)
            {
                p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] +
                               p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4;
            }
        }
        set_bnd(0, p);
    }
    
    for (int i = 1; i <= N; ++i)
    {
        for (int j = 1; j <= N; ++j)
        {
            u[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
            v[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
        }
    }
    set_bnd(1, u);
    set_bnd(2, v);
}

void FluidSolver::set_bnd(int b, std::vector<double>& x)
{
    for (int i = 1; i <= N; ++i)
    {
        x[IX(0, i)] = (b == 1) ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N + 1, i)] = (b == 1) ? -x[IX(N, i)] : x[IX(N, i)];
        x[IX(i, 0)] = (b == 2) ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N + 1)] = (b == 2) ? -x[IX(i, N)] : x[IX(i, N)];
    }
    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N + 1)] = 0.5 * (x[IX(1, N + 1)] + x[IX(0, N)]);
    x[IX(N + 1, 0)] = 0.5 * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
    x[IX(N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

void FluidSolver::step(double uqff_g)
{
    // Add UQFF gravity-like force as body force in v (vertical direction)
    for (int i = 1; i <= N; ++i)
    {
        for (int j = 1; j <= N; ++j)
        {
            v[IX(i, j)] += dt_ns * uqff_g; // Integrate UQFF acceleration into velocity
        }
    }

    diffuse(1, u_prev, u, visc);
    diffuse(2, v_prev, v, visc);
    project(u_prev, v_prev, u, v);
    advect(1, u, u_prev);
    advect(2, v, v_prev);
    project(u, v, u_prev, v_prev);
}

void FluidSolver::add_jet_force(double force)
{
    // Add force in the center as a jet (simulating SCm expulsion)
    for (int i = N / 4; i <= 3 * N / 4; ++i)
    {
        v[IX(i, N / 2)] += force;
    }
}

void FluidSolver::print_velocity_field()
{
    std::cout << "Velocity field (magnitude):" << std::endl;
    for (int j = N; j >= 1; --j)
    { // Print top to bottom
        for (int i = 1; i <= N; ++i)
        {
            double mag = std::sqrt(u[IX(i, j)] * u[IX(i, j)] + v[IX(i, j)] * v[IX(i, j)]);
            char sym = (mag > 1.0) ? '#' : (mag > 0.5) ? '+' : (mag > 0.1) ? '.' : ' ';
            std::cout << sym;
        }
        std::cout << std::endl;
    }
}

} // namespace UQFFCore
