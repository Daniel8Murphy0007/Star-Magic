#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include <vector>

class FluidSolver
{
public:
    std::vector<double> u, v, u_prev, v_prev, dens, dens_prev;
    int size;
    double dt, visc, diff;

    FluidSolver();
    void add_source(std::vector<double> &x, std::vector<double> &s);
    void diffuse(int b, std::vector<double> &x, std::vector<double> &x0, double diff);
    void advect(int b, std::vector<double> &d, std::vector<double> &d0);
    void project(std::vector<double> &u, std::vector<double> &v, std::vector<double> &p, std::vector<double> &div);
    void set_bnd(int b, std::vector<double> &x);
    void step(double uqff_g = 0.0);
    void add_jet_force(double force);
    void print_velocity_field();
};

#endif // FLUID_SOLVER_H
