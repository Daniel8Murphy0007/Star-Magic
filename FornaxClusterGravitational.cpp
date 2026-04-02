// STANDALONE_FORNAXCLUSTERGRAVITATIONAL
#include "FornaxClusterGravitational.h"
#include <iostream>

FornaxClusterGravitational::FornaxClusterGravitational() {
    // Populate key Fornax members
    galaxies = {
        {"NGC 1399", 2.0e12*M_sun, 0, 0, 0, 0, 0, 0},  // cD galaxy
        {"NGC 1316", 5.0e11*M_sun, 1.0e3*kpc, 0, 0, 50e3, 0, 0},
        {"NGC 1404", 4.0e11*M_sun,-0.5e3*kpc, 0, 0,-30e3, 0, 0},
        {"NGC 1365", 3.0e11*M_sun, 2.0e3*kpc, 1.5e3*kpc, 0, 20e3, 10e3, 0},
        {"NGC 1340", 2.0e11*M_sun,-1.0e3*kpc,-0.5e3*kpc, 0,-20e3, 5e3, 0},
        {"NGC 1380", 1.5e11*M_sun, 0.5e3*kpc,-1.0e3*kpc, 0, 10e3,-15e3, 0}
    };
}

void FornaxClusterGravitational::compute_forces() {
    int N = (int)galaxies.size();
    for (int i = 0; i < N; i++) {
        double fx = 0, fy = 0, fz = 0;
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            double dx = galaxies[j].x - galaxies[i].x;
            double dy = galaxies[j].y - galaxies[i].y;
            double dz = galaxies[j].z - galaxies[i].z;
            double r3 = std::pow(dx*dx+dy*dy+dz*dz, 1.5) + 1e30;
            double f  = G * galaxies[i].mass * galaxies[j].mass / r3;
            fx += f*dx; fy += f*dy; fz += f*dz;
        }
        galaxies[i].vx += fx/galaxies[i].mass * time_step;
        galaxies[i].vy += fy/galaxies[i].mass * time_step;
        galaxies[i].vz += fz/galaxies[i].mass * time_step;
    }
}

void FornaxClusterGravitational::update_positions(double dt) {
    for (auto& g : galaxies) {
        g.x += g.vx * dt;
        g.y += g.vy * dt;
        g.z += g.vz * dt;
    }
}

// UQFF-modified cluster gravity:
// g = G*M/r^2 * (1 + rho_SCm/rho_UA) * (1 + f_TRZ)
double FornaxClusterGravitational::g_cluster_UQFF(double r) const {
    double g0 = G * M_cluster / (r * r);
    return g0 * (1.0 + rho_SCm / rho_UA) * (1.0 + f_TRZ);
}

// sigma_v^2 = G*M_cluster / (2*R_cluster) [virial theorem]
double FornaxClusterGravitational::velocity_dispersion() const {
    return std::sqrt(G * M_cluster / (2.0 * R_cluster));
}

// r_tidal = r_orbit * (m_gal / (3*M_cluster))^(1/3)
double FornaxClusterGravitational::tidal_radius(double m_gal, double r_orbit) const {
    return r_orbit * std::cbrt(m_gal / (3.0 * M_cluster));
}

std::string FornaxClusterGravitational::primary_equation() const {
    return "g_cluster_UQFF = G*M_cluster/r^2 * (1+rho_SCm/rho_UA) * (1+f_TRZ); sigma_v = sqrt(G*M/(2*R))";
}

void FornaxClusterGravitational::self_update() { curr_t += time_step; compute_forces(); update_positions(time_step); }
void FornaxClusterGravitational::self_expand() { R_cluster *= 1.001; }
void FornaxClusterGravitational::simulate(int num_steps) {
    for (int i = 0; i < num_steps; i++) {
        self_update();
        std::cout << "t=" << curr_t/3.156e7 << " yr  galaxies=" << galaxies.size()
                  << "  sigma_v=" << velocity_dispersion()/1e3 << " km/s\n";
    }
}

#ifdef STANDALONE_FORNAXCLUSTERGRAVITATIONAL
int main() {
    FornaxClusterGravitational fc;
    std::cout << "Fornax Cluster Gravitational Simulation\n";
    std::cout << fc.primary_equation() << "\n";
    std::cout << "sigma_v = " << fc.velocity_dispersion()/1e3 << " km/s\n";
    std::cout << "r_tidal(NGC1316) = " << fc.tidal_radius(5e11*fc.M_sun,1e3*fc.kpc)/fc.kpc << " kpc\n";
    fc.simulate(3);
    return 0;
}
#endif
