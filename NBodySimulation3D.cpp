// STANDALONE_NBODYSIMULATION3D
#include "NBodySimulation3D.h"
#include <iostream>
#include <numeric>

void NBodySimulation3D::add_particle(double m, double x, double y, double z,
                         double vx, double vy, double vz) {
    particles.push_back({m, x, y, z, vx, vy, vz, 0, 0, 0});
}

void NBodySimulation3D::compute_accelerations() {
    for (auto& p : particles) { p.ax = p.ay = p.az = 0.0; }
    int N = (int)particles.size();
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            double dz = particles[j].z - particles[i].z;
            double r2 = dx*dx + dy*dy + dz*dz + softening*softening;
            double r3 = r2 * std::sqrt(r2);
            double f_ij = G / r3;
            particles[i].ax += f_ij * particles[j].mass * dx;
            particles[i].ay += f_ij * particles[j].mass * dy;
            particles[i].az += f_ij * particles[j].mass * dz;
            particles[j].ax -= f_ij * particles[i].mass * dx;
            particles[j].ay -= f_ij * particles[i].mass * dy;
            particles[j].az -= f_ij * particles[i].mass * dz;
        }
    }
}

void NBodySimulation3D::step_euler() {
    compute_accelerations();
    for (auto& p : particles) {
        p.x += p.vx * dt; p.y += p.vy * dt; p.z += p.vz * dt;
        p.vx += p.ax * dt; p.vy += p.ay * dt; p.vz += p.az * dt;
    }
    curr_t += dt;
}

void NBodySimulation3D::step_leapfrog() {
    // Half-kick
    for (auto& p : particles) {
        p.vx += 0.5*p.ax*dt; p.vy += 0.5*p.ay*dt; p.vz += 0.5*p.az*dt;
    }
    // Drift
    for (auto& p : particles) {
        p.x += p.vx*dt; p.y += p.vy*dt; p.z += p.vz*dt;
    }
    compute_accelerations();
    // Half-kick again
    for (auto& p : particles) {
        p.vx += 0.5*p.ax*dt; p.vy += 0.5*p.ay*dt; p.vz += 0.5*p.az*dt;
    }
    curr_t += dt;
}

double NBodySimulation3D::total_energy() const {
    double KE = 0, PE = 0;
    for (const auto& p : particles)
        KE += 0.5 * p.mass * (p.vx*p.vx + p.vy*p.vy + p.vz*p.vz);
    int N = (int)particles.size();
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++) {
            double dx = particles[j].x-particles[i].x, dy = particles[j].y-particles[i].y,
                   dz = particles[j].z-particles[i].z;
            double r = std::sqrt(dx*dx+dy*dy+dz*dz+softening*softening);
            PE -= G * particles[i].mass * particles[j].mass / r;
        }
    return KE + PE;
}

std::string NBodySimulation3D::primary_equation() const {
    return "a_i = sum_j G*m_j*(r_j-r_i)/(|r_ij|^2+eps^2)^(3/2); r+=v*dt; v+=a*dt (Euler) / KDK (Leapfrog)";
}

void NBodySimulation3D::self_update() { step_leapfrog(); }
void NBodySimulation3D::self_expand() { dt *= 1.01; }
void NBodySimulation3D::simulate(int num_steps) {
    for (int i = 0; i < num_steps; i++) {
        step_leapfrog();
        std::cout << "t=" << curr_t/3.156e7 << " yr  E=" << total_energy() << " J\n";
    }
}

#ifdef STANDALONE_NBODYSIMULATION3D
int main() {
    NBodySimulation3D sim;
    // 2-body Sun-Jupiter analogue
    sim.add_particle(1.989e30, 0, 0, 0, 0, 0, 0);
    sim.add_particle(1.898e27, 5.2*1.496e11, 0, 0, 0, 13.1e3, 0);
    std::cout << "3D N-Body Simulation\n";
    std::cout << sim.primary_equation() << "\n";
    sim.simulate(5);
    return 0;
}
#endif
