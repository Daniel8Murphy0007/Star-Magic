#ifndef N_BODY_SIMULATION_3D_H
#define N_BODY_SIMULATION_3D_H
// STANDALONE_NBODYSIMULATION3D
// PAPER_691: 3D N-Body Gravitational Simulation Framework
// Euler integrator for stellar/galactic scale particle dynamics
// Source: grok_share_ba508f76c8e.txt (#108) | Session 174
#include <vector>
#include <string>
#include <cmath>

struct Particle3D {
    double mass;
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;  // acceleration (computed each step)
};

class NBodySimulation3D {
public:

    // UQFF Universal Constants
    static constexpr double G        = 6.6743e-11;  // m^3 kg^-1 s^-2
    static constexpr double c        = 3.0e8;       // m/s
    static constexpr double hbar     = 1.0546e-34;  // J*s
    static constexpr double mu_0     = 1.2566e-6;   // H/m
    static constexpr double k_B      = 1.3806e-23;  // J/K
    static constexpr double M_sun    = 1.989e30;    // kg
    static constexpr double kpc      = 3.086e19;    // m (1 kpc)
    static constexpr double Mpc      = 3.086e22;    // m (1 Mpc)
    static constexpr double rho_UA   = 7.09e-36;    // J/m^3 Universal Aether
    static constexpr double rho_SCm  = 7.09e-37;    // J/m^3 Schwarzschild condensate
    static constexpr double f_TRZ    = 0.1;         // time-reversal zone factor
    static constexpr double kappa    = 5.0e-4;      // UQFF calibration day^-1
    static constexpr double SSq      = 0.57;        // superstring quenching
    static constexpr double mu_J     = 3.38e23;     // J*m magnetic string coupling
    static constexpr double Lambda   = 1.1e-52;     // m^-2 cosmological constant
    static constexpr double H_0      = 2.269e-18;   // s^-1 Hubble constant
    static constexpr double t_H      = 4.35e17;     // s Hubble time

    std::vector<Particle3D> particles;
    double dt         = 3.156e10;   // time step ~1000 yr
    double softening  = 1.0 * kpc; // gravitational softening
    mutable double curr_t = 0.0;

    NBodySimulation3D() = default;

    void add_particle(double m, double x, double y, double z,
                      double vx, double vy, double vz);

    // Compute pairwise gravitational accelerations
    // a_i = sum_j G*m_j*(r_j - r_i) / (|r_j-r_i|^2 + eps^2)^(3/2)
    void compute_accelerations();

    // Euler step: r(t+dt) = r(t) + v(t)*dt, v(t+dt) = v(t) + a(t)*dt
    void step_euler();

    // Leapfrog (kick-drift-kick) for better energy conservation
    void step_leapfrog();

    // Total kinetic + potential energy
    double total_energy() const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // N_BODY_SIMULATION_3D_H
