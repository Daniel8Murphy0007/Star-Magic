#ifndef FORNAX_CLUSTER_GRAVITATIONAL_H
#define FORNAX_CLUSTER_GRAVITATIONAL_H
// STANDALONE_FORNAXCLUSTERGRAVITATIONAL
// PAPER_690: Fornax Cluster Gravitational N-Body Dynamics
// 58-member galaxy cluster at 19.5 Mpc, dominated by NGC 1399 (cD galaxy)
// Source: grok_share_ba508f76c8e.txt | Session 174
#include <vector>
#include <string>
#include <cmath>

struct FornaxGalaxy {
    std::string name;
    double mass;       // kg
    double x, y, z;   // position m
    double vx, vy, vz;// velocity m/s
};

class FornaxClusterGravitational {
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

    // Cluster parameters (Fornax Cluster)
    double M_cluster  = 7.0e13 * M_sun;  // total cluster mass
    double R_cluster  = 1.0e3  * kpc;    // virial radius ~1 Mpc
    double distance   = 19.5   * Mpc;    // distance to Earth
    int    N_galaxies = 58;              // confirmed members

    // UQFF cluster correction
    // g_cluster_UQFF = G*M_cluster/r^2 * (1 + rho_SCm/rho_UA) * (1 + f_TRZ)

    std::vector<FornaxGalaxy> galaxies;
    double time_step = 3.156e14;  // 10 Myr steps
    mutable double curr_t = 0.0;

    FornaxClusterGravitational();  // populates member galaxies

    // Gravitational force between two galaxies
    void   compute_forces();
    void   update_positions(double dt);

    // UQFF-modified cluster gravitational acceleration
    double g_cluster_UQFF(double r) const;

    // Velocity dispersion (virial theorem)
    // sigma_v^2 = G * M_cluster / (2 * R_cluster)
    double velocity_dispersion() const;

    // Tidal radius for member galaxy
    double tidal_radius(double m_gal, double r_orbit) const;

    std::string primary_equation() const;

    // Self-expand/self-update/self-simulate
    virtual void self_update();
    virtual void self_expand();
    virtual void simulate(int num_steps = 100);

};

#endif // FORNAX_CLUSTER_GRAVITATIONAL_H
