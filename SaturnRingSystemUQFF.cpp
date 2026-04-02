// STANDALONE_SATURNRINGSYSTEMUQFF
#include "SaturnRingSystemUQFF.h"
#include <iostream>
#include <fstream>

// Solar gravitational pull at Saturn orbit, with Hubble + TRZ correction
// g_solar = G*M_Sun/r_orbit^2 * (1 + H_0*t) * (1 + f_TRZ)
double SaturnRingSystemUQFF::g_solar(double t) const {
    double g_base = G * M_Sun / (r_orbit * r_orbit);
    return g_base * (1.0 + H_0 * t) * (1.0 + f_TRZ);
}

// Saturn self-gravity at equatorial surface
// g_self = G*M_Saturn/r_eq^2  => ~10.44 m/s^2
double SaturnRingSystemUQFF::g_self() const {
    return G * M_Saturn / (r_eq * r_eq);
}

// Ring tidal acceleration from main ring system
// T_ring = G*M_ring/r_ring^2 ~ 2.043e-7 m/s^2
double SaturnRingSystemUQFF::T_ring() const {
    return G * M_ring / (r_ring * r_ring);
}

// Wind dynamic pressure converted to acceleration
// a_wind ~ rho_atm * v_wind^2 scaled ~ 2.5e-7 m/s^2
double SaturnRingSystemUQFF::a_wind() const {
    double P_wind = rho_atm * v_wind * v_wind;
    // Scale by 1/M_Saturn * V_eff (effective column volume ~ 1e18 m^3)
    double V_eff = 1.0e18;
    return P_wind * V_eff / M_Saturn;
}

// Electromagnetic aether term with [UA]/[SCm] correction
// a_EM = q_p*(v x B)*(1 + rho_UA/rho_SCm)*1e-12
//      = q_p*v_charged*B_Sat / m_proton * 11 * 1e-12
double SaturnRingSystemUQFF::a_EM() const {
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_charged * B_Sat / m_p;
    double correction = 1.0 + rho_UA / rho_SCm;  // = 11
    return a_bare * correction * 1.0e-12;
}

// Full Saturn MUGE
// g_Saturn(t) = g_solar(t) + g_self + T_ring + a_wind + a_EM
// Numerical example (t=0): ~ 10.44 m/s^2 (surface gravity)
double SaturnRingSystemUQFF::g_Saturn(double t) const {
    return g_solar(t) + g_self() + T_ring() + a_wind() + a_EM();
}

std::string SaturnRingSystemUQFF::primary_equation() const {
    return "g_Saturn(t) = G*M_Sun/r_orbit^2*(1+H0*t)*(1+f_TRZ) + G*M_Sat/r_eq^2 + G*M_ring/r_ring^2 + rho*v_wind^2/M_eff + q*(vxB)*(1+rho_UA/rho_SCm)*1e-12";
}

void SaturnRingSystemUQFF::self_update() {
    curr_t += time_step;
    r_orbit *= 1.0 + 1.0e-10 * time_step; // very slow orbital drift
}

void SaturnRingSystemUQFF::self_expand() {
    r_ring *= 1.001;
    M_ring *= 0.9999; // ring mass slowly dissipates
}

void SaturnRingSystemUQFF::simulate(int num_steps) {
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        double g = g_Saturn(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e7 << " yr"
                  << "  g_Saturn=" << g << " m/s^2\n";
    }
}

#ifdef STANDALONE_SATURNRINGSYSTEMUQFF
int main() {
    SaturnRingSystemUQFF sat;
    std::cout << "Saturn Ring System UQFF Simulation\n";
    std::cout << sat.primary_equation() << "\n\n";
    std::cout << "g_Saturn(t=0) = " << sat.g_Saturn(0.0) << " m/s^2\n";
    sat.simulate(5);
    return 0;
}
#endif
