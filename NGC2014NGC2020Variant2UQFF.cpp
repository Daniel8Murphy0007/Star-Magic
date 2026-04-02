// STANDALONE_NGC2014NGC2020VARIANT2UQFF
#include "NGC2014NGC2020Variant2UQFF.h"
#include <iostream>

// WR mass decreases as wind ejects material
double NGC2014NGC2020Variant2UQFF::M_WR_t(double t) const {
    return M_WR * (1.0 + M_0_WR * std::exp(-t / tau_SF_WR));
}

// WR radiation pressure on cone gas
// a_rad = L_WR / (4*pi*r_cone^2 * c) * (rho_WR/m_H)
double NGC2014NGC2020Variant2UQFF::a_WR_radiation() const {
    const double m_H = 1.67e-27;
    double flux = L_WR / (4.0 * M_PI * r_cone * r_cone * c);
    return flux * (rho_WR / m_H);
}

double NGC2014NGC2020Variant2UQFF::a_EM_WR() const {
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_gas * B_WR / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}

// WR cone MUGE variant 2
double NGC2014NGC2020Variant2UQFF::g_WRcone(double t) const {
    double M  = M_WR_t(t);
    double g_n = G * M / (r_cone * r_cone);
    double g_core = g_n * (1.0 + H_0 * t) * (1.0 + f_TRZ);
    return g_core + a_WR_radiation() + a_EM_WR();
}

std::string NGC2014NGC2020Variant2UQFF::primary_equation() const {
    return "g_WRcone(t)=G*M_WR(t)/r^2*(1+H0*t)*(1+f_TRZ)+L_WR/(4pi*r^2*c)*(rho_WR/m_H)+q*v*B/m*(1+rho_UA/rho_SCm)*1e-12";
}

void NGC2014NGC2020Variant2UQFF::self_update() {
    curr_t += time_step;
}

void NGC2014NGC2020Variant2UQFF::self_expand() {
    r_cone  *= 1.002;
    rho_WR  *= 0.998;
}

void NGC2014NGC2020Variant2UQFF::simulate(int num_steps) {
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        double g = g_WRcone(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\n";
    }
}

#ifdef STANDALONE_NGC2014NGC2020VARIANT2UQFF
int main() {
    NGC2014NGC2020Variant2UQFF wrc;
    std::cout << "NGC 2014/2020 WR-variant 2 UQFF\n";
    std::cout << wrc.primary_equation() << "\n\n";
    std::cout << "g(2.5 Myr) = " << wrc.g_WRcone(7.89e13) << " m/s^2\n";
    wrc.simulate(3);
    return 0;
}
#endif
