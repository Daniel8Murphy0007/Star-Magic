// STANDALONE_NGC2525BARREDSPIRAL2UQFF
#include "NGC2525BarredSpiral2UQFF.h"
#include <iostream>

double NGC2525BarredSpiral2UQFF::M_tot() const {
    return M_stellar + M_DM;
}

double NGC2525BarredSpiral2UQFF::H_z() const {
    return H_0 * std::sqrt(0.3 * std::pow(1.0 + z_ngc, 3) + 0.7);
}

// SN Ia feedback: energy impulse decaying over ~10 yr
// a_SN = L_SN / (c * M_tot * r_gal) * exp(-t/tau_SN)
double NGC2525BarredSpiral2UQFF::a_SN(double t) const {
    double a_peak = L_SN / (c * M_tot() * r_gal);
    return a_peak * std::exp(-t / tau_SN);
}

// EM bar streaming with aether correction
double NGC2525BarredSpiral2UQFF::a_EM_bar() const {
    const double m_p = 1.673e-27;
    double a_bare = q_p * v_bar * B_bar / m_p;
    return a_bare * (1.0 + rho_UA / rho_SCm) * 1.0e-12;
}

// Full MUGE v2
double NGC2525BarredSpiral2UQFF::g_NGC2525v2(double t) const {
    double g_grav = G * M_tot() / (r_gal * r_gal);
    double g_core = g_grav * (1.0 + H_z() * t) * (1.0 + f_TRZ);
    return g_core + a_SN(t) + a_EM_bar();
}

std::string NGC2525BarredSpiral2UQFF::primary_equation() const {
    return "g_NGC2525v2(t)=G*(M_s+M_DM)/r^2*(1+H_z*t)*(1+f_TRZ)+L_SN/(c*M*r)*exp(-t/tau_SN)+q*v_bar*B/(m_p)*(1+rho_UA/rho_SCm)*1e-12";
}

void NGC2525BarredSpiral2UQFF::self_update() {
    curr_t += time_step;
}

void NGC2525BarredSpiral2UQFF::self_expand() {
    r_gal *= 1.0001;
    M_DM  *= 1.00005;
}

void NGC2525BarredSpiral2UQFF::simulate(int num_steps) {
    for (int step = 0; step < num_steps; step++) {
        double t = curr_t + step * time_step;
        double g = g_NGC2525v2(t);
        std::cout << "Step " << step << "  t=" << t / 3.156e13
                  << " Myr  g=" << g << " m/s^2\n";
    }
}

#ifdef STANDALONE_NGC2525BARREDSPIRAL2UQFF
int main() {
    NGC2525BarredSpiral2UQFF n;
    std::cout << "NGC 2525 Barred Spiral v2 UQFF\n";
    std::cout << n.primary_equation() << "\n\n";
    std::cout << "g(0) = " << n.g_NGC2525v2(0.0) << " m/s^2\n";
    n.simulate(3);
    return 0;
}
#endif
