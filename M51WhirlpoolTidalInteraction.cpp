// STANDALONE_M51WHIRLPOOLTIDALINTERACTION
#include "M51WhirlpoolTidalInteraction.h"
#include <iostream>

double M51WhirlpoolTidalInteraction::tidal_force(double R_d) const {
    return 2.0 * G * M_comp * R_d / std::pow(d_sep, 3);
}

double M51WhirlpoolTidalInteraction::pitch_angle() const {
    // pitch ~ arctan(F_tidal / v_circ^2 / R) approximated
    double v_circ = std::sqrt(G * M_main / R_spiral);
    double F_t    = tidal_force(R_spiral);
    return std::atan(F_t / (v_circ * v_circ / R_spiral));
}

double M51WhirlpoolTidalInteraction::g_M51_UQFF(double r) const {
    double M_tot = M_main + M_comp;
    double H_z   = H_0 * std::sqrt(0.3 * std::pow(1.0 + z_M51, 3) + 0.7);
    double g0    = G * M_tot / (r * r);
    double g_ism = rho_ISM * 1.0e51 * g0;  // ISM fluid correction
    return g0 * (1.0 + H_z) + g_ism * (rho_SCm / rho_UA);
}

// SFE_UQFF: star formation efficiency enhanced by tidal interaction
// SFE = SFR * (1 + F_tidal/g_grav) * (1 + f_TRZ)
double M51WhirlpoolTidalInteraction::SFE_UQFF() const {
    double g_grav = G * M_main / (R_spiral * R_spiral);
    double F_t    = tidal_force(R_spiral);
    return SFR * (1.0 + F_t / g_grav) * (1.0 + f_TRZ);
}

std::string M51WhirlpoolTidalInteraction::primary_equation() const {
    return "g_M51 = G*(M_main+M_comp)/r^2*(1+H(z)) + rho_ISM*V*g*(rho_SCm/rho_UA); F_tidal=2*G*M_comp*R/d^3";
}

void M51WhirlpoolTidalInteraction::self_update() { curr_t += time_step; d_sep *= 0.99999; }
void M51WhirlpoolTidalInteraction::self_expand() { R_spiral *= 1.0001; }
void M51WhirlpoolTidalInteraction::simulate(int num_steps) {
    for (int i = 0; i < num_steps; i++) {
        self_update();
        std::cout << "t=" << curr_t/3.156e7 << " yr  pitch=" << pitch_angle()*180/M_PI
                  << " deg  SFE=" << SFE_UQFF() << " M_sun/yr\n";
    }
}

#ifdef STANDALONE_M51WHIRLPOOLTIDALINTERACTION
int main() {
    M51WhirlpoolTidalInteraction m51;
    std::cout << "M51 Whirlpool Tidal Interaction\n";
    std::cout << m51.primary_equation() << "\n";
    std::cout << "Tidal force at spiral arm: " << m51.tidal_force(m51.R_spiral) << " N/kg\n";
    std::cout << "Pitch angle: " << m51.pitch_angle()*180/M_PI << " deg\n";
    m51.simulate(3);
    return 0;
}
#endif
