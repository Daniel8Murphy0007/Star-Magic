// STANDALONE_SOMBREROGALAXYM104NGC4594
#include "SombreroGalaxyM104NGC4594.h"
#include <iostream>

// sigma^2 = G*M_bulge / (5 * R_bulge) [virial estimate]
double SombreroGalaxyM104NGC4594::sigma_bulge() const {
    return std::sqrt(G * M_bulge / (5.0 * R_bulge));
}

// v_c(R) = sqrt(G*(M_bulge+M_disk)/R) with UQFF oscillation correction
double SombreroGalaxyM104NGC4594::v_circular(double R) const {
    double M_enc = M_bulge + M_disk * (R / R_disk < 1.0 ? R / R_disk : 1.0);
    double v0    = std::sqrt(G * M_enc / R);
    return v0 * (1.0 + 0.5 * rho_SCm / rho_UA);
}

// Phi_UQFF: bulge (Hernquist) + disk (log) + cosmological
// Phi ~ -G*M_bulge/(r+a) - G*M_disk/(2*R_d)*log(R^2+h^2) + Lambda*c^2*r^2/6
double SombreroGalaxyM104NGC4594::Phi_UQFF(double r) const {
    double a     = R_bulge;
    double Phi_b = -G * M_bulge / (r + a);
    double R     = r;  // edge-on projection r ~= R
    double Phi_d = -G * M_disk / (2.0 * R_disk) * std::log(R * R + h_disk * h_disk);
    double Phi_L = Lambda * c * c * r * r / 6.0;
    return Phi_b + Phi_d + Phi_L;
}

// Dust lane wavefunction: psi = exp(-z^2/(2*h^2)) * cos(omega_i*t)
double SombreroGalaxyM104NGC4594::psi_dust_lane(double z, double t) const {
    return std::exp(-z * z / (2.0 * h_disk * h_disk)) * std::cos(5.0e-18 * t);
}

std::string SombreroGalaxyM104NGC4594::primary_equation() const {
    return "v_c = sqrt(G*M_enc/R)*(1+rho_SCm/rho_UA/2); Phi = -G*M_b/(r+a) - G*M_d*log(R^2+h^2)/(2*R_d) + Lambda*c^2*r^2/6";
}

void SombreroGalaxyM104NGC4594::self_update() { curr_t += time_step; }
void SombreroGalaxyM104NGC4594::self_expand() { R_disk *= 1.0001; }
void SombreroGalaxyM104NGC4594::simulate(int num_steps) {
    for (int i = 0; i < num_steps; i++) {
        self_update();
        std::cout << "t=" << curr_t/3.156e7 << " yr  sigma=" << sigma_bulge()/1e3
                  << " km/s  v_c(R_disk)=" << v_circular(R_disk)/1e3 << " km/s\n";
    }
}

#ifdef STANDALONE_SOMBREROGALAXYM104NGC4594
int main() {
    SombreroGalaxyM104NGC4594 sg;
    std::cout << "Sombrero Galaxy M104/NGC4594\n";
    std::cout << sg.primary_equation() << "\n";
    std::cout << "sigma_bulge = " << sg.sigma_bulge()/1e3 << " km/s\n";
    std::cout << "v_circ(10 kpc) = " << sg.v_circular(10e3*sg.kpc)/1e3 << " km/s\n";
    sg.simulate(3);
    return 0;
}
#endif
