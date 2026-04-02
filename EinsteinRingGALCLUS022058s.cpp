// STANDALONE_EINSTEINRINGGALCLUS022058S
#include "EinsteinRingGALCLUS022058s.h"
#include <iostream>

// R_E = sqrt(4*G*M*D_LS / (c^2 * D_L * D_S)) * D_L
double EinsteinRingGALCLUS022058s::einstein_radius() const {
    return D_L * std::sqrt(4.0 * G * M_lens * D_LS / (c * c * D_L * D_S));
}

// u = theta/theta_E (dimensionless impact parameter)
double EinsteinRingGALCLUS022058s::magnification(double theta_arcsec) const {
    double theta_E = theta_E_arcsec * (M_PI / (180.0 * 3600.0));  // rad
    double theta   = theta_arcsec   * (M_PI / (180.0 * 3600.0));
    double u       = theta / theta_E;
    return (u * u + 2.0) / (u * std::sqrt(u * u + 4.0));
}

// alpha_UQFF = 4*G*M/(c^2*r) * (1 + rho_SCm/rho_UA)
double EinsteinRingGALCLUS022058s::deflection_angle_UQFF(double r) const {
    double alpha_GR = 4.0 * G * M_lens / (c * c * r);
    return alpha_GR * (1.0 + rho_SCm / rho_UA);
}

// Sigma_crit = c^2*D_S / (4*pi*G*D_L*D_LS)
double EinsteinRingGALCLUS022058s::convergence(double Sigma) const {
    double Sigma_crit = c * c * D_S / (4.0 * M_PI * G * D_L * D_LS);
    return Sigma / Sigma_crit;
}

std::string EinsteinRingGALCLUS022058s::primary_equation() const {
    return "R_E=sqrt(4*G*M*D_LS/(c^2*D_L*D_S))*D_L; alpha_UQFF=4*G*M/(c^2*r)*(1+rho_SCm/rho_UA)";
}

void EinsteinRingGALCLUS022058s::self_update() { curr_t += time_step; }
void EinsteinRingGALCLUS022058s::self_expand() { theta_E_arcsec *= 1.0001; }
void EinsteinRingGALCLUS022058s::simulate(int num_steps) {
    for (int i = 0; i < num_steps; i++) {
        self_update();
        std::cout << "R_E = " << einstein_radius()/kpc << " kpc  mu(10") = "
                  << magnification(10.0) << "\n";
    }
}

#ifdef STANDALONE_EINSTEINRINGGALCLUS022058S
int main() {
    EinsteinRingGALCLUS022058s ring;
    std::cout << "Einstein Ring GAL-CLUS-022058s\n";
    std::cout << ring.primary_equation() << "\n";
    std::cout << "R_E = " << ring.einstein_radius()/ring.kpc << " kpc\n";
    std::cout << "Magnification(5") = " << ring.magnification(5.0) << "\n";
    ring.simulate(3);
    return 0;
}
#endif
