// STANDALONE_FORNAXCONSTELLATIONUHDF
#include "FornaxConstellationUHDF.h"
#include <iostream>

double FornaxConstellationUHDF::H_z(double zz) const {
    return H_0 * std::sqrt(omega_m * std::pow(1.0 + zz, 3) + omega_L);
}

// dV/dz = c/H(z) * D_c^2 * area_sr
// D_c ~ integral c/H(z')dz' (comoving distance) — approximated as (c/H_0)*z/sqrt(omega_m)
double FornaxConstellationUHDF::dV_dz(double z) const {
    double D_c = (c / H_0) * z / std::sqrt(omega_m * std::pow(1.0 + z, 3) + omega_L);
    return (c / H_z(z)) * D_c * D_c * area_sq_deg;
}

double FornaxConstellationUHDF::N_UQFF(double z) const {
    return N_galaxies * std::pow(1.0 + z, -1.5)
           * (1.0 + rho_SCm / rho_UA) * (1.0 + f_TRZ);
}

double FornaxConstellationUHDF::schechter_phi(double L, double L_star, double phi_star, double alpha) const {
    double x = L / L_star;
    return phi_star * std::pow(x, alpha) * std::exp(-x) / L_star;
}

std::string FornaxConstellationUHDF::primary_equation() const {
    return "N_UQFF(z)=N_obs*(1+z)^(-1.5)*(1+rho_SCm/rho_UA)*(1+f_TRZ); dV/dz=c/H(z)*D_c^2*dOmega";
}

void FornaxConstellationUHDF::self_update() { curr_t += time_step; }
void FornaxConstellationUHDF::self_expand() { z_max += 0.1; }
void FornaxConstellationUHDF::simulate(int num_steps) {
    for (int i = 0; i < num_steps; i++) {
        self_update();
        double z_curr = z_min + (z_max - z_min) * (i + 1.0) / num_steps;
        std::cout << "z=" << z_curr << "  N_UQFF=" << N_UQFF(z_curr) << "  dV=" << dV_dz(z_curr) << " m^3\n";
    }
}

#ifdef STANDALONE_FORNAXCONSTELLATIONUHDF
int main() {
    FornaxConstellationUHDF hudf;
    std::cout << "Fornax UHDF Galaxy Statistics\n";
    std::cout << hudf.primary_equation() << "\n";
    std::cout << "N_UQFF(z=1) = " << hudf.N_UQFF(1.0) << "\n";
    std::cout << "N_UQFF(z=6) = " << hudf.N_UQFF(6.0) << "\n";
    hudf.simulate(5);
    return 0;
}
#endif
