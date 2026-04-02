// STANDALONE_AGNJETDYNAMICSBLANDFORDZNAJEK
#include "AGNJetDynamicsBlandfordZnajek.h"
#include <iostream>

// Gravitational radius r_g = G*M_BH/c^2
double AGNJetDynamicsBlandfordZnajek::r_g() const {
    return G * M_BH / (c * c);
}

// P_BZ: Blandford-Znajek jet power
// P_BZ ~ (1/4)(kappa_BZ)*phi_B^2*Omega_H^2/c
// Simplified: P_BZ = 0.044 * a^2 * B^2 * r_g^2 * c
double AGNJetDynamicsBlandfordZnajek::P_BZ() const {
    double rg = r_g();
    double kappa_BZ = 0.044;
    return kappa_BZ * a_spin * a_spin * B_field * B_field * rg * rg * c;
}

double AGNJetDynamicsBlandfordZnajek::lorentz_factor(double v) const {
    double beta = v / c;
    return 1.0 / std::sqrt(1.0 - beta * beta);
}

// Poynting flux: S = (E x B)/mu_0 ~ c*B^2/(4*pi*mu_0) approximated
double AGNJetDynamicsBlandfordZnajek::poynting_flux() const {
    return c * B_field * B_field / mu_0;
}

// Hoop stress: sigma = B_toroidal^2 / mu_0 (collimation force)
double AGNJetDynamicsBlandfordZnajek::hoop_stress(double B_toroidal) const {
    return B_toroidal * B_toroidal / mu_0;
}

// Jet luminosity (synchrotron + inverse Compton)
// L_jet ~ P_BZ * gamma_jet (bulk kinetic conversion)
double AGNJetDynamicsBlandfordZnajek::jet_luminosity() const {
    return P_BZ() * gamma_jet;
}

// UQFF-modulated jet acceleration
// Aether suppression: (1 - rho_SCm/rho_UA) and TRZ factor
double AGNJetDynamicsBlandfordZnajek::g_jet_UQFF() const {
    double p = P_BZ();
    double aether_factor = 1.0 - rho_SCm / rho_UA;
    double trz_factor    = 1.0 - f_TRZ;
    return p * aether_factor * trz_factor;
}

std::string AGNJetDynamicsBlandfordZnajek::primary_equation() const {
    return "P_BZ = kappa_BZ * a^2 * B^2 * r_g^2 * c; g_jet_UQFF = P_BZ*(1-rho_SCm/rho_UA)*(1-f_TRZ)";
}

void AGNJetDynamicsBlandfordZnajek::self_update() {
    curr_t += time_step;
    jet_length *= 1.001;  // jet propagates outward
}

void AGNJetDynamicsBlandfordZnajek::self_expand() {
    B_field *= 0.999;  // field weakens as jet expands
    gamma_jet *= 1.001;
}

void AGNJetDynamicsBlandfordZnajek::simulate(int num_steps) {
    for (int i = 0; i < num_steps; i++) {
        self_update();
        std::cout << "t=" << curr_t/3.156e7 << " yr  P_BZ=" << P_BZ()
                  << " W  L_jet=" << jet_luminosity() << " W\n";
    }
}

#ifdef STANDALONE_AGNJETDYNAMICSBLANDFORDZNAJEK
int main() {
    AGNJetDynamicsBlandfordZnajek agn;
    std::cout << "AGN Jet Dynamics (BZ Mechanism)\n";
    std::cout << agn.primary_equation() << "\n";
    std::cout << "P_BZ = " << agn.P_BZ() << " W\n";
    std::cout << "r_g  = " << agn.r_g() << " m\n";
    std::cout << "L_jet= " << agn.jet_luminosity() << " W\n";
    agn.simulate(3);
    return 0;
}
#endif
