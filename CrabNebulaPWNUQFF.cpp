// STANDALONE_CRABNEBULAPWNUQFF
#include "CrabNebulaPWNUQFF.h"
#include <iostream>

// L_sd = I * |Omega * dOmega/dt|
// dOmega/dt from magnetic dipole braking: -B^2*R_ns^6*Omega^3/(6*I*c^3)
double CrabNebulaPWNUQFF::spin_down_luminosity() const {
    double R_ns     = 1.0e4;   // neutron star radius ~10 km
    double dOmega   = -B_pulsar*B_pulsar * std::pow(R_ns, 6) * std::pow(Omega_psr, 3)
                      / (6.0 * I_psr * c * c * c);
    return -I_psr * Omega_psr * dOmega;
}

// t_sync = 9*m_e^3*c^5 / (4*r_e^2*c*B^2*gamma_e)  [simplified]
double CrabNebulaPWNUQFF::synchrotron_cooling(double gamma_e) const {
    double m_e  = 9.109e-31;
    double r_e  = 2.818e-15;
    return 9.0 * m_e*m_e*m_e * c*c*c*c*c
           / (4.0 * r_e*r_e * c * B_pulsar*B_pulsar * gamma_e);
}

// v_SNR_UQFF: Sedov-Taylor with UQFF suppression
// v = sqrt(2 * E_SN * (1-f_TRZ) / M_ejecta) * (1 + rho_SCm/rho_UA)
double CrabNebulaPWNUQFF::v_SNR_UQFF() const {
    double v0 = std::sqrt(2.0 * E_SN * (1.0 - f_TRZ) / M_ejecta);
    return v0 * (1.0 + rho_SCm / rho_UA);
}

// B_energy = B^2/(2*mu_0) [magnetic energy density in PWN]
double CrabNebulaPWNUQFF::B_energy_PWN() const {
    double B_PWN = 5.0e-8;  // typical PWN field ~50 nT
    return B_PWN * B_PWN / (2.0 * mu_0);
}

std::string CrabNebulaPWNUQFF::primary_equation() const {
    return "L_sd=I*|Omega*dOmega_dt|; v_SNR=sqrt(2*E_SN*(1-f_TRZ)/M_ej)*(1+rho_SCm/rho_UA)";
}

void CrabNebulaPWNUQFF::self_update() { curr_t += time_step; R_SNR += v_SNR_UQFF() * time_step; }
void CrabNebulaPWNUQFF::self_expand() { rho_PWN *= 0.99; }
void CrabNebulaPWNUQFF::simulate(int num_steps) {
    for (int i = 0; i < num_steps; i++) {
        self_update();
        std::cout << "t=" << (age_SNR+curr_t)/3.156e7 << " yr  R_SNR=" << R_SNR/3.086e16
                  << " ly  L_sd=" << spin_down_luminosity() << " W\n";
    }
}

#ifdef STANDALONE_CRABNEBULAPWNUQFF
int main() {
    CrabNebulaPWNUQFF crab;
    std::cout << "Crab Nebula PWN UQFF\n";
    std::cout << crab.primary_equation() << "\n";
    std::cout << "L_sd = " << crab.spin_down_luminosity() << " W\n";
    std::cout << "v_SNR_UQFF = " << crab.v_SNR_UQFF()/1e3 << " km/s\n";
    crab.simulate(3);
    return 0;
}
#endif
