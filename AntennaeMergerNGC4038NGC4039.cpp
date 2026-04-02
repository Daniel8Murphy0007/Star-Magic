// STANDALONE_ANTENNAEMERGERNGC4038NGC4039
#include "AntennaeMergerNGC4038NGC4039.h"
#include <iostream>

double AntennaeMergerNGC4038NGC4039::E_binding() const {
    return -G * M_4038 * M_4039 / d_sep;
}

// t_dyn = 1.17 * v_rel^3 / (G^2 * M_tot * rho_avg * ln_Lambda)
double AntennaeMergerNGC4038NGC4039::t_dynamical_friction() const {
    double M_tot  = M_4038 + M_4039;
    double rho_avg = M_tot / (4.0/3.0 * M_PI * std::pow(d_sep/2.0, 3));
    double ln_L   = 10.0;  // Coulomb logarithm
    return 1.17 * std::pow(v_rel, 3) / (G * G * M_tot * rho_avg * ln_L);
}

double AntennaeMergerNGC4038NGC4039::g_merge_UQFF(double r) const {
    double M_tot = M_4038 + M_4039;
    double g0    = G * M_tot / (r * r);
    double shock = rho_shock * 1.0e51 * g0;  // shock fluid term
    return g0 * (1.0 + f_TRZ) + shock * (rho_SCm / rho_UA);
}

// SFR_UQFF enhanced by shock compression
double AntennaeMergerNGC4038NGC4039::SFR_UQFF() const {
    return SFR_burst * (rho_shock / rho_UA) * (1.0 + rho_SCm / rho_UA);
}

std::string AntennaeMergerNGC4038NGC4039::primary_equation() const {
    return "g_merge = G*(M1+M2)/r^2*(1+f_TRZ) + rho_shock*V*g*(rho_SCm/rho_UA); t_df=1.17*v^3/(G^2*M*rho*ln_L)";
}

void AntennaeMergerNGC4038NGC4039::self_update() { curr_t += time_step; d_sep -= v_rel * time_step * 0.01; }
void AntennaeMergerNGC4038NGC4039::self_expand() { rho_shock *= 1.001; }
void AntennaeMergerNGC4038NGC4039::simulate(int num_steps) {
    for (int i = 0; i < num_steps; i++) {
        self_update();
        std::cout << "t=" << curr_t/3.156e7 << " yr  d_sep=" << d_sep/kpc
                  << " kpc  SFR=" << SFR_UQFF() << " M_sun/yr\n";
    }
}

#ifdef STANDALONE_ANTENNAEMERGERNGC4038NGC4039
int main() {
    AntennaeMergerNGC4038NGC4039 ant;
    std::cout << "Antennae Galaxies Merger\n";
    std::cout << ant.primary_equation() << "\n";
    std::cout << "E_binding = " << ant.E_binding() << " J\n";
    std::cout << "t_df = " << ant.t_dynamical_friction()/3.156e7 << " yr\n";
    ant.simulate(3);
    return 0;
}
#endif
