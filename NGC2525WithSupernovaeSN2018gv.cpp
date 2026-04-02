// STANDALONE_NGC2525WITHSUPERNOVAESN2018GV
#include "NGC2525WithSupernovaeSN2018gv.h"
#include <iostream>

double NGC2525WithSupernovaeSN2018gv::absolute_magnitude(double delta_m15) const {
    return -19.3 + 0.74 * (delta_m15 - 1.1);
}

// Simplified Arnett light curve: Gaussian rise + exponential decline
double NGC2525WithSupernovaeSN2018gv::L_light_curve(double t) const {
    double rise  = L_SN_peak * std::exp(-0.5 * std::pow((t - t_rise) / (t_rise/3.0), 2));
    double decay = (t > t_rise) ? std::exp(-(t - t_rise) / t_decline) : 1.0;
    return rise * decay;
}

double NGC2525WithSupernovaeSN2018gv::L_SN_UQFF(double t) const {
    return L_light_curve(t) * (1.0 + rho_SCm / rho_UA) * (1.0 - f_TRZ);
}

// Bar + disk rotation: v_c = sqrt(G*M/R + G*M_bar*(1-exp(-R/R_bar))/R)
double NGC2525WithSupernovaeSN2018gv::v_circ_barred(double R) const {
    double v_disk = std::sqrt(G * M_galaxy / R);
    double v_bar  = std::sqrt(G * M_galaxy * 0.3 * (1.0 - std::exp(-R / R_bar)) / R);
    return std::sqrt(v_disk * v_disk + v_bar * v_bar);
}

std::string NGC2525WithSupernovaeSN2018gv::primary_equation() const {
    return "L_SN(t)=L_peak*exp(-0.5*((t-t_r)/sigma)^2)*exp(-(t-t_r)/t_d); L_UQFF=L_SN*(1+rho_SCm/rho_UA)*(1-f_TRZ)";
}

void NGC2525WithSupernovaeSN2018gv::self_update() { curr_t += time_step; }
void NGC2525WithSupernovaeSN2018gv::self_expand() { SFR_gal *= 1.001; }
void NGC2525WithSupernovaeSN2018gv::simulate(int num_steps) {
    for (int i = 0; i < num_steps; i++) {
        self_update();
        std::cout << "t=" << curr_t/86400 << " days  L_UQFF=" << L_SN_UQFF(curr_t) << " W\n";
    }
}

#ifdef STANDALONE_NGC2525WITHSUPERNOVAESN2018GV
int main() {
    NGC2525WithSupernovaeSN2018gv ngc2525;
    std::cout << "NGC 2525 with SN 2018gv\n";
    std::cout << ngc2525.primary_equation() << "\n";
    std::cout << "L_peak UQFF = " << ngc2525.L_SN_UQFF(ngc2525.t_rise) << " W\n";
    ngc2525.simulate(10);  // 10 day increments
    return 0;
}
#endif
