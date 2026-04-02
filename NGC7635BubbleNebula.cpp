// STANDALONE_NGC7635BUBBLENEBULA
#include "NGC7635BubbleNebula.h"
#include <iostream>

double NGC7635BubbleNebula::L_wind() const {
    return 0.5 * M_dot_wind * v_wind * v_wind;
}

// R_b(t) = 0.88 * (L_w / rho_ISM)^(1/5) * t^(3/5) [Weaver 1977]
double NGC7635BubbleNebula::R_bubble_expansion(double t) const {
    return 0.88 * std::pow(L_wind() / rho_ISM, 0.2) * std::pow(t, 0.6);
}

// UQFF-enhanced wind velocity
double NGC7635BubbleNebula::v_wind_UQFF() const {
    return v_wind * (1.0 + f_TRZ) * std::sqrt(rho_UA / rho_SCm);
}

// Gamma_eff = 5/3 (monatomic ideal gas) -> compression = 4
double NGC7635BubbleNebula::shock_compression_ratio() const {
    double gamma = 5.0/3.0;
    return (gamma + 1.0) / (gamma - 1.0);
}

std::string NGC7635BubbleNebula::primary_equation() const {
    return "R_b(t)=0.88*(L_w/rho_ISM)^(1/5)*t^(3/5); v_UQFF=v_wind*(1+f_TRZ)*sqrt(rho_UA/rho_SCm)";
}

void NGC7635BubbleNebula::self_update() {
    curr_t += time_step;
    R_bubble = R_bubble_expansion(t_age + curr_t);
}
void NGC7635BubbleNebula::self_expand() { rho_ISM *= 0.999; }
void NGC7635BubbleNebula::simulate(int num_steps) {
    for (int i = 0; i < num_steps; i++) {
        self_update();
        std::cout << "t=" << (t_age+curr_t)/3.156e7 << " yr  R=" << R_bubble/3.086e16
                  << " ly  v_UQFF=" << v_wind_UQFF()/1e3 << " km/s\n";
    }
}

#ifdef STANDALONE_NGC7635BUBBLENEBULA
int main() {
    NGC7635BubbleNebula ngc7635;
    std::cout << "NGC 7635 Bubble Nebula\n";
    std::cout << ngc7635.primary_equation() << "\n";
    std::cout << "L_wind = " << ngc7635.L_wind() << " W\n";
    std::cout << "R(100 kyr) = " << ngc7635.R_bubble_expansion(ngc7635.t_age)/3.086e16 << " ly\n";
    ngc7635.simulate(3);
    return 0;
}
#endif
