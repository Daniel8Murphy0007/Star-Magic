// STANDALONE_UQFFEQUATIONMATHEMATICALDERIVATION
#include "UQFFEquationMathematicalDerivation.h"
#include <iostream>

// Ug1_i = mu_dipole * B_i (B_i = B0 / (i+1) as layers weaken)
double UQFFEquationMathematicalDerivation::Ug1(int layer, double r) const {
    double mu_dipole = I_coil * A_coil * omega_mag;
    double B_i       = 1.0e-10 / (layer + 1);  // layer-dependent B
    return mu_dipole * B_i;
}

// Ug2_i = (mu_0 * H_aether / (i+1))^2 / (2*mu_0)
double UQFFEquationMathematicalDerivation::Ug2(int layer) const {
    double B_super = mu_0 * H_aether / (layer + 1);
    return B_super * B_super / (2.0 * mu_0);
}

double UQFFEquationMathematicalDerivation::Ug3(int layer, double r, double M_ext) const {
    double w = 1.0 / (layer + 1.0);  // layer weight ~1/i
    return w * G * M_ext / (r * r);
}

// Ug4_i = (k_4 * E_react_0 * exp(-kappa * t)) / N_DIM
double UQFFEquationMathematicalDerivation::Ug4(int layer, double t) const {
    (void)layer;
    return k_4 * E_react_0 * std::exp(-kappa / 86400.0 * t) / N_DIM;
}

// Full sum over 26 layers
double UQFFEquationMathematicalDerivation::g_UQFF_26D(double r, double t, double M_ext) const {
    double g_total = 0.0;
    for (int i = 0; i < N_DIM; i++) {
        g_total += Ug1(i, r) + Ug2(i) + Ug3(i, r, M_ext) + Ug4(i, t);
    }
    double U_i_term = (rho_SCm / rho_UA) * omega_mag * std::cos(M_PI * 0.0) * (1.0 + f_TRZ);
    double g_lambda = Lambda * c * c / 3.0;
    return g_total + U_i_term + g_lambda;
}

// V_UQFF = -G*M/(r) * (1 + rho_SCm/rho_UA) * (1 - f_TRZ)
double UQFFEquationMathematicalDerivation::V_UQFF(double r, double M) const {
    return -G * M / r * (1.0 + rho_SCm / rho_UA) * (1.0 - f_TRZ);
}

// psi(r,t) = A * cos(k*r - omega*t) (real part of plane wave solution)
double UQFFEquationMathematicalDerivation::psi_gravity_wave(double r, double t, double E) const {
    double m = M_sun;  // reference mass
    double k = std::sqrt(2.0 * m * std::abs(E)) / hbar;
    double omega_gw = E / hbar;
    return std::cos(k * r - omega_gw * t);
}

std::string UQFFEquationMathematicalDerivation::primary_equation() const {
    return "g_UQFF_26D(r,t) = sum_{i=1}^{26}[Ug1_i+Ug2_i+Ug3_i+Ug4_i] + U_i + Lambda*c^2/3; V_UQFF=-G*M/r*(1+rho_SCm/rho_UA)*(1-f_TRZ)";
}

void UQFFEquationMathematicalDerivation::self_update() { curr_t += time_step; }
void UQFFEquationMathematicalDerivation::self_expand() { E_react_0 *= 0.999; }
void UQFFEquationMathematicalDerivation::simulate(int num_steps) {
    double r_test = 10.0 * kpc;
    for (int i = 0; i < num_steps; i++) {
        self_update();
        double g = g_UQFF_26D(r_test, curr_t);
        std::cout << "t=" << curr_t/3.156e7 << " yr  g_26D=" << g << " m/s^2\n";
    }
}

#ifdef STANDALONE_UQFFEQUATIONMATHEMATICALDERIVATION
int main() {
    UQFFEquationMathematicalDerivation uqff;
    std::cout << "UQFF Equation Mathematical Derivation (26D)\n";
    std::cout << uqff.primary_equation() << "\n";
    std::cout << "g_26D(10 kpc, 0) = " << uqff.g_UQFF_26D(10.0*uqff.kpc, 0) << " m/s^2\n";
    uqff.simulate(3);
    return 0;
}
#endif
