// InertiaUQFFModule.h
// UQFF for Quantum Waves, Inertial Operator, Universal Inertia,
// Bosonic Energy (Doc 43.d).
// Equations: wavefunction ψ, twist phase β·sin(ωt), Î inertial operator,
//            pseudo-monopole B, universal inertia, bosonic energy,
//            magnetic Hamiltonian, E_wave (~1.17e-105 J for n=4).
// Three-leg proofset: cons~1, vac ratio 1.683e-97, Q scale 3.333e-23.
// Copyright - Daniel T. Murphy.

#ifndef INERTIA_UQFF_MODULE_H
#define INERTIA_UQFF_MODULE_H

#include <map>
#include <string>
#include <complex>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>

enum class InertiaSystemType {
    QUANTUM_WAVES,
    INERTIAL_OPERATOR,
    UNIVERSAL_INERTIA,
    BOSONIC_ENERGY,
    GENERIC
};

class InertiaUQFFModule {
private:
    std::map<std::string, double> variables;
    InertiaSystemType current_system;

    double computeVacDensityRatio();
    double computeQuantumScalingFactor();
    double computeConservation(double E_in, double E_out);

public:
    InertiaUQFFModule(InertiaSystemType sys = InertiaSystemType::GENERIC);

    void setSystem(InertiaSystemType sys);

    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core equations
    std::complex<double> computeWaveFunction(double r, double theta, double phi, double t);
    double computeTwistPhase(double t);

    // Î = lambda_I * (d/dt + i*omega_m * r_vec · grad) ψ  — scalar approximation
    double computeInertialOperator(double psi, double t);

    // Pseudo-monopole B field
    double computePseudoMonopoleB(double r);

    // Universal inertia: lambda_I * (rho_SCm/rho_UA) * omega_i * cos(pi*t_n) * (1+F_RZ)
    double computeUniversalInertia(double t, double t_n);

    // Bosonic energy: 0.5*m*omega_r^2*x^2 + hbar*omega_r*(n+0.5)
    double computeBosonicEnergy(double x, int n);

    // Magnetic Hamiltonian: -mu*B
    double computeMagneticHamiltonian(double mu_val, double B);

    // E_wave: compound scaling ~ 1.17e-105 J for n_levels=4
    double computeEwave(int n_levels);

    // Three-leg proofset
    double computeThreeLegLeg1(double E_in, double E_out);
    double computeThreeLegLeg2();
    double computeThreeLegLeg3();

    double computeUQFF(double t);

    std::string getEquationText();
    std::string getSolutions(double t);

    void printVariables();
};

#endif // INERTIA_UQFF_MODULE_H
