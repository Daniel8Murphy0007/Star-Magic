
// UQFFNebulaTriadicModule.h
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration) for Master Gravity Compressed, Resonance, and Buoyancy Equations in NGC 3596, NGC 1961, NGC 5335, NGC 2014, NGC 2020.
// This module can be plugged into a base program (e.g., 'nebula_triadic_sim.cpp') by including this header and linking the .cpp.
// Usage in base: #include "UQFFNebulaTriadicModule.h"
// UQFFNebulaTriadicModule mod; mod.computeMasterEquations(system); mod.updateVariable("F_rel", {new_real, new_imag});
// All variables are stored in a std::map for dynamic addition/subtraction/update, using complex<double> for real/imaginary components.
// Nothing is negligible: Includes all terms - Ug compressed, resonance U_r/U_dp, U_Bi buoyancy, DPM creation, SMBH dynamics, Triadic UQFF frame, gas nebula integration.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Integral approximated as integrand * x2 (quadratic root); imag parts small and not fully scaled; LENR dominant due to low ω_0; x2 from quadratic solver approx; F_rel from LEP 1998; [SSq], t_n placeholders.
// Multi-system params: NGC3596 M=1e41 kg r=1e21 m; NGC1961 M=2e41 kg r=2e21 m; NGC5335 M=3e41 kg r=3e21 m; NGC2014 M=4e41 kg r=4e21 m; NGC2020 M=5e41 kg r=5e21 m.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 23, 2025.

#ifndef UQFF_NEBULA_TRIADIC_MODULE_H
#define UQFF_NEBULA_TRIADIC_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

using cdouble = std::complex<double>;

class UQFFNebulaTriadicModule {
private:
    std::map<std::string, cdouble> variables;
    cdouble computeIntegrand(double t, const std::string& system);
    cdouble computeDPM_resonance(const std::string& system);
    cdouble computeX2(const std::string& system);
    cdouble computeQuadraticRoot(cdouble a, cdouble b, cdouble c);
    cdouble computeLENRTerm(const std::string& system);
    double computeG(double t, const std::string& system);
    cdouble computeQ_wave(double t, const std::string& system);
    cdouble computeUb1(const string& system);
    cdouble computeUi(double t, const std::string& system);
    void setSystemParams(const std::string& system);
    cdouble computeGravityCompressed(const std::string& system);
    cdouble computeResonanceUr(int U_dp, int U_r, const std::string& system);
    cdouble computeBuoyancyUbi(const std::string& system);
    cdouble computeGasNebulaIntegration(const std::string& system);

public:
    // Constructor: Initialize all variables with multi-system defaults
    UQFFNebulaTriadicModule();

    // Dynamic variable operations (complex)
    void updateVariable(const std::string& name, cdouble value);
    void addToVariable(const std::string& name, cdouble delta);
    void subtractFromVariable(const std::string& name, cdouble delta);

    // Core computation: Master Equations (F_U_Bi_i, Gravity Compressed, Resonance, Buoyancy) for system (approx integral)
    cdouble computeMasterEquations(const std::string& system, double t);

    // Sub-equations
    cdouble computeCompressed(const std::string& system, double t);  // Integrand
    cdouble computeResonant(const std::string& system);
    cdouble computeBuoyancy(const std::string& system);
    cdouble computeSuperconductive(const std::string& system, double t);
    double computeCompressedG(const std::string& system, double t);  // g(r,t)

    // Output descriptive text of the equation
    std::string getEquationText(const std::string& system);

    // Print all current variables (for debugging/updates)
    void printVariables();
};

#endif // UQFF_NEBULA_TRIADIC_MODULE_H

