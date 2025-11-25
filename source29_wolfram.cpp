// source29_wolfram.cpp
// Wolfram companion file for Sombrero Galaxy (M104)
// Extracts all physics terms from source29.cpp into PhysicsTerm classes
// Astronomical system: Sombrero Galaxy M104 (M=1e11 M_sun, r=2.36e20 m, z=0.0063, M_BH=1e9 M_sun)
// Key: Prominent dust lane, superconductivity correction (1-B/B_crit), massive SMBH, Virgo Cluster member

#include "PhysicsTerm.h"
#include "PhysicsTermRegistry.h"
#include <cmath>
#include <memory>
#include <complex>

// Forward declare delegation
void registerWolframTerms_source28(PhysicsTermRegistry& registry);

// ========================================
// Self-Expanding Framework Classes (2.0)
// ========================================

// Dynamic Vacuum Fluctuation Term
class DynamicVacuumTerm : public PhysicsTerm {
private:
    double amplitude;
    double frequency;

public:
    DynamicVacuumTerm(double amp = 1e-10, double freq = 1e-15)
        : amplitude(amp), frequency(freq) {}

    double compute(double t) const override {
        double rho_vac = 7.09e-36; // J/m^3
        return amplitude * rho_vac * sin(frequency * t);
    }

    std::string getName() const override { return "DynamicVacuum_Sombrero"; }
    std::string getDescription() const override {
        return "Time-varying vacuum energy for Sombrero galaxy";
    }
    std::string getCategory() const override { return "SelfExpanding"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"version", 2.0}, {"amplitude", amplitude}, {"frequency", frequency}};
    }
};

// Quantum Coupling Term (non-local entanglement)
class QuantumCouplingTerm : public PhysicsTerm {
private:
    double coupling_strength;

public:
    QuantumCouplingTerm(double strength = 1e-40)
        : coupling_strength(strength) {}

    double compute(double t) const override {
        double hbar = 1.0546e-34; // J·s
        double M = 1.989e41;      // 1e11 M_sun in kg
        double r = 2.36e20;       // 25k ly in m
        return coupling_strength * (hbar * hbar) / (M * r * r) * cos(t / 1e6);
    }

    std::string getName() const override { return "QuantumCoupling_Sombrero"; }
    std::string getDescription() const override {
        return "Non-local quantum effects for Sombrero galaxy";
    }
    std::string getCategory() const override { return "SelfExpanding"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"version", 2.0}, {"coupling_strength", coupling_strength}};
    }
};

// ========================================
// Sombrero Physics Terms
// ========================================

// Term 1: Base gravity with Hz, superconductivity correction, and f_TRZ
class SombreroBaseGravityTerm : public PhysicsTerm {
private:
    double G;       // 6.674e-11 m^3/kg·s^2
    double M;       // 1.989e41 kg (1e11 M_sun)
    double r;       // 2.36e20 m (25k ly)
    double H0;      // 70 km/s/Mpc
    double z;       // 0.0063 (Virgo Cluster)
    double Omega_m; // 0.3
    double Omega_Lambda; // 0.7
    double B;       // 1e-10 T
    double B_crit;  // 1e11 T (converted from 10^15 G)
    double f_TRZ;   // 0.1 time-reversal factor
    double Mpc_to_m; // 3.086e22 m

public:
    SombreroBaseGravityTerm()
        : G(6.674e-11), M(1.989e41), r(2.36e20), H0(70.0), z(0.0063),
          Omega_m(0.3), Omega_Lambda(0.7), B(1e-10), B_crit(1e11),
          f_TRZ(0.1), Mpc_to_m(3.086e22) {}

    double compute(double t) const override {
        // H(z) in s^-1
        double Hz_kms = H0 * sqrt(Omega_m * pow(1.0 + z, 3) + Omega_Lambda);
        double Hz = (Hz_kms * 1e3) / Mpc_to_m;
        
        double expansion = 1.0 + Hz * t;
        double sc_correction = 1.0 - (B / B_crit); // Superconductivity
        double tr_factor = 1.0 + f_TRZ;
        return (G * M / (r * r)) * expansion * sc_correction * tr_factor;
    }

    std::string getName() const override { return "Sombrero_BaseGravity"; }
    std::string getDescription() const override {
        return "Base gravity with Hz, SC, TR: g·(1+Hz·t)·(1-B/B_c)·(1+f_TRZ)";
    }
    std::string getCategory() const override { return "CoreGravity"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"M_kg", M}, {"r_m", r}, {"z", z}, {"B_crit_T", B_crit}, {"f_TRZ", f_TRZ}};
    }
};

// Term 2: Supermassive black hole (1e9 M_sun)
class SombreroSMBHTerm : public PhysicsTerm {
private:
    double G;       // 6.674e-11 m^3/kg·s^2
    double M_BH;    // 1e9 M_sun = 1.989e39 kg
    double r_BH;    // 1e15 m (core scale)

public:
    SombreroSMBHTerm()
        : G(6.674e-11), M_BH(1e9 * 1.989e30), r_BH(1e15) {}

    double compute(double t) const override {
        return (G * M_BH) / (r_BH * r_BH);
    }

    std::string getName() const override { return "Sombrero_SMBH"; }
    std::string getDescription() const override {
        return "Supermassive black hole: G·M_BH/r_BH² (M_BH=1e9 M_sun)";
    }
    std::string getCategory() const override { return "BlackHole"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"M_BH_kg", M_BH}, {"M_BH_Msun", M_BH / 1.989e30}, {"r_BH_m", r_BH}};
    }
};

// Term 3: UQFF Unification Ug sum (Ug1 + Ug4, Ug2/Ug3=0)
class SombreroUQFFUnificationTerm : public PhysicsTerm {
private:
    double G, M, r;
    double f_sc; // 1.0 superconductive factor

public:
    SombreroUQFFUnificationTerm()
        : G(6.674e-11), M(1.989e41), r(2.36e20), f_sc(1.0) {}

    double compute(double t) const override {
        double Ug1 = (G * M) / (r * r);
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double Ug4 = Ug1 * f_sc;
        return Ug1 + Ug2 + Ug3 + Ug4;
    }

    std::string getName() const override { return "Sombrero_UQFF_Unification"; }
    std::string getDescription() const override {
        return "UQFF unification Ug = Ug1 + Ug4 (Ug2/Ug3=0)";
    }
    std::string getCategory() const override { return "UQFF"; }
};

// Term 4: Cosmological constant Lambda
class SombreroCosmologicalConstantTerm : public PhysicsTerm {
private:
    double Lambda;  // 1.1e-52 m^-2
    double c_light; // 3e8 m/s

public:
    SombreroCosmologicalConstantTerm()
        : Lambda(1.1e-52), c_light(3e8) {}

    double compute(double t) const override {
        return (Lambda * c_light * c_light) / 3.0;
    }

    std::string getName() const override { return "Sombrero_CosmologicalConstant"; }
    std::string getDescription() const override {
        return "Cosmological constant Λc²/3 acceleration";
    }
    std::string getCategory() const override { return "DarkEnergy"; }
};

// Term 5: Quantum uncertainty Heisenberg
class SombreroQuantumUncertaintyTerm : public PhysicsTerm {
private:
    double hbar;          // 1.0546e-34 J·s
    double Delta_x;       // 1e20 m
    double Delta_p;       // 1.0546e-54 kg·m/s
    double integral_psi;  // 1.0 (ground state)
    double t_Hubble;      // 4.35e17 s

public:
    SombreroQuantumUncertaintyTerm()
        : hbar(1.0546e-34), Delta_x(1e20), Delta_p(1.0546e-54),
          integral_psi(1.0), t_Hubble(4.35e17) {}

    double compute(double t) const override {
        double unc = sqrt(Delta_x * Delta_p);
        return (hbar / unc) * integral_psi * (2 * M_PI / t_Hubble);
    }

    std::string getName() const override { return "Sombrero_QuantumUncertainty"; }
    std::string getDescription() const override {
        return "Quantum uncertainty: (ℏ/√(Δx·Δp))·∫|ψ|²·(2π/t_H)";
    }
    std::string getCategory() const override { return "QuantumCorrection"; }
};

// Term 6: Electromagnetic Lorentz force (q v × B)
class SombreroElectromagneticTerm : public PhysicsTerm {
private:
    double q_charge;     // 1.602e-19 C
    double v_orbit;      // 2e5 m/s
    double B;            // 1e-10 T
    double proton_mass;  // 1.673e-27 kg
    double rho_vac_UA;   // 7.09e-36 J/m^3
    double rho_vac_SCm;  // 7.09e-37 J/m^3
    double scale_macro;  // 1e-12

public:
    SombreroElectromagneticTerm()
        : q_charge(1.602e-19), v_orbit(2e5), B(1e-10), proton_mass(1.673e-27),
          rho_vac_UA(7.09e-36), rho_vac_SCm(7.09e-37), scale_macro(1e-12) {}

    double compute(double t) const override {
        double em_base = (q_charge * v_orbit * B) / proton_mass;
        double corr_UA = 1.0 + (rho_vac_UA / rho_vac_SCm); // ~11
        return em_base * corr_UA * scale_macro;
    }

    std::string getName() const override { return "Sombrero_Electromagnetic"; }
    std::string getDescription() const override {
        return "EM Lorentz force: (q·v·B/m_p)·(1+ρ_UA/ρ_SCm)·scale";
    }
    std::string getCategory() const override { return "Electromagnetic"; }
};

// Term 7: Fluid density coupling (dust lane ISM)
class SombreroFluidDensityTerm : public PhysicsTerm {
private:
    double rho_fluid; // 1e-21 kg/m^3
    double r;         // 2.36e20 m
    double G, M;

public:
    SombreroFluidDensityTerm()
        : rho_fluid(1e-21), r(2.36e20), G(6.674e-11), M(1.989e41) {}

    double compute(double t) const override {
        double V = (4.0 / 3.0) * M_PI * r * r * r;
        double g_base = (G * M) / (r * r);
        return rho_fluid * V * g_base;
    }

    std::string getName() const override { return "Sombrero_FluidDensity"; }
    std::string getDescription() const override {
        return "Dust lane ISM fluid coupling: ρ_fluid·V·g";
    }
    std::string getCategory() const override { return "FluidDynamics"; }
};

// Term 8: Resonant oscillatory waves (Aether-mediated)
class SombreroOscillatoryWaveTerm : public PhysicsTerm {
private:
    double A_osc;     // 1e-15 m/s^2
    double k_osc;     // 1e-21 m^-1
    double omega_osc; // 1e-15 s^-1
    double x_pos;     // 0.0 m (central)
    double exp_factor; // 2π/13.8 Gyr

public:
    SombreroOscillatoryWaveTerm()
        : A_osc(1e-15), k_osc(1e-21), omega_osc(1e-15),
          x_pos(0.0), exp_factor(2 * M_PI / 13.8) {}

    double compute(double t) const override {
        double cos_term = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        std::complex<double> exp_complex(A_osc * exp(std::complex<double>(0, arg)));
        double real_exp = exp_complex.real();
        return cos_term + exp_factor * real_exp;
    }

    std::string getName() const override { return "Sombrero_OscillatoryWave"; }
    std::string getDescription() const override {
        return "Resonant oscillatory wave: 2A·cos(kx)cos(ωt) + (2π/13.8)·A·Re[e^(i(kx-ωt))]";
    }
    std::string getCategory() const override { return "WavePhenomena"; }
};

// Term 9: Dark matter halo with perturbations
class SombreroDarkMatterPerturbationTerm : public PhysicsTerm {
private:
    double M_visible;   // 0.8 * 1.989e41 kg
    double M_DM;        // 0.2 * 1.989e41 kg (lower DM fraction for bulge-dominant)
    double delta_rho;   // 1e-25 kg/m^3
    double rho;         // 1e-20 kg/m^3
    double G, M, r;

public:
    SombreroDarkMatterPerturbationTerm()
        : M_visible(0.8 * 1.989e41), M_DM(0.2 * 1.989e41),
          delta_rho(1e-25), rho(1e-20),
          G(6.674e-11), M(1.989e41), r(2.36e20) {}

    double compute(double t) const override {
        double pert = delta_rho / rho;
        double curv = 3 * G * M / (r * r * r);
        return (M_visible + M_DM) * (pert + curv);
    }

    std::string getName() const override { return "Sombrero_DarkMatterPerturbation"; }
    std::string getDescription() const override {
        return "Dark matter halo: (M_vis+M_DM)·(δρ/ρ + 3GM/r³)";
    }
    std::string getCategory() const override { return "DarkMatter"; }
};

// Term 10: Dust friction/drag (prominent dust lane)
class SombreroDustFrictionTerm : public PhysicsTerm {
private:
    double rho_dust;  // 1e-22 kg/m^3
    double v_orbit;   // 2e5 m/s
    double rho_mass;  // 1e-15 kg/m^3
    double scale_macro; // 1e-12

public:
    SombreroDustFrictionTerm()
        : rho_dust(1e-22), v_orbit(2e5), rho_mass(1e-15), scale_macro(1e-12) {}

    double compute(double t) const override {
        double force_dust = rho_dust * (v_orbit * v_orbit);
        return (force_dust / rho_mass) * scale_macro;
    }

    std::string getName() const override { return "Sombrero_DustFriction"; }
    std::string getDescription() const override {
        return "Dust lane friction/drag: (ρ_dust·v²/ρ_mass)·scale";
    }
    std::string getCategory() const override { return "DustFriction"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"rho_dust", rho_dust}, {"v_orbit_m/s", v_orbit}, {"friction_force", rho_dust * v_orbit * v_orbit}};
    }
};

// Term 11: Superconductivity quantum field correction
class SombreroSuperconductivityTerm : public PhysicsTerm {
private:
    double G, M, r;
    double B;      // 1e-10 T
    double B_crit; // 1e11 T

public:
    SombreroSuperconductivityTerm()
        : G(6.674e-11), M(1.989e41), r(2.36e20), B(1e-10), B_crit(1e11) {}

    double compute(double t) const override {
        double g_base = (G * M) / (r * r);
        double sc_reduction = B / B_crit;
        // Returns the correction term: -g * (B/B_crit)
        return -g_base * sc_reduction;
    }

    std::string getName() const override { return "Sombrero_Superconductivity"; }
    std::string getDescription() const override {
        return "Superconductivity correction: -g·(B/B_crit) quantum field effect";
    }
    std::string getCategory() const override { return "Superconductivity"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"B_T", B}, {"B_crit_T", B_crit}, {"B/B_crit", B / B_crit}};
    }
};

// ========================================
// Registration Function
// ========================================

void registerWolframTerms_source29(PhysicsTermRegistry& registry) {
    // First, inherit all 278 terms from source28 (Andromeda)
    registerWolframTerms_source28(registry);

    // Register 2 self-expanding framework terms
    registry.registerTerm(std::make_unique<DynamicVacuumTerm>());
    registry.registerTerm(std::make_unique<QuantumCouplingTerm>());

    // Register 11 Sombrero physics terms
    registry.registerTerm(std::make_unique<SombreroBaseGravityTerm>());
    registry.registerTerm(std::make_unique<SombreroSMBHTerm>());
    registry.registerTerm(std::make_unique<SombreroUQFFUnificationTerm>());
    registry.registerTerm(std::make_unique<SombreroCosmologicalConstantTerm>());
    registry.registerTerm(std::make_unique<SombreroQuantumUncertaintyTerm>());
    registry.registerTerm(std::make_unique<SombreroElectromagneticTerm>());
    registry.registerTerm(std::make_unique<SombreroFluidDensityTerm>());
    registry.registerTerm(std::make_unique<SombreroOscillatoryWaveTerm>());
    registry.registerTerm(std::make_unique<SombreroDarkMatterPerturbationTerm>());
    registry.registerTerm(std::make_unique<SombreroDustFrictionTerm>());
    registry.registerTerm(std::make_unique<SombreroSuperconductivityTerm>());

    // Total: 278 (inherited) + 13 (new) = 291 classes
}
