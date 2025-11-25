// source25_wolfram.cpp
// Wolfram companion file for NGC 1275 (Perseus A - Magnetic Monster)
// Extracts all physics terms from source25.cpp into PhysicsTerm classes
// Astronomical system: NGC 1275 (M=1e11 M_sun, r=1.893e21 m, M_BH=8e8 M_sun, z=0.0176)
// Key: AGN with B(t) decay, F(t) filament support, black hole, cooling flows

#include "PhysicsTerm.h"
#include "PhysicsTermRegistry.h"
#include <cmath>
#include <memory>

// Forward declare delegation
void registerWolframTerms_source24(PhysicsTermRegistry& registry);

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

    std::string getName() const override { return "DynamicVacuum_NGC1275"; }
    std::string getDescription() const override {
        return "Time-varying vacuum energy for NGC 1275 AGN";
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
        double r = 1.893e21;      // 200k ly in m
        return coupling_strength * (hbar * hbar) / (M * r * r) * cos(t / 1e6);
    }

    std::string getName() const override { return "QuantumCoupling_NGC1275"; }
    std::string getDescription() const override {
        return "Non-local quantum effects for active galactic nucleus";
    }
    std::string getCategory() const override { return "SelfExpanding"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"version", 2.0}, {"coupling_strength", coupling_strength}};
    }
};

// ========================================
// NGC 1275 Physics Terms
// ========================================

// Term 1: Base gravity with Hz, B(t), F(t) corrections
class NGC1275BaseGravityTerm : public PhysicsTerm {
private:
    double G;       // 6.674e-11 m^3/kg·s^2
    double M;       // 1.989e41 kg (1e11 M_sun)
    double r;       // 1.893e21 m (200k ly)
    double Hz;      // 2.20e-18 s^-1 (z=0.0176)
    double B0;      // 5e-9 T
    double tau_B;   // 3.156e15 s (100 Myr)
    double B_crit;  // 4.4e13 T
    double F0;      // 0.1 (filament factor)
    double tau_fil; // 3.156e15 s (100 Myr)

public:
    NGC1275BaseGravityTerm()
        : G(6.674e-11), M(1.989e41), r(1.893e21), Hz(2.20e-18),
          B0(5e-9), tau_B(3.156e15), B_crit(4.4e13),
          F0(0.1), tau_fil(3.156e15) {}

    double compute(double t) const override {
        double B_t = B0 * exp(-t / tau_B);
        double F_t = F0 * exp(-t / tau_fil);
        double ug1 = (G * M) / (r * r);
        double corr_H = 1 + Hz * t;
        double corr_B = 1 - B_t / B_crit;
        double corr_F = 1 + F_t;
        return ug1 * corr_H * corr_B * corr_F;
    }

    std::string getName() const override { return "NGC1275_BaseGravity"; }
    std::string getDescription() const override {
        return "Base gravity with Hz, B(t) decay, F(t) filament: g·(1+Hz·t)·(1-B(t)/B_c)·(1+F(t))";
    }
    std::string getCategory() const override { return "CoreGravity"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"M_kg", M}, {"r_m", r}, {"Hz_s-1", Hz}, {"B0_T", B0}, 
                {"tau_B_s", tau_B}, {"F0", F0}, {"tau_fil_s", tau_fil}};
    }
};

// Term 2: Black hole contribution
class NGC1275BlackHoleTerm : public PhysicsTerm {
private:
    double G;       // 6.674e-11 m^3/kg·s^2
    double M_BH;    // 1.592e39 kg (8e8 M_sun)
    double r_BH;    // 2.5e12 m

public:
    NGC1275BlackHoleTerm()
        : G(6.674e-11), M_BH(1.592e39), r_BH(2.5e12) {}

    double compute(double t) const override {
        return (G * M_BH) / (r_BH * r_BH);
    }

    std::string getName() const override { return "NGC1275_BlackHole"; }
    std::string getDescription() const override {
        return "Central SMBH acceleration: G·M_BH/r_BH²";
    }
    std::string getCategory() const override { return "BlackHole"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"M_BH_kg", M_BH}, {"M_BH_Msun", 8e8}, {"r_BH_m", r_BH}};
    }
};

// Term 3: UQFF Unification with f_TRZ, B(t), F(t)
class NGC1275UQFFUnificationTerm : public PhysicsTerm {
private:
    double G, M, r;
    double f_TRZ;   // 0.01
    double B0, tau_B, B_crit;
    double F0, tau_fil;

public:
    NGC1275UQFFUnificationTerm()
        : G(6.674e-11), M(1.989e41), r(1.893e21), f_TRZ(0.01),
          B0(5e-9), tau_B(3.156e15), B_crit(4.4e13),
          F0(0.1), tau_fil(3.156e15) {}

    double compute(double t) const override {
        double B_t = B0 * exp(-t / tau_B);
        double F_t = F0 * exp(-t / tau_fil);
        double Ug1 = (G * M) / (r * r);
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - B_t / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ) * (1 + F_t);
    }

    std::string getName() const override { return "NGC1275_UQFF_Unification"; }
    std::string getDescription() const override {
        return "UQFF unification Ug = (Ug1+Ug2+Ug3+Ug4)*(1+f_TRZ)*(1+F(t))";
    }
    std::string getCategory() const override { return "UQFF"; }
};

// Term 4: Cosmological constant Lambda
class NGC1275CosmologicalConstantTerm : public PhysicsTerm {
private:
    double Lambda;  // 1.1056e-52 m^-2
    double c_light; // 2.998e8 m/s

public:
    NGC1275CosmologicalConstantTerm()
        : Lambda(1.1056e-52), c_light(2.998e8) {}

    double compute(double t) const override {
        return (Lambda * c_light * c_light) / 3.0;
    }

    std::string getName() const override { return "NGC1275_CosmologicalConstant"; }
    std::string getDescription() const override {
        return "Cosmological constant Λc²/3 acceleration";
    }
    std::string getCategory() const override { return "DarkEnergy"; }
};

// Term 5: Scaled electromagnetic with UA vacuum and B(t)
class NGC1275ElectromagneticTerm : public PhysicsTerm {
private:
    double q_charge;     // 1.602e-19 C
    double gas_v;        // 5e5 m/s
    double B0;           // 5e-9 T
    double tau_B;        // 3.156e15 s
    double proton_mass;  // 1.673e-27 kg
    double rho_vac_UA;   // 6.3e-10 kg/m^3
    double rho_vac_SCm;  // 1e-26 kg/m^3
    double scale_EM;     // 1e-12

public:
    NGC1275ElectromagneticTerm()
        : q_charge(1.602e-19), gas_v(5e5), B0(5e-9), tau_B(3.156e15),
          proton_mass(1.673e-27), rho_vac_UA(6.3e-10), rho_vac_SCm(1e-26), scale_EM(1e-12) {}

    double compute(double t) const override {
        double B_t = B0 * exp(-t / tau_B);
        double cross_vB = gas_v * B_t;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
        return (em_base * corr_UA) * scale_EM;
    }

    std::string getName() const override { return "NGC1275_Electromagnetic"; }
    std::string getDescription() const override {
        return "Scaled EM with UA vacuum and B(t): (q*v×B(t)/m_p)*(1+ρ_UA/ρ_SCm)*scale";
    }
    std::string getCategory() const override { return "Electromagnetic"; }
};

// Term 6: Quantum uncertainty
class NGC1275QuantumUncertaintyTerm : public PhysicsTerm {
private:
    double hbar;          // 1.0546e-34 J·s
    double delta_x;       // 1e21 m
    double delta_p;       // 1e-20 kg·m/s
    double integral_psi;  // 1e-5
    double t_Hubble;      // 4.35e17 s

public:
    NGC1275QuantumUncertaintyTerm()
        : hbar(1.0546e-34), delta_x(1e21), delta_p(1e-20),
          integral_psi(1e-5), t_Hubble(4.35e17) {}

    double compute(double t) const override {
        double sqrt_unc = sqrt(delta_x * delta_p);
        return (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);
    }

    std::string getName() const override { return "NGC1275_QuantumUncertainty"; }
    std::string getDescription() const override {
        return "Quantum uncertainty: (ℏ/√(Δx·Δp))·∫|ψ|²·(2π/t_H)";
    }
    std::string getCategory() const override { return "QuantumCorrection"; }
};

// Term 7: Fluid density coupling
class NGC1275FluidDensityTerm : public PhysicsTerm {
private:
    double rho_fluid; // 1e-21 kg/m^3
    double r;         // 1.893e21 m
    double G, M;

public:
    NGC1275FluidDensityTerm()
        : rho_fluid(1e-21), r(1.893e21), G(6.674e-11), M(1.989e41) {}

    double compute(double t) const override {
        double ug1 = (G * M) / (r * r);
        double V = (4.0 / 3.0) * M_PI * r * r * r;
        return (rho_fluid * V * ug1) / M;
    }

    std::string getName() const override { return "NGC1275_FluidDensity"; }
    std::string getDescription() const override {
        return "Fluid density coupling: (ρ_fluid·V·g)/M";
    }
    std::string getCategory() const override { return "FluidDynamics"; }
};

// Term 8: Oscillatory wave components
class NGC1275OscillatoryWaveTerm : public PhysicsTerm {
private:
    double A_osc;       // 1e-15 m/s^2
    double k_osc;       // 1e-22 m^-1
    double omega_osc;   // 1e-15 s^-1
    double x_pos;       // 1e21 m
    double t_Hubble_gyr; // 4.35e17 s (13.8 Gyr)

public:
    NGC1275OscillatoryWaveTerm()
        : A_osc(1e-15), k_osc(1e-22), omega_osc(1e-15),
          x_pos(1e21), t_Hubble_gyr(4.35e17) {}

    double compute(double t) const override {
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        return term_osc1 + term_osc2;
    }

    std::string getName() const override { return "NGC1275_OscillatoryWave"; }
    std::string getDescription() const override {
        return "Oscillatory wave: 2A·cos(kx)cos(ωt) + (2π/t_H)·A·cos(kx-ωt)";
    }
    std::string getCategory() const override { return "WavePhenomena"; }
};

// Term 9: Dark matter + density perturbations
class NGC1275DarkMatterPerturbationTerm : public PhysicsTerm {
private:
    double G, M, r;
    double M_DM_factor;         // 5.0
    double delta_rho_over_rho;  // 1e-5

public:
    NGC1275DarkMatterPerturbationTerm()
        : G(6.674e-11), M(1.989e41), r(1.893e21),
          M_DM_factor(5.0), delta_rho_over_rho(1e-5) {}

    double compute(double t) const override {
        double M_dm = M * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M / (r * r * r);
        double term_dm_force_like = (M + M_dm) * (pert1 + pert2);
        return term_dm_force_like / M;
    }

    std::string getName() const override { return "NGC1275_DarkMatterPerturbation"; }
    std::string getDescription() const override {
        return "Dark matter + density perturbation: (M+M_DM)·(δρ/ρ + 3GM/r³)/M";
    }
    std::string getCategory() const override { return "DarkMatter"; }
};

// Term 10: Cooling flow term (pressure/density for acceleration)
class NGC1275CoolingFlowTerm : public PhysicsTerm {
private:
    double rho_cool;  // 1e-20 kg/m^3
    double v_cool;    // 3e3 m/s
    double rho_fluid; // 1e-21 kg/m^3

public:
    NGC1275CoolingFlowTerm()
        : rho_cool(1e-20), v_cool(3e3), rho_fluid(1e-21) {}

    double compute(double t) const override {
        double cool_pressure = rho_cool * v_cool * v_cool;
        return cool_pressure / rho_fluid;
    }

    std::string getName() const override { return "NGC1275_CoolingFlow"; }
    std::string getDescription() const override {
        return "Cooling flow: ρ_cool·v_cool²/ρ_fluid for acceleration";
    }
    std::string getCategory() const override { return "CoolingFlow"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"rho_cool", rho_cool}, {"v_cool_m/s", v_cool}, {"cool_pressure_Pa", rho_cool * v_cool * v_cool}};
    }
};

// Term 11: Magnetic field B(t) decay - direct contribution
class NGC1275MagneticDecayTerm : public PhysicsTerm {
private:
    double B0;       // 5e-9 T
    double tau_B;    // 3.156e15 s (100 Myr)
    double B_crit;   // 4.4e13 T
    double G, M, r;

public:
    NGC1275MagneticDecayTerm()
        : B0(5e-9), tau_B(3.156e15), B_crit(4.4e13),
          G(6.674e-11), M(1.989e41), r(1.893e21) {}

    double compute(double t) const override {
        double B_t = B0 * exp(-t / tau_B);
        double ug1 = (G * M) / (r * r);
        // B field reduces gravity by (1 - B/B_crit), so decay increases it
        return -ug1 * (B_t / B_crit);
    }

    std::string getName() const override { return "NGC1275_MagneticDecay"; }
    std::string getDescription() const override {
        return "Magnetic field decay B(t) = B0·e^(-t/τ_B) affects gravity: -g·(B(t)/B_c)";
    }
    std::string getCategory() const override { return "MagneticField"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"B0_T", B0}, {"tau_B_s", tau_B}, {"tau_B_Myr", tau_B / 3.156e13}, {"B_crit_T", B_crit}};
    }
};

// Term 12: Filament support F(t) - direct contribution
class NGC1275FilamentSupportTerm : public PhysicsTerm {
private:
    double F0;       // 0.1
    double tau_fil;  // 3.156e15 s (100 Myr)
    double G, M, r;

public:
    NGC1275FilamentSupportTerm()
        : F0(0.1), tau_fil(3.156e15),
          G(6.674e-11), M(1.989e41), r(1.893e21) {}

    double compute(double t) const override {
        double F_t = F0 * exp(-t / tau_fil);
        double ug1 = (G * M) / (r * r);
        // Filament increases gravity by (1+F(t))
        return ug1 * F_t;
    }

    std::string getName() const override { return "NGC1275_FilamentSupport"; }
    std::string getDescription() const override {
        return "Filament support F(t) = F0·e^(-t/τ_fil) enhances gravity: g·F(t)";
    }
    std::string getCategory() const override { return "FilamentDynamics"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"F0", F0}, {"tau_fil_s", tau_fil}, {"tau_fil_Myr", tau_fil / 3.156e13}};
    }
};

// ========================================
// Registration Function
// ========================================

void registerWolframTerms_source25(PhysicsTermRegistry& registry) {
    // First, inherit all 227 terms from source24 (Horsehead Nebula)
    registerWolframTerms_source24(registry);

    // Register 2 self-expanding framework terms
    registry.registerTerm(std::make_unique<DynamicVacuumTerm>());
    registry.registerTerm(std::make_unique<QuantumCouplingTerm>());

    // Register 12 NGC 1275 physics terms
    registry.registerTerm(std::make_unique<NGC1275BaseGravityTerm>());
    registry.registerTerm(std::make_unique<NGC1275BlackHoleTerm>());
    registry.registerTerm(std::make_unique<NGC1275UQFFUnificationTerm>());
    registry.registerTerm(std::make_unique<NGC1275CosmologicalConstantTerm>());
    registry.registerTerm(std::make_unique<NGC1275ElectromagneticTerm>());
    registry.registerTerm(std::make_unique<NGC1275QuantumUncertaintyTerm>());
    registry.registerTerm(std::make_unique<NGC1275FluidDensityTerm>());
    registry.registerTerm(std::make_unique<NGC1275OscillatoryWaveTerm>());
    registry.registerTerm(std::make_unique<NGC1275DarkMatterPerturbationTerm>());
    registry.registerTerm(std::make_unique<NGC1275CoolingFlowTerm>());
    registry.registerTerm(std::make_unique<NGC1275MagneticDecayTerm>());
    registry.registerTerm(std::make_unique<NGC1275FilamentSupportTerm>());

    // Total: 227 (inherited) + 14 (new) = 241 classes
}
