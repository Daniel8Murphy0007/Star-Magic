// source20_wolfram.cpp
// Wolfram companion file for NGC 2525 barred spiral galaxy (z=0.016)
// Extracts all physics terms from source20.cpp into PhysicsTerm classes
// Astronomical system: NGC 2525 (M=1.0000225e10 M_sun, r=2.836e20 m, z=0.016)
// Key: Central SMBH (M_BH=2.25e7 M_sun), supernova mass loss M_SN(t)

#include "PhysicsTerm.h"
#include "PhysicsTermRegistry.h"
#include <cmath>
#include <memory>

// Forward declare delegation
void registerWolframTerms_source19(PhysicsTermRegistry& registry);

// ========================================
// Self-Expanding Framework Classes (2.0)
// ========================================

// Dynamic Vacuum Fluctuation Term
class DynamicVacuumTerm : public PhysicsTerm {
private:
    double rho_vac_base;
    double time_evolution_factor;

public:
    DynamicVacuumTerm(double rho_base = 1e-26, double time_factor = 1e-18)
        : rho_vac_base(rho_base), time_evolution_factor(time_factor) {}

    double compute(double t) const override {
        double rho_vac_t = rho_vac_base * (1 + time_evolution_factor * t);
        return rho_vac_t * 8.0 * M_PI * 6.674e-11 / (3.0 * pow(2.998e8, 2));
    }

    std::string getName() const override { return "DynamicVacuum_NGC2525"; }
    std::string getDescription() const override {
        return "Time-varying vacuum energy density (ρ_vac(t)) for NGC 2525 galaxy";
    }
    std::string getCategory() const override { return "SelfExpanding"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"version", 2.0}, {"rho_vac_base", rho_vac_base}, {"time_factor", time_evolution_factor}};
    }
};

// Quantum Coupling Term (non-local entanglement)
class QuantumCouplingTerm : public PhysicsTerm {
private:
    double coupling_strength;
    double coherence_length;

public:
    QuantumCouplingTerm(double strength = 1e-40, double length = 1e15)
        : coupling_strength(strength), coherence_length(length) {}

    double compute(double t) const override {
        double phase = 2.0 * M_PI * t / (365.25 * 86400); // Annual phase
        return coupling_strength * exp(-coherence_length / 1e20) * cos(phase);
    }

    std::string getName() const override { return "QuantumCoupling_NGC2525"; }
    std::string getDescription() const override {
        return "Non-local quantum entanglement correction for NGC 2525";
    }
    std::string getCategory() const override { return "SelfExpanding"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"version", 2.0}, {"coupling_strength", coupling_strength}, {"coherence_length", coherence_length}};
    }
};

// ========================================
// NGC 2525 Galaxy Physics Terms
// ========================================

// Term 1: Base gravity with H(z), B corrections
class NGC2525BaseGravityTerm : public PhysicsTerm {
private:
    double ug1_base; // ~G*M/r^2 = 4.187e-11 m/s^2
    double Hz;       // 2.19e-18 s^-1 for z=0.016
    double B;        // 1e-5 T
    double B_crit;   // 4.4e13 T

public:
    NGC2525BaseGravityTerm()
        : ug1_base(4.187e-11), Hz(2.19e-18), B(1e-5), B_crit(4.4e13) {}

    double compute(double t) const override {
        double corr_H = 1 + Hz * t;
        double corr_B = 1 - B / B_crit;
        return ug1_base * corr_H * corr_B;
    }

    std::string getName() const override { return "NGC2525_BaseGravity"; }
    std::string getDescription() const override {
        return "Galaxy base gravity with Hubble H(z=0.016) and B-field corrections";
    }
    std::string getCategory() const override { return "CoreGravity"; }
};

// Term 2: Black hole gravitational influence (central SMBH)
class NGC2525BlackHoleTerm : public PhysicsTerm {
private:
    double g_BH; // G*M_BH/r_BH^2 (cached)
    double G;
    double M_BH;   // 2.25e7 M_sun = 4.48e37 kg
    double r_BH;   // 1.496e11 m (1 AU)

public:
    NGC2525BlackHoleTerm()
        : G(6.674e-11), M_BH(4.48e37), r_BH(1.496e11) {
        g_BH = G * M_BH / (r_BH * r_BH);
    }

    double compute(double t) const override {
        return g_BH; // Static BH influence
    }

    std::string getName() const override { return "NGC2525_BlackHole"; }
    std::string getDescription() const override {
        return "Central SMBH gravitational influence (M_BH=2.25e7 M_sun, r_BH=1 AU)";
    }
    std::string getCategory() const override { return "BlackHolePhysics"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"M_BH_kg", M_BH}, {"r_BH_m", r_BH}, {"g_BH", g_BH}};
    }
};

// Term 3: UQFF unification with f_TRZ
class NGC2525UQFFUnificationTerm : public PhysicsTerm {
private:
    double ug1_base;
    double f_TRZ;  // 0.01
    double B, B_crit;

public:
    NGC2525UQFFUnificationTerm()
        : ug1_base(4.187e-11), f_TRZ(0.01), B(1e-5), B_crit(4.4e13) {}

    double compute(double t) const override {
        double Ug1 = ug1_base;
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - B / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ);
    }

    std::string getName() const override { return "NGC2525_UQFF_Unification"; }
    std::string getDescription() const override {
        return "UQFF unification term Ug = (Ug1+Ug2+Ug3+Ug4)*(1+f_TRZ) for galaxy";
    }
    std::string getCategory() const override { return "UQFF"; }
};

// Term 4: Cosmological constant Λ
class NGC2525CosmologicalConstantTerm : public PhysicsTerm {
private:
    double Lambda;  // 1.1056e-52 m^-2
    double c_light; // 2.998e8 m/s

public:
    NGC2525CosmologicalConstantTerm()
        : Lambda(1.1056e-52), c_light(2.998e8) {}

    double compute(double t) const override {
        return (Lambda * c_light * c_light) / 3.0;
    }

    std::string getName() const override { return "NGC2525_CosmologicalConstant"; }
    std::string getDescription() const override {
        return "Cosmological constant Λc²/3 acceleration for galaxy";
    }
    std::string getCategory() const override { return "DarkEnergy"; }
};

// Term 5: Scaled electromagnetic with UA vacuum
class NGC2525ElectromagneticTerm : public PhysicsTerm {
private:
    double q_charge;     // 1.602e-19 C
    double gas_v;        // 5e5 m/s
    double B;            // 1e-5 T
    double proton_mass;  // 1.673e-27 kg
    double rho_vac_UA;   // 6.3e-10 kg/m^3
    double rho_vac_SCm;  // 1e-26 kg/m^3
    double scale_EM;     // 1e-12

public:
    NGC2525ElectromagneticTerm()
        : q_charge(1.602e-19), gas_v(5e5), B(1e-5), proton_mass(1.673e-27),
          rho_vac_UA(6.3e-10), rho_vac_SCm(1e-26), scale_EM(1e-12) {}

    double compute(double t) const override {
        double cross_vB = gas_v * B;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
        return (em_base * corr_UA) * scale_EM;
    }

    std::string getName() const override { return "NGC2525_Electromagnetic"; }
    std::string getDescription() const override {
        return "Scaled EM with UA vacuum correction: (q*v×B/m_p)*(1+ρ_UA/ρ_SCm)*scale";
    }
    std::string getCategory() const override { return "Electromagnetic"; }
};

// Term 6: Quantum uncertainty
class NGC2525QuantumUncertaintyTerm : public PhysicsTerm {
private:
    double hbar;          // 1.055e-34 J·s
    double delta_x;       // 1e15 m
    double delta_p;       // 1e-20 kg·m/s
    double integral_psi;  // 1e-5
    double t_Hubble;      // 4.35e17 s

public:
    NGC2525QuantumUncertaintyTerm()
        : hbar(1.055e-34), delta_x(1e15), delta_p(1e-20),
          integral_psi(1e-5), t_Hubble(4.35e17) {}

    double compute(double t) const override {
        double sqrt_unc = sqrt(delta_x * delta_p);
        return (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);
    }

    std::string getName() const override { return "NGC2525_QuantumUncertainty"; }
    std::string getDescription() const override {
        return "Quantum uncertainty: (ℏ/√(Δx·Δp))·∫|ψ|²·(2π/t_H)";
    }
    std::string getCategory() const override { return "QuantumCorrection"; }
};

// Term 7: Fluid density coupling
class NGC2525FluidDensityTerm : public PhysicsTerm {
private:
    double rho_fluid; // 1e-21 kg/m^3
    double r;         // 2.836e20 m
    double ug1_base;  // 4.187e-11 m/s^2
    double M;         // 1.989e40 kg

public:
    NGC2525FluidDensityTerm()
        : rho_fluid(1e-21), r(2.836e20), ug1_base(4.187e-11), M(1.989e40) {}

    double compute(double t) const override {
        double V = (4.0 / 3.0) * M_PI * r * r * r;
        return (rho_fluid * V * ug1_base) / M;
    }

    std::string getName() const override { return "NGC2525_FluidDensity"; }
    std::string getDescription() const override {
        return "Fluid density coupling: (ρ_fluid·V·g_base)/M";
    }
    std::string getCategory() const override { return "FluidDynamics"; }
};

// Term 8: Oscillatory wave components
class NGC2525OscillatoryWaveTerm : public PhysicsTerm {
private:
    double A_osc;       // 1e-15 m/s^2
    double k_osc;       // 1e-18 m^-1
    double omega_osc;   // 1e-15 s^-1
    double x_pos;       // 1e20 m
    double t_Hubble_gyr; // 13.8 Gyr = 4.35e17 s

public:
    NGC2525OscillatoryWaveTerm()
        : A_osc(1e-15), k_osc(1e-18), omega_osc(1e-15),
          x_pos(1e20), t_Hubble_gyr(4.35e17) {}

    double compute(double t) const override {
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        return term_osc1 + term_osc2;
    }

    std::string getName() const override { return "NGC2525_OscillatoryWave"; }
    std::string getDescription() const override {
        return "Oscillatory wave: 2A·cos(kx)cos(ωt) + (2π/t_H)·A·cos(kx-ωt)";
    }
    std::string getCategory() const override { return "WavePhenomena"; }
};

// Term 9: Dark matter + density perturbations
class NGC2525DarkMatterPerturbationTerm : public PhysicsTerm {
private:
    double G;                   // 6.674e-11 m^3/kg·s^2
    double M;                   // 1.989e40 kg
    double r;                   // 2.836e20 m
    double M_DM_factor;         // 5.0
    double delta_rho_over_rho;  // 1e-5

public:
    NGC2525DarkMatterPerturbationTerm()
        : G(6.674e-11), M(1.989e40), r(2.836e20),
          M_DM_factor(5.0), delta_rho_over_rho(1e-5) {}

    double compute(double t) const override {
        double M_dm = M * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M / (r * r * r);
        double term_dm_force_like = (M + M_dm) * (pert1 + pert2);
        return term_dm_force_like / M;
    }

    std::string getName() const override { return "NGC2525_DarkMatterPerturbation"; }
    std::string getDescription() const override {
        return "Dark matter + density perturbation: (M+M_DM)·(δρ/ρ + 3GM/r³)/M";
    }
    std::string getCategory() const override { return "DarkMatter"; }
};

// Term 10: Supernova mass loss (time-dependent)
class NGC2525SupernovaMassLossTerm : public PhysicsTerm {
private:
    double G;       // 6.674e-11 m^3/kg·s^2
    double M_SN0;   // 2.79e30 kg (1.4 M_sun)
    double tau_SN;  // 3.156e7 s (1 year)
    double r;       // 2.836e20 m

public:
    NGC2525SupernovaMassLossTerm()
        : G(6.674e-11), M_SN0(2.79e30), tau_SN(3.156e7), r(2.836e20) {}

    double compute(double t) const override {
        double M_SN_t = M_SN0 * exp(-t / tau_SN);
        return -(G * M_SN_t) / (r * r); // Negative acceleration (mass loss)
    }

    std::string getName() const override { return "NGC2525_SupernovaMassLoss"; }
    std::string getDescription() const override {
        return "Supernova mass loss: -G·M_SN(t)/r² with M_SN(t)=M_SN0·e^(-t/τ_SN)";
    }
    std::string getCategory() const override { return "StellarEvolution"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"M_SN0_kg", M_SN0}, {"tau_SN_s", tau_SN}, {"decay_time_yr", tau_SN / 3.156e7}};
    }
};

// ========================================
// Registration Function
// ========================================

void registerWolframTerms_source20(PhysicsTermRegistry& registry) {
    // First, inherit all 169 terms from source19 (Einstein Ring)
    registerWolframTerms_source19(registry);

    // Register 2 self-expanding framework terms
    registry.registerTerm(std::make_unique<DynamicVacuumTerm>());
    registry.registerTerm(std::make_unique<QuantumCouplingTerm>());

    // Register 10 NGC 2525 galaxy physics terms
    registry.registerTerm(std::make_unique<NGC2525BaseGravityTerm>());
    registry.registerTerm(std::make_unique<NGC2525BlackHoleTerm>());
    registry.registerTerm(std::make_unique<NGC2525UQFFUnificationTerm>());
    registry.registerTerm(std::make_unique<NGC2525CosmologicalConstantTerm>());
    registry.registerTerm(std::make_unique<NGC2525ElectromagneticTerm>());
    registry.registerTerm(std::make_unique<NGC2525QuantumUncertaintyTerm>());
    registry.registerTerm(std::make_unique<NGC2525FluidDensityTerm>());
    registry.registerTerm(std::make_unique<NGC2525OscillatoryWaveTerm>());
    registry.registerTerm(std::make_unique<NGC2525DarkMatterPerturbationTerm>());
    registry.registerTerm(std::make_unique<NGC2525SupernovaMassLossTerm>());

    // Total: 169 (inherited) + 12 (new) = 181 classes
}
