// source23_wolfram.cpp
// Wolfram companion file for Antennae Galaxies (NGC 4038 & NGC 4039)
// Extracts all physics terms from source23.cpp into PhysicsTerm classes
// Astronomical system: NGC 4038/4039 (M0=2e11 M_sun, r=2.838e20 m, z=0.0105)
// Key: Interacting galaxy merger with M(t) star formation, I(t) interaction factor

#include "PhysicsTerm.h"
#include "PhysicsTermRegistry.h"
#include <cmath>
#include <memory>

// Forward declare delegation
void registerWolframTerms_source22(PhysicsTermRegistry& registry);

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

    std::string getName() const override { return "DynamicVacuum_Antennae"; }
    std::string getDescription() const override {
        return "Time-varying vacuum energy for Antennae Galaxies merger";
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
        double M = 3.978e41;      // 2e11 M_sun in kg
        double r = 2.838e20;      // 30k ly in m
        return coupling_strength * (hbar * hbar) / (M * r * r) * cos(t / 1e6);
    }

    std::string getName() const override { return "QuantumCoupling_Antennae"; }
    std::string getDescription() const override {
        return "Non-local quantum effects for galaxy merger";
    }
    std::string getCategory() const override { return "SelfExpanding"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"version", 2.0}, {"coupling_strength", coupling_strength}};
    }
};

// ========================================
// Antennae Galaxies Physics Terms
// ========================================

// Term 1: Base gravity with M(t) star formation, H(z), B, I(t) corrections
class AntennaeBaseGravityTerm : public PhysicsTerm {
private:
    double G;       // 6.674e-11 m^3/kg·s^2
    double M0;      // 3.978e41 kg (2e11 M_sun)
    double r;       // 2.838e20 m (30k ly)
    double Hz;      // 2.19e-18 s^-1 (z=0.0105)
    double B;       // 1e-5 T
    double B_crit;  // 4.4e13 T
    double SFR_factor; // 20 / 2e11 = 1e-10
    double tau_SF;  // 3.156e15 s (100 Myr)
    double I0;      // 0.1 (interaction factor)
    double tau_merger; // 1.262e16 s (400 Myr)

public:
    AntennaeBaseGravityTerm()
        : G(6.674e-11), M0(3.978e41), r(2.838e20), Hz(2.19e-18),
          B(1e-5), B_crit(4.4e13), SFR_factor(1e-10), tau_SF(3.156e15),
          I0(0.1), tau_merger(1.262e16) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double M_t = M0 * (1 + M_dot);
        double I_t = I0 * exp(-t / tau_merger);
        double ug1_t = (G * M_t) / (r * r);
        double corr_H = 1 + Hz * t;
        double corr_B = 1 - B / B_crit;
        double corr_I = 1 + I_t;
        return ug1_t * corr_H * corr_B * corr_I;
    }

    std::string getName() const override { return "Antennae_BaseGravity"; }
    std::string getDescription() const override {
        return "Base gravity with M(t) star formation, H(z), B, and I(t) merger interaction";
    }
    std::string getCategory() const override { return "CoreGravity"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"M0_kg", M0}, {"r_m", r}, {"Hz_s-1", Hz}, {"SFR_factor", SFR_factor}, 
                {"tau_SF_s", tau_SF}, {"I0", I0}, {"tau_merger_s", tau_merger}};
    }
};

// Term 2: UQFF Unification with f_TRZ and I(t)
class AntennaeUQFFUnificationTerm : public PhysicsTerm {
private:
    double G, M0, r;
    double f_TRZ;   // 0.01
    double B, B_crit;
    double SFR_factor, tau_SF;
    double I0, tau_merger;

public:
    AntennaeUQFFUnificationTerm()
        : G(6.674e-11), M0(3.978e41), r(2.838e20), f_TRZ(0.01),
          B(1e-5), B_crit(4.4e13), SFR_factor(1e-10), tau_SF(3.156e15),
          I0(0.1), tau_merger(1.262e16) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double M_t = M0 * (1 + M_dot);
        double I_t = I0 * exp(-t / tau_merger);
        double Ug1 = (G * M_t) / (r * r);
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - B / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ) * (1 + I_t);
    }

    std::string getName() const override { return "Antennae_UQFF_Unification"; }
    std::string getDescription() const override {
        return "UQFF unification Ug = (Ug1+Ug2+Ug3+Ug4)*(1+f_TRZ)*(1+I(t))";
    }
    std::string getCategory() const override { return "UQFF"; }
};

// Term 3: Cosmological constant Lambda
class AntennaeCosmologicalConstantTerm : public PhysicsTerm {
private:
    double Lambda;  // 1.1056e-52 m^-2
    double c_light; // 2.998e8 m/s

public:
    AntennaeCosmologicalConstantTerm()
        : Lambda(1.1056e-52), c_light(2.998e8) {}

    double compute(double t) const override {
        return (Lambda * c_light * c_light) / 3.0;
    }

    std::string getName() const override { return "Antennae_CosmologicalConstant"; }
    std::string getDescription() const override {
        return "Cosmological constant Λc²/3 acceleration";
    }
    std::string getCategory() const override { return "DarkEnergy"; }
};

// Term 4: Scaled electromagnetic with UA vacuum
class AntennaeElectromagneticTerm : public PhysicsTerm {
private:
    double q_charge;     // 1.602e-19 C
    double gas_v;        // 5e5 m/s
    double B;            // 1e-5 T
    double proton_mass;  // 1.673e-27 kg
    double rho_vac_UA;   // 6.3e-10 kg/m^3
    double rho_vac_SCm;  // 1e-26 kg/m^3
    double scale_EM;     // 1e-12

public:
    AntennaeElectromagneticTerm()
        : q_charge(1.602e-19), gas_v(5e5), B(1e-5), proton_mass(1.673e-27),
          rho_vac_UA(6.3e-10), rho_vac_SCm(1e-26), scale_EM(1e-12) {}

    double compute(double t) const override {
        double cross_vB = gas_v * B;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
        return (em_base * corr_UA) * scale_EM;
    }

    std::string getName() const override { return "Antennae_Electromagnetic"; }
    std::string getDescription() const override {
        return "Scaled EM with UA vacuum: (q*v×B/m_p)*(1+ρ_UA/ρ_SCm)*scale";
    }
    std::string getCategory() const override { return "Electromagnetic"; }
};

// Term 5: Quantum uncertainty
class AntennaeQuantumUncertaintyTerm : public PhysicsTerm {
private:
    double hbar;          // 1.0546e-34 J·s
    double delta_x;       // 1e20 m
    double delta_p;       // 1e-20 kg·m/s
    double integral_psi;  // 1e-5
    double t_Hubble;      // 4.35e17 s

public:
    AntennaeQuantumUncertaintyTerm()
        : hbar(1.0546e-34), delta_x(1e20), delta_p(1e-20),
          integral_psi(1e-5), t_Hubble(4.35e17) {}

    double compute(double t) const override {
        double sqrt_unc = sqrt(delta_x * delta_p);
        return (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);
    }

    std::string getName() const override { return "Antennae_QuantumUncertainty"; }
    std::string getDescription() const override {
        return "Quantum uncertainty: (ℏ/√(Δx·Δp))·∫|ψ|²·(2π/t_H)";
    }
    std::string getCategory() const override { return "QuantumCorrection"; }
};

// Term 6: Fluid density coupling
class AntennaeFluidDensityTerm : public PhysicsTerm {
private:
    double rho_fluid; // 1e-21 kg/m^3
    double r;         // 2.838e20 m
    double G, M0;
    double SFR_factor, tau_SF;

public:
    AntennaeFluidDensityTerm()
        : rho_fluid(1e-21), r(2.838e20), G(6.674e-11),
          M0(3.978e41), SFR_factor(1e-10), tau_SF(3.156e15) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double M_t = M0 * (1 + M_dot);
        double ug1_t = (G * M_t) / (r * r);
        double V = (4.0 / 3.0) * M_PI * r * r * r;
        return (rho_fluid * V * ug1_t) / M_t;
    }

    std::string getName() const override { return "Antennae_FluidDensity"; }
    std::string getDescription() const override {
        return "Fluid density coupling: (ρ_fluid·V·g)/M(t)";
    }
    std::string getCategory() const override { return "FluidDynamics"; }
};

// Term 7: Oscillatory wave components
class AntennaeOscillatoryWaveTerm : public PhysicsTerm {
private:
    double A_osc;       // 1e-15 m/s^2
    double k_osc;       // 1e-21 m^-1
    double omega_osc;   // 1e-15 s^-1
    double x_pos;       // 1e20 m
    double t_Hubble_gyr; // 4.35e17 s (13.8 Gyr)

public:
    AntennaeOscillatoryWaveTerm()
        : A_osc(1e-15), k_osc(1e-21), omega_osc(1e-15),
          x_pos(1e20), t_Hubble_gyr(4.35e17) {}

    double compute(double t) const override {
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        return term_osc1 + term_osc2;
    }

    std::string getName() const override { return "Antennae_OscillatoryWave"; }
    std::string getDescription() const override {
        return "Oscillatory wave: 2A·cos(kx)cos(ωt) + (2π/t_H)·A·cos(kx-ωt)";
    }
    std::string getCategory() const override { return "WavePhenomena"; }
};

// Term 8: Dark matter + density perturbations
class AntennaeDarkMatterPerturbationTerm : public PhysicsTerm {
private:
    double G, M0, r;
    double M_DM_factor;         // 5.0
    double delta_rho_over_rho;  // 1e-5
    double SFR_factor, tau_SF;

public:
    AntennaeDarkMatterPerturbationTerm()
        : G(6.674e-11), M0(3.978e41), r(2.838e20),
          M_DM_factor(5.0), delta_rho_over_rho(1e-5),
          SFR_factor(1e-10), tau_SF(3.156e15) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double M_t = M0 * (1 + M_dot);
        double M_dm = M_t * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M_t / (r * r * r);
        double term_dm_force_like = (M_t + M_dm) * (pert1 + pert2);
        return term_dm_force_like / M_t;
    }

    std::string getName() const override { return "Antennae_DarkMatterPerturbation"; }
    std::string getDescription() const override {
        return "Dark matter + density perturbation: (M+M_DM)·(δρ/ρ + 3GM/r³)/M";
    }
    std::string getCategory() const override { return "DarkMatter"; }
};

// Term 9: Stellar feedback (pressure/density for acceleration)
class AntennaeStellarFeedbackTerm : public PhysicsTerm {
private:
    double rho_wind;   // 1e-21 kg/m^3
    double v_wind;     // 2e6 m/s
    double rho_fluid;  // 1e-21 kg/m^3

public:
    AntennaeStellarFeedbackTerm()
        : rho_wind(1e-21), v_wind(2e6), rho_fluid(1e-21) {}

    double compute(double t) const override {
        double wind_pressure = rho_wind * v_wind * v_wind;
        return wind_pressure / rho_fluid;
    }

    std::string getName() const override { return "Antennae_StellarFeedback"; }
    std::string getDescription() const override {
        return "Stellar feedback: ρ_wind·v_wind²/ρ_fluid for acceleration";
    }
    std::string getCategory() const override { return "StellarFeedback"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"rho_wind", rho_wind}, {"v_wind_m/s", v_wind}, {"wind_pressure_Pa", rho_wind * v_wind * v_wind}};
    }
};

// ========================================
// Registration Function
// ========================================

void registerWolframTerms_source23(PhysicsTermRegistry& registry) {
    // First, inherit all 204 terms from source22 (Bubble Nebula)
    registerWolframTerms_source22(registry);

    // Register 2 self-expanding framework terms
    registry.registerTerm(std::make_unique<DynamicVacuumTerm>());
    registry.registerTerm(std::make_unique<QuantumCouplingTerm>());

    // Register 9 Antennae Galaxies physics terms
    registry.registerTerm(std::make_unique<AntennaeBaseGravityTerm>());
    registry.registerTerm(std::make_unique<AntennaeUQFFUnificationTerm>());
    registry.registerTerm(std::make_unique<AntennaeCosmologicalConstantTerm>());
    registry.registerTerm(std::make_unique<AntennaeElectromagneticTerm>());
    registry.registerTerm(std::make_unique<AntennaeQuantumUncertaintyTerm>());
    registry.registerTerm(std::make_unique<AntennaeFluidDensityTerm>());
    registry.registerTerm(std::make_unique<AntennaeOscillatoryWaveTerm>());
    registry.registerTerm(std::make_unique<AntennaeDarkMatterPerturbationTerm>());
    registry.registerTerm(std::make_unique<AntennaeStellarFeedbackTerm>());

    // Total: 204 (inherited) + 11 (new) = 215 classes
}
