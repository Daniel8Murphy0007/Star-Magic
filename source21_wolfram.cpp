// source21_wolfram.cpp
// Wolfram companion file for NGC 3603 Extreme Star Cluster
// Extracts all physics terms from source21.cpp into PhysicsTerm classes
// Astronomical system: NGC 3603 (M0=400,000 M_sun, r=8.998e15 m, z=0.0071)
// Key: Young massive star cluster, stellar wind feedback, cavity pressure P(t)

#include "PhysicsTerm.h"
#include "PhysicsTermRegistry.h"
#include <cmath>
#include <memory>

// Forward declare delegation
void registerWolframTerms_source20(PhysicsTermRegistry& registry);

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

    std::string getName() const override { return "DynamicVacuum_NGC3603"; }
    std::string getDescription() const override {
        return "Time-varying vacuum energy for NGC 3603 star cluster";
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
        double M = 7.96e35;       // 400,000 M_sun in kg
        double r = 8.998e15;      // 9.5 ly in m
        return coupling_strength * (hbar * hbar) / (M * r * r) * cos(t / 1e6);
    }

    std::string getName() const override { return "QuantumCoupling_NGC3603"; }
    std::string getDescription() const override {
        return "Non-local quantum effects for star cluster";
    }
    std::string getCategory() const override { return "SelfExpanding"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"version", 2.0}, {"coupling_strength", coupling_strength}};
    }
};

// ========================================
// NGC 3603 Star Cluster Physics Terms
// ========================================

// Term 1: Base gravity with mass growth M(t), H0, B corrections
class NGC3603BaseGravityTerm : public PhysicsTerm {
private:
    double G;       // 6.674e-11 m^3/kg·s^2
    double M0;      // 7.96e35 kg (400,000 M_sun)
    double r;       // 8.998e15 m (9.5 ly)
    double H0;      // 2.268e-18 s^-1
    double B;       // 1e-5 T
    double B_crit;  // 4.4e13 T
    double M_dot_factor; // 0.01 (dimensionless)
    double tau_SF;  // 3.156e13 s (1 Myr)

public:
    NGC3603BaseGravityTerm()
        : G(6.674e-11), M0(7.96e35), r(8.998e15), H0(2.268e-18),
          B(1e-5), B_crit(4.4e13), M_dot_factor(0.01), tau_SF(3.156e13) {}

    double compute(double t) const override {
        double M_t = M0 * (1 + M_dot_factor * (1 - exp(-t / tau_SF)));
        double ug1_t = (G * M_t) / (r * r);
        double corr_H = 1 + H0 * t;
        double corr_B = 1 - B / B_crit;
        return ug1_t * corr_H * corr_B;
    }

    std::string getName() const override { return "NGC3603_BaseGravity"; }
    std::string getDescription() const override {
        return "Base gravity with mass growth M(t), Hubble H0, and B-field corrections";
    }
    std::string getCategory() const override { return "CoreGravity"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"M0_kg", M0}, {"r_m", r}, {"H0_s-1", H0}, {"tau_SF_s", tau_SF}};
    }
};

// Term 2: UQFF Unification with f_TRZ
class NGC3603UQFFUnificationTerm : public PhysicsTerm {
private:
    double G;
    double M0;
    double r;
    double f_TRZ;   // 0.01
    double B, B_crit;
    double M_dot_factor;
    double tau_SF;

public:
    NGC3603UQFFUnificationTerm()
        : G(6.674e-11), M0(7.96e35), r(8.998e15), f_TRZ(0.01),
          B(1e-5), B_crit(4.4e13), M_dot_factor(0.01), tau_SF(3.156e13) {}

    double compute(double t) const override {
        double M_t = M0 * (1 + M_dot_factor * (1 - exp(-t / tau_SF)));
        double ug1_t = (G * M_t) / (r * r);
        double Ug1 = ug1_t;
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - B / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ);
    }

    std::string getName() const override { return "NGC3603_UQFF_Unification"; }
    std::string getDescription() const override {
        return "UQFF unification Ug = (Ug1+Ug2+Ug3+Ug4)*(1+f_TRZ) for star cluster";
    }
    std::string getCategory() const override { return "UQFF"; }
};

// Term 3: Cosmological constant Lambda
class NGC3603CosmologicalConstantTerm : public PhysicsTerm {
private:
    double Lambda;  // 1.1056e-52 m^-2
    double c_light; // 2.998e8 m/s

public:
    NGC3603CosmologicalConstantTerm()
        : Lambda(1.1056e-52), c_light(2.998e8) {}

    double compute(double t) const override {
        return (Lambda * c_light * c_light) / 3.0;
    }

    std::string getName() const override { return "NGC3603_CosmologicalConstant"; }
    std::string getDescription() const override {
        return "Cosmological constant Λc²/3 acceleration";
    }
    std::string getCategory() const override { return "DarkEnergy"; }
};

// Term 4: Scaled electromagnetic with UA vacuum
class NGC3603ElectromagneticTerm : public PhysicsTerm {
private:
    double q_charge;     // 1.602e-19 C
    double gas_v;        // 5e5 m/s
    double B;            // 1e-5 T
    double proton_mass;  // 1.673e-27 kg
    double rho_vac_UA;   // 6.3e-10 kg/m^3
    double rho_vac_SCm;  // 1e-26 kg/m^3
    double scale_EM;     // 1e-12

public:
    NGC3603ElectromagneticTerm()
        : q_charge(1.602e-19), gas_v(5e5), B(1e-5), proton_mass(1.673e-27),
          rho_vac_UA(6.3e-10), rho_vac_SCm(1e-26), scale_EM(1e-12) {}

    double compute(double t) const override {
        double cross_vB = gas_v * B;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
        return (em_base * corr_UA) * scale_EM;
    }

    std::string getName() const override { return "NGC3603_Electromagnetic"; }
    std::string getDescription() const override {
        return "Scaled EM with UA vacuum: (q*v×B/m_p)*(1+ρ_UA/ρ_SCm)*scale";
    }
    std::string getCategory() const override { return "Electromagnetic"; }
};

// Term 5: Quantum uncertainty
class NGC3603QuantumUncertaintyTerm : public PhysicsTerm {
private:
    double hbar;          // 1.0546e-34 J·s
    double delta_x;       // 1e15 m
    double delta_p;       // 1e-20 kg·m/s
    double integral_psi;  // 1e-5
    double t_Hubble;      // 4.35e17 s

public:
    NGC3603QuantumUncertaintyTerm()
        : hbar(1.0546e-34), delta_x(1e15), delta_p(1e-20),
          integral_psi(1e-5), t_Hubble(4.35e17) {}

    double compute(double t) const override {
        double sqrt_unc = sqrt(delta_x * delta_p);
        return (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);
    }

    std::string getName() const override { return "NGC3603_QuantumUncertainty"; }
    std::string getDescription() const override {
        return "Quantum uncertainty: (ℏ/√(Δx·Δp))·∫|ψ|²·(2π/t_H)";
    }
    std::string getCategory() const override { return "QuantumCorrection"; }
};

// Term 6: Fluid density coupling
class NGC3603FluidDensityTerm : public PhysicsTerm {
private:
    double rho_fluid; // 1e-21 kg/m^3
    double r;         // 8.998e15 m
    double G;
    double M0;
    double M_dot_factor;
    double tau_SF;

public:
    NGC3603FluidDensityTerm()
        : rho_fluid(1e-21), r(8.998e15), G(6.674e-11),
          M0(7.96e35), M_dot_factor(0.01), tau_SF(3.156e13) {}

    double compute(double t) const override {
        double M_t = M0 * (1 + M_dot_factor * (1 - exp(-t / tau_SF)));
        double ug1_t = (G * M_t) / (r * r);
        double V = (4.0 / 3.0) * M_PI * r * r * r;
        return (rho_fluid * V * ug1_t) / M_t;
    }

    std::string getName() const override { return "NGC3603_FluidDensity"; }
    std::string getDescription() const override {
        return "Fluid density coupling: (ρ_fluid·V·g)/M(t)";
    }
    std::string getCategory() const override { return "FluidDynamics"; }
};

// Term 7: Oscillatory wave components
class NGC3603OscillatoryWaveTerm : public PhysicsTerm {
private:
    double A_osc;       // 1e-15 m/s^2
    double k_osc;       // 1e-18 m^-1
    double omega_osc;   // 1e-15 s^-1
    double x_pos;       // 1e15 m
    double t_Hubble_gyr; // 4.35e17 s (13.8 Gyr)

public:
    NGC3603OscillatoryWaveTerm()
        : A_osc(1e-15), k_osc(1e-18), omega_osc(1e-15),
          x_pos(1e15), t_Hubble_gyr(4.35e17) {}

    double compute(double t) const override {
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        return term_osc1 + term_osc2;
    }

    std::string getName() const override { return "NGC3603_OscillatoryWave"; }
    std::string getDescription() const override {
        return "Oscillatory wave: 2A·cos(kx)cos(ωt) + (2π/t_H)·A·cos(kx-ωt)";
    }
    std::string getCategory() const override { return "WavePhenomena"; }
};

// Term 8: Dark matter + density perturbations
class NGC3603DarkMatterPerturbationTerm : public PhysicsTerm {
private:
    double G;
    double M0;
    double r;
    double M_DM_factor;         // 5.0
    double delta_rho_over_rho;  // 1e-5
    double M_dot_factor;
    double tau_SF;

public:
    NGC3603DarkMatterPerturbationTerm()
        : G(6.674e-11), M0(7.96e35), r(8.998e15),
          M_DM_factor(5.0), delta_rho_over_rho(1e-5),
          M_dot_factor(0.01), tau_SF(3.156e13) {}

    double compute(double t) const override {
        double M_t = M0 * (1 + M_dot_factor * (1 - exp(-t / tau_SF)));
        double M_dm = M_t * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M_t / (r * r * r);
        double term_dm_force_like = (M_t + M_dm) * (pert1 + pert2);
        return term_dm_force_like / M_t;
    }

    std::string getName() const override { return "NGC3603_DarkMatterPerturbation"; }
    std::string getDescription() const override {
        return "Dark matter + density perturbation: (M+M_DM)·(δρ/ρ + 3GM/r³)/M";
    }
    std::string getCategory() const override { return "DarkMatter"; }
};

// Term 9: Stellar wind feedback (pressure/density for acceleration)
class NGC3603StellarWindTerm : public PhysicsTerm {
private:
    double rho_wind;   // 1e-20 kg/m^3
    double v_wind;     // 2e6 m/s
    double rho_fluid;  // 1e-21 kg/m^3

public:
    NGC3603StellarWindTerm()
        : rho_wind(1e-20), v_wind(2e6), rho_fluid(1e-21) {}

    double compute(double t) const override {
        double wind_pressure = rho_wind * v_wind * v_wind;
        return wind_pressure / rho_fluid;
    }

    std::string getName() const override { return "NGC3603_StellarWind"; }
    std::string getDescription() const override {
        return "Stellar wind feedback: ρ_wind·v_wind²/ρ_fluid for acceleration";
    }
    std::string getCategory() const override { return "StellarFeedback"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"rho_wind", rho_wind}, {"v_wind_m/s", v_wind}, {"wind_pressure_Pa", rho_wind * v_wind * v_wind}};
    }
};

// Term 10: Cavity pressure P(t) / rho_fluid
class NGC3603CavityPressureTerm : public PhysicsTerm {
private:
    double P0;         // 4e-8 Pa (initial pressure)
    double tau_exp;    // 3.156e13 s (1 Myr expansion timescale)
    double rho_fluid;  // 1e-21 kg/m^3

public:
    NGC3603CavityPressureTerm()
        : P0(4e-8), tau_exp(3.156e13), rho_fluid(1e-21) {}

    double compute(double t) const override {
        double P_t = P0 * exp(-t / tau_exp);
        return P_t / rho_fluid;
    }

    std::string getName() const override { return "NGC3603_CavityPressure"; }
    std::string getDescription() const override {
        return "Cavity pressure acceleration: P(t)/ρ_fluid with P(t)=P0·e^(-t/τ_exp)";
    }
    std::string getCategory() const override { return "GasDynamics"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"P0_Pa", P0}, {"tau_exp_s", tau_exp}, {"tau_exp_Myr", tau_exp / 3.156e13}};
    }
};

// ========================================
// Registration Function
// ========================================

void registerWolframTerms_source21(PhysicsTermRegistry& registry) {
    // First, inherit all 181 terms from source20 (NGC 2525)
    registerWolframTerms_source20(registry);

    // Register 2 self-expanding framework terms
    registry.registerTerm(std::make_unique<DynamicVacuumTerm>());
    registry.registerTerm(std::make_unique<QuantumCouplingTerm>());

    // Register 10 NGC 3603 star cluster physics terms
    registry.registerTerm(std::make_unique<NGC3603BaseGravityTerm>());
    registry.registerTerm(std::make_unique<NGC3603UQFFUnificationTerm>());
    registry.registerTerm(std::make_unique<NGC3603CosmologicalConstantTerm>());
    registry.registerTerm(std::make_unique<NGC3603ElectromagneticTerm>());
    registry.registerTerm(std::make_unique<NGC3603QuantumUncertaintyTerm>());
    registry.registerTerm(std::make_unique<NGC3603FluidDensityTerm>());
    registry.registerTerm(std::make_unique<NGC3603OscillatoryWaveTerm>());
    registry.registerTerm(std::make_unique<NGC3603DarkMatterPerturbationTerm>());
    registry.registerTerm(std::make_unique<NGC3603StellarWindTerm>());
    registry.registerTerm(std::make_unique<NGC3603CavityPressureTerm>());

    // Total: 181 (inherited) + 12 (new) = 193 classes
}
