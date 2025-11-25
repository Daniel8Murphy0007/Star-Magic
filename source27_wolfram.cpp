// source27_wolfram.cpp
// Wolfram companion file for NGC 1792 ("The Stellar Forge")
// Extracts all physics terms from source27.cpp into PhysicsTerm classes
// Astronomical system: NGC 1792 Starburst Galaxy (M0=1e10 M_sun, r=7.569e20 m, z=0.0095)
// Key: Intense star formation with M(t) growth, supernova feedback

#include "PhysicsTerm.h"
#include "PhysicsTermRegistry.h"
#include <cmath>
#include <memory>

// Forward declare delegation
void registerWolframTerms_source26(PhysicsTermRegistry& registry);

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

    std::string getName() const override { return "DynamicVacuum_NGC1792"; }
    std::string getDescription() const override {
        return "Time-varying vacuum energy for NGC 1792 starburst galaxy";
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
        double M = 1.989e40;      // 1e10 M_sun in kg
        double r = 7.569e20;      // 80k ly in m
        return coupling_strength * (hbar * hbar) / (M * r * r) * cos(t / 1e6);
    }

    std::string getName() const override { return "QuantumCoupling_NGC1792"; }
    std::string getDescription() const override {
        return "Non-local quantum effects for starburst galaxy";
    }
    std::string getCategory() const override { return "SelfExpanding"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"version", 2.0}, {"coupling_strength", coupling_strength}};
    }
};

// ========================================
// NGC 1792 Physics Terms
// ========================================

// Term 1: Base gravity with M(t), Hz, B corrections
class NGC1792BaseGravityTerm : public PhysicsTerm {
private:
    double G;       // 6.674e-11 m^3/kg·s^2
    double M0;      // 1.989e40 kg (1e10 M_sun)
    double r;       // 7.569e20 m (80k ly)
    double Hz;      // 2.19e-18 s^-1 (z=0.0095)
    double B;       // 1e-5 T
    double B_crit;  // 4.4e13 T
    double SFR_factor; // 10/1e10 = 1e-9
    double tau_SF;  // 3.156e15 s (100 Myr)

public:
    NGC1792BaseGravityTerm()
        : G(6.674e-11), M0(1.989e40), r(7.569e20), Hz(2.19e-18),
          B(1e-5), B_crit(4.4e13), SFR_factor(1e-9), tau_SF(3.156e15) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double M_t = M0 * (1 + M_dot);
        double ug1_t = (G * M_t) / (r * r);
        double corr_H = 1 + Hz * t;
        double corr_B = 1 - B / B_crit;
        return ug1_t * corr_H * corr_B;
    }

    std::string getName() const override { return "NGC1792_BaseGravity"; }
    std::string getDescription() const override {
        return "Base gravity with M(t), Hz, B: g·(1+Hz·t)·(1-B/B_c)";
    }
    std::string getCategory() const override { return "CoreGravity"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"M0_kg", M0}, {"r_m", r}, {"Hz_s-1", Hz}, {"z", 0.0095}, 
                {"SFR_factor", SFR_factor}, {"tau_SF_s", tau_SF}};
    }
};

// Term 2: UQFF Unification with f_TRZ and M(t)
class NGC1792UQFFUnificationTerm : public PhysicsTerm {
private:
    double G, M0, r;
    double f_TRZ;   // 0.01
    double B, B_crit;
    double SFR_factor, tau_SF;

public:
    NGC1792UQFFUnificationTerm()
        : G(6.674e-11), M0(1.989e40), r(7.569e20), f_TRZ(0.01),
          B(1e-5), B_crit(4.4e13), SFR_factor(1e-9), tau_SF(3.156e15) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double M_t = M0 * (1 + M_dot);
        double Ug1 = (G * M_t) / (r * r);
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - B / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ);
    }

    std::string getName() const override { return "NGC1792_UQFF_Unification"; }
    std::string getDescription() const override {
        return "UQFF unification Ug = (Ug1+Ug2+Ug3+Ug4)*(1+f_TRZ)";
    }
    std::string getCategory() const override { return "UQFF"; }
};

// Term 3: Cosmological constant Lambda
class NGC1792CosmologicalConstantTerm : public PhysicsTerm {
private:
    double Lambda;  // 1.1056e-52 m^-2
    double c_light; // 2.998e8 m/s

public:
    NGC1792CosmologicalConstantTerm()
        : Lambda(1.1056e-52), c_light(2.998e8) {}

    double compute(double t) const override {
        return (Lambda * c_light * c_light) / 3.0;
    }

    std::string getName() const override { return "NGC1792_CosmologicalConstant"; }
    std::string getDescription() const override {
        return "Cosmological constant Λc²/3 acceleration";
    }
    std::string getCategory() const override { return "DarkEnergy"; }
};

// Term 4: Scaled electromagnetic with UA vacuum
class NGC1792ElectromagneticTerm : public PhysicsTerm {
private:
    double q_charge;     // 1.602e-19 C
    double gas_v;        // 5e5 m/s
    double B;            // 1e-5 T
    double proton_mass;  // 1.673e-27 kg
    double rho_vac_UA;   // 6.3e-10 kg/m^3
    double rho_vac_SCm;  // 1e-26 kg/m^3
    double scale_EM;     // 1e-12

public:
    NGC1792ElectromagneticTerm()
        : q_charge(1.602e-19), gas_v(5e5), B(1e-5), proton_mass(1.673e-27),
          rho_vac_UA(6.3e-10), rho_vac_SCm(1e-26), scale_EM(1e-12) {}

    double compute(double t) const override {
        double cross_vB = gas_v * B;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
        return (em_base * corr_UA) * scale_EM;
    }

    std::string getName() const override { return "NGC1792_Electromagnetic"; }
    std::string getDescription() const override {
        return "Scaled EM with UA vacuum: (q*v×B/m_p)*(1+ρ_UA/ρ_SCm)*scale";
    }
    std::string getCategory() const override { return "Electromagnetic"; }
};

// Term 5: Quantum uncertainty
class NGC1792QuantumUncertaintyTerm : public PhysicsTerm {
private:
    double hbar;          // 1.0546e-34 J·s
    double delta_x;       // 1e20 m
    double delta_p;       // 1e-20 kg·m/s
    double integral_psi;  // 1e-5
    double t_Hubble;      // 4.35e17 s

public:
    NGC1792QuantumUncertaintyTerm()
        : hbar(1.0546e-34), delta_x(1e20), delta_p(1e-20),
          integral_psi(1e-5), t_Hubble(4.35e17) {}

    double compute(double t) const override {
        double sqrt_unc = sqrt(delta_x * delta_p);
        return (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);
    }

    std::string getName() const override { return "NGC1792_QuantumUncertainty"; }
    std::string getDescription() const override {
        return "Quantum uncertainty: (ℏ/√(Δx·Δp))·∫|ψ|²·(2π/t_H)";
    }
    std::string getCategory() const override { return "QuantumCorrection"; }
};

// Term 6: Fluid density coupling with M(t)
class NGC1792FluidDensityTerm : public PhysicsTerm {
private:
    double rho_fluid; // 1e-21 kg/m^3
    double r;         // 7.569e20 m
    double G, M0;
    double SFR_factor, tau_SF;

public:
    NGC1792FluidDensityTerm()
        : rho_fluid(1e-21), r(7.569e20), G(6.674e-11), M0(1.989e40),
          SFR_factor(1e-9), tau_SF(3.156e15) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double M_t = M0 * (1 + M_dot);
        double ug1_t = (G * M_t) / (r * r);
        double V = (4.0 / 3.0) * M_PI * r * r * r;
        return (rho_fluid * V * ug1_t) / M_t;
    }

    std::string getName() const override { return "NGC1792_FluidDensity"; }
    std::string getDescription() const override {
        return "Fluid density coupling: (ρ_fluid·V·g)/M(t)";
    }
    std::string getCategory() const override { return "FluidDynamics"; }
};

// Term 7: Oscillatory wave components
class NGC1792OscillatoryWaveTerm : public PhysicsTerm {
private:
    double A_osc;       // 1e-15 m/s^2
    double k_osc;       // 1e-21 m^-1
    double omega_osc;   // 1e-15 s^-1
    double x_pos;       // 1e20 m
    double t_Hubble_gyr; // 4.35e17 s (13.8 Gyr)

public:
    NGC1792OscillatoryWaveTerm()
        : A_osc(1e-15), k_osc(1e-21), omega_osc(1e-15),
          x_pos(1e20), t_Hubble_gyr(4.35e17) {}

    double compute(double t) const override {
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        return term_osc1 + term_osc2;
    }

    std::string getName() const override { return "NGC1792_OscillatoryWave"; }
    std::string getDescription() const override {
        return "Oscillatory wave: 2A·cos(kx)cos(ωt) + (2π/t_H)·A·cos(kx-ωt)";
    }
    std::string getCategory() const override { return "WavePhenomena"; }
};

// Term 8: Dark matter + density perturbations with M(t)
class NGC1792DarkMatterPerturbationTerm : public PhysicsTerm {
private:
    double G, M0, r;
    double M_DM_factor;         // 5.0
    double delta_rho_over_rho;  // 1e-5
    double SFR_factor, tau_SF;

public:
    NGC1792DarkMatterPerturbationTerm()
        : G(6.674e-11), M0(1.989e40), r(7.569e20),
          M_DM_factor(5.0), delta_rho_over_rho(1e-5),
          SFR_factor(1e-9), tau_SF(3.156e15) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double M_t = M0 * (1 + M_dot);
        double M_dm = M_t * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M_t / (r * r * r);
        double term_dm_force_like = (M_t + M_dm) * (pert1 + pert2);
        return term_dm_force_like / M_t;
    }

    std::string getName() const override { return "NGC1792_DarkMatterPerturbation"; }
    std::string getDescription() const override {
        return "Dark matter + density perturbation: (M(t)+M_DM)·(δρ/ρ + 3GM(t)/r³)/M(t)";
    }
    std::string getCategory() const override { return "DarkMatter"; }
};

// Term 9: Supernova feedback (pressure/density for acceleration)
class NGC1792SupernovaFeedbackTerm : public PhysicsTerm {
private:
    double rho_wind;  // 1e-21 kg/m^3
    double v_wind;    // 2e6 m/s
    double rho_fluid; // 1e-21 kg/m^3

public:
    NGC1792SupernovaFeedbackTerm()
        : rho_wind(1e-21), v_wind(2e6), rho_fluid(1e-21) {}

    double compute(double t) const override {
        double wind_pressure = rho_wind * v_wind * v_wind;
        return wind_pressure / rho_fluid;
    }

    std::string getName() const override { return "NGC1792_SupernovaFeedback"; }
    std::string getDescription() const override {
        return "Supernova feedback: ρ_wind·v_wind²/ρ_fluid for acceleration";
    }
    std::string getCategory() const override { return "SupernovaFeedback"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"rho_wind", rho_wind}, {"v_wind_m/s", v_wind}, {"wind_pressure_Pa", rho_wind * v_wind * v_wind}};
    }
};

// Term 10: Star formation M(t) - direct contribution
class NGC1792StarFormationTerm : public PhysicsTerm {
private:
    double G, M0, r;
    double SFR_factor;  // 1e-9
    double tau_SF;      // 3.156e15 s (100 Myr)

public:
    NGC1792StarFormationTerm()
        : G(6.674e-11), M0(1.989e40), r(7.569e20),
          SFR_factor(1e-9), tau_SF(3.156e15) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double delta_M = M0 * M_dot;
        // Mass growth increases gravity
        return (G * delta_M) / (r * r);
    }

    std::string getName() const override { return "NGC1792_StarFormation"; }
    std::string getDescription() const override {
        return "Star formation M(t) = M0·(1+SFR·e^(-t/τ)) increases gravity: G·ΔM/r²";
    }
    std::string getCategory() const override { return "StarFormation"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"SFR_factor", SFR_factor}, {"tau_SF_s", tau_SF}, {"tau_SF_Myr", tau_SF / 3.156e13}};
    }
};

// ========================================
// Registration Function
// ========================================

void registerWolframTerms_source27(PhysicsTermRegistry& registry) {
    // First, inherit all 254 terms from source26 (HUDF Galaxies)
    registerWolframTerms_source26(registry);

    // Register 2 self-expanding framework terms
    registry.registerTerm(std::make_unique<DynamicVacuumTerm>());
    registry.registerTerm(std::make_unique<QuantumCouplingTerm>());

    // Register 10 NGC 1792 physics terms
    registry.registerTerm(std::make_unique<NGC1792BaseGravityTerm>());
    registry.registerTerm(std::make_unique<NGC1792UQFFUnificationTerm>());
    registry.registerTerm(std::make_unique<NGC1792CosmologicalConstantTerm>());
    registry.registerTerm(std::make_unique<NGC1792ElectromagneticTerm>());
    registry.registerTerm(std::make_unique<NGC1792QuantumUncertaintyTerm>());
    registry.registerTerm(std::make_unique<NGC1792FluidDensityTerm>());
    registry.registerTerm(std::make_unique<NGC1792OscillatoryWaveTerm>());
    registry.registerTerm(std::make_unique<NGC1792DarkMatterPerturbationTerm>());
    registry.registerTerm(std::make_unique<NGC1792SupernovaFeedbackTerm>());
    registry.registerTerm(std::make_unique<NGC1792StarFormationTerm>());

    // Total: 254 (inherited) + 12 (new) = 266 classes
}
