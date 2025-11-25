// source26_wolfram.cpp
// Wolfram companion file for Hubble Ultra Deep Field (HUDF) Galaxies
// Extracts all physics terms from source26.cpp into PhysicsTerm classes
// Astronomical system: HUDF Galaxies Field (M0=1e12 M_sun, r=1.23e27 m, z_avg=3.5)
// Key: Cosmic deep field with M(t) star formation, I(t) galaxy interactions

#include "PhysicsTerm.h"
#include "PhysicsTermRegistry.h"
#include <cmath>
#include <memory>

// Forward declare delegation
void registerWolframTerms_source25(PhysicsTermRegistry& registry);

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

    std::string getName() const override { return "DynamicVacuum_HUDF"; }
    std::string getDescription() const override {
        return "Time-varying vacuum energy for HUDF cosmic field";
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
        double M = 1.989e42;      // 1e12 M_sun in kg
        double r = 1.23e27;       // 1.3e11 ly in m
        return coupling_strength * (hbar * hbar) / (M * r * r) * cos(t / 1e6);
    }

    std::string getName() const override { return "QuantumCoupling_HUDF"; }
    std::string getDescription() const override {
        return "Non-local quantum effects for cosmic deep field";
    }
    std::string getCategory() const override { return "SelfExpanding"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"version", 2.0}, {"coupling_strength", coupling_strength}};
    }
};

// ========================================
// HUDF Galaxies Physics Terms
// ========================================

// Term 1: Base gravity with M(t), Hz, B, I(t) corrections
class HUDFBaseGravityTerm : public PhysicsTerm {
private:
    double G;       // 6.674e-11 m^3/kg·s^2
    double M0;      // 1.989e42 kg (1e12 M_sun)
    double r;       // 1.23e27 m (1.3e11 ly)
    double Hz;      // 2.5e-18 s^-1 (z=3.5)
    double B;       // 1e-10 T
    double B_crit;  // 4.4e13 T
    double SFR_factor; // 1.0
    double tau_SF;  // 3.156e16 s (1 Gyr)
    double I0;      // 0.05
    double tau_inter; // 3.156e16 s (1 Gyr)

public:
    HUDFBaseGravityTerm()
        : G(6.674e-11), M0(1.989e42), r(1.23e27), Hz(2.5e-18),
          B(1e-10), B_crit(4.4e13), SFR_factor(1.0), tau_SF(3.156e16),
          I0(0.05), tau_inter(3.156e16) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double M_t = M0 * (1 + M_dot);
        double I_t = I0 * exp(-t / tau_inter);
        double ug1_t = (G * M_t) / (r * r);
        double corr_H = 1 + Hz * t;
        double corr_B = 1 - B / B_crit;
        double corr_I = 1 + I_t;
        return ug1_t * corr_H * corr_B * corr_I;
    }

    std::string getName() const override { return "HUDF_BaseGravity"; }
    std::string getDescription() const override {
        return "Base gravity with M(t), Hz, B, I(t): g·(1+Hz·t)·(1-B/B_c)·(1+I(t))";
    }
    std::string getCategory() const override { return "CoreGravity"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"M0_kg", M0}, {"r_m", r}, {"Hz_s-1", Hz}, {"z_avg", 3.5}, 
                {"SFR_factor", SFR_factor}, {"I0", I0}};
    }
};

// Term 2: UQFF Unification with f_TRZ, M(t), I(t)
class HUDFUQFFUnificationTerm : public PhysicsTerm {
private:
    double G, M0, r;
    double f_TRZ;   // 0.01
    double B, B_crit;
    double SFR_factor, tau_SF;
    double I0, tau_inter;

public:
    HUDFUQFFUnificationTerm()
        : G(6.674e-11), M0(1.989e42), r(1.23e27), f_TRZ(0.01),
          B(1e-10), B_crit(4.4e13), SFR_factor(1.0), tau_SF(3.156e16),
          I0(0.05), tau_inter(3.156e16) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double M_t = M0 * (1 + M_dot);
        double I_t = I0 * exp(-t / tau_inter);
        double Ug1 = (G * M_t) / (r * r);
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - B / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ) * (1 + I_t);
    }

    std::string getName() const override { return "HUDF_UQFF_Unification"; }
    std::string getDescription() const override {
        return "UQFF unification Ug = (Ug1+Ug2+Ug3+Ug4)*(1+f_TRZ)*(1+I(t))";
    }
    std::string getCategory() const override { return "UQFF"; }
};

// Term 3: Cosmological constant Lambda
class HUDFCosmologicalConstantTerm : public PhysicsTerm {
private:
    double Lambda;  // 1.1056e-52 m^-2
    double c_light; // 2.998e8 m/s

public:
    HUDFCosmologicalConstantTerm()
        : Lambda(1.1056e-52), c_light(2.998e8) {}

    double compute(double t) const override {
        return (Lambda * c_light * c_light) / 3.0;
    }

    std::string getName() const override { return "HUDF_CosmologicalConstant"; }
    std::string getDescription() const override {
        return "Cosmological constant Λc²/3 acceleration";
    }
    std::string getCategory() const override { return "DarkEnergy"; }
};

// Term 4: Scaled electromagnetic with UA vacuum
class HUDFElectromagneticTerm : public PhysicsTerm {
private:
    double q_charge;     // 1.602e-19 C
    double gas_v;        // 5e5 m/s
    double B;            // 1e-10 T
    double proton_mass;  // 1.673e-27 kg
    double rho_vac_UA;   // 6.3e-10 kg/m^3
    double rho_vac_SCm;  // 1e-26 kg/m^3
    double scale_EM;     // 1e-12

public:
    HUDFElectromagneticTerm()
        : q_charge(1.602e-19), gas_v(5e5), B(1e-10), proton_mass(1.673e-27),
          rho_vac_UA(6.3e-10), rho_vac_SCm(1e-26), scale_EM(1e-12) {}

    double compute(double t) const override {
        double cross_vB = gas_v * B;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
        return (em_base * corr_UA) * scale_EM;
    }

    std::string getName() const override { return "HUDF_Electromagnetic"; }
    std::string getDescription() const override {
        return "Scaled EM with UA vacuum: (q*v×B/m_p)*(1+ρ_UA/ρ_SCm)*scale";
    }
    std::string getCategory() const override { return "Electromagnetic"; }
};

// Term 5: Quantum uncertainty
class HUDFQuantumUncertaintyTerm : public PhysicsTerm {
private:
    double hbar;          // 1.0546e-34 J·s
    double delta_x;       // 1e27 m
    double delta_p;       // 1e-20 kg·m/s
    double integral_psi;  // 1e-5
    double t_Hubble;      // 4.35e17 s

public:
    HUDFQuantumUncertaintyTerm()
        : hbar(1.0546e-34), delta_x(1e27), delta_p(1e-20),
          integral_psi(1e-5), t_Hubble(4.35e17) {}

    double compute(double t) const override {
        double sqrt_unc = sqrt(delta_x * delta_p);
        return (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);
    }

    std::string getName() const override { return "HUDF_QuantumUncertainty"; }
    std::string getDescription() const override {
        return "Quantum uncertainty: (ℏ/√(Δx·Δp))·∫|ψ|²·(2π/t_H)";
    }
    std::string getCategory() const override { return "QuantumCorrection"; }
};

// Term 6: Fluid density coupling
class HUDFFluidDensityTerm : public PhysicsTerm {
private:
    double rho_fluid; // 1e-21 kg/m^3
    double r;         // 1.23e27 m
    double G, M0;
    double SFR_factor, tau_SF;

public:
    HUDFFluidDensityTerm()
        : rho_fluid(1e-21), r(1.23e27), G(6.674e-11), M0(1.989e42),
          SFR_factor(1.0), tau_SF(3.156e16) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double M_t = M0 * (1 + M_dot);
        double ug1_t = (G * M_t) / (r * r);
        double V = (4.0 / 3.0) * M_PI * r * r * r;
        return (rho_fluid * V * ug1_t) / M_t;
    }

    std::string getName() const override { return "HUDF_FluidDensity"; }
    std::string getDescription() const override {
        return "Fluid density coupling: (ρ_fluid·V·g)/M(t)";
    }
    std::string getCategory() const override { return "FluidDynamics"; }
};

// Term 7: Oscillatory wave components
class HUDFOscillatoryWaveTerm : public PhysicsTerm {
private:
    double A_osc;       // 1e-15 m/s^2
    double k_osc;       // 1e-28 m^-1
    double omega_osc;   // 1e-15 s^-1
    double x_pos;       // 1e27 m
    double t_Hubble_gyr; // 4.35e17 s (13.8 Gyr)

public:
    HUDFOscillatoryWaveTerm()
        : A_osc(1e-15), k_osc(1e-28), omega_osc(1e-15),
          x_pos(1e27), t_Hubble_gyr(4.35e17) {}

    double compute(double t) const override {
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        return term_osc1 + term_osc2;
    }

    std::string getName() const override { return "HUDF_OscillatoryWave"; }
    std::string getDescription() const override {
        return "Oscillatory wave: 2A·cos(kx)cos(ωt) + (2π/t_H)·A·cos(kx-ωt)";
    }
    std::string getCategory() const override { return "WavePhenomena"; }
};

// Term 8: Dark matter + density perturbations with M(t)
class HUDFDarkMatterPerturbationTerm : public PhysicsTerm {
private:
    double G, M0, r;
    double M_DM_factor;         // 5.0
    double delta_rho_over_rho;  // 1e-5
    double SFR_factor, tau_SF;

public:
    HUDFDarkMatterPerturbationTerm()
        : G(6.674e-11), M0(1.989e42), r(1.23e27),
          M_DM_factor(5.0), delta_rho_over_rho(1e-5),
          SFR_factor(1.0), tau_SF(3.156e16) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double M_t = M0 * (1 + M_dot);
        double M_dm = M_t * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M_t / (r * r * r);
        double term_dm_force_like = (M_t + M_dm) * (pert1 + pert2);
        return term_dm_force_like / M_t;
    }

    std::string getName() const override { return "HUDF_DarkMatterPerturbation"; }
    std::string getDescription() const override {
        return "Dark matter + density perturbation: (M(t)+M_DM)·(δρ/ρ + 3GM(t)/r³)/M(t)";
    }
    std::string getCategory() const override { return "DarkMatter"; }
};

// Term 9: Merger feedback term (pressure/density for acceleration)
class HUDFMergerFeedbackTerm : public PhysicsTerm {
private:
    double rho_wind;  // 1e-22 kg/m^3
    double v_wind;    // 1e6 m/s
    double rho_fluid; // 1e-21 kg/m^3

public:
    HUDFMergerFeedbackTerm()
        : rho_wind(1e-22), v_wind(1e6), rho_fluid(1e-21) {}

    double compute(double t) const override {
        double wind_pressure = rho_wind * v_wind * v_wind;
        return wind_pressure / rho_fluid;
    }

    std::string getName() const override { return "HUDF_MergerFeedback"; }
    std::string getDescription() const override {
        return "Merger feedback: ρ_wind·v_wind²/ρ_fluid for acceleration";
    }
    std::string getCategory() const override { return "MergerFeedback"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"rho_wind", rho_wind}, {"v_wind_m/s", v_wind}, {"wind_pressure_Pa", rho_wind * v_wind * v_wind}};
    }
};

// Term 10: Star formation M(t) - direct contribution
class HUDFStarFormationTerm : public PhysicsTerm {
private:
    double G, M0, r;
    double SFR_factor;  // 1.0
    double tau_SF;      // 3.156e16 s (1 Gyr)

public:
    HUDFStarFormationTerm()
        : G(6.674e-11), M0(1.989e42), r(1.23e27),
          SFR_factor(1.0), tau_SF(3.156e16) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double delta_M = M0 * M_dot;
        // Mass growth increases gravity
        return (G * delta_M) / (r * r);
    }

    std::string getName() const override { return "HUDF_StarFormation"; }
    std::string getDescription() const override {
        return "Star formation M(t) = M0·(1+SFR·e^(-t/τ)) increases gravity: G·ΔM/r²";
    }
    std::string getCategory() const override { return "StarFormation"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"SFR_factor", SFR_factor}, {"tau_SF_s", tau_SF}, {"tau_SF_Gyr", tau_SF / 3.156e16}};
    }
};

// Term 11: Galaxy interaction I(t) - direct contribution
class HUDFGalaxyInteractionTerm : public PhysicsTerm {
private:
    double I0;         // 0.05
    double tau_inter;  // 3.156e16 s (1 Gyr)
    double G, M0, r;
    double SFR_factor, tau_SF;

public:
    HUDFGalaxyInteractionTerm()
        : I0(0.05), tau_inter(3.156e16),
          G(6.674e-11), M0(1.989e42), r(1.23e27),
          SFR_factor(1.0), tau_SF(3.156e16) {}

    double compute(double t) const override {
        double M_dot = SFR_factor * exp(-t / tau_SF);
        double M_t = M0 * (1 + M_dot);
        double I_t = I0 * exp(-t / tau_inter);
        double ug1_t = (G * M_t) / (r * r);
        // Interaction enhances gravity by (1+I(t))
        return ug1_t * I_t;
    }

    std::string getName() const override { return "HUDF_GalaxyInteraction"; }
    std::string getDescription() const override {
        return "Galaxy interaction I(t) = I0·e^(-t/τ) enhances gravity: g·I(t)";
    }
    std::string getCategory() const override { return "GalaxyInteraction"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"I0", I0}, {"tau_inter_s", tau_inter}, {"tau_inter_Gyr", tau_inter / 3.156e16}};
    }
};

// ========================================
// Registration Function
// ========================================

void registerWolframTerms_source26(PhysicsTermRegistry& registry) {
    // First, inherit all 241 terms from source25 (NGC 1275)
    registerWolframTerms_source25(registry);

    // Register 2 self-expanding framework terms
    registry.registerTerm(std::make_unique<DynamicVacuumTerm>());
    registry.registerTerm(std::make_unique<QuantumCouplingTerm>());

    // Register 11 HUDF Galaxies physics terms
    registry.registerTerm(std::make_unique<HUDFBaseGravityTerm>());
    registry.registerTerm(std::make_unique<HUDFUQFFUnificationTerm>());
    registry.registerTerm(std::make_unique<HUDFCosmologicalConstantTerm>());
    registry.registerTerm(std::make_unique<HUDFElectromagneticTerm>());
    registry.registerTerm(std::make_unique<HUDFQuantumUncertaintyTerm>());
    registry.registerTerm(std::make_unique<HUDFFluidDensityTerm>());
    registry.registerTerm(std::make_unique<HUDFOscillatoryWaveTerm>());
    registry.registerTerm(std::make_unique<HUDFDarkMatterPerturbationTerm>());
    registry.registerTerm(std::make_unique<HUDFMergerFeedbackTerm>());
    registry.registerTerm(std::make_unique<HUDFStarFormationTerm>());
    registry.registerTerm(std::make_unique<HUDFGalaxyInteractionTerm>());

    // Total: 241 (inherited) + 13 (new) = 254 classes
}
