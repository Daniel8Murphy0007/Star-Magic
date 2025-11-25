// source24_wolfram.cpp
// Wolfram companion file for Horsehead Nebula (Barnard 33)
// Extracts all physics terms from source24.cpp into PhysicsTerm classes
// Astronomical system: Horsehead Nebula - Barnard 33 (M0=1000 M_sun, r=2.365e16 m, B=1e-6 T)
// Key: Dark nebula with erosion factor E(t) = E_0·(1-e^(-t/τ_erosion))

#include "PhysicsTerm.h"
#include "PhysicsTermRegistry.h"
#include <cmath>
#include <memory>

// Forward declare delegation
void registerWolframTerms_source23(PhysicsTermRegistry& registry);

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

    std::string getName() const override { return "DynamicVacuum_Horsehead"; }
    std::string getDescription() const override {
        return "Time-varying vacuum energy for Horsehead Nebula dark cloud";
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
        double M = 1.989e33;      // 1000 M_sun in kg
        double r = 2.365e16;      // 2.5 ly in m
        return coupling_strength * (hbar * hbar) / (M * r * r) * cos(t / 1e6);
    }

    std::string getName() const override { return "QuantumCoupling_Horsehead"; }
    std::string getDescription() const override {
        return "Non-local quantum effects for dark nebula";
    }
    std::string getCategory() const override { return "SelfExpanding"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"version", 2.0}, {"coupling_strength", coupling_strength}};
    }
};

// ========================================
// Horsehead Nebula Physics Terms
// ========================================

// Term 1: Base gravity with H0, B, E(t) erosion corrections
class HorseheadBaseGravityTerm : public PhysicsTerm {
private:
    double G;       // 6.674e-11 m^3/kg·s^2
    double M;       // 1.989e33 kg (1000 M_sun)
    double r;       // 2.365e16 m (2.5 ly)
    double H0;      // 2.3e-18 s^-1
    double B;       // 1e-6 T
    double B_crit;  // 4.4e13 T
    double E_0;     // 0.1 (erosion factor)
    double tau_erosion; // 1.578e14 s (5 Myr)

public:
    HorseheadBaseGravityTerm()
        : G(6.674e-11), M(1.989e33), r(2.365e16), H0(2.3e-18),
          B(1e-6), B_crit(4.4e13), E_0(0.1), tau_erosion(1.578e14) {}

    double compute(double t) const override {
        double E_t = E_0 * (1 - exp(-t / tau_erosion));
        double ug1 = (G * M) / (r * r);
        double corr_H = 1 + H0 * t;
        double corr_B = 1 - B / B_crit;
        double corr_E = 1 - E_t;
        return ug1 * corr_H * corr_B * corr_E;
    }

    std::string getName() const override { return "Horsehead_BaseGravity"; }
    std::string getDescription() const override {
        return "Base gravity with H0, B, and E(t) erosion: g·(1+H0·t)·(1-B/B_c)·(1-E(t))";
    }
    std::string getCategory() const override { return "CoreGravity"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"M_kg", M}, {"r_m", r}, {"H0_s-1", H0}, {"E_0", E_0}, {"tau_erosion_s", tau_erosion}};
    }
};

// Term 2: UQFF Unification with f_TRZ and E(t)
class HorseheadUQFFUnificationTerm : public PhysicsTerm {
private:
    double G, M, r;
    double f_TRZ;   // 0.01
    double B, B_crit;
    double E_0, tau_erosion;

public:
    HorseheadUQFFUnificationTerm()
        : G(6.674e-11), M(1.989e33), r(2.365e16), f_TRZ(0.01),
          B(1e-6), B_crit(4.4e13), E_0(0.1), tau_erosion(1.578e14) {}

    double compute(double t) const override {
        double E_t = E_0 * (1 - exp(-t / tau_erosion));
        double Ug1 = (G * M) / (r * r);
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - B / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ) * (1 - E_t);
    }

    std::string getName() const override { return "Horsehead_UQFF_Unification"; }
    std::string getDescription() const override {
        return "UQFF unification Ug = (Ug1+Ug2+Ug3+Ug4)*(1+f_TRZ)*(1-E(t))";
    }
    std::string getCategory() const override { return "UQFF"; }
};

// Term 3: Cosmological constant Lambda
class HorseheadCosmologicalConstantTerm : public PhysicsTerm {
private:
    double Lambda;  // 1.1056e-52 m^-2
    double c_light; // 2.998e8 m/s

public:
    HorseheadCosmologicalConstantTerm()
        : Lambda(1.1056e-52), c_light(2.998e8) {}

    double compute(double t) const override {
        return (Lambda * c_light * c_light) / 3.0;
    }

    std::string getName() const override { return "Horsehead_CosmologicalConstant"; }
    std::string getDescription() const override {
        return "Cosmological constant Λc²/3 acceleration";
    }
    std::string getCategory() const override { return "DarkEnergy"; }
};

// Term 4: Scaled electromagnetic with UA vacuum
class HorseheadElectromagneticTerm : public PhysicsTerm {
private:
    double q_charge;     // 1.602e-19 C
    double gas_v;        // 5e5 m/s
    double B;            // 1e-6 T
    double proton_mass;  // 1.673e-27 kg
    double rho_vac_UA;   // 6.3e-10 kg/m^3
    double rho_vac_SCm;  // 1e-26 kg/m^3
    double scale_EM;     // 1e-12

public:
    HorseheadElectromagneticTerm()
        : q_charge(1.602e-19), gas_v(5e5), B(1e-6), proton_mass(1.673e-27),
          rho_vac_UA(6.3e-10), rho_vac_SCm(1e-26), scale_EM(1e-12) {}

    double compute(double t) const override {
        double cross_vB = gas_v * B;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
        return (em_base * corr_UA) * scale_EM;
    }

    std::string getName() const override { return "Horsehead_Electromagnetic"; }
    std::string getDescription() const override {
        return "Scaled EM with UA vacuum: (q*v×B/m_p)*(1+ρ_UA/ρ_SCm)*scale";
    }
    std::string getCategory() const override { return "Electromagnetic"; }
};

// Term 5: Quantum uncertainty
class HorseheadQuantumUncertaintyTerm : public PhysicsTerm {
private:
    double hbar;          // 1.0546e-34 J·s
    double delta_x;       // 1e16 m
    double delta_p;       // 1e-20 kg·m/s
    double integral_psi;  // 1e-5
    double t_Hubble;      // 4.35e17 s

public:
    HorseheadQuantumUncertaintyTerm()
        : hbar(1.0546e-34), delta_x(1e16), delta_p(1e-20),
          integral_psi(1e-5), t_Hubble(4.35e17) {}

    double compute(double t) const override {
        double sqrt_unc = sqrt(delta_x * delta_p);
        return (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);
    }

    std::string getName() const override { return "Horsehead_QuantumUncertainty"; }
    std::string getDescription() const override {
        return "Quantum uncertainty: (ℏ/√(Δx·Δp))·∫|ψ|²·(2π/t_H)";
    }
    std::string getCategory() const override { return "QuantumCorrection"; }
};

// Term 6: Fluid density coupling
class HorseheadFluidDensityTerm : public PhysicsTerm {
private:
    double rho_fluid; // 1e-21 kg/m^3
    double r;         // 2.365e16 m
    double G, M;

public:
    HorseheadFluidDensityTerm()
        : rho_fluid(1e-21), r(2.365e16), G(6.674e-11), M(1.989e33) {}

    double compute(double t) const override {
        double ug1 = (G * M) / (r * r);
        double V = (4.0 / 3.0) * M_PI * r * r * r;
        return (rho_fluid * V * ug1) / M;
    }

    std::string getName() const override { return "Horsehead_FluidDensity"; }
    std::string getDescription() const override {
        return "Fluid density coupling: (ρ_fluid·V·g)/M";
    }
    std::string getCategory() const override { return "FluidDynamics"; }
};

// Term 7: Oscillatory wave components
class HorseheadOscillatoryWaveTerm : public PhysicsTerm {
private:
    double A_osc;       // 1e-15 m/s^2
    double k_osc;       // 1e-17 m^-1
    double omega_osc;   // 1e-15 s^-1
    double x_pos;       // 1e16 m
    double t_Hubble_gyr; // 4.35e17 s (13.8 Gyr)

public:
    HorseheadOscillatoryWaveTerm()
        : A_osc(1e-15), k_osc(1e-17), omega_osc(1e-15),
          x_pos(1e16), t_Hubble_gyr(4.35e17) {}

    double compute(double t) const override {
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        return term_osc1 + term_osc2;
    }

    std::string getName() const override { return "Horsehead_OscillatoryWave"; }
    std::string getDescription() const override {
        return "Oscillatory wave: 2A·cos(kx)cos(ωt) + (2π/t_H)·A·cos(kx-ωt)";
    }
    std::string getCategory() const override { return "WavePhenomena"; }
};

// Term 8: Dark matter + density perturbations
class HorseheadDarkMatterPerturbationTerm : public PhysicsTerm {
private:
    double G, M, r;
    double M_DM_factor;         // 5.0
    double delta_rho_over_rho;  // 1e-5

public:
    HorseheadDarkMatterPerturbationTerm()
        : G(6.674e-11), M(1.989e33), r(2.365e16),
          M_DM_factor(5.0), delta_rho_over_rho(1e-5) {}

    double compute(double t) const override {
        double M_dm = M * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M / (r * r * r);
        double term_dm_force_like = (M + M_dm) * (pert1 + pert2);
        return term_dm_force_like / M;
    }

    std::string getName() const override { return "Horsehead_DarkMatterPerturbation"; }
    std::string getDescription() const override {
        return "Dark matter + density perturbation: (M+M_DM)·(δρ/ρ + 3GM/r³)/M";
    }
    std::string getCategory() const override { return "DarkMatter"; }
};

// Term 9: Stellar wind feedback (pressure/density for acceleration)
class HorseheadStellarWindTerm : public PhysicsTerm {
private:
    double rho_wind;   // 1e-21 kg/m^3
    double v_wind;     // 2e6 m/s
    double rho_fluid;  // 1e-21 kg/m^3

public:
    HorseheadStellarWindTerm()
        : rho_wind(1e-21), v_wind(2e6), rho_fluid(1e-21) {}

    double compute(double t) const override {
        double wind_pressure = rho_wind * v_wind * v_wind;
        return wind_pressure / rho_fluid;
    }

    std::string getName() const override { return "Horsehead_StellarWind"; }
    std::string getDescription() const override {
        return "Stellar wind feedback: ρ_wind·v_wind²/ρ_fluid for acceleration";
    }
    std::string getCategory() const override { return "StellarFeedback"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"rho_wind", rho_wind}, {"v_wind_m/s", v_wind}, {"wind_pressure_Pa", rho_wind * v_wind * v_wind}};
    }
};

// Term 10: Erosion Factor E(t) - direct contribution
class HorseheadErosionTerm : public PhysicsTerm {
private:
    double E_0;          // 0.1
    double tau_erosion;  // 1.578e14 s (5 Myr)
    double G, M, r;

public:
    HorseheadErosionTerm()
        : E_0(0.1), tau_erosion(1.578e14),
          G(6.674e-11), M(1.989e33), r(2.365e16) {}

    double compute(double t) const override {
        double E_t = E_0 * (1 - exp(-t / tau_erosion));
        double ug1 = (G * M) / (r * r);
        // Erosion reduces gravity by factor (1-E(t))
        return -ug1 * E_t;
    }

    std::string getName() const override { return "Horsehead_ErosionFactor"; }
    std::string getDescription() const override {
        return "Erosion E(t) = E_0·(1-e^(-t/τ)) reduces gravity: -g·E(t)";
    }
    std::string getCategory() const override { return "ErosionDynamics"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"E_0", E_0}, {"tau_erosion_s", tau_erosion}, {"tau_erosion_Myr", tau_erosion / 3.156e13}};
    }
};

// ========================================
// Registration Function
// ========================================

void registerWolframTerms_source24(PhysicsTermRegistry& registry) {
    // First, inherit all 215 terms from source23 (Antennae Galaxies)
    registerWolframTerms_source23(registry);

    // Register 2 self-expanding framework terms
    registry.registerTerm(std::make_unique<DynamicVacuumTerm>());
    registry.registerTerm(std::make_unique<QuantumCouplingTerm>());

    // Register 10 Horsehead Nebula physics terms
    registry.registerTerm(std::make_unique<HorseheadBaseGravityTerm>());
    registry.registerTerm(std::make_unique<HorseheadUQFFUnificationTerm>());
    registry.registerTerm(std::make_unique<HorseheadCosmologicalConstantTerm>());
    registry.registerTerm(std::make_unique<HorseheadElectromagneticTerm>());
    registry.registerTerm(std::make_unique<HorseheadQuantumUncertaintyTerm>());
    registry.registerTerm(std::make_unique<HorseheadFluidDensityTerm>());
    registry.registerTerm(std::make_unique<HorseheadOscillatoryWaveTerm>());
    registry.registerTerm(std::make_unique<HorseheadDarkMatterPerturbationTerm>());
    registry.registerTerm(std::make_unique<HorseheadStellarWindTerm>());
    registry.registerTerm(std::make_unique<HorseheadErosionTerm>());

    // Total: 215 (inherited) + 12 (new) = 227 classes
}
