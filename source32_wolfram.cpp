// source32_wolfram.cpp - Wolfram Companion for Crab Nebula UQFF Module
// Auto-generated wolfram integration terms for pulsar-driven supernova remnant
// System: Crab Nebula - M1 in Taurus
// Mass: M=4.6 M_sun = 9.149e30 kg (ejected stellar envelope + pulsar)
// Radius: r(t)=r0+v_exp·t, r0=5.2e16 m (5.5 ly initial), v_exp=1.5e6 m/s (expansion velocity)
// Age: t=971 years since 1054 AD supernova
// Current radius (t=971 yr): r≈9.8e16 m (10.4 ly)
// Redshift: z=0.0015 (local, 6500 ly)
// Pulsar: P_pulsar=5e31 W (relativistic wind), period=33 ms (PSR B0531+21)
// Magnetic field: B=1e-8 T (nebula average, synchrotron radiation)
// Shock velocity: v_shock=1.5e6 m/s (expansion rate)
// Filament density: ρ_fluid=1e-21 kg/m³
// Superconductivity: (1-B/B_crit), B_crit=1e11 T
// DM fraction: 0% (M_DM=0, visible mass only)
// Unique physics: r(t) time-dependent expansion, pulsar wind pressure, magnetic Lorentz force
//
// Copyright - Daniel T. Murphy, Wolfram integration Oct 09, 2025
// Part of Star-Magic UQFF framework - 317 + 13 = 330 total PhysicsTerm classes

#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <memory>
#include <complex>

// ===========================================================================================
// PHYSICS TERM BASE CLASS (Shared across all wolfram companions)
// ===========================================================================================

class PhysicsTerm {
public:
    virtual ~PhysicsTerm() {}
    virtual double compute(double t) const = 0;
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
    virtual std::string getCategory() const = 0;
    virtual std::map<std::string, double> getMetadata() const { return {}; }
};

// ===========================================================================================
// SELF-EXPANDING FRAMEWORK TERMS (2 universal terms in all wolfram companions)
// ===========================================================================================

class DynamicVacuumTerm : public PhysicsTerm {
private:
    double amplitude;
    double frequency;
    double rho_vac_UA;
public:
    DynamicVacuumTerm(double amp = 1e-10, double freq = 1e-15, double rho = 7.09e-36)
        : amplitude(amp), frequency(freq), rho_vac_UA(rho) {}
    
    double compute(double t) const override {
        return amplitude * rho_vac_UA * std::sin(frequency * t);
    }
    
    std::string getName() const override { return "DynamicVacuum"; }
    std::string getDescription() const override { 
        return "Time-varying vacuum energy density oscillation"; 
    }
    std::string getCategory() const override { return "vacuum"; }
};

class QuantumCouplingTerm : public PhysicsTerm {
private:
    double coupling_strength;
    double hbar;
    double M;
    double r0, v_exp;
public:
    QuantumCouplingTerm(double strength = 1e-40, double mass = 9.149e30, 
                        double radius_init = 5.2e16, double v = 1.5e6)
        : coupling_strength(strength), hbar(1.0546e-34), M(mass), 
          r0(radius_init), v_exp(v) {}
    
    double compute(double t) const override {
        double r = r0 + v_exp * t;
        return coupling_strength * (hbar * hbar) / (M * r * r) * std::cos(t / 1e6);
    }
    
    std::string getName() const override { return "QuantumCoupling"; }
    std::string getDescription() const override { 
        return "Non-local quantum entanglement effects with r(t) evolution"; 
    }
    std::string getCategory() const override { return "quantum"; }
};

// ===========================================================================================
// CRAB NEBULA PHYSICS TERMS (11 system-specific terms)
// ===========================================================================================

// 1. Base Gravity with Time-Dependent Radius r(t)=r0+v_exp·t
class CrabBaseGravityTerm : public PhysicsTerm {
private:
    double G, M, r0, v_exp;
    double H0, Omega_m, Omega_Lambda, z, Mpc_to_m;
    double f_TRZ, B, B_crit;
public:
    CrabBaseGravityTerm()
        : G(6.6743e-11), M(9.149e30), r0(5.2e16), v_exp(1.5e6),
          H0(70.0), Omega_m(0.3), Omega_Lambda(0.7), z(0.0015),
          Mpc_to_m(3.086e22), f_TRZ(0.1), B(1e-8), B_crit(1e11) {}
    
    double compute(double t) const override {
        // Time-dependent radius: r(t) = r0 + v_exp·t
        double r = r0 + v_exp * t;
        
        // H(z) evolution: Hz = H0·sqrt(Omega_m·(1+z)³ + Omega_Lambda) in s⁻¹
        double Hz_kms = H0 * std::sqrt(Omega_m * std::pow(1.0 + z, 3) + Omega_Lambda);
        double Hz = (Hz_kms * 1e3) / Mpc_to_m;
        
        // Base gravity with expansion, superconductivity, time-reversal
        double expansion = 1.0 + Hz * t;
        double sc_correction = 1.0 - (B / B_crit);
        double tr_factor = 1.0 + f_TRZ;
        
        return (G * M / (r * r)) * expansion * sc_correction * tr_factor;
    }
    
    std::string getName() const override { return "CrabBaseGravity"; }
    std::string getDescription() const override {
        return "Base Newtonian gravity with time-dependent radius r(t)=r0+v_exp·t, "
               "H(z) expansion, superconductivity (1-B/B_crit), time-reversal correction";
    }
    std::string getCategory() const override { return "gravity"; }
    std::map<std::string, double> getMetadata() const override {
        return {{"M", M}, {"r0", r0}, {"v_exp", v_exp}, {"age_yr", 971}, {"z", z}};
    }
};

// 2. Pulsar Wind Pressure (Dominant outward force)
class CrabPulsarWindTerm : public PhysicsTerm {
private:
    double P_pulsar, pi, v_shock, c, rho_fluid, scale_macro;
    double r0, v_exp;
public:
    CrabPulsarWindTerm()
        : P_pulsar(5e31), pi(3.141592653589793), v_shock(1.5e6),
          c(3e8), rho_fluid(1e-21), scale_macro(1e-12),
          r0(5.2e16), v_exp(1.5e6) {}
    
    double compute(double t) const override {
        double r = r0 + v_exp * t;
        // Wind pressure: P = (L_pulsar / 4πr²) · (1 + v_shock/c)
        double pressure = (P_pulsar / (4 * pi * r * r)) * (1.0 + v_shock / c);
        // Acceleration: a = P/ρ
        return (pressure / rho_fluid) * scale_macro;
    }
    
    std::string getName() const override { return "CrabPulsarWind"; }
    std::string getDescription() const override {
        return "Relativistic pulsar wind pressure acceleration: "
               "a_wind = [P_pulsar/(4πr²)·(1+v_shock/c)]/ρ·scale with P=5e31 W";
    }
    std::string getCategory() const override { return "pulsarwind"; }
};

// 3. Magnetic Lorentz Force (Synchrotron radiation electrons)
class CrabMagneticLorentzTerm : public PhysicsTerm {
private:
    double q, v_shock, B, m_e, scale_macro;
public:
    CrabMagneticLorentzTerm()
        : q(1.602e-19), v_shock(1.5e6), B(1e-8),
          m_e(9.11e-31), scale_macro(1e-12) {}
    
    double compute(double t) const override {
        // Lorentz force on electrons: F = q·v×B
        // Acceleration: a = F/m_e = (q·v·B)/m_e
        return (q * v_shock * B / m_e) * scale_macro;
    }
    
    std::string getName() const override { return "CrabMagneticLorentz"; }
    std::string getDescription() const override {
        return "Lorentz force on relativistic electrons in nebula magnetic field: "
               "a_mag = (q·v_shock·B/m_e)·scale with B=1e-8 T (synchrotron source)";
    }
    std::string getCategory() const override { return "magnetic"; }
};

// 4. UQFF Unification with r(t) (Ug1 + Ug4 only, Ug2/Ug3=0 for remnant)
class CrabUQFFUnificationTerm : public PhysicsTerm {
private:
    double G, M, r0, v_exp, f_sc;
public:
    CrabUQFFUnificationTerm()
        : G(6.6743e-11), M(9.149e30), r0(5.2e16), v_exp(1.5e6), f_sc(1.0) {}
    
    double compute(double t) const override {
        double r = r0 + v_exp * t;
        double Ug1 = (G * M) / (r * r);
        double Ug4 = Ug1 * f_sc;
        return Ug1 + Ug4;  // Ug2=Ug3=0 (no d²Φ/dt², no moon)
    }
    
    std::string getName() const override { return "CrabUQFFUnification"; }
    std::string getDescription() const override {
        return "UQFF subterms with r(t): Ug1 + Ug4 = G·M/r(t)²·(1+f_sc), Ug2/Ug3=0";
    }
    std::string getCategory() const override { return "unified"; }
};

// 5. Cosmological Constant (Lambda c²/3)
class CrabCosmologicalConstantTerm : public PhysicsTerm {
private:
    double Lambda, c;
public:
    CrabCosmologicalConstantTerm()
        : Lambda(1.1e-52), c(3e8) {}
    
    double compute(double t) const override {
        return Lambda * (c * c) / 3.0;
    }
    
    std::string getName() const override { return "CrabCosmologicalConstant"; }
    std::string getDescription() const override {
        return "Dark energy contribution: Λc²/3 with Λ=1.1e-52 m⁻²";
    }
    std::string getCategory() const override { return "vacuum"; }
};

// 6. Quantum Uncertainty (Heisenberg for particle quantum effects)
class CrabQuantumUncertaintyTerm : public PhysicsTerm {
private:
    double hbar, Delta_x, pi, t_Hubble;
public:
    CrabQuantumUncertaintyTerm()
        : hbar(1.0546e-34), Delta_x(1e-10), pi(3.141592653589793),
          t_Hubble(13.8e9 * 3.156e7) {}
    
    double compute(double t) const override {
        double Delta_p = hbar / Delta_x;  // Heisenberg
        double unc = std::sqrt(Delta_x * Delta_p);
        // Quantum term: (ℏ/√(Δx·Δp))·∫|ψ|²·(2π/t_H) with integral≈1 (ground state)
        return (hbar / unc) * 1.0 * (2 * pi / t_Hubble);
    }
    
    std::string getName() const override { return "CrabQuantumUncertainty"; }
    std::string getDescription() const override {
        return "Quantum uncertainty for nebular particles: "
               "(ℏ/√(Δx·Δp))·<ψ|H|ψ>·(2π/t_H) with Δx=1e-10 m";
    }
    std::string getCategory() const override { return "quantum"; }
};

// 7. Electromagnetic Lorentz Force (q·v×B for shock physics)
class CrabElectromagneticTerm : public PhysicsTerm {
private:
    double q, v_shock, B, m_proton, rho_UA, rho_SCm, scale_macro;
public:
    CrabElectromagneticTerm()
        : q(1.602e-19), v_shock(1.5e6), B(1e-8), m_proton(1.673e-27),
          rho_UA(7.09e-36), rho_SCm(7.09e-37), scale_macro(1e-12) {}
    
    double compute(double t) const override {
        double em_base = (q * v_shock * B) / m_proton;
        double ua_scm_ratio = 1.0 + (rho_UA / rho_SCm);  // ≈10
        return em_base * ua_scm_ratio * scale_macro;
    }
    
    std::string getName() const override { return "CrabElectromagnetic"; }
    std::string getDescription() const override {
        return "Lorentz force on shock-accelerated ions: (q·v_shock·B/m_p)·(1+ρ_UA/ρ_SCm)·scale "
               "with v_shock=1.5e6 m/s (expansion/shock velocity)";
    }
    std::string getCategory() const override { return "electromagnetic"; }
};

// 8. Fluid Density (Nebular filament ISM)
class CrabFluidDensityTerm : public PhysicsTerm {
private:
    double rho_fluid, V, G, M, r0, v_exp;
public:
    CrabFluidDensityTerm()
        : rho_fluid(1e-21), V(1e3), G(6.6743e-11), M(9.149e30),
          r0(5.2e16), v_exp(1.5e6) {}
    
    double compute(double t) const override {
        double r = r0 + v_exp * t;
        double g_base = (G * M) / (r * r);
        return rho_fluid * V * g_base;
    }
    
    std::string getName() const override { return "CrabFluidDensity"; }
    std::string getDescription() const override {
        return "Nebular filament density-volume-gravity coupling: ρ_fluid·V·g(r(t)) with ρ=1e-21 kg/m³";
    }
    std::string getCategory() const override { return "fluid"; }
};

// 9. Oscillatory Wave (Aether-mediated wisp dynamics)
class CrabOscillatoryWaveTerm : public PhysicsTerm {
private:
    double A, k, omega, x, pi;
public:
    CrabOscillatoryWaveTerm()
        : A(1e-10), k(1e20), omega(1e15), x(0.0), pi(3.141592653589793) {}
    
    double compute(double t) const override {
        // Resonant wave: 2A·cos(kx)cos(ωt) + (2π/13.8)·A·Re[e^(i(kx-ωt))]
        double cos_term = 2 * A * std::cos(k * x) * std::cos(omega * t);
        std::complex<double> exp_i(0, k * x - omega * t);
        double real_exp = A * std::exp(exp_i).real();
        double exp_factor = (2 * pi / 13.8);  // Gyr normalization
        return cos_term + exp_factor * real_exp;
    }
    
    std::string getName() const override { return "CrabOscillatoryWave"; }
    std::string getDescription() const override {
        return "Aether-mediated resonant wave for wisp dynamics: "
               "2A·cos(kx)cos(ωt) + (2π/13.8)·A·Re[e^(i(kx-ωt))] (synchrotron wisps)";
    }
    std::string getCategory() const override { return "wave"; }
};

// 10. Dark Matter Perturbation with r(t) (M_DM=0, visible mass only)
class CrabDarkMatterPerturbationTerm : public PhysicsTerm {
private:
    double M_visible, M_DM, delta_rho, rho, G, M, r0, v_exp;
public:
    CrabDarkMatterPerturbationTerm()
        : M_visible(9.149e30), M_DM(0.0), delta_rho(1e-22), rho(1e-21),
          G(6.6743e-11), M(9.149e30), r0(5.2e16), v_exp(1.5e6) {}
    
    double compute(double t) const override {
        double r = r0 + v_exp * t;
        double pert = delta_rho / rho;  // ≈0.1
        double curv = 3 * G * M / (r * r * r);
        return (M_visible + M_DM) * (pert + curv);
    }
    
    std::string getName() const override { return "CrabDarkMatterPerturbation"; }
    std::string getDescription() const override {
        return "Visible mass perturbation with r(t) curvature (M_DM=0): "
               "M_vis·(δρ/ρ + 3GM/r(t)³)";
    }
    std::string getCategory() const override { return "darkmatter"; }
};

// 11. Superconductivity Correction with r(t) (Quantum field effect near pulsar)
class CrabSuperconductivityTerm : public PhysicsTerm {
private:
    double G, M, r0, v_exp, B, B_crit;
public:
    CrabSuperconductivityTerm()
        : G(6.6743e-11), M(9.149e30), r0(5.2e16), v_exp(1.5e6),
          B(1e-8), B_crit(1e11) {}
    
    double compute(double t) const override {
        double r = r0 + v_exp * t;
        double g_base = (G * M) / (r * r);
        // Superconductivity reduces gravity: -g·(B/B_crit)
        return -g_base * (B / B_crit);
    }
    
    std::string getName() const override { return "CrabSuperconductivity"; }
    std::string getDescription() const override {
        return "Quantum field superconductivity correction with r(t): -g(r(t))·(B/B_crit) "
               "with B=1e-8 T, B_crit=1e11 T";
    }
    std::string getCategory() const override { return "superconductivity"; }
};

// ===========================================================================================
// WOLFRAM REGISTRATION FUNCTION (Delegation chain: source31 → source32)
// ===========================================================================================

// Forward declaration of previous registration (source31_wolfram.cpp)
void registerWolframTerms_source31(std::vector<std::unique_ptr<PhysicsTerm>>& terms);

void registerWolframTerms_source32(std::vector<std::unique_ptr<PhysicsTerm>>& terms) {
    // Inherit all previous terms from source31 (317 classes)
    registerWolframTerms_source31(terms);
    
    // Add source32 self-expanding framework terms (2)
    terms.push_back(std::make_unique<DynamicVacuumTerm>());
    terms.push_back(std::make_unique<QuantumCouplingTerm>());
    
    // Add Crab Nebula physics terms (11)
    terms.push_back(std::make_unique<CrabBaseGravityTerm>());
    terms.push_back(std::make_unique<CrabPulsarWindTerm>());
    terms.push_back(std::make_unique<CrabMagneticLorentzTerm>());
    terms.push_back(std::make_unique<CrabUQFFUnificationTerm>());
    terms.push_back(std::make_unique<CrabCosmologicalConstantTerm>());
    terms.push_back(std::make_unique<CrabQuantumUncertaintyTerm>());
    terms.push_back(std::make_unique<CrabElectromagneticTerm>());
    terms.push_back(std::make_unique<CrabFluidDensityTerm>());
    terms.push_back(std::make_unique<CrabOscillatoryWaveTerm>());
    terms.push_back(std::make_unique<CrabDarkMatterPerturbationTerm>());
    terms.push_back(std::make_unique<CrabSuperconductivityTerm>());
    
    // Total: 317 (source31) + 2 (framework) + 11 (Crab) = 330 PhysicsTerm classes
}

// ===========================================================================================
// PHYSICS SUMMARY FOR CRAB NEBULA
// ===========================================================================================
/*
Crab Nebula (M1) - Pulsar-Driven Supernova Remnant with Time-Dependent Expansion

SYSTEM PARAMETERS:
- Mass: M=4.6 M_sun = 9.149e30 kg (ejected stellar envelope + pulsar)
- Initial radius: r0=5.2e16 m (5.5 light-years at discovery)
- Expansion velocity: v_exp=1.5e6 m/s (0.5% speed of light)
- Age: t=971 years (since 1054 AD supernova)
- Current radius: r(t=971 yr) ≈ 9.8e16 m (10.4 ly)
- Redshift: z=0.0015 (local, distance ≈6500 ly in Taurus)
- Pulsar: PSR B0531+21, P_pulsar=5e31 W, rotation period=33 ms
- Magnetic field: B=1e-8 T (nebula average, drives synchrotron emission)
- Shock velocity: v_shock=1.5e6 m/s (equals expansion rate)
- Filament density: ρ_fluid=1e-21 kg/m³

UNIQUE PHYSICS:
1. **Time-Dependent Radius**: r(t) = r0 + v_exp·t
   - Radius grows linearly with time due to supersonic expansion
   - All gravity terms scale as 1/r(t)², weakening with age
   - At t=971 yr: r≈9.8e16 m (88% larger than initial)
   
2. **Pulsar Wind Pressure**: a_wind = [P_pulsar/(4πr²)·(1+v_shock/c)]/ρ
   - Relativistic wind from rapidly rotating neutron star
   - Dominant outward force: a_wind ≈ 1.48e6 m/s² (overwhelms gravity)
   - Drives nebula expansion and synchrotron emission
   - Decreases as 1/r² as nebula expands
   
3. **Magnetic Lorentz Force**: a_mag = (q·v_shock·B/m_e)·scale
   - Acts on relativistic electrons in nebula magnetic field
   - Produces synchrotron radiation (optical/X-ray wisps)
   - B=1e-8 T average (locally much stronger near pulsar)
   
4. **Wisp Dynamics**: Oscillatory waves from moving magnetic knots
   - Aether-mediated resonant waves: 2A·cos(kx)cos(ωt) + complex exp
   - Visible as bright wisps moving outward at ~c/2

OBSERVATIONAL CONTEXT:
- SN 1054 observed by Chinese/Arab astronomers (guest star visible in daylight)
- First recognized supernova remnant associated with historical event
- Central pulsar discovered 1968 (first optical pulsar)
- Synchrotron nebula powered by pulsar wind (luminosity 75,000 L_sun)
- Famous wisp features visible in time-lapse Hubble images
- X-ray/gamma-ray source from high-energy electrons

GRAVITY TERMS (11 total):
1. BaseGravity: (G·M/r(t)²)·(1+Hz·t)·(1-B/B_c)·(1+f_TRZ) with r(t) expansion
2. PulsarWind: [P_pulsar/(4πr²)·(1+v_shock/c)]/ρ·scale (DOMINANT outward force)
3. MagneticLorentz: (q·v_shock·B/m_e)·scale (electron acceleration)
4. UQFFUnification: Ug1 + Ug4 = G·M/r(t)²·(1+f_sc) with r(t)
5. CosmologicalConstant: Λc²/3 (dark energy)
6. QuantumUncertainty: (ℏ/√(Δx·Δp))·<ψ|H|ψ>·(2π/t_H) (particle quantum)
7. Electromagnetic: (q·v_shock·B/m_p)·(1+ρ_UA/ρ_SCm)·scale (shock ions)
8. FluidDensity: ρ_fluid·V·g(r(t)) (filament ISM)
9. OscillatoryWave: 2A·cos(kx)cos(ωt) + Aether wave (wisp oscillations)
10. DarkMatterPerturbation: M_vis·(δρ/ρ + 3GM/r(t)³) (no DM)
11. Superconductivity: -g(r(t))·(B/B_crit) (quantum field with r(t))

PHYSICS CATEGORIES:
- gravity, pulsarwind (NEW), magnetic (NEW for electrons), unified, vacuum, quantum,
  electromagnetic, fluid, wave, darkmatter, superconductivity

NEW CATEGORIES INTRODUCED:
- pulsarwind: Relativistic pulsar wind pressure
- magnetic: Magnetic Lorentz force (distinct from EM term which uses proton mass)

TOTAL ACCUMULATED CLASSES: 330
- source21-31: 317 classes
- source32 framework: 2 classes (DynamicVacuum, QuantumCoupling)
- source32 Crab: 11 classes
- Delegation: source31 → source32 → [future source33]

Copyright - Daniel T. Murphy, Wolfram integration analyzed Oct 09, 2025
*/
