// source45_wolfram.cpp - Spiral Galaxies & Supernovae UQFF Wolfram Companion
// Wolfram Language integration for SpiralSupernovaeUQFFModule (Spiral Galaxy Evolution + SN)
// Spiral galaxy: M=1.989e41 kg (~1e11 M_sun, typical spiral), r=9.258e20 m (~30 kpc, 98,000 ly)
// Rotation curve: Ω_p~20 km/s/kpc (pattern speed), v_rot~200 km/s (flat rotation curve)
// Supernova: L_SN=1e36 W (Type Ia/II luminosity), H0=73 km/s/Mpc, z up to 1.5 (cosmological)
// Physics: Spiral density waves T_spiral, supernova SN_term, dark matter rotation curve, Lambda(z)
// UQFF replaces SM gravity with spiral wave + SN-driven dynamics
// Copyright - Daniel T. Murphy, analyzed Oct 09, 2025

#include <cmath>
#include <complex>
#include <string>
#include <map>
#include <memory>
#include <vector>

// ===========================================================================================
// WOLFRAM PHYSICS TERM CLASSES (PhysicsTerm Interface)
// ===========================================================================================

class PhysicsTerm {
public:
    virtual ~PhysicsTerm() {}
    virtual double compute(double t) const = 0;
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
    virtual std::string getCategory() const = 0;
};

// Framework terms (inherited from source44)
class DynamicVacuumTerm : public PhysicsTerm {
private:
    double amplitude, frequency;
public:
    DynamicVacuumTerm(double amp = 1e-10, double freq = 1e-15)
        : amplitude(amp), frequency(freq) {}
    
    double compute(double t) const override {
        double rho_vac = 7.09e-36;  // J/m^3
        return amplitude * rho_vac * std::sin(frequency * t);
    }
    
    std::string getName() const override { return "DynamicVacuum"; }
    std::string getDescription() const override { 
        return "Galactic vacuum energy: a_vac=amp·rho_vac·sin(freq·t)";
    }
    std::string getCategory() const override { return "vacuum"; }
};

class QuantumCouplingTerm : public PhysicsTerm {
private:
    double coupling_strength;
    double M, r;  // Spiral galaxy mass and radius
public:
    QuantumCouplingTerm(double strength = 1e-40, double mass = 1.989e41, double radius = 9.258e20)
        : coupling_strength(strength), M(mass), r(radius) {}
    
    double compute(double t) const override {
        double hbar = 1.0546e-34;  // J·s
        return coupling_strength * (hbar * hbar) / (M * r * r) * std::cos(t / 1e6);
    }
    
    std::string getName() const override { return "QuantumCoupling"; }
    std::string getDescription() const override {
        return "Galactic quantum coupling: a_q=(coupling·ℏ²)/(M·r²)·cos(t/1e6)";
    }
    std::string getCategory() const override { return "quantum"; }
};

// Spiral galaxy-specific terms
class SpiralDensityWaveTerm : public PhysicsTerm {
private:
    double Omega_p, r, k_wave, m_arms, G, M;
public:
    SpiralDensityWaveTerm()
        : Omega_p(6.48e-16), r(9.258e20), k_wave(1e-20), m_arms(2.0),  // Omega_p ~ 20 km/s/kpc in rad/s
          G(6.674e-11), M(1.989e41) {}
    
    double compute(double t) const override {
        // Spiral density wave: T_spiral = (G·M·Omega_p²·cos(m·θ))/r² where θ = k·r - Omega_p·t
        // Simplified: T_spiral ~ Omega_p² · cos(k·r - Omega_p·t)
        double theta = k_wave * r - Omega_p * t;
        double T_spiral = (G * M * Omega_p * Omega_p * std::cos(m_arms * theta)) / (r * r);
        return T_spiral;
    }
    
    std::string getName() const override { return "SpiralDensityWave"; }
    std::string getDescription() const override {
        return "Spiral density wave: T_spiral=(G·M·Ω_p²·cos(m·θ))/r² with Ω_p~20 km/s/kpc, m=2 arms";
    }
    std::string getCategory() const override { return "spiral"; }
};

class RotationCurveFlatTerm : public PhysicsTerm {
private:
    double v_rot, r;
public:
    RotationCurveFlatTerm()
        : v_rot(2e5), r(9.258e20) {}  // v_rot ~ 200 km/s (flat rotation curve)
    
    double compute(double t) const override {
        // Flat rotation curve: a_rot = v²/r (centripetal from dark matter halo)
        return (v_rot * v_rot) / r;
    }
    
    std::string getName() const override { return "RotationCurveFlat"; }
    std::string getDescription() const override {
        return "Flat rotation curve: a_rot=v²/r with v~200 km/s (dark matter evidence)";
    }
    std::string getCategory() const override { return "rotation"; }
};

class SupernovaLuminosityTerm : public PhysicsTerm {
private:
    double L_SN, c, r;
public:
    SupernovaLuminosityTerm()
        : L_SN(1e36), c(3e8), r(9.258e20) {}  // L_SN ~ 10³⁶ W (Type Ia/II peak)
    
    double compute(double t) const override {
        // Supernova radiation pressure: P_SN = L_SN/(4πr²c)
        double pi = 3.141592653589793;
        double P_SN = L_SN / (4.0 * pi * r * r * c);
        
        // Acceleration (assuming ISM density ρ~1e-21 kg/m³)
        double rho_ISM = 1e-21;
        return P_SN / rho_ISM;
    }
    
    std::string getName() const override { return "SupernovaLuminosity"; }
    std::string getDescription() const override {
        return "Supernova radiation: a_SN=P_SN/ρ with P=L_SN/(4πr²c), L~1e36 W (Type Ia/II)";
    }
    std::string getCategory() const override { return "supernova"; }
};

class SupernovaShockwaveTerm : public PhysicsTerm {
private:
    double E_SN, M_ejecta, v_shock;
public:
    SupernovaShockwaveTerm()
        : E_SN(1e44), M_ejecta(1.4 * 1.989e30), v_shock(1e7) {}  // E_SN ~ 10^51 erg = 1e44 J, v~1e4 km/s
    
    double compute(double t) const override {
        // Shockwave kinetic energy: E_k = (1/2)·M_ejecta·v²
        // Acceleration from shock: a_shock ~ E_SN / (M_ejecta·r) (energy deposition rate)
        double r = 9.258e20;
        return E_SN / (M_ejecta * r);
    }
    
    std::string getName() const override { return "SupernovaShockwave"; }
    std::string getDescription() const override {
        return "SN shockwave: a_shock=E_SN/(M_ejecta·r) with E~1e44 J (10^51 erg), v~1e4 km/s";
    }
    std::string getCategory() const override { return "supernova"; }
};

class DarkMatterHaloGalacticTerm : public PhysicsTerm {
private:
    double f_DM, G, M, r;
public:
    DarkMatterHaloGalacticTerm()
        : f_DM(0.85), G(6.674e-11), M(1.989e41), r(9.258e20) {}  // 85% DM (typical spiral)
    
    double compute(double t) const override {
        // a_DM = (f_DM·G·M)/r² (dark matter halo)
        return (f_DM * G * M) / (r * r);
    }
    
    std::string getName() const override { return "DarkMatterHaloGalactic"; }
    std::string getDescription() const override {
        return "DM halo (85%): a_DM=(f_DM·G·M)/r² enabling flat rotation curve";
    }
    std::string getCategory() const override { return "dark_matter"; }
};

class CosmologicalLambdaRedshiftTerm : public PhysicsTerm {
private:
    double Lambda, c, H0;
public:
    CosmologicalLambdaRedshiftTerm()
        : Lambda(1.11e-52), c(3e8), H0(2.37e-18) {}  // H0 = 73 km/s/Mpc in SI
    
    double compute(double t) const override {
        // Lambda(z) with redshift dependence: a_Λ(z) = (Λ·c²)/3 · (1+z)³
        // For z up to 1.5 (distant SNe): simplified as a_Λ = (Λ·c²)/3
        return (Lambda * c * c) / 3.0;
    }
    
    std::string getName() const override { return "CosmologicalLambdaRedshift"; }
    std::string getDescription() const override {
        return "Lambda(z): a_Λ=(Λ·c²)/3 for SN cosmology (z up to 1.5)";
    }
    std::string getCategory() const override { return "cosmology"; }
};

class ISMFluidDynamicsTerm : public PhysicsTerm {
private:
    double rho_ISM, V_sys, G, M, r;
public:
    ISMFluidDynamicsTerm()
        : rho_ISM(1e-21), V_sys(3.32e63), G(6.674e-11),  // V ~ (4/3)πr³ for spiral disk
          M(1.989e41), r(9.258e20) {}
    
    double compute(double t) const override {
        // ISM fluid: a_fluid = rho_ISM · V_sys · (G·M)/r²
        double g_base = (G * M) / (r * r);
        return rho_ISM * V_sys * g_base;
    }
    
    std::string getName() const override { return "ISMFluidDynamics"; }
    std::string getDescription() const override {
        return "ISM fluid: a_fluid=rho_ISM·V·g_base with rho~1e-21 kg/m³ (interstellar medium)";
    }
    std::string getCategory() const override { return "fluid"; }
};

class SpiralResonanceTerm : public PhysicsTerm {
private:
    double A, k, omega_spiral, x, pi;
public:
    SpiralResonanceTerm()
        : A(1e-10), k(1e-20), omega_spiral(6.48e-16 * 2.0 * 3.141592653589793),  // ω ~ Omega_p
          x(0.0), pi(3.141592653589793) {}
    
    double compute(double t) const override {
        // Standing wave component (spiral resonances)
        double cos_term = 2.0 * A * std::cos(k * x) * std::cos(omega_spiral * t);
        
        // Traveling wave component
        std::complex<double> i_unit(0.0, 1.0);
        std::complex<double> exp_arg = i_unit * (k * x - omega_spiral * t);
        std::complex<double> exp_term = A * std::exp(exp_arg);
        double real_exp = exp_term.real();
        
        // UQFF cosmological constant factor (2π/13.8)
        double exp_factor = (2.0 * pi) / 13.8;
        
        return cos_term + exp_factor * real_exp;
    }
    
    std::string getName() const override { return "SpiralResonance"; }
    std::string getDescription() const override {
        return "Spiral resonance: 2A·cos(kx)·cos(ωt)+(2π/13.8)·A·Re[exp] at pattern speed Ω_p";
    }
    std::string getCategory() const override { return "resonance"; }
};

// ===========================================================================================
// WOLFRAM TERM REGISTRY (Accumulation: 464 + 10 = 474 terms)
// ===========================================================================================

void registerWolframTerms_source45(std::vector<std::unique_ptr<PhysicsTerm>>& terms) {
    // Inherit 464 terms from source44_wolfram.cpp (delegate to previous registration)
    extern void registerWolframTerms_source44(std::vector<std::unique_ptr<PhysicsTerm>>&);
    registerWolframTerms_source44(terms);
    
    // Add 10 new Spiral & Supernovae terms (465-474)
    terms.push_back(std::make_unique<DynamicVacuumTerm>());                            // 465
    terms.push_back(std::make_unique<QuantumCouplingTerm>(1e-40, 1.989e41, 9.258e20)); // 466 (Spiral params)
    terms.push_back(std::make_unique<SpiralDensityWaveTerm>());                        // 467
    terms.push_back(std::make_unique<RotationCurveFlatTerm>());                        // 468
    terms.push_back(std::make_unique<SupernovaLuminosityTerm>());                      // 469
    terms.push_back(std::make_unique<SupernovaShockwaveTerm>());                       // 470
    terms.push_back(std::make_unique<DarkMatterHaloGalacticTerm>());                   // 471
    terms.push_back(std::make_unique<CosmologicalLambdaRedshiftTerm>());               // 472
    terms.push_back(std::make_unique<ISMFluidDynamicsTerm>());                         // 473
    terms.push_back(std::make_unique<SpiralResonanceTerm>());                          // 474
}

// ===========================================================================================
// PHYSICS SUMMARY
// ===========================================================================================
/*
SPIRAL GALAXIES & SUPERNOVAE UQFF MODULE

SYSTEM: Spiral Galaxy (Milky Way-like) + Type Ia/II Supernovae
- Mass: M = 1.989e41 kg (~1e11 M_sun, typical spiral galaxy)
- Radius: r = 9.258e20 m (~30 kpc, 98,000 light-years disk)
- Rotation velocity: v_rot ~ 200 km/s (flat rotation curve)
- Pattern speed: Ω_p ~ 20 km/s/kpc = 6.48e-16 rad/s (spiral arm rotation)
- Number of arms: m = 2 (typical grand design spiral)
- Hubble constant: H0 = 73 km/s/Mpc = 2.37e-18 s^-1
- Supernova luminosity: L_SN ~ 1e36 W (10⁴³ erg/s at peak)
- SN energy: E_SN ~ 1e44 J (10⁵¹ erg, canonical explosion energy)
- Redshift range: z = 0 to 1.5 (local to cosmological SNe)

DOMINANT PHYSICS:
1. Spiral density wave: T_spiral=(G·M·Ω_p²·cos(m·θ))/r² with m=2 arms, Ω_p~20 km/s/kpc
2. Flat rotation curve: a_rot=v²/r with v~200 km/s (dark matter signature)
3. Supernova luminosity: a_SN=P_SN/ρ with P=L_SN/(4πr²c), L~1e36 W (Type Ia/II)
4. SN shockwave: a_shock=E_SN/(M_ejecta·r) with E~1e44 J, v~1e4 km/s
5. Dark matter halo (85%): a_DM=(f_DM·G·M)/r² enables flat rotation
6. Cosmological Lambda(z): a_Λ=(Λ·c²)/3 for SN cosmology (z up to 1.5)
7. ISM fluid dynamics: a_fluid=rho_ISM·V·g_base with rho~1e-21 kg/m³
8. Spiral resonance: Standing/traveling waves at pattern speed Ω_p

KEY FEATURES:
- Flat rotation curve: Smoking gun for dark matter (~85% of mass)
- Spiral arms: Density waves propagating at Ω_p (NOT material arms)
- Supernovae rate: ~1-3 per century for Milky Way-like galaxy
- Type Ia SNe: Standard candles for cosmology (discovered dark energy acceleration)
- Type II SNe: Core-collapse of massive stars (>8 M_sun)

OBSERVATIONAL DATA:
- Rotation curve measured via H I 21 cm line (neutral hydrogen)
- Spiral arms traced by H II regions, OB associations, molecular clouds
- Type Ia SNe: Perlmutter, Schmidt, Riess (1998, Nobel Prize 2011)
- Hubble diagram: z vs distance modulus (redshift-distance relation)

UQFF PARADIGM:
- No SM gravity - replaced by spiral wave + dark matter halo dynamics
- Flat rotation curve emergent from UQFF dark matter term (85% f_DM)
- Spiral density waves = UQFF resonant structures (NOT material arms)
- Supernovae = Local energy injection modifying UQFF fluid terms
- Cosmological Lambda(z) = Aether expansion (dark energy) from UQFF

10 PhysicsTerm classes (framework + spiral + SN + DM + cosmology)
Total accumulated: 474 classes
*/
