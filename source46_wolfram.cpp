// source46_wolfram.cpp - NGC 6302 Butterfly Nebula UQFF Wolfram Companion
// Wolfram Language integration for NGC6302UQFFModule (Planetary Nebula Evolution)
// Planetary nebula: M=3.98e30 kg (~2 M_sun central star + ejecta), r=9.46e15 m (~0.3 light-years)
// Stellar wind: v_wind=1e5 m/s (~100 km/s), t_eject=2000 yr (ejection time)
// z=0.00095 (distance ~3,800 light-years in Scorpius), rho~1e-20 kg/m³ (nebular density)
// Physics: Stellar wind W_shock, ionized ejecta fluid, resonant oscillatory, DM fraction ~0.85
// UQFF replaces SM gravity with stellar wind-driven shock dynamics
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

// Framework terms (inherited from source45)
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
        return "Planetary nebula vacuum energy: a_vac=amp·rho_vac·sin(freq·t)";
    }
    std::string getCategory() const override { return "vacuum"; }
};

class QuantumCouplingTerm : public PhysicsTerm {
private:
    double coupling_strength;
    double M, r;  // NGC 6302 mass and radius
public:
    QuantumCouplingTerm(double strength = 1e-40, double mass = 3.98e30, double radius = 9.46e15)
        : coupling_strength(strength), M(mass), r(radius) {}
    
    double compute(double t) const override {
        double hbar = 1.0546e-34;  // J·s
        return coupling_strength * (hbar * hbar) / (M * r * r) * std::cos(t / 1e6);
    }
    
    std::string getName() const override { return "QuantumCoupling"; }
    std::string getDescription() const override {
        return "PN quantum coupling: a_q=(coupling·ℏ²)/(M·r²)·cos(t/1e6)";
    }
    std::string getCategory() const override { return "quantum"; }
};

// NGC 6302-specific terms (planetary nebula)
class StellarWindShockTerm : public PhysicsTerm {
private:
    double rho, v_wind, t_eject;
    double seconds_per_year;
public:
    StellarWindShockTerm()
        : rho(1e-20), v_wind(1e5), t_eject(2000.0),
          seconds_per_year(3.154e7) {}  // v_wind ~ 100 km/s
    
    double compute(double t) const override {
        // W_shock = rho · v_wind² · (1 + t/t_eject) (time-dependent shock acceleration)
        double t_years = t / seconds_per_year;
        double time_factor = 1.0 + (t_years / t_eject);
        return rho * v_wind * v_wind * time_factor;
    }
    
    std::string getName() const override { return "StellarWindShock"; }
    std::string getDescription() const override {
        return "Stellar wind shock: W_shock=ρ·v_wind²·(1+t/t_eject) with v~100 km/s, t_eject=2000 yr";
    }
    std::string getCategory() const override { return "stellar_wind"; }
};

class BipolarOutflowTerm : public PhysicsTerm {
private:
    double M_dot, v_wind, r;
public:
    BipolarOutflowTerm()
        : M_dot(1e-6 * 1.989e30 / 3.154e7), v_wind(1e5), r(9.46e15) {}  // Ṁ ~ 1e-6 M_sun/yr
    
    double compute(double t) const override {
        // Bipolar outflow momentum: a_outflow = (Ṁ · v_wind) / (M · r)
        double M = 3.98e30;
        return (M_dot * v_wind) / (M * r);
    }
    
    std::string getName() const override { return "BipolarOutflow"; }
    std::string getDescription() const override {
        return "Bipolar outflow: a_outflow=(Ṁ·v_wind)/(M·r) with Ṁ~1e-6 M_sun/yr (hourglass shape)";
    }
    std::string getCategory() const override { return "stellar_wind"; }
};

class IonizedEjectaFluidTerm : public PhysicsTerm {
private:
    double rho_ejecta, V_sys, G, M, r;
public:
    IonizedEjectaFluidTerm()
        : rho_ejecta(1e-20), V_sys(3.54e48), G(6.674e-11),  // V ~ (4/3)πr³
          M(3.98e30), r(9.46e15) {}
    
    double compute(double t) const override {
        // Ionized ejecta fluid: a_fluid = rho_ejecta · V_sys · (G·M)/r²
        double g_base = (G * M) / (r * r);
        return rho_ejecta * V_sys * g_base;
    }
    
    std::string getName() const override { return "IonizedEjectaFluid"; }
    std::string getDescription() const override {
        return "Ionized ejecta fluid: a_fluid=rho_ejecta·V·g_base with rho~1e-20 kg/m³";
    }
    std::string getCategory() const override { return "fluid"; }
};

class CentralStarRadiationTerm : public PhysicsTerm {
private:
    double L_star, c, r;
public:
    CentralStarRadiationTerm()
        : L_star(1e4 * 3.828e26), c(3e8), r(9.46e15) {}  // L ~ 10,000 L_sun (hot white dwarf)
    
    double compute(double t) const override {
        // Radiation pressure from central star: P_rad = L/(4πr²c)
        double pi = 3.141592653589793;
        double P_rad = L_star / (4.0 * pi * r * r * c);
        
        // Acceleration (assuming ejecta density ρ~1e-20 kg/m³)
        double rho_ejecta = 1e-20;
        return P_rad / rho_ejecta;
    }
    
    std::string getName() const override { return "CentralStarRadiation"; }
    std::string getDescription() const override {
        return "Central star radiation: a_rad=P_rad/ρ with P=L/(4πr²c), L~10,000 L_sun (hot WD)";
    }
    std::string getCategory() const override { return "radiation"; }
};

class DarkMatterPlanetaryNebulaTerm : public PhysicsTerm {
private:
    double f_DM, G, M, r;
public:
    DarkMatterPlanetaryNebulaTerm()
        : f_DM(0.85), G(6.674e-11), M(3.98e30), r(9.46e15) {}  // 85% DM fraction
    
    double compute(double t) const override {
        // a_DM = (f_DM·G·M)/r² (dark matter contribution)
        return (f_DM * G * M) / (r * r);
    }
    
    std::string getName() const override { return "DarkMatterPlanetaryNebula"; }
    std::string getDescription() const override {
        return "DM (~85%): a_DM=(f_DM·G·M)/r² with f_DM=0.85";
    }
    std::string getCategory() const override { return "dark_matter"; }
};

class PlanetaryNebulaResonanceTerm : public PhysicsTerm {
private:
    double A, k, omega_PN, x, pi;
public:
    PlanetaryNebulaResonanceTerm()
        : A(1e-10), k(1e-15), omega_PN(1e-8 * 2.0 * 3.141592653589793),  // k~1/r_PN, low freq
          x(0.0), pi(3.141592653589793) {}
    
    double compute(double t) const override {
        // Standing wave component (PN shell oscillations)
        double cos_term = 2.0 * A * std::cos(k * x) * std::cos(omega_PN * t);
        
        // Traveling wave component
        std::complex<double> i_unit(0.0, 1.0);
        std::complex<double> exp_arg = i_unit * (k * x - omega_PN * t);
        std::complex<double> exp_term = A * std::exp(exp_arg);
        double real_exp = exp_term.real();
        
        // UQFF cosmological constant factor (2π/13.8)
        double exp_factor = (2.0 * pi) / 13.8;
        
        return cos_term + exp_factor * real_exp;
    }
    
    std::string getName() const override { return "PlanetaryNebulaResonance"; }
    std::string getDescription() const override {
        return "PN shell resonance: 2A·cos(kx)·cos(ωt)+(2π/13.8)·A·Re[exp] at PN freq";
    }
    std::string getCategory() const override { return "resonance"; }
};

class LorentzEjectaTerm : public PhysicsTerm {
private:
    double q, v, B;
public:
    LorentzEjectaTerm()
        : q(1.602e-19), v(1e5), B(1e-7) {}  // v~wind velocity, B~0.1 microGauss (weak PN field)
    
    double compute(double t) const override {
        // a_Lorentz = q·|v×B| (electromagnetic force on ionized ejecta)
        // Assuming perpendicular v and B
        return q * v * B;
    }
    
    std::string getName() const override { return "LorentzEjecta"; }
    std::string getDescription() const override {
        return "Electromagnetic: a_L=q·|v×B| from PN B-field (~0.1 µG) and wind (~100 km/s)";
    }
    std::string getCategory() const override { return "electromagnetic"; }
};

class QuantumIntegralPNTerm : public PhysicsTerm {
private:
    double hbar, c, r;
public:
    QuantumIntegralPNTerm()
        : hbar(1.0546e-34), c(3e8), r(9.46e15) {}
    
    double compute(double t) const override {
        // Quantum integral (Casimir-like): a_q_int = (ℏ·c)/(r³)
        return (hbar * c) / (r * r * r);
    }
    
    std::string getName() const override { return "QuantumIntegralPN"; }
    std::string getDescription() const override {
        return "Quantum vacuum integral: a_q_int=(ℏ·c)/r³ (Casimir at PN scale)";
    }
    std::string getCategory() const override { return "quantum"; }
};

// ===========================================================================================
// WOLFRAM TERM REGISTRY (Accumulation: 474 + 10 = 484 terms)
// ===========================================================================================

void registerWolframTerms_source46(std::vector<std::unique_ptr<PhysicsTerm>>& terms) {
    // Inherit 474 terms from source45_wolfram.cpp (delegate to previous registration)
    extern void registerWolframTerms_source45(std::vector<std::unique_ptr<PhysicsTerm>>&);
    registerWolframTerms_source45(terms);
    
    // Add 10 new NGC 6302 terms (475-484)
    terms.push_back(std::make_unique<DynamicVacuumTerm>());                            // 475
    terms.push_back(std::make_unique<QuantumCouplingTerm>(1e-40, 3.98e30, 9.46e15));  // 476 (NGC 6302 params)
    terms.push_back(std::make_unique<StellarWindShockTerm>());                         // 477
    terms.push_back(std::make_unique<BipolarOutflowTerm>());                           // 478
    terms.push_back(std::make_unique<IonizedEjectaFluidTerm>());                       // 479
    terms.push_back(std::make_unique<CentralStarRadiationTerm>());                     // 480
    terms.push_back(std::make_unique<DarkMatterPlanetaryNebulaTerm>());                // 481
    terms.push_back(std::make_unique<PlanetaryNebulaResonanceTerm>());                 // 482
    terms.push_back(std::make_unique<LorentzEjectaTerm>());                            // 483
    terms.push_back(std::make_unique<QuantumIntegralPNTerm>());                        // 484
}

// ===========================================================================================
// PHYSICS SUMMARY
// ===========================================================================================
/*
NGC 6302 BUTTERFLY NEBULA UQFF MODULE (Planetary Nebula Evolution)

SYSTEM: NGC 6302 (Bug Nebula, Butterfly Nebula)
- Mass: M = 3.98e30 kg (~2 M_sun, central star + ejecta)
- Radius: r = 9.46e15 m (~0.3 light-years, 1 parsec, nebular shell)
- Central star: Hot white dwarf, T_eff ~ 250,000 K (one of hottest known)
- Luminosity: L_star ~ 10,000 L_sun = 3.828e30 W (extreme UV emission)
- Stellar wind: v_wind ~ 100 km/s = 1e5 m/s
- Mass loss rate: Ṁ ~ 1e-6 M_sun/yr (past, now ceased)
- Ejection time: t_eject ~ 2000 years ago (bipolar outflow event)
- Distance: z = 0.00095 (~3,800 light-years in Scorpius)
- Ejecta density: rho ~ 1e-20 kg/m³ (ionized hydrogen)
- Magnetic field: B ~ 0.1 microGauss = 1e-7 T (weak PN field)

DOMINANT PHYSICS:
1. Stellar wind shock: W_shock=ρ·v_wind²·(1+t/t_eject) with v~100 km/s, time-dependent
2. Bipolar outflow: a_outflow=(Ṁ·v_wind)/(M·r) creates hourglass/butterfly morphology
3. Ionized ejecta fluid: a_fluid=rho_ejecta·V·g_base (H, He II ionization)
4. Central star radiation: a_rad=P_rad/ρ with L~10,000 L_sun (hot white dwarf drives ionization)
5. Dark matter (85%): a_DM=(f_DM·G·M)/r² (typical PN DM fraction)
6. PN shell resonance: Standing/traveling waves at PN frequencies
7. Lorentz force: q·|v×B| from weak B-field and fast wind
8. Quantum integral: a_q_int=(ℏ·c)/r³ (Casimir at 0.3 ly scale)

KEY FEATURES:
- **Butterfly/hourglass shape:** Bipolar outflow from rotating central star + equatorial disk
- **Extreme chemistry:** Richest known PN in heavy elements (Fe, Ni, Cr detected)
- **Hot central star:** T ~ 250,000 K (among hottest planetary nebula nuclei)
- **Fast winds:** v_wind up to 600 km/s in polar lobes (fastest PN winds)
- **Dust lanes:** Optically thick equatorial torus blocks central star in visible
- **Age:** Ejected ~2000 years ago (very young planetary nebula)

OBSERVATIONAL DATA:
- Discovered by James Dunlop (1826)
- Catalogued as NGC 6302 by John Dreyer (1888)
- HST images reveal intricate substructure (knots, jets, filaments)
- Spitzer infrared: Cold dust (~100 K) in equatorial disk
- Chemical abundances: [Fe/H] ~ +0.5 (iron-rich, unusual for PN)

UQFF PARADIGM:
- No SM gravity - replaced by stellar wind-driven shock dynamics
- Bipolar morphology emergent from UQFF fluid + wind terms
- Central star radiation drives ionization cascades (H II, He II, O III forbidden lines)
- Dark matter (85%) provides gravitational confinement of ejecta
- Stellar wind shock W_shock time-dependent: Peaks during ejection, decays over millennia

10 PhysicsTerm classes (framework + stellar wind + bipolar + radiation + DM)
Total accumulated: 484 classes
*/
