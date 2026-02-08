// Wolfram-Enhanced MUGE Resonance Component Classes from source5.cpp
// Generated: November 29, 2025 (Phase 3)
// Source: UQFF Module 5 - Resonance MUGE Terms
// Total Classes: 14 (13 resonance components + 1 wrapper)
// Status: PHASE 3 COMPLETE - Resonance MUGE modularization

#include <cmath>
#include <string>
#include <map>
#include <memory>

// Constants
const double PI = 3.141592653589793;
const double c_light = 3.0e8; // Speed of light (m/s)

// ============================================================================
// MUGE RESONANCE COMPONENT CLASSES (13 Classes)
// Modular breakdown of compute_resonance_MUGE (Source5.cpp lines 641-712)
// ============================================================================

// CLASS 31: Source5ResonanceADPMTerm - Base DPM Acceleration
class Source5ResonanceADPMTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Extract parameters with SGR1745 defaults
        double I = params.count("I") ? params.at("I") : 1e45;
        double A = params.count("A") ? params.at("A") : 7e22;
        double omega1 = params.count("omega1") ? params.at("omega1") : 1e-8;
        double omega2 = params.count("omega2") ? params.at("omega2") : 5e-9;
        double fDPM = params.count("fDPM") ? params.at("fDPM") : 1e12;
        double Evac_neb = params.count("Evac_neb") ? params.at("Evac_neb") : 7.09e-36;
        double c_res = params.count("c_res") ? params.at("c_res") : 3e8;
        double Vsys = params.count("Vsys") ? params.at("Vsys") : 1e56;

        // FDPM = I * A * (omega1 - omega2)
        double FDPM = I * A * (omega1 - omega2);

        // aDPM = FDPM * fDPM * Evac_neb * c_res * Vsys
        return FDPM * fDPM * Evac_neb * c_res * Vsys;
    }

    std::string getName() const override { return "Source5_ResonanceADPM"; }

    std::string getDescription() const override {
        return "Base DPM acceleration: aDPM = FDPM*fDPM*Evac_neb*c_res*Vsys where FDPM=I*A*(omega1-omega2)";
    }

    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 32: Source5ResonanceATHzTerm - THz Frequency Contribution
class Source5ResonanceATHzTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        double aDPM = params.count("aDPM") ? params.at("aDPM") : 0.0;
        double fTHz = params.count("fTHz") ? params.at("fTHz") : 1e12;
        double Evac_neb = params.count("Evac_neb") ? params.at("Evac_neb") : 7.09e-36;
        double vexp = params.count("vexp") ? params.at("vexp") : 1e6;
        double Evac_ISM = params.count("Evac_ISM") ? params.at("Evac_ISM") : 7.09e-37;
        double c_res = params.count("c_res") ? params.at("c_res") : 3e8;

        // aTHz = fTHz * Evac_neb * vexp * aDPM / (Evac_ISM * c_res)
        return (Evac_ISM * c_res != 0.0) ? 
               (fTHz * Evac_neb * vexp * aDPM) / (Evac_ISM * c_res) : 0.0;
    }

    std::string getName() const override { return "Source5_ResonanceATHz"; }

    std::string getDescription() const override {
        return "THz frequency: aTHz = fTHz*Evac_neb*vexp*aDPM/(Evac_ISM*c_res)";
    }

    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 33: Source5ResonanceAvacDiffTerm - Vacuum Energy Differential
class Source5ResonanceAvacDiffTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        double aDPM = params.count("aDPM") ? params.at("aDPM") : 0.0;
        double Delta_Evac = params.count("Delta_Evac") ? params.at("Delta_Evac") : 6.381e-36;
        double vexp = params.count("vexp") ? params.at("vexp") : 1e6;
        double Evac_neb = params.count("Evac_neb") ? params.at("Evac_neb") : 7.09e-36;
        double c_res = params.count("c_res") ? params.at("c_res") : 3e8;

        // avac_diff = Delta_Evac * vexp^2 * aDPM / (Evac_neb * c_res^2)
        return (Evac_neb * c_res * c_res != 0.0) ? 
               (Delta_Evac * vexp * vexp * aDPM) / (Evac_neb * c_res * c_res) : 0.0;
    }

    std::string getName() const override { return "Source5_ResonanceAvacDiff"; }

    std::string getDescription() const override {
        return "Vacuum energy differential: avac_diff = Delta_Evac*vexp^2*aDPM/(Evac_neb*c_res^2)";
    }

    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 34: Source5ResonanceASuperFreqTerm - Superconductive Frequency
class Source5ResonanceASuperFreqTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        double aDPM = params.count("aDPM") ? params.at("aDPM") : 0.0;
        double Fsuper = params.count("Fsuper") ? params.at("Fsuper") : 6.287e-19;
        double fTHz = params.count("fTHz") ? params.at("fTHz") : 1e12;
        double Evac_neb = params.count("Evac_neb") ? params.at("Evac_neb") : 7.09e-36;
        double c_res = params.count("c_res") ? params.at("c_res") : 3e8;

        // asuper_freq = Fsuper * fTHz * aDPM / (Evac_neb * c_res)
        return (Evac_neb * c_res != 0.0) ? 
               (Fsuper * fTHz * aDPM) / (Evac_neb * c_res) : 0.0;
    }

    std::string getName() const override { return "Source5_ResonanceASuperFreq"; }

    std::string getDescription() const override {
        return "Superconductive frequency: asuper_freq = Fsuper*fTHz*aDPM/(Evac_neb*c_res)";
    }

    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 35: Source5ResonanceAAetherResTerm - Aether Resonance
class Source5ResonanceAAetherResTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        double aDPM = params.count("aDPM") ? params.at("aDPM") : 0.0;
        double UA_SCM = params.count("UA_SCM") ? params.at("UA_SCM") : 1e-10;
        double omega_i = params.count("omega_i") ? params.at("omega_i") : 1e-6;
        double fTHz = params.count("fTHz") ? params.at("fTHz") : 1e12;
        double fTRZ = params.count("fTRZ") ? params.at("fTRZ") : 1e-15;

        // aaether_res = UA_SCM * omega_i * fTHz * aDPM * (1 + fTRZ)
        return UA_SCM * omega_i * fTHz * aDPM * (1.0 + fTRZ);
    }

    std::string getName() const override { return "Source5_ResonanceAAetherRes"; }

    std::string getDescription() const override {
        return "Aether resonance: aaether_res = UA_SCM*omega_i*fTHz*aDPM*(1+fTRZ)";
    }

    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 36: Source5ResonanceUg4iTerm - Reactor Efficiency Contribution
class Source5ResonanceUg4iTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        double aDPM = params.count("aDPM") ? params.at("aDPM") : 0.0;
        double k4_res = params.count("k4_res") ? params.at("k4_res") : 2.0;
        double freact = params.count("freact") ? params.at("freact") : 1e-3;
        double Evac_neb = params.count("Evac_neb") ? params.at("Evac_neb") : 7.09e-36;
        double c_res = params.count("c_res") ? params.at("c_res") : 3e8;
        double t_sys = params.count("t_age") ? params.at("t_age") : t;

        // Ereact = 1046 * exp(-0.0005*t)
        double Ereact = 1046.0 * std::exp(-0.0005 * t_sys);

        // Ug4i = k4_res * Ereact * freact * aDPM / Evac_neb * c_res
        return (Evac_neb != 0.0) ? 
               (k4_res * Ereact * freact * aDPM * c_res / Evac_neb) : 0.0;
    }

    std::string getName() const override { return "Source5_ResonanceUg4i"; }

    std::string getDescription() const override {
        return "Reactor efficiency: Ug4i = k4_res*Ereact*freact*aDPM*c_res/Evac_neb where Ereact=1046*exp(-0.0005*t)";
    }

    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 37: Source5ResonanceAQuantumFreqTerm - Quantum Frequency
class Source5ResonanceAQuantumFreqTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        double aDPM = params.count("aDPM") ? params.at("aDPM") : 0.0;
        double fquantum = params.count("fquantum") ? params.at("fquantum") : 1e12;
        double Evac_neb = params.count("Evac_neb") ? params.at("Evac_neb") : 7.09e-36;
        double Evac_ISM = params.count("Evac_ISM") ? params.at("Evac_ISM") : 7.09e-37;
        double c_res = params.count("c_res") ? params.at("c_res") : 3e8;

        // aquantum_freq = fquantum * Evac_neb * aDPM / (Evac_ISM * c_res)
        return (Evac_ISM * c_res != 0.0) ? 
               (fquantum * Evac_neb * aDPM) / (Evac_ISM * c_res) : 0.0;
    }

    std::string getName() const override { return "Source5_ResonanceAQuantumFreq"; }

    std::string getDescription() const override {
        return "Quantum frequency: aquantum_freq = fquantum*Evac_neb*aDPM/(Evac_ISM*c_res)";
    }

    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 38: Source5ResonanceAAetherFreqTerm - Aether Frequency
class Source5ResonanceAAetherFreqTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        double aDPM = params.count("aDPM") ? params.at("aDPM") : 0.0;
        double fAether = params.count("fAether") ? params.at("fAether") : 1e12;
        double Evac_neb = params.count("Evac_neb") ? params.at("Evac_neb") : 7.09e-36;
        double Evac_ISM = params.count("Evac_ISM") ? params.at("Evac_ISM") : 7.09e-37;
        double c_res = params.count("c_res") ? params.at("c_res") : 3e8;

        // aAether_freq = fAether * Evac_neb * aDPM / (Evac_ISM * c_res)
        return (Evac_ISM * c_res != 0.0) ? 
               (fAether * Evac_neb * aDPM) / (Evac_ISM * c_res) : 0.0;
    }

    std::string getName() const override { return "Source5_ResonanceAAetherFreq"; }

    std::string getDescription() const override {
        return "Aether frequency: aAether_freq = fAether*Evac_neb*aDPM/(Evac_ISM*c_res)";
    }

    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 39: Source5ResonanceAFluidFreqTerm - Fluid Frequency
class Source5ResonanceAFluidFreqTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        double ffluid = params.count("ffluid") ? params.at("ffluid") : 3.465e-8;
        double Evac_neb = params.count("Evac_neb") ? params.at("Evac_neb") : 7.09e-36;
        double Vsys = params.count("Vsys") ? params.at("Vsys") : 1e56;
        double Evac_ISM = params.count("Evac_ISM") ? params.at("Evac_ISM") : 7.09e-37;
        double c_res = params.count("c_res") ? params.at("c_res") : 3e8;

        // afluid_freq = ffluid * Evac_neb * Vsys / (Evac_ISM * c_res)
        return (Evac_ISM * c_res != 0.0) ? 
               (ffluid * Evac_neb * Vsys) / (Evac_ISM * c_res) : 0.0;
    }

    std::string getName() const override { return "Source5_ResonanceAFluidFreq"; }

    std::string getDescription() const override {
        return "Fluid frequency: afluid_freq = ffluid*Evac_neb*Vsys/(Evac_ISM*c_res)";
    }

    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 40: Source5ResonanceOscTermTerm - Oscillation Term (Placeholder)
class Source5ResonanceOscTermTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        return 0.0;  // Placeholder for future oscillation physics
    }

    std::string getName() const override { return "Source5_ResonanceOscTerm"; }

    std::string getDescription() const override {
        return "Oscillation term: placeholder for future resonant oscillation physics";
    }

    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 41: Source5ResonanceAExpFreqTerm - Expansion Frequency
class Source5ResonanceAExpFreqTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        double aDPM = params.count("aDPM") ? params.at("aDPM") : 0.0;
        double Evac_neb = params.count("Evac_neb") ? params.at("Evac_neb") : 7.09e-36;
        double Evac_ISM = params.count("Evac_ISM") ? params.at("Evac_ISM") : 7.09e-37;
        double c_res = params.count("c_res") ? params.at("c_res") : 3e8;
        double t_sys = params.count("t_age") ? params.at("t_age") : t;
        double z = params.count("z") ? params.at("z") : 0.01;  // Redshift

        // H_z = H0 * sqrt(Omega_m * (1+z)^3 + Omega_Lambda)
        // Simplified: H_z ~ H0 (1+z)^1.5
        double H0 = 2.269e-18;  // Hubble constant
        double H_z = H0 * std::pow(1.0 + z, 1.5);

        // fexp = 2*PI*H_z*t
        double fexp = 2.0 * PI * H_z * t_sys;

        // aexp_freq = fexp * Evac_neb * aDPM / (Evac_ISM * c_res)
        return (Evac_ISM * c_res != 0.0) ? 
               (fexp * Evac_neb * aDPM) / (Evac_ISM * c_res) : 0.0;
    }

    std::string getName() const override { return "Source5_ResonanceAExpFreq"; }

    std::string getDescription() const override {
        return "Expansion frequency: aexp_freq = fexp*Evac_neb*aDPM/(Evac_ISM*c_res) where fexp=2*PI*H_z*t";
    }

    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 42: Source5ResonanceFTRZTerm - TRZ Frequency
class Source5ResonanceFTRZTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        double fTRZ = params.count("fTRZ") ? params.at("fTRZ") : 1e-15;
        return fTRZ;  // Direct contribution
    }

    std::string getName() const override { return "Source5_ResonanceFTRZ"; }

    std::string getDescription() const override {
        return "TRZ frequency: direct contribution fTRZ (transition redshift zone)";
    }

    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// CLASS 43: Source5ResonanceAWormholeTerm - Wormhole Contribution
class Source5ResonanceAWormholeTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        double r = params.count("r") ? params.at("r") : 1e13;
        double b = params.count("b_wormhole") ? params.at("b_wormhole") : 1.0;
        double f_worm = params.count("f_worm") ? params.at("f_worm") : 1.0;
        double Evac_neb = params.count("Evac_neb") ? params.at("Evac_neb") : 7.09e-36;

        // a_wormhole = f_worm * Evac_neb * (1 / (b^2 + r^2))
        return f_worm * Evac_neb * (1.0 / (b * b + r * r));
    }

    std::string getName() const override { return "Source5_ResonanceAWormhole"; }

    std::string getDescription() const override {
        return "Wormhole contribution: a_wormhole = f_worm*Evac_neb/(b^2+r^2)";
    }

    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// ============================================================================
// COMBINED RESONANCE MUGE TERM (Wrapper)
// ============================================================================

class Source5ResonanceMUGETerm : public PhysicsTerm {
private:
    std::unique_ptr<Source5ResonanceADPMTerm> aDPM_term;
    std::unique_ptr<Source5ResonanceATHzTerm> aTHz_term;
    std::unique_ptr<Source5ResonanceAvacDiffTerm> avac_diff_term;
    std::unique_ptr<Source5ResonanceASuperFreqTerm> asuper_freq_term;
    std::unique_ptr<Source5ResonanceAAetherResTerm> aaether_res_term;
    std::unique_ptr<Source5ResonanceUg4iTerm> Ug4i_term;
    std::unique_ptr<Source5ResonanceAQuantumFreqTerm> aquantum_freq_term;
    std::unique_ptr<Source5ResonanceAAetherFreqTerm> aAether_freq_term;
    std::unique_ptr<Source5ResonanceAFluidFreqTerm> afluid_freq_term;
    std::unique_ptr<Source5ResonanceOscTermTerm> Osc_term_term;
    std::unique_ptr<Source5ResonanceAExpFreqTerm> aexp_freq_term;
    std::unique_ptr<Source5ResonanceFTRZTerm> fTRZ_term;
    std::unique_ptr<Source5ResonanceAWormholeTerm> a_wormhole_term;
    
public:
    Source5ResonanceMUGETerm()
        : aDPM_term(std::make_unique<Source5ResonanceADPMTerm>()),
          aTHz_term(std::make_unique<Source5ResonanceATHzTerm>()),
          avac_diff_term(std::make_unique<Source5ResonanceAvacDiffTerm>()),
          asuper_freq_term(std::make_unique<Source5ResonanceASuperFreqTerm>()),
          aaether_res_term(std::make_unique<Source5ResonanceAAetherResTerm>()),
          Ug4i_term(std::make_unique<Source5ResonanceUg4iTerm>()),
          aquantum_freq_term(std::make_unique<Source5ResonanceAQuantumFreqTerm>()),
          aAether_freq_term(std::make_unique<Source5ResonanceAAetherFreqTerm>()),
          afluid_freq_term(std::make_unique<Source5ResonanceAFluidFreqTerm>()),
          Osc_term_term(std::make_unique<Source5ResonanceOscTermTerm>()),
          aexp_freq_term(std::make_unique<Source5ResonanceAExpFreqTerm>()),
          fTRZ_term(std::make_unique<Source5ResonanceFTRZTerm>()),
          a_wormhole_term(std::make_unique<Source5ResonanceAWormholeTerm>()) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Compute aDPM first (base term)
        double aDPM = aDPM_term->compute(t, params);
        
        // Create augmented params with aDPM
        std::map<std::string, double> aug_params = params;
        aug_params["aDPM"] = aDPM;
        
        // Compute all remaining terms
        double aTHz = aTHz_term->compute(t, aug_params);
        double avac_diff = avac_diff_term->compute(t, aug_params);
        double asuper_freq = asuper_freq_term->compute(t, aug_params);
        double aaether_res = aaether_res_term->compute(t, aug_params);
        double Ug4i = Ug4i_term->compute(t, aug_params);
        double aquantum_freq = aquantum_freq_term->compute(t, aug_params);
        double aAether_freq = aAether_freq_term->compute(t, aug_params);
        double afluid_freq = afluid_freq_term->compute(t, aug_params);
        double Osc_term = Osc_term_term->compute(t, aug_params);
        double aexp_freq = aexp_freq_term->compute(t, aug_params);
        double fTRZ = fTRZ_term->compute(t, aug_params);
        double a_wormhole = a_wormhole_term->compute(t, aug_params);
        
        // Sum all components
        return aDPM + aTHz + avac_diff + asuper_freq + aaether_res + Ug4i + 
               aquantum_freq + aAether_freq + afluid_freq + Osc_term + aexp_freq + fTRZ + a_wormhole;
    }
    
    std::string getName() const override { return "Source5_ResonanceMUGE_Full"; }
    
    std::string getDescription() const override {
        return "Full Resonance MUGE: sum of all 13 resonance terms (aDPM + aTHz + avac_diff + ... + a_wormhole)";
    }
    
    bool validate(const std::map<std::string, double>&) const override { return true; }
};

// ============================================================================
// REGISTRATION FUNCTION
// ============================================================================

void registerResonanceTerms_source5(PhysicsTermRegistry& registry) {
    // Individual Component Classes (13)
    registry.registerPhysicsTerm("Source5_ResonanceADPM", std::make_unique<Source5ResonanceADPMTerm>(), "wolfram_resonance");
    registry.registerPhysicsTerm("Source5_ResonanceATHz", std::make_unique<Source5ResonanceATHzTerm>(), "wolfram_resonance");
    registry.registerPhysicsTerm("Source5_ResonanceAvacDiff", std::make_unique<Source5ResonanceAvacDiffTerm>(), "wolfram_resonance");
    registry.registerPhysicsTerm("Source5_ResonanceASuperFreq", std::make_unique<Source5ResonanceASuperFreqTerm>(), "wolfram_resonance");
    registry.registerPhysicsTerm("Source5_ResonanceAAetherRes", std::make_unique<Source5ResonanceAAetherResTerm>(), "wolfram_resonance");
    registry.registerPhysicsTerm("Source5_ResonanceUg4i", std::make_unique<Source5ResonanceUg4iTerm>(), "wolfram_resonance");
    registry.registerPhysicsTerm("Source5_ResonanceAQuantumFreq", std::make_unique<Source5ResonanceAQuantumFreqTerm>(), "wolfram_resonance");
    registry.registerPhysicsTerm("Source5_ResonanceAAetherFreq", std::make_unique<Source5ResonanceAAetherFreqTerm>(), "wolfram_resonance");
    registry.registerPhysicsTerm("Source5_ResonanceAFluidFreq", std::make_unique<Source5ResonanceAFluidFreqTerm>(), "wolfram_resonance");
    registry.registerPhysicsTerm("Source5_ResonanceOscTerm", std::make_unique<Source5ResonanceOscTermTerm>(), "wolfram_resonance");
    registry.registerPhysicsTerm("Source5_ResonanceAExpFreq", std::make_unique<Source5ResonanceAExpFreqTerm>(), "wolfram_resonance");
    registry.registerPhysicsTerm("Source5_ResonanceFTRZ", std::make_unique<Source5ResonanceFTRZTerm>(), "wolfram_resonance");
    registry.registerPhysicsTerm("Source5_ResonanceAWormhole", std::make_unique<Source5ResonanceAWormholeTerm>(), "wolfram_resonance");
    
    // Combined Wrapper (1)
    registry.registerPhysicsTerm("Source5_ResonanceMUGE_Full", std::make_unique<Source5ResonanceMUGETerm>(), "wolfram_resonance");
}

// ============================================================================
// TOTAL: 14 CLASSES (13 components + 1 wrapper)
// - aDPM: Base DPM acceleration
// - aTHz: THz frequency contribution
// - avac_diff: Vacuum energy differential
// - asuper_freq: Superconductive frequency
// - aaether_res: Aether resonance
// - Ug4i: Reactor efficiency contribution
// - aquantum_freq: Quantum frequency
// - aAether_freq: Aether frequency
// - afluid_freq: Fluid frequency
// - Osc_term: Oscillation placeholder
// - aexp_freq: Expansion frequency
// - fTRZ: TRZ frequency
// - a_wormhole: Wormhole contribution
// - ResonanceMUGE_Full: Complete equation wrapper
// ============================================================================
