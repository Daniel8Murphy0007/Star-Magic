// Wolfram-Enhanced Resonance MUGE Physics Terms from source7.cpp
// Generated: November 30, 2025
// Modularization: PHASE 4 - Resonance file (source4/5 pattern)
// Total Classes: 13 Resonance MUGE terms
//
// RESONANCE MUGE PHYSICS (13 terms):
// - ResonanceADPMTerm: Differential Plasmotic Motion resonance (vortex dynamics)
// - ResonanceATHzTerm: Terahertz frequency coupling (expansion velocity modulation)
// - ResonanceAvacDiffTerm: Vacuum energy difference (nebula vs ISM)
// - ResonanceASuperFreqTerm: Superconductor frequency (Cooper pair resonance)
// - ResonanceAAetherResTerm: Aether resonance (dark energy replacement)
// - ResonanceUg4iTerm: Integrated Ug4 term (cosmological decay + feedback)
// - ResonanceAQuantumFreqTerm: Quantum frequency (Planck-scale oscillations)
// - ResonanceAAetherFreqTerm: Aether frequency (vacuum energy squared)
// - ResonanceAFluidFreqTerm: Fluid frequency (Navier-Stokes coupling)
// - ResonanceOscTerm: Oscillatory component (approximated to 0 in current model)
// - ResonanceAExpFreqTerm: Expansion frequency (H(z) × t / 2π)
// - ResonanceFTRZTerm: Traversable wormhole scaling factor
// - ResonanceWormholeTerm: Wormhole metric correction (f(r) = 1 - b/r)
//
// PHYSICS VALIDITY:
// - All terms physically motivated from UQFF (Unified Quantum Field Framework)
// - Resonance MUGE dominates at cosmological scales (vacuum energy, expansion frequencies)
//
// DEPENDENCIES:
// - MUGESystem structure (from source7_wolfram.cpp)
// - ResonanceParams structure (from source7_wolfram.cpp)
// - Constants: G, c, PI, H0
//
// AUTHOR: Daniel T. Murphy
// COPYRIGHT: Analyzed Oct 10, 2025; Modularized Nov 30, 2025

#include <cmath>
#include <string>
#include <stdexcept>
#include <map>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// PHYSICS STRUCTURES (from source7_wolfram.cpp)
// ============================================================================

struct MUGESystem
{
    std::string name;
    double I;              // Current (A)
    double A;              // Area (m²)
    double omega1;         // Angular velocity 1 (rad/s)
    double omega2;         // Angular velocity 2 (rad/s)
    double Vsys;           // System volume (m³)
    double vexp;           // Expansion velocity (m/s)
    double t;              // Time (s)
    double z;              // Redshift
    double ffluid;         // Fluid frequency (Hz)
    double M;              // Mass (kg)
    double r;              // Radius (m)
    double B;              // Magnetic field (T)
    double Bcrit;          // Critical magnetic field (T)
    double rho_fluid;      // Fluid density (kg/m³)
    double g_local;        // Local gravity (m/s²)
    double M_DM;           // Dark matter mass (kg)
    double delta_rho_rho;  // Density contrast (dimensionless)
};

struct ResonanceParams
{
    double fDPM = 1e12;           // DPM frequency (Hz)
    double fTHz = 1e12;           // THz frequency (Hz)
    double Evac_neb = 7.09e-36;   // Vacuum energy nebula (J/m³)
    double Evac_ISM = 7.09e-37;   // Vacuum energy ISM (J/m³)
    double Delta_Evac = 6.381e-36;// Vacuum energy difference (J/m³)
    double Fsuper = 6.287e-19;    // Superconductor force (dimensionless)
    double UA_SCM = 10;           // UA scaling (dimensionless)
    double omega_i = 1e-8;        // Angular frequency (rad/s)
    double k4_res = 1.0;          // k4 resonance coupling (dimensionless)
    double freact = 1e10;         // Reactor frequency (Hz)
    double fquantum = 1.445e-17;  // Quantum frequency (Hz)
    double fAether = 1.576e-35;   // Aether frequency (Hz)
    double fosc = 4.57e14;        // Oscillation frequency (Hz)
    double fTRZ = 0.1;            // TRZ factor (dimensionless)
    double c_res = 3e8;           // Speed of light (m/s)
};

// ============================================================================
// CONSTANTS
// ============================================================================

const double PI = M_PI;
const double H0 = 2.269e-18;     // Hubble constant (s⁻¹) = 70 km/s/Mpc

// ============================================================================
// BASE PHYSICS TERM CLASS
// ============================================================================

class PhysicsTerm
{
protected:
    std::string description;
    std::map<std::string, double> parameters;

public:
    virtual ~PhysicsTerm() = default;
    virtual double compute() const = 0;
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const { return description; }
    virtual bool validate() const { return true; }
};

// ============================================================================
// RESONANCE MUGE TERMS (13 classes)
// ============================================================================

class ResonanceADPMTerm : public PhysicsTerm
{
private:
    MUGESystem system;
    ResonanceParams resonance;

public:
    ResonanceADPMTerm(const MUGESystem& sys, const ResonanceParams& res) 
        : system(sys), resonance(res)
    {
        description = "aDPM: Differential Plasmotic Motion = FDPM * fDPM * Evac_neb * c * Vsys (vortex dynamics)";
    }

    double compute() const override
    {
        // FDPM = I * A * (omega1 - omega2)
        double FDPM = system.I * system.A * (system.omega1 - system.omega2);
        return FDPM * resonance.fDPM * resonance.Evac_neb * resonance.c_res * system.Vsys;
    }

    std::string getName() const override { return "ResonanceADPM"; }
};

class ResonanceATHzTerm : public PhysicsTerm
{
private:
    double aDPM;
    MUGESystem system;
    ResonanceParams resonance;

public:
    ResonanceATHzTerm(double aDPM_val, const MUGESystem& sys, const ResonanceParams& res)
        : aDPM(aDPM_val), system(sys), resonance(res)
    {
        description = "aTHz: Terahertz frequency coupling = aDPM * fTHz * vexp / c";
    }

    double compute() const override
    {
        return aDPM * resonance.fTHz * system.vexp / resonance.c_res;
    }

    std::string getName() const override { return "ResonanceATHz"; }
};

class ResonanceAvacDiffTerm : public PhysicsTerm
{
private:
    double aDPM;
    MUGESystem system;
    ResonanceParams resonance;

public:
    ResonanceAvacDiffTerm(double aDPM_val, const MUGESystem& sys, const ResonanceParams& res)
        : aDPM(aDPM_val), system(sys), resonance(res)
    {
        description = "avac_diff: Vacuum energy difference = aDPM * Delta_Evac * vexp / c";
    }

    double compute() const override
    {
        return aDPM * resonance.Delta_Evac * system.vexp / resonance.c_res;
    }

    std::string getName() const override { return "ResonanceAvacDiff"; }
};

class ResonanceASuperFreqTerm : public PhysicsTerm
{
private:
    double aDPM;
    ResonanceParams resonance;

public:
    ResonanceASuperFreqTerm(double aDPM_val, const ResonanceParams& res)
        : aDPM(aDPM_val), resonance(res)
    {
        description = "asuper_freq: Superconductor frequency = aDPM * Fsuper * UA_SCM * omega_i (Cooper pair resonance)";
    }

    double compute() const override
    {
        return aDPM * resonance.Fsuper * resonance.UA_SCM * resonance.omega_i;
    }

    std::string getName() const override { return "ResonanceASuperFreq"; }
};

class ResonanceAAetherResTerm : public PhysicsTerm
{
private:
    double aDPM;
    ResonanceParams resonance;

public:
    ResonanceAAetherResTerm(double aDPM_val, const ResonanceParams& res)
        : aDPM(aDPM_val), resonance(res)
    {
        description = "aaether_res: Aether resonance = aDPM * k4_res * Evac_neb * freact (dark energy replacement)";
    }

    double compute() const override
    {
        return aDPM * resonance.k4_res * resonance.Evac_neb * resonance.freact;
    }

    std::string getName() const override { return "ResonanceAAetherRes"; }
};

class ResonanceUg4iTerm : public PhysicsTerm
{
private:
    double aDPM;
    MUGESystem system;
    ResonanceParams resonance;

public:
    ResonanceUg4iTerm(double aDPM_val, const MUGESystem& sys, const ResonanceParams& res)
        : aDPM(aDPM_val), system(sys), resonance(res)
    {
        description = "Ug4i: Integrated Ug4 term = k4 * Evac_ISM * omega_i * t (cosmological decay)";
    }

    double compute() const override
    {
        return resonance.k4_res * resonance.Evac_ISM * resonance.omega_i * system.t;
    }

    std::string getName() const override { return "ResonanceUg4i"; }
};

class ResonanceAQuantumFreqTerm : public PhysicsTerm
{
private:
    double aDPM;
    ResonanceParams resonance;

public:
    ResonanceAQuantumFreqTerm(double aDPM_val, const ResonanceParams& res)
        : aDPM(aDPM_val), resonance(res)
    {
        description = "aquantum_freq: Quantum frequency = aDPM * fquantum * Evac_neb² (Planck-scale oscillations)";
    }

    double compute() const override
    {
        return aDPM * resonance.fquantum * resonance.Evac_neb * resonance.Evac_neb;
    }

    std::string getName() const override { return "ResonanceAQuantumFreq"; }
};

class ResonanceAAetherFreqTerm : public PhysicsTerm
{
private:
    double aDPM;
    ResonanceParams resonance;

public:
    ResonanceAAetherFreqTerm(double aDPM_val, const ResonanceParams& res)
        : aDPM(aDPM_val), resonance(res)
    {
        description = "aAether_freq: Aether frequency = aDPM * fAether * Evac_neb² (vacuum energy squared)";
    }

    double compute() const override
    {
        return aDPM * resonance.fAether * resonance.Evac_neb * resonance.Evac_neb;
    }

    std::string getName() const override { return "ResonanceAAetherFreq"; }
};

class ResonanceAFluidFreqTerm : public PhysicsTerm
{
private:
    MUGESystem system;
    ResonanceParams resonance;

public:
    ResonanceAFluidFreqTerm(const MUGESystem& sys, const ResonanceParams& res)
        : system(sys), resonance(res)
    {
        description = "afluid_freq: Fluid frequency = ffluid * Vsys * omega_i (Navier-Stokes coupling)";
    }

    double compute() const override
    {
        return system.ffluid * system.Vsys * resonance.omega_i;
    }

    std::string getName() const override { return "ResonanceAFluidFreq"; }
};

class ResonanceOscTerm : public PhysicsTerm
{
public:
    ResonanceOscTerm()
    {
        description = "Osc_term: Oscillatory component (approximated to 0 in current model)";
    }

    double compute() const override
    {
        return 0.0;
    }

    std::string getName() const override { return "ResonanceOsc"; }
};

class ResonanceAExpFreqTerm : public PhysicsTerm
{
private:
    double aDPM;
    MUGESystem system;
    ResonanceParams resonance;
    double H_z;

public:
    ResonanceAExpFreqTerm(double aDPM_val, const MUGESystem& sys, const ResonanceParams& res, double Hz = H0)
        : aDPM(aDPM_val), system(sys), resonance(res), H_z(Hz)
    {
        description = "aexp_freq: Expansion frequency = aDPM * H(z) * t / (2π) (cosmological time evolution)";
    }

    double compute() const override
    {
        return aDPM * H_z * system.t / (2.0 * PI);
    }

    std::string getName() const override { return "ResonanceAExpFreq"; }
};

class ResonanceFTRZTerm : public PhysicsTerm
{
private:
    ResonanceParams resonance;

public:
    ResonanceFTRZTerm(const ResonanceParams& res) : resonance(res)
    {
        description = "fTRZ: Traversable wormhole scaling factor (dimensionless)";
    }

    double compute() const override
    {
        return resonance.fTRZ;
    }

    std::string getName() const override { return "ResonanceFTRZ"; }
};

class ResonanceWormholeTerm : public PhysicsTerm
{
private:
    double r;
    double b;
    double f_worm;
    double Evac_neb;

public:
    ResonanceWormholeTerm(double radius, double b_param = 1.0, double f_w = 1.0, double Evac = 7.09e-36)
        : r(radius), b(b_param), f_worm(f_w), Evac_neb(Evac)
    {
        description = "a_wormhole: Wormhole metric correction f(r) = 1 - b/r (Morris-Thorne traversable wormhole)";
    }

    double compute() const override
    {
        if (r <= b)
            throw std::runtime_error("ResonanceWormholeTerm: r <= b (inside throat)");
        double f_r = 1.0 - b / r;
        return f_worm * f_r * Evac_neb;
    }

    std::string getName() const override { return "ResonanceWormhole"; }
};
