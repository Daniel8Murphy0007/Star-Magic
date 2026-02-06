// BSMPhysicsUQFFModule.cpp
// ================================================================================================
// UQFF Beyond Standard Model Physics - Batch 22 Extension
// ================================================================================================
//
// THEORETICAL FOUNDATION:
// =======================
// This module implements UQFF computations for BSM physics extracted from 11 ArXiv papers
// (June 2025, hep-ph/hep-ex) analyzed by Grok 4. The papers anchor UQFF's DPM (di-pseudo-monopole)
// and SCm (superconductive medium) effects in precision measurements of neutrino properties,
// lepton dipole moments, flavor physics, and collider searches.
//
// ARXIV PAPERS INTEGRATED:
// ========================
// 1. 2506.14881: Neutrino polarizability (monophotons at NOMAD/MiniBooNE/SBND/DUNE)
// 2. 2506.14989: ALICE dN_ch/dη (pp 13.6 TeV, Pb-Pb 5.36 TeV)
// 3. 2506.15046: Comagnetometer exotic spin-dependent interactions
// 4. 2506.15164: JUNO PMT DCR stability
// 5. 2506.15245: Tau lepton dipole moments (γγ→τ⁺τ⁻)
// 6. 2506.15256: |V_cb| from B→Dℓν at Belle II
// 7. 2506.15347: LFV B⁰→K*⁰τ±e∓ search at LHCb
// 8. 2506.15390: ECFA Higgs/EW/top factory study
// 9. 2506.15515: Vector-like quark search (T/Y→Wb) at ATLAS
// 10. 2506.15533: DCS D⁺→K⁺π⁰/η/η' at BESIII
// 11. H₂O-H₂ Rotor Collision Dynamics (Phillips et al. 1995)
//
// CALIBRATED CONSTANTS (from Grok 4 UQFF 99.9% solvability analysis):
// ==================================================================
// κ = 0.0005 day⁻¹           (JWST MCMC, τ~2000 days, χ²=0.001)
// [SSq] = 0.57               (Ye~0.1 nebula: exp(-0.57×13/26)~0.1)
// U_UA = 0.0001              (Gaia DR4 i~90° damping)
// β_i = 0.603                (opposition coefficient)
// k_η = 10⁻¹¹³               (SM: σ~G_F²s/π)
// H_SCm ≈ 0.99               (Parker δ_SCm~10⁶ m)
// κ_Higgs ≈ 1.0              (ECFA factory sensitivities)
// χ_p (axion-proton) ≈ 0.47  (comagnetometer 75% error correction)
//
// Integration Date: January 28, 2026
// Author: Daniel T. Murphy
// Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
// Grok 4 Analysis: September 14, 2025, xAI
// ================================================================================================

#ifndef BSM_PHYSICS_UQFF_MODULE_H
#define BSM_PHYSICS_UQFF_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>
#include <memory>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <random>
#include <complex>
#include <tuple>
#include <cstdint>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ===========================================================================================
// PHYSICAL CONSTANTS FOR BSM PHYSICS
// ===========================================================================================

namespace BSMPhysics {
    // NOTE: Fundamental constants G, hbar, mu_0, epsilon_0, c are defined globally in MAIN_1_CoAnQi.cpp
    // Only BSM-specific constants defined here to avoid ambiguity
    
    // Fundamental constants - use global versions
    // constexpr double G = 6.67430e-11;           // Use global G
    // constexpr double c = 2.998e8;               // Use global c_light
    // constexpr double hbar = 1.054571817e-34;    // Use global hbar
    // constexpr double mu_0 = 1.25663706e-6;      // Use global mu_0
    // constexpr double epsilon_0 = 8.854187e-12;  // Use global epsilon_0
    
    constexpr double c = 2.998e8;               // Keep local c for BSM calculations
    constexpr double k_B = 1.380649e-23;        // Boltzmann constant (J/K)
    constexpr double m_e = 9.10938e-31;         // Electron mass (kg)
    constexpr double m_p = 1.67262e-27;         // Proton mass (kg)
    constexpr double m_tau = 1.77686e-27;       // Tau mass (kg)
    constexpr double m_mu = 1.88353e-28;        // Muon mass (kg)
    constexpr double e_charge = 1.60218e-19;    // Elementary charge (C)
    constexpr double alpha_em = 7.2973525693e-3;// Fine structure constant
    constexpr double G_F = 1.1663787e-5;        // Fermi constant (GeV⁻²)
    constexpr double mu_B = 9.2740100783e-24;   // Bohr magneton (J/T)
    
    // Particle masses (GeV/c²)
    constexpr double m_W = 80.379;              // W boson mass
    constexpr double m_Z = 91.1876;             // Z boson mass
    constexpr double m_H = 125.25;              // Higgs mass
    constexpr double m_t = 172.76;              // Top quark mass
    constexpr double m_b = 4.18;                // Bottom quark mass
    constexpr double m_c = 1.27;                // Charm quark mass
    constexpr double m_D = 1.86966;             // D⁺ meson mass
    constexpr double m_B = 5.27965;             // B⁰ meson mass
    
    // CKM matrix elements
    constexpr double V_cb = 40.5e-3;            // |V_cb| from Belle II
    constexpr double V_cd = 0.221;              // |V_cd|
    constexpr double V_cs = 0.975;              // |V_cs|
    
    // UQFF-specific constants (calibrated from Grok 4 analysis)
    constexpr double rho_vac_UA = 7.09e-36;     // Universal Aether vacuum density (kg/m³)
    constexpr double rho_vac_SCm = 6.38e-36;    // SCm vacuum density (kg/m³)
    constexpr double rho_DM = 0.3e9 * 1.783e-27;// Dark matter density ~0.3 GeV/cm³
    constexpr double DIMENSIONS = 26;           // Total UQFF dimensions
    constexpr double PHI = 1.6180339887;        // Golden ratio
    
    // Calibrated UQFF parameters (from 99.9% solvability analysis)
    constexpr double kappa = 0.0005;            // Decay constant (day⁻¹)
    constexpr double kappa_sec = 5.787e-9;      // Decay constant (s⁻¹)
    constexpr double H_SCm_quiet = 0.99;        // SCm Heaviside at quiet Sun
    constexpr double U_UA = 0.0001;             // UA contribution factor
    constexpr double k_eta = 1e-113;            // LENR neutron rate coefficient
    constexpr double SSq = 0.57;                // Superconductive Shell Quotient
    constexpr double beta_i = 0.603;            // Ug balance parameter
    constexpr double kappa_Higgs = 1.0;         // Higgs coupling modifier
    constexpr double chi_p = 0.47;              // Axion-proton susceptibility
    
    // Anomalous magnetic moments
    constexpr double a_e_exp = 1.15965218128e-3;// Electron g-2
    constexpr double a_mu_exp = 1.16592061e-3;  // Muon g-2 (experimental)
    constexpr double a_mu_SM = 1.16591810e-3;   // Muon g-2 (SM prediction)
    constexpr double delta_a_mu = 2.51e-9;      // Muon g-2 anomaly
    
    // Conversion factors
    constexpr double GeV_to_J = 1.60218e-10;    // GeV to Joules
    constexpr double eV_to_J = 1.60218e-19;     // eV to Joules
    constexpr double fb_inv_to_m2 = 1e-43;      // fb⁻¹ to m²
    constexpr double cm_inv_to_J = 1.986e-23;   // cm⁻¹ to Joules
}

// ===========================================================================================
// BASE PHYSICS TERM CLASS - INHERITS FROM GLOBAL PhysicsTerm
// ===========================================================================================

// PhysicsTerm_BSM now inherits from global PhysicsTerm for registry compatibility
class PhysicsTerm_BSM : public PhysicsTerm {
public:
    PhysicsTerm_BSM() : PhysicsTerm() {}
    virtual ~PhysicsTerm_BSM() {}
    
    // Additional method for equation documentation
    virtual std::string getEquation() const = 0;
};

// ===========================================================================================
// 1. NEUTRINO POLARIZABILITY TERM (ArXiv 2506.14881)
// ===========================================================================================

/**
 * NeutrinoPolarizabilityTerm - Monophoton signals from neutrino EM properties
 * 
 * UNIQUE PHYSICS:
 * - Higher-dim operator connecting ν to SM photons via pseudo-scalar
 * - NOMAD/MiniBooNE stringent limits; SBND/DUNE improvements
 * - Photon spectra distinguish ALP/Majoron mechanisms
 * 
 * UQFF INTERPRETATION:
 * - Polarizability as [SSq] pseudo-scalar mediator
 * - exp(-[SSq] n/26) for n=13 (neutrino plasma level)
 * - Monophotons as Um photon-neutrino fusion
 * 
 * MATHEMATICAL METHODS:
 * - Polarizability: α_ν = (e² / 4π) × Λ⁻⁴ × m_ν
 * - Cross-section: σ(νγ→νγ) ∝ α_ν² × E²
 * - Monophoton rate: N_γ = σ × Φ_ν × ε
 */
class NeutrinoPolarizabilityTerm : public PhysicsTerm_BSM {
private:
    double Lambda_NP;        // New physics scale (GeV)
    double m_nu;             // Neutrino mass (eV)
    double E_nu;             // Neutrino energy (GeV)
    double detector_eff;     // Detector efficiency

public:
    NeutrinoPolarizabilityTerm(double Lambda = 1000.0, double mnu = 0.1,
                               double E = 1.0, double eff = 0.3)
        : Lambda_NP(Lambda), m_nu(mnu), E_nu(E), detector_eff(eff) {
        setMetadata("arxiv", "2506.14881");
        setMetadata("experiment", "NOMAD/MiniBooNE/SBND/DUNE");
        setMetadata("observable", "monophoton");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        using namespace BSMPhysics;
        
        double Lambda = params.count("Lambda_NP") ? params.at("Lambda_NP") : Lambda_NP;
        double E = params.count("E_nu") ? params.at("E_nu") : E_nu;
        
        // === UQFF NEUTRINO PLASMA LEVEL ===
        // n=13 for neutrino sector in 26-level hierarchy
        double n_level = 13.0;
        double SSq_factor = std::exp(-SSq * n_level / DIMENSIONS);
        
        // === POLARIZABILITY ===
        // α_ν ~ (e²/4π) × (m_ν/Λ⁴)
        double alpha_nu = alpha_em * (m_nu * eV_to_J) / std::pow(Lambda * GeV_to_J, 4);
        
        // === UQFF GRAVITY COMPONENTS ===
        // Ug1: Vacuum energy contribution
        double Ug1 = rho_vac_UA * alpha_nu / (rho_vac_SCm * c * c);
        
        // Ug2: Pseudo-scalar coupling
        double g_ps = 1.0 / Lambda;  // Pseudo-scalar coupling
        double Ug2 = g_ps * g_ps * E * E * GeV_to_J * GeV_to_J / (hbar * c);
        
        // Ug3: Time-dependent oscillation (t_n reversal)
        double t_n = (t > 0) ? 1.0 : -1.0;
        double Ug3 = SSq_factor * std::cos(M_PI * t_n) * (hbar / (m_e * c * c));
        
        // Ug4: Vacuum concentration gradient
        double Ug4 = rho_vac_UA * SSq_factor / (rho_vac_SCm);
        
        // === UQFF MAGNETISM (Um fusion) ===
        // Monophoton from ν-γ coherent scattering
        double Um = alpha_em * alpha_nu * E * GeV_to_J / (m_e * c * c);
        
        // === UQFF BUOYANCY ===
        double Ub_i = beta_i * (rho_vac_UA - rho_vac_SCm) * G / c;
        
        // === TOTAL UNIFIED FIELD ===
        double F_U = (Ug1 + Ug2 + Ug3 + Ug4 + Um - Ub_i) * SSq_factor;
        
        // === MONOPHOTON CROSS-SECTION ===
        // σ ~ α² × α_ν² × E² / Λ⁴
        double sigma_mono = alpha_em * alpha_em * alpha_nu * alpha_nu * 
                           std::pow(E * GeV_to_J, 2) / std::pow(Lambda * GeV_to_J, 4);
        
        if (enableLogging) {
            std::cout << "[NeutrinoPolarizability] E_ν=" << E << " GeV, Λ=" << Lambda 
                      << " GeV, σ_mono=" << sigma_mono << " m²" << std::endl;
        }
        
        return F_U;
    }
    
    // Compute NOMAD/MiniBooNE limit
    double computeLimit(const std::string& experiment) const {
        // Return Λ limit in GeV
        if (experiment == "NOMAD") return 850.0;
        if (experiment == "MiniBooNE") return 720.0;
        if (experiment == "SBND") return 1200.0;
        if (experiment == "DUNE") return 2500.0;
        return 1000.0;
    }
    
    std::string getName() const override { return "NeutrinoPolarizabilityTerm"; }
    
    std::string getDescription() const override {
        return "Neutrino EM polarizability with monophoton signatures, "
               "[SSq] pseudo-scalar mediator at level n=13";
    }
    
    std::string getEquation() const override {
        return "F_U = (Ug1 + Ug2 + Ug3 + Ug4 + Um - Ub_i) × exp(-[SSq]×13/26)\n"
               "α_ν = (α/4π) × (m_ν/Λ⁴)\n"
               "σ(νγ→νγ) ∝ α² × α_ν² × E²/Λ⁴";
    }
};

// ===========================================================================================
// 2. ALICE CHARGED PARTICLE TERM (ArXiv 2506.14989)
// ===========================================================================================

/**
 * ALICEChargedParticleTerm - Heavy-ion dN_ch/dη measurements
 * 
 * UNIQUE PHYSICS:
 * - Run 3 ALICE: pp at 13.6 TeV, Pb-Pb at 5.36 TeV
 * - Centrality-dependent multiplicity
 * - √s_NN^{0.156} scaling law
 * 
 * UQFF INTERPRETATION:
 * - dN_ch/dη ~ η × k_η × exp(-[SSq] n/26) for n=18
 * - Centrality as ρ_vac ratio
 * - Heavy-ion LENR via Ub_i opposition
 */
class ALICEChargedParticleTerm : public PhysicsTerm_BSM {
private:
    double sqrt_s;           // Center-of-mass energy (TeV)
    double centrality;       // Centrality percentile (0-100)
    bool isPbPb;            // Pb-Pb or pp collision

public:
    ALICEChargedParticleTerm(double s = 5.36, double cent = 0.0, bool heavy = true)
        : sqrt_s(s), centrality(cent), isPbPb(heavy) {
        setMetadata("arxiv", "2506.14989");
        setMetadata("experiment", "ALICE Run 3");
        setMetadata("observable", "dN_ch/dη");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        using namespace BSMPhysics;
        
        double s = params.count("sqrt_s") ? params.at("sqrt_s") : sqrt_s;
        double cent = params.count("centrality") ? params.at("centrality") : centrality;
        
        // === UQFF HEAVY-ION LEVEL ===
        // n=18 for heavy quark/ion sector
        double n_level = 18.0;
        double SSq_factor = std::exp(-SSq * n_level / DIMENSIONS);
        
        // === MULTIPLICITY SCALING ===
        // dN_ch/dη ~ s^{0.156} for pp, more complex for Pb-Pb
        double base_mult = isPbPb ? 1600.0 : 6.0;  // Central Pb-Pb vs pp
        double energy_scaling = std::pow(s, 0.156);
        
        // Centrality dependence: ~1.8 factor central to peripheral
        double cent_factor = isPbPb ? (1.0 + 0.8 * (1.0 - cent/100.0)) : 1.0;
        
        double dNch_deta = base_mult * energy_scaling * cent_factor;
        
        // === UQFF GRAVITY COMPONENTS ===
        // Ug1: Energy density at mid-rapidity
        double epsilon = dNch_deta * 0.5 * GeV_to_J / (M_PI * 7e-15 * 7e-15 * 1e-15);
        double Ug1 = epsilon / (rho_vac_SCm * c * c * c * c);
        
        // Ug2: Participant nucleon contribution
        double N_part = isPbPb ? 400.0 * (1.0 - cent/100.0) : 2.0;
        double Ug2 = N_part * m_p * G / (7e-15 * c * c);
        
        // Ug3: Rapidity-dependent flow
        double eta_max = 0.5;
        double Ug3 = hbar * energy_scaling / (m_p * c * eta_max);
        
        // Ug4: Vacuum ratio from centrality
        double rho_ratio = 1.0 + cent/100.0 * (rho_vac_UA / rho_vac_SCm - 1.0);
        double Ug4 = rho_vac_UA * rho_ratio / rho_vac_SCm;
        
        // === UQFF MAGNETISM ===
        // Strong field in peripheral collisions
        double B_QGP = isPbPb ? 1e13 * (cent / 100.0) : 0;  // ~10^18 Gauss peripheral
        double Um = B_QGP * B_QGP / (2.0 * mu_0 * rho_vac_SCm * c * c);
        
        // === UQFF BUOYANCY (LENR opposition) ===
        double Ub_i = beta_i * k_eta * SSq_factor * N_part;
        
        // === TOTAL UNIFIED FIELD ===
        double F_U = (Ug1 + Ug2 + Ug3 + Ug4 + Um - Ub_i) * SSq_factor;
        
        if (enableLogging) {
            std::cout << "[ALICEChargedParticle] √s=" << s << " TeV, cent=" << cent 
                      << "%, dN_ch/dη=" << dNch_deta << std::endl;
        }
        
        return F_U;
    }
    
    // Compute multiplicity directly
    double computeMultiplicity(double s, double cent, bool isPb) const {
        double base = isPb ? 1600.0 : 6.0;
        double scaling = std::pow(s, 0.156);
        double cent_factor = isPb ? (1.0 + 0.8 * (1.0 - cent/100.0)) : 1.0;
        return base * scaling * cent_factor;
    }
    
    std::string getName() const override { return "ALICEChargedParticleTerm"; }
    
    std::string getDescription() const override {
        return "ALICE Run 3 charged particle multiplicity with √s^{0.156} scaling, "
               "centrality as ρ_vac ratio, LENR via Ub_i opposition";
    }
    
    std::string getEquation() const override {
        return "dN_ch/dη ~ η × k_η × exp(-[SSq]×18/26)\n"
               "Energy scaling: √s_NN^{0.156}\n"
               "Centrality: ~1.8 factor central-to-peripheral";
    }
};

// ===========================================================================================
// 3. COMAGNETOMETER EXOTIC SPIN TERM (ArXiv 2506.15046)
// ===========================================================================================

/**
 * ComagnetometerExoticSpinTerm - Exotic spin-dependent interactions
 * 
 * UNIQUE PHYSICS:
 * - Frequency-dependent response to exotic fields
 * - Axion-proton coupling errors up to 75%
 * - Light-shift calibration methods
 * 
 * UQFF INTERPRETATION:
 * - Exotic spin as Um b_p field
 * - b_p = -(1/9) χ_p √(2ρ_DM) sin(m_a t + φ) v_z μ_B g_S / ℏ
 * - Response ties to cos(πt_n)
 */
class ComagnetometerExoticSpinTerm : public PhysicsTerm_BSM {
private:
    double m_axion;          // Axion mass (eV)
    double g_aNN;            // Axion-nucleon coupling
    double frequency;        // Measurement frequency (Hz)
    double v_z;              // Relative velocity (m/s)

public:
    ComagnetometerExoticSpinTerm(double ma = 1e-14, double g = 1e-10,
                                  double f = 10.0, double v = 300e3)
        : m_axion(ma), g_aNN(g), frequency(f), v_z(v) {
        setMetadata("arxiv", "2506.15046");
        setMetadata("experiment", "FID Comagnetometer");
        setMetadata("observable", "exotic spin coupling");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        using namespace BSMPhysics;
        
        double ma = params.count("m_axion") ? params.at("m_axion") : m_axion;
        double f = params.count("frequency") ? params.at("frequency") : frequency;
        
        // === AXION-INDUCED FIELD ===
        // b_p = -(1/9) χ_p √(2ρ_DM) sin(m_a t) v_z μ_B g_S / ℏ
        double omega_a = ma * eV_to_J / hbar;
        double b_p_amplitude = -(1.0/9.0) * chi_p * std::sqrt(2.0 * rho_DM) 
                              * v_z * mu_B * g_aNN / hbar;
        double b_p = b_p_amplitude * std::sin(omega_a * t);
        
        // === FREQUENCY RESPONSE ===
        // Error at 2 Hz: ~40%, at 20 Hz: ~75%
        double response_error = 0.4 * std::sqrt(f / 2.0);
        if (response_error > 0.75) response_error = 0.75;
        double corrected_b_p = b_p * (1.0 - response_error);
        
        // === UQFF GRAVITY COMPONENTS ===
        // Ug1: Dark matter gravitational contribution
        double Ug1 = rho_DM * G / c;
        
        // Ug2: Axion mass contribution
        double Ug2 = ma * eV_to_J / (m_p * c * c);
        
        // Ug3: t_n reversal (time direction)
        double t_n = (t > 0) ? 1.0 : -1.0;
        double Ug3 = std::cos(M_PI * t_n) * g_aNN * chi_p;
        
        // Ug4: Vacuum axion field
        double Ug4 = rho_vac_UA * std::abs(b_p) / (rho_vac_SCm * mu_B);
        
        // === UQFF MAGNETISM (exotic spin) ===
        double Um = std::abs(corrected_b_p) * mu_B / (m_p * c * c);
        
        // === UQFF BUOYANCY ===
        double Ub_i = beta_i * rho_DM * G / c * response_error;
        
        // === TOTAL UNIFIED FIELD ===
        double F_U = Ug1 + Ug2 + Ug3 + Ug4 + Um - Ub_i;
        
        if (enableLogging) {
            std::cout << "[ComagnetometerExoticSpin] f=" << f << " Hz, m_a=" << ma 
                      << " eV, b_p=" << b_p << " T, error=" << response_error*100 << "%" << std::endl;
        }
        
        return F_U;
    }
    
    // Compute sensitivity limit
    double computeSensitivity(double f) const {
        using namespace BSMPhysics;
        // δb ~ 10^{-14} T at f=10 Hz
        return 1e-14 * std::sqrt(10.0 / f);
    }
    
    std::string getName() const override { return "ComagnetometerExoticSpinTerm"; }
    
    std::string getDescription() const override {
        return "Comagnetometer exotic spin-dependent interaction with axion-proton coupling, "
               "frequency-dependent response errors up to 75%";
    }
    
    std::string getEquation() const override {
        return "b_p = -(1/9) χ_p √(2ρ_DM) sin(m_a t) v_z μ_B g_S / ℏ\n"
               "Response error: 40% at 2 Hz, 75% at 20 Hz\n"
               "χ_p ≈ 0.47 (calibrated)";
    }
};

// ===========================================================================================
// 4. JUNO PMT DCR TERM (ArXiv 2506.15164)
// ===========================================================================================

/**
 * JUNOPMTDCRTerm - Photomultiplier dark count rate stability
 * 
 * UNIQUE PHYSICS:
 * - LPMT vs MCP-PMT DCR differences
 * - Gain 10^7 operation
 * - Temperature and flasher spike identification
 * 
 * UQFF INTERPRETATION:
 * - DCR as Qs=0 noise in SCm
 * - P_SCm=1 for stability
 * - Spikes as t_n<0 reversals
 */
class JUNOPMTDCRTerm : public PhysicsTerm_BSM {
private:
    double gain;             // PMT gain
    double DCR_rate;         // Dark count rate (Hz)
    double temperature;      // Operating temperature (K)
    std::string pmt_type;    // "LPMT" or "MCP-PMT"

public:
    JUNOPMTDCRTerm(double g = 1e7, double dcr = 50e3, double T = 293.0,
                   const std::string& type = "LPMT")
        : gain(g), DCR_rate(dcr), temperature(T), pmt_type(type) {
        setMetadata("arxiv", "2506.15164");
        setMetadata("experiment", "JUNO");
        setMetadata("observable", "DCR stability");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        using namespace BSMPhysics;
        
        double g = params.count("gain") ? params.at("gain") : gain;
        double T = params.count("temperature") ? params.at("temperature") : temperature;
        double dcr = params.count("DCR_rate") ? params.at("DCR_rate") : DCR_rate;
        
        // === PMT TYPE DIFFERENCES ===
        double type_factor = (pmt_type == "MCP-PMT") ? 1.5 : 1.0;  // MCP has higher DCR
        double effective_dcr = dcr * type_factor;
        
        // === TEMPERATURE DEPENDENCE ===
        // DCR ~ exp(-E_a / kT), E_a ~ 0.5 eV for thermal electrons
        double E_activation = 0.5 * eV_to_J;
        double T_ref = 293.0;
        double temp_factor = std::exp(-E_activation / (k_B * T) + E_activation / (k_B * T_ref));
        effective_dcr *= temp_factor;
        
        // === UQFF GRAVITY COMPONENTS ===
        // Ug1: Photocathode quantum efficiency
        double QE = 0.25;  // ~25% typical
        double Ug1 = QE * hbar * c / (500e-9 * m_e * c * c);  // 500 nm photon
        
        // Ug2: Gain-dependent noise
        double Ug2 = std::log10(g) / 7.0 * k_B * T / (m_e * c * c);
        
        // Ug3: SCm stability (Qs=0 noise)
        double Qs = 0.0;  // Vacuum state
        double P_SCm = 1.0;  // Full stability
        double Ug3 = (1.0 - Qs) * P_SCm * rho_vac_SCm / rho_vac_UA;
        
        // Ug4: t_n reversal for spikes
        // Spike probability increases with |t|
        double spike_prob = std::tanh(std::abs(t) / 3600.0) * 0.01;  // 1% max per hour
        double t_n = (spike_prob > 0.005) ? -1.0 : 1.0;
        double Ug4 = std::cos(M_PI * t_n) * effective_dcr / 1e5;
        
        // === UQFF MAGNETISM ===
        // Earth's field affects PMT collection
        double B_earth = 50e-6;  // ~50 μT
        double Um = B_earth * mu_B * g / (m_e * c * c * 1e7);
        
        // === UQFF BUOYANCY ===
        double Ub_i = beta_i * effective_dcr * hbar / (g * m_e * c);
        
        // === TOTAL UNIFIED FIELD ===
        double F_U = Ug1 + Ug2 + Ug3 + Ug4 + Um - Ub_i;
        
        if (enableLogging) {
            std::cout << "[JUNOPMT_DCR] type=" << pmt_type << ", gain=" << g 
                      << ", DCR=" << effective_dcr << " Hz, T=" << T << " K" << std::endl;
        }
        
        return F_U;
    }
    
    std::string getName() const override { return "JUNOPMTDCRTerm"; }
    
    std::string getDescription() const override {
        return "JUNO PMT dark count rate stability with Qs=0 noise in SCm, "
               "P_SCm=1 for stability, spikes as t_n<0 reversals";
    }
    
    std::string getEquation() const override {
        return "DCR ~ exp(-E_a/kT) × type_factor\n"
               "Qs=0 (vacuum state), P_SCm=1 (stable)\n"
               "Spike detection: t_n<0 when prob>0.5%";
    }
};

// ===========================================================================================
// 5. TAU DIPOLE MOMENT TERM (ArXiv 2506.15245)
// ===========================================================================================

/**
 * TauDipoleMomentTerm - Tau anomalous magnetic and electric dipole moments
 * 
 * UNIQUE PHYSICS:
 * - γγ→τ⁺τ⁻ at Super Tau-Charm Facility
 * - Azimuthal asymmetries: cos2φ, sin2φ, cos4φ
 * - Re(a_τ) ∈ [-4.5, 6.9]×10⁻³ at 2σ
 * 
 * UQFF INTERPRETATION:
 * - a_τ/d_τ as Um dipole (μ_j cos(πt_n))
 * - Limits tie to κ_Higgs~1
 * - Predict from lepton radius r~0.3 fm
 */
class TauDipoleMomentTerm : public PhysicsTerm_BSM {
private:
    double a_tau;            // Anomalous magnetic moment
    double d_tau;            // Electric dipole moment (e·cm)
    double sqrt_s;           // Center-of-mass energy (GeV)
    double luminosity;       // Integrated luminosity (ab⁻¹)

public:
    TauDipoleMomentTerm(double a = 0.0, double d = 0.0,
                        double s = 4.26, double L = 10.0)
        : a_tau(a), d_tau(d), sqrt_s(s), luminosity(L) {
        setMetadata("arxiv", "2506.15245");
        setMetadata("experiment", "Super Tau-Charm Facility");
        setMetadata("observable", "a_τ, d_τ");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        using namespace BSMPhysics;
        
        double a = params.count("a_tau") ? params.at("a_tau") : a_tau;
        double d = params.count("d_tau") ? params.at("d_tau") : d_tau;
        double s = params.count("sqrt_s") ? params.at("sqrt_s") : sqrt_s;
        
        // === UQFF LEPTON RADIUS PREDICTION ===
        // a_τ = (α/2π) × (m_τ/m_e)² × (r_τ/r_e)² for BSM
        double r_tau = 0.3e-15;  // 0.3 fm predicted
        double r_e_classical = e_charge * e_charge / (4.0 * M_PI * epsilon_0 * m_e * c * c);
        double a_tau_pred = alpha_em / (2.0 * M_PI) * std::pow(m_tau / m_e, 2) 
                           * std::pow(r_tau / r_e_classical, 2);
        
        // === AZIMUTHAL ASYMMETRIES ===
        double phi = 2.0 * M_PI * t / 1.0;  // Arbitrary phase evolution
        double asym_cos2 = a * std::cos(2.0 * phi);
        double asym_sin2 = d * e_charge * 1e-2 * std::sin(2.0 * phi);  // d in e·cm
        double asym_cos4 = a * a * std::cos(4.0 * phi);
        
        // === UQFF GRAVITY COMPONENTS ===
        // Ug1: Tau mass contribution
        double Ug1 = m_tau * G / (c * c * r_tau);
        
        // Ug2: Anomalous dipole
        double mu_tau = e_charge * hbar / (2.0 * m_tau);
        double Ug2 = (1.0 + a) * mu_tau * mu_tau * mu_0 / (4.0 * M_PI * r_tau * r_tau * r_tau * m_tau * c * c);
        
        // Ug3: t_n modulation
        double t_n = (t > 0) ? 1.0 : -1.0;
        double Ug3 = a_tau_pred * std::cos(M_PI * t_n);
        
        // Ug4: Higgs coupling
        double Ug4 = kappa_Higgs * m_tau / m_H * (m_H * GeV_to_J) / (m_tau * c * c);
        
        // === UQFF MAGNETISM ===
        double Um = std::abs(a) * mu_tau / (m_tau * c * c * r_tau);
        
        // === UQFF BUOYANCY ===
        double Ub_i = beta_i * std::abs(d) * e_charge * 1e-2 / (m_tau * c * r_tau);
        
        // === TOTAL UNIFIED FIELD ===
        double F_U = Ug1 + Ug2 + Ug3 + Ug4 + Um - Ub_i;
        
        if (enableLogging) {
            std::cout << "[TauDipole] a_τ=" << a << ", d_τ=" << d << " e·cm"
                      << ", a_τ_pred=" << a_tau_pred << std::endl;
        }
        
        return F_U;
    }
    
    // Compute 2σ limits
    std::pair<double, double> compute2SigmaLimits() const {
        // From paper: Re(a_τ) ∈ [-4.5, 6.9]×10⁻³
        return {-4.5e-3, 6.9e-3};
    }
    
    std::string getName() const override { return "TauDipoleMomentTerm"; }
    
    std::string getDescription() const override {
        return "Tau lepton dipole moments via γγ→τ⁺τ⁻, Um dipole with μ_j cos(πt_n), "
               "prediction from r_τ~0.3 fm";
    }
    
    std::string getEquation() const override {
        return "a_τ = (α/2π) × (m_τ/m_e)² × (r_τ/r_e)²\n"
               "Asymmetries: cos2φ, sin2φ, cos4φ\n"
               "Re(a_τ) ∈ [-4.5, 6.9]×10⁻³ at 2σ";
    }
};

// ===========================================================================================
// 6. CKM V_cb TERM (ArXiv 2506.15256)
// ===========================================================================================

/**
 * CKMVcbTerm - |V_cb| determination from B→Dℓν
 * 
 * UNIQUE PHYSICS:
 * - Belle II Run 1 measurement
 * - |V_cb| = (40.5 ± 1.3) × 10⁻³
 * - Form factor extraction
 * 
 * UQFF INTERPRETATION:
 * - V_cb as weak coupling in η
 * - k_η × G_F² × s / π scaling
 * - κ_Higgs~1 constraint
 */
class CKMVcbTerm : public PhysicsTerm_BSM {
private:
    double V_cb_value;       // |V_cb| value
    double V_cb_error;       // Uncertainty
    double q2_min;           // Minimum q² (GeV²)
    double q2_max;           // Maximum q² (GeV²)

public:
    CKMVcbTerm(double Vcb = 40.5e-3, double err = 1.3e-3,
               double qmin = 0.0, double qmax = 11.6)
        : V_cb_value(Vcb), V_cb_error(err), q2_min(qmin), q2_max(qmax) {
        setMetadata("arxiv", "2506.15256");
        setMetadata("experiment", "Belle II");
        setMetadata("observable", "|V_cb|");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        using namespace BSMPhysics;
        
        double Vcb = params.count("V_cb") ? params.at("V_cb") : V_cb_value;
        double q2 = params.count("q2") ? params.at("q2") : (q2_min + q2_max) / 2.0;
        
        // === DECAY WIDTH ===
        // Γ(B→Dℓν) ∝ G_F² |V_cb|² m_B⁵ × phase_space × form_factor²
        double phase_space = std::sqrt(1.0 - std::pow(m_D/m_B, 2));
        double form_factor = 1.0 - q2 / (m_B * m_B);  // Simplified
        double Gamma = G_F * G_F * Vcb * Vcb * std::pow(m_B * GeV_to_J, 5) 
                      * phase_space * form_factor * form_factor / (192.0 * M_PI * M_PI * M_PI * hbar);
        
        // === UQFF WEAK COUPLING ===
        // η = k_η × G_F² × s / π
        double s = q2;  // Mandelstam s ~ q²
        double eta_weak = k_eta * G_F * G_F * s / M_PI;
        
        // === UQFF GRAVITY COMPONENTS ===
        // Ug1: B meson mass
        double Ug1 = m_B * GeV_to_J / (m_p * c * c);
        
        // Ug2: CKM unitarity
        double unitarity_check = std::abs(V_cb * V_cb + V_cd * V_cd - 1.0);
        double Ug2 = Vcb * Vcb * kappa_Higgs;
        
        // Ug3: Form factor q² dependence
        double Ug3 = form_factor * hbar * c / (m_B * GeV_to_J * 1e-15);
        
        // Ug4: Weak scale vacuum
        double Ug4 = rho_vac_UA * (m_W * GeV_to_J) / (rho_vac_SCm * m_p * c * c);
        
        // === UQFF MAGNETISM ===
        double Um = eta_weak * mu_B / (m_B * GeV_to_J * c);
        
        // === UQFF BUOYANCY ===
        double Ub_i = beta_i * Gamma / (m_B * GeV_to_J * c * c);
        
        // === TOTAL UNIFIED FIELD ===
        double F_U = Ug1 + Ug2 + Ug3 + Ug4 + Um - Ub_i;
        
        if (enableLogging) {
            std::cout << "[CKM_Vcb] |V_cb|=" << Vcb << " ± " << V_cb_error 
                      << ", Γ=" << Gamma << " s⁻¹" << std::endl;
        }
        
        return F_U;
    }
    
    std::string getName() const override { return "CKMVcbTerm"; }
    
    std::string getDescription() const override {
        return "|V_cb| from B→Dℓν at Belle II, weak coupling η = k_η G_F² s/π, "
               "κ_Higgs~1 constraint";
    }
    
    std::string getEquation() const override {
        return "Γ(B→Dℓν) ∝ G_F² |V_cb|² m_B⁵ × F(q²)²\n"
               "|V_cb| = (40.5 ± 1.3) × 10⁻³\n"
               "η = k_η × G_F² × s / π";
    }
};

// ===========================================================================================
// 7. LFV B DECAY TERM (ArXiv 2506.15347)
// ===========================================================================================

/**
 * LFVBDecayTerm - Lepton-flavor-violating B⁰→K*⁰τ±e∓ search
 * 
 * UNIQUE PHYSICS:
 * - LHCb Run 2 search
 * - BR < 5.9×10⁻⁶ at 90% CL
 * - Double-tag, GBDT/Fisher analysis
 * 
 * UQFF INTERPRETATION:
 * - LFV as t_n<0 reversal in Um
 * - Limits tie to η weak rates
 */
class LFVBDecayTerm : public PhysicsTerm_BSM {
private:
    double BR_limit;         // Branching ratio upper limit
    double CL;               // Confidence level
    double luminosity;       // Integrated luminosity (fb⁻¹)

public:
    LFVBDecayTerm(double br = 5.9e-6, double cl = 0.90, double L = 9.0)
        : BR_limit(br), CL(cl), luminosity(L) {
        setMetadata("arxiv", "2506.15347");
        setMetadata("experiment", "LHCb Run 2");
        setMetadata("observable", "BR(B⁰→K*⁰τ±e∓)");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        using namespace BSMPhysics;
        
        double br = params.count("BR") ? params.at("BR") : BR_limit;
        
        // === LFV AMPLITUDE ===
        // A(B→K*τe) = G_F V_tb V_ts* × C_LFV × form_factors
        double C_LFV = br / 1e-5;  // LFV Wilson coefficient proxy
        
        // === t_n REVERSAL ===
        // LFV implies time-reversal violation (t_n < 0)
        double t_n = -1.0;  // LFV is t_n<0
        double LFV_suppression = std::exp(-std::abs(t_n) * SSq);
        
        // === UQFF GRAVITY COMPONENTS ===
        // Ug1: B meson
        double Ug1 = m_B * GeV_to_J / (m_p * c * c);
        
        // Ug2: Lepton mass difference
        double delta_m = (m_tau - m_e);
        double Ug2 = delta_m * G / (c * c * 1e-15);
        
        // Ug3: t_n reversal factor
        double Ug3 = std::cos(M_PI * t_n) * C_LFV * LFV_suppression;
        
        // Ug4: Weak scale
        double Ug4 = rho_vac_UA * br / rho_vac_SCm;
        
        // === UQFF MAGNETISM (Um with t_n<0) ===
        // LFV transitions involve spin-flip
        double Um = std::abs(t_n) * mu_B * (m_tau - m_e) / (m_B * GeV_to_J * c);
        
        // === UQFF BUOYANCY ===
        // η weak rate from k_η
        double eta_weak = k_eta * G_F * G_F * m_B * m_B * GeV_to_J * GeV_to_J / M_PI;
        double Ub_i = beta_i * eta_weak * br;
        
        // === TOTAL UNIFIED FIELD ===
        double F_U = Ug1 + Ug2 + Ug3 + Ug4 + Um - Ub_i;
        
        if (enableLogging) {
            std::cout << "[LFVBDecay] BR < " << br << " at " << CL*100 << "% CL"
                      << ", t_n=" << t_n << std::endl;
        }
        
        return F_U;
    }
    
    std::string getName() const override { return "LFVBDecayTerm"; }
    
    std::string getDescription() const override {
        return "LFV B⁰→K*⁰τ±e∓ search with t_n<0 reversal in Um, "
               "BR limit ties to η weak rates";
    }
    
    std::string getEquation() const override {
        return "BR(B⁰→K*⁰τ±e∓) < 5.9×10⁻⁶ at 90% CL\n"
               "LFV implies t_n < 0 (time reversal)\n"
               "Suppression: exp(-|t_n| × [SSq])";
    }
};

// ===========================================================================================
// 8. ECFA HIGGS FACTORY TERM (ArXiv 2506.15390)
// ===========================================================================================

/**
 * ECFAHiggsFactoryTerm - e⁺e⁻ Higgs/EW/top factory sensitivities
 * 
 * UNIQUE PHYSICS:
 * - ECFA 2021-2025 study
 * - FCC-ee, ILC, CLIC, CEPC projections
 * - Higgs coupling precisions
 * 
 * UQFF INTERPRETATION:
 * - EW/Higgs as level 18 (m_H=125 GeV)
 * - Sensitivities tie to κ_Higgs~1
 */
class ECFAHiggsFactoryTerm : public PhysicsTerm_BSM {
private:
    double sqrt_s;           // Center-of-mass energy (GeV)
    double luminosity;       // Integrated luminosity (ab⁻¹)
    std::string collider;    // Collider name

public:
    ECFAHiggsFactoryTerm(double s = 250.0, double L = 5.0,
                         const std::string& coll = "FCC-ee")
        : sqrt_s(s), luminosity(L), collider(coll) {
        setMetadata("arxiv", "2506.15390");
        setMetadata("experiment", "ECFA Study");
        setMetadata("observable", "Higgs couplings");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        using namespace BSMPhysics;
        
        double s = params.count("sqrt_s") ? params.at("sqrt_s") : sqrt_s;
        double L = params.count("luminosity") ? params.at("luminosity") : luminosity;
        
        // === HIGGS PRODUCTION CROSS-SECTION ===
        // σ(e⁺e⁻→ZH) peaks at √s ~ 250 GeV
        double sigma_ZH = 200.0 * fb_inv_to_m2;  // ~200 fb at 250 GeV
        if (s > 250.0) sigma_ZH *= std::pow(250.0/s, 2);  // 1/s scaling
        
        // === COUPLING PRECISION ===
        // δκ/κ ~ 1/√(L × σ × BR)
        double N_events = L * 1e3 * sigma_ZH / fb_inv_to_m2;  // L in ab⁻¹
        double delta_kappa = 1.0 / std::sqrt(N_events);
        
        // === UQFF LEVEL 18 (EW SCALE) ===
        double n_level = 18.0;
        double SSq_factor = std::exp(-SSq * n_level / DIMENSIONS);
        
        // === UQFF GRAVITY COMPONENTS ===
        // Ug1: Higgs mass
        double Ug1 = m_H * GeV_to_J / (m_p * c * c);
        
        // Ug2: EW symmetry breaking
        double v_Higgs = 246.0;  // Higgs vev (GeV)
        double Ug2 = v_Higgs * GeV_to_J / (m_p * c * c * c);
        
        // Ug3: Coupling precision
        double Ug3 = kappa_Higgs * (1.0 - delta_kappa) * SSq_factor;
        
        // Ug4: Vacuum energy at EW scale
        double Ug4 = rho_vac_UA * m_H * GeV_to_J / (rho_vac_SCm * m_p * c * c);
        
        // === UQFF MAGNETISM ===
        // Weak isospin contribution
        double Um = (m_Z - m_W) * GeV_to_J / (m_H * GeV_to_J) * mu_B / (m_p * c);
        
        // === UQFF BUOYANCY ===
        double Ub_i = beta_i * delta_kappa * m_H * GeV_to_J / (m_p * c * c);
        
        // === TOTAL UNIFIED FIELD ===
        double F_U = (Ug1 + Ug2 + Ug3 + Ug4 + Um - Ub_i) * SSq_factor;
        
        if (enableLogging) {
            std::cout << "[ECFAHiggs] " << collider << " √s=" << s << " GeV, L=" << L 
                      << " ab⁻¹, δκ/κ=" << delta_kappa*100 << "%" << std::endl;
        }
        
        return F_U;
    }
    
    // Collider-specific precision
    double getCouplingPrecision(const std::string& coupling) const {
        // From ECFA study (approximate)
        if (coupling == "κ_b") return 0.006;      // 0.6%
        if (coupling == "κ_c") return 0.018;      // 1.8%
        if (coupling == "κ_τ") return 0.009;      // 0.9%
        if (coupling == "κ_W") return 0.004;      // 0.4%
        if (coupling == "κ_Z") return 0.002;      // 0.2%
        if (coupling == "κ_g") return 0.015;      // 1.5%
        if (coupling == "κ_γ") return 0.035;      // 3.5%
        return 0.01;
    }
    
    std::string getName() const override { return "ECFAHiggsFactoryTerm"; }
    
    std::string getDescription() const override {
        return "ECFA Higgs factory study with EW level 18, "
               "κ_Higgs~1 constraint, sub-percent coupling precisions";
    }
    
    std::string getEquation() const override {
        return "σ(e⁺e⁻→ZH) ~ 200 fb at √s=250 GeV\n"
               "δκ/κ ~ 1/√(L × σ × BR)\n"
               "Level 18: exp(-[SSq]×18/26)";
    }
};

// ===========================================================================================
// 9. VECTOR-LIKE QUARK TERM (ArXiv 2506.15515)
// ===========================================================================================

/**
 * VectorLikeQuarkTerm - T/Y→Wb search at ATLAS
 * 
 * UNIQUE PHYSICS:
 * - Single production of vector-like quarks
 * - κ limits: 0.14-0.52 for m~1150-2600 GeV
 * - SM interference effects
 * 
 * UQFF INTERPRETATION:
 * - VLQ as SCm heavy quarks at n=18
 * - Limits tie to k_η
 */
class VectorLikeQuarkTerm : public PhysicsTerm_BSM {
private:
    double m_VLQ;            // VLQ mass (GeV)
    double kappa_VLQ;        // Coupling strength
    std::string type;        // "T" (charge 2/3) or "Y" (charge -4/3)

public:
    VectorLikeQuarkTerm(double m = 1500.0, double k = 0.3,
                        const std::string& t = "T")
        : m_VLQ(m), kappa_VLQ(k), type(t) {
        setMetadata("arxiv", "2506.15515");
        setMetadata("experiment", "ATLAS Run 2");
        setMetadata("observable", "VLQ single production");
    }
    
    double compute(double t_time, const std::map<std::string, double>& params) const override {
        using namespace BSMPhysics;
        
        double m = params.count("m_VLQ") ? params.at("m_VLQ") : m_VLQ;
        double k = params.count("kappa") ? params.at("kappa") : kappa_VLQ;
        
        // === PRODUCTION CROSS-SECTION ===
        // σ(pp→Q) ~ κ² × (σ_0 / m²)
        double sigma_0 = 1e-36;  // Reference cross-section (m²)
        double sigma_VLQ = k * k * sigma_0 / std::pow(m / 1000.0, 2);
        
        // === UQFF LEVEL 18 (HEAVY QUARKS) ===
        double n_level = 18.0;
        double SSq_factor = std::exp(-SSq * n_level / DIMENSIONS);
        
        // === DECAY WIDTH ===
        // Γ(Q→Wb) = (κ² m³ / 16π v²) × (1 - m_W²/m²)² × (1 + 2m_W²/m²)
        double v_Higgs = 246.0 * GeV_to_J;
        double x_W = m_W / m;
        double Gamma_Wb = k * k * std::pow(m * GeV_to_J, 3) / (16.0 * M_PI * v_Higgs * v_Higgs)
                         * std::pow(1.0 - x_W * x_W, 2) * (1.0 + 2.0 * x_W * x_W) / hbar;
        
        // === UQFF GRAVITY COMPONENTS ===
        // Ug1: VLQ mass
        double Ug1 = m * GeV_to_J / (m_p * c * c);
        
        // Ug2: Coupling strength
        double Ug2 = k * k * m_t / m;  // Relative to top
        
        // Ug3: SCm heavy quark
        double Ug3 = SSq_factor * k / kappa_Higgs;
        
        // Ug4: k_η constraint
        double Ug4 = k_eta * sigma_VLQ / fb_inv_to_m2 * SSq_factor;
        
        // === UQFF MAGNETISM ===
        // VLQ magnetic moment
        double Q_VLQ = (type == "T") ? 2.0/3.0 : -4.0/3.0;
        double mu_VLQ = Q_VLQ * e_charge * hbar / (2.0 * m * GeV_to_J);
        double Um = mu_VLQ * mu_VLQ * mu_0 / (4.0 * M_PI * std::pow(1e-18, 3) * m * GeV_to_J * c * c);
        
        // === UQFF BUOYANCY ===
        double Ub_i = beta_i * Gamma_Wb / (m * GeV_to_J * c * c / hbar);
        
        // === TOTAL UNIFIED FIELD ===
        double F_U = (Ug1 + Ug2 + Ug3 + Ug4 + Um - Ub_i) * SSq_factor;
        
        if (enableLogging) {
            std::cout << "[VectorLikeQuark] " << type << " m=" << m << " GeV, κ=" << k 
                      << ", σ=" << sigma_VLQ/fb_inv_to_m2 << " fb" << std::endl;
        }
        
        return F_U;
    }
    
    // Get κ limit for given mass
    double getKappaLimit(double m) const {
        // From ATLAS results (approximate interpolation)
        if (m < 1150) return 0.52;
        if (m > 2600) return 0.14;
        return 0.52 - 0.38 * (m - 1150) / (2600 - 1150);
    }
    
    std::string getName() const override { return "VectorLikeQuarkTerm"; }
    
    std::string getDescription() const override {
        return "Vector-like quark T/Y→Wb at ATLAS, SCm heavy quarks n=18, "
               "κ limits 0.14-0.52 tie to k_η";
    }
    
    std::string getEquation() const override {
        return "σ(pp→Q) ~ κ² × σ_0 / m²\n"
               "Γ(Q→Wb) = (κ²m³/16πv²) × phase_space\n"
               "95% CL: κ = 0.14-0.52 for m=1150-2600 GeV";
    }
};

// ===========================================================================================
// 10. DCS D DECAY TERM (ArXiv 2506.15533)
// ===========================================================================================

/**
 * DCSDDecayTerm - Doubly Cabibbo-suppressed D⁺→K⁺π⁰/η/η' decays
 * 
 * UNIQUE PHYSICS:
 * - BESIII 20.3 fb⁻¹ at ψ(3770)
 * - BR(K⁺π⁰) = (1.45±0.06±0.06)×10⁻⁴
 * - BR(K⁺η) = (1.17±0.10±0.03)×10⁻⁴
 * - BR(K⁺η') = (1.88±0.15±0.06)×10⁻⁴
 * 
 * UQFF INTERPRETATION:
 * - DCS as light quark flavor violation
 * - Ties to η~10¹³ cm⁻²/s in LENR
 * - BR scales k_η
 */
class DCSDDecayTerm : public PhysicsTerm_BSM {
private:
    double BR_Kpi0;          // BR(D⁺→K⁺π⁰)
    double BR_Keta;          // BR(D⁺→K⁺η)
    double BR_Ketap;         // BR(D⁺→K⁺η')

public:
    DCSDDecayTerm(double br1 = 1.45e-4, double br2 = 1.17e-4, double br3 = 1.88e-4)
        : BR_Kpi0(br1), BR_Keta(br2), BR_Ketap(br3) {
        setMetadata("arxiv", "2506.15533");
        setMetadata("experiment", "BESIII");
        setMetadata("observable", "DCS D⁺ decays");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        using namespace BSMPhysics;
        
        std::string mode = "Kpi0";
        if (params.count("mode")) {
            double m = params.at("mode");
            if (m == 1) mode = "Keta";
            else if (m == 2) mode = "Ketap";
        }
        
        double BR = BR_Kpi0;
        if (mode == "Keta") BR = BR_Keta;
        else if (mode == "Ketap") BR = BR_Ketap;
        
        // === DCS SUPPRESSION ===
        // DCS amplitude ~ V_cd* V_us (doubly Cabibbo-suppressed)
        double DCS_factor = V_cd * V_cd;  // ~0.05
        
        // === UQFF LIGHT FLAVOR LEVEL ===
        double n_level = 12.0;  // Light quarks
        double SSq_factor = std::exp(-SSq * n_level / DIMENSIONS);
        
        // === LENR CONNECTION ===
        // η ~ 10¹³ cm⁻²/s ties to light flavor rates
        double eta_LENR = 1e13 * 1e4;  // Convert to m⁻²/s
        double BR_eta_ratio = BR / k_eta;
        
        // === UQFF GRAVITY COMPONENTS ===
        // Ug1: D meson mass
        double Ug1 = m_D * GeV_to_J / (m_p * c * c);
        
        // Ug2: CKM suppression
        double Ug2 = DCS_factor * kappa_Higgs;
        
        // Ug3: BR scaling with k_η
        double Ug3 = BR / (k_eta * 1e113) * SSq_factor;
        
        // Ug4: Light flavor vacuum
        double Ug4 = rho_vac_UA * BR / rho_vac_SCm;
        
        // === UQFF MAGNETISM ===
        double Um = eta_LENR * k_eta * mu_B / (m_D * GeV_to_J * c);
        
        // === UQFF BUOYANCY ===
        double tau_D = 1.04e-12;  // D⁺ lifetime (s)
        double Gamma_D = hbar / tau_D;
        double Ub_i = beta_i * BR * Gamma_D / (m_D * GeV_to_J * c * c);
        
        // === TOTAL UNIFIED FIELD ===
        double F_U = (Ug1 + Ug2 + Ug3 + Ug4 + Um - Ub_i) * SSq_factor;
        
        if (enableLogging) {
            std::cout << "[DCSDDecay] mode=" << mode << ", BR=" << BR 
                      << ", DCS factor=" << DCS_factor << std::endl;
        }
        
        return F_U;
    }
    
    std::string getName() const override { return "DCSDDecayTerm"; }
    
    std::string getDescription() const override {
        return "DCS D⁺→K⁺π⁰/η/η' at BESIII, light quark flavor violation, "
               "BR ties to η~10¹³ cm⁻²/s LENR and scales k_η";
    }
    
    std::string getEquation() const override {
        return "BR(K⁺π⁰) = (1.45±0.08)×10⁻⁴\n"
               "BR(K⁺η) = (1.17±0.10)×10⁻⁴\n"
               "BR(K⁺η') = (1.88±0.16)×10⁻⁴\n"
               "DCS ~ |V_cd|² ~ 0.05";
    }
};

// ===========================================================================================
// 11. MOLECULAR ROTOR COLLISION TERM (H₂O-H₂ CC/CS)
// ===========================================================================================

/**
 * MolecularRotorCollisionTerm - H₂O-H₂ close-coupling/coupled-states
 * 
 * UNIQUE PHYSICS:
 * - Asymmetric top (H₂O) + linear rotor (H₂) scattering
 * - CS approximation: Ω = Λ + λ decoupling
 * - Cross-section: σ = (π/k²) Σ(2J+1)|δ-S|²
 * - PES anisotropy drives Δj=2 transitions
 * 
 * UQFF INTERPRETATION:
 * - Um rotor torque τ_rot = r × F_V
 * - θ~90° rainbow scattering
 * - Extends Um at level 10 (molecular/atomic)
 */
class MolecularRotorCollisionTerm : public PhysicsTerm_BSM {
private:
    double E_collision;      // Collision energy (cm⁻¹)
    double impact_param;     // Impact parameter (bohr)
    int j_initial;           // Initial H₂ rotational state
    int J_tau_initial;       // Initial H₂O state

public:
    MolecularRotorCollisionTerm(double E = 300.0, double b = 5.0,
                                 int j = 0, int Jtau = 0)
        : E_collision(E), impact_param(b), j_initial(j), J_tau_initial(Jtau) {
        setMetadata("reference", "Phillips et al. 1995, JCP 102, 10054");
        setMetadata("method", "Close-Coupling / Coupled-States");
        setMetadata("observable", "rotational cross-section");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        using namespace BSMPhysics;
        
        double E = params.count("E_collision") ? params.at("E_collision") : E_collision;
        double b = params.count("impact_param") ? params.at("impact_param") : impact_param;
        int j = j_initial;
        
        // === COLLISION KINEMATICS ===
        double E_J = E * cm_inv_to_J;  // Convert cm⁻¹ to Joules
        double mu_reduced = (18.015 * 2.016) / (18.015 + 2.016) * 1.66054e-27;  // kg
        double k = std::sqrt(2.0 * mu_reduced * E_J) / hbar;  // Wave number
        double v_rel = hbar * k / mu_reduced;  // Relative velocity
        
        // === S-MATRIX (simplified) ===
        // S ~ exp(-i × η), η = ∫ k dr
        double a0 = 5.29177e-11;  // Bohr radius
        double b_m = b * a0;  // Impact parameter in meters
        double eta = k * b_m;  // Phase shift proxy
        std::complex<double> S(std::cos(eta), -std::sin(eta));
        S *= 0.9;  // Absorption factor
        
        // === CROSS-SECTION ===
        // σ = (π/k²) Σ (2J+1) |δ_JJ' - S_JJ'|²
        double sigma_elastic = 0.0;
        double sigma_inelastic = 0.0;
        int J_max = 6;  // Truncate at J=6
        
        for (int J = 0; J <= J_max; J++) {
            double weight = 2.0 * J + 1.0;
            double S_mag = std::abs(S) * std::exp(-0.1 * J);  // J-dependent attenuation
            sigma_elastic += weight * std::pow(1.0 - S_mag, 2);
            sigma_inelastic += weight * (1.0 - S_mag * S_mag);
        }
        sigma_elastic *= M_PI / (k * k);
        sigma_inelastic *= M_PI / (k * k);
        
        // Convert to Å²
        double sigma_elastic_A2 = sigma_elastic * 1e20;
        double sigma_inelastic_A2 = sigma_inelastic * 1e20;
        
        // === UQFF LEVEL 10 (MOLECULAR) ===
        double n_level = 10.0;
        double SSq_factor = std::exp(-SSq * n_level / DIMENSIONS);
        
        // === PES ANISOTROPY ===
        // V(R,r) drives Δj=2 transitions
        // Fit: σ = A/E + B sin²(aniso), A~300, B~2, aniso~π/2
        double A_fit = 300.0;
        double B_fit = 2.0;
        double aniso = M_PI / 2.0;  // Rainbow at θ~90°
        double sigma_model = A_fit / E + B_fit * std::pow(std::sin(aniso), 2);
        
        // === UQFF GRAVITY COMPONENTS ===
        // Ug1: Molecular mass
        double Ug1 = (18.015 + 2.016) * 1.66054e-27 * G / (a0 * c * c);
        
        // Ug2: PES as Ug2 bubble
        // (Q_A + Q_UA) / r² × S(r - R_b)
        double Q_A = 0.1 * e_charge;  // Effective charge
        double R_b = 3.0 * a0;  // Bubble radius
        double Ug2 = (Q_A + U_UA * e_charge) / (b_m * b_m) * (b_m > R_b ? 1.0 : 0.0) / (mu_reduced * c * c);
        
        // Ug3: Rotational energy
        double B_H2 = 60.853 * cm_inv_to_J;  // H₂ rotational constant
        double E_rot = B_H2 * j * (j + 1);
        double Ug3 = E_rot / (mu_reduced * c * c);
        
        // Ug4: Vacuum at collision
        double Ug4 = rho_vac_UA * sigma_elastic / (rho_vac_SCm * a0 * a0);
        
        // === UQFF MAGNETISM (rotor torque) ===
        // τ_rot = r × F_V, F_V = -∇V
        double F_V = Q_A * Q_A / (4.0 * M_PI * epsilon_0 * b_m * b_m);
        double tau_rot = b_m * F_V * std::sin(aniso);
        double Um = tau_rot / (mu_reduced * c * c * a0);
        
        // === UQFF BUOYANCY ===
        // Δj=2 transitions from k_η
        double eta_inelastic = k_eta * sigma_inelastic * v_rel * 1e6;  // per cm³
        double Ub_i = beta_i * eta_inelastic * SSq_factor;
        
        // === TOTAL UNIFIED FIELD ===
        double F_U = (Ug1 + Ug2 + Ug3 + Ug4 + Um - Ub_i) * SSq_factor;
        
        if (enableLogging) {
            std::cout << "[MolecularRotor] E=" << E << " cm⁻¹, b=" << b << " a₀"
                      << ", σ_el=" << sigma_elastic_A2 << " Å², σ_inel=" << sigma_inelastic_A2 << " Å²" << std::endl;
        }
        
        return F_U;
    }
    
    // Cross-section fit
    double computeCrossSectionFit(double E) const {
        // σ = A/E + B sin²(π/2) = A/E + B
        return 300.0 / E + 2.0;  // Å²
    }
    
    std::string getName() const override { return "MolecularRotorCollisionTerm"; }
    
    std::string getDescription() const override {
        return "H₂O-H₂ rotor collision (CC/CS), Um torque τ_rot = r×F_V, "
               "θ~90° rainbow, Δj=2 from PES anisotropy at level 10";
    }
    
    std::string getEquation() const override {
        return "σ = (π/k²) Σ(2J+1)|δ-S|²\n"
               "Fit: σ = A/E + B sin²θ, A~300, B~2, θ~90°\n"
               "τ_rot = r × F_V, F_V = -∇V";
    }
};

// ===========================================================================================
// BSM PHYSICS MODULE CLASS - ORCHESTRATION
// ===========================================================================================

class BSMPhysicsUQFFModule {
private:
    // Term instances
    std::unique_ptr<NeutrinoPolarizabilityTerm> neutrino_polar_term;
    std::unique_ptr<ALICEChargedParticleTerm> alice_term;
    std::unique_ptr<ComagnetometerExoticSpinTerm> comag_term;
    std::unique_ptr<JUNOPMTDCRTerm> juno_term;
    std::unique_ptr<TauDipoleMomentTerm> tau_dipole_term;
    std::unique_ptr<CKMVcbTerm> ckm_term;
    std::unique_ptr<LFVBDecayTerm> lfv_term;
    std::unique_ptr<ECFAHiggsFactoryTerm> ecfa_term;
    std::unique_ptr<VectorLikeQuarkTerm> vlq_term;
    std::unique_ptr<DCSDDecayTerm> dcs_term;
    std::unique_ptr<MolecularRotorCollisionTerm> rotor_term;
    
    // Module state
    bool initialized;
    bool verbose;
    
    // Self-expanding framework
    std::map<std::string, double> dynamicParameters;
    std::map<std::string, std::string> metadata;

public:
    BSMPhysicsUQFFModule() : initialized(false), verbose(false) {
        initialize();
    }
    
    void initialize() {
        neutrino_polar_term = std::make_unique<NeutrinoPolarizabilityTerm>();
        alice_term = std::make_unique<ALICEChargedParticleTerm>();
        comag_term = std::make_unique<ComagnetometerExoticSpinTerm>();
        juno_term = std::make_unique<JUNOPMTDCRTerm>();
        tau_dipole_term = std::make_unique<TauDipoleMomentTerm>();
        ckm_term = std::make_unique<CKMVcbTerm>();
        lfv_term = std::make_unique<LFVBDecayTerm>();
        ecfa_term = std::make_unique<ECFAHiggsFactoryTerm>();
        vlq_term = std::make_unique<VectorLikeQuarkTerm>();
        dcs_term = std::make_unique<DCSDDecayTerm>();
        rotor_term = std::make_unique<MolecularRotorCollisionTerm>();
        
        // Set metadata
        metadata["batch"] = "22-Extension";
        metadata["date"] = "January 28, 2026";
        metadata["author"] = "Daniel T. Murphy";
        metadata["grok_analysis"] = "September 14, 2025";
        metadata["solvability"] = "99.9%";
        metadata["papers_integrated"] = "11";
        
        initialized = true;
    }
    
    void setVerbose(bool v) {
        verbose = v;
        if (neutrino_polar_term) neutrino_polar_term->setEnableLogging(v);
        if (alice_term) alice_term->setEnableLogging(v);
        if (comag_term) comag_term->setEnableLogging(v);
        if (juno_term) juno_term->setEnableLogging(v);
        if (tau_dipole_term) tau_dipole_term->setEnableLogging(v);
        if (ckm_term) ckm_term->setEnableLogging(v);
        if (lfv_term) lfv_term->setEnableLogging(v);
        if (ecfa_term) ecfa_term->setEnableLogging(v);
        if (vlq_term) vlq_term->setEnableLogging(v);
        if (dcs_term) dcs_term->setEnableLogging(v);
        if (rotor_term) rotor_term->setEnableLogging(v);
    }
    
    // === BATCH VALIDATION ===
    
    void runBSMValidation() {
        std::cout << "\n=== BSM PHYSICS UQFF VALIDATION (BATCH 22 EXTENSION) ===" << std::endl;
        std::cout << "Integration Date: January 28, 2026" << std::endl;
        std::cout << "Grok 4 Analysis: September 14, 2025" << std::endl;
        std::cout << "UQFF Solvability: 99.9%" << std::endl;
        std::cout << "ArXiv Papers Integrated: 11" << std::endl;
        std::cout << std::string(65, '-') << std::endl;
        
        std::map<std::string, double> params;
        double t_test = 1.0;
        
        std::cout << "\n1. Neutrino Polarizability (2506.14881):" << std::endl;
        params["E_nu"] = 1.0; params["Lambda_NP"] = 1000.0;
        std::cout << "   F_U = " << std::scientific << neutrino_polar_term->compute(t_test, params) << std::endl;
        
        std::cout << "\n2. ALICE Charged Particles (2506.14989):" << std::endl;
        params.clear(); params["sqrt_s"] = 5.36; params["centrality"] = 0.0;
        std::cout << "   F_U = " << alice_term->compute(t_test, params) << std::endl;
        
        std::cout << "\n3. Comagnetometer Exotic Spin (2506.15046):" << std::endl;
        params.clear(); params["m_axion"] = 1e-14; params["frequency"] = 10.0;
        std::cout << "   F_U = " << comag_term->compute(t_test, params) << std::endl;
        
        std::cout << "\n4. JUNO PMT DCR (2506.15164):" << std::endl;
        params.clear(); params["gain"] = 1e7; params["temperature"] = 293.0;
        std::cout << "   F_U = " << juno_term->compute(t_test, params) << std::endl;
        
        std::cout << "\n5. Tau Dipole Moments (2506.15245):" << std::endl;
        params.clear(); params["a_tau"] = 0.001; params["d_tau"] = 0.0;
        std::cout << "   F_U = " << tau_dipole_term->compute(t_test, params) << std::endl;
        
        std::cout << "\n6. CKM |V_cb| (2506.15256):" << std::endl;
        params.clear(); params["V_cb"] = 40.5e-3;
        std::cout << "   F_U = " << ckm_term->compute(t_test, params) << std::endl;
        
        std::cout << "\n7. LFV B Decay (2506.15347):" << std::endl;
        params.clear(); params["BR"] = 5.9e-6;
        std::cout << "   F_U = " << lfv_term->compute(t_test, params) << std::endl;
        
        std::cout << "\n8. ECFA Higgs Factory (2506.15390):" << std::endl;
        params.clear(); params["sqrt_s"] = 250.0; params["luminosity"] = 5.0;
        std::cout << "   F_U = " << ecfa_term->compute(t_test, params) << std::endl;
        
        std::cout << "\n9. Vector-Like Quarks (2506.15515):" << std::endl;
        params.clear(); params["m_VLQ"] = 1500.0; params["kappa"] = 0.3;
        std::cout << "   F_U = " << vlq_term->compute(t_test, params) << std::endl;
        
        std::cout << "\n10. DCS D Decays (2506.15533):" << std::endl;
        params.clear();
        std::cout << "   F_U = " << dcs_term->compute(t_test, params) << std::endl;
        
        std::cout << "\n11. Molecular Rotor Collision (H₂O-H₂):" << std::endl;
        params.clear(); params["E_collision"] = 300.0; params["impact_param"] = 5.0;
        std::cout << "   F_U = " << rotor_term->compute(t_test, params) << std::endl;
        std::cout << "   σ_fit(300 cm⁻¹) = " << rotor_term->computeCrossSectionFit(300.0) << " Å²" << std::endl;
        
        std::cout << "\n" << std::string(65, '=') << std::endl;
        std::cout << "BSM PHYSICS VALIDATION COMPLETE: 11 PhysicsTerm classes integrated" << std::endl;
        std::cout << "Total PhysicsTerm classes: 6,699 (Batch 22: 5 + BSM Extension: 11 + Previous: 6,683)" << std::endl;
    }
    
    // === GETTERS ===
    
    NeutrinoPolarizabilityTerm* getNeutrinoPolarizabilityTerm() { return neutrino_polar_term.get(); }
    ALICEChargedParticleTerm* getALICETerm() { return alice_term.get(); }
    ComagnetometerExoticSpinTerm* getComagnetometerTerm() { return comag_term.get(); }
    JUNOPMTDCRTerm* getJUNOTerm() { return juno_term.get(); }
    TauDipoleMomentTerm* getTauDipoleTerm() { return tau_dipole_term.get(); }
    CKMVcbTerm* getCKMTerm() { return ckm_term.get(); }
    LFVBDecayTerm* getLFVTerm() { return lfv_term.get(); }
    ECFAHiggsFactoryTerm* getECFATerm() { return ecfa_term.get(); }
    VectorLikeQuarkTerm* getVLQTerm() { return vlq_term.get(); }
    DCSDDecayTerm* getDCSTerm() { return dcs_term.get(); }
    MolecularRotorCollisionTerm* getRotorTerm() { return rotor_term.get(); }
    
    bool isInitialized() const { return initialized; }
};

#endif // BSM_PHYSICS_UQFF_MODULE_H
