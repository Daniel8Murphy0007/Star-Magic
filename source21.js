// ================================================================================================
// Header: NGC3603.js
//
// Description: SELF-EXPANDING JavaScript Module for NGC 3603 Extreme Star Cluster Class
//              This is the eleventh module in a series of 500+ code files for the Universal Quantum
//              Field Framework (UQFF) simulations, focusing on young massive star cluster evolution
//              and gravity equations derived from Hubble datasets, high-energy lab simulations, and
//              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
//
// Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for NGC 3603 evolution.
//          Includes ALL terms: base gravity with mass growth M(t), cosmic expansion (H_0), magnetic
//          correction (static B), UQFF Ug components with f_TRZ, Lambda, quantum uncertainty, scaled
//          EM with [UA], fluid dynamics, oscillatory waves, DM/density perturbations, stellar wind
//          feedback (pressure / density for acc), and cavity pressure P(t) / rho_fluid.
//          Supports dynamic variable updates for all parameters.
//
// Integration: Designed for inclusion in base program 'index.js'.
//              Instantiate class in main: const ngc3603 = new NGC3603();
//              Compute: const g = ngc3603.compute_g_NGC3603(t);
//
// Key Features:
//   - Default values from UQFF document: M0 = 400,000 Msun, r = 8.998e15 m (9.5 ly), tau_SF = 1 Myr,
//     rho_wind = 1e-20 kg/m^3, v_wind = 2e6 m/s, P0 = 4e-8 (Pa), tau_exp = 1 Myr, B = 1e-5 T.
//   - Units handled: Msun to kg, ly to m; pressure terms as P / rho_fluid for acceleration.
//   - Setter methods for updates: setVariable(new_val) or addToVariable(delta)/subtractFromVariable(delta).
//   - Computes g_NGC3603(r, t) with every term explicitly included.
//
// Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
// Date: October 08, 2025
//
// Enhanced: November 04, 2025 - Added self-expanding capabilities
// Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
// ================================================================================================

// ===========================================================================================
// SELF-EXPANDING FRAMEWORK: Dynamic Physics Term System
// ===========================================================================================

class PhysicsTerm {
    constructor() {
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;
    }

    compute(t, params) {
        // Base implementation - override in subclasses
        return 0.0;
    }

    getName() {
        return "BasePhysicsTerm";
    }

    getDescription() {
        return "Base physics term";
    }

    validate(params) {
        return true;
    }

    // Self-expanding methods
    registerDynamicTerm(term) {
        if (this.enableDynamicTerms) {
            this.dynamicTerms.push(term);
            if (this.enableLogging) {
                console.log(`Registered dynamic term: ${term.getName()}`);
            }
        }
    }

    setDynamicParameter(key, value) {
        this.dynamicParameters.set(key, value);
    }

    getDynamicParameter(key, defaultValue = 0.0) {
        return this.dynamicParameters.get(key) || defaultValue;
    }

    setEnableLogging(enable) {
        this.enableLogging = enable;
    }

    setLearningRate(rate) {
        this.learningRate = rate;
    }

    exportState(filename) {
        // Export current state for cross-module communication
        const state = {
            name: this.getName(),
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            metadata: Object.fromEntries(this.metadata),
            enableDynamicTerms: this.enableDynamicTerms,
            learningRate: this.learningRate
        };

        if (this.enableLogging) {
            console.log(`Exporting state for ${this.getName()} to ${filename}`);
        }

        // In a real implementation, this would write to a file
        return state;
    }
}

class DynamicVacuumTerm extends PhysicsTerm {
    constructor(amplitude = 1e-10, frequency = 1e-15) {
        super();
        this.amplitude = amplitude;
        this.frequency = frequency;
    }

    compute(t, params) {
        const rho_vac = params.get('rho_vac_UA') || 7.09e-36;
        return this.amplitude * rho_vac * Math.sin(this.frequency * t);
    }

    getName() {
        return "DynamicVacuum";
    }

    getDescription() {
        return "Time-varying vacuum energy";
    }
}

class QuantumCouplingTerm extends PhysicsTerm {
    constructor(coupling_strength = 1e-40) {
        super();
        this.coupling_strength = coupling_strength;
    }

    compute(t, params) {
        const hbar = params.get('hbar') || 1.0546e-34;
        const M = params.get('M') || 1.989e30;
        const r = params.get('r') || 1e4;
        return this.coupling_strength * (hbar * hbar) / (M * r * r) * Math.cos(t / 1e6);
    }

    getName() {
        return "QuantumCoupling";
    }

    getDescription() {
        return "Non-local quantum effects";
    }
}

// ===========================================================================================
// ENHANCED CLASS WITH SELF-EXPANDING CAPABILITIES
// ===========================================================================================

class NGC3603 extends PhysicsTerm {
    constructor() {
        super();
        this.metadata.set('enhanced', 'true');
        this.metadata.set('version', '2.0-Enhanced');
        this.metadata.set('system', 'NGC_3603');

        this.initializeDefaults();
    }

    // Initialization method (called in constructor)
    initializeDefaults() {
        // Core parameters (original UQFF - preserved)
        this.G = 6.6743e-11;              // Gravitational constant
        this.M0 = 400000.0 * 1.989e30;    // Initial mass (kg) - 400,000 Msun
        this.r = 9.5 * 9.461e15;          // Radius (m) - 9.5 ly
        this.H0 = 2.184e-18;              // Hubble constant (s^-1)
        this.B = 1e-5;                    // Static magnetic field (T)
        this.B_crit = 1e11;               // Critical B field (T)
        this.Lambda = 1.1e-52;            // Cosmological constant
        this.c_light = 3e8;               // Speed of light
        this.q_charge = 1.602e-19;        // Charge (proton)
        this.gas_v = 1e5;                 // Gas velocity for EM (m/s)
        this.f_TRZ = 0.1;                 // Time-reversal factor
        this.M_dot_factor = 1.0;          // Star formation factor (dimensionless)
        this.tau_SF = 1e6 * 3.156e7;      // Star formation timescale (s)
        this.rho_wind = 1e-20;            // Wind density (kg/m^3)
        this.v_wind = 2e6;                // Wind velocity (m/s)
        this.rho_fluid = 1e-20;           // Fluid density (kg/m^3)
        this.P0 = 4e-8;                   // Initial pressure (Pa)
        this.tau_exp = 1e6 * 3.156e7;     // Expansion timescale (s)
        this.rho_vac_UA = 7.09e-36;       // UA vacuum density (J/m^3)
        this.rho_vac_SCm = 7.09e-37;      // SCm vacuum density (J/m^3)
        this.scale_EM = 1e-12;            // EM scaling factor
        this.proton_mass = 1.673e-27;     // Proton mass for EM acceleration

        // Additional parameters for full inclusion of terms
        this.hbar = 1.0546e-34;           // Reduced Planck's constant
        this.t_Hubble = 13.8e9 * 3.156e7; // Hubble time (s)
        this.t_Hubble_gyr = 13.8;         // Hubble time in Gyr
        this.delta_x = 1e-10;             // Position uncertainty (m)
        this.delta_p = this.hbar / this.delta_x; // Momentum uncertainty (kg m/s)
        this.integral_psi = 1.0;          // Wavefunction integral approximation
        this.A_osc = 1e-10;               // Oscillatory amplitude (m/s^2)
        this.k_osc = 1.0 / this.r;        // Wave number (1/m)
        this.omega_osc = 2 * Math.PI / (this.r / this.c_light); // Angular frequency (rad/s)
        this.x_pos = this.r;              // Position for oscillation (m)
        this.M_DM_factor = 0.1;           // Dark matter mass fraction
        this.delta_rho_over_rho = 1e-5;   // Density perturbation fraction

        this.updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    updateCache() {
        this.ug1_base = (this.G * this.M0) / (this.r * this.r);
    }

    // Universal setter for any variable (by name, for flexibility)
    setVariable(varName, newValue) {
        switch(varName) {
            case "G": this.G = newValue; break;
            case "M0": this.M0 = newValue; break;
            case "r": this.r = newValue; break;
            case "H0": this.H0 = newValue; break;
            case "B": this.B = newValue; break;
            case "B_crit": this.B_crit = newValue; break;
            case "Lambda": this.Lambda = newValue; break;
            case "c_light": this.c_light = newValue; break;
            case "q_charge": this.q_charge = newValue; break;
            case "gas_v": this.gas_v = newValue; break;
            case "f_TRZ": this.f_TRZ = newValue; break;
            case "M_dot_factor": this.M_dot_factor = newValue; break;
            case "tau_SF": this.tau_SF = newValue; break;
            case "rho_wind": this.rho_wind = newValue; break;
            case "v_wind": this.v_wind = newValue; break;
            case "rho_fluid": this.rho_fluid = newValue; break;
            case "P0": this.P0 = newValue; break;
            case "tau_exp": this.tau_exp = newValue; break;
            case "rho_vac_UA": this.rho_vac_UA = newValue; break;
            case "rho_vac_SCm": this.rho_vac_SCm = newValue; break;
            case "scale_EM": this.scale_EM = newValue; break;
            case "proton_mass": this.proton_mass = newValue; break;
            // Full terms
            case "hbar": this.hbar = newValue; break;
            case "t_Hubble": this.t_Hubble = newValue; break;
            case "t_Hubble_gyr": this.t_Hubble_gyr = newValue; break;
            case "delta_x": this.delta_x = newValue; break;
            case "delta_p": this.delta_p = newValue; break;
            case "integral_psi": this.integral_psi = newValue; break;
            case "A_osc": this.A_osc = newValue; break;
            case "k_osc": this.k_osc = newValue; break;
            case "omega_osc": this.omega_osc = newValue; break;
            case "x_pos": this.x_pos = newValue; break;
            case "M_DM_factor": this.M_DM_factor = newValue; break;
            case "delta_rho_over_rho": this.delta_rho_over_rho = newValue; break;
            default:
                console.error(`Error: Unknown variable '${varName}'.`);
                return false;
        }
        this.updateCache();
        return true;
    }

    // Addition method for variables
    addToVariable(varName, delta) {
        return this.setVariable(varName, this.getVariable(varName) + delta);
    }

    // Subtraction method for variables
    subtractFromVariable(varName, delta) {
        return this.addToVariable(varName, -delta);
    }

    // Getter for any variable (helper for add/subtract)
    getVariable(varName) {
        switch(varName) {
            case "G": return this.G;
            case "M0": return this.M0;
            case "r": return this.r;
            case "H0": return this.H0;
            case "B": return this.B;
            case "B_crit": return this.B_crit;
            case "Lambda": return this.Lambda;
            case "c_light": return this.c_light;
            case "q_charge": return this.q_charge;
            case "gas_v": return this.gas_v;
            case "f_TRZ": return this.f_TRZ;
            case "M_dot_factor": return this.M_dot_factor;
            case "tau_SF": return this.tau_SF;
            case "rho_wind": return this.rho_wind;
            case "v_wind": return this.v_wind;
            case "rho_fluid": return this.rho_fluid;
            case "P0": return this.P0;
            case "tau_exp": return this.tau_exp;
            case "rho_vac_UA": return this.rho_vac_UA;
            case "rho_vac_SCm": return this.rho_vac_SCm;
            case "scale_EM": return this.scale_EM;
            case "proton_mass": return this.proton_mass;
            // Full terms
            case "hbar": return this.hbar;
            case "t_Hubble": return this.t_Hubble;
            case "t_Hubble_gyr": return this.t_Hubble_gyr;
            case "delta_x": return this.delta_x;
            case "delta_p": return this.delta_p;
            case "integral_psi": return this.integral_psi;
            case "A_osc": return this.A_osc;
            case "k_osc": return this.k_osc;
            case "omega_osc": return this.omega_osc;
            case "x_pos": return this.x_pos;
            case "M_DM_factor": return this.M_DM_factor;
            case "delta_rho_over_rho": return this.delta_rho_over_rho;
            default:
                console.error(`Error: Unknown variable '${varName}'.`);
                return 0.0;
        }
    }

    // M(t) computation - mass growth over time
    M_t(t) {
        const M_dot = this.M_dot_factor * Math.exp(-t / this.tau_SF);
        return this.M0 * (1 + M_dot);
    }

    // P(t) computation - pressure evolution over time
    P_t(t) {
        return this.P0 * Math.exp(-t / this.tau_exp);
    }

    // Ug terms computation
    compute_Ug(Mt) {
        const Ug1 = (this.G * Mt) / (this.r * this.r);
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const corr_B = 1 - this.B / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ);
    }

    // Volume computation for fluid
    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    // Main MUGE computation (includes ALL terms)
    compute_g_NGC3603(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        const Mt = this.M_t(t);
        const Pt = this.P_t(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);

        // Term 1: Base + H0 + B corrections
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - this.B / this.B_crit;
        const term1 = ug1_t * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ
        const term2 = this.compute_Ug(Mt);

        // Term 3: Lambda
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled EM with UA
        const cross_vB = this.gas_v * this.B;  // Magnitude, assuming perpendicular
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        const term4 = (em_base * corr_UA) * this.scale_EM;

        // Quantum uncertainty term
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Fluid term (effective acceleration)
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * ug1_t) / Mt;

        // Oscillatory terms (real parts)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term (converted to acceleration)
        const M_dm = Mt * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / Mt;

        // Stellar wind feedback term (pressure / density for acceleration)
        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        const term_wind = wind_pressure / this.rho_fluid;

        // Cavity pressure term (P(t) / rho_fluid for acceleration)
        const term_pressure = Pt / this.rho_fluid;

        // Total g_NGC3603 (all terms summed)
        const total_g = term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind + term_pressure;

        // Add dynamic terms if enabled
        let dynamic_contribution = 0.0;
        if (this.enableDynamicTerms) {
            const params = new Map([
                ['rho_vac_UA', this.rho_vac_UA],
                ['hbar', this.hbar],
                ['M', Mt],
                ['r', this.r]
            ]);

            for (const term of this.dynamicTerms) {
                dynamic_contribution += term.compute(t, params);
            }
        }

        return total_g + dynamic_contribution;
    }

    // Debug/Output method (for transparency in base program)
    printParameters() {
        console.log("NGC 3603 Parameters:");
        console.log(`G: ${this.G}, M0: ${this.M0}, r: ${this.r}`);
        console.log(`H0: ${this.H0}, B: ${this.B}, B_crit: ${this.B_crit}`);
        console.log(`f_TRZ: ${this.f_TRZ}, M_dot_factor: ${this.M_dot_factor}, tau_SF: ${this.tau_SF}`);
        console.log(`rho_fluid: ${this.rho_fluid}, rho_wind: ${this.rho_wind}, v_wind: ${this.v_wind}`);
        console.log(`P0: ${this.P0}, tau_exp: ${this.tau_exp}`);
        console.log(`gas_v: ${this.gas_v}, M_DM_factor: ${this.M_DM_factor}`);
        console.log(`A_osc: ${this.A_osc}, delta_rho_over_rho: ${this.delta_rho_over_rho}`);
        console.log(`ug1_base: ${this.ug1_base}`);
    }

    // Example computation at t=500k years (for testing)
    exampleAt500kYears() {
        const t_example = 5e5 * 3.156e7;
        return this.compute_g_NGC3603(t_example);
    }

    // Get comprehensive diagnostics
    getDiagnostics(t) {
        const Mt = this.M_t(t);
        const Pt = this.P_t(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);

        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - this.B / this.B_crit;
        const term1 = ug1_t * corr_H * corr_B;

        const term2 = this.compute_Ug(Mt);
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        const cross_vB = this.gas_v * this.B;
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        const term4 = (em_base * corr_UA) * this.scale_EM;

        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * ug1_t) / Mt;

        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        const M_dm = Mt * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / Mt;

        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        const term_wind = wind_pressure / this.rho_fluid;

        const term_pressure = Pt / this.rho_fluid;

        return {
            time: t,
            mass: Mt,
            pressure: Pt,
            components: {
                base_gravity: term1,
                universal_gravity: term2,
                dark_energy: term3,
                electromagnetic: term4,
                quantum_uncertainty: term_q,
                fluid_dynamics: term_fluid,
                oscillatory_waves: term_osc,
                dark_matter: term_DM,
                stellar_wind: term_wind,
                cavity_pressure: term_pressure
            }
        };
    }
}

module.exports = { NGC3603 };