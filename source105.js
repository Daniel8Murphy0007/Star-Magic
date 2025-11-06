// GalacticBlackHoleModule.js
// JavaScript implementation of the Mass of the Galactic Black Hole (M_bh) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes M_bh=8.15e36 kg ≈4.1e6 M_sun; scales M_bh / d_g in Universal Buoyancy U_bi and Ug4.
// Converted from source105.cpp
// Variables in Map; example for Sun at t=0, t_n=0.
// Approximations: cos(π t_n)=1; (1 + ε_sw ρ_vac,sw)≈1; α=0.001 s^-1; f_feedback=0.1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

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
        throw new Error('compute() must be implemented by subclass');
    }

    getName() {
        throw new Error('getName() must be implemented by subclass');
    }

    getDescription() {
        throw new Error('getDescription() must be implemented by subclass');
    }

    validate(params) {
        return true;
    }
}

class DynamicVacuumTerm extends PhysicsTerm {
    constructor(amp = 1e-10, freq = 1e-15) {
        super();
        this.amplitude = amp;
        this.frequency = freq;
    }

    compute(t, params) {
        const rho_vac = params.has('rho_vac_UA') ? params.get('rho_vac_UA') : 7.09e-36;
        return this.amplitude * rho_vac * Math.sin(this.frequency * t);
    }

    getName() {
        return 'DynamicVacuum';
    }

    getDescription() {
        return 'Time-varying vacuum energy';
    }
}

class QuantumCouplingTerm extends PhysicsTerm {
    constructor(strength = 1e-40) {
        super();
        this.coupling_strength = strength;
    }

    compute(t, params) {
        const hbar = params.has('hbar') ? params.get('hbar') : 1.0546e-34;
        const M = params.has('M') ? params.get('M') : 1.989e30;
        const r = params.has('r') ? params.get('r') : 1e4;
        return this.coupling_strength * (hbar * hbar) / (M * r * r) * Math.cos(t / 1e6);
    }

    getName() {
        return 'QuantumCoupling';
    }

    getDescription() {
        return 'Non-local quantum effects';
    }
}

// ===========================================================================================
// ENHANCED CLASS WITH SELF-EXPANDING CAPABILITIES
// ===========================================================================================

class GalacticBlackHoleModule {
    constructor() {
        // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;

        this.metadata.set('enhanced', 'true');
        this.metadata.set('version', '2.0-Enhanced');

        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();

        // Universal constants
        this.variables.set('M_sun', 1.989e30);              // kg
        this.variables.set('M_bh', 8.15e36);                // kg (Sgr A*)

        // Shared params for terms
        this.variables.set('beta_1', 0.6);                  // Unitless
        this.variables.set('U_g1', 1.39e26);                // J/m^3
        this.variables.set('Omega_g', 7.3e-16);             // rad/s
        this.variables.set('d_g', 2.55e20);                 // m
        this.variables.set('epsilon_sw', 0.001);            // Unitless
        this.variables.set('rho_vac_sw', 8e-21);            // J/m^3
        this.variables.set('U_UA', 1.0);                    // Normalized
        this.variables.set('t_n', 0.0);                     // s
        this.variables.set('pi', Math.PI);

        // Ug4 params
        this.variables.set('k_4', 1.0);                     // Unitless
        this.variables.set('rho_vac_SCm', 7.09e-37);        // J/m^3
        this.variables.set('alpha', 0.001);                 // s^-1
        this.variables.set('f_feedback', 0.1);              // Unitless
    }

    // ========== DYNAMIC VARIABLE OPERATIONS ==========

    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
        } else {
            if (this.enableLogging) {
                console.log(`Variable '${name}' not found. Adding with value ${value}`);
            }
            this.variables.set(name, value);
        }
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
        } else {
            if (this.enableLogging) {
                console.log(`Variable '${name}' not found. Adding with delta ${delta}`);
            }
            this.variables.set(name, delta);
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    // ========== CORE COMPUTATIONS ==========

    /**
     * Compute M_bh (kg)
     * @returns {number} Black hole mass in kg
     */
    computeM_bh() {
        return this.variables.get('M_bh');
    }

    /**
     * M_bh in M_sun
     * @returns {number} Black hole mass in solar masses
     */
    computeM_bhInMsun() {
        return this.computeM_bh() / this.variables.get('M_sun');
    }

    /**
     * M_bh / d_g (kg/m)
     * @returns {number} Mass per unit distance
     */
    computeMbhOverDg() {
        return this.computeM_bh() / this.variables.get('d_g');
    }

    /**
     * U_b1 example (J/m^3)
     * Universal Buoyancy with SMBH scaling
     * @returns {number} U_b1 in J/m^3
     */
    computeU_b1() {
        const beta_1 = this.variables.get('beta_1');
        const U_g1 = this.variables.get('U_g1');
        const Omega_g = this.variables.get('Omega_g');
        const mbh_over_dg = this.computeMbhOverDg();
        const swirl_factor = 1.0 + this.variables.get('epsilon_sw') * this.variables.get('rho_vac_sw');
        const U_UA = this.variables.get('U_UA');
        const cos_term = Math.cos(this.variables.get('pi') * this.variables.get('t_n'));
        return -beta_1 * U_g1 * Omega_g * mbh_over_dg * swirl_factor * U_UA * cos_term;
    }

    /**
     * U_g4 example (J/m^3)
     * Star-BH interaction component with SMBH
     * @returns {number} U_g4 in J/m^3
     */
    computeU_g4() {
        const k_4 = this.variables.get('k_4');
        const rho_vac_SCm = this.variables.get('rho_vac_SCm');
        const mbh_over_dg = this.computeMbhOverDg();
        const exp_term = Math.exp(-this.variables.get('alpha') * this.variables.get('t_n'));
        const cos_term = Math.cos(this.variables.get('pi') * this.variables.get('t_n'));
        const feedback_factor = 1.0 + this.variables.get('f_feedback');
        return k_4 * (rho_vac_SCm * this.computeM_bh() / this.variables.get('d_g')) * exp_term * cos_term * feedback_factor;
    }

    // ========== SELF-EXPANDING FRAMEWORK METHODS ==========

    registerDynamicTerm(term) {
        if (term instanceof PhysicsTerm) {
            this.dynamicTerms.push(term);
            if (this.enableLogging) {
                console.log(`Registered dynamic term: ${term.getName()}`);
            }
        } else {
            throw new Error('Term must be instance of PhysicsTerm');
        }
    }

    setDynamicParameter(name, value) {
        this.dynamicParameters.set(name, value);
        if (this.enableLogging) {
            console.log(`Set dynamic parameter ${name} = ${value}`);
        }
    }

    getDynamicParameter(name) {
        return this.dynamicParameters.get(name);
    }

    exportState() {
        const state = {
            variables: Object.fromEntries(this.variables),
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            metadata: Object.fromEntries(this.metadata),
            enableDynamicTerms: this.enableDynamicTerms,
            enableLogging: this.enableLogging,
            learningRate: this.learningRate
        };
        return JSON.stringify(state, null, 2);
    }

    importState(stateJson) {
        const state = JSON.parse(stateJson);
        this.variables = new Map(Object.entries(state.variables));
        this.dynamicParameters = new Map(Object.entries(state.dynamicParameters));
        this.metadata = new Map(Object.entries(state.metadata));
        this.enableDynamicTerms = state.enableDynamicTerms;
        this.enableLogging = state.enableLogging;
        this.learningRate = state.learningRate;
    }

    setEnableLogging(enabled) {
        this.enableLogging = enabled;
    }

    setLearningRate(rate) {
        this.learningRate = rate;
    }

    // ========== OUTPUT AND DISPLAY METHODS ==========

    getEquationText() {
        return `U_bi = -β_i U_gi Ω_g (M_bh / d_g) (1 + ε_sw ρ_vac,sw) U_UA cos(π t_n)
U_g4 = k_4 (ρ_vac,[SCm] M_bh / d_g) e^{-α t} cos(π t_n) (1 + f_feedback)
Where M_bh = 8.15e36 kg ≈4.1e6 M_sun (Sgr A*).
M_bh / d_g ≈3.20e16 kg/m;
Example U_b1 ≈ -1.94e27 J/m³; U_g4 ≈2.50e-20 J/m³ (t_n=0).
Role: Scales SMBH gravity in buoyancy/Ug4; drives galactic dynamics/mergers.
UQFF: Central mass for star formation/nebulae; resolves parsec problem.`;
    }

    printVariables() {
        console.log('\nCurrent Variables:');
        for (const [name, value] of this.variables) {
            console.log(`  ${name} = ${value.toExponential(4)}`);
        }
    }

    printBlackHoleProperties() {
        console.log('\n=== Galactic Black Hole Properties ===');
        console.log(`  M_bh = ${this.computeM_bh().toExponential(4)} kg`);
        console.log(`  M_bh = ${this.computeM_bhInMsun().toExponential(4)} M_sun`);
        console.log(`  M_bh / d_g = ${this.computeMbhOverDg().toExponential(4)} kg/m`);
        console.log(`  U_b1 = ${this.computeU_b1().toExponential(4)} J/m³`);
        console.log(`  U_g4 = ${this.computeU_g4().toExponential(4)} J/m³`);
    }

    printTimeEvolution(times) {
        console.log('\n=== Time Evolution of SMBH Contributions ===');
        for (const t of times) {
            this.updateVariable('t_n', t);
            console.log(`\nTime t_n = ${t.toExponential(2)} s:`);
            console.log(`  U_b1 = ${this.computeU_b1().toExponential(4)} J/m³`);
            console.log(`  U_g4 = ${this.computeU_g4().toExponential(4)} J/m³`);
            console.log(`  exp(-α t_n) = ${Math.exp(-this.variables.get('alpha') * t).toExponential(4)}`);
            console.log(`  cos(π t_n) = ${Math.cos(Math.PI * t).toFixed(6)}`);
        }
        // Reset to t_n = 0
        this.updateVariable('t_n', 0.0);
    }

    printModuleInfo() {
        console.log('\n=== GalacticBlackHoleModule Info ===');
        console.log(`  Module Version:     ${this.metadata.get('version')}`);
        console.log(`  Enhanced:           ${this.metadata.get('enhanced')}`);
        console.log(`  Variables:          ${this.variables.size}`);
        console.log(`  M_bh:               ${this.variables.get('M_bh').toExponential(4)} kg`);
        console.log(`  d_g:                ${this.variables.get('d_g').toExponential(4)} m`);
        console.log(`  α:                  ${this.variables.get('alpha').toExponential(4)} s^-1`);
        console.log(`  Dynamic Terms:      ${this.dynamicTerms.length}`);
        console.log(`  Dynamic Parameters: ${this.dynamicParameters.size}`);
        console.log(`  Logging Enabled:    ${this.enableLogging}`);
        console.log(`  Learning Rate:      ${this.learningRate}`);
    }
}

// ===========================================================================================
// EXPORTS
// ===========================================================================================

module.exports = {
    PhysicsTerm,
    DynamicVacuumTerm,
    QuantumCouplingTerm,
    GalacticBlackHoleModule
};

// Example usage:
// const { GalacticBlackHoleModule } = require('./source105.js');
// const mod = new GalacticBlackHoleModule();
// console.log(`M_bh ≈ ${mod.computeM_bhInMsun().toExponential(2)} M_sun`);
// console.log(`U_b1 = ${mod.computeU_b1().toExponential(4)} J/m³`);
// console.log(mod.getEquationText());
// mod.updateVariable('M_bh', 9e36);
// mod.printVariables();
