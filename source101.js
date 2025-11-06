// HeliosphereThicknessModule.js
// JavaScript implementation of the Heliosphere Thickness Factor (H_SCm) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes H_SCm ≈1 (unitless) and its scaling in Universal Gravity U_g2 term.
// Variables in Map; example for Sun at t=0, t_n=0, r=R_b=1.496e13 m.
// Approximations: S(r - R_b)=1; δ_sw v_sw=5001; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.
// Converted to JavaScript: Nov 6, 2025

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
        throw new Error('PhysicsTerm.compute() must be implemented by subclass');
    }

    getName() {
        throw new Error('PhysicsTerm.getName() must be implemented by subclass');
    }

    getDescription() {
        throw new Error('PhysicsTerm.getDescription() must be implemented by subclass');
    }

    validate(params) {
        return true;
    }
}

class DynamicVacuumTerm extends PhysicsTerm {
    constructor(amplitude = 1e-10, frequency = 1e-15) {
        super();
        this.amplitude = amplitude;
        this.frequency = frequency;
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
    constructor(coupling_strength = 1e-40) {
        super();
        this.coupling_strength = coupling_strength;
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

class HeliosphereThicknessModule {
    constructor() {
        // Core variables using Map
        this.variables = new Map();
        
        // Self-expanding framework members
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;
        
        // Set metadata
        this.metadata.set('enhanced', 'true');
        this.metadata.set('version', '2.0-Enhanced');
        
        // Initialize universal constants
        this.variables.set('H_SCm', 1.0);              // Unitless ≈1
        this.variables.set('k_2', 1.2);                // Coupling
        this.variables.set('rho_vac_UA', 7.09e-36);    // J/m³
        this.variables.set('rho_vac_SCm', 7.09e-37);   // J/m³
        this.variables.set('M_s', 1.989e30);           // kg (Sun)
        this.variables.set('r', 1.496e13);             // m (R_b)
        this.variables.set('R_b', 1.496e13);           // m
        this.variables.set('delta_sw', 0.01);          // Unitless
        this.variables.set('v_sw', 5e5);               // m/s
        this.variables.set('E_react', 1e46);           // J
        this.variables.set('S_r_Rb', 1.0);             // Step function
        this.variables.set('pi', Math.PI);
        this.variables.set('t_n', 0.0);                // s
        
        // Derived variables
        this.variables.set('rho_sum', this.variables.get('rho_vac_UA') + this.variables.get('rho_vac_SCm'));
        this.variables.set('swirl_factor', 1.0 + this.variables.get('delta_sw') * this.variables.get('v_sw'));
    }

    // ========== DYNAMIC VARIABLE OPERATIONS ==========
    
    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
            this._updateDerivedVariables(name);
        } else {
            console.warn(`Variable '${name}' not found. Adding with value ${value}`);
            this.variables.set(name, value);
        }
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
            this._updateDerivedVariables(name);
        } else {
            console.warn(`Variable '${name}' not found. Adding with delta ${delta}`);
            this.variables.set(name, delta);
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    _updateDerivedVariables(name) {
        // Update derived variables when dependencies change
        if (name === 'rho_vac_UA' || name === 'rho_vac_SCm') {
            this.variables.set('rho_sum', 
                this.variables.get('rho_vac_UA') + this.variables.get('rho_vac_SCm'));
        } else if (name === 'delta_sw' || name === 'v_sw') {
            this.variables.set('swirl_factor', 
                1.0 + this.variables.get('delta_sw') * this.variables.get('v_sw'));
        }
    }

    // ========== CORE COMPUTATIONS ==========
    
    /**
     * Compute H_SCm ≈1 (unitless heliosphere thickness factor)
     * @returns {number} H_SCm value
     */
    computeH_SCm() {
        return this.variables.get('H_SCm');
    }

    /**
     * Compute U_g2 with H_SCm scaling
     * U_g2 = k_2 * [(ρ_vac,UA + ρ_vac,SCm) M_s / r²] * S(r - R_b) * (1 + δ_sw v_sw) * H_SCm * E_react
     * @param {number} t - Time (s)
     * @param {number} t_n - Reference time (s)
     * @returns {number} U_g2 in J/m³
     */
    computeU_g2(t, t_n) {
        const k_2 = this.variables.get('k_2');
        const rho_sum = this.variables.get('rho_sum');
        const M_s = this.variables.get('M_s');
        const r = this.variables.get('r');
        const S_r_Rb = this.variables.get('S_r_Rb');
        const swirl_factor = this.variables.get('swirl_factor');
        const H_SCm = this.computeH_SCm();
        const E_react = this.variables.get('E_react');
        
        // Apply dynamic terms if enabled
        let dynamicContribution = 0.0;
        if (this.enableDynamicTerms && this.dynamicTerms.length > 0) {
            for (const term of this.dynamicTerms) {
                dynamicContribution += term.compute(t, this.variables);
            }
        }
        
        // Core U_g2 calculation
        const baseU_g2 = k_2 * (rho_sum * M_s / (r * r)) * S_r_Rb * swirl_factor * H_SCm * E_react;
        
        return baseU_g2 + dynamicContribution;
    }

    /**
     * Compute U_g2 without H_SCm variation (H=1 fixed)
     * @param {number} t - Time (s)
     * @param {number} t_n - Reference time (s)
     * @returns {number} U_g2 with H_SCm=1 in J/m³
     */
    computeU_g2_no_H(t, t_n) {
        const orig_H = this.variables.get('H_SCm');
        this.variables.set('H_SCm', 1.0);
        const result = this.computeU_g2(t, t_n);
        this.variables.set('H_SCm', orig_H);
        return result;
    }

    // ========== SELF-EXPANDING FRAMEWORK METHODS ==========
    
    /**
     * Register a dynamic physics term
     * @param {PhysicsTerm} term - Physics term to add
     */
    registerDynamicTerm(term) {
        if (term instanceof PhysicsTerm) {
            this.dynamicTerms.push(term);
            if (this.enableLogging) {
                console.log(`Registered dynamic term: ${term.getName()}`);
            }
        } else {
            throw new Error('Term must be an instance of PhysicsTerm');
        }
    }

    /**
     * Set a dynamic parameter
     * @param {string} name - Parameter name
     * @param {number} value - Parameter value
     */
    setDynamicParameter(name, value) {
        this.dynamicParameters.set(name, value);
        if (this.enableLogging) {
            console.log(`Set dynamic parameter: ${name} = ${value}`);
        }
    }

    /**
     * Get a dynamic parameter
     * @param {string} name - Parameter name
     * @returns {number|undefined} Parameter value
     */
    getDynamicParameter(name) {
        return this.dynamicParameters.get(name);
    }

    /**
     * Export state for persistence or cross-module communication
     * @returns {object} State object
     */
    exportState() {
        return {
            variables: Object.fromEntries(this.variables),
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            metadata: Object.fromEntries(this.metadata),
            enableDynamicTerms: this.enableDynamicTerms,
            enableLogging: this.enableLogging,
            learningRate: this.learningRate,
            dynamicTermCount: this.dynamicTerms.length
        };
    }

    /**
     * Import state from exported data
     * @param {object} state - State object
     */
    importState(state) {
        if (state.variables) {
            this.variables = new Map(Object.entries(state.variables));
        }
        if (state.dynamicParameters) {
            this.dynamicParameters = new Map(Object.entries(state.dynamicParameters));
        }
        if (state.metadata) {
            this.metadata = new Map(Object.entries(state.metadata));
        }
        if (typeof state.enableDynamicTerms !== 'undefined') {
            this.enableDynamicTerms = state.enableDynamicTerms;
        }
        if (typeof state.enableLogging !== 'undefined') {
            this.enableLogging = state.enableLogging;
        }
        if (typeof state.learningRate !== 'undefined') {
            this.learningRate = state.learningRate;
        }
    }

    /**
     * Set learning rate for auto-optimization
     * @param {number} rate - Learning rate
     */
    setLearningRate(rate) {
        this.learningRate = rate;
    }

    /**
     * Enable or disable logging
     * @param {boolean} enabled - Logging state
     */
    setEnableLogging(enabled) {
        this.enableLogging = enabled;
    }

    // ========== OUTPUT METHODS ==========
    
    /**
     * Get equation description text
     * @returns {string} Equation description
     */
    getEquationText() {
        return `U_g2 = k_2 * [(ρ_vac,UA + ρ_vac,SCm) M_s / r²] * S(r - R_b) * (1 + δ_sw v_sw) * H_SCm * E_react
Where H_SCm ≈1 (unitless heliosphere thickness factor);
Scales outer field bubble gravity for heliopause extent (~120 AU).
Example r=R_b=1.496e13 m, t=0: U_g2 ≈1.18e53 J/m³ (H=1);
If H_SCm=1.1: ≈1.30e53 J/m³ (+10%).
Role: Adjusts [SCm] influence in heliosphere; minimal but flexible for boundary variations.
UQFF: Models solar wind dominance; key for nebular/heliospheric dynamics.`;
    }

    /**
     * Print all current variables
     */
    printVariables() {
        console.log('\n=== HeliosphereThicknessModule Variables ===');
        for (const [key, value] of this.variables) {
            console.log(`  ${key} = ${value.toExponential(4)}`);
        }
    }

    /**
     * Print U_g2 comparison with and without H_SCm variation
     * @param {number} t - Time (s)
     * @param {number} t_n - Reference time (s)
     */
    printU_g2Comparison(t = 0.0, t_n = 0.0) {
        const u_g2_with_H = this.computeU_g2(t, t_n);
        const u_g2_no_H = this.computeU_g2_no_H(t, t_n);
        const H_SCm = this.computeH_SCm();
        
        console.log('\n=== U_g2 Comparison ===');
        console.log(`  H_SCm = ${H_SCm.toFixed(4)}`);
        console.log(`  U_g2 (with H_SCm) = ${u_g2_with_H.toExponential(4)} J/m³`);
        console.log(`  U_g2 (H=1) = ${u_g2_no_H.toExponential(4)} J/m³`);
        console.log(`  Difference = ${((u_g2_with_H - u_g2_no_H) / u_g2_no_H * 100).toFixed(2)}%`);
    }

    /**
     * Print component breakdown
     * @param {number} t - Time (s)
     * @param {number} t_n - Reference time (s)
     */
    printComponentBreakdown(t = 0.0, t_n = 0.0) {
        const k_2 = this.variables.get('k_2');
        const rho_sum = this.variables.get('rho_sum');
        const M_s = this.variables.get('M_s');
        const r = this.variables.get('r');
        const S_r_Rb = this.variables.get('S_r_Rb');
        const swirl_factor = this.variables.get('swirl_factor');
        const H_SCm = this.computeH_SCm();
        const E_react = this.variables.get('E_react');
        
        const base_term = rho_sum * M_s / (r * r);
        const scaled_term = base_term * S_r_Rb * swirl_factor;
        const with_H = scaled_term * H_SCm;
        const final_U_g2 = k_2 * with_H * E_react;
        
        console.log('\n=== U_g2 Component Breakdown ===');
        console.log(`  k_2 = ${k_2.toFixed(2)}`);
        console.log(`  ρ_sum = ${rho_sum.toExponential(4)} J/m³`);
        console.log(`  M_s = ${M_s.toExponential(4)} kg`);
        console.log(`  r = ${r.toExponential(4)} m`);
        console.log(`  Base term (ρ_sum * M_s / r²) = ${base_term.toExponential(4)} J/m`);
        console.log(`  S(r-R_b) = ${S_r_Rb.toFixed(1)}`);
        console.log(`  Swirl factor (1 + δ_sw v_sw) = ${swirl_factor.toExponential(4)}`);
        console.log(`  Scaled term = ${scaled_term.toExponential(4)} J/m`);
        console.log(`  H_SCm = ${H_SCm.toFixed(4)}`);
        console.log(`  With H_SCm = ${with_H.toExponential(4)} J/m`);
        console.log(`  E_react = ${E_react.toExponential(4)} J`);
        console.log(`  Final U_g2 = ${final_U_g2.toExponential(4)} J/m³`);
    }

    /**
     * Print module information
     */
    printModuleInfo() {
        console.log('\n=== HeliosphereThicknessModule Info ===');
        console.log(`  Module Version:     ${this.metadata.get('version')}`);
        console.log(`  Enhanced:           ${this.metadata.get('enhanced')}`);
        console.log(`  Variables:          ${this.variables.size}`);
        console.log(`  Dynamic Terms:      ${this.dynamicTerms.length}`);
        console.log(`  Dynamic Parameters: ${this.dynamicParameters.size}`);
        console.log(`  Logging Enabled:    ${this.enableLogging}`);
        console.log(`  Learning Rate:      ${this.learningRate}`);
    }
}

// Export classes
module.exports = {
    HeliosphereThicknessModule,
    PhysicsTerm,
    DynamicVacuumTerm,
    QuantumCouplingTerm
};

// Example usage:
// const { HeliosphereThicknessModule } = require('./source101.js');
// const mod = new HeliosphereThicknessModule();
// const h = mod.computeH_SCm();
// console.log(`H_SCm ≈ ${h}`);
// const u_g2 = mod.computeU_g2(0.0, 0.0);
// console.log(`U_g2 = ${u_g2} J/m³`);
// console.log(mod.getEquationText());
// mod.updateVariable('H_SCm', 1.1);
// mod.printVariables();
// Sample: H_SCm=1; U_g2 ≈1.18e53 J/m³; +10% for H=1.1.
