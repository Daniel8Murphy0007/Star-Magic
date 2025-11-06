// OuterFieldBubbleModule.js
// JavaScript conversion of OuterFieldBubbleModule.cpp
// Modular implementation of the Radius of the Outer Field Bubble (R_b) in the UQFF framework.
// This module computes R_b=1.496e13 m (100 AU); defines S(r - R_b) step function in Universal Gravity U_g2 term.
// Variables in Map; example for Sun at t=0; S=1 for r >= R_b, 0 otherwise.
// Approximations: S step=1 at r=R_b; δ_sw v_sw=5001; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025. Converted Nov 6, 2025.

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
        throw new Error('PhysicsTerm.compute() must be implemented by derived class');
    }

    getName() {
        throw new Error('PhysicsTerm.getName() must be implemented by derived class');
    }

    getDescription() {
        throw new Error('PhysicsTerm.getDescription() must be implemented by derived class');
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

class OuterFieldBubbleModule {
    constructor() {
        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();
        
        // Universal constants
        this.variables.set('R_b', 1.496e13);                    // m (100 AU)
        this.variables.set('AU_to_m', 1.496e11);                // m/AU
        this.variables.set('k_2', 1.2);                         // Coupling
        this.variables.set('rho_vac_UA', 7.09e-36);             // J/m^3
        this.variables.set('rho_vac_SCm', 7.09e-37);            // J/m^3
        this.variables.set('M_s', 1.989e30);                    // kg
        this.variables.set('r', 1.496e13);                      // m (default = R_b)
        this.variables.set('delta_sw', 0.01);                   // Unitless
        this.variables.set('v_sw', 5e5);                        // m/s
        this.variables.set('H_SCm', 1.0);                       // Unitless
        this.variables.set('E_react', 1e46);                    // J

        // Derived
        this.variables.set('rho_sum', this.variables.get('rho_vac_UA') + this.variables.get('rho_vac_SCm'));
        this.variables.set('swirl_factor', 1.0 + this.variables.get('delta_sw') * this.variables.get('v_sw'));

        // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;

        this.metadata.set('enhanced', 'true');
        this.metadata.set('version', '2.0-Enhanced');
    }

    // ========== DYNAMIC VARIABLE OPERATIONS ==========

    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
            if (name === 'rho_vac_UA' || name === 'rho_vac_SCm') {
                this.variables.set('rho_sum', this.variables.get('rho_vac_UA') + this.variables.get('rho_vac_SCm'));
            } else if (name === 'delta_sw' || name === 'v_sw') {
                this.variables.set('swirl_factor', 1.0 + this.variables.get('delta_sw') * this.variables.get('v_sw'));
            }
        } else {
            console.warn(`Variable '${name}' not found. Adding with value ${value}`);
            this.variables.set(name, value);
        }
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
            if (name === 'rho_vac_UA' || name === 'rho_vac_SCm') {
                this.variables.set('rho_sum', this.variables.get('rho_vac_UA') + this.variables.get('rho_vac_SCm'));
            } else if (name === 'delta_sw' || name === 'v_sw') {
                this.variables.set('swirl_factor', 1.0 + this.variables.get('delta_sw') * this.variables.get('v_sw'));
            }
        } else {
            console.warn(`Variable '${name}' not found. Adding with delta ${delta}`);
            this.variables.set(name, delta);
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    // ========== CORE COMPUTATIONS ==========

    computeR_b() {
        return this.variables.get('R_b');
    }

    computeR_bInAU() {
        return this.computeR_b() / this.variables.get('AU_to_m');
    }

    computeS_r_Rb(r) {
        return (r >= this.computeR_b()) ? 1.0 : 0.0;
    }

    computeU_g2(r) {
        this.variables.set('r', r);
        const k_2 = this.variables.get('k_2');
        const rho_sum = this.variables.get('rho_sum');
        const M_s = this.variables.get('M_s');
        const s_step = this.computeS_r_Rb(r);
        const swirl_factor = this.variables.get('swirl_factor');
        const h_scm = this.variables.get('H_SCm');
        const e_react = this.variables.get('E_react');
        return k_2 * (rho_sum * M_s / (r * r)) * s_step * swirl_factor * h_scm * e_react;
    }

    // ========== OUTPUT AND DEBUGGING ==========

    getEquationText() {
        return `U_g2 = k_2 * [(ρ_vac,[UA] + ρ_vac,[SCm]) M_s / r^2] * S(r - R_b) * (1 + δ_sw v_sw) * H_SCm * E_react
Where R_b = 1.496e13 m (100 AU, outer bubble radius);
S(r - R_b) = 1 (r >= R_b), 0 otherwise (step function).
Example r=R_b: U_g2 ≈1.18e53 J/m³; r < R_b (e.g., 1 AU): U_g2=0.
Role: Defines external gravity boundary (~heliopause); activates U_g2 beyond R_b.
UQFF: Separates internal/external fields; models heliosphere/nebular extent.`;
    }

    printVariables() {
        console.log('Current Variables:');
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential(3)}`);
        }
    }

    printU_g2Comparison(r_inside = 1.496e11, r_boundary = 1.496e13, r_outside = 1.5e13) {
        console.log('U_g2 Comparison across R_b boundary:');
        console.log(`  Inside (r = ${r_inside.toExponential(2)} m, ${(r_inside/this.variables.get('AU_to_m')).toFixed(1)} AU):`);
        console.log(`    S(r - R_b) = ${this.computeS_r_Rb(r_inside)}, U_g2 = ${this.computeU_g2(r_inside).toExponential(3)} J/m³`);
        console.log(`  Boundary (r = ${r_boundary.toExponential(2)} m, ${(r_boundary/this.variables.get('AU_to_m')).toFixed(1)} AU):`);
        console.log(`    S(r - R_b) = ${this.computeS_r_Rb(r_boundary)}, U_g2 = ${this.computeU_g2(r_boundary).toExponential(3)} J/m³`);
        console.log(`  Outside (r = ${r_outside.toExponential(2)} m, ${(r_outside/this.variables.get('AU_to_m')).toFixed(1)} AU):`);
        console.log(`    S(r - R_b) = ${this.computeS_r_Rb(r_outside)}, U_g2 = ${this.computeU_g2(r_outside).toExponential(3)} J/m³`);
    }

    // ========== SELF-EXPANDING FRAMEWORK INTERFACE ==========

    registerDynamicTerm(term) {
        if (term instanceof PhysicsTerm) {
            this.dynamicTerms.push(term);
            if (this.enableLogging) {
                console.log(`Registered dynamic term: ${term.getName()}`);
            }
            return true;
        }
        console.error('Invalid physics term');
        return false;
    }

    setDynamicParameter(name, value) {
        this.dynamicParameters.set(name, value);
        if (this.enableLogging) {
            console.log(`Set dynamic parameter ${name} = ${value}`);
        }
    }

    getDynamicParameter(name, defaultValue = 0.0) {
        return this.dynamicParameters.has(name) ? this.dynamicParameters.get(name) : defaultValue;
    }

    computeDynamicContribution(t) {
        if (!this.enableDynamicTerms || this.dynamicTerms.length === 0) {
            return 0.0;
        }

        let total = 0.0;
        const params = new Map([...this.variables, ...this.dynamicParameters]);
        
        for (const term of this.dynamicTerms) {
            if (term.validate(params)) {
                total += term.compute(t, params);
            } else if (this.enableLogging) {
                console.warn(`Term ${term.getName()} failed validation`);
            }
        }
        
        return total;
    }

    exportState() {
        const state = {
            variables: Array.from(this.variables.entries()),
            dynamicParameters: Array.from(this.dynamicParameters.entries()),
            metadata: Array.from(this.metadata.entries()),
            enableDynamicTerms: this.enableDynamicTerms,
            enableLogging: this.enableLogging,
            learningRate: this.learningRate,
            dynamicTermCount: this.dynamicTerms.length
        };
        return JSON.stringify(state, null, 2);
    }

    importState(jsonString) {
        try {
            const state = JSON.parse(jsonString);
            this.variables = new Map(state.variables);
            this.dynamicParameters = new Map(state.dynamicParameters);
            this.metadata = new Map(state.metadata);
            this.enableDynamicTerms = state.enableDynamicTerms;
            this.enableLogging = state.enableLogging;
            this.learningRate = state.learningRate;
            if (this.enableLogging) {
                console.log('State imported successfully');
            }
            return true;
        } catch (error) {
            console.error('Failed to import state:', error);
            return false;
        }
    }

    setEnableLogging(enable) {
        this.enableLogging = enable;
    }

    setLearningRate(rate) {
        this.learningRate = rate;
    }

    getMetadata(key) {
        return this.metadata.get(key);
    }

    setMetadata(key, value) {
        this.metadata.set(key, value);
    }
}

// Export for Node.js
if (typeof module !== 'undefined' && module.exports) {
    module.exports = { OuterFieldBubbleModule, PhysicsTerm, DynamicVacuumTerm, QuantumCouplingTerm };
}

// Example usage:
// const { OuterFieldBubbleModule } = require('./source110.js');
// const mod = new OuterFieldBubbleModule();
// const rb = mod.computeR_b();
// console.log('R_b =', rb, 'm (', mod.computeR_bInAU(), 'AU)');
// const u_g2 = mod.computeU_g2(1.5e13);  // r > R_b
// console.log('U_g2 (r=1.5e13 m) =', u_g2, 'J/m³');
// console.log(mod.getEquationText());
// mod.updateVariable('R_b', 2e13);
// mod.printVariables();
// Sample: R_b=1.496e13 m (100 AU); U_g2≈1.18e53 J/m³ (r>=R_b); 0 inside.
