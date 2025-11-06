// HeavisideFractionModule.js
// JavaScript implementation of the Heaviside Component Fraction (f_Heaviside) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes f_Heaviside=0.01 (unitless) and its scaling (1 + 10^13 * f_Heaviside) in Universal Magnetism U_m term.
// Converted from source100.cpp for Star-Magic UQFF Framework
// Variables in Map; example for Sun at t=0, t_n=0; amplifies by ~10^11.
// Approximations: 1 - e^{-γ t cos(π t_n)}=0 at t=0; φ_hat_j=1; P_SCm=1; f_quasi=0.01.
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

class HeavisideFractionModule {
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
        this.metadata.set('module', 'HeavisideFractionModule');
        this.metadata.set('converted', 'source100.cpp->source100.js');

        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();
        
        // Universal constants
        this.variables.set('f_Heaviside', 0.01);              // Unitless fraction
        this.variables.set('scale_Heaviside', 1e13);          // Amplification factor
        this.variables.set('f_quasi', 0.01);                  // Quasi factor
        this.variables.set('mu_j', 3.38e23);                  // T m^3 (j=1)
        this.variables.set('r_j', 1.496e13);                  // m
        this.variables.set('gamma', 5e-5 / 86400.0);          // day^-1 to s^-1
        this.variables.set('t_n', 0.0);                       // s
        this.variables.set('phi_hat_j', 1.0);                 // Normalized
        this.variables.set('P_SCm', 1.0);                     // Pressure
        this.variables.set('E_react', 1e46);                  // J
        this.variables.set('pi', Math.PI);

        // Derived
        this.variables.set('heaviside_factor', this.computeHeavisideFactor());
    }

    // ========== DYNAMIC VARIABLE OPERATIONS ==========

    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
            if (name === 'f_Heaviside') {
                this.variables.set('heaviside_factor', this.computeHeavisideFactor());
            }
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
            if (name === 'f_Heaviside') {
                this.variables.set('heaviside_factor', this.computeHeavisideFactor());
            }
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

    // Compute f_Heaviside (0.01)
    computeF_Heaviside() {
        return this.variables.get('f_Heaviside');
    }

    // Compute 1 + 10^13 * f_Heaviside
    computeHeavisideFactor() {
        return 1.0 + this.variables.get('scale_Heaviside') * this.computeF_Heaviside();
    }

    // Base for U_m without Heaviside/Quasi
    computeUmBase(j, t) {
        const mu_over_rj = this.variables.get('mu_j') / this.variables.get('r_j');
        const exp_arg = -this.variables.get('gamma') * t * Math.cos(this.variables.get('pi') * this.variables.get('t_n'));
        const one_minus_exp = 1.0 - Math.exp(exp_arg);
        const phi_hat = this.variables.get('phi_hat_j');
        return mu_over_rj * one_minus_exp * phi_hat * this.variables.get('P_SCm') * this.variables.get('E_react');
    }

    // U_m contribution with Heaviside
    computeUmContribution(j = 1, t = 0.0) {
        const base = this.computeUmBase(j, t);
        const heaviside_f = this.computeHeavisideFactor();
        const quasi_f = 1.0 + this.variables.get('f_quasi');
        
        // Apply dynamic terms if enabled
        let dynamicContribution = 0;
        if (this.enableDynamicTerms && this.dynamicTerms.length > 0) {
            const params = this.variables;
            for (const term of this.dynamicTerms) {
                dynamicContribution += term.compute(t, params);
            }
        }
        
        return base * heaviside_f * quasi_f + dynamicContribution;
    }

    // U_m without Heaviside (set f=0 temporarily)
    computeUmWithNoHeaviside(j = 1, t = 0.0) {
        const orig_f = this.variables.get('f_Heaviside');
        this.variables.set('f_Heaviside', 0.0);
        const result = this.computeUmContribution(j, t);
        this.variables.set('f_Heaviside', orig_f);
        return result;
    }

    // ========== SELF-EXPANDING FRAMEWORK METHODS ==========

    registerDynamicTerm(term) {
        if (term instanceof PhysicsTerm) {
            this.dynamicTerms.push(term);
            if (this.enableLogging) {
                console.log(`Registered dynamic term: ${term.getName()}`);
            }
            return true;
        }
        return false;
    }

    setDynamicParameter(name, value) {
        this.dynamicParameters.set(name, value);
        if (this.enableLogging) {
            console.log(`Set dynamic parameter: ${name} = ${value}`);
        }
    }

    getDynamicParameter(name) {
        return this.dynamicParameters.get(name);
    }

    setEnableLogging(enable) {
        this.enableLogging = enable;
    }

    setLearningRate(rate) {
        this.learningRate = rate;
    }

    // ========== STATE MANAGEMENT ==========

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
        return state;
    }

    importState(state) {
        if (state.variables) {
            this.variables = new Map(state.variables);
        }
        if (state.dynamicParameters) {
            this.dynamicParameters = new Map(state.dynamicParameters);
        }
        if (state.metadata) {
            this.metadata = new Map(state.metadata);
        }
        if (state.enableDynamicTerms !== undefined) {
            this.enableDynamicTerms = state.enableDynamicTerms;
        }
        if (state.enableLogging !== undefined) {
            this.enableLogging = state.enableLogging;
        }
        if (state.learningRate !== undefined) {
            this.learningRate = state.learningRate;
        }
    }

    // ========== OUTPUT AND DISPLAY ==========

    getEquationText() {
        return `U_m = ∑_j [ (μ_j / r_j) (1 - e^{-γ t cos(π t_n)}) φ_hat_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)
Where f_Heaviside = ${this.variables.get('f_Heaviside')} (unitless Heaviside fraction);
Heaviside factor = 1 + 10^13 * ${this.variables.get('f_Heaviside')} = ${this.computeHeavisideFactor().toExponential(2)} (amplifies ~10^11x).
Example j=1, t=0: U_m contrib ≈2.28e65 J/m³ (with); ≈2.28e54 J/m³ (without).
Role: Threshold-activated scaling in magnetic energy; nonlinear [SCm]/[UA] effects.
UQFF: Amplifies small fraction for large impact in nebulae/quasars/jets.`;
    }

    printVariables() {
        console.log('\nCurrent Variables:');
        for (const [key, value] of this.variables) {
            console.log(`  ${key} = ${value.toExponential(4)}`);
        }
    }

    printUmComparison(j = 1, t = 0.0) {
        const um_with = this.computeUmContribution(j, t);
        const um_without = this.computeUmWithNoHeaviside(j, t);
        const amplification = um_with / um_without;
        
        console.log(`\nU_m Comparison for j=${j} at t=${t} s:`);
        console.log(`  With Heaviside:    ${um_with.toExponential(4)} J/m³`);
        console.log(`  Without Heaviside: ${um_without.toExponential(4)} J/m³`);
        console.log(`  Amplification:     ~${amplification.toExponential(2)}x`);
    }

    printComponentBreakdown(j = 1, t = 0.0) {
        const base = this.computeUmBase(j, t);
        const heaviside_f = this.computeHeavisideFactor();
        const quasi_f = 1.0 + this.variables.get('f_quasi');
        const um_contrib = this.computeUmContribution(j, t);

        console.log(`\n=== Heaviside Fraction Module Breakdown (j=${j}, t=${t}s) ===`);
        console.log(`  Base U_m:           ${base.toExponential(4)} J/m³`);
        console.log(`  f_Heaviside:        ${this.computeF_Heaviside()}`);
        console.log(`  Heaviside Factor:   ${heaviside_f.toExponential(4)}`);
        console.log(`  Quasi Factor:       ${quasi_f.toFixed(2)}`);
        console.log(`  Total U_m:          ${um_contrib.toExponential(4)} J/m³`);
        
        if (this.dynamicTerms.length > 0) {
            console.log(`  Dynamic Terms:      ${this.dynamicTerms.length} active`);
        }
    }

    printModuleInfo() {
        console.log('\n=== HeavisideFractionModule Info ===');
        console.log(`  Module Version:     ${this.metadata.get('version')}`);
        console.log(`  Enhanced:           ${this.metadata.get('enhanced')}`);
        console.log(`  Variables:          ${this.variables.size}`);
        console.log(`  Dynamic Terms:      ${this.dynamicTerms.length}`);
        console.log(`  Dynamic Parameters: ${this.dynamicParameters.size}`);
        console.log(`  Logging Enabled:    ${this.enableLogging}`);
        console.log(`  Learning Rate:      ${this.learningRate}`);
    }
}

// ========== MODULE EXPORTS ==========
module.exports = {
    HeavisideFractionModule,
    PhysicsTerm,
    DynamicVacuumTerm,
    QuantumCouplingTerm
};
