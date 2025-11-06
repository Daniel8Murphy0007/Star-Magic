// NegativeTimeModule.js
// JavaScript implementation of the Negative Time Factor (t_n) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes t_n = t - t_0 (s or days, allows t_n < 0); used in cos(π t_n) for oscillations and exp(-γ t cos(π t_n)) for growth/decay.
// Converted from source106.cpp
// Variables in Map; defaults t_0=0, t=0 (t_n=0); example for U_m term with t_n negative.
// Approximations: cos even function; γ=5e-5 day^-1; at t_n=-1, exp term negative (growth).
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

class NegativeTimeModule {
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
        this.variables.set('t_0', 0.0);                     // Reference time (s/days)
        this.variables.set('t', 0.0);                       // Current time
        this.variables.set('gamma', 5e-5);                  // day^-1 (example)
        this.variables.set('pi', Math.PI);
        this.variables.set('mu_over_rj', 2.26e10);          // T m^2 (example)
        this.variables.set('P_SCm', 1.0);                   // Normalized
        this.variables.set('E_react', 1e46);                // J
        this.variables.set('heaviside_f', 1e11 + 1.0);      // 1 + 10^13 * 0.01
        this.variables.set('quasi_f', 1.01);                // 1 + 0.01
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
     * Compute t_n = t - t_0 (s/days)
     * @param {number} t - Current time
     * @returns {number} Normalized time t_n (can be negative)
     */
    computeT_n(t) {
        this.variables.set('t', t);
        return t - this.variables.get('t_0');
    }

    /**
     * Compute cos(π t_n)
     * @param {number} t - Current time
     * @returns {number} Cosine of π*t_n
     */
    computeCosPiTn(t) {
        const t_n = this.computeT_n(t);
        return Math.cos(this.variables.get('pi') * t_n);
    }

    /**
     * Compute exp(-γ t cos(π t_n))
     * @param {number} gamma - Decay/growth rate
     * @param {number} t - Current time
     * @returns {number} Exponential term
     */
    computeExpTerm(gamma, t) {
        const cos_pi_tn = this.computeCosPiTn(t);
        const arg = -gamma * t * cos_pi_tn;
        return Math.exp(arg);
    }

    /**
     * Compute 1 - exp(-γ t cos(π t_n))
     * Used in U_m and other field contributions
     * @param {number} gamma - Decay/growth rate
     * @param {number} t - Current time
     * @returns {number} One minus exponential term
     */
    computeOneMinusExp(gamma, t) {
        return 1.0 - this.computeExpTerm(gamma, t);
    }

    /**
     * Simplified U_m example contribution
     * @param {number} t - Current time
     * @param {number} mu_over_rj - Magnetic moment over distance ratio
     * @returns {number} U_m contribution in J/m^3
     */
    computeUmExample(t, mu_over_rj = null) {
        if (mu_over_rj === null) {
            mu_over_rj = this.variables.get('mu_over_rj');
        }
        const gamma = this.variables.get('gamma');
        const one_minus_exp = this.computeOneMinusExp(gamma, t);
        const phi_hat = 1.0;
        const p_scm = this.variables.get('P_SCm');
        const e_react = this.variables.get('E_react');
        const heaviside_f = this.variables.get('heaviside_f');
        const quasi_f = this.variables.get('quasi_f');
        return (mu_over_rj * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
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
        return `t_n = t - t_0 (s/days, allows t_n < 0 for time-reversal);
Used in: cos(π t_n) for oscillations; exp(-γ t cos(π t_n)) for decay/growth.
In U_m: ... (1 - exp(-γ t cos(π t_n))) ...;
Negative t_n: e.g., t_n=-1 ⇒ cos(-π)=-1 ⇒ exp(γ t) >1 (growth, negentropic).
Example t=1000 days, γ=5e-5 day^-1, t_0=0: 1-exp ≈0.049, U_m ≈1.12e66 J/m³.
t_n=-1000: same (cos even); t_n=-1: 1-exp ≈ -0.051 (growth phase).
Role: Models cyclic/TRZ dynamics; forward/reverse time in nebulae/mergers/jets.`;
    }

    printVariables() {
        console.log('\nCurrent Variables:');
        for (const [name, value] of this.variables) {
            console.log(`  ${name} = ${value.toExponential(4)}`);
        }
    }

    printTnEffects(t, gamma = null) {
        if (gamma === null) {
            gamma = this.variables.get('gamma');
        }

        // Positive t_n example
        const t_n_pos = this.computeT_n(t);
        const cos_pos = this.computeCosPiTn(t);
        const exp_pos = this.computeExpTerm(gamma, t);
        const one_minus_pos = this.computeOneMinusExp(gamma, t);
        const um_pos = this.computeUmExample(t);

        // Negative t_n: adjust t_0 to make t_n negative
        const orig_t0 = this.variables.get('t_0');
        this.variables.set('t_0', t + 1.0);  // t_n = t - (t+1) = -1
        const t_n_neg = this.computeT_n(t);
        const cos_neg = this.computeCosPiTn(t);
        const exp_neg = this.computeExpTerm(gamma, t);
        const one_minus_neg = this.computeOneMinusExp(gamma, t);
        const um_neg = this.computeUmExample(t);

        this.variables.set('t_0', orig_t0);  // Restore

        console.log(`\nt_n Effects at t=${t} (γ=${gamma}):`);
        console.log(`  Positive t_n (${t_n_pos.toFixed(1)}): cos(π t_n)=${cos_pos.toFixed(6)}, 1-exp=${one_minus_pos.toFixed(6)}, U_m≈${um_pos.toExponential(4)} J/m³`);
        console.log(`  Negative t_n (${t_n_neg.toFixed(1)}): cos(π t_n)=${cos_neg.toFixed(6)}, 1-exp=${one_minus_neg.toFixed(6)}, U_m≈${um_neg.toExponential(4)} J/m³`);
    }

    printTimeEvolution(times, gamma = null) {
        if (gamma === null) {
            gamma = this.variables.get('gamma');
        }

        console.log('\n=== Time Evolution with Negative Time Factor ===');
        for (const t of times) {
            const t_n = this.computeT_n(t);
            const cos_term = this.computeCosPiTn(t);
            const exp_term = this.computeExpTerm(gamma, t);
            const one_minus = this.computeOneMinusExp(gamma, t);
            console.log(`\nTime t = ${t.toExponential(2)} days:`);
            console.log(`  t_n = ${t_n.toExponential(2)}`);
            console.log(`  cos(π t_n) = ${cos_term.toFixed(6)}`);
            console.log(`  exp(-γ t cos(π t_n)) = ${exp_term.toExponential(4)}`);
            console.log(`  1 - exp = ${one_minus.toFixed(6)}`);
        }
    }

    printModuleInfo() {
        console.log('\n=== NegativeTimeModule Info ===');
        console.log(`  Module Version:     ${this.metadata.get('version')}`);
        console.log(`  Enhanced:           ${this.metadata.get('enhanced')}`);
        console.log(`  Variables:          ${this.variables.size}`);
        console.log(`  t_0:                ${this.variables.get('t_0').toExponential(4)} days`);
        console.log(`  γ:                  ${this.variables.get('gamma').toExponential(4)} day^-1`);
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
    NegativeTimeModule
};

// Example usage:
// const { NegativeTimeModule } = require('./source106.js');
// const mod = new NegativeTimeModule();
// const t = 1000.0;  // days
// mod.printTnEffects(t);
// console.log(mod.getEquationText());
// mod.updateVariable('t_0', 500.0);
// mod.printVariables();
