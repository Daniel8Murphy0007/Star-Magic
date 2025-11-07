// ScmReactivityDecayModule.js
// JavaScript implementation of the [SCm] Reactivity Decay Rate (κ) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes κ=0.0005 day⁻¹ (~5.8e-6 s⁻¹); used in E_react = 10^46 * exp(-κ t) for decay in U_m, U_bi, etc.
// Variables in Map; example for Sun at t=0 (E_react=1e46); t=2000 days: ~3.68e45.
// Approximations: t in days; timescale ~5.5 years; integrates into U_m example.
// Converted from C++ source113.cpp
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

class ScmReactivityDecayModule {
    constructor() {
        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();
        
        // Universal constants - [SCm] Reactivity Decay Rate
        this.variables.set('kappa_day', 0.0005);                // day⁻¹ (reactivity decay rate)
        this.variables.set('day_to_s', 86400.0);                // s/day conversion
        this.variables.set('E_react_base', 1e46);               // J (base reactive energy)
        this.variables.set('t_day', 0.0);                       // days (current time)
        this.variables.set('mu_over_rj', 2.26e10);              // T/m (magnetic moment ratio)
        this.variables.set('P_SCm', 1.0);                       // Normalized penetration factor
        this.variables.set('heaviside_f', 1e11 + 1.0);          // 1 + 10^13 × 0.01
        this.variables.set('quasi_f', 1.01);                    // 1 + 0.01
        this.variables.set('one_minus_exp', 1.0);               // At t=0 (placeholder for full computation)

        // Derived
        this.variables.set('kappa_s', this.computeKappa_s());   // s⁻¹

        // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;
        
        this.metadata.set('enhanced', 'true');
        this.metadata.set('version', '2.0-Enhanced');
        this.metadata.set('module', 'ScmReactivityDecayModule');
        this.metadata.set('description', '[SCm] Reactivity Decay Rate for UQFF E_react term');
    }

    // ========== DYNAMIC VARIABLE OPERATIONS ==========

    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
            if (name === 'kappa_day') {
                this.variables.set('kappa_s', this.computeKappa_s());
            }
        } else {
            if (this.enableLogging) {
                console.warn(`Variable '${name}' not found. Adding with value ${value}`);
            }
            this.variables.set(name, value);
        }
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
            if (name === 'kappa_day') {
                this.variables.set('kappa_s', this.computeKappa_s());
            }
        } else {
            if (this.enableLogging) {
                console.warn(`Variable '${name}' not found. Adding with delta ${delta}`);
            }
            this.variables.set(name, delta);
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    getVariable(name) {
        return this.variables.get(name);
    }

    // ========== CORE COMPUTATIONS ==========

    /**
     * Compute κ (day⁻¹) - reactivity decay rate
     * Returns 0.0005 day⁻¹ by default
     */
    computeKappa_day() {
        return this.variables.get('kappa_day');
    }

    /**
     * Compute κ in s⁻¹
     * Formula: κ_day / 86400
     * Returns ~5.8e-6 s⁻¹
     */
    computeKappa_s() {
        return this.computeKappa_day() / this.variables.get('day_to_s');
    }

    /**
     * Compute E_react with decay
     * Formula: E_react = 10^46 × exp(-κ t)
     * @param {number} t_day - Time in days
     * @returns {number} E_react in Joules
     */
    computeE_react(t_day) {
        this.variables.set('t_day', t_day);
        const kappa = this.computeKappa_day();
        const arg = -kappa * t_day;
        return this.variables.get('E_react_base') * Math.exp(arg);
    }

    /**
     * Compute decay timescale (1/κ in days)
     * Returns ~2000 days (~5.5 years)
     */
    computeTimescale() {
        return 1.0 / this.computeKappa_day();
    }

    /**
     * Compute timescale in years
     */
    computeTimescaleYears() {
        return this.computeTimescale() / 365.25;
    }

    /**
     * Compute decay fraction at given time
     * Formula: exp(-κ t)
     */
    computeDecayFraction(t_day) {
        const kappa = this.computeKappa_day();
        return Math.exp(-kappa * t_day);
    }

    /**
     * Simplified U_m example with E_react
     * Formula: (μ/r × (1-exp) × φ_hat) × P_SCm × E_react × Heaviside × quasi
     * @param {number} t_day - Time in days
     * @returns {number} U_m in J/m³
     */
    computeUmExample(t_day) {
        const e_react = this.computeE_react(t_day);
        const one_minus_exp = this.variables.get('one_minus_exp');  // Placeholder; full would compute from source111
        const phi_hat = 1.0;
        const p_scm = this.variables.get('P_SCm');
        const heaviside_f = this.variables.get('heaviside_f');
        const quasi_f = this.variables.get('quasi_f');
        
        return (this.variables.get('mu_over_rj') * one_minus_exp * phi_hat) * 
               p_scm * e_react * heaviside_f * quasi_f;
    }

    /**
     * Print decay effects at given time
     */
    printDecayEffects(t_day = 2000.0) {
        const e_react = this.computeE_react(t_day);
        const um_ex = this.computeUmExample(t_day);
        const fraction = e_react / this.variables.get('E_react_base');
        
        console.log(`\n=== [SCm] Reactivity Decay Effects at t=${t_day.toFixed(1)} days ===`);
        console.log(`  κ (day⁻¹):         ${this.computeKappa_day().toExponential(3)}`);
        console.log(`  κ (s⁻¹):           ${this.computeKappa_s().toExponential(3)}`);
        console.log(`  Timescale:         ${this.computeTimescale().toFixed(1)} days (~${this.computeTimescaleYears().toFixed(2)} years)`);
        console.log(`  E_react:           ${e_react.toExponential(3)} J`);
        console.log(`  Decay fraction:    ${fraction.toFixed(4)} (${(fraction * 100).toFixed(2)}% of initial)`);
        console.log(`  U_m example:       ${um_ex.toExponential(3)} J/m³`);
    }

    // ========== OUTPUT AND DEBUGGING ==========

    getEquationText() {
        return `E_react = 10^46 × exp(-κ t) (t days); κ=0.0005 day⁻¹ (~5.8e-6 s⁻¹, timescale ~5.5 years).
In U_m, U_bi, U_i, U_gi: ... × E_react × ... (decays [SCm] reactivity).
Example t=0: E_react=1e46 J; t=2000 days: ~3.68e45 J (~36.8%).
U_m (t=0): ≈2.28e65 J/m³; t=2000: ≈8.39e64 J/m³.
Role: Gradual [SCm]-[UA] interaction loss; temporal evolution in jets/nebulae/mergers.
UQFF: Models reactivity decay; energy dissipation over cosmic time.`;
    }

    printVariables() {
        console.log('\n=== ScmReactivityDecayModule Variables ===');
        for (const [key, value] of this.variables) {
            console.log(`  ${key.padEnd(20)} = ${typeof value === 'number' ? value.toExponential(4) : value}`);
        }
    }

    // ========== SELF-EXPANDING FRAMEWORK: DYNAMIC TERMS ==========

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

    computeDynamicTerms(t) {
        if (!this.enableDynamicTerms || this.dynamicTerms.length === 0) {
            return 0.0;
        }
        
        let total = 0.0;
        for (const term of this.dynamicTerms) {
            if (term.validate(this.variables)) {
                total += term.compute(t, this.variables);
            }
        }
        return total;
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
        return JSON.stringify(state, null, 2);
    }

    importState(stateJson) {
        try {
            const state = JSON.parse(stateJson);
            this.variables = new Map(state.variables);
            this.dynamicParameters = new Map(state.dynamicParameters);
            this.metadata = new Map(state.metadata);
            this.enableDynamicTerms = state.enableDynamicTerms;
            this.enableLogging = state.enableLogging;
            this.learningRate = state.learningRate;
            return true;
        } catch (error) {
            console.error('Failed to import state:', error);
            return false;
        }
    }
}

// ===========================================================================================
// MODULE EXPORTS
// ===========================================================================================

module.exports = {
    ScmReactivityDecayModule,
    PhysicsTerm,
    DynamicVacuumTerm,
    QuantumCouplingTerm
};
