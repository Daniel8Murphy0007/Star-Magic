// ReciprocationDecayModule.js
// JavaScript conversion of ReciprocationDecayModule.cpp
// Modular implementation of the Reciprocation Decay Rate (γ) in the UQFF framework.
// This module computes γ=0.00005 day⁻¹ (~5.8e-10 s⁻¹); used in exp(-γ t cos(π t_n)) for U_m decay.
// Variables in Map; example for t=1000 days, t_n=0; 1-exp ≈0.049.
// Approximations: cos(π t_n)=1; timescale ~55 years; μ_j / r_j=2.26e10 T m².
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

class ReciprocationDecayModule {
    constructor() {
        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();
        
        // Universal constants
        this.variables.set('gamma_day', 0.00005);               // day⁻¹
        this.variables.set('day_to_s', 86400.0);                // s/day
        this.variables.set('t_n', 0.0);                         // days
        this.variables.set('t_day', 0.0);                       // days
        this.variables.set('pi', 3.141592653589793);
        this.variables.set('mu_over_rj', 2.26e10);              // T m²
        this.variables.set('P_SCm', 1.0);                       // Normalized
        this.variables.set('E_react', 1e46);                    // J
        this.variables.set('heaviside_f', 1e11 + 1.0);          // 1 + 10^13 * 0.01
        this.variables.set('quasi_f', 1.01);                    // 1 + 0.01

        // Derived
        this.variables.set('gamma_s', this.computeGamma_s());

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
            if (name === 'gamma_day') {
                this.variables.set('gamma_s', this.computeGamma_s());
            }
        } else {
            console.warn(`Variable '${name}' not found. Adding with value ${value}`);
            this.variables.set(name, value);
        }
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
            if (name === 'gamma_day') {
                this.variables.set('gamma_s', this.computeGamma_s());
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

    computeGamma_day() {
        return this.variables.get('gamma_day');
    }

    computeGamma_s() {
        return this.computeGamma_day() / this.variables.get('day_to_s');
    }

    computeCosPiTn(t_n) {
        this.variables.set('t_n', t_n);
        return Math.cos(this.variables.get('pi') * t_n);
    }

    computeExpTerm(t_day, t_n) {
        this.variables.set('t_day', t_day);
        const cos_pi_tn = this.computeCosPiTn(t_n);
        const arg = -this.computeGamma_day() * t_day * cos_pi_tn;
        return Math.exp(arg);
    }

    computeOneMinusExp(t_day, t_n) {
        return 1.0 - this.computeExpTerm(t_day, t_n);
    }

    computeUmExample(t_day, t_n, mu_over_rj = null) {
        if (mu_over_rj === null) {
            mu_over_rj = this.variables.get('mu_over_rj');
        }
        const one_minus_exp = this.computeOneMinusExp(t_day, t_n);
        const phi_hat = 1.0;
        const p_scm = this.variables.get('P_SCm');
        const e_react = this.variables.get('E_react');
        const heaviside_f = this.variables.get('heaviside_f');
        const quasi_f = this.variables.get('quasi_f');
        return (mu_over_rj * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
    }

    // ========== OUTPUT AND DEBUGGING ==========

    getEquationText() {
        return `γ = 0.00005 day⁻¹ (~5.8e-10 s⁻¹; timescale ~55 years);
In U_m: ... (1 - exp(-γ t cos(π t_n))) ... (t days, reciprocating decay/growth).
Negative cos(π t_n): exp(+γ t) >1 (growth, negentropic TRZ).
Example t=1000 days, t_n=0: 1-exp ≈0.049, U_m ≈1.12e66 J/m³.
UQFF: Slow decay for magnetic strings; cyclic via cos(π t_n) in jets/nebulae/mergers.`;
    }

    printVariables() {
        console.log('Current Variables:');
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential(3)}`);
        }
    }

    printDecayEffects(t_day = 1000.0, t_n = 0.0) {
        const cos_pi = this.computeCosPiTn(t_n);
        const exp_val = this.computeExpTerm(t_day, t_n);
        const one_minus = this.computeOneMinusExp(t_day, t_n);
        const um_ex = this.computeUmExample(t_day, t_n);
        console.log(`Decay Effects at t=${t_day} days, t_n=${t_n}:`);
        console.log(`cos(π t_n) = ${cos_pi.toFixed(6)}`);
        console.log(`exp(-γ t cos(π t_n)) = ${exp_val.toExponential(6)}`);
        console.log(`1 - exp(...) = ${one_minus.toExponential(6)}`);
        console.log(`U_m example contrib = ${um_ex.toExponential(3)} J/m³`);
    }

    computeTimescale() {
        // Timescale is approximately 1/γ in days
        return 1.0 / this.computeGamma_day();
    }

    computeTimescaleYears() {
        return this.computeTimescale() / 365.25;
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
    module.exports = { ReciprocationDecayModule, PhysicsTerm, DynamicVacuumTerm, QuantumCouplingTerm };
}

// Example usage:
// const { ReciprocationDecayModule } = require('./source111.js');
// const mod = new ReciprocationDecayModule();
// const gamma = mod.computeGamma_day();
// console.log('γ =', gamma, 'day⁻¹');
// mod.printDecayEffects(1000.0, 0.0);
// console.log(mod.getEquationText());
// mod.updateVariable('gamma_day', 0.0001);
// mod.printVariables();
// Sample: γ=5e-5 day⁻¹; t=1000 days: 1-exp≈0.049; U_m≈1.12e66 J/m³.
