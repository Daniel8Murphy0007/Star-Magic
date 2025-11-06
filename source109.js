// QuasiLongitudinalModule.js
// JavaScript conversion of QuasiLongitudinalModule.cpp
// Modular implementation of the Quasi-Longitudinal Wave Factor (f_quasi) in the UQFF framework.
// This module computes f_quasi=0.01 (unitless) and its scaling (1 + f_quasi) in Universal Magnetism U_m term.
// Variables in Map; example for Sun at t=0, t_n=0; minor 1% increase in U_m.
// Approximations: 1 - e^{-γ t cos(π t_n)}=0 at t=0; φ_hat_j=1; P_SCm=1; f_Heaviside=0.01 (1 + 10^13 f=1e11+1).
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

class QuasiLongitudinalModule {
    constructor() {
        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();
        
        // Universal constants
        this.variables.set('f_quasi', 0.01);                    // Unitless fraction
        this.variables.set('mu_j', 3.38e23);                    // T·m^3 (j=1)
        this.variables.set('r_j', 1.496e13);                    // m
        this.variables.set('gamma', 5e-5 / 86400.0);            // s^-1 (0.00005 day^-1)
        this.variables.set('t_n', 0.0);                         // s
        this.variables.set('phi_hat_j', 1.0);                   // Normalized
        this.variables.set('P_SCm', 1.0);                       // Pressure
        this.variables.set('E_react', 1e46);                    // J
        this.variables.set('f_Heaviside', 0.01);                // For Heaviside
        this.variables.set('scale_Heaviside', 1e13);            // Amplification
        this.variables.set('pi', 3.141592653589793);

        // Derived
        this.variables.set('quasi_factor', this.computeQuasiFactor());
        this.variables.set('heaviside_factor', 1.0 + this.variables.get('scale_Heaviside') * this.variables.get('f_Heaviside'));

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
            if (name === 'f_quasi') {
                this.variables.set('quasi_factor', this.computeQuasiFactor());
            }
        } else {
            console.warn(`Variable '${name}' not found. Adding with value ${value}`);
            this.variables.set(name, value);
        }
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
            if (name === 'f_quasi') {
                this.variables.set('quasi_factor', this.computeQuasiFactor());
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

    computeF_quasi() {
        return this.variables.get('f_quasi');
    }

    computeQuasiFactor() {
        return 1.0 + this.computeF_quasi();
    }

    computeUmBase(j, t) {
        const mu_over_rj = this.variables.get('mu_j') / this.variables.get('r_j');
        const exp_arg = -this.variables.get('gamma') * t * Math.cos(this.variables.get('pi') * this.variables.get('t_n'));
        const one_minus_exp = 1.0 - Math.exp(exp_arg);
        const phi_hat = this.variables.get('phi_hat_j');
        return mu_over_rj * one_minus_exp * phi_hat * this.variables.get('P_SCm') * this.variables.get('E_react');
    }

    computeUmContribution(j, t) {
        const base = this.computeUmBase(j, t);
        const quasi_f = this.computeQuasiFactor();
        const heaviside_f = this.variables.get('heaviside_factor');
        return base * heaviside_f * quasi_f;
    }

    computeUmWithNoQuasi(j, t) {
        const orig_f = this.variables.get('f_quasi');
        this.variables.set('f_quasi', 0.0);
        const result = this.computeUmContribution(j, t);
        this.variables.set('f_quasi', orig_f);
        return result;
    }

    // ========== OUTPUT AND DEBUGGING ==========

    getEquationText() {
        return `U_m = Σ_j [ (μ_j / r_j) (1 - e^{-γ t cos(π t_n)}) φ_hat_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)
Where f_quasi = 0.01 (unitless quasi-longitudinal wave factor);
Quasi factor = 1 + 0.01 = 1.01 (1% increase).
Example j=1, t=0: U_m contrib ≈2.28e65 J/m³ (with); ≈2.26e65 J/m³ (without; -1%).
Role: Minor scaling for quasi-longitudinal waves in magnetic strings; subtle [SCm]/[UA] wave effects.
UQFF: Enhances wave propagation in jets/nebulae; small but cumulative in dynamics.`;
    }

    printVariables() {
        console.log('Current Variables:');
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential(3)}`);
        }
    }

    printUmComparison(j = 1, t = 0.0) {
        const um_with = this.computeUmContribution(j, t);
        const um_without = this.computeUmWithNoQuasi(j, t);
        const percent_increase = ((um_with - um_without) / um_without) * 100.0;
        console.log(`U_m Comparison for j=${j} at t=${t} s:`);
        console.log(`With quasi: ${um_with.toExponential(3)} J/m³`);
        console.log(`Without quasi: ${um_without.toExponential(3)} J/m³`);
        console.log(`Increase: +${percent_increase.toFixed(1)}%`);
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
    module.exports = { QuasiLongitudinalModule, PhysicsTerm, DynamicVacuumTerm, QuantumCouplingTerm };
}

// Example usage:
// const { QuasiLongitudinalModule } = require('./source109.js');
// const mod = new QuasiLongitudinalModule();
// const quasi_f = mod.computeQuasiFactor();
// console.log('Quasi Factor =', quasi_f);
// mod.printUmComparison(1, 0.0);
// console.log(mod.getEquationText());
// mod.updateVariable('f_quasi', 0.02);
// mod.printVariables();
// Sample: Factor=1.01; U_m with=2.28e65 J/m³ (+1% vs without).
