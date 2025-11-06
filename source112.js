// ScmPenetrationModule.js
// JavaScript implementation of the [SCm] Penetration Factor (P_SCm) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes P_SCm ≈1 (unitless for Sun, ~1e-3 for planets); scales P_SCm in Universal Magnetism U_m term.
// Variables in Map; example for Sun at t=0, t_n=0; full penetration for plasma cores.
// Approximations: 1 - e^{-γ t cos(π t_n)}=0 at t=0; φ_hat_j=1; μ_j / r_j=2.26e10 T m⁻¹.
// Converted from C++ source112.cpp
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

class ScmPenetrationModule {
    constructor() {
        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();
        
        // Universal constants - [SCm] Penetration Factor
        this.variables.set('P_SCm', 1.0);                      // Unitless ≈1 for Sun (full plasma penetration)
        this.variables.set('P_SCm_planet', 1e-3);              // For planets (solid core, 1000× reduction)
        this.variables.set('mu_j', 3.38e23);                   // T·m³ (magnetic moment)
        this.variables.set('r_j', 1.496e13);                   // m (solar radius scale)
        this.variables.set('gamma', 5e-5 / 86400.0);           // s⁻¹ (decay rate from source111)
        this.variables.set('t_n', 0.0);                        // Normalized time for reciprocation
        this.variables.set('phi_hat_j', 1.0);                  // Normalized field direction
        this.variables.set('E_react', 1e46);                   // J (reactive energy)
        this.variables.set('f_Heaviside', 0.01);               // Unitless (Heaviside factor from source109)
        this.variables.set('f_quasi', 0.01);                   // Unitless (quasi factor from source109)
        this.variables.set('pi', Math.PI);
        this.variables.set('scale_Heaviside', 1e13);           // Amplification scale

        // Derived
        this.variables.set('heaviside_factor', 
            1.0 + this.variables.get('scale_Heaviside') * this.variables.get('f_Heaviside'));

        // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;
        
        this.metadata.set('enhanced', 'true');
        this.metadata.set('version', '2.0-Enhanced');
        this.metadata.set('module', 'ScmPenetrationModule');
        this.metadata.set('description', '[SCm] Penetration Factor for UQFF Universal Magnetism');
    }

    // ========== DYNAMIC VARIABLE OPERATIONS ==========

    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
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
     * Compute P_SCm penetration factor
     * Returns ≈1 for Sun (unitless), ~1e-3 for planets
     */
    computeP_SCm() {
        return this.variables.get('P_SCm');
    }

    /**
     * Compute base U_m without additional factors
     * Formula: (μ_j / r_j) × (1 - e^{-γ t cos(π t_n)}) × φ_hat_j × P_SCm × E_react
     */
    computeUmBase(t) {
        const mu_over_rj = this.variables.get('mu_j') / this.variables.get('r_j');
        const gamma = this.variables.get('gamma');
        const t_n = this.variables.get('t_n');
        const pi = this.variables.get('pi');
        const exp_arg = -gamma * t * Math.cos(pi * t_n);
        const one_minus_exp = 1.0 - Math.exp(exp_arg);
        const phi_hat = this.variables.get('phi_hat_j');
        const p_scm = this.computeP_SCm();
        const e_react = this.variables.get('E_react');
        
        return mu_over_rj * one_minus_exp * phi_hat * p_scm * e_react;
    }

    /**
     * Compute U_m contribution with P_SCm, Heaviside, and quasi factors
     * Formula: U_m_base × (1 + 10^13 × f_Heaviside) × (1 + f_quasi)
     * Returns J/m³
     */
    computeUmContribution(t) {
        const base = this.computeUmBase(t);
        const heaviside_f = this.variables.get('heaviside_factor');
        const quasi_f = 1.0 + this.variables.get('f_quasi');
        return base * heaviside_f * quasi_f;
    }

    /**
     * Compute U_m for planet (P_SCm=1e-3)
     * Temporarily sets P_SCm to planetary value
     */
    computeUmPlanet(t) {
        const orig_p = this.variables.get('P_SCm');
        this.variables.set('P_SCm', this.variables.get('P_SCm_planet'));
        const result = this.computeUmContribution(t);
        this.variables.set('P_SCm', orig_p);
        return result;
    }

    /**
     * Compute stellar vs planetary scaling ratio
     */
    computeScalingRatio() {
        return this.variables.get('P_SCm') / this.variables.get('P_SCm_planet');
    }

    /**
     * Print comparison between stellar and planetary U_m
     */
    printComparison(t) {
        const um_star = this.computeUmContribution(t);
        const um_planet = this.computeUmPlanet(t);
        const ratio = um_star / um_planet;
        
        console.log(`\n=== [SCm] Penetration Comparison at t=${t.toExponential(3)} s ===`);
        console.log(`  P_SCm (stellar):   ${this.variables.get('P_SCm').toExponential(3)}`);
        console.log(`  P_SCm (planetary): ${this.variables.get('P_SCm_planet').toExponential(3)}`);
        console.log(`  U_m (stellar):     ${um_star.toExponential(3)} J/m³`);
        console.log(`  U_m (planetary):   ${um_planet.toExponential(3)} J/m³`);
        console.log(`  Scaling ratio:     ${ratio.toExponential(3)} (${ratio.toFixed(0)}×)`);
    }

    // ========== OUTPUT AND DEBUGGING ==========

    getEquationText() {
        return `U_m = [ (μ_j / r_j) (1 - e^{-γ t cos(π t_n)}) φ_hat_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)
Where P_SCm ≈1 (unitless [SCm] penetration factor; ~1e-3 for planets).
Scales magnetic energy for [SCm] interior interaction.
Example Sun t=0: U_m ≈2.28e65 J/m³ (P_SCm=1);
Planet: ≈2.28e62 J/m³ (P_SCm=1e-3, -3 orders).
Role: Full for stellar plasma, reduced for solid cores; [SCm] influence on strings.
UQFF: Models penetration in jets/nebulae/formation; massless [SCm] dynamics.`;
    }

    printVariables() {
        console.log('\n=== ScmPenetrationModule Variables ===');
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
    ScmPenetrationModule,
    PhysicsTerm,
    DynamicVacuumTerm,
    QuantumCouplingTerm
};
