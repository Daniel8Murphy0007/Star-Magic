// source114.js
// SolarCycleFrequencyModule - Solar Cycle Frequency (ω_c) UQFF Module
// Converted from source114.cpp - Maintains all self-expanding dynamics
// 
// Physics: ω_c = 2π / 3.96e8 s⁻¹ ≈ 1.59e-8 rad/s (period ~12.55 years)
// Models solar cycle periodicity; used in sin(ω_c t) for μ_j in U_m
// μ_j = (10³ + 0.4 sin(ω_c t)) × 3.38e20 T·m³ (cyclic magnetic variation)
// Near 11-year Hale solar cycle; affects jets/nebulae/stellar formation
//
// Self-Expanding Framework: Dynamic variables, physics terms, state management
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025

// ===========================================================================================
// SELF-EXPANDING FRAMEWORK: Dynamic Physics Term System
// ===========================================================================================

class PhysicsTerm {
    constructor() {
        this.dynamicParameters = new Map();
        this.metadata = new Map();
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
    constructor(couplingStrength = 1e-40) {
        super();
        this.couplingStrength = couplingStrength;
    }

    compute(t, params) {
        const hbar = params.has('hbar') ? params.get('hbar') : 1.0546e-34;
        const M = params.has('M') ? params.get('M') : 1.989e30;
        const r = params.has('r') ? params.get('r') : 1e4;
        return this.couplingStrength * (hbar * hbar) / (M * r * r) * Math.cos(t / 1e6);
    }

    getName() {
        return 'QuantumCoupling';
    }

    getDescription() {
        return 'Non-local quantum effects';
    }
}

// ===========================================================================================
// ENHANCED SOLAR CYCLE FREQUENCY MODULE WITH SELF-EXPANDING CAPABILITIES
// ===========================================================================================

class SolarCycleFrequencyModule {
    constructor() {
        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();
        
        // Universal constants
        this.variables.set('pi', 3.141592653589793);
        this.variables.set('period', 3.96e8);              // s (~12.55 years)
        this.variables.set('base_mu', 3.38e20);            // T·m³
        this.variables.set('B_j', 1e3);                    // Base T
        this.variables.set('t', 0.0);                      // s

        // Derived (computed)
        this.variables.set('omega_c', this.computeOmega_c());

        // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;

        // Framework metadata
        this.metadata.set('enhanced', 'true');
        this.metadata.set('version', '2.0-Enhanced');
        this.metadata.set('module', 'SolarCycleFrequencyModule');
        this.metadata.set('physics', 'Solar cycle frequency ω_c for periodic magnetic variations');
    }

    // ========== CORE COMPUTATIONS (Original UQFF - Preserved) ==========

    /**
     * Compute ω_c = 2π / period
     * @returns {number} Angular frequency in rad/s
     */
    computeOmega_c() {
        return 2.0 * this.variables.get('pi') / this.variables.get('period');
    }

    /**
     * Compute sin(ω_c t)
     * @param {number} t - Time in seconds
     * @returns {number} sin(ω_c t)
     */
    computeSinOmegaCT(t) {
        this.variables.set('t', t);
        return Math.sin(this.variables.get('omega_c') * t);
    }

    /**
     * Example μ_j = (10³ + 0.4 sin(ω_c t)) × 3.38e20
     * Models cyclic magnetic variation with solar cycle
     * @param {number} t - Time in seconds
     * @returns {number} Magnetic moment in T·m³
     */
    computeMuJExample(t) {
        const sin_omega = this.computeSinOmegaCT(t);
        const b_j = this.variables.get('B_j') + 0.4 * sin_omega;
        return b_j * this.variables.get('base_mu');
    }

    /**
     * Compute period in years
     * @returns {number} Period in years
     */
    computePeriodYears() {
        return this.variables.get('period') / (365.25 * 86400);
    }

    /**
     * Compute frequency in Hz
     * @returns {number} Frequency in Hz
     */
    computeFrequencyHz() {
        return this.variables.get('omega_c') / (2.0 * this.variables.get('pi'));
    }

    /**
     * Compute magnetic field variation at time t
     * @param {number} t - Time in seconds
     * @returns {number} B_j(t) in Tesla
     */
    computeBJVariation(t) {
        const sin_omega = this.computeSinOmegaCT(t);
        return this.variables.get('B_j') + 0.4 * sin_omega;
    }

    /**
     * Compute percentage change in μ_j from baseline
     * @param {number} t - Time in seconds
     * @returns {number} Percentage change
     */
    computeMuJPercentChange(t) {
        const mu_t = this.computeMuJExample(t);
        const mu_0 = this.computeMuJExample(0);
        return ((mu_t - mu_0) / mu_0) * 100;
    }

    // ========== DYNAMIC VARIABLE OPERATIONS ==========

    /**
     * Update a variable value
     * @param {string} name - Variable name
     * @param {number} value - New value
     */
    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
            if (name === 'period') {
                this.variables.set('omega_c', this.computeOmega_c());
            }
        } else {
            console.warn(`Variable '${name}' not found. Adding with value ${value}`);
            this.variables.set(name, value);
        }
    }

    /**
     * Add delta to variable
     * @param {string} name - Variable name
     * @param {number} delta - Delta to add
     */
    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
            if (name === 'period') {
                this.variables.set('omega_c', this.computeOmega_c());
            }
        } else {
            console.warn(`Variable '${name}' not found. Adding with delta ${delta}`);
            this.variables.set(name, delta);
        }
    }

    /**
     * Subtract delta from variable
     * @param {string} name - Variable name
     * @param {number} delta - Delta to subtract
     */
    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    // ========== SELF-EXPANDING FRAMEWORK METHODS ==========

    /**
     * Register a dynamic physics term
     * @param {PhysicsTerm} term - Physics term to add
     */
    registerDynamicTerm(term) {
        this.dynamicTerms.push(term);
        if (this.enableLogging) {
            console.log(`Registered dynamic term: ${term.getName()}`);
        }
    }

    /**
     * Compute contribution from all dynamic terms
     * @param {number} t - Time parameter
     * @returns {number} Sum of dynamic term contributions
     */
    computeDynamicTerms(t) {
        if (!this.enableDynamicTerms || this.dynamicTerms.length === 0) {
            return 0;
        }

        let total = 0;
        const params = new Map(this.variables);
        
        for (const term of this.dynamicTerms) {
            total += term.compute(t, params);
        }

        return total;
    }

    /**
     * Set a dynamic parameter
     * @param {string} name - Parameter name
     * @param {number} value - Parameter value
     */
    setDynamicParameter(name, value) {
        this.dynamicParameters.set(name, value);
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
     * Set learning rate for optimization
     * @param {number} rate - Learning rate
     */
    setLearningRate(rate) {
        this.learningRate = rate;
    }

    /**
     * Enable or disable logging
     * @param {boolean} enable - Enable logging
     */
    setEnableLogging(enable) {
        this.enableLogging = enable;
    }

    /**
     * Export module state
     * @returns {string} JSON state
     */
    exportState() {
        return JSON.stringify({
            variables: Array.from(this.variables.entries()),
            dynamicParameters: Array.from(this.dynamicParameters.entries()),
            metadata: Array.from(this.metadata.entries()),
            enableDynamicTerms: this.enableDynamicTerms,
            enableLogging: this.enableLogging,
            learningRate: this.learningRate
        });
    }

    /**
     * Import module state
     * @param {string} stateJson - JSON state
     * @returns {boolean} Success
     */
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

    // ========== OUTPUT AND DESCRIPTIVE METHODS ==========

    /**
     * Get equation text
     * @returns {string} Equation description
     */
    getEquationText() {
        return `ω_c = 2π / 3.96e8 s⁻¹ ≈ 1.59e-8 rad/s (period ~12.55 yr, near 11-yr solar cycle);
In U_m: μ_j = (10³ + 0.4 sin(ω_c t)) × 3.38e20 T·m³ (cyclic magnetic variation).
In U_g3: ... cos(ω_s t π) ... (ω_s Sun rotation, but ω_c for cycle).
Example t=0: sin=0 → μ_j=3.38e23 T·m³;
t=3.14e7 s (~1 yr): sin≈0.477 → μ_j≈3.381e23 T·m³ (+0.019%).
Role: Models solar cycle periodicity; magnetic activity in strings/fields.
UQFF: Cyclic effects in jets/nebulae/formation; near 11-yr Hale cycle.`;
    }

    /**
     * Print all current variables
     */
    printVariables() {
        console.log('Current Variables:');
        for (const [key, value] of this.variables) {
            console.log(`  ${key} = ${value.toExponential(3)}`);
        }
    }

    /**
     * Print solar cycle effects at time t
     * @param {number} t - Time in seconds
     */
    printSolarCycleEffects(t) {
        console.log(`\n=== Solar Cycle Effects at t = ${t.toExponential(2)} s ===`);
        console.log(`ω_c (rad/s): ${this.variables.get('omega_c').toExponential(3)}`);
        console.log(`Period (years): ${this.computePeriodYears().toFixed(2)}`);
        console.log(`Frequency (Hz): ${this.computeFrequencyHz().toExponential(3)}`);
        console.log(`sin(ω_c t): ${this.computeSinOmegaCT(t).toFixed(6)}`);
        console.log(`B_j(t) (T): ${this.computeBJVariation(t).toExponential(3)}`);
        console.log(`μ_j(t) (T·m³): ${this.computeMuJExample(t).toExponential(3)}`);
        console.log(`Δμ_j from baseline: ${this.computeMuJPercentChange(t).toFixed(6)}%`);
    }
}

// ===========================================================================================
// MODULE EXPORTS
// ===========================================================================================

module.exports = {
    SolarCycleFrequencyModule,
    DynamicVacuumTerm,
    QuantumCouplingTerm,
    PhysicsTerm
};
