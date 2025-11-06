// UgIndexModule.js
// JavaScript implementation of the Index for Discrete Universal Gravity Ranges (i) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module uses i=1 to 4 to label Ug1-Ug4; computes sum_{i=1}^4 k_i * U_gi for F_U contribution.
// Variables in Map; defaults for Sun at t=0; i labels: 1=Internal Dipole, 2=Outer Bubble, 3=Magnetic Disk, 4=Star-BH.
// Approximations: k_i from coupling; sum ≈1.42e53 J/m³ (Ug2 dominant).
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

class UgIndexModule {
    constructor() {
        // Core variables using Map
        this.variables = new Map();
        
        // Coupling constants (unitless) - k1 to k4
        this.k_values = [1.5, 1.2, 1.8, 1.0];
        
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
        
        // U_gi defaults (J/m³, Sun t=0)
        this.variables.set('U_g1', 1.39e26);      // Internal Dipole
        this.variables.set('U_g2', 1.18e53);      // Outer Field Bubble
        this.variables.set('U_g3', 1.8e49);       // Magnetic Strings Disk
        this.variables.set('U_g4', 2.50e-20);     // Star-Black Hole Interactions
        
        // Shared params (placeholders)
        this.variables.set('t_n', 0.0);           // s
        this.variables.set('pi', Math.PI);
    }

    // ========== DYNAMIC VARIABLE OPERATIONS ==========
    
    updateVariable(name, value) {
        this.variables.set(name, value);
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
        } else {
            console.warn(`Variable '${name}' not found. Adding with delta ${delta}`);
            this.variables.set(name, delta);
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    // ========== CORE COMPUTATIONS ==========
    
    /**
     * Get the range of index i (1 to 4)
     * @returns {number} Maximum index value (4)
     */
    getIndexRange() {
        return 4;  // i=1 to 4
    }

    /**
     * Compute U_gi for a given index i
     * @param {number} i - Index (1-4): 1=Internal Dipole, 2=Outer Bubble, 3=Magnetic Disk, 4=Star-BH
     * @returns {number} U_gi in J/m³
     */
    computeU_gi(i) {
        const key = `U_g${i}`;
        if (this.variables.has(key)) {
            return this.variables.get(key);
        }
        console.error(`U_g${i} not found. Returning 0.`);
        return 0.0;
    }

    /**
     * Compute coupling constant k_i for a given index i
     * @param {number} i - Index (1-4)
     * @returns {number} k_i (unitless)
     */
    computeK_i(i) {
        if (i < 1 || i > 4) {
            console.error(`Invalid i: ${i}. Using k1.`);
            return this.k_values[0];
        }
        return this.k_values[i - 1];  // Convert to 0-based index
    }

    /**
     * Compute k_i * U_gi for a given index i
     * @param {number} i - Index (1-4)
     * @returns {number} k_i * U_gi in J/m³
     */
    computeKUgi(i) {
        return this.computeK_i(i) * this.computeU_gi(i);
    }

    /**
     * Compute sum of k_i * U_gi over index range
     * Σ_{i=i_min}^{i_max} k_i * U_gi for F_U contribution
     * @param {number} i_min - Minimum index (default 1)
     * @param {number} i_max - Maximum index (default 4)
     * @returns {number} Sum in J/m³
     */
    computeSumKUgi(i_min = 1, i_max = 4) {
        let sum = 0.0;
        
        // Add base contributions
        for (let i = i_min; i <= i_max; i++) {
            sum += this.computeKUgi(i);
        }
        
        // Apply dynamic terms if enabled
        if (this.enableDynamicTerms && this.dynamicTerms.length > 0) {
            const t = this.variables.get('t_n') || 0.0;
            for (const term of this.dynamicTerms) {
                sum += term.compute(t, this.variables);
            }
        }
        
        return sum;
    }

    /**
     * Get label for index i
     * @param {number} i - Index (1-4)
     * @returns {string} Label
     */
    getIndexLabel(i) {
        const labels = {
            1: 'Internal Dipole',
            2: 'Outer Field Bubble',
            3: 'Magnetic Strings Disk',
            4: 'Star-Black Hole'
        };
        return labels[i] || 'Unknown';
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
            k_values: this.k_values.slice(),
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
        if (state.k_values) {
            this.k_values = state.k_values.slice();
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
        return `F_U = Σ_{i=1}^4 [k_i * U_gi(r,t,M_s,ρ_s,T_s,B_s,ρ_vac,[SCm],ρ_vac,[UA],t_n) - ρ_i * ... ] + other terms
i (dimensionless integer): Labels Ug ranges; i=1: Internal Dipole, i=2: Outer Bubble,
i=3: Magnetic Disk, i=4: Star-BH.
Discretizes gravity for summation; enables scale-specific modeling.
Example Sun t=0: Σ k_i U_gi ≈1.42e53 J/m³ (Ug2 dominant).
Role: Structures Ug contributions; extensible for more ranges.`;
    }

    /**
     * Print all current variables
     */
    printVariables() {
        console.log('\n=== UgIndexModule Variables ===');
        for (const [key, value] of this.variables) {
            console.log(`  ${key} = ${value.toExponential(4)}`);
        }
        console.log(`  k_values: [k1=${this.k_values[0]}, k2=${this.k_values[1]}, k3=${this.k_values[2]}, k4=${this.k_values[3]}]`);
    }

    /**
     * Print breakdown by index i
     */
    printIndexBreakdown() {
        console.log('\n=== Ug Index Breakdown (i=1 to 4) ===');
        for (let i = 1; i <= 4; i++) {
            const ugi = this.computeU_gi(i);
            const ki = this.computeK_i(i);
            const kugi = this.computeKUgi(i);
            const label = this.getIndexLabel(i);
            console.log(`  i=${i} (${label}):`);
            console.log(`    U_g${i} = ${ugi.toExponential(4)} J/m³`);
            console.log(`    k${i} = ${ki.toFixed(2)}`);
            console.log(`    k_i * U_gi = ${kugi.toExponential(4)} J/m³`);
        }
        const sum = this.computeSumKUgi();
        console.log(`\n  Total Σ k_i * U_gi = ${sum.toExponential(4)} J/m³`);
    }

    /**
     * Print component contributions
     */
    printComponentContributions() {
        console.log('\n=== Component Contributions ===');
        const sum = this.computeSumKUgi();
        for (let i = 1; i <= 4; i++) {
            const kugi = this.computeKUgi(i);
            const percentage = (kugi / sum * 100).toFixed(2);
            const label = this.getIndexLabel(i);
            console.log(`  ${label} (i=${i}): ${percentage}%`);
        }
    }

    /**
     * Print module information
     */
    printModuleInfo() {
        console.log('\n=== UgIndexModule Info ===');
        console.log(`  Module Version:     ${this.metadata.get('version')}`);
        console.log(`  Enhanced:           ${this.metadata.get('enhanced')}`);
        console.log(`  Variables:          ${this.variables.size}`);
        console.log(`  Index Range:        1 to ${this.getIndexRange()}`);
        console.log(`  Coupling Constants: k=[${this.k_values.join(', ')}]`);
        console.log(`  Dynamic Terms:      ${this.dynamicTerms.length}`);
        console.log(`  Dynamic Parameters: ${this.dynamicParameters.size}`);
        console.log(`  Logging Enabled:    ${this.enableLogging}`);
        console.log(`  Learning Rate:      ${this.learningRate}`);
    }
}

// Export classes
module.exports = {
    UgIndexModule,
    PhysicsTerm,
    DynamicVacuumTerm,
    QuantumCouplingTerm
};

// Example usage:
// const { UgIndexModule } = require('./source102.js');
// const mod = new UgIndexModule();
// const sum = mod.computeSumKUgi();
// console.log(`Σ k_i U_gi = ${sum} J/m³`);
// mod.printIndexBreakdown();
// console.log(mod.getEquationText());
// mod.updateVariable('U_g3', 2e49);
// mod.printVariables();
// Sample: Sum ≈1.42e53 J/m³; i structures 4 Ug ranges.
