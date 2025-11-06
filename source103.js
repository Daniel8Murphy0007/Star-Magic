// source103.js - InertiaCouplingModule
// JavaScript conversion of InertiaCouplingModule.cpp for Star-Magic UQFF Framework
// Computes λ_i=1.0 (unitless, uniform for i=1-4) and scales U_i in F_U: -λ_i [ρ_i U_i E_react]
// Variables: λ (lambda), ρ_vac_SCm, ρ_vac_UA, ω_s, f_TRZ, E_react, t_n, α_decay
// Example: Sun at t=0, t_n=0: U_i ≈1.38e-47 J/m³, contrib ≈ -0.138 J/m³ per i
// Sum over i=1 to 4: ≈ -0.552 J/m³
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025
// Converted to JavaScript: November 6, 2025

'use strict';

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
        const rho_vac = params.get('rho_vac_UA') || 7.09e-36;
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
        const hbar = params.get('hbar') || 1.0546e-34;
        const M = params.get('M') || 1.989e30;
        const r = params.get('r') || 1e4;
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

class InertiaCouplingModule {
    constructor() {
        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();
        
        // Universal constants for Sun at t=0, level 13
        this.variables.set('lambda', 1.0);                      // Uniform λ_i (unitless)
        this.variables.set('rho_vac_SCm', 7.09e-37);            // J/m³
        this.variables.set('rho_vac_UA', 7.09e-36);             // J/m³
        this.variables.set('omega_s', 2.5e-6);                  // rad/s (Sun rotation)
        this.variables.set('f_TRZ', 0.1);                       // Unitless
        this.variables.set('E_react', 1e46);                    // J
        this.variables.set('pi', Math.PI);
        this.variables.set('t_n', 0.0);                         // s
        this.variables.set('alpha_decay', 0.0005);              // For E_react exp decay
        
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
    
    // Compute λ_i (uniform 1.0)
    computeLambda_i(i) {
        return this.variables.get('lambda');
    }

    // Compute U_i for index i at time t (J/m³)
    computeU_i(i, t) {
        const lambda_i = this.computeLambda_i(i);
        const rho_sc = this.variables.get('rho_vac_SCm');
        const rho_ua = this.variables.get('rho_vac_UA');
        const omega_s_t = this.variables.get('omega_s');
        const cos_term = Math.cos(this.variables.get('pi') * this.variables.get('t_n'));
        const trz_factor = 1.0 + this.variables.get('f_TRZ');
        
        return lambda_i * rho_sc * rho_ua * omega_s_t * cos_term * trz_factor;
    }

    // Compute inertia term -λ_i U_i E_react (J/m³)
    computeInertiaTerm(i, t) {
        const u_i = this.computeU_i(i, t);
        const alpha = this.variables.get('alpha_decay');
        const e_react = this.variables.get('E_react') * Math.exp(-alpha * t);
        
        return -this.computeLambda_i(i) * u_i * e_react;
    }

    // Sum over i=1 to 4 (total F_U contribution, J/m³)
    computeSumInertiaTerms(t) {
        let sum = 0.0;
        for (let i = 1; i <= 4; i++) {
            sum += this.computeInertiaTerm(i, t);
        }
        
        // Add dynamic terms if enabled
        if (this.enableDynamicTerms && this.dynamicTerms.length > 0) {
            const params = this.variables;
            for (const term of this.dynamicTerms) {
                sum += term.compute(t, params);
            }
        }
        
        return sum;
    }

    // Get index label
    getIndexLabel(i) {
        const labels = {
            1: 'Ug1 (Internal Dipole)',
            2: 'Ug2 (Outer Field Bubble)',
            3: 'Ug3 (Magnetic Strings Disk)',
            4: 'Ug4 (Star-BH Interactions)'
        };
        return labels[i] || `i=${i}`;
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
        return {
            variables: Array.from(this.variables.entries()),
            dynamicParameters: Array.from(this.dynamicParameters.entries()),
            metadata: Array.from(this.metadata.entries()),
            enableDynamicTerms: this.enableDynamicTerms,
            learningRate: this.learningRate,
            numDynamicTerms: this.dynamicTerms.length
        };
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
        if (typeof state.enableDynamicTerms !== 'undefined') {
            this.enableDynamicTerms = state.enableDynamicTerms;
        }
        if (typeof state.learningRate !== 'undefined') {
            this.learningRate = state.learningRate;
        }
    }

    // ========== OUTPUT AND DISPLAY METHODS ==========
    
    getEquationText() {
        return `F_U = ... - λ_i [ρ_i * U_i * E_react] + ...
U_i = λ_i * ρ_vac,[SCm] * ρ_vac,[UA] * ω_s(t) * cos(π t_n) * (1 + f_TRZ)
Where λ_i = 1.0 (unitless, uniform for i=1-4: Ug1-Ug4);
E_react = 1e46 * e^{-α t} (α=5e-4);
Example Sun t=0, t_n=0: U_i ≈1.38e-47 J/m³; -λ_i U_i E_react ≈ -0.138 J/m³ (per i).
Role: Scales resistive inertia; uniform baseline opposition to dynamics.
UQFF: Consistent across scales; aids stability in interiors/disks/mergers.`;
    }

    printVariables() {
        console.log('\n=== InertiaCouplingModule Variables ===');
        for (const [name, value] of this.variables.entries()) {
            console.log(`  ${name} = ${value.toExponential(4)}`);
        }
    }

    printInertiaBreakdown(t = 0.0) {
        console.log(`\n=== Inertia Breakdown at t=${t} s ===`);
        for (let i = 1; i <= 4; i++) {
            const u_i = this.computeU_i(i, t);
            const term = this.computeInertiaTerm(i, t);
            console.log(`  ${this.getIndexLabel(i)}: U_i = ${u_i.toExponential(4)} J/m³, Term = ${term.toExponential(4)} J/m³`);
        }
        const sum = this.computeSumInertiaTerms(t);
        console.log(`  Total Σ Terms = ${sum.toExponential(4)} J/m³`);
    }

    printComponentContributions(t = 0.0) {
        console.log(`\n=== Component Contributions at t=${t} s ===`);
        const terms = [];
        let totalAbs = 0;
        
        for (let i = 1; i <= 4; i++) {
            const term = this.computeInertiaTerm(i, t);
            terms.push({ index: i, value: term, label: this.getIndexLabel(i) });
            totalAbs += Math.abs(term);
        }
        
        for (const item of terms) {
            const percent = totalAbs > 0 ? (Math.abs(item.value) / totalAbs * 100) : 0;
            console.log(`  ${item.label}: ${percent.toFixed(2)}%`);
        }
    }

    printModuleInfo() {
        console.log('\n=== InertiaCouplingModule Info ===');
        console.log(`  Module Version:     ${this.metadata.get('version')}`);
        console.log(`  Enhanced:           ${this.metadata.get('enhanced')}`);
        console.log(`  Variables:          ${this.variables.size}`);
        console.log(`  λ (lambda):         ${this.variables.get('lambda')}`);
        console.log(`  Dynamic Terms:      ${this.dynamicTerms.length}`);
        console.log(`  Dynamic Parameters: ${this.dynamicParameters.size}`);
        console.log(`  Logging Enabled:    ${this.enableLogging}`);
        console.log(`  Learning Rate:      ${this.learningRate}`);
    }
}

// ===========================================================================================
// MODULE EXPORTS
// ===========================================================================================

module.exports = {
    InertiaCouplingModule,
    PhysicsTerm,
    DynamicVacuumTerm,
    QuantumCouplingTerm
};

// ===========================================================================================
// USAGE EXAMPLES (Comment out when integrating)
// ===========================================================================================

/*
// Example 1: Basic usage
const module = new InertiaCouplingModule();
const sum = module.computeSumInertiaTerms(0.0);
console.log(`Sum Inertia Terms = ${sum.toExponential(4)} J/m³`);
module.printInertiaBreakdown(0.0);
console.log(module.getEquationText());

// Example 2: Dynamic variable updates
module.updateVariable('lambda', 1.1);
module.updateVariable('E_react', 2e46);
module.printVariables();

// Example 3: State management
const state = module.exportState();
const newModule = new InertiaCouplingModule();
newModule.importState(state);

// Example 4: Time evolution
for (let t = 0; t <= 1e8; t += 1e7) {
    const sum = module.computeSumInertiaTerms(t);
    console.log(`t=${t}s: Σ = ${sum.toExponential(4)} J/m³`);
}
*/
