// source104.js - MagneticMomentModule
// JavaScript conversion of MagneticMomentModule.cpp for Star-Magic UQFF Framework
// Computes μ_j = (10³ + 0.4 sin(ω_c t)) * 3.38e20 T·m³; scales μ_j / r_j in Universal Magnetism U_m and Ug3
// Variables: base_mu, ω_c, r_j, γ, t_n, φ_hat_j, P_SCm, E_react, f_Heaviside, f_quasi, k3
// Example: j=1 at t=0: μ_j ≈3.38e23 T·m³; U_m contrib ≈2.28e65 J/m³; Ug3 ≈1.8e49 J/m³
// Role: Quantifies string magnetic strength; drives Um/Ug3 for jets/disks/nebulae
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

class MagneticMomentModule {
    constructor() {
        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();
        
        // Universal constants
        this.variables.set('base_mu', 3.38e20);                 // T·m³ (definition); example uses 3.38e23
        this.variables.set('omega_c', 2.5e-6);                  // rad/s
        this.variables.set('r_j', 1.496e13);                    // m (for j=1, AU distance)
        this.variables.set('gamma', 5e-5 / 86400.0);            // s⁻¹ (0.00005 day⁻¹)
        this.variables.set('t_n', 0.0);                         // s
        this.variables.set('phi_hat_j', 1.0);                   // Normalized
        this.variables.set('P_SCm', 1.0);                       // Pressure
        this.variables.set('E_react', 1e46);                    // J
        this.variables.set('f_Heaviside', 0.01);                // Unitless
        this.variables.set('f_quasi', 0.01);                    // Unitless
        this.variables.set('k3', 1.8);                          // Coupling for Ug3
        this.variables.set('pi', Math.PI);
        
        // Derived defaults
        this.variables.set('B_j', 1e3);                         // Base T
        this.variables.set('scale_Heaviside', 1e13);            // Amplification
        
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
        this.variables.set(name, value);
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
    
    // Compute μ_j(t) - Magnetic moment of j-th string (T·m³)
    computeMu_j(j, t) {
        const sin_term = Math.sin(this.variables.get('omega_c') * t);
        const b_j = this.variables.get('B_j') + 0.4 * sin_term;  // T
        return b_j * this.variables.get('base_mu');  // T·m³
    }

    // Compute B_j(t) - Base magnetic field (T)
    computeB_j(t) {
        return this.variables.get('B_j') + 0.4 * Math.sin(this.variables.get('omega_c') * t);
    }

    // Example U_m contribution for j-th string (J/m³, simplified)
    computeUmContrib(j, t) {
        const mu_j = this.computeMu_j(j, t);
        const r_j = this.variables.get('r_j');
        const gamma = this.variables.get('gamma');
        const t_n = this.variables.get('t_n');
        const pi = this.variables.get('pi');
        
        const exp_arg = -gamma * t * Math.cos(pi * t_n);
        const one_minus_exp = 1.0 - Math.exp(exp_arg);
        const phi_hat = this.variables.get('phi_hat_j');
        const heaviside_f = 1.0 + this.variables.get('scale_Heaviside') * this.variables.get('f_Heaviside');
        const quasi_f = 1.0 + this.variables.get('f_quasi');
        
        return (mu_j / r_j * one_minus_exp * phi_hat) * 
               this.variables.get('P_SCm') * 
               this.variables.get('E_react') * 
               heaviside_f * 
               quasi_f;
    }

    // Example Ug3 contribution (J/m³)
    computeUg3Contrib(t) {
        const b_j = this.computeB_j(t);
        const omega_c = this.variables.get('omega_c');
        const pi = this.variables.get('pi');
        const cos_term = Math.cos(omega_c * t * pi);  // Approximation
        const p_core = 1.0;
        const e_react = this.variables.get('E_react');
        
        return this.variables.get('k3') * b_j * cos_term * p_core * e_react;
    }

    // Get string label
    getStringLabel(j) {
        const labels = {
            1: 'Magnetic String 1 (Primary)',
            2: 'Magnetic String 2 (Secondary)',
            3: 'Magnetic String 3 (Tertiary)',
            4: 'Magnetic String 4 (Quaternary)'
        };
        return labels[j] || `String j=${j}`;
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
        return `μ_j = (10³ + 0.4 sin(ω_c t)) * 3.38e20 T·m³
Where ω_c=2.5e-6 rad/s; units T·m³ (magnetic dipole strength).
In U_m: μ_j [μ_j / r_j * (1 - e^{-γ t cos(π t_n)}) φ_hat_j ] P_SCm E_react (1 + 10¹³ f_Heaviside) (1 + f_quasi)
In Ug3: k3 * μ_j B_j cos(ω_s t π) P_core E_react; B_j = 10³ + 0.4 sin(ω_c t) T.
Example j=1, t=0: μ_j ≈3.38e23 T·m³; U_m contrib ≈2.28e65 J/m³; Ug3 ≈1.8e49 J/m³.
Role: Quantifies string magnetic strength; drives Um/Ug3 for jets/disks/nebulae.`;
    }

    printVariables() {
        console.log('\n=== MagneticMomentModule Variables ===');
        for (const [name, value] of this.variables.entries()) {
            console.log(`  ${name} = ${value.toExponential(4)}`);
        }
    }

    printMomentContributions(j = 1, t = 0.0) {
        const mu = this.computeMu_j(j, t);
        const b = this.computeB_j(t);
        const um = this.computeUmContrib(j, t);
        const ug3 = this.computeUg3Contrib(t);
        
        console.log(`\n=== Magnetic Moment ${this.getStringLabel(j)} at t=${t} s ===`);
        console.log(`  μ_j = ${mu.toExponential(4)} T·m³`);
        console.log(`  B_j = ${b.toExponential(4)} T`);
        console.log(`  U_m contrib = ${um.toExponential(4)} J/m³`);
        console.log(`  Ug3 contrib = ${ug3.toExponential(4)} J/m³`);
    }

    printTimeEvolution(j = 1, times = [0, 1e6, 1e7, 1e8]) {
        console.log(`\n=== Time Evolution for ${this.getStringLabel(j)} ===`);
        for (const t of times) {
            const mu = this.computeMu_j(j, t);
            const b = this.computeB_j(t);
            console.log(`  t=${t.toExponential(1)}s: μ_j=${mu.toExponential(4)} T·m³, B_j=${b.toExponential(4)} T`);
        }
    }

    printModuleInfo() {
        console.log('\n=== MagneticMomentModule Info ===');
        console.log(`  Module Version:     ${this.metadata.get('version')}`);
        console.log(`  Enhanced:           ${this.metadata.get('enhanced')}`);
        console.log(`  Variables:          ${this.variables.size}`);
        console.log(`  base_mu:            ${this.variables.get('base_mu').toExponential(4)} T·m³`);
        console.log(`  ω_c:                ${this.variables.get('omega_c').toExponential(4)} rad/s`);
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
    MagneticMomentModule,
    PhysicsTerm,
    DynamicVacuumTerm,
    QuantumCouplingTerm
};

// ===========================================================================================
// USAGE EXAMPLES (Comment out when integrating)
// ===========================================================================================

/*
// Example 1: Basic usage
const module = new MagneticMomentModule();
const t = 0.0;
module.printMomentContributions(1, t);
console.log(module.getEquationText());

// Example 2: Dynamic variable updates
module.updateVariable('base_mu', 4e20);
module.updateVariable('B_j', 2e3);
module.printVariables();

// Example 3: Time evolution
module.printTimeEvolution(1, [0, 1e6, 5e6, 1e7]);

// Example 4: State management
const state = module.exportState();
const newModule = new MagneticMomentModule();
newModule.importState(state);
*/
