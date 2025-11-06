// source107.js
// PiConstantModule - Mathematical Constant Pi (π) in UQFF Framework
// Converted from source107.cpp - Maintains all self-expanding dynamics
// 
// Physics: π ≈ 3.141592653589793 (unitless mathematical constant)
// Role: Defines periodicity in oscillations; C=2π r; trig args (sin/cos with 2π cycle)
// Applications: 
//   - μ_j = (10³ + 0.4 sin(ω_c t)) × 3.38e20 Tï¿½m³, ω_c = 2π / period
//   - U_g1: ... cos(π t_n) ... (time-reversal oscillations)
//   - Solar cycles, rotations in nebulae/quasars
// Example: t=0, t_n=0: sin(ω_c t)=0 → μ_j=3.38e23 Tï¿½m³; cos(π t_n)=1
//
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025

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

class PiConstantModule {
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

        // Mathematical constants
        this.variables.set('pi', 3.141592653589793);           // Unitless
        this.variables.set('t_n', 0.0);                        // Negative time factor (days)
        this.variables.set('t', 0.0);                          // Time (s)
        this.variables.set('period', 3.96e8);                  // s (example solar cycle ~12.5 years)
        this.variables.set('omega_c', 2.0 * this.variables.get('pi') / 3.96e8);  // rad/s
        this.variables.set('base_mu', 3.38e20);                // Tï¿½m³
        this.variables.set('B_j', 1e3);                        // Base T
    }

    // ========== DYNAMIC VARIABLE OPERATIONS ==========

    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
            if (name === 'period') {
                this.variables.set('omega_c', 2.0 * this.variables.get('pi') / value);
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
            const currentValue = this.variables.get(name);
            this.variables.set(name, currentValue + delta);
            if (name === 'period') {
                this.variables.set('omega_c', 2.0 * this.variables.get('pi') / this.variables.get(name));
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

    // Compute π ≈ 3.141592653589793
    computePi() {
        return this.variables.get('pi');
    }

    // Compute cos(π t_n)
    computeCosPiTn(t_n) {
        this.variables.set('t_n', t_n);
        return Math.cos(this.computePi() * t_n);
    }

    // Compute sin(ω_c t), ω_c = 2π / period
    computeSinOmegaCT(t) {
        this.variables.set('t', t);
        return Math.sin(this.variables.get('omega_c') * t);
    }

    // Example μ_j = (10³ + 0.4 sin(ω_c t)) × 3.38e20 Tï¿½m³
    computeMuJExample(t) {
        const sin_omega = this.computeSinOmegaCT(t);
        const b_j = this.variables.get('B_j') + 0.4 * sin_omega;
        return b_j * this.variables.get('base_mu');
    }

    // Example cos(π t_n) in U_g1
    computeUg1CosTerm(t_n) {
        return this.computeCosPiTn(t_n);
    }

    // Compute tan(π/4) for validation (should be 1.0)
    computeTanPiOver4() {
        return Math.tan(this.variables.get('pi') / 4.0);
    }

    // Compute 2π (full circle in radians)
    computeTwoPi() {
        return 2.0 * this.variables.get('pi');
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
            console.log(`Set dynamic parameter ${name} = ${value}`);
        }
    }

    getDynamicParameter(name, defaultValue = 0.0) {
        return this.dynamicParameters.has(name) ? this.dynamicParameters.get(name) : defaultValue;
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

            if (this.enableLogging) {
                console.log('State imported successfully');
            }
            return true;
        } catch (error) {
            console.error('Failed to import state:', error.message);
            return false;
        }
    }

    // ========== OUTPUT AND DEBUGGING ==========

    getEquationText() {
        return `π ≈ 3.141592653589793 (unitless mathematical constant).
Role: Defines periodicity in oscillations; C=2π r; trig args (sin/cos with 2π cycle).
In U_m: μ_j = (10³ + 0.4 sin(ω_c t)) × 3.38e20; ω_c = 2π / period.
In U_g1: ... cos(π t_n) ... (time-reversal oscillations).
Example t=0, t_n=0: sin(ω_c t)=0 → μ_j=3.38e23 Tï¿½m³; cos(π t_n)=1.
UQFF: Ensures cyclic/TRZ dynamics; solar cycles, rotations in nebulae/quasars.`;
    }

    printVariables() {
        console.log('\n=== PiConstantModule Variables ===');
        for (const [key, value] of this.variables.entries()) {
            console.log(`  ${key} = ${value.toExponential(6)}`);
        }
    }

    printPiApplications(t, t_n) {
        console.log('\n=== Pi Applications in UQFF ===');
        console.log(`Time: t = ${t.toExponential(3)} s, t_n = ${t_n.toExponential(3)} days`);
        console.log(`\n1. Core Constant:`);
        console.log(`   π = ${this.computePi().toFixed(15)}`);
        console.log(`   2π = ${this.computeTwoPi().toFixed(15)}`);
        console.log(`   tan(π/4) = ${this.computeTanPiOver4().toFixed(6)} (validation: should be 1.0)`);

        console.log(`\n2. Oscillatory Functions:`);
        const cos_pi_tn = this.computeCosPiTn(t_n);
        console.log(`   cos(π t_n) = ${cos_pi_tn.toFixed(6)}`);
        const sin_omega = this.computeSinOmegaCT(t);
        console.log(`   sin(ω_c t) = ${sin_omega.toFixed(6)}`);
        console.log(`   ω_c = ${this.variables.get('omega_c').toExponential(6)} rad/s`);
        console.log(`   period = ${this.variables.get('period').toExponential(3)} s`);

        console.log(`\n3. Physical Applications:`);
        const mu_j = this.computeMuJExample(t);
        console.log(`   μ_j(t) = ${mu_j.toExponential(4)} Tï¿½m³`);
        console.log(`   B_j = ${this.variables.get('B_j').toExponential(3)} T`);
        
        const ug1_cos = this.computeUg1CosTerm(t_n);
        console.log(`   U_g1 cos term = ${ug1_cos.toFixed(6)}`);

        if (this.dynamicTerms.length > 0) {
            const dynamic_contrib = this.computeDynamicTerms(t);
            console.log(`\n4. Dynamic Terms:`);
            console.log(`   Total dynamic contribution = ${dynamic_contrib.toExponential(4)}`);
            console.log(`   Number of terms = ${this.dynamicTerms.length}`);
        }
    }

    getModuleInfo() {
        console.log('\n=== PiConstantModule Info ===');
        console.log(`  Module Version:     ${this.metadata.get('version')}`);
        console.log(`  Enhanced:           ${this.metadata.get('enhanced')}`);
        console.log(`  Variables:          ${this.variables.size}`);
        console.log(`  π:                  ${this.variables.get('pi').toFixed(15)}`);
        console.log(`  ω_c:                ${this.variables.get('omega_c').toExponential(6)} rad/s`);
        console.log(`  period:             ${this.variables.get('period').toExponential(3)} s`);
        console.log(`  Dynamic Terms:      ${this.dynamicTerms.length}`);
        console.log(`  Dynamic Parameters: ${this.dynamicParameters.size}`);
        console.log(`  Logging Enabled:    ${this.enableLogging}`);
        console.log(`  Learning Rate:      ${this.learningRate}`);
    }
}

// Export for use in other modules
module.exports = {
    PiConstantModule,
    PhysicsTerm,
    DynamicVacuumTerm,
    QuantumCouplingTerm
};

// Example usage (commented out for module use):
/*
const { PiConstantModule } = require('./source107.js');

const mod = new PiConstantModule();

// Basic π computations
console.log('\n=== Pi Constant Validation ===');
console.log(`π = ${mod.computePi()}`);
console.log(`2π = ${mod.computeTwoPi()}`);
console.log(`tan(π/4) = ${mod.computeTanPiOver4()}`);

// Oscillatory applications
console.log('\n=== Oscillatory Applications ===');
console.log(`cos(π × 0) = ${mod.computeCosPiTn(0.0)}`);
console.log(`cos(π × 1) = ${mod.computeCosPiTn(1.0)}`);
console.log(`sin(ω_c × 0) = ${mod.computeSinOmegaCT(0.0)}`);

// Physical example: μ_j
const mu_j = mod.computeMuJExample(0.0);
console.log(`\nμ_j (t=0) = ${mu_j.toExponential(4)} Tï¿½m³`);

// Full applications display
mod.printPiApplications(0.0, 0.0);
mod.printPiApplications(3.96e8, 1.0);

// Dynamic terms
mod.registerDynamicTerm(new DynamicVacuumTerm());
mod.printPiApplications(1e6, 0.0);

// Equation text
console.log('\n' + mod.getEquationText());

// Module info
mod.getModuleInfo();
*/
