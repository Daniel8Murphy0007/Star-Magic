// source108.js
// CorePenetrationModule - Planetary Core Penetration Factor (P_core) in UQFF Framework
// Converted from source108.cpp - Maintains all self-expanding dynamics
// 
// Physics: P_core ≈1 for Sun (unitless), ~1e-3 for planets
// Role: Scales P_core in Universal Gravity U_g3 term (magnetic strings disk)
// Formula: U_g3 = k_3 × μ_j B_j(r,θ,t,ρ_vac,[SCm]) × cos(ω_s(t) t π) × P_core × E_react
// Applications:
//   - Core penetration: Full for stellar plasma (P_core=1), reduced for solid cores (P_core~1e-3)
//   - Adjusts magnetic disk gravity for core [SCm] influence
//   - Models nebulae/star formation/planetary disks
// Example: Sun t=0: U_g3 ≈1.8e49 J/m³ (P_core=1); Planet: ≈1.8e46 J/m³ (P_core=1e-3)
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

class CorePenetrationModule {
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
        this.variables.set('P_core', 1.0);                      // Unitless ≈1 for Sun
        this.variables.set('k_3', 1.8);                         // Coupling constant
        this.variables.set('B_j', 1e3);                         // T (base magnetic field)
        this.variables.set('omega_s', 2.5e-6);                  // rad/s (solar rotation frequency)
        this.variables.set('P_core_planet', 1e-3);              // For planets (3 orders lower)
        this.variables.set('E_react', 1e46);                    // J (reaction energy)
        this.variables.set('pi', 3.141592653589793);
        this.variables.set('t', 0.0);                           // s (time)
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
            const currentValue = this.variables.get(name);
            this.variables.set(name, currentValue + delta);
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

    // Compute P_core ≈1 for Sun
    computeP_core() {
        return this.variables.get('P_core');
    }

    // Compute U_g3 with P_core: U_g3 = k_3 × B_j × cos(ω_s t π) × P_core × E_react
    computeU_g3(t) {
        this.variables.set('t', t);
        const k_3 = this.variables.get('k_3');
        const b_j = this.variables.get('B_j');
        const omega_s = this.variables.get('omega_s');
        const pi = this.variables.get('pi');
        const cos_term = Math.cos(omega_s * t * pi);
        const p_core = this.computeP_core();
        const e_react = this.variables.get('E_react');
        
        return k_3 * b_j * cos_term * p_core * e_react;
    }

    // U_g3 for planet (P_core=1e-3)
    computeU_g3_planet(t) {
        const orig_p = this.variables.get('P_core');
        this.variables.set('P_core', this.variables.get('P_core_planet'));
        const result = this.computeU_g3(t);
        this.variables.set('P_core', orig_p);
        return result;
    }

    // Compute scaling factor between stellar and planetary cores
    computeScalingFactor() {
        return this.variables.get('P_core') / this.variables.get('P_core_planet');
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
        return `U_g3 = k_3 × μ_j B_j(r,θ,t,ρ_vac,[SCm]) × cos(ω_s(t) t π) × P_core × E_react
Where P_core ≈1 (unitless for Sun, ~1e-3 for planets; core penetration).
Scales magnetic disk gravity for core [SCm] influence.
Example Sun t=0: U_g3 ≈1.8e49 J/m³ (P_core=1);
Planet: ≈1.8e46 J/m³ (P_core=1e-3, -3 orders).
Role: Adjusts core interactions; full for stellar plasma, reduced for solid cores.
UQFF: Models penetration in nebulae/star formation/disks.`;
    }

    printVariables() {
        console.log('\n=== CorePenetrationModule Variables ===');
        for (const [key, value] of this.variables.entries()) {
            console.log(`  ${key} = ${value.toExponential(6)}`);
        }
    }

    printCorePenetration(t) {
        console.log('\n=== Core Penetration Analysis ===');
        console.log(`Time: t = ${t.toExponential(3)} s`);
        
        console.log(`\n1. Core Penetration Factors:`);
        const p_core_sun = this.variables.get('P_core');
        const p_core_planet = this.variables.get('P_core_planet');
        console.log(`   P_core (Sun):    ${p_core_sun.toExponential(3)} (full plasma penetration)`);
        console.log(`   P_core (Planet): ${p_core_planet.toExponential(3)} (solid core, 3 orders lower)`);
        console.log(`   Scaling Factor:  ${this.computeScalingFactor().toExponential(3)}x`);

        console.log(`\n2. U_g3 Magnetic Strings Disk Energy:`);
        const u_g3_sun = this.computeU_g3(t);
        const u_g3_planet = this.computeU_g3_planet(t);
        console.log(`   U_g3 (Sun):      ${u_g3_sun.toExponential(4)} J/m³`);
        console.log(`   U_g3 (Planet):   ${u_g3_planet.toExponential(4)} J/m³`);
        console.log(`   Ratio:           ${(u_g3_sun / u_g3_planet).toExponential(3)}x`);

        console.log(`\n3. Physical Parameters:`);
        console.log(`   k_3:             ${this.variables.get('k_3').toExponential(3)} (coupling)`);
        console.log(`   B_j:             ${this.variables.get('B_j').toExponential(3)} T`);
        console.log(`   ω_s:             ${this.variables.get('omega_s').toExponential(3)} rad/s`);
        console.log(`   E_react:         ${this.variables.get('E_react').toExponential(3)} J`);
        
        const omega_s = this.variables.get('omega_s');
        const pi = this.variables.get('pi');
        const cos_term = Math.cos(omega_s * t * pi);
        console.log(`   cos(ω_s t π):    ${cos_term.toFixed(6)}`);

        if (this.dynamicTerms.length > 0) {
            const dynamic_contrib = this.computeDynamicTerms(t);
            console.log(`\n4. Dynamic Terms:`);
            console.log(`   Total dynamic contribution = ${dynamic_contrib.toExponential(4)}`);
            console.log(`   Number of terms = ${this.dynamicTerms.length}`);
        }
    }

    getModuleInfo() {
        console.log('\n=== CorePenetrationModule Info ===');
        console.log(`  Module Version:     ${this.metadata.get('version')}`);
        console.log(`  Enhanced:           ${this.metadata.get('enhanced')}`);
        console.log(`  Variables:          ${this.variables.size}`);
        console.log(`  P_core (Sun):       ${this.variables.get('P_core').toExponential(3)}`);
        console.log(`  P_core (Planet):    ${this.variables.get('P_core_planet').toExponential(3)}`);
        console.log(`  k_3:                ${this.variables.get('k_3').toExponential(3)}`);
        console.log(`  Dynamic Terms:      ${this.dynamicTerms.length}`);
        console.log(`  Dynamic Parameters: ${this.dynamicParameters.size}`);
        console.log(`  Logging Enabled:    ${this.enableLogging}`);
        console.log(`  Learning Rate:      ${this.learningRate}`);
    }
}

// Export for use in other modules
module.exports = {
    CorePenetrationModule,
    PhysicsTerm,
    DynamicVacuumTerm,
    QuantumCouplingTerm
};

// Example usage (commented out for module use):
/*
const { CorePenetrationModule } = require('./source108.js');

const mod = new CorePenetrationModule();

// Basic P_core computations
console.log('\n=== Core Penetration Factor Validation ===');
console.log(`P_core (Sun) = ${mod.computeP_core()}`);

// U_g3 calculations
console.log('\n=== U_g3 Magnetic Strings Disk Energy ===');
const u_g3_sun = mod.computeU_g3(0.0);
console.log(`U_g3 (Sun, t=0) = ${u_g3_sun.toExponential(4)} J/m³`);

const u_g3_planet = mod.computeU_g3_planet(0.0);
console.log(`U_g3 (Planet, t=0) = ${u_g3_planet.toExponential(4)} J/m³`);

console.log(`Scaling factor: ${(u_g3_sun / u_g3_planet).toExponential(3)}x`);

// Full core penetration analysis
mod.printCorePenetration(0.0);
mod.printCorePenetration(1e6);

// Dynamic terms
mod.registerDynamicTerm(new DynamicVacuumTerm());
mod.printCorePenetration(1e6);

// Equation text
console.log('\n' + mod.getEquationText());

// Module info
mod.getModuleInfo();
*/
