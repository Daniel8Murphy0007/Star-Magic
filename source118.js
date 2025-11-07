// StellarRotationModule.js
// JavaScript implementation of the Stellar/Planetary Rotation Rate (ω_s) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ω_s=2.5e-6 rad/s (~29-day Sun period); scales ω_s(t) in U_g3 cos(ω_s t π) and U_i ω_s cos(π t_n).
// Variables in Map; example for Sun at t=0, t_n=0; U_g3 ≈1.8e49 J/m³, U_i ≈1.38e-47 J/m³.
// Approximations: cos(π t_n)=1; f_TRZ=0.1; λ_i=1.0; ρ_vac sum=7.80e-36 J/m³.
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
        throw new Error("compute() must be implemented by subclass");
    }

    getName() {
        throw new Error("getName() must be implemented by subclass");
    }

    getDescription() {
        throw new Error("getDescription() must be implemented by subclass");
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
        return "DynamicVacuum";
    }

    getDescription() {
        return "Time-varying vacuum energy";
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
        return "QuantumCoupling";
    }

    getDescription() {
        return "Non-local quantum effects";
    }
}

// ===========================================================================================
// ENHANCED CLASS WITH SELF-EXPANDING CAPABILITIES
// ===========================================================================================

class StellarRotationModule {
    constructor() {
        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();
        
        // Universal constants
        this.variables.set('omega_s', 2.5e-6);              // rad/s
        this.variables.set('k_3', 1.8);                     // Coupling U_g3
        this.variables.set('B_j', 1e3);                     // T
        this.variables.set('P_core', 1.0);                  // Unitless
        this.variables.set('E_react', 1e46);                // J
        this.variables.set('lambda_i', 1.0);                // Unitless U_i
        this.variables.set('rho_vac_SCm', 7.09e-37);        // J/m^3
        this.variables.set('rho_vac_UA', 7.09e-36);         // J/m^3
        this.variables.set('f_TRZ', 0.1);                   // Unitless
        this.variables.set('pi', Math.PI);
        this.variables.set('t', 0.0);                       // s
        this.variables.set('t_n', 0.0);                     // s
        
        // Derived
        this.variables.set('rho_sum', 
            this.variables.get('rho_vac_SCm') + this.variables.get('rho_vac_UA'));
        this.variables.set('day_to_s', 86400.0);            // s/day
        
        // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;
        
        this.metadata.set('enhanced', 'true');
        this.metadata.set('version', '2.0-Enhanced');
        this.metadata.set('module', 'StellarRotationModule');
        this.metadata.set('parameter', 'omega_s');
        this.metadata.set('description', 'Stellar/Planetary Rotation Rate (ω_s = 2.5e-6 rad/s)');
    }

    // ========== VARIABLE MANAGEMENT ==========
    
    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
            if (name === 'rho_vac_SCm' || name === 'rho_vac_UA') {
                this.variables.set('rho_sum', 
                    this.variables.get('rho_vac_SCm') + this.variables.get('rho_vac_UA'));
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
            this.variables.set(name, this.variables.get(name) + delta);
            if (name === 'rho_vac_SCm' || name === 'rho_vac_UA') {
                this.variables.set('rho_sum', 
                    this.variables.get('rho_vac_SCm') + this.variables.get('rho_vac_UA'));
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
    
    computeOmega_s() {
        return this.variables.get('omega_s');
    }

    computeOmega_s_t(t) {
        this.variables.set('t', t);
        return this.computeOmega_s();  // Constant for Sun (no time dependence)
    }

    computePeriod_days() {
        const period_s = 2.0 * Math.PI / this.computeOmega_s();
        return period_s / this.variables.get('day_to_s');
    }

    computeU_g3(t) {
        const k_3 = this.variables.get('k_3');
        const b_j = this.variables.get('B_j');
        const omega_s = this.variables.get('omega_s');
        const pi = this.variables.get('pi');
        const cos_term = Math.cos(omega_s * t * pi);
        const p_core = this.variables.get('P_core');
        const e_react = this.variables.get('E_react');
        
        let result = k_3 * b_j * cos_term * p_core * e_react;
        
        // Add dynamic terms if enabled
        if (this.enableDynamicTerms && this.dynamicTerms.length > 0) {
            for (const term of this.dynamicTerms) {
                result += term.compute(t, this.variables);
            }
        }
        
        return result;
    }

    computeU_i(t, t_n) {
        const lambda_i = this.variables.get('lambda_i');
        const rho_sc = this.variables.get('rho_vac_SCm');
        const rho_ua = this.variables.get('rho_vac_UA');
        const omega_s_t = this.computeOmega_s_t(t);
        const pi = this.variables.get('pi');
        const cos_pi_tn = Math.cos(pi * t_n);
        const trz_factor = 1.0 + this.variables.get('f_TRZ');
        
        return lambda_i * rho_sc * rho_ua * omega_s_t * cos_pi_tn * trz_factor;
    }

    computeRotationScaling(omega_factor, t) {
        const original_omega = this.variables.get('omega_s');
        this.variables.set('omega_s', original_omega * omega_factor);
        
        const U_g3 = this.computeU_g3(t);
        const U_i = this.computeU_i(t, 0.0);
        const period_days = this.computePeriod_days();
        
        // Restore original
        this.variables.set('omega_s', original_omega);
        
        return {
            omega_factor,
            U_g3,
            U_i,
            period_days
        };
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

    getDynamicParameter(name) {
        return this.dynamicParameters.get(name);
    }

    exportState() {
        return {
            variables: Object.fromEntries(this.variables),
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            metadata: Object.fromEntries(this.metadata),
            enableDynamicTerms: this.enableDynamicTerms,
            dynamicTermCount: this.dynamicTerms.length
        };
    }

    setEnableLogging(enabled) {
        this.enableLogging = enabled;
    }

    setLearningRate(rate) {
        this.learningRate = rate;
    }

    // ========== OUTPUT AND DIAGNOSTICS ==========
    
    getEquationText() {
        return `U_g3 = k_3 * B_j * cos(ω_s(t) t π) * P_core * E_react
U_i = λ_i * ρ_vac,[SCm] * ρ_vac,[UA] * ω_s(t) * cos(π t_n) * (1 + f_TRZ)
Where ω_s = 2.5e-6 rad/s (~29-day Sun equatorial rotation);
Scales rotational oscillations/inertia.
Example t=0, t_n=0: U_g3 ≈1.8e49 J/m³; U_i ≈1.38e-47 J/m³.
Role: Introduces spin in disk gravity/inertia; stellar/planetary dynamics.
UQFF: Rotational effects in nebulae/disks/formation/mergers.`;
    }

    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential(3)}`);
        }
    }

    printRotationEffects(t = 0.0, t_n = 0.0) {
        console.log("\n=== Stellar/Planetary Rotation Rate (ω_s) Effects ===");
        console.log(`ω_s = ${this.computeOmega_s().toExponential(3)} rad/s`);
        console.log(`Period = ${this.computePeriod_days().toFixed(2)} days`);
        console.log(`\nAt t = ${t.toExponential(2)} s, t_n = ${t_n}:`);
        console.log(`U_g3 = ${this.computeU_g3(t).toExponential(3)} J/m³`);
        console.log(`U_i = ${this.computeU_i(t, t_n).toExponential(3)} J/m³`);
        
        // Show rotation scaling
        console.log("\n=== Rotation Rate Scaling ===");
        const scalings = [
            { name: "Slow rotator (0.5×)", factor: 0.5 },
            { name: "Solar (1.0×)", factor: 1.0 },
            { name: "Fast rotator (2.0×)", factor: 2.0 }
        ];
        
        for (const scaling of scalings) {
            const result = this.computeRotationScaling(scaling.factor, t);
            console.log(`${scaling.name}: Period = ${result.period_days.toFixed(2)} days, U_g3 = ${result.U_g3.toExponential(2)} J/m³`);
        }
    }
}

// ===========================================================================================
// EXPORTS
// ===========================================================================================

module.exports = {
    StellarRotationModule,
    PhysicsTerm,
    DynamicVacuumTerm,
    QuantumCouplingTerm
};

// Example usage:
// const { StellarRotationModule } = require('./source118.js');
// const mod = new StellarRotationModule();
// const omega = mod.computeOmega_s();
// console.log(`ω_s = ${omega} rad/s (~${mod.computePeriod_days()} days)`);
// const u_g3 = mod.computeU_g3(0.0);
// console.log(`U_g3 = ${u_g3} J/m³`);
// mod.printRotationEffects();
