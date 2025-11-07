// source117.js
// StellarMassModule - Stellar/Planetary Mass (M_s) in Universal Quantum Field Superconductive Framework (UQFF)
// JavaScript conversion from source117.cpp
// Computes M_s=1.989e30 kg (1 M_sun for Sun); scales M_s / r^2 in Universal Gravity U_g1 and U_g2 terms
// Variables in Map; example for Sun at r=1.496e13 m; U_g2 ≈1.18e53 J/m³
// Approximations: S(r - R_b)=1; (1 + δ_sw v_sw)=5001; H_SCm=1; E_react=1e46
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025
// Enhanced with self-expanding framework - November 2025

// ===========================================================================================
// SELF-EXPANDING FRAMEWORK: Dynamic Physics Term System
// ===========================================================================================

class PhysicsTerm {
    constructor() {
        this.dynamicParameters = new Map();
        this.metadata = new Map();
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
// ENHANCED CLASS WITH SELF-EXPANDING CAPABILITIES
// ===========================================================================================

class StellarMassModule {
    constructor() {
        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();

        // Universal constants
        this.variables.set('M_s', 1.989e30);                    // kg (Sun)
        this.variables.set('M_sun', 1.989e30);                  // kg
        this.variables.set('k_1', 1.5);                         // Coupling for U_g1
        this.variables.set('k_2', 1.2);                         // Coupling for U_g2
        this.variables.set('rho_vac_UA', 7.09e-36);             // J/m³
        this.variables.set('rho_vac_SCm', 7.09e-37);            // J/m³
        this.variables.set('r', 1.496e13);                      // m (example R_b)
        this.variables.set('R_b', 1.496e13);                    // m
        this.variables.set('S_r_Rb', 1.0);                      // Step
        this.variables.set('delta_sw', 0.01);                   // Unitless
        this.variables.set('v_sw', 5e5);                        // m/s
        this.variables.set('H_SCm', 1.0);                       // Unitless
        this.variables.set('E_react', 1e46);                    // J

        // Derived
        this.variables.set('rho_sum', this.variables.get('rho_vac_UA') + this.variables.get('rho_vac_SCm'));
        this.variables.set('swirl_factor', 1.0 + this.variables.get('delta_sw') * this.variables.get('v_sw'));

        // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;

        this.metadata.set('enhanced', 'true');
        this.metadata.set('version', '2.0-Enhanced');
        this.metadata.set('module', 'StellarMassModule');
        this.metadata.set('converted', 'November 2025');
    }

    // ========== VARIABLE MANAGEMENT ==========

    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
            this._updateDerivedVariables(name);
        } else {
            console.warn(`Variable '${name}' not found. Adding with value ${value}`);
            this.variables.set(name, value);
        }
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
            this._updateDerivedVariables(name);
        } else {
            console.warn(`Variable '${name}' not found. Adding with delta ${delta}`);
            this.variables.set(name, delta);
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    _updateDerivedVariables(name) {
        if (name === 'rho_vac_UA' || name === 'rho_vac_SCm') {
            this.variables.set('rho_sum', this.variables.get('rho_vac_UA') + this.variables.get('rho_vac_SCm'));
        } else if (name === 'delta_sw' || name === 'v_sw') {
            this.variables.set('swirl_factor', 1.0 + this.variables.get('delta_sw') * this.variables.get('v_sw'));
        }
    }

    // ========== CORE COMPUTATIONS ==========

    // Compute M_s (kg)
    computeM_s() {
        return this.variables.get('M_s');
    }

    // M_s in M_sun
    computeM_sInMsun() {
        return this.computeM_s() / this.variables.get('M_sun');
    }

    // M_s / r² (kg/m²)
    computeM_sOverR2(r) {
        this.variables.set('r', r);
        if (r === 0.0) return 0.0;
        return this.computeM_s() / (r * r);
    }

    // U_g1 example (internal, simplified)
    computeU_g1(r) {
        const k_1 = this.variables.get('k_1');
        const rho_sum = this.variables.get('rho_sum');
        const m_over_r2 = this.computeM_sOverR2(r);
        const e_react = this.variables.get('E_react');
        return k_1 * rho_sum * m_over_r2 * e_react;  // Simplified
    }

    // U_g2 example (outer bubble)
    computeU_g2(r) {
        this.variables.set('r', r);
        const k_2 = this.variables.get('k_2');
        const rho_sum = this.variables.get('rho_sum');
        const R_b = this.variables.get('R_b');
        const s_step = (r >= R_b) ? 1.0 : 0.0;
        const swirl_factor = this.variables.get('swirl_factor');
        const h_scm = this.variables.get('H_SCm');
        const e_react = this.variables.get('E_react');
        return k_2 * rho_sum * this.computeM_sOverR2(r) * s_step * swirl_factor * h_scm * e_react;
    }

    // Compute U_g1 and U_g2 ratio (comparison between internal vs external gravity)
    computeGravityRatio(r) {
        const u_g1 = this.computeU_g1(r);
        const u_g2 = this.computeU_g2(r);
        if (u_g2 === 0.0) return 0.0;
        return u_g1 / u_g2;
    }

    // Compute mass scaling effects at different masses
    computeMassScaling(mass_factor, r) {
        const original_M_s = this.variables.get('M_s');
        const original_U_g2 = this.computeU_g2(r);
        
        this.updateVariable('M_s', original_M_s * mass_factor);
        const scaled_U_g2 = this.computeU_g2(r);
        
        // Restore original
        this.updateVariable('M_s', original_M_s);
        
        return {
            mass_factor: mass_factor,
            original_M_s: original_M_s,
            scaled_M_s: original_M_s * mass_factor,
            original_U_g2: original_U_g2,
            scaled_U_g2: scaled_U_g2,
            scaling_ratio: scaled_U_g2 / original_U_g2
        };
    }

    // Print stellar mass effects at a given radius
    printStellarMassEffects(r) {
        console.log('\n=== Stellar Mass Effects ===');
        console.log(`M_s: ${this.computeM_s().toExponential(3)} kg`);
        console.log(`M_s: ${this.computeM_sInMsun().toFixed(2)} M_sun`);
        console.log(`r: ${r.toExponential(3)} m`);
        console.log(`M_s / r²: ${this.computeM_sOverR2(r).toExponential(3)} kg/m²`);
        console.log(`U_g1 (internal): ${this.computeU_g1(r).toExponential(3)} J/m³`);
        console.log(`U_g2 (external): ${this.computeU_g2(r).toExponential(3)} J/m³`);
        console.log(`U_g1/U_g2 ratio: ${this.computeGravityRatio(r).toFixed(3)}`);
        console.log('============================\n');
    }

    // ========== SELF-EXPANDING FRAMEWORK METHODS ==========

    registerDynamicTerm(term) {
        if (term instanceof PhysicsTerm) {
            this.dynamicTerms.push(term);
            if (this.enableLogging) {
                console.log(`Registered dynamic term: ${term.getName()}`);
            }
        }
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

    computeDynamicTerms(t) {
        if (!this.enableDynamicTerms) return 0.0;

        let total = 0.0;
        for (const term of this.dynamicTerms) {
            if (term.validate(this.variables)) {
                total += term.compute(t, this.variables);
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

    importState(stateJson) {
        const state = JSON.parse(stateJson);
        this.variables = new Map(state.variables);
        this.dynamicParameters = new Map(state.dynamicParameters);
        this.metadata = new Map(state.metadata);
        this.enableDynamicTerms = state.enableDynamicTerms;
        this.enableLogging = state.enableLogging;
        this.learningRate = state.learningRate;
    }

    setEnableLogging(enable) {
        this.enableLogging = enable;
    }

    setLearningRate(rate) {
        this.learningRate = rate;
    }

    // ========== OUTPUT METHODS ==========

    getEquationText() {
        return `U_g1 = k_1 * ρ_vac,[UA/SCm] * (M_s / r²) * ... E_react (internal dipole);
U_g2 = k_2 * ρ_vac,[UA/SCm] * (M_s / r²) * S(r - R_b) * (1 + δ_sw v_sw) * H_SCm * E_react (outer bubble).
Where M_s = 1.989e30 kg (1 M_sun for Sun).
Scales gravity by mass; M_s / r² ≈8.89e3 kg/m² at r=1.496e13 m.
Example U_g2 (r=R_b): ≈1.18e53 J/m³.
Role: Central mass drives internal/external gravity; stellar/planetary dynamics.
UQFF: Mass-dependent fields for nebulae/formation/mergers.`;
    }

    printVariables() {
        console.log('\nCurrent Variables:');
        for (const [key, value] of this.variables.entries()) {
            console.log(`${key} = ${value.toExponential(3)}`);
        }
    }

    getModuleInfo() {
        return {
            module: 'StellarMassModule',
            version: this.metadata.get('version'),
            enhanced: this.metadata.get('enhanced'),
            variableCount: this.variables.size,
            M_s_kg: this.computeM_s(),
            M_s_Msun: this.computeM_sInMsun(),
            dynamicTerms: this.dynamicTerms.length,
            dynamicParameters: this.dynamicParameters.size,
            loggingEnabled: this.enableLogging,
            learningRate: this.learningRate
        };
    }
}

// ===========================================================================================
// MODULE EXPORTS
// ===========================================================================================

module.exports = {
    StellarMassModule,
    PhysicsTerm,
    DynamicVacuumTerm,
    QuantumCouplingTerm
};
