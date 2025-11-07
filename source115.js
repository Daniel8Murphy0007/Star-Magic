// source115.js
// SolarWindModulationModule - Solar Wind Modulation Factor (δ_sw) UQFF Module
// JavaScript conversion from Source115.cpp
// Maintains self-expanding framework with dynamic physics terms and state management
//
// Core Physics:
// - δ_sw = 0.01 (unitless solar wind modulation factor)
// - Modulation Factor: 1 + δ_sw × v_sw
// - v_sw = 5×10⁵ m/s (solar wind velocity)
// - Amplification: ~5001× at standard parameters
//
// U_g2 Computation:
// U_g2 = k_2 × [(ρ_vac,UA + ρ_vac,SCm) × M_s / r²] × S(r - R_b) × (1 + δ_sw × v_sw) × H_SCm × E_react
//
// Where:
// - k_2 = 1.2 (coupling constant)
// - ρ_vac,UA = 7.09×10⁻³⁶ J/m³ (Universal Ambiance vacuum density)
// - ρ_vac,SCm = 7.09×10⁻³⁷ J/m³ (Superconductive Medium vacuum density)
// - M_s = 1.989×10³⁰ kg (solar mass)
// - r = distance from Sun (m)
// - R_b = 1.496×10¹³ m (100 AU, heliopause boundary)
// - S(r - R_b) = step function (1 if r ≥ R_b, else 0)
// - H_SCm = 1.0 (Superconductive Medium factor)
// - E_react = 10⁴⁶ J (reactive energy)
//
// Example at r = R_b = 1.496×10¹³ m:
// - With modulation: U_g2 ≈ 1.18×10⁵³ J/m³
// - Without modulation (δ_sw=0): U_g2 ≈ 2.36×10⁴⁹ J/m³
// - Amplification ratio: ~5000×
//
// Physical Significance:
// - Models heliosphere dynamics beyond heliopause
// - Solar wind momentum/pressure enhances external gravity
// - Critical for nebular dynamics, star formation, galactic interactions
// - Explains wind influence on interstellar medium
//
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025
// Enhanced: JavaScript conversion with self-expanding framework

// ===========================================================================================
// SELF-EXPANDING FRAMEWORK: Dynamic Physics Term System
// ===========================================================================================

class PhysicsTerm {
    constructor() {
        this.dynamicParameters = new Map();
        this.metadata = new Map();
    }

    compute(t, params) {
        throw new Error("PhysicsTerm.compute() must be implemented by subclass");
    }

    getName() {
        throw new Error("PhysicsTerm.getName() must be implemented by subclass");
    }

    getDescription() {
        throw new Error("PhysicsTerm.getDescription() must be implemented by subclass");
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
        return "Time-varying vacuum energy contribution to solar wind modulation";
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
        return "Non-local quantum effects in heliosphere";
    }
}

// ===========================================================================================
// ENHANCED CLASS WITH SELF-EXPANDING CAPABILITIES
// ===========================================================================================

class SolarWindModulationModule {
    constructor() {
        // Core variables (Map-based for dynamic management)
        this.variables = new Map();

        // Universal constants
        this.variables.set('delta_sw', 0.01);           // Unitless modulation factor
        this.variables.set('v_sw', 5e5);                // m/s (solar wind velocity)
        this.variables.set('k_2', 1.2);                 // Coupling constant
        this.variables.set('rho_vac_UA', 7.09e-36);     // J/m³ (UA vacuum density)
        this.variables.set('rho_vac_SCm', 7.09e-37);    // J/m³ (SCm vacuum density)
        this.variables.set('M_s', 1.989e30);            // kg (solar mass)
        this.variables.set('r', 1.496e13);              // m (distance, default = R_b)
        this.variables.set('R_b', 1.496e13);            // m (heliopause at 100 AU)
        this.variables.set('S_r_Rb', 1.0);              // Step function value
        this.variables.set('H_SCm', 1.0);               // Superconductive Medium factor
        this.variables.set('E_react', 1e46);            // J (reactive energy)

        // Derived variables
        this.variables.set('rho_sum', 
            this.variables.get('rho_vac_UA') + this.variables.get('rho_vac_SCm'));
        this.variables.set('modulation_factor', this.computeModulationFactor());

        // Self-expanding framework
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;

        this.metadata.set('enhanced', 'true');
        this.metadata.set('version', '2.0-Enhanced');
        this.metadata.set('module', 'SolarWindModulationModule');
        this.metadata.set('description', 'Solar wind modulation factor δ_sw and U_g2 computation');
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

    _updateDerivedVariables(changedName) {
        // Update modulation_factor if delta_sw or v_sw changes
        if (changedName === 'delta_sw' || changedName === 'v_sw') {
            this.variables.set('modulation_factor', this.computeModulationFactor());
        }
        // Update rho_sum if vacuum densities change
        else if (changedName === 'rho_vac_UA' || changedName === 'rho_vac_SCm') {
            this.variables.set('rho_sum', 
                this.variables.get('rho_vac_UA') + this.variables.get('rho_vac_SCm'));
        }
    }

    // ========== CORE COMPUTATIONS ==========

    computeDelta_sw() {
        return this.variables.get('delta_sw');
    }

    computeModulationFactor() {
        const delta_sw = this.variables.get('delta_sw');
        const v_sw = this.variables.get('v_sw');
        return 1.0 + delta_sw * v_sw;
    }

    computeU_g2(r) {
        this.variables.set('r', r);
        
        const k_2 = this.variables.get('k_2');
        const rho_sum = this.variables.get('rho_sum');
        const M_s = this.variables.get('M_s');
        const R_b = this.variables.get('R_b');
        const H_SCm = this.variables.get('H_SCm');
        const E_react = this.variables.get('E_react');
        
        // Step function: S(r - R_b) = 1 if r >= R_b, else 0
        const s_step = (r >= R_b) ? 1.0 : 0.0;
        const mod_factor = this.computeModulationFactor();
        
        const U_g2 = k_2 * (rho_sum * M_s / (r * r)) * s_step * mod_factor * H_SCm * E_react;
        
        // Add dynamic terms if enabled
        let dynamicContribution = 0;
        if (this.enableDynamicTerms && this.dynamicTerms.length > 0) {
            const params = new Map([
                ['r', r],
                ['M_s', M_s],
                ['rho_vac_UA', this.variables.get('rho_vac_UA')],
                ['delta_sw', this.variables.get('delta_sw')],
                ['v_sw', this.variables.get('v_sw')]
            ]);
            
            for (const term of this.dynamicTerms) {
                if (term.validate(params)) {
                    dynamicContribution += term.compute(0, params);
                }
            }
        }
        
        return U_g2 + dynamicContribution;
    }

    computeU_g2_no_mod(r) {
        const original_delta = this.variables.get('delta_sw');
        this.variables.set('delta_sw', 0.0);
        this.variables.set('modulation_factor', this.computeModulationFactor());
        
        const result = this.computeU_g2(r);
        
        this.variables.set('delta_sw', original_delta);
        this.variables.set('modulation_factor', this.computeModulationFactor());
        
        return result;
    }

    computeAmplificationRatio(r) {
        const with_mod = this.computeU_g2(r);
        const without_mod = this.computeU_g2_no_mod(r);
        return without_mod !== 0 ? with_mod / without_mod : 0;
    }

    computeModulationPercentage() {
        const mod_factor = this.computeModulationFactor();
        return ((mod_factor - 1.0) / 1.0) * 100;
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

    setEnableLogging(enabled) {
        this.enableLogging = enabled;
    }

    setLearningRate(rate) {
        this.learningRate = rate;
    }

    exportState() {
        const state = {
            variables: Object.fromEntries(this.variables),
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            metadata: Object.fromEntries(this.metadata),
            settings: {
                enableDynamicTerms: this.enableDynamicTerms,
                enableLogging: this.enableLogging,
                learningRate: this.learningRate
            },
            dynamicTermCount: this.dynamicTerms.length
        };
        return JSON.stringify(state, null, 2);
    }

    importState(stateJson) {
        const state = JSON.parse(stateJson);
        
        this.variables = new Map(Object.entries(state.variables));
        this.dynamicParameters = new Map(Object.entries(state.dynamicParameters || {}));
        this.metadata = new Map(Object.entries(state.metadata || {}));
        
        if (state.settings) {
            this.enableDynamicTerms = state.settings.enableDynamicTerms;
            this.enableLogging = state.settings.enableLogging;
            this.learningRate = state.settings.learningRate;
        }
    }

    // ========== OUTPUT AND DEBUGGING ==========

    getEquationText() {
        return `U_g2 = k_2 × [(ρ_vac,UA + ρ_vac,SCm) × M_s / r²] × S(r - R_b) × (1 + δ_sw × v_sw) × H_SCm × E_react

Where δ_sw = ${this.variables.get('delta_sw')} (unitless solar wind modulation factor)
Modulation Factor = 1 + δ_sw × v_sw = ${this.computeModulationFactor().toExponential(3)}
v_sw = ${this.variables.get('v_sw').toExponential(2)} m/s (solar wind velocity)

Example at r = R_b = ${this.variables.get('R_b').toExponential(3)} m:
- With modulation: U_g2 ≈ ${this.computeU_g2(this.variables.get('R_b')).toExponential(3)} J/m³
- Without modulation (δ_sw=0): U_g2 ≈ ${this.computeU_g2_no_mod(this.variables.get('R_b')).toExponential(3)} J/m³
- Amplification ratio: ~${this.computeAmplificationRatio(this.variables.get('R_b')).toFixed(0)}×

Role: Enhances external gravity via solar wind momentum/pressure beyond heliopause (r ≥ R_b)
UQFF: Models heliosphere dynamics; wind influence on nebular/star formation, galactic interactions`;
    }

    printVariables() {
        console.log("\n=== SolarWindModulationModule Variables ===");
        for (const [name, value] of this.variables) {
            console.log(`  ${name} = ${typeof value === 'number' ? value.toExponential(3) : value}`);
        }
        console.log(`\n  Derived:`);
        console.log(`  Modulation Factor = ${this.computeModulationFactor().toExponential(3)}`);
        console.log(`  Modulation % = ${this.computeModulationPercentage().toFixed(1)}%`);
    }

    printSolarWindEffects(r) {
        const delta_sw = this.computeDelta_sw();
        const mod_factor = this.computeModulationFactor();
        const U_g2_with = this.computeU_g2(r);
        const U_g2_without = this.computeU_g2_no_mod(r);
        const amplification = this.computeAmplificationRatio(r);
        const mod_percentage = this.computeModulationPercentage();

        console.log("\n=== Solar Wind Modulation Effects ===");
        console.log(`  Distance r:              ${r.toExponential(3)} m`);
        console.log(`  Heliopause R_b:          ${this.variables.get('R_b').toExponential(3)} m`);
        console.log(`  Solar wind velocity:     ${this.variables.get('v_sw').toExponential(2)} m/s`);
        console.log(`  δ_sw (modulation):       ${delta_sw}`);
        console.log(`  Modulation factor:       ${mod_factor.toExponential(3)}`);
        console.log(`  Modulation percentage:   ${mod_percentage.toFixed(1)}%`);
        console.log(`  U_g2 (with modulation):  ${U_g2_with.toExponential(3)} J/m³`);
        console.log(`  U_g2 (no modulation):    ${U_g2_without.toExponential(3)} J/m³`);
        console.log(`  Amplification ratio:     ${amplification.toFixed(0)}×`);
        console.log(`  Step function S(r-R_b):  ${(r >= this.variables.get('R_b')) ? 1 : 0}`);
    }
}

// Export for use in other modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        SolarWindModulationModule,
        PhysicsTerm,
        DynamicVacuumTerm,
        QuantumCouplingTerm
    };
}
