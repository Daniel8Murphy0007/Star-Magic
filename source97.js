// Source97UQFFModule.js
// JavaScript implementation of the Feedback Factor Module in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes f_feedback=0.1 for ΔM_BH=1 dex (10x mass increase); scales (1 + f_feedback) in U_g4 term.
// Dynamic variables in Map; supports updateVariable, addToVariable, subtractFromVariable operations.
// Self-expanding framework with PhysicsTerm classes, dynamic parameters, metadata, and state management.
// Example: const mod = new Source97UQFFModule(); mod.computeU_g4(0.0, 0.0); mod.updateVariable("f_feedback", new_value);
// Watermark: Copyright - Daniel T. Murphy, analyzed Nov 6, 2025.

class PhysicsTerm {
    constructor() {
        this.dynamicParameters = new Map();
        this.metadata = new Map();
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

    setDynamicParameter(name, value) {
        this.dynamicParameters.set(name, value);
        if (this.enableLogging) {
            console.log(`PhysicsTerm ${this.getName()}: Set ${name} = ${value}`);
        }
    }

    getDynamicParameter(name, defaultValue = 0) {
        return this.dynamicParameters.get(name) ?? defaultValue;
    }
}

class DynamicVacuumTerm extends PhysicsTerm {
    constructor(amplitude = 1e-10, frequency = 1e-15) {
        super();
        this.amplitude = amplitude;
        this.frequency = frequency;
        this.metadata.set("type", "vacuum_energy");
        this.metadata.set("description", "Time-varying vacuum energy contribution");
    }

    compute(t, params) {
        const rho_vac = params.get("rho_vac_UA") ?? 7.09e-36;
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
        this.metadata.set("type", "quantum_coupling");
        this.metadata.set("description", "Non-local quantum coupling effects");
    }

    compute(t, params) {
        const hbar = params.get("hbar") ?? 1.0546e-34;
        const M = params.get("M") ?? 1.989e30;
        const r = params.get("r") ?? 1e4;
        return this.coupling_strength * (hbar * hbar) / (M * r * r) * Math.cos(t / 1e6);
    }

    getName() {
        return "QuantumCoupling";
    }

    getDescription() {
        return "Non-local quantum effects";
    }
}

class Source97UQFFModule {
    constructor() {
        // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;

        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();

        // Initialize framework defaults
        this.metadata.set("enhanced", "true");
        this.metadata.set("version", "2.0-Enhanced");
        this.metadata.set("module", "FeedbackFactorModule");
        this.metadata.set("description", "Feedback factor (f_feedback) in UQFF for black hole mass scaling");

        // Universal constants
        this.variables.set("f_feedback", 0.1);                  // Unitless, for ΔM_BH=1 dex
        this.variables.set("delta_M_BH_dex", 1.0);              // 1 dex = factor 10
        this.variables.set("M_bh_initial", 8.15e36);            // kg (Sgr A*)
        this.variables.set("k_4", 1.0);                         // Coupling for Ug4
        this.variables.set("rho_vac_SCm", 7.09e-37);            // J/m³
        this.variables.set("d_g", 2.55e20);                     // m
        this.variables.set("alpha", 0.001 / 86400.0);           // day^-1 to s^-1
        this.variables.set("pi", Math.PI);
        this.variables.set("t", 0.0);                           // s
        this.variables.set("t_n", 0.0);                         // s

        // Compute derived values
        this.computeM_bh_final();
    }

    // ========== DYNAMIC VARIABLE OPERATIONS ==========
    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
            if (this.enableLogging) {
                console.log(`Source97UQFFModule: Updated ${name} = ${value}`);
            }
        } else {
            console.warn(`Variable '${name}' not found. Adding with value ${value}`);
            this.variables.set(name, value);
        }

        // Recalculate dependent values
        if (name === "delta_M_BH_dex" || name === "M_bh_initial") {
            this.computeM_bh_final();
        }
    }

    addToVariable(name, delta) {
        const current = this.variables.get(name);
        if (current !== undefined) {
            this.variables.set(name, current + delta);
            if (this.enableLogging) {
                console.log(`Source97UQFFModule: Added ${delta} to ${name}, new value = ${current + delta}`);
            }
        } else {
            console.warn(`Variable '${name}' not found. Adding with delta ${delta}`);
            this.variables.set(name, delta);
        }

        // Recalculate dependent values
        if (name === "delta_M_BH_dex" || name === "M_bh_initial") {
            this.computeM_bh_final();
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    // ========== CORE COMPUTATIONS ==========
    computeF_feedback() {
        return this.variables.get("f_feedback");
    }

    computeDeltaM_BH() {
        return this.variables.get("delta_M_BH_dex");
    }

    computeM_bh_final() {
        const factor = Math.pow(10.0, this.computeDeltaM_BH());
        const initial = this.variables.get("M_bh_initial");
        const final_mass = initial * factor;
        this.variables.set("M_bh_final", final_mass);
        return final_mass;
    }

    computeU_g4(t, t_n) {
        const k_4 = this.variables.get("k_4");
        const rho_vac_SCm = this.variables.get("rho_vac_SCm");
        const M_bh = this.computeM_bh_final();  // Use final mass for feedback scenario
        const d_g = this.variables.get("d_g");
        const alpha = this.variables.get("alpha");
        const pi = this.variables.get("pi");
        const f_feedback = this.computeF_feedback();

        const exp_term = Math.exp(-alpha * t);
        const cos_term = Math.cos(pi * t_n);
        const feedback_factor = 1.0 + f_feedback;

        return k_4 * (rho_vac_SCm * M_bh / d_g) * exp_term * cos_term * feedback_factor;
    }

    computeU_g4_no_feedback(t, t_n) {
        const orig_f = this.variables.get("f_feedback");
        this.variables.set("f_feedback", 0.0);
        const result = this.computeU_g4(t, t_n);
        this.variables.set("f_feedback", orig_f);  // Restore
        return result;
    }

    // ========== SELF-EXPANDING FRAMEWORK METHODS ==========
    registerDynamicTerm(term) {
        this.dynamicTerms.push(term);
        if (this.enableLogging) {
            console.log(`Source97UQFFModule: Registered dynamic term ${term.getName()}`);
        }
    }

    setDynamicParameter(name, value) {
        this.dynamicParameters.set(name, value);
        if (this.enableLogging) {
            console.log(`Source97UQFFModule: Set dynamic parameter ${name} = ${value}`);
        }
    }

    getDynamicParameter(name, defaultValue = 0) {
        return this.dynamicParameters.get(name) ?? defaultValue;
    }

    computeDynamicTerms(t, params) {
        let total = 0;
        for (const term of this.dynamicTerms) {
            if (this.enableDynamicTerms) {
                total += term.compute(t, params);
            }
        }
        return total;
    }

    // ========== OUTPUT METHODS ==========
    getEquationText() {
        return "U_g4 = k_4 * (ρ_vac,[SCm] M_bh / d_g) * e^{-α t} * cos(π t_n) * (1 + f_feedback)\n" +
               "Where f_feedback = 0.1 (unitless, for ΔM_BH = 1 dex = 10x mass increase);\n" +
               "ΔM_BH =1 dex → M_bh_final = 10 * M_bh_initial (8.15e36 kg → 8.15e37 kg).\n" +
               "Example t=0, t_n=0: U_g4 ≈2.75e-20 J/m³ (with); ≈2.50e-20 J/m³ (without; +10%).\n" +
               "Role: Scales feedback in star-BH interactions; regulates AGN, mergers, star formation.\n" +
               "UQFF: Enhances energy density for 10x M_BH; resolves final parsec, quasar jets.";
    }

    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential(2)}`);
        }
    }

    printU_g4_comparison(t = 0.0, t_n = 0.0) {
        const u_with = this.computeU_g4(t, t_n);
        const u_without = this.computeU_g4_no_feedback(t, t_n);
        const delta_percent = ((u_with - u_without) / u_without) * 100.0;

        console.log(`U_g4 Comparison at t=${t} s, t_n=${t_n} s:`);
        console.log(`With feedback: ${u_with.toExponential(2)} J/m³`);
        console.log(`Without feedback: ${u_without.toExponential(2)} J/m³`);
        console.log(`Difference: +${delta_percent.toFixed(1)}%`);
    }

    printFeedbackAnalysis() {
        console.log("=== FEEDBACK FACTOR ANALYSIS ===");
        console.log(`f_feedback = ${this.computeF_feedback()} (unitless)`);
        console.log(`ΔM_BH = ${this.computeDeltaM_BH()} dex`);
        console.log(`M_bh_initial = ${this.variables.get("M_bh_initial").toExponential(2)} kg`);
        console.log(`M_bh_final = ${this.computeM_bh_final().toExponential(2)} kg`);
        console.log(`Mass increase factor = ${Math.pow(10, this.computeDeltaM_BH()).toFixed(0)}x`);
        console.log("");
        this.printU_g4_comparison(0.0, 0.0);
        console.log("");
        console.log("Physical Role:");
        console.log("- Scales feedback in star-black hole interactions");
        console.log("- Regulates AGN activity, galaxy mergers, star formation");
        console.log("- UQFF enhancement for 10x black hole mass scenarios");
        console.log("- Resolves final parsec problem and quasar jet dynamics");
    }

    // ========== STATE MANAGEMENT ==========
    exportState() {
        const state = {
            variables: Object.fromEntries(this.variables),
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            metadata: Object.fromEntries(this.metadata),
            enableDynamicTerms: this.enableDynamicTerms,
            enableLogging: this.enableLogging,
            learningRate: this.learningRate,
            dynamicTermsCount: this.dynamicTerms.length,
            timestamp: new Date().toISOString()
        };
        return JSON.stringify(state, null, 2);
    }

    importState(jsonString) {
        try {
            const state = JSON.parse(jsonString);
            this.variables = new Map(Object.entries(state.variables || {}));
            this.dynamicParameters = new Map(Object.entries(state.dynamicParameters || {}));
            this.metadata = new Map(Object.entries(state.metadata || {}));
            this.enableDynamicTerms = state.enableDynamicTerms ?? true;
            this.enableLogging = state.enableLogging ?? false;
            this.learningRate = state.learningRate ?? 0.001;

            if (this.enableLogging) {
                console.log("Source97UQFFModule: State imported successfully");
            }
        } catch (error) {
            console.error("Source97UQFFModule: Failed to import state:", error.message);
        }
    }

    // ========== UTILITY METHODS ==========
    setEnableLogging(enable) {
        this.enableLogging = enable;
    }

    setLearningRate(rate) {
        this.learningRate = rate;
    }

    getMetadata(key) {
        return this.metadata.get(key);
    }

    setMetadata(key, value) {
        this.metadata.set(key, value);
    }
}

// Example usage:
// const mod = new Source97UQFFModule();
// console.log("M_bh_final =", mod.computeM_bh_final().toExponential(2), "kg");
// mod.printU_g4_comparison(0.0, 0.0);
// console.log(mod.getEquationText());
// mod.updateVariable("delta_M_BH_dex", 2.0);  // 100x mass
// mod.printVariables();

// Export for use in other modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = { Source97UQFFModule };
}

module.exports = { Source97UQFFModule };