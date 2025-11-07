// source123.js - TimeReversalZoneModule
// JavaScript implementation of Time-Reversal Zone Factor (f_TRZ) in UQFF
// Computes f_TRZ=0.1 (unitless); scales (1 + f_TRZ) in Universal Inertia U_i
// Converted from source123.cpp
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

class TimeReversalZoneModule {
    constructor() {
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;
        
        this.metadata.set("enhanced", "true");
        this.metadata.set("version", "2.0-Enhanced");
        
        this.variables = new Map();
        // Time-reversal zone parameters
        this.variables.set("f_TRZ", 0.1);                   // Unitless TRZ factor
        this.variables.set("lambda_i", 1.0);                // Coupling
        this.variables.set("rho_vac_SCm", 7.09e-37);        // J/m³
        this.variables.set("rho_vac_UA", 7.09e-36);         // J/m³
        this.variables.set("omega_s", 2.5e-6);              // rad/s
        this.variables.set("pi", Math.PI);
        this.variables.set("t", 0.0);                       // s
        this.variables.set("t_n", 0.0);                     // s
        
        // Derived
        this.variables.set("rho_product", 
            this.variables.get("rho_vac_SCm") * this.variables.get("rho_vac_UA"));
        this.variables.set("trz_factor", this.computeTRZFactor());
    }

    updateVariable(name, value) {
        this.variables.set(name, value);
        if (name === "f_TRZ") {
            this.variables.set("trz_factor", this.computeTRZFactor());
        } else if (name === "rho_vac_SCm" || name === "rho_vac_UA") {
            this.variables.set("rho_product", 
                this.variables.get("rho_vac_SCm") * this.variables.get("rho_vac_UA"));
        }
    }

    addToVariable(name, delta) {
        const current = this.variables.get(name) || 0;
        this.updateVariable(name, current + delta);
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    computeF_TRZ() {
        return this.variables.get("f_TRZ");
    }

    computeTRZFactor() {
        return 1.0 + this.computeF_TRZ();
    }

    computeU_i_base(t, t_n) {
        const lambda_i = this.variables.get("lambda_i");
        const rho_product = this.variables.get("rho_product");
        const omega_s_t = this.variables.get("omega_s");
        const pi = this.variables.get("pi");
        const cos_pi_tn = Math.cos(pi * t_n);
        return lambda_i * rho_product * omega_s_t * cos_pi_tn;
    }

    computeU_i(t, t_n) {
        this.variables.set("t", t);
        const base = this.computeU_i_base(t, t_n);
        const trz_f = this.computeTRZFactor();
        return base * trz_f;
    }

    computeU_i_no_TRZ(t, t_n) {
        const orig_f = this.variables.get("f_TRZ");
        this.variables.set("f_TRZ", 0.0);
        const result = this.computeU_i(t, t_n);
        this.variables.set("f_TRZ", orig_f);
        return result;
    }

    getEquationText() {
        return "U_i = λ_i * ρ_vac,[SCm] * ρ_vac,[UA] * ω_s(t) * cos(π t_n) * (1 + f_TRZ)\n" +
               "Where f_TRZ = 0.1 (unitless time-reversal zone factor; +10% negentropic enhancement)\n" +
               "TRZ: Regions for time-reversal/negentropy (COP>1, vacuum extraction).\n" +
               "Example Sun t=0, t_n=0: U_i ≈1.38e-47 J/m³ (with); ≈1.25e-47 J/m³ (without; -9.1%).\n" +
               "In F_U: -∇ ρ_i U_i E_react (resistive, TRZ-boosted).\n" +
               "Role: Stabilizes via negentropy; TRZ in nebulae/formation/mergers/biology.\n" +
               "UQFF: Integrates pondermotive force/time asymmetry; Aether superfluid effects.";
    }

    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables.entries()) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }

    printUiComparison(t = 0.0, t_n = 0.0) {
        const u_i_with = this.computeU_i(t, t_n);
        const u_i_without = this.computeU_i_no_TRZ(t, t_n);
        const percent_increase = ((u_i_with - u_i_without) / u_i_without) * 100.0;
        console.log(`U_i Comparison at t=${t} s, t_n=${t_n}:`);
        console.log(`With TRZ: ${u_i_with.toExponential()} J/m³`);
        console.log(`Without TRZ: ${u_i_without.toExponential()} J/m³`);
        console.log(`Increase: +${percent_increase.toFixed(1)}%`);
    }

    exportState() {
        return {
            variables: Object.fromEntries(this.variables),
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            metadata: Object.fromEntries(this.metadata)
        };
    }
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { TimeReversalZoneModule };
}
