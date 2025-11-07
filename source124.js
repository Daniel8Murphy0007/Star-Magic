// source124.js - Ug1DefectModule
// JavaScript implementation of Ug1 Defect Factor (δ_def) in UQFF
// Computes δ_def = 0.01 * sin(0.001 t); scales (1 + δ_def) in Universal Gravity U_g1
// Converted from source124.cpp
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

class Ug1DefectModule {
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
        // U_g1 defect parameters
        this.variables.set("amplitude", 0.01);              // Unitless
        this.variables.set("freq", 0.001);                  // day⁻¹
        this.variables.set("k_1", 1.5);                     // Coupling
        this.variables.set("mu_s", 3.38e23);                // T²m³
        this.variables.set("M_s", 1.989e30);                // kg
        this.variables.set("alpha", 0.001);                 // day⁻¹
        this.variables.set("t_n", 0.0);                     // days
        this.variables.set("pi", Math.PI);
        this.variables.set("t_day", 0.0);                   // days
        this.variables.set("r", 1.496e11);                  // m (Earth-Sun example)
        
        // Derived
        this.variables.set("period_days", 2.0 * Math.PI / this.variables.get("freq"));
    }

    updateVariable(name, value) {
        this.variables.set(name, value);
        if (name === "freq") {
            this.variables.set("period_days", 2.0 * Math.PI / value);
        }
    }

    addToVariable(name, delta) {
        const current = this.variables.get(name) || 0;
        this.updateVariable(name, current + delta);
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    computeDelta_def(t_day) {
        this.variables.set("t_day", t_day);
        const amplitude = this.variables.get("amplitude");
        const freq = this.variables.get("freq");
        return amplitude * Math.sin(freq * t_day);
    }

    computeU_g1(t_day, r) {
        this.variables.set("r", r);
        const k_1 = this.variables.get("k_1");
        const mu_s = this.variables.get("mu_s");
        const M_s = this.variables.get("M_s");
        const grad_ms_r = M_s / (r * r);  // Approx ∇(M_s / r) = M_s / r²
        const alpha = this.variables.get("alpha");
        const exp_term = Math.exp(-alpha * t_day);
        const pi = this.variables.get("pi");
        const t_n = this.variables.get("t_n");
        const cos_tn = Math.cos(pi * t_n);
        const defect_factor = 1.0 + this.computeDelta_def(t_day);
        return k_1 * mu_s * grad_ms_r * exp_term * cos_tn * defect_factor;
    }

    computePeriod_years() {
        return this.variables.get("period_days") / 365.25;
    }

    getEquationText() {
        return "U_g1 = k_1 * μ_s * ∇(M_s / r) * e^{-α t} * cos(π t_n) * (1 + δ_def)\n" +
               "Where δ_def = 0.01 * sin(0.001 t) (unitless, t days; period ~17.22 yr).\n" +
               "Small oscillatory defect (~±1%) in internal dipole gravity.\n" +
               "Example t=0, r=1.496e11 m: δ_def=0, U_g1 ≈4.51e31 J/m³\n" +
               "t=1570.8 days: δ_def=0.01, U_g1 ≈4.56e31 J/m³ (+1.1%).\n" +
               "Role: Time-dependent perturbations; internal dynamics/[SCm] variations.\n" +
               "UQFF: Cyclic defects in stellar gravity; for formation/nebular stability.";
    }

    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables.entries()) {
            console.log(`${key} = ${value.toExponential()}`);
        }
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
    module.exports = { Ug1DefectModule };
}
