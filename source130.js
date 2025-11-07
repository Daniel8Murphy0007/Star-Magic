// source130.js - UniversalInertiaVacuumModule
// JavaScript implementation of Vacuum Energy Density of Universal Inertia (ρ_vac,Ui) in UQFF
// Computes ρ_vac,Ui = 2.84e-36 J/m³ (Sun, level 13); reference scale for U_i
// Converted from source130.cpp
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

class UniversalInertiaVacuumModule130 {
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
        this.variables.set("rho_vac_Ui", 2.84e-36);         // J/m³
        this.variables.set("lambda_i", 1.0);                // Unitless
        this.variables.set("rho_vac_SCm", 7.09e-37);        // J/m³
        this.variables.set("rho_vac_UA", 7.09e-36);         // J/m³
        this.variables.set("omega_s", 2.5e-6);              // rad/s
        this.variables.set("f_TRZ", 0.1);                   // Unitless
        this.variables.set("pi", Math.PI);
        this.variables.set("t", 0.0);                       // s
        this.variables.set("t_n", 0.0);                     // s
        
        this.variables.set("rho_product", this.variables.get("rho_vac_SCm") * this.variables.get("rho_vac_UA"));
    }

    updateVariable(name, value) {
        this.variables.set(name, value);
        if (name === "rho_vac_SCm" || name === "rho_vac_UA") {
            this.variables.set("rho_product", this.variables.get("rho_vac_SCm") * this.variables.get("rho_vac_UA"));
        }
    }

    addToVariable(name, delta) {
        const current = this.variables.get(name) || 0;
        this.updateVariable(name, current + delta);
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    computeRho_vac_Ui() {
        return this.variables.get("rho_vac_Ui");
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
        const trz_factor = 1.0 + this.variables.get("f_TRZ");
        return base * trz_factor;
    }

    getEquationText() {
        return "U_i = λ_i * ρ_vac,[SCm] * ρ_vac,[UA] * ω_s(t) * cos(π t_n) * (1 + f_TRZ)\n" +
               "ρ_vac,Ui = 2.84e-36 J/m³ (Sun level 13, inertia vacuum scale; not direct in eq.).\n" +
               "Provides reference for U_i magnitude; inertial resistance from [SCm]/[UA].\n" +
               "Example t=0, t_n=0: U_i ≈1.38e-47 J/m³ (consistent scale with ρ_vac,Ui).\n" +
               "In F_U: -∑ λ_i U_i E_react (resistive inertia).\n" +
               "Role: Quantifies vacuum inertia energy; opposes dynamics in nebulae/formation.\n" +
               "UQFF: Small-scale reference for cosmic inertia; [SCm]-[UA] resistance.";
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
    module.exports = { UniversalInertiaVacuumModule130 };
}
