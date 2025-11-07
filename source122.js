// source122.js - SurfaceTemperatureModule
// JavaScript implementation of Surface Temperature (T_s) in UQFF
// Computes T_s=5778 K (Sun); scales magnetic field B_j by T_s/T_s_ref
// Converted from source122.cpp
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

class SurfaceTemperatureModule {
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
        // Surface temperature parameters
        this.variables.set("T_s", 5778.0);              // K (Sun effective)
        this.variables.set("T_s_ref", 5778.0);          // K (reference)
        this.variables.set("k_3", 1.8);                 // Coupling
        this.variables.set("B_ref", 1e3);               // T (base magnetic field)
        this.variables.set("omega_s", 2.5e-6);          // rad/s (solar cycle)
        this.variables.set("P_core", 1.0);              // Unitless
        this.variables.set("E_react", 1e46);            // J
        this.variables.set("pi", Math.PI);
        this.variables.set("t", 0.0);                   // s
    }

    updateVariable(name, value) {
        this.variables.set(name, value);
    }

    addToVariable(name, delta) {
        const current = this.variables.get(name) || 0;
        this.variables.set(name, current + delta);
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    computeT_s() {
        return this.variables.get("T_s");
    }

    computeB_j_hypothetical(t, T_s) {
        this.variables.set("t", t);
        const B_ref = this.variables.get("B_ref");
        const omega_s = this.variables.get("omega_s");
        const T_s_ref = this.variables.get("T_s_ref");
        const base_b = B_ref + 0.4 * Math.sin(omega_s * t);
        return base_b * (T_s / T_s_ref);
    }

    computeU_g3_example(t, T_s) {
        const k_3 = this.variables.get("k_3");
        const b_j = this.computeB_j_hypothetical(t, T_s);
        const omega_s = this.variables.get("omega_s");
        const pi = this.variables.get("pi");
        const P_core = this.variables.get("P_core");
        const E_react = this.variables.get("E_react");
        const cos_term = Math.cos(omega_s * t * pi);
        return k_3 * b_j * cos_term * P_core * E_react;
    }

    getEquationText() {
        return "B_j ≈ (10³ + 0.4 sin(ω_s t)) * (T_s / T_s_ref) T (hypothetical)\n" +
               "U_g3 = k_3 * B_j * cos(ω_s t π) * P_core * E_react\n" +
               "Where T_s = 5778 K (Sun effective photosphere; °C=5505).\n" +
               "T_s_ref=5778 K; scales string fields by temperature.\n" +
               "Example t=0, T_s=5778 K: B_j≈1e3 T, U_g3≈1.8e49 J/m³\n" +
               "T_s=10000 K: B_j≈1730 T, U_g3≈3.11e49 J/m³ (+73%).\n" +
               "Role: Thermal baseline for magnetic strength; variability in U_g3/disks.\n" +
               "UQFF: Temperature-dependent fields; extensible for radiation/formation.";
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
    module.exports = { SurfaceTemperatureModule };
}
