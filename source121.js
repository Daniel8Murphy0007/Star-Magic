// source121.js - SurfaceMagneticFieldModule
// JavaScript implementation of Surface Magnetic Field (B_s) in UQFF
// Computes B_s range [1e-4, 0.4] T for Sun; influences U_g3 magnetic strings
// Converted from source121.cpp
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

class SurfaceMagneticFieldModule {
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
        // Surface magnetic field parameters
        this.variables.set("B_s_min", 1e-4);                    // T (quiet Sun)
        this.variables.set("B_s_max", 0.4);                     // T (sunspot max)
        this.variables.set("B_ref", 0.4);                       // T (reference field)
        this.variables.set("B_s", 1e-4);                        // T (current field)
        this.variables.set("k_3", 1.0);                         // U_g3 coupling
        this.variables.set("omega_s", 2 * Math.PI / (11 * 365.25 * 24 * 3600)); // rad/s (solar cycle)
        this.variables.set("phi", 0.0);                         // Phase
        this.variables.set("P_core", 1.0);                      // Core penetration
        this.variables.set("E_react", 1e46);                    // J (reactive energy)
        this.variables.set("M_s", 1.989e30);                    // kg
        this.variables.set("r", 6.96e8);                        // m
    }

    updateVariable(name, value) {
        this.variables.set(name, value);
    }

    computeB_j(t) {
        const B_s = this.variables.get("B_s");
        const B_ref = this.variables.get("B_ref");
        const omega_s = this.variables.get("omega_s");
        const phi = this.variables.get("phi");
        return (B_s / B_ref) * Math.cos(omega_s * t + phi);
    }

    computeU_g3_example(t) {
        const k_3 = this.variables.get("k_3");
        const M_s = this.variables.get("M_s");
        const r = this.variables.get("r");
        const B_j = this.computeB_j(t);
        const P_core = this.variables.get("P_core");
        const E_react = this.variables.get("E_react");
        return k_3 * (M_s / (r * r)) * B_j * P_core * E_react;
    }

    getEquationText() {
        return "B_j = (B_s / B_ref) * cos(ω_s t + φ)\n" +
               "U_g3 = k_3 * (M_s / r²) * B_j * P_core * E_react\n" +
               "Where B_s ∈ [1e-4, 0.4] T (quiet to active Sun), B_ref=0.4 T.\n" +
               "B_s=1e-4 T → U_g3≈4.5e45 J/m³; B_s=0.4 T → U_g3≈1.8e49 J/m³.\n" +
               "Role: Surface magnetic field modulates magnetic string energy density.\n" +
               "UQFF: Solar cycle variation in U_g3 via B_s oscillation.";
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
    module.exports = { SurfaceMagneticFieldModule };
}
