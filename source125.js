// source125.js - Ug3DiskVectorModule
// JavaScript implementation of Unit Vector in Ug3 Disk Plane (φ̂_j) in UQFF
// Computes φ̂_j (unit vector, magnitude=1; [cos θ_j, sin θ_j, 0]); scales in U_m
// Converted from source125.cpp
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

class Ug3DiskVectorModule {
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
        // Ug3 disk vector parameters
        this.variables.set("theta_j", 0.0);                 // rad (azimuthal angle)
        this.variables.set("mu_j", 3.38e23);                // T²m³ (j=1)
        this.variables.set("r_j", 1.496e13);                // m
        this.variables.set("gamma", 5e-5 / 86400.0);        // s⁻¹
        this.variables.set("t_n", 0.0);                     // s
        this.variables.set("P_SCm", 1.0);                   // Pressure
        this.variables.set("E_react", 1e46);                // J
        this.variables.set("f_Heaviside", 0.01);            // Unitless
        this.variables.set("f_quasi", 0.01);                // Unitless
        this.variables.set("pi", Math.PI);
        
        // Derived
        this.variables.set("scale_Heaviside", 1e13);
        this.variables.set("heaviside_factor", 
            1.0 + this.variables.get("scale_Heaviside") * this.variables.get("f_Heaviside"));
    }

    updateVariable(name, value) {
        this.variables.set(name, value);
        if (name === "f_Heaviside" || name === "scale_Heaviside") {
            this.variables.set("heaviside_factor",
                1.0 + this.variables.get("scale_Heaviside") * this.variables.get("f_Heaviside"));
        }
    }

    addToVariable(name, delta) {
        const current = this.variables.get(name) || 0;
        this.updateVariable(name, current + delta);
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    computePhiHat_j(j) {
        const theta = this.variables.get("theta_j");
        return [Math.cos(theta), Math.sin(theta), 0.0];
    }

    computePhiHatMagnitude(j) {
        const phi = this.computePhiHat_j(j);
        return Math.sqrt(phi[0]*phi[0] + phi[1]*phi[1] + phi[2]*phi[2]);  // =1
    }

    computeUmBase(t) {
        const mu_over_rj = this.variables.get("mu_j") / this.variables.get("r_j");
        const gamma = this.variables.get("gamma");
        const pi = this.variables.get("pi");
        const t_n = this.variables.get("t_n");
        const exp_arg = -gamma * t * Math.cos(pi * t_n);
        const one_minus_exp = 1.0 - Math.exp(exp_arg);
        const phi_mag = this.computePhiHatMagnitude(1);  // =1
        const p_scm = this.variables.get("P_SCm");
        const e_react = this.variables.get("E_react");
        return mu_over_rj * one_minus_exp * phi_mag * p_scm * e_react;
    }

    computeUmContribution(t, j) {
        const base = this.computeUmBase(t);
        const heaviside_f = this.variables.get("heaviside_factor");
        const quasi_f = 1.0 + this.variables.get("f_quasi");
        return base * heaviside_f * quasi_f;
    }

    getEquationText() {
        return "U_m = Σ_j [ (μ_j / r_j) (1 - e^{-γ t cos(π t_n)}) φ̂_j ] P_SCm E_react (1 + 10¹³ f_Heaviside) (1 + f_quasi)\n" +
               "Where φ̂_j = [cos θ_j, sin θ_j, 0] (unit vector in Ug3 disk plane, |φ̂_j|=1)\n" +
               "Specifies azimuthal direction for j-th string in disk (e.g., galactic plane).\n" +
               "Example j=1, θ_j=0, t=0: φ̂_j=[1,0,0], U_m ≈2.28e65 J/m³ (mag=1).\n" +
               "Role: Directional geometry for magnetic contributions in disks/nebulae.\n" +
               "UQFF: Vector orientation in U_m/U_g3; collimation in jets/disks/formation.";
    }

    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables.entries()) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }

    printVectorAndUm(j = 1, t = 0.0) {
        const phi = this.computePhiHat_j(j);
        const mag = this.computePhiHatMagnitude(j);
        const um = this.computeUmContribution(t, j);
        const theta_j = this.variables.get("theta_j");
        console.log(`φ̂_${j} at θ_j=${theta_j} rad, t=${t} s:`);
        console.log(`φ̂_j = [${phi[0].toExponential()}, ${phi[1].toExponential()}, ${phi[2].toExponential()}] (mag=${mag})`);
        console.log(`U_m contrib = ${um.toExponential()} J/m³`);
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
    module.exports = { Ug3DiskVectorModule };
}
