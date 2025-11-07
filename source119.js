// source119.js - StepFunctionModule
// JavaScript implementation of the Step Function (S) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes S(r - R_b) = 1 for r > R_b, 0 otherwise; activates U_g2 outside outer field bubble.
// Converted from Source119.cpp
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

class StepFunctionModule {
    constructor() {
        // Self-expanding framework members
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;

        // Initialize metadata
        this.metadata.set("enhanced", "true");
        this.metadata.set("version", "2.0-Enhanced");

        // Core variables
        this.variables = new Map();
        
        // Universal constants
        this.variables.set("R_b", 1.496e13);                    // m (100 AU)
        this.variables.set("k_2", 1.2);                         // Coupling
        this.variables.set("rho_vac_UA", 7.09e-36);             // J/m^3
        this.variables.set("rho_vac_SCm", 7.09e-37);            // J/m^3
        this.variables.set("M_s", 1.989e30);                    // kg
        this.variables.set("r", 1.496e13);                      // m (default = R_b)
        this.variables.set("delta_sw", 0.01);                   // Unitless
        this.variables.set("v_sw", 5e5);                        // m/s
        this.variables.set("H_SCm", 1.0);                       // Unitless
        this.variables.set("E_react", 1e46);                    // J

        // Derived
        this.variables.set("rho_sum", this.variables.get("rho_vac_UA") + this.variables.get("rho_vac_SCm"));
        this.variables.set("swirl_factor", 1.0 + this.variables.get("delta_sw") * this.variables.get("v_sw"));
    }

    // Update variable
    updateVariable(name, value) {
        this.variables.set(name, value);
        if (name === "rho_vac_UA" || name === "rho_vac_SCm") {
            this.variables.set("rho_sum", this.variables.get("rho_vac_UA") + this.variables.get("rho_vac_SCm"));
        } else if (name === "delta_sw" || name === "v_sw") {
            this.variables.set("swirl_factor", 1.0 + this.variables.get("delta_sw") * this.variables.get("v_sw"));
        }
    }

    // Add delta to variable
    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
            if (name === "rho_vac_UA" || name === "rho_vac_SCm") {
                this.variables.set("rho_sum", this.variables.get("rho_vac_UA") + this.variables.get("rho_vac_SCm"));
            } else if (name === "delta_sw" || name === "v_sw") {
                this.variables.set("swirl_factor", 1.0 + this.variables.get("delta_sw") * this.variables.get("v_sw"));
            }
        } else {
            this.variables.set(name, delta);
        }
    }

    // Subtract delta from variable
    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    // Compute S(r - R_b): 1 if r > R_b, 0 otherwise (treat = as 1)
    computeS_r_Rb(r) {
        return (r >= this.variables.get("R_b")) ? 1.0 : 0.0;
    }

    // Compute U_g2 with S(r - R_b)
    computeU_g2(r) {
        this.variables.set("r", r);
        const k_2 = this.variables.get("k_2");
        const rho_sum = this.variables.get("rho_sum");
        const M_s = this.variables.get("M_s");
        const s_step = this.computeS_r_Rb(r);
        const swirl_factor = this.variables.get("swirl_factor");
        const h_scm = this.variables.get("H_SCm");
        const e_react = this.variables.get("E_react");
        
        return k_2 * (rho_sum * M_s / (r * r)) * s_step * swirl_factor * h_scm * e_react;
    }

    // Get equation text
    getEquationText() {
        return "U_g2 = k_2 * [(ρ_vac,[UA] + ρ_vac,[SCm]) M_s / r^2] * S(r - R_b) * (1 + δ_sw v_sw) * H_SCm * E_react\n" +
               "Where S(r - R_b) = 1 (r > R_b), 0 otherwise (Heaviside step; =1 at boundary).\n" +
               "Defines outer bubble activation beyond R_b=1.496e13 m (100 AU).\n" +
               "Example r=1.496e13 m: S=1, U_g2 ≈1.18e53 J/m³;\n" +
               "r=1e11 m: S=0, U_g2=0; r=1e14 m: S=1, U_g2≈1.18e51 J/m³.\n" +
               "Role: Sharp transition internal/external gravity; heliopause-like boundary.\n" +
               "UQFF: Separates regimes for heliodynamics/nebular formation.";
    }

    // Print all current variables
    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }

    // Self-expanding framework methods
    registerDynamicTerm(term) {
        if (this.enableDynamicTerms) {
            this.dynamicTerms.push(term);
            if (this.enableLogging) {
                console.log(`Registered dynamic term: ${term.name || 'unnamed'}`);
            }
        }
    }

    setDynamicParameter(name, value) {
        this.dynamicParameters.set(name, value);
        if (this.enableLogging) {
            console.log(`Set dynamic parameter: ${name} = ${value}`);
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
        return {
            variables: Object.fromEntries(this.variables),
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            metadata: Object.fromEntries(this.metadata),
            enableDynamicTerms: this.enableDynamicTerms,
            learningRate: this.learningRate
        };
    }
}

// Export module
if (typeof module !== 'undefined' && module.exports) {
    module.exports = { StepFunctionModule };
}
