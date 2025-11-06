class Source93UQFFModule {
    constructor() {
        // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;

        // Initialize metadata
        this.metadata.set("enhanced", "true");
        this.metadata.set("version", "2.0-Enhanced");

        // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
        this.variables = new Map();

        // Universal constants
        this.variables.set("epsilon_sw", 0.001);                // Buoyancy modulation (unitless)
        this.variables.set("rho_vac_sw", 8e-21);                // J/m³ (solar wind energy density)
        this.variables.set("beta_1", 0.6);                      // From buoyancy coupling
        this.variables.set("U_g1", 1.39e26);                    // J/m³ (Ug1 example)
        this.variables.set("Omega_g", 7.3e-16);                 // rad/s
        this.variables.set("M_bh", 8.15e36);                    // kg
        this.variables.set("d_g", 2.55e20);                     // m
        this.variables.set("E_react", 1.0);                     // Normalized
        this.variables.set("U_UA", 1.0);                        // Universal Aether factor
        this.variables.set("t_n", 0.0);                         // s
        this.variables.set("pi", Math.PI);

        // Derived defaults
        this.variables.set("modulation_factor", this.computeModulationFactor());
    }

    // ========== SELF-EXPANDING FRAMEWORK METHODS ==========
    registerDynamicTerm(term) {
        if (this.enableDynamicTerms) {
            this.dynamicTerms.push(term);
            if (this.enableLogging) {
                console.log(`Registered dynamic term: ${term.getName()}`);
            }
        }
    }

    setDynamicParameter(key, value) {
        this.dynamicParameters.set(key, value);
        if (this.enableLogging) {
            console.log(`Set dynamic parameter ${key} = ${value}`);
        }
    }

    getDynamicParameter(key) {
        return this.dynamicParameters.get(key);
    }

    setEnableLogging(enable) {
        this.enableLogging = enable;
    }

    setLearningRate(rate) {
        this.learningRate = rate;
    }

    exportState(filename) {
        const state = {
            variables: Object.fromEntries(this.variables),
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            metadata: Object.fromEntries(this.metadata),
            enableDynamicTerms: this.enableDynamicTerms,
            enableLogging: this.enableLogging,
            learningRate: this.learningRate
        };
        // In a real implementation, this would write to a file
        if (this.enableLogging) {
            console.log(`Exporting state to ${filename}`);
        }
        return state;
    }

    // ========== DYNAMIC VARIABLE OPERATIONS ==========
    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
            if (name === "epsilon_sw" || name === "rho_vac_sw") {
                this.variables.set("modulation_factor", this.computeModulationFactor());
            }
        } else {
            console.warn(`Variable '${name}' not found. Adding with value ${value}`);
            this.variables.set(name, value);
        }
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            const currentValue = this.variables.get(name);
            this.variables.set(name, currentValue + delta);
            if (name === "epsilon_sw" || name === "rho_vac_sw") {
                this.variables.set("modulation_factor", this.computeModulationFactor());
            }
        } else {
            console.warn(`Variable '${name}' not found. Adding with delta ${delta}`);
            this.variables.set(name, delta);
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    // ========== CORE COMPUTATIONS ==========
    computeEpsilon_sw() {
        return this.variables.get("epsilon_sw");
    }

    computeModulationFactor() {
        const epsilon_sw = this.variables.get("epsilon_sw");
        const rho_vac_sw = this.variables.get("rho_vac_sw");
        return 1.0 + epsilon_sw * rho_vac_sw;
    }

    computeU_b1() {
        const beta_1 = this.variables.get("beta_1");
        const U_g1 = this.variables.get("U_g1");
        const Omega_g = this.variables.get("Omega_g");
        const M_bh_over_d_g = this.variables.get("M_bh") / this.variables.get("d_g");
        const E_react = this.variables.get("E_react");
        const mod_factor = this.computeModulationFactor();
        const U_UA = this.variables.get("U_UA");
        const pi = this.variables.get("pi");
        const t_n = this.variables.get("t_n");
        const cos_term = Math.cos(pi * t_n);

        return -beta_1 * U_g1 * Omega_g * M_bh_over_d_g * E_react * mod_factor * U_UA * cos_term;
    }

    // ========== OUTPUT METHODS ==========
    getEquationText() {
        return "Modulation Factor = 1 + ε_sw * ρ_vac,sw\n" +
               "Where ε_sw = 0.001 (unitless); ρ_vac,sw = 8e-21 J/m³.\n" +
               "In U_bi: ... * (1 + ε_sw * ρ_vac,sw) * ... ≈1 (negligible correction ~8e-24).\n" +
               "U_b1 = -β₁ U_g1 Ω_g (M_bh / d_g) * modulation * U_UA * cos(π t_n)\n" +
               "≈ -1.94e27 J/m³ (at t_n=0, Sun params; modulation ≈1).\n" +
               "Role: Minor solar wind density effect on buoyancy; stabilizes heliosphere/nebulae.\n" +
               "UQFF: Scales counterforce to Ug; negligible but flexible for variations.";
    }

    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential(6)}`);
        }
    }

    // ========== COMPUTATION METHODS ==========
    computeSolarWindModulation() {
        return {
            epsilon_sw: this.computeEpsilon_sw(),
            modulation_factor: this.computeModulationFactor(),
            u_b1: this.computeU_b1(),
            equation: this.getEquationText()
        };
    }
}

module.exports = Source93UQFFModule;