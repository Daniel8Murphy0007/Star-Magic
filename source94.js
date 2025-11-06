class Source94UQFFModule {
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

        // Coupling constants (unitless)
        this.k_values = [1.5, 1.2, 1.8, 1.0];       // k1=1.5, k2=1.2, k3=1.8, k4=1.0

        // U_gi defaults (example from Sun at t=0, J/m^3)
        this.variables.set("U_g1", 1.39e26);         // Internal Dipole
        this.variables.set("U_g2", 1.18e53);         // Outer Field Bubble
        this.variables.set("U_g3", 1.8e49);          // Magnetic Strings Disk
        this.variables.set("U_g4", 2.50e-20);        // Star-Black Hole Interactions

        // Shared params (placeholders)
        this.variables.set("mu_s", 1.0);             // Magnetic moment
        this.variables.set("M_s", 1.989e30);         // Stellar mass kg
        this.variables.set("r", 1e11);               // m
        this.variables.set("alpha", 1e-10);          // Decay rate s^-1
        this.variables.set("t_n", 0.0);              // s
        this.variables.set("pi", Math.PI);
        this.variables.set("delta_def", 0.0);        // Deformation
        this.variables.set("rho_vac_UA", 7.09e-36);  // J/m^3
        this.variables.set("rho_vac_SCm", 7.09e-37); // J/m^3
        this.variables.set("S_r_Rb", 1.0);           // Step function
        this.variables.set("delta_sw", 0.0);         // Swirl deformation
        this.variables.set("v_sw", 0.0);             // Solar wind velocity
        this.variables.set("H_SCm", 1.0);            // Heaviside SCm
        this.variables.set("E_react", 1.0);          // Reactive energy
        this.variables.set("M_bh", 8.15e36);         // kg
        this.variables.set("d_g", 2.55e20);          // m
        this.variables.set("f_feedback", 0.0);       // Feedback factor
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
            if (this.enableLogging) {
                console.log(`Updated variable ${name} = ${value}`);
            }
        } else {
            if (this.enableLogging) {
                console.log(`Variable '${name}' not found. Adding with value ${value}`);
            }
            this.variables.set(name, value);
        }
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            const currentValue = this.variables.get(name);
            this.variables.set(name, currentValue + delta);
            if (this.enableLogging) {
                console.log(`Added ${delta} to ${name}, new value = ${currentValue + delta}`);
            }
        } else {
            if (this.enableLogging) {
                console.log(`Variable '${name}' not found. Adding with delta ${delta}`);
            }
            this.variables.set(name, delta);
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    // ========== CORE COMPUTATIONS ==========
    computeK_i(i) {
        if (i < 1 || i > 4) {
            console.warn(`Invalid i: ${i}. Using k1.`);
            return this.k_values[0];
        }
        return this.k_values[i-1];
    }

    computeU_gi(i) {
        const key = `U_g${i}`;
        if (this.variables.has(key)) {
            return this.variables.get(key);
        }
        console.warn(`U_g${i} not defined. Returning 0.`);
        return 0.0;
    }

    computeK_Ugi(i) {
        return this.computeK_i(i) * this.computeU_gi(i);
    }

    computeAllK_Ugi() {
        const k_ugi = [];
        for (let i = 1; i <= 4; i++) {
            k_ugi.push(this.computeK_Ugi(i));
        }
        return k_ugi;
    }

    computeSumK_Ugi() {
        const all = this.computeAllK_Ugi();
        return all.reduce((sum, val) => sum + val, 0);
    }

    // ========== OUTPUT METHODS ==========
    getEquationText() {
        return "F_U = ∑ [k_i * U_gi(r,t,M_s,μ_s,T_s,B_s,ρ_vac,[SCm],ρ_vac,[UA],t_n) - β_i * ... ] + other terms\n" +
               "k_i (unitless): k1=1.5 (Ug1 Internal Dipole), k2=1.2 (Ug2 Outer Bubble),\n" +
               "k3=1.8 (Ug3 Magnetic Disk), k4=1.0 (Ug4 Star-BH).\n" +
               "U_g1 = k1 * μ_s √(M_s/r) e^{-α t} cos(π t_n) (1+δ_def);\n" +
               "U_g2 = k2 * (ρ_UA + ρ_SCm) M_s / r^2 * S(r-R_b) (1+δ_sw v_sw) H_SCm E_react;\n" +
               "U_g3 = k3 * (ρ_SCm + ρ_UA) Ω_g M_s / d_g * cos(π t_n) (1 + f_mag_str);\n" +
               "U_g4 = k4 * ρ_SCm M_bh / d_g * e^{-α t} cos(π t_n) (1 + f_feedback).\n" +
               "Example Sun t=0: ∑ k_i U_gi ≈1.42e53 J/m³ (Ug2 dominant).\n" +
               "Role: Scales Ug strengths; normalizes contributions in F_U.";
    }

    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential(6)}`);
        }
        console.log(`k_i values: k1=${this.k_values[0]}, k2=${this.k_values[1]}, k3=${this.k_values[2]}, k4=${this.k_values[3]}`);
    }

    printK_Ugi() {
        const all = this.computeAllK_Ugi();
        console.log("Scaled Ug Terms k_i * U_gi (J/m³):");
        for (let i = 1; i <= 4; i++) {
            console.log(`k${i} * U_g${i} = ${all[i-1].toExponential(6)}`);
        }
        console.log(`Sum ∑ k_i U_gi = ${this.computeSumK_Ugi().toExponential(6)}`);
    }

    // ========== COMPUTATION METHODS ==========
    computeUgCoupling() {
        return {
            k_values: this.k_values.slice(),
            allK_Ugi: this.computeAllK_Ugi(),
            sumK_Ugi: this.computeSumK_Ugi(),
            equation: this.getEquationText()
        };
    }
}

module.exports = Source94UQFFModule;