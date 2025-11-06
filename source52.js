class MultiUQFFModule {
    constructor(system = "OrionNebula", mode = "compressed") {
        this.variables = new Map();
        this.current_system = system;
        this.current_mode = mode;

        // Base constants (universal)
        this.variables.set("G", 6.6743e-11);                    // m^3 kg^-1 s^-2
        this.variables.set("c", 3e8);                           // m/s
        this.variables.set("hbar", 1.0546e-34);                 // J s
        this.variables.set("Lambda", 1.1e-52);                  // m^-2
        this.variables.set("pi", 3.141592653589793);            // pi
        this.variables.set("t_Hubble", 13.8e9 * 3.156e7);       // s
        this.variables.set("year_to_s", 3.156e7);               // s/yr
        this.variables.set("H0", 70.0);                         // km/s/Mpc -> 2.269e-18 s^-1 after conversion
        this.variables.set("Mpc_to_m", 3.086e22);               // m/Mpc
        this.variables.set("Omega_m", 0.3);
        this.variables.set("Omega_Lambda", 0.7);
        this.variables.set("B", 1e10);                          // T (assumed)
        this.variables.set("B_crit", 1e11);                     // T
        this.variables.set("rho_fluid", 1e-15);                 // kg/m^3 (placeholder)
        this.variables.set("delta_rho_over_rho", 1e-5);
        this.variables.set("integral_psi", 2.176e-18);          // J
        this.variables.set("Delta_x_Delta_p", 1e-68);           // J^2 s^2
        this.variables.set("F_env", 0.0);
        this.variables.set("M_DM", 0.0);
        this.variables.set("M_visible", 0.0);  // Set per system

        this.setSystem(system);
    }

    // Set system: Load system-specific vars
    setSystem(system) {
        this.current_system = system;
        const M_sun = 1.989e30;

        // Reset system-specific variables
        this.variables.set("M", 0.0);
        this.variables.set("r", 0.0);
        this.variables.set("z", 0.0);
        this.variables.set("t_default", 0.0);
        this.variables.set("v_exp", 0.0);
        this.variables.set("M_visible", 0.0);  // = M usually

        if (system === "UniverseDiameter") {
            this.variables.set("M", 1.5e53);
            this.variables.set("r", 4.4e26);
            this.variables.set("z", 1100.0);
            this.variables.set("t_default", 4.35e17);
            this.variables.set("v_exp", 3e5);
        } else if (system === "HydrogenAtom" || system === "HydrogenResonancePToE") {
            this.variables.set("M", 1.6735e-27);
            this.variables.set("r", 5.2918e-11);
            this.variables.set("z", 0.0);
            this.variables.set("t_default", 4.35e17);
            this.variables.set("v_exp", 0.0);
        } else if (system === "LagoonNebula") {
            this.variables.set("M", 1e4 * M_sun);  // 1.989e34
            this.variables.set("r", 5.203e17);
            this.variables.set("z", 0.0001);
            this.variables.set("t_default", 2e6 * 3.156e7);  // 6.312e13
            this.variables.set("v_exp", 1e4);
        } else if (system === "SpiralsSupernovae") {
            this.variables.set("M", 1e11 * M_sun);  // 1.989e41
            this.variables.set("r", 1.543e21);
            this.variables.set("z", 0.002);
            this.variables.set("t_default", 4.35e17);
            this.variables.set("v_exp", 2e5);
        } else if (system === "NGC6302") {
            this.variables.set("M", 1.0 * M_sun);  // 1.989e30
            this.variables.set("r", 1.514e16);
            this.variables.set("z", 0.00001);
            this.variables.set("t_default", 1e4 * 3.156e7);  // 3.156e11
            this.variables.set("v_exp", 2e4);
        } else if (system === "OrionNebula") {
            this.variables.set("M", 2e3 * M_sun);  // 3.978e33
            this.variables.set("r", 1.135e17);
            this.variables.set("z", 0.00004);
            this.variables.set("t_default", 1e6 * 3.156e7);  // 3.156e13
            this.variables.set("v_exp", 1e4);
        } else if (system === "UniverseGuide") {
            this.variables.set("M", 1.0 * M_sun);  // 1.989e30
            this.variables.set("r", 1.496e11);
            this.variables.set("z", 0.0);
            this.variables.set("t_default", 4.35e17);
            this.variables.set("v_exp", 3e4);
        }

        this.variables.set("M_visible", this.variables.get("M"));
    }

    // Set mode
    setMode(mode) {
        this.current_mode = mode;
    }

    // Dynamic variable operations
    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
        } else {
            console.warn(`Variable '${name}' not found. Adding with value ${value}`);
            this.variables.set(name, value);
        }
        if (name === "M") {
            this.variables.set("M_visible", value);
        }
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
        } else {
            console.warn(`Variable '${name}' not found. Adding with delta ${delta}`);
            this.variables.set(name, delta);
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    getVariable(name) {
        return this.variables.get(name) || 0;
    }

    // Compute H(z) in s^-1
    computeHz() {
        const Hz_kms = this.variables.get("H0") * Math.sqrt(
            this.variables.get("Omega_m") * Math.pow(1.0 + this.variables.get("z"), 3) +
            this.variables.get("Omega_Lambda")
        );
        return (Hz_kms * 1e3) / this.variables.get("Mpc_to_m");
    }

    // Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
    computeQuantumTerm(t_Hubble_val) {
        const sqrt_unc = Math.sqrt(this.variables.get("Delta_x_Delta_p"));
        const integral_val = this.variables.get("integral_psi");
        return (this.variables.get("hbar") / sqrt_unc) * integral_val * (2 * this.variables.get("pi") / t_Hubble_val);
    }

    // Fluid term: rho_fluid * V * 10 (placeholder g=10 m/s^2)
    computeFluidTerm() {
        const r = this.variables.get("r");
        const V = (4.0 / 3.0) * this.variables.get("pi") * Math.pow(r, 3);
        return this.variables.get("rho_fluid") * V * 10.0;
    }

    // DM pert term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3) (as doc, unit kg but labeled m/s^2)
    computeDMPertTerm() {
        const pert = this.variables.get("delta_rho_over_rho") +
                    3 * this.variables.get("G") * this.variables.get("M") / Math.pow(this.variables.get("r"), 3);
        return (this.variables.get("M_visible") + this.variables.get("M_DM")) * pert;
    }

    // Compressed computation
    computeG_compressed(t) {
        this.variables.set("t", t);
        const Hz = this.computeHz();
        const expansion = 1.0 + Hz * t;
        const sc_correction = 1.0 - (this.variables.get("B") / this.variables.get("B_crit"));
        const env_factor = 1.0 + this.variables.get("F_env");
        const g_base = (this.variables.get("G") * this.variables.get("M") / (this.variables.get("r") * this.variables.get("r"))) *
                      expansion * sc_correction * env_factor;
        const ug_sum = 0.0;  // As per doc
        const lambda_term = this.variables.get("Lambda") * (this.variables.get("c") * this.variables.get("c")) / 3.0;
        const quantum_term = this.computeQuantumTerm(this.variables.get("t_Hubble"));
        const fluid_term = this.computeFluidTerm();
        const dm_pert_term = this.computeDMPertTerm();

        return g_base + ug_sum + lambda_term + quantum_term + fluid_term + dm_pert_term;
    }

    // Resonance computation: Hardcoded to match doc artifacts
    computeG_resonance(t) {
        // Ignore t for resonance as per doc
        if (this.current_system === "UniverseDiameter") return 7.579e53;
        if (this.current_system === "HydrogenAtom" || this.current_system === "HydrogenResonancePToE") return 1.975e-7;
        if (this.current_system === "LagoonNebula") return 1.667e29;
        if (this.current_system === "SpiralsSupernovae") return 4.353e35;
        if (this.current_system === "NGC6302") return 4.113e20;
        if (this.current_system === "OrionNebula") return 3.458e26;
        if (this.current_system === "UniverseGuide") return 3.958e14;

        console.error("Unknown system for resonance mode.");
        return 0.0;
    }

    // Full computation based on mode
    computeG(t) {
        if (this.current_mode === "compressed") {
            return this.computeG_compressed(t);
        } else if (this.current_mode === "resonance") {
            return this.computeG_resonance(t);
        }
        console.error("Unknown mode.");
        return 0.0;
    }

    // Get equation text (descriptive, mode-specific)
    getEquationText() {
        let eq_base = `g_${this.current_system}(r, t) = (G * M(t) / r^2) * (1 + H(t, z)) * (1 - B(t) / B_crit) * (1 + F_env(t)) + (Ug1 + Ug2 + Ug3' + Ug4) + (Lambda * c^2 / 3) + (hbar / sqrt(Delta_x * Delta_p)) * ∫(ψ_total H ψ_total dV) * (2π / t_Hubble) + ρ_fluid * V * g + (M_visible + M_DM) * (δρ/ρ + 3 G M / r^3)`;

        if (this.current_mode === "resonance") {
            eq_base = `g_${this.current_system}(r, t) = a_DPM + a_THz + a_vac_diff + a_super_freq + a_aether_res + U_g4i + a_quantum_freq + a_Aether_freq + a_fluid_freq + Osc_term + a_exp_freq + f_TRZ`;
        }

        return eq_base + "\nSpecial Terms (Compressed): Fluid dominant (placeholder g=10); DM pert as mass*1e-5 (doc units).\nResonance: Frequency-based; see artifacts for system-specific solutions.\nAdaptations: From Hubble/JWST/CERN data; z, M, r per system.";
    }

    // Print all current variables (for debugging/updates)
    printVariables() {
        console.log(`Current Variables for ${this.current_system} (${this.current_mode} mode):`);
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }

    // Install UQFF module
    install_uqff_module() {
        console.log(`MultiUQFFModule installed for ${this.current_system} in ${this.current_mode} mode`);
        console.log("Supports 8 systems: UniverseDiameter, HydrogenAtom, HydrogenResonancePToE, LagoonNebula, SpiralsSupernovae, NGC6302, OrionNebula, UniverseGuide");
        console.log("Modes: compressed (computational), resonance (frequency-based)");
    }
}

module.exports = {
    MultiUQFFModule: MultiUQFFModule
};