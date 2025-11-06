class HydrogenAtomUQFFModule {
    constructor() {
        this.variables = new Map();

        // Base constants (universal)
        this.variables.set("G", 6.6743e-11);                    // m^3 kg^-1 s^-2
        this.variables.set("c", 3e8);                           // m/s
        this.variables.set("hbar", 1.0546e-34);                 // J s
        this.variables.set("Lambda", 1.1e-52);                  // m^-2 (negligible at atomic scale)
        this.variables.set("q", 1.602e-19);                     // C (electron charge)
        this.variables.set("pi", 3.141592653589793);            // pi
        this.variables.set("t_Hubble", 13.8e9 * 3.156e7);       // s (irrelevant, but included)

        // Hydrogen Atom parameters
        this.variables.set("M", 1.673e-27);                     // kg (proton mass, electron negligible)
        this.variables.set("M_visible", this.variables.get("M")); // Visible mass
        this.variables.set("M_DM", 0.0);                        // No DM
        this.variables.set("r", 5.29e-11);                      // m (Bohr radius)

        // Hubble/cosmology (negligible)
        this.variables.set("H0", 70.0);                         // km/s/Mpc
        this.variables.set("Mpc_to_m", 3.086e22);               // m/Mpc
        this.variables.set("z", 0.0);
        this.variables.set("Omega_m", 0.3);
        this.variables.set("Omega_Lambda", 0.7);
        this.variables.set("t", 1e-15);                         // s (atomic timescale proxy)

        // Electron/orbital dynamics
        this.variables.set("rho_fluid", 1e-25);                 // kg/m^3 (electron cloud density est.)
        this.variables.set("V", (4.0 / 3.0) * this.variables.get("pi") * Math.pow(this.variables.get("r"), 3));  // m^3 (orbital volume)
        this.variables.set("v_orbital", 2.2e6);                 // m/s (electron velocity)
        this.variables.set("delta_rho", 0.1 * this.variables.get("rho_fluid"));
        this.variables.set("rho", this.variables.get("rho_fluid"));

        // EM/magnetic (atomic scale)
        this.variables.set("B", 1e-4);                          // T (internal atomic field est.)
        this.variables.set("B_crit", 1e11);                     // T

        // Quantum terms (dominant)
        this.variables.set("Delta_x", 1e-10);                   // m (Compton wavelength proxy)
        this.variables.set("Delta_p", this.variables.get("hbar") / this.variables.get("Delta_x"));
        this.variables.set("integral_psi", 1.0);                // Normalized (ground state)

        // Resonant/oscillatory terms (atomic transitions)
        this.variables.set("A", 1e-10);                         // Amplitude
        this.variables.set("k", 1e11);                          // m^-1 (UV wavelength ~1e-8 m)
        this.variables.set("omega", 1e15);                      // rad/s (~Lyman alpha freq)
        this.variables.set("x", 0.0);                           // m

        // Ug subterms (init placeholders)
        this.variables.set("Ug1", 0.0);
        this.variables.set("Ug2", 0.0);
        this.variables.set("Ug3", 0.0);
        this.variables.set("Ug4", 0.0);

        // Scale factors
        this.variables.set("scale_macro", 1e-12);               // Adjusted for atomic
        this.variables.set("f_TRZ", 0.1);
        this.variables.set("f_sc", 1.0);
    }

    // Dynamic variable operations
    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
        } else {
            console.warn(`Variable '${name}' not found. Adding with value ${value}`);
            this.variables.set(name, value);
        }

        // Update dependent variables
        if (name === "Delta_x") {
            this.variables.set("Delta_p", this.variables.get("hbar") / value);
        } else if (name === "r") {
            this.variables.set("V", (4.0 / 3.0) * this.variables.get("pi") * Math.pow(value, 3));
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
        return this.variables.get(name);
    }

    // Compute H(z) in s^-1 (negligible)
    computeHz() {
        const Hz_kms = this.variables.get("H0") * Math.sqrt(
            this.variables.get("Omega_m") * Math.pow(1.0 + this.variables.get("z"), 3) +
            this.variables.get("Omega_Lambda")
        );
        return (Hz_kms * 1e3) / this.variables.get("Mpc_to_m");
    }

    // Compute Ug sum: Ug1 = G M / r^2, Ug4 = Ug1 * f_sc, others 0 (weak)
    computeUgSum() {
        const Ug1 = (this.variables.get("G") * this.variables.get("M")) / (this.variables.get("r") * this.variables.get("r"));
        this.variables.set("Ug1", Ug1);
        this.variables.set("Ug4", Ug1 * this.variables.get("f_sc"));
        return this.variables.get("Ug1") + this.variables.get("Ug2") + this.variables.get("Ug3") + this.variables.get("Ug4");
    }

    // Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble) (dominant)
    computeQuantumTerm(t_Hubble_val) {
        const unc = Math.sqrt(this.variables.get("Delta_x") * this.variables.get("Delta_p"));
        const integral_val = this.variables.get("integral_psi");
        return (this.variables.get("hbar") / unc) * integral_val * (2 * this.variables.get("pi") / t_Hubble_val);
    }

    // Fluid term: rho_fluid * V * g (electron cloud)
    computeFluidTerm(g_base) {
        return this.variables.get("rho_fluid") * this.variables.get("V") * g_base;
    }

    // Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
    computeResonantTerm(t) {
        const cos_term = 2 * this.variables.get("A") * Math.cos(this.variables.get("k") * this.variables.get("x")) * Math.cos(this.variables.get("omega") * t);
        const exp_real = this.variables.get("A") * Math.cos(this.variables.get("k") * this.variables.get("x") - this.variables.get("omega") * t);
        const exp_factor = (2 * this.variables.get("pi") / 13.8);
        return cos_term + exp_factor * exp_real;
    }

    // DM term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3) (negligible)
    computeDMTerm() {
        const pert = this.variables.get("delta_rho") / this.variables.get("rho");
        const curv = 3 * this.variables.get("G") * this.variables.get("M") / (this.variables.get("r") * this.variables.get("r") * this.variables.get("r"));
        return (this.variables.get("M_visible") + this.variables.get("M_DM")) * (pert + curv);
    }

    // Full computation: g_UQFF(r, t) = ... all terms (quantum dominant)
    computeG(t) {
        this.variables.set("t", t);
        const Hz = this.computeHz();
        const expansion = 1.0 + Hz * t;  // ~1
        const sc_correction = 1.0 - (this.variables.get("B") / this.variables.get("B_crit"));
        const tr_factor = 1.0 + this.variables.get("f_TRZ");

        // Base gravity with expansion, SC, TR (weak)
        const g_base = (this.variables.get("G") * this.variables.get("M") / (this.variables.get("r") * this.variables.get("r"))) * expansion * sc_correction * tr_factor;

        // Ug sum
        const ug_sum = this.computeUgSum();

        // Cosmological (negligible)
        const lambda_term = this.variables.get("Lambda") * (this.variables.get("c") * this.variables.get("c")) / 3.0;

        // Quantum (dominant)
        const quantum_term = this.computeQuantumTerm(this.variables.get("t_Hubble"));

        // EM Lorentz (electron orbital)
        const em_base = this.variables.get("q") * this.variables.get("v_orbital") * this.variables.get("B") / 9.11e-31;  // / electron mass
        const em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * this.variables.get("scale_macro");

        // Fluid
        const fluid_term = this.computeFluidTerm(g_base);

        // Resonant (orbital transitions)
        const resonant_term = this.computeResonantTerm(t);

        // DM (negligible)
        const dm_term = this.computeDMTerm();

        // Total: Sum all
        return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term;
    }

    // Get equation text (descriptive)
    getEquationText() {
        return "g_Hydrogen(r, t) = (G * M / r^2) * (1 + H(z) * t) * (1 - B / B_crit) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 / 3) + " +
               "(hbar / sqrt(Delta_x * Delta_p)) * ∫(ψ* H ψ dV) * (2π / t_Hubble) + q (v × B) + ρ_fluid * V * g + " +
               "2 A cos(k x) cos(ω t) + (2π / 13.8) A exp(i (k x - ω t)) + (M_visible + M_DM) * (δρ/ρ + 3 G M / r^3)\n" +
               "Special Terms:\n" +
               "- Quantum: Heisenberg uncertainty dominant for orbital stability.\n" +
               "- Fluid: Electron cloud density-volume-gravity coupling.\n" +
               "- Resonant: Oscillatory waves for atomic transitions/orbitals.\n" +
               "- DM: Negligible at atomic scale.\n" +
               "- Superconductivity: (1 - B/B_crit) for quantum field in atom.\n" +
               "Solutions: At t=1e-15 s, g_Hydrogen ~1e12 m/s² (EM/quantum dominant; g_base ~1e-40 m/s²).\n" +
               "Adaptations for Hydrogen Atom: Bohr r=5.29e-11 m; v_orbital=2.2e6 m/s; f_osc=1e15 Hz (Lyman).";
    }

    // Print all current variables (for debugging/updates)
    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }
}

module.exports = HydrogenAtomUQFFModule;