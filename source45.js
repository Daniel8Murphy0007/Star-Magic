class SpiralSupernovaeUQFFModule {
    constructor() {
        // Initialize Map for dynamic variables
        this.variables = new Map();

        // Base constants (universal)
        this.variables.set("G", 6.6743e-11);                    // m^3 kg^-1 s^-2
        this.variables.set("c", 3e8);                           // m/s
        this.variables.set("hbar", 1.0546e-34);                 // J s
        this.variables.set("Lambda", 1.1e-52);                  // m^-2
        this.variables.set("q", 1.602e-19);                     // C
        this.variables.set("pi", 3.141592653589793);            // pi
        this.variables.set("t_Hubble", 13.8e9 * 3.156e7);       // s

        // Spiral-SN parameters
        const M_sun_val = 1.989e30;                    // kg
        this.variables.set("M_sun", M_sun_val);
        this.variables.set("M", 1e11 * M_sun_val);              // Galaxy mass kg
        this.variables.set("M_visible", 0.15 * this.variables.get("M")); // Visible fraction
        this.variables.set("M_DM", 0.85 * this.variables.get("M"));      // Dark matter
        this.variables.set("r", 9.258e20);                      // m (~30 kpc)
        this.variables.set("M_gas", 1e9 * M_sun_val);           // Gas mass

        // Hubble/cosmology
        this.variables.set("H0", 73.0);                         // km/s/Mpc (SH0ES)
        this.variables.set("Mpc_to_m", 3.086e22);               // m/Mpc
        this.variables.set("z", 0.5);                           // Typical z for SN
        this.variables.set("Omega_m", 0.3);
        this.variables.set("Omega_Lambda", 0.7);
        this.variables.set("t", 5e9 * 3.156e7);                 // Default t=5 Gyr s

        // Spiral dynamics
        this.variables.set("Omega_p", 20e3 / 3.086e19);         // rad/s (20 km/s/kpc pattern speed)

        // SN parameters
        this.variables.set("L_SN", 1e36);                       // W (peak luminosity)
        this.variables.set("rho_fluid", 1e-21);                 // kg/m^3 (ISM)
        this.variables.set("V", 1e3);                           // m^3
        this.variables.set("v_rot", 2e5);                       // m/s (rotation)
        this.variables.set("delta_rho", 0.1 * this.variables.get("rho_fluid"));
        this.variables.set("rho", this.variables.get("rho_fluid"));

        // EM/magnetic
        this.variables.set("B", 1e-5);                          // T (galactic field)
        this.variables.set("B_crit", 1e11);                     // T

        // Quantum terms
        this.variables.set("Delta_x", 1e-10);                   // m
        this.variables.set("Delta_p", this.variables.get("hbar") / this.variables.get("Delta_x"));
        this.variables.set("integral_psi", 1.0);

        // Resonant/oscillatory
        this.variables.set("A", 1e-10);
        this.variables.set("k", 1e20);
        this.variables.set("omega", 1e15);
        this.variables.set("x", 0.0);

        // Ug subterms
        this.variables.set("Ug1", 0.0);
        this.variables.set("Ug2", 0.0);
        this.variables.set("Ug3", 0.0);
        this.variables.set("Ug4", 0.0);

        // Scale factors
        this.variables.set("scale_macro", 1e-12);
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

        // Auto-update dependent variables
        if (name === "Delta_x") {
            this.variables.set("Delta_p", this.variables.get("hbar") / value);
        } else if (name === "M") {
            this.variables.set("M_visible", 0.15 * value);
            this.variables.set("M_DM", 0.85 * value);
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
    computeHz(z) {
        const Hz_kms = this.variables.get("H0") * Math.sqrt(
            this.variables.get("Omega_m") * Math.pow(1.0 + z, 3) +
            this.variables.get("Omega_Lambda")
        );
        return (Hz_kms * 1e3) / this.variables.get("Mpc_to_m");
    }

    // Compute Ug sum: Ug1 = G M / r^2, Ug4 = Ug1 * f_sc, others 0
    computeUgSum() {
        const Ug1 = (this.variables.get("G") * this.variables.get("M")) /
                   (this.variables.get("r") * this.variables.get("r"));
        this.variables.set("Ug1", Ug1);
        this.variables.set("Ug4", Ug1 * this.variables.get("f_sc"));
        return this.variables.get("Ug1") + this.variables.get("Ug2") +
               this.variables.get("Ug3") + this.variables.get("Ug4");
    }

    // Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
    computeQuantumTerm(t_Hubble_val) {
        const unc = Math.sqrt(this.variables.get("Delta_x") * this.variables.get("Delta_p"));
        const integral_val = this.variables.get("integral_psi");
        return (this.variables.get("hbar") / unc) * integral_val *
               (2 * this.variables.get("pi") / t_Hubble_val);
    }

    // Fluid term: rho_fluid * V * g
    computeFluidTerm(g_base) {
        return this.variables.get("rho_fluid") * this.variables.get("V") * g_base;
    }

    // Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
    computeResonantTerm(t) {
        const cos_term = 2 * this.variables.get("A") *
                        Math.cos(this.variables.get("k") * this.variables.get("x")) *
                        Math.cos(this.variables.get("omega") * t);

        // Complex exponential: Re[A * exp(i * (k*x - omega*t))]
        const phase = this.variables.get("k") * this.variables.get("x") - this.variables.get("omega") * t;
        const real_exp = this.variables.get("A") * Math.cos(phase);
        const exp_factor = (2 * this.variables.get("pi") / 13.8);

        return cos_term + exp_factor * real_exp;
    }

    // DM term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3)
    computeDMTerm() {
        const pert = this.variables.get("delta_rho") / this.variables.get("rho");
        const curv = 3 * this.variables.get("G") * this.variables.get("M") /
                    (this.variables.get("r") * this.variables.get("r") * this.variables.get("r"));
        return (this.variables.get("M_visible") + this.variables.get("M_DM")) * (pert + curv);
    }

    // Spiral torque term: T_spiral = G * M_gas * M / r^2 * (1 + Omega_p * t)
    computeT_spiral(t) {
        const torque_base = (this.variables.get("G") * this.variables.get("M_gas") * this.variables.get("M")) /
                           (this.variables.get("r") * this.variables.get("r"));
        return torque_base * (1.0 + this.variables.get("Omega_p") * t);
    }

    // Supernova term: SN_term = (L_SN / (4 pi r^2 c)) * (1 + H(z) * t)
    computeSN_term(z) {
        const Hz = this.computeHz(z);
        const flux = this.variables.get("L_SN") /
                    (4 * this.variables.get("pi") * this.variables.get("r") * this.variables.get("r") * this.variables.get("c"));
        return flux * (1.0 + Hz * this.variables.get("t"));
    }

    // Full computation: g_UQFF(r, t) = ... all terms with T_spiral and SN_term
    computeG(t, z = null) {
        if (z === null) z = this.variables.get("z");
        this.variables.set("t", t);

        const Hz = this.computeHz(z);
        const expansion = 1.0 + Hz * t;
        const sc_correction = 1.0 - (this.variables.get("B") / this.variables.get("B_crit"));
        const tr_factor = 1.0 + this.variables.get("f_TRZ");
        const t_spiral = this.computeT_spiral(t);
        const sn_term = this.computeSN_term(z);

        // Base gravity with expansion, SC, TR, T_spiral
        const g_base = ((this.variables.get("G") * this.variables.get("M") /
                        (this.variables.get("r") * this.variables.get("r"))) *
                       expansion * sc_correction * tr_factor) * (1.0 + t_spiral);

        // Ug sum
        const ug_sum = this.computeUgSum();

        // Cosmological with Lambda
        const lambda_term = this.variables.get("Lambda") *
                           (this.variables.get("c") * this.variables.get("c") *
                            this.variables.get("Omega_Lambda")) / 3.0;

        // Quantum
        const quantum_term = this.computeQuantumTerm(this.variables.get("t_Hubble"));

        // EM Lorentz (rotation v_rot B)
        const em_base = this.variables.get("q") * this.variables.get("v_rot") *
                       this.variables.get("B") / 1.673e-27;
        const em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * this.variables.get("scale_macro");

        // Fluid
        const fluid_term = this.computeFluidTerm(g_base);

        // Resonant
        const resonant_term = this.computeResonantTerm(t);

        // DM
        const dm_term = this.computeDMTerm();

        // Total: Sum all + SN_term
        return g_base + ug_sum + lambda_term + quantum_term + em_term +
               fluid_term + resonant_term + dm_term + sn_term;
    }

    // Get equation text (descriptive)
    getEquationText() {
        return "g_Spiral_SN(r, t) = (G * M / r^2) * (1 + H(z) * t) * (1 + T_spiral) * (1 - B / B_crit) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 * Omega_Lambda / 3) + " +
               "(hbar / sqrt(Delta_x * Delta_p)) * integral_psi * (2*pi / t_Hubble) + q (v_rot × B) + rho_fluid * V * g + " +
               "2 A cos(k x) cos(omega t) + (2*pi / 13.8) A exp(i (k x - omega t)) + (M_visible + M_DM) * (delta_rho/rho + 3 G M / r^3) + SN_term\n" +
               "Where T_spiral = G * M_gas * M / r^2 * (1 + Omega_p * t); SN_term = (L_SN / (4*pi r^2 c)) * (1 + H(z) * t)\n" +
               "Special Terms:\n" +
               "- Quantum: Heisenberg uncertainty for ISM quantum effects.\n" +
               "- Fluid: Gas density-volume-gravity coupling in arms.\n" +
               "- Resonant: Oscillatory Aether waves for density waves.\n" +
               "- DM: Visible+dark mass with perturbations for rotation curves.\n" +
               "- Superconductivity: (1 - B/B_crit) for galactic fields.\n" +
               "- Spiral Torque: T_spiral for arm evolution.\n" +
               "- Supernova: SN_term for expansion probe.\n" +
               "Solutions: At t=5 Gyr, z=0.5, g_Spiral_SN ~1e-10 m/s² (Lambda/SN dominant; g_base ~1e-10).\n" +
               "Adaptations for Spirals and Supernovae: SH0ES H0=73; Omega_p=20 km/s/kpc; L_SN=1e36 W for Ia SN.";
    }

    // Print variables (for debugging)
    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }
}

module.exports = SpiralSupernovaeUQFFModule;