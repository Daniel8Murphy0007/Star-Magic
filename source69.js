class UQFFCompressionModule {
    constructor() {
        this.variables = new Map();
        this.current_system = "General";
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;
        this.metadata = new Map();
        this.metadata.set("enhanced", "true");
        this.metadata.set("version", "2.0-Enhanced");

        // Universal constants
        this.variables.set("G", 6.6743e-11);                    // m^3 kg^-1 s^-2
        this.variables.set("c", 3e8);                           // m/s
        this.variables.set("hbar", 1.0546e-34);                 // J s
        this.variables.set("Lambda", 1.1e-52);                  // m^-2
        this.variables.set("q", 1.602e-19);                     // C
        this.variables.set("pi", 3.141592653589793);            // pi
        this.variables.set("t_Hubble", 13.8e9 * 3.156e7);       // s
        this.variables.set("year_to_s", 3.156e7);               // s/yr
        this.variables.set("H0", 67.15);                        // km/s/Mpc
        this.variables.set("Mpc_to_m", 3.086e22);               // m/Mpc
        this.variables.set("Omega_m", 0.3);
        this.variables.set("Omega_Lambda", 0.7);

        // General system params (overridden by setSystem)
        this.variables.set("M", 1e33);                          // kg (placeholder)
        this.variables.set("M0", this.variables.get("M"));
        this.variables.set("SFR", 1e30);                        // kg/yr (Msun/yr equiv)
        this.variables.set("r", 1e17);                          // m
        this.variables.set("z", 0.001);                         // Redshift (general)
        this.variables.set("M_visible", 0.15 * this.variables.get("M"));
        this.variables.set("M_DM", 0.85 * this.variables.get("M"));
        this.variables.set("t", 1e6 * this.variables.get("year_to_s"));  // Default 1 Myr

        // Dynamics
        this.variables.set("rho_fluid", 1e-20);                 // kg/m^3
        this.variables.set("V", 1e3);                           // m^3
        this.variables.set("B", 1e-5);                          // T
        this.variables.set("B_crit", 1e11);                     // T
        this.variables.set("Delta_x", 1e-10);                   // m
        this.variables.set("Delta_p", this.variables.get("hbar") / this.variables.get("Delta_x"));
        this.variables.set("integral_psi", 1.0);                // Normalized

        // Wave/oscillatory
        this.variables.set("A", 1e-10);
        this.variables.set("k", 1e20);
        this.variables.set("omega", 1e15);                      // rad/s
        this.variables.set("x", 0.0);
        this.variables.set("v", 1e3);                           // m/s for v x B

        // Ug subterms (initial)
        this.variables.set("Ug1", 0.0);
        this.variables.set("Ug2", 0.0);
        this.variables.set("Ug3", 0.0);                         // Will be Ug3'
        this.variables.set("Ug4", 0.0);

        // Environmental sub-terms (for F_env)
        this.variables.set("F_wind", 0.0);
        this.variables.set("F_rad", 0.0);
        this.variables.set("F_SN", 0.0);
        this.variables.set("F_BH", 0.0);
        this.variables.set("F_erode", 0.0);
        this.variables.set("F_lensing", 0.0);
        this.variables.set("F_mag", 0.0);
        this.variables.set("F_decay", 0.0);
        this.variables.set("F_coll", 0.0);
        this.variables.set("F_evo", 0.0);
        this.variables.set("F_merge", 0.0);
        this.variables.set("F_sf", 0.0);

        // External for Ug3'
        this.variables.set("M_ext", 0.0);
        this.variables.set("r_ext", 1e18);                      // m

        // Scales
        this.variables.set("scale_macro", 1e-12);
        this.variables.set("f_TRZ", 0.1);
        this.variables.set("f_sc", 1.0);
        this.variables.set("delta_rho", 0.1 * this.variables.get("rho_fluid"));
        this.variables.set("rho", this.variables.get("rho_fluid"));
    }

    setSystem(sys_name) {
        this.current_system = sys_name;
        if (sys_name === "Magnetar") {
            this.variables.set("M", 2e30);  // kg
            this.variables.set("r", 1e4);   // m
            this.variables.set("B", 1e11);  // T
            this.variables.set("M_ext", 4e6 * 1.989e30);  // Sgr A* mass
            this.variables.set("r_ext", 8e19);  // ~kpc
            this.variables.set("F_mag", 1e40);  // Magnetic energy
            this.variables.set("F_decay", Math.exp(-this.variables.get("t") / 1e6));  // Outburst decay
        } else if (sys_name === "SagittariusA") {
            this.variables.set("M", 4e6 * 1.989e30);  // kg
            this.variables.set("r", 1e10);  // m (event horizon scale)
            this.variables.set("z", 0.00034);
            this.variables.set("F_BH", 1e42);  // GW term approx
        } else if (sys_name === "Pillars") {
            this.variables.set("M", 1e33);
            this.variables.set("r", 1e17);
            this.variables.set("F_erode", 0.1 * (this.variables.get("t") / (3e5 * this.variables.get("year_to_s"))));
            this.variables.set("F_wind", this.variables.get("rho_fluid") * Math.pow(1e4, 2));  // v_wind ~10 km/s
        }
        // Auto-update dependents
        this.variables.set("M_visible", 0.15 * this.variables.get("M"));
        this.variables.set("M_DM", 0.85 * this.variables.get("M"));
        this.variables.set("M0", this.variables.get("M"));
        this.variables.set("Delta_p", this.variables.get("hbar") / this.variables.get("Delta_x"));
    }

    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
        } else {
            console.log(`Variable '${name}' not found. Adding.`);
            this.variables.set(name, value);
        }
        if (name === "Delta_x") {
            this.variables.set("Delta_p", this.variables.get("hbar") / value);
        } else if (name === "M") {
            this.variables.set("M_visible", 0.15 * value);
            this.variables.set("M_DM", 0.85 * value);
            this.variables.set("M0", value);
        }
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
        } else {
            this.variables.set(name, delta);
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    computeHtz(z_val) {
        const Hz_kms = this.variables.get("H0") * Math.sqrt(
            this.variables.get("Omega_m") * Math.pow(1.0 + z_val, 3) + this.variables.get("Omega_Lambda")
        );
        return (Hz_kms * 1e3) / this.variables.get("Mpc_to_m");
    }

    computeFenv(t) {
        const f_wind = this.variables.get("F_wind");
        const f_erode = this.variables.get("F_erode");
        const f_lensing = this.variables.get("F_lensing");
        const f_mag = this.variables.get("F_mag");
        const f_decay = this.variables.get("F_decay");
        const f_coll = this.variables.get("F_coll");
        const f_evo = this.variables.get("F_evo");
        const f_merge = this.variables.get("F_merge");
        const f_sf = this.variables.get("F_sf");
        const f_sn = this.variables.get("F_SN");
        const f_rad = this.variables.get("F_rad");
        const f_bh = this.variables.get("F_BH");
        return f_wind + f_erode + f_lensing + f_mag + f_decay + f_coll + f_evo + f_merge + f_sf + f_sn + f_rad + f_bh;
    }

    computeUg3prime() {
        return (this.variables.get("G") * this.variables.get("M_ext")) / (this.variables.get("r_ext") * this.variables.get("r_ext"));
    }

    computePsiTotal(t) {
        const mag_term = this.variables.get("q") * this.variables.get("v") * this.variables.get("B");
        const standing = 2 * this.variables.get("A") * Math.cos(this.variables.get("k") * this.variables.get("x")) * Math.cos(this.variables.get("omega") * t);
        const exp_term_real = this.variables.get("A") * Math.cos(this.variables.get("k") * this.variables.get("x") - this.variables.get("omega") * t);
        const quantum_wave = (2 * this.variables.get("pi") / 13.8) * exp_term_real;
        return mag_term + standing + quantum_wave;
    }

    computeQuantumTerm(t_Hubble_val) {
        const unc = Math.sqrt(this.variables.get("Delta_x") * this.variables.get("Delta_p"));
        const psi_tot = this.computePsiTotal(this.variables.get("t"));
        return (this.variables.get("hbar") / unc) * this.variables.get("integral_psi") * (2 * this.variables.get("pi") / t_Hubble_val) * psi_tot;
    }

    computeFluidTerm(g_base) {
        return this.variables.get("rho_fluid") * this.variables.get("V") * g_base;
    }

    computeDMTerm() {
        const pert = this.variables.get("delta_rho") / this.variables.get("rho");
        const curv = 3 * this.variables.get("G") * this.variables.get("M") / (this.variables.get("r") * this.variables.get("r") * this.variables.get("r"));
        return (this.variables.get("M_visible") + this.variables.get("M_DM")) * (pert + curv);
    }

    computeUgSum() {
        const Ug1 = (this.variables.get("G") * this.variables.get("M")) / (this.variables.get("r") * this.variables.get("r"));
        this.variables.set("Ug1", Ug1);
        this.variables.set("Ug3", this.computeUg3prime());
        this.variables.set("Ug4", Ug1 * this.variables.get("f_sc"));
        return this.variables.get("Ug1") + this.variables.get("Ug2") + this.variables.get("Ug3") + this.variables.get("Ug4");
    }

    computeMsfFactor(t) {
        const t_yr = t / this.variables.get("year_to_s");
        return (this.variables.get("SFR") * t_yr) / this.variables.get("M0");
    }

    computeG(t) {
        this.variables.set("t", t);
        const Htz = this.computeHtz(this.variables.get("z"));
        const expansion = 1.0 + Htz * t;
        const sc_correction = 1.0 - (this.variables.get("B") / this.variables.get("B_crit"));
        const msf_factor = this.computeMsfFactor(t);
        const m_factor = 1.0 + msf_factor;
        const f_env = this.computeFenv(t);

        // Base gravity with corrections and M(t)
        const g_base = (this.variables.get("G") * this.variables.get("M") * m_factor / (this.variables.get("r") * this.variables.get("r"))) * expansion * sc_correction * (1.0 + f_env);

        // Ug sum
        const ug_sum = this.computeUgSum();

        // Cosmological
        const lambda_term = this.variables.get("Lambda") * (this.variables.get("c") * this.variables.get("c")) / 3.0;

        // Quantum
        const quantum_term = this.computeQuantumTerm(this.variables.get("t_Hubble"));

        // Fluid
        const fluid_term = this.computeFluidTerm(g_base);

        // DM
        const dm_term = this.computeDMTerm();

        // Total
        return g_base + ug_sum + lambda_term + quantum_term + fluid_term + dm_term;
    }

    getEquationText() {
        return "g_UQFF(r, t) = (G * M(t) / r^2) * (1 + H(t, z)) * (1 - B(t) / B_crit) * (1 + F_env(t)) + (Ug1 + Ug2 + Ug3' + Ug4) + (Lambda * c^2 / 3) + " +
               "(hbar / sqrt(Delta_x * Delta_p)) * ∫(psi_total * H * psi_total dV) * (2π / t_Hubble) + ρ_fluid * V * g + " +
               "(M_visible + M_DM) * (δρ/ρ + 3 G M / r^3)\n" +
               "Where: H(t, z) = H0 * sqrt(Ωm (1+z)^3 + ΩΛ); M(t) = M * (1 + M_sf(t)); M_sf(t) = (SFR * t_yr) / M0;\n" +
               "F_env(t) = ∑ F_i(t) [winds, erosion, lensing, mag, decay, coll, evo, merge, sf, SN, rad, BH];\n" +
               "Ug3' = G M_ext / r_ext^2; psi_total = q(v × B) + 2A cos(kx) cos(ωt) + (2π/13.8) A Re[exp(i(kx - ωt))];\n" +
               "Compression Advancements: Unified expansion, modular env effects, consolidated waves/gravity terms for 19+ systems.\n" +
               "Adaptations: setSystem('Magnetar') for SGR 1745-2900; etc. Solutions: g ~1e-10 to 1e-12 m/s² typical.";
    }

    printVariables() {
        console.log("System:", this.current_system);
        console.log("Current Variables:");
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }
}

module.exports = { UQFFCompressionModule };