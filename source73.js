class NGC1300EnhancedUQFFModule {
    constructor() {
        // Self-expanding framework members
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;
        this.metadata.set("enhanced", "true");
        this.metadata.set("version", "2.0-Enhanced");

        // Core variables (NGC 1300 parameters)
        this.variables = new Map();

        // Universal constants
        this.variables.set("G", 6.6743e-11);                    // m^3 kg^-1 s^-2
        this.variables.set("c", 3e8);                           // m/s
        this.variables.set("hbar", 1.0546e-34);                 // J s
        this.variables.set("Lambda", 1.1e-52);                  // m^-2
        this.variables.set("q", 1.602e-19);                     // C
        this.variables.set("pi", 3.141592653589793);            // pi
        this.variables.set("t_Hubble", 13.8e9 * 3.156e7);       // s
        this.variables.set("year_to_s", 3.156e7);               // s/yr
        this.variables.set("H0", 70.0);                         // km/s/Mpc
        this.variables.set("Mpc_to_m", 3.086e22);               // m/Mpc
        this.variables.set("Omega_m", 0.3);
        this.variables.set("Omega_Lambda", 0.7);

        const M_sun_val = 1.989e30;                             // kg
        const kpc_val = 3.086e19;                               // m

        // NGC 1300 parameters
        this.variables.set("M_visible", 7e10 * M_sun_val);      // kg
        this.variables.set("M_DM", 3e10 * M_sun_val);           // kg
        this.variables.set("M", this.variables.get("M_visible") + this.variables.get("M_DM"));  // Total initial
        this.variables.set("M0", this.variables.get("M"));
        this.variables.set("SFR", 1 * M_sun_val / this.variables.get("year_to_s"));    // kg/s
        this.variables.set("r", 11.79e3 * kpc_val);             // m
        this.variables.set("z", 0.005);                         // Redshift
        this.variables.set("v_arm", 200e3);                     // m/s (gas velocity)
        this.variables.set("t", 1e9 * this.variables.get("year_to_s"));  // Default t=1 Gyr s

        // Dynamics
        this.variables.set("rho_fluid", 1e-21);                 // kg/m^3
        this.variables.set("V", 1e50);                          // m^3
        this.variables.set("B", 1e-5);                          // T
        this.variables.set("B_crit", 1e11);                     // T
        this.variables.set("Delta_x", 1e-10);                   // m
        this.variables.set("Delta_p", this.variables.get("hbar") / this.variables.get("Delta_x"));
        this.variables.set("integral_psi", 1.0);                // Normalized

        // Wave/oscillatory for spiral arms
        this.variables.set("A", 1e-10);
        this.variables.set("k", 1e20);
        this.variables.set("omega", 1e-15);                     // rad/s for density waves
        this.variables.set("x", 0.0);
        this.variables.set("v", this.variables.get("v_arm"));   // m/s
        this.variables.set("sigma", 1e3 * kpc_val);             // m for Gaussian

        // Ug subterms & Ui
        this.variables.set("Ug1", 0.0);                         // Dipole
        this.variables.set("Ug2", 0.0);                         // Superconductor
        this.variables.set("Ug3", 0.0);                         // External
        this.variables.set("Ug4", 0.0);                         // Reaction
        this.variables.set("Ui", 0.0);
        this.variables.set("mu_0", 4 * this.variables.get("pi") * 1e-7); // H/m
        this.variables.set("rho_vac_SCm", 7.09e-37);            // J/m^3
        this.variables.set("rho_vac_UA", 7.09e-36);             // J/m^3
        this.variables.set("lambda_I", 1.0);
        this.variables.set("omega_i", 1e-8);                    // rad/s
        this.variables.set("t_n", 0.0);
        this.variables.set("F_RZ", 0.01);
        this.variables.set("k_4", 1.0);
        this.variables.set("k_SF", 1e-10);                      // N/Msun, adjusted to m/s^2
        this.variables.set("omega_spin", 1e-4);                 // rad/s (bar rotation)
        this.variables.set("I_dipole", 1e20);                   // A
        this.variables.set("A_dipole", 1e15);                   // m^2
        this.variables.set("H_aether", 1e-6);                   // A/m
        this.variables.set("delta_rho_over_rho", 1e-5);

        // Scales
        this.variables.set("scale_macro", 1e-12);
        this.variables.set("f_TRZ", 0.1);
        this.variables.set("f_sc", 1.0);
        this.variables.set("v_r", 1e3);                         // m/s radial velocity
        this.variables.set("rho", this.variables.get("rho_fluid"));
    }

    // Dynamic variable operations
    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
        } else {
            if (this.enableLogging) {
                console.log(`Variable '${name}' not found. Adding.`);
            }
            this.variables.set(name, value);
        }

        // Update dependent variables
        if (name === "Delta_x") {
            this.variables.set("Delta_p", this.variables.get("hbar") / value);
        } else if (name === "M") {
            this.variables.set("M_visible", 0.7 * value);
            this.variables.set("M_DM", 0.3 * value);
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

    getVariable(name) {
        return this.variables.get(name);
    }

    // Compute H(t, z)
    computeHtz(z_val) {
        const Hz_kms = this.variables.get("H0") * Math.sqrt(
            this.variables.get("Omega_m") * Math.pow(1.0 + z_val, 3) +
            this.variables.get("Omega_Lambda")
        );
        return (Hz_kms * 1e3) / this.variables.get("Mpc_to_m");
    }

    // M(t) star formation factor
    computeMsfFactor(t) {
        return this.variables.get("SFR") * t / this.variables.get("M0");
    }

    // r(t)
    computeRt(t) {
        return this.variables.get("r") + this.variables.get("v_r") * t;
    }

    // F_env(t) - environmental forces
    computeFenv(t) {
        const F_bar = 0.1 * (this.variables.get("G") * this.variables.get("M")) /
                      (this.variables.get("r") * this.variables.get("r"));  // Bar funneling
        const F_SF = this.variables.get("k_SF") * this.variables.get("SFR") / 1.989e30;  // Normalize to m/s^2
        const F_wave = this.variables.get("rho_fluid") * Math.pow(this.variables.get("v_arm"), 2);  // Density wave
        return F_bar + F_SF + F_wave;
    }

    // Ug1: dipole
    computeUg1(t) {
        const mu_dipole = this.variables.get("I_dipole") * this.variables.get("A_dipole") * this.variables.get("omega_spin");
        return mu_dipole * this.variables.get("B");
    }

    // Ug2: superconductor
    computeUg2(t) {
        const B_super = this.variables.get("mu_0") * this.variables.get("H_aether");
        return (B_super * B_super) / (2 * this.variables.get("mu_0"));
    }

    // Ug3': external (bar as external)
    computeUg3prime(t) {
        const M_bar = 0.2 * this.variables.get("M");  // Bar mass fraction
        const r_bar = 0.3 * this.variables.get("r");  // Bar radius
        return (this.variables.get("G") * M_bar) / (r_bar * r_bar);
    }

    // Ug4: reaction
    computeUg4(t) {
        const E_react = 1e46 * Math.exp(-0.0005 * t);
        return this.variables.get("k_4") * E_react;
    }

    // Ui
    computeUi(t) {
        return this.variables.get("lambda_I") *
               (this.variables.get("rho_vac_SCm") / this.variables.get("rho_vac_UA")) *
               this.variables.get("omega_i") *
               Math.cos(this.variables.get("pi") * this.variables.get("t_n")) *
               (1 + this.variables.get("F_RZ"));
    }

    // Psi integral (simplified spiral wave function)
    computePsiIntegral(r, t) {
        const A = this.variables.get("A");
        const m = 2.0;  // m-mode for spiral
        const omega = this.variables.get("omega");
        const sigma = this.variables.get("sigma");

        // Complex exponential for spiral pattern
        const real_part = A * Math.exp(-r*r / (2 * sigma * sigma)) * Math.cos(m * 0 - omega * t);
        const imag_part = A * Math.exp(-r*r / (2 * sigma * sigma)) * Math.sin(m * 0 - omega * t);

        return real_part * real_part + imag_part * imag_part;  // |psi|^2
    }

    // Quantum term
    computeQuantumTerm(t_Hubble_val, r) {
        const unc = Math.sqrt(this.variables.get("Delta_x") * this.variables.get("Delta_p"));
        const psi_int = this.computePsiIntegral(r, this.variables.get("t"));
        return (this.variables.get("hbar") / unc) *
               this.variables.get("integral_psi") *
               (2 * this.variables.get("pi") / t_Hubble_val) *
               psi_int;
    }

    // Fluid term
    computeFluidTerm(g_base) {
        return this.variables.get("rho_fluid") * this.variables.get("V") * g_base;
    }

    // DM term
    computeDMTerm(r) {
        const pert = this.variables.get("delta_rho_over_rho");
        const curv = 3 * this.variables.get("G") * this.variables.get("M") / (r * r * r);
        return (this.variables.get("M_visible") + this.variables.get("M_DM")) * (pert + curv);
    }

    // Ug sum
    computeUgSum(r) {
        const Ug_base = (this.variables.get("G") * this.variables.get("M")) / (r * r);
        this.variables.set("Ug1", this.computeUg1(this.variables.get("t")));
        this.variables.set("Ug2", this.computeUg2(this.variables.get("t")));
        this.variables.set("Ug3", this.computeUg3prime(this.variables.get("t")));
        this.variables.set("Ug4", this.computeUg4(this.variables.get("t")));
        return Ug_base + this.variables.get("Ug1") + this.variables.get("Ug2") +
               this.variables.get("Ug3") + this.variables.get("Ug4");
    }

    // Full g_NGC1300 computation
    computeG(t, r) {
        this.variables.set("t", t);
        const msf_factor = this.computeMsfFactor(t);
        const m_factor = 1.0 + msf_factor;
        const Hz = this.computeHtz(this.variables.get("z"));
        const expansion = 1.0 + Hz * t;
        const sc_correction = 1.0 - (this.variables.get("B") / this.variables.get("B_crit"));
        const f_env = this.computeFenv(t);
        const tr_factor = 1.0 + this.variables.get("f_TRZ");

        // Base gravity
        const g_base = (this.variables.get("G") * this.variables.get("M") * m_factor / (r * r)) *
                      expansion * sc_correction * (1.0 + f_env) * tr_factor;

        // Ug sum (includes base? Adjust: Ug sum without base)
        const ug_sum = this.computeUgSum(r) - g_base;  // Subtract to avoid double-count

        // Cosmological
        const lambda_term = this.variables.get("Lambda") * (this.variables.get("c") * this.variables.get("c")) / 3.0;

        // Ui
        const ui_term = this.computeUi(t);

        // Quantum
        const quantum_term = this.computeQuantumTerm(this.variables.get("t_Hubble"), r);

        // Fluid
        const fluid_term = this.computeFluidTerm(g_base);

        // DM
        const dm_term = this.computeDMTerm(r);

        // Total
        return g_base + ug_sum + lambda_term + ui_term + quantum_term + fluid_term + dm_term;
    }

    // Equation description
    getEquationText() {
        return "g_NGC1300(r, t) = (G * M(t) / r(t)^2) * (1 + H(t, z)) * (1 - B(t) / B_crit) * (1 + F_env(t)) + " +
               "(U_g1 + U_g2 + U_g3' + U_g4) + U_i + (Lambda * c^2 / 3) + " +
               "(hbar / sqrt(?x * ?p)) * ? (?_total * H * ?_total dV) * (2? / t_Hubble) + " +
               "?_fluid * V * g + (M_visible + M_DM) * (??/? + 3 G M / r^3)\n" +
               "Where: M(t) = M * (1 + M_SF(t)); M_SF(t) = SFR * t; r(t) = r0 + v_r t;\n" +
               "H(t, z) = H0 * sqrt(?m (1+z)^3 + ??); F_env(t) = F_bar + F_SF + F_wave;\n" +
               "F_bar = 0.1 G M / r^2; F_wave = ? v_arm^2; U_g1 = ?_dipole * B; U_g2 = B_super^2 / (2 ?0);\n" +
               "U_g3' = G M_bar / r_bar^2; U_g4 = k4 * E_react(t); U_i = ?_I * (?_SCm/?_UA) * ?_i * cos(? t_n) * (1 + F_RZ);\n" +
               "?_total = A exp(-r^2/(2?^2)) exp(i(m? - ? t)) + bar terms; Insights: Attractive (g_base, Ug1, Ug3') vs. Repulsive (U_g2, ?) advance UQFF.\n" +
               "Adaptations: Hubble ACS 2004 data; SFR=1 Msun/yr; M=1e11 Msun. Solutions: g ~2e36 m/s² at t=1 Gyr (DM/fluid dominant).";
    }

    // Print all current variables (for debugging)
    printVariables() {
        console.log('NGC 1300 Variables:');
        for (const [key, value] of this.variables.entries()) {
            console.log(`  ${key} = ${typeof value === 'number' ? value.toExponential(3) : value}`);
        }
    }

    // Dynamic parameter update method
    updateParameter(paramName, newValue) {
        if (this.hasOwnProperty(paramName)) {
            this[paramName] = newValue;
            if (this.updateCache) {
                this.updateCache();
            }
            return true;
        }
        return this.updateVariable(paramName, newValue) || false;
    }

    // Dynamic method expansion
    expand(methodName, methodFunction) {
        if (typeof methodFunction === 'function') {
            this[methodName] = methodFunction;
            return true;
        }
        return false;
    }

    // Install NGC 1300 UQFF module (for compatibility)
    install_uqff_module() {
        console.log("NGC1300EnhancedUQFFModule installed for NGC 1300 Barred Spiral Galaxy Evolution");
        console.log("Models gravitational dynamics with bar-driven gas funneling, spiral arm density waves, star formation, dust lanes, and dark matter");
        console.log("Supports dynamic variable management and comprehensive UQFF physics framework");
    }
}

module.exports = { NGC1300EnhancedUQFFModule };