// Source71 NGC1316UQFFModule
// JavaScript implementation of the NGC1316UQFFModule for NGC 1316 (Hubble Spies Cosmic Dust Bunnies) gravitational dynamics
// Models NGC 1316's gravitational dynamics, incorporating merger history, tidal forces, star cluster disruption, dust lanes, AGN jets/radio lobes, and dark matter.
// Supports dynamic variable management and F_env(t) with tidal and cluster terms.

class NGC1316UQFFModule {
    constructor() {
        this.variables = new Map();
        this.initializeConstants();
    }

    // Initialize universal constants and NGC 1316-specific parameters
    initializeConstants() {
        // Universal constants
        this.variables.set("G", 6.6743e-11);                    // m^3 kg^-1 s^-2
        this.variables.set("c", 3e8);                           // m/s
        this.variables.set("hbar", 1.0546e-34);                 // J s
        this.variables.set("Lambda", 1.1e-52);                  // m^-2
        this.variables.set("q", 1.602e-19);                     // C
        this.variables.set("pi", Math.PI);                      // pi
        this.variables.set("t_Hubble", 13.8e9 * 3.156e7);       // s
        this.variables.set("year_to_s", 3.156e7);               // s/yr
        this.variables.set("H0", 70.0);                         // km/s/Mpc
        this.variables.set("Mpc_to_m", 3.086e22);               // m/Mpc
        this.variables.set("Omega_m", 0.3);
        this.variables.set("Omega_Lambda", 0.7);

        const M_sun_val = 1.989e30;                             // kg
        const kpc_val = 3.086e19;                               // m

        // NGC 1316 parameters
        this.variables.set("M_visible", 3.5e11 * M_sun_val);    // kg
        this.variables.set("M_DM", 1.5e11 * M_sun_val);         // kg
        this.variables.set("M", this.variables.get("M_visible") + this.variables.get("M_DM"));  // Total initial
        this.variables.set("M0", this.variables.get("M"));
        this.variables.set("M_spiral", 1e10 * M_sun_val);       // kg (merger progenitor)
        this.variables.set("d_spiral", 50e3 * kpc_val);         // m
        this.variables.set("M_BH", 1e8 * M_sun_val);            // kg (AGN BH)
        this.variables.set("M_cluster", 1e6 * M_sun_val);       // kg
        this.variables.set("r", 46e3 * kpc_val);                // m
        this.variables.set("z", 0.005);                         // Redshift
        this.variables.set("tau_merge", 1e9 * this.variables.get("year_to_s"));  // s
        this.variables.set("t", 2e9 * this.variables.get("year_to_s"));  // Default t=2 Gyr s

        // Dynamics
        this.variables.set("rho_dust", 1e-21);                  // kg/m^3
        this.variables.set("V", 1e51);                          // m^3
        this.variables.set("B", 1e-4);                          // T (AGN)
        this.variables.set("B_crit", 1e11);                     // T
        this.variables.set("Delta_x", 1e-10);                   // m
        this.variables.set("Delta_p", this.variables.get("hbar") / this.variables.get("Delta_x"));
        this.variables.set("integral_psi", 1.0);
        this.variables.set("A", 1e-10);
        this.variables.set("k", 1e20);
        this.variables.set("omega", 1e-16);
        this.variables.set("x", 0.0);
        this.variables.set("v", 1e3);
        this.variables.set("sigma", 6.172e22);
        this.variables.set("Ug1", 0.0);
        this.variables.set("Ug2", 0.0);
        this.variables.set("Ug3", 0.0);
        this.variables.set("Ug4", 0.0);
        this.variables.set("Ui", 0.0);
        this.variables.set("mu_0", 1.2566370614359173e-6);
        this.variables.set("rho_vac_SCm", 7.09e-37);
        this.variables.set("rho_vac_UA", 7.09e-36);
        this.variables.set("lambda_I", 1.0);
        this.variables.set("omega_i", 1e-8);
        this.variables.set("t_n", 0.0);
        this.variables.set("F_RZ", 1e-2);
        this.variables.set("k_4", 1.0);
        this.variables.set("k_cluster", 1e-12);
        this.variables.set("omega_spin", 1e-3);
        this.variables.set("I_dipole", 1e20);
        this.variables.set("A_dipole", 1e15);
        this.variables.set("H_aether", 1e-5);
        this.variables.set("delta_rho_over_rho", 1e-5);
        this.variables.set("scale_macro", 1e-12);
        this.variables.set("f_TRZ", 1e-1);
        this.variables.set("f_sc", 1.0);
        this.variables.set("v_r", 1e3);
        this.variables.set("rho", 1e-21);
        this.variables.set("integral_psi", 1.0);                // Normalized

        // Wave/oscillatory for dust lanes
        this.variables.set("A", 1e-10);
        this.variables.set("k", 1e20);
        this.variables.set("omega", 1e-16);                     // rad/s for dust dynamics
        this.variables.set("x", 0.0);
        this.variables.set("v", 1e3);                           // m/s
        this.variables.set("sigma", 2e3 * kpc_val);             // m for Gaussian

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
        this.variables.set("k_cluster", 1e-12);                 // N/Msun, adjusted to m/s^2
        this.variables.set("omega_spin", 1e-3);                 // rad/s BH spin
        this.variables.set("I_dipole", 1e20);                   // A
        this.variables.set("A_dipole", 1e15);                   // m^2
        this.variables.set("H_aether", 1e-5);                   // A/m
        this.variables.set("delta_rho_over_rho", 1e-5);

        // Scales
        this.variables.set("scale_macro", 1e-12);
        this.variables.set("f_TRZ", 0.1);
        this.variables.set("f_sc", 1.0);
        this.variables.set("v_r", 1e3);                         // m/s radial velocity
        this.variables.set("rho", this.variables.get("rho_dust"));
    }

    // Dynamic variable operations
    updateVariable(name, value) {
        this.variables.set(name, value);

        // Update dependent variables
        if (name === "Delta_x") {
            this.variables.set("Delta_p", this.variables.get("hbar") / value);
        } else if (name === "M_visible" || name === "M_DM") {
            const newM = this.variables.get("M_visible") + this.variables.get("M_DM");
            this.variables.set("M", newM);
            this.variables.set("M0", newM);
        }
    }

    addToVariable(name, delta) {
        const current = this.variables.get(name) || 0;
        this.variables.set(name, current + delta);
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

    // M_merge(t)
    computeMmerge(t) {
        return 1e10 * 1.989e30 * Math.exp(-t / this.variables.get("tau_merge"));
    }

    // r(t)
    computeRt(t) {
        return this.variables.get("r") + this.variables.get("v_r") * t;
    }

    // F_env(t)
    computeFenv(t) {
        const F_tidal = (this.variables.get("G") * this.variables.get("M_spiral")) /
                       (this.variables.get("d_spiral") * this.variables.get("d_spiral"));
        const F_cluster = this.variables.get("k_cluster") * (this.variables.get("M_cluster") / 1.989e30);  // Normalize to m/s^2
        return F_tidal + F_cluster;
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

    // Ug3': external
    computeUg3prime(t) {
        return (this.variables.get("G") * this.variables.get("M_spiral")) /
               (this.variables.get("d_spiral") * this.variables.get("d_spiral"));
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

    // Psi integral (simplified dust lanes)
    computePsiIntegral(r, t) {
        const A = this.variables.get("A");
        const m = 2.0;
        const omega = this.variables.get("omega");
        const sigma = this.variables.get("sigma");

        // Simplified complex exponential for dust lane dynamics
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

    // Fluid with rho_dust
    computeFluidTerm(g_base) {
        return this.variables.get("rho_dust") * this.variables.get("V") * g_base;
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

    // Core computation: g_NGC1316(r, t)
    computeG(t, r) {
        this.variables.set("t", t);

        const m_merge = this.computeMmerge(t);
        const m_factor = 1.0 + m_merge / this.variables.get("M0");
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

    // Get equation text
    getEquationText() {
        return "g_NGC1316(r, t) = (G * M(t) / r(t)^2) * (1 + H(t, z)) * (1 - B(t) / B_crit) * (1 + F_env(t)) + " +
               "(U_g1 + U_g2 + U_g3' + U_g4) + U_i + (Lambda * c^2 / 3) + " +
               "(hbar / sqrt(?x * ?p)) * ? (?_total * H * ?_total dV) * (2? / t_Hubble) + " +
               "?_dust * V * g + (M_visible + M_DM) * (??/? + 3 G M / r^3)\\n" +
               "Where: M(t) = M * (1 + M_merge(t)); M_merge(t) = 1e10 Msun * exp(-t/?); r(t) = r0 + v_r t;\\n" +
               "H(t, z) = H0 * sqrt(?m (1+z)^3 + ??); F_env(t) = F_tidal + F_cluster;\\n" +
               "F_tidal = G M_spiral / d^2; F_cluster = k_cluster * M_cluster; U_g1 = ?_dipole * B; U_g2 = B_super^2 / (2 ?0);\\n" +
               "U_g3' = G M_spiral / d^2; U_g4 = k4 * E_react(t); U_i = ?_I * (?_SCm/?_UA) * ?_i * cos(? t_n) * (1 + F_RZ);\\n" +
               "?_total = A exp(-r^2/(2?^2)) exp(i(m? - ? t)) + BH terms; Insights: Attractive (g_base, Ug1, Ug3') vs. Repulsive (U_g1, U_g2, ?) advance UQFF.\\n" +
               "Adaptations: Hubble ACS 2003 data; M=5e11 Msun; rho_dust=1e-21 kg/m^3. Solutions: g ~2e37 m/s² at t=2 Gyr (DM/dust dominant).";
    }

    // Print variables (for debugging)
    printVariables() {
        console.log("NGC 1316 Variables:");
        for (const [key, value] of this.variables) {
            if (typeof value === 'number' && !isNaN(value)) {
                console.log(`${key} = ${value.toExponential()}`);
            } else {
                console.log(`${key} = ${value} (non-numeric or NaN)`);
            }
        }
    }

    // Install NGC 1316 UQFF module (for compatibility)
    install_uqff_module() {
        console.log("NGC1316UQFFModule installed for NGC 1316 (Hubble Spies Cosmic Dust Bunnies) gravitational dynamics");
        console.log("Models NGC 1316's merger history, tidal forces, star cluster disruption, dust lanes, AGN jets, and dark matter");
        console.log("Supports dynamic variable management and F_env(t) with tidal and cluster terms");
    }
}

module.exports = { NGC1316UQFFModule };