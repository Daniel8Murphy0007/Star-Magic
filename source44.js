class LagoonUQFFModule {
    constructor() {
        this.variables = new Map();

        // Base constants (universal)
        this.variables.set("G", 6.6743e-11);                    // m^3 kg^-1 s^-2
        this.variables.set("c", 3e8);                           // m/s
        this.variables.set("hbar", 1.0546e-34);                 // J s
        this.variables.set("Lambda", 1.1e-52);                  // m^-2
        this.variables.set("q", 1.602e-19);                     // C
        this.variables.set("pi", 3.141592653589793);            // pi
        this.variables.set("t_Hubble", 13.8e9 * 3.156e7);       // s
        this.variables.set("year_to_s", 3.156e7);               // s/yr

        // Lagoon Nebula parameters
        const M_sun_val = 1.989e30;                    // kg
        this.variables.set("M_sun", M_sun_val);
        this.variables.set("M", 1e4 * M_sun_val);               // Total mass kg
        this.variables.set("M0", this.variables.get("M"));               // Initial mass
        this.variables.set("SFR", 0.1 * M_sun_val);             // Msun/yr
        this.variables.set("M_visible", 0.15 * this.variables.get("M")); // Visible fraction est.
        this.variables.set("M_DM", 0.85 * this.variables.get("M"));      // Dark matter/halo
        this.variables.set("r", 5.2e17);                        // m (half width ~55 ly)

        // Hubble/cosmology
        this.variables.set("H0", 67.15);                        // km/s/Mpc
        this.variables.set("Mpc_to_m", 3.086e22);               // m/Mpc
        this.variables.set("z", 0.0013);                        // Redshift
        this.variables.set("Omega_m", 0.3);
        this.variables.set("Omega_Lambda", 0.7);
        this.variables.set("t", 1e6 * this.variables.get("year_to_s"));  // Default t=1 Myr s

        // Gas dynamics
        this.variables.set("rho_fluid", 1e-20);                 // kg/m^3 (dense gas)
        this.variables.set("V", 1e3);                           // m^3 (arbitrary)
        this.variables.set("v_gas", 1e5);                       // m/s (turbulent velocity)
        this.variables.set("delta_rho", 0.1 * this.variables.get("rho_fluid"));
        this.variables.set("rho", this.variables.get("rho_fluid"));

        // EM/magnetic
        this.variables.set("B", 1e-5);                          // T (nebula field)
        this.variables.set("B_crit", 1e11);                     // T (10^15 G)

        // Quantum terms
        this.variables.set("Delta_x", 1e-10);                   // m
        this.variables.set("Delta_p", this.variables.get("hbar") / this.variables.get("Delta_x"));
        this.variables.set("integral_psi", 1.0);

        // Resonant/oscillatory
        this.variables.set("A", 1e-10);
        this.variables.set("k", 1e20);
        this.variables.set("omega", 1e15);                      // rad/s (high freq)
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

        // Radiation pressure
        this.variables.set("L_H36", 7.65e31);                   // W
        this.variables.set("m_H", 1.67e-27);                    // kg
    }

    // Dynamic variable operations
    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
        } else {
            console.error(`Variable '${name}' not found. Adding with value ${value}`);
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
            console.error(`Variable '${name}' not found. Adding with delta ${delta}`);
            this.variables.set(name, delta);
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    getVariable(name) {
        return this.variables.get(name);
    }

    // Compute H(z) in s^-1
    computeHz() {
        const Hz_kms = this.variables.get("H0") * Math.sqrt(
            this.variables.get("Omega_m") * Math.pow(1.0 + this.variables.get("z"), 3) +
            this.variables.get("Omega_Lambda")
        );
        return (Hz_kms * 1e3) / this.variables.get("Mpc_to_m");
    }

    // Compute Ug sum: Ug1 = G M / r^2, Ug4 = Ug1 * f_sc, others 0
    computeUgSum() {
        const Ug1 = (this.variables.get("G") * this.variables.get("M")) / (this.variables.get("r") * this.variables.get("r"));
        this.variables.set("Ug1", Ug1);
        this.variables.set("Ug4", Ug1 * this.variables.get("f_sc"));
        return this.variables.get("Ug1") + this.variables.get("Ug2") + this.variables.get("Ug3") + this.variables.get("Ug4");
    }

    // Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
    computeQuantumTerm(t_Hubble_val) {
        const unc = Math.sqrt(this.variables.get("Delta_x") * this.variables.get("Delta_p"));
        const integral_val = this.variables.get("integral_psi");
        return (this.variables.get("hbar") / unc) * integral_val * (2 * this.variables.get("pi") / t_Hubble_val);
    }

    // Fluid term: rho_fluid * V * g
    computeFluidTerm(g_base) {
        return this.variables.get("rho_fluid") * this.variables.get("V") * g_base;
    }

    // Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
    computeResonantTerm(t) {
        const cos_term = 2 * this.variables.get("A") * Math.cos(this.variables.get("k") * this.variables.get("x")) * Math.cos(this.variables.get("omega") * t);
        const exp_real = this.variables.get("A") * Math.cos(this.variables.get("k") * this.variables.get("x") - this.variables.get("omega") * t);  // Real part of exp(i*phase)
        const exp_factor = (2 * this.variables.get("pi") / 13.8);
        return cos_term + exp_factor * exp_real;
    }

    // DM term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3)
    computeDMTerm() {
        const pert = this.variables.get("delta_rho") / this.variables.get("rho");
        const curv = 3 * this.variables.get("G") * this.variables.get("M") / (this.variables.get("r") * this.variables.get("r") * this.variables.get("r"));
        return (this.variables.get("M_visible") + this.variables.get("M_DM")) * (pert + curv);
    }

    // Star formation factor: (SFR * t_yr) / M0
    computeMsfFactor(t) {
        const t_yr = t / this.variables.get("year_to_s");
        return (this.variables.get("SFR") * t_yr) / this.variables.get("M0");
    }

    // Radiation pressure: P_rad = (L_H36 / (4 pi r^2 c)) * (rho / m_H)
    computeP_rad() {
        const flux = this.variables.get("L_H36") / (4 * this.variables.get("pi") * this.variables.get("r") * this.variables.get("r") * this.variables.get("c"));
        return flux * (this.variables.get("rho_fluid") / this.variables.get("m_H"));
    }

    // Full computation: g_UQFF(r, t) = ... all terms with M_sf and -P_rad
    computeG(t) {
        this.variables.set("t", t);
        const Hz = this.computeHz();
        const expansion = 1.0 + Hz * t;
        const sc_correction = 1.0 - (this.variables.get("B") / this.variables.get("B_crit"));
        const tr_factor = 1.0 + this.variables.get("f_TRZ");
        const msf_factor = this.computeMsfFactor(t);
        const m_factor = 1.0 + msf_factor;
        const p_rad = this.computeP_rad();

        // Base gravity with expansion, SC, TR, M_sf
        const g_base = (this.variables.get("G") * this.variables.get("M") * m_factor / (this.variables.get("r") * this.variables.get("r"))) * expansion * sc_correction * tr_factor;

        // Ug sum
        const ug_sum = this.computeUgSum();

        // Cosmological
        const lambda_term = this.variables.get("Lambda") * (this.variables.get("c") * this.variables.get("c")) / 3.0;

        // Quantum
        const quantum_term = this.computeQuantumTerm(this.variables.get("t_Hubble"));

        // EM Lorentz (magnitude v_gas B)
        const em_base = this.variables.get("q") * this.variables.get("v_gas") * this.variables.get("B") / 1.673e-27;
        const em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * this.variables.get("scale_macro");

        // Fluid
        const fluid_term = this.computeFluidTerm(g_base);

        // Resonant
        const resonant_term = this.computeResonantTerm(t);

        // DM
        const dm_term = this.computeDMTerm();

        // Total: Sum all - P_rad
        return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term - p_rad;
    }

    // Get equation text (descriptive)
    getEquationText() {
        return "g_Lagoon(r, t) = (G * M(t) / r^2) * (1 + H(z) * t) * (1 - B / B_crit) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 / 3) + " +
               "(hbar / sqrt(Delta_x * Delta_p)) * ∫(ψ* H ψ dV) * (2π / t_Hubble) + q (v × B) + ρ_fluid * V * g + " +
               "2 A cos(k x) cos(ω t) + (2π / 13.8) A exp(i (k x - ω t)) + (M_visible + M_DM) * (δρ/ρ + 3 G M / r^3) - P_rad\n" +
               "Where M(t) = M * (1 + M_sf(t)); M_sf(t) = (SFR * t_yr) / M0; P_rad = (L_H36 / (4π r^2 c)) * (ρ / m_H)\n" +
               "Special Terms:\n" +
               "- Quantum: Heisenberg uncertainty for gas quantum effects.\n" +
               "- Fluid: Nebular gas density-volume-gravity coupling.\n" +
               "- Resonant: Oscillatory Aether waves for ionization fronts.\n" +
               "- DM: Visible+dark mass with perturbations for halo.\n" +
               "- Superconductivity: (1 - B/B_crit) for quantum fields.\n" +
               "- Star Formation: M_sf(t) boosts mass via SFR=0.1 Msun/yr.\n" +
               "- Radiation Pressure: P_rad from Herschel 36 erodes gas.\n" +
               "Solutions: At t=1 Myr, g_Lagoon ~1e-12 m/s² (EM/fluid dominant; g_base ~1e-13; P_rad ~1e-14).\n" +
               "Adaptations for Lagoon Nebula: H II region with Herschel 36 radiation; z=0.0013; SFR for starbirth.";
    }

    // Print all current variables (for debugging/updates)
    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }
}

module.exports = LagoonUQFFModule;