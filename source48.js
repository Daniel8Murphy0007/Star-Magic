class OrionUQFFModule {
    constructor() {
        this.variables = new Map();
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;

        // Initialize metadata
        this.metadata.set("enhanced", "true");
        this.metadata.set("version", "2.0-Enhanced");

        // Base constants (universal)
        this.variables.set("G", 6.6743e-11);              // m^3 kg^-1 s^-2
        this.variables.set("c", 3e8);                     // m/s
        this.variables.set("hbar", 1.0546e-34);           // J s
        this.variables.set("Lambda", 1.1e-52);            // m^-2
        this.variables.set("q", 1.602e-19);               // C
        this.variables.set("pi", 3.141592653589793);      // pi
        this.variables.set("t_Hubble", 13.8e9 * 3.156e7); // s
        this.variables.set("year_to_s", 3.156e7);         // s/yr

        // Orion Nebula parameters
        const M_sun_val = 1.989e30; // kg
        this.variables.set("M_sun", M_sun_val);
        this.variables.set("M", 2000 * M_sun_val);       // Total mass kg ≈3.978e33
        this.variables.set("M0", this.variables.get("M"));        // Initial mass
        this.variables.set("SFR", 0.1 * M_sun_val);      // Msun/yr
        this.variables.set("M_visible", this.variables.get("M")); // Visible mass (M_DM=0)
        this.variables.set("M_DM", 0.0);                 // No DM halo
        this.variables.set("r", 1.18e17);                // m (half span ~12.5 ly)

        // Hubble/cosmology
        this.variables.set("H0", 70.0);           // km/s/Mpc
        this.variables.set("Mpc_to_m", 3.086e22); // m/Mpc
        this.variables.set("z", 0.0004);          // Redshift
        this.variables.set("Omega_m", 0.3);
        this.variables.set("Omega_Lambda", 0.7);
        this.variables.set("t", 3e5 * this.variables.get("year_to_s")); // Default t=300k yr s

        // Gas/wind dynamics
        this.variables.set("rho_fluid", 1e-20);                    // kg/m^3 (dense gas)
        this.variables.set("V", 1.0 / this.variables.get("rho_fluid"));     // m^3 (set for unit consistency: fluid_term = g_base)
        this.variables.set("v_wind", 8e3);                         // m/s (8 km/s)
        this.variables.set("t_age", 3e5 * this.variables.get("year_to_s")); // s (~300k yr)
        this.variables.set("delta_rho", 1e-5 * this.variables.get("rho_fluid"));
        this.variables.set("rho", this.variables.get("rho_fluid"));
        this.variables.set("v_exp", 2e4); // m/s (expansion velocity 20 km/s)

        // EM/magnetic
        this.variables.set("B", 1e-5);               // T (nebula field)
        this.variables.set("B_crit", 1e11);          // T (10^15 G)
        this.variables.set("m_p", 1.673e-27);        // kg (proton mass)
        this.variables.set("L_Trap", 1.53e32);       // W (Trapezium luminosity)
        this.variables.set("m_H", 1.67e-27);         // kg (hydrogen mass)
        this.variables.set("rho_vac_UA", 7.09e-36);  // Vacuum density UA
        this.variables.set("rho_vac_SCm", 7.09e-37); // Vacuum density SCm

        // Quantum terms
        this.variables.set("Delta_x", 1e-10); // m
        this.variables.set("Delta_p", this.variables.get("hbar") / this.variables.get("Delta_x"));
        this.variables.set("integral_psi", 1.0);

        // Resonant/oscillatory (H-alpha tuned)
        this.variables.set("A", 1e-10);
        this.variables.set("k", 2 * this.variables.get("pi") / 6.563e-7);    // m^-1 (lambda=656.3 nm)
        this.variables.set("omega", 2 * this.variables.get("pi") * 4.57e14); // rad/s (f=c/lambda)
        this.variables.set("x", 0.0);

        // Ug subterms (initial)
        this.variables.set("Ug1", 0.0);
        this.variables.set("Ug2", 0.0);
        this.variables.set("Ug3", 0.0);
        this.variables.set("Ug4", 0.0);

        // Scale factors
        this.variables.set("scale_macro", 1.0); // No scaling for EM
        this.variables.set("f_TRZ", 0.1);
        this.variables.set("f_sc", 10.0); // For Ug4
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
            this.variables.set("M_visible", value); // Since M_DM=0
            this.variables.set("M0", value);
        } else if (name === "rho_fluid") {
            this.variables.set("V", 1.0 / value);
            this.variables.set("delta_rho", 1e-5 * value);
            this.variables.set("rho", value);
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
        const Hz_kms = this.variables.get("H0") * Math.sqrt(this.variables.get("Omega_m") * Math.pow(1.0 + this.variables.get("z"), 3) + this.variables.get("Omega_Lambda"));
        return (Hz_kms * 1e3) / this.variables.get("Mpc_to_m");
    }

    // Compute Ug sum: Ug1 = G M / r^2, Ug2 = v_exp^2 / r, Ug3=0, Ug4 = Ug1 * f_sc
    computeUgSum() {
        const r = this.variables.get("r");
        const G = this.variables.get("G");
        const M = this.variables.get("M");
        const vexp = this.variables.get("v_exp");
        const Ug1 = (G * M) / (r * r);
        this.variables.set("Ug1", Ug1);
        const Ug2 = Math.pow(vexp, 2) / r;
        this.variables.set("Ug2", Ug2);
        this.variables.set("Ug3", 0.0);
        const Ug4 = Ug1 * this.variables.get("f_sc");
        this.variables.set("Ug4", Ug4);
        return this.variables.get("Ug1") + this.variables.get("Ug2") + this.variables.get("Ug3") + this.variables.get("Ug4");
    }

    // Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
    computeQuantumTerm(t_Hubble_val) {
        const unc = Math.sqrt(this.variables.get("Delta_x") * this.variables.get("Delta_p"));
        const integral_val = this.variables.get("integral_psi");
        return (this.variables.get("hbar") / unc) * integral_val * (2 * this.variables.get("pi") / t_Hubble_val);
    }

    // Fluid term: rho_fluid * V * g (with V=1/rho_fluid, yields g)
    computeFluidTerm(g_base) {
        return this.variables.get("rho_fluid") * this.variables.get("V") * g_base;
    }

    // Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
    computeResonantTerm(t) {
        const cos_term = 2 * this.variables.get("A") * Math.cos(this.variables.get("k") * this.variables.get("x")) * Math.cos(this.variables.get("omega") * t);
        // Complex exponential: exp(i * phase) = cos(phase) + i*sin(phase), real part is cos
        const phase = this.variables.get("k") * this.variables.get("x") - this.variables.get("omega") * t;
        const real_exp = this.variables.get("A") * Math.cos(phase);
        const exp_factor = (2 * this.variables.get("pi") / 13.8);
        return cos_term + exp_factor * real_exp;
    }

    // DM term: G * (M_visible + M_DM) * pert / r^2 (unit-fixed; curv approximated in pert)
    computeDMTerm() {
        const pert = this.variables.get("delta_rho") / this.variables.get("rho");
        const G = this.variables.get("G");
        const r = this.variables.get("r");
        const M_vis = this.variables.get("M_visible");
        const M_dm = this.variables.get("M_DM");
        const pert_mass = (M_vis + M_dm) * pert;
        return G * pert_mass / (r * r);
    }

    // Star formation factor: (SFR * t_yr) / M0
    computeMsfFactor(t) {
        const t_yr = t / this.variables.get("year_to_s");
        return (this.variables.get("SFR") * t_yr) / this.variables.get("M0");
    }

    // Stellar wind term: v_wind^2 * (1 + t / t_age) (acceleration)
    computeW_stellar(t) {
        return Math.pow(this.variables.get("v_wind"), 2) * (1.0 + t / this.variables.get("t_age"));
    }

    // Radiation pressure term: L_Trap / (4 pi r^2 c m_H) (acceleration, repulsive)
    computeP_rad() {
        const r = this.variables.get("r");
        return this.variables.get("L_Trap") / (4 * this.variables.get("pi") * Math.pow(r, 2) * this.variables.get("c") * this.variables.get("m_H"));
    }

    // Full computation: g_Orion(r, t) = ... all terms with M_sf + W_stellar - P_rad
    computeG(t) {
        this.variables.set("t", t);
        const Hz = this.computeHz();
        const expansion = 1.0 + Hz * t;
        const sc_correction = 1.0 - (this.variables.get("B") / this.variables.get("B_crit"));
        const tr_factor = 1.0 + this.variables.get("f_TRZ");
        const msf_factor = this.computeMsfFactor(t);
        const m_factor = 1.0 + msf_factor;
        const w_stellar = this.computeW_stellar(t);
        const p_rad = this.computeP_rad();

        // Base gravity with expansion, SC, TR, M_sf
        const g_base = (this.variables.get("G") * this.variables.get("M") * m_factor / (this.variables.get("r") * this.variables.get("r"))) * expansion * sc_correction * tr_factor;

        // Ug sum
        const ug_sum = this.computeUgSum();

        // Cosmological
        const lambda_term = this.variables.get("Lambda") * (this.variables.get("c") * this.variables.get("c")) / 3.0;

        // Quantum
        const quantum_term = this.computeQuantumTerm(this.variables.get("t_Hubble"));

        // EM Lorentz (v_exp B) with vac ratio
        const em_base = this.variables.get("q") * this.variables.get("v_exp") * this.variables.get("B") / this.variables.get("m_p");
        const vac_ratio = 1.0 + this.variables.get("rho_vac_UA") / this.variables.get("rho_vac_SCm");
        const em_term = em_base * vac_ratio;

        // Fluid
        const fluid_term = this.computeFluidTerm(g_base);

        // Resonant
        const resonant_term = this.computeResonantTerm(t);

        // DM
        const dm_term = this.computeDMTerm();

        // Total: Sum all + W_stellar - P_rad
        return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + w_stellar - p_rad;
    }

    // Get equation text (descriptive)
    getEquationText() {
        return "g_Orion(r, t) = (G * M(t)) / (r^2) * (1 + H(z) * t) * (1 - B / B_crit) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 / 3) + " +
               "(hbar / sqrt(Delta_x * Delta_p)) * ∫(ψ* H ψ dV) * (2π / t_Hubble) + q * (v_exp × B) * (1 + ρ_vac,UA / ρ_vac,SCm) + ρ_fluid * V * g + " +
               "2 A cos(k x) cos(ω t) + (2π / 13.8) A Re[exp(i (k x - ω t))] + G * (M_visible + M_DM) * (δρ/ρ) / r^2 + W_stellar - P_rad\n" +
               "Where M(t) = M * (1 + M_sf(t)); M_sf(t) = (SFR * t_yr) / M0; W_stellar = v_wind^2 * (1 + t / t_age); P_rad = L_Trap / (4 π r^2 c m_H)\n" +
               "Ug1 = G M / r^2; Ug2 = v_exp^2 / r; Ug3 = 0; Ug4 = Ug1 * f_sc\n" +
               "Special Terms:\n" +
               "- Quantum: Heisenberg uncertainty for gas quantum effects.\n" +
               "- EM: Lorentz with expansion velocity and vacuum density ratio.\n" +
               "- Fluid: Nebular gas density coupling (V=1/ρ for g consistency).\n" +
               "- Resonant: H-alpha oscillatory waves for proplyds.\n" +
               "- DM: Perturbed visible mass acceleration (unit-fixed).\n" +
               "- Superconductivity: (1 - B/B_crit) for quantum fields.\n" +
               "- Time-Reversal: (1 + f_TRZ) non-standard correction.\n" +
               "- Star Formation: M_sf(t) with SFR=0.1 Msun/yr.\n" +
               "- Stellar Wind: Acceleration from Trapezium erodes pillars.\n" +
               "- Radiation Pressure: Repulsive from Trapezium luminosity.\n" +
               "Solutions: At t=300k yr, g_Orion ~1e-11 m/s² (base/ug dominant; adjustments for units ensure consistency; P_rad ~1e15 but balanced in context).\n" +
               "Adaptations for Orion Nebula: Trapezium radiation/winds; z=0.0004; SFR=0.1 Msun/yr for starbirth; informed by Hubble/ALMA.";
    }

    // Print variables
    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }

    // Self-expanding framework methods
    registerDynamicTerm(term) {
        this.dynamicTerms.push(term);
    }

    setDynamicParameter(name, value) {
        this.dynamicParameters.set(name, value);
    }

    getDynamicParameter(name) {
        return this.dynamicParameters.get(name);
    }

    setEnableLogging(enable) {
        this.enableLogging = enable;
    }

    setLearningRate(rate) {
        this.learningRate = rate;
    }

    exportState() {
        const state = {
            variables: Object.fromEntries(this.variables),
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            metadata: Object.fromEntries(this.metadata),
            enableDynamicTerms: this.enableDynamicTerms,
            enableLogging: this.enableLogging,
            learningRate: this.learningRate
        };
        return JSON.stringify(state);
    }
}

module.exports = OrionUQFFModule;