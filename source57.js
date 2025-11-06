class MultiCompressedUQFFModule {
    constructor(system = "MagnetarSGR1745") {
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
        this.variables.set("B", 1e-5);                          // T (default)
        this.variables.set("B_crit", 1e11);                     // T
        this.variables.set("rho_fluid", 1e-20);                 // kg/m^3 (default)
        this.variables.set("delta_rho_over_rho", 1e-5);
        this.variables.set("integral_psi_total", 1.0);          // Combined waves
        this.variables.set("Delta_x_Delta_p", 1e-68);           // J^2 s^2
        this.variables.set("M_DM", 0.0);                        // Default no DM
        this.variables.set("M_visible", 0.0);                   // Set per system
        this.variables.set("M_ext", 0.0);                       // For Ug3'
        this.variables.set("r_ext", 0.0);
        this.variables.set("f_sc", 10.0);

        // Set initial system
        this.setSystem(system);
    }

    // Set system: Load system-specific vars
    setSystem(system) {
        this.current_system = system;
        const M_sun = 1.989e30;
        this.variables.set("M", 0.0);
        this.variables.set("r", 0.0);
        this.variables.set("z", 0.0);
        this.variables.set("t_default", 0.0);
        this.variables.set("SFR", 0.0);
        this.variables.set("M0", 0.0);
        this.variables.set("M_visible", 0.0);
        this.variables.set("M_ext", 0.0);
        this.variables.set("r_ext", 0.0);
        this.variables.set("v_wind", 0.0);  // For F_env

        if (system === "MagnetarSGR1745") {
            this.variables.set("M", 2.8 * M_sun);               // kg
            this.variables.set("r", 1e4);                       // m
            this.variables.set("z", 0.026);
            this.variables.set("t_default", 1e3 * this.variables.get("year_to_s"));     // 1 kyr
            this.variables.set("SFR", 0.0);
            this.variables.set("M0", this.variables.get("M"));
            this.variables.set("M_visible", this.variables.get("M"));
            this.variables.set("M_ext", 4e6 * M_sun);           // Sgr A* M_BH
            this.variables.set("r_ext", 8e9);                   // m (distance)
            this.variables.set("v_wind", 1e5);                  // m/s
        } else if (system === "SagittariusA") {
            this.variables.set("M", 4e6 * M_sun);
            this.variables.set("r", 1e10);                      // m (event horizon scale)
            this.variables.set("z", 0.0);
            this.variables.set("t_default", 1e6 * this.variables.get("year_to_s"));     // 1 Myr
            this.variables.set("SFR", 0.0);
            this.variables.set("M0", this.variables.get("M"));
            this.variables.set("M_visible", this.variables.get("M"));
            this.variables.set("M_ext", 0.0);
            this.variables.set("r_ext", 0.0);
            this.variables.set("v_wind", 1e8);                  // m/s (relativistic)
        } else if (system === "TapestryStarbirth" || system === "Westerlund2") {
            this.variables.set("M", 1e4 * M_sun);
            this.variables.set("r", 1e18);                      // m (~10 pc)
            this.variables.set("z", 0.001);
            this.variables.set("t_default", 5e6 * this.variables.get("year_to_s"));     // 5 Myr
            this.variables.set("SFR", 0.1 * M_sun);
            this.variables.set("M0", this.variables.get("M"));
            this.variables.set("M_visible", this.variables.get("M"));
            this.variables.set("M_ext", 0.0);
            this.variables.set("r_ext", 0.0);
            this.variables.set("v_wind", 1e3);                  // m/s
        } else if (system === "PillarsCreation") {
            this.variables.set("M", 800 * M_sun);
            this.variables.set("r", 3e17);                      // m (~3 ly)
            this.variables.set("z", 0.0018);
            this.variables.set("t_default", 2e6 * this.variables.get("year_to_s"));     // 2 Myr
            this.variables.set("SFR", 0.1 * M_sun);
            this.variables.set("M0", this.variables.get("M"));
            this.variables.set("M_visible", this.variables.get("M"));
            this.variables.set("M_ext", 0.0);
            this.variables.set("r_ext", 0.0);
            this.variables.set("v_wind", 1e4);                  // m/s
        } else if (system === "RingsRelativity") {
            this.variables.set("M", 1e11 * M_sun);              // Galaxy mass
            this.variables.set("r", 1e21);                      // m (~100 kpc)
            this.variables.set("z", 0.5);
            this.variables.set("t_default", 1e10 * this.variables.get("year_to_s"));    // 10 Gyr
            this.variables.set("SFR", 0.0);
            this.variables.set("M0", this.variables.get("M"));
            this.variables.set("M_visible", this.variables.get("M"));
            this.variables.set("M_ext", 0.0);
            this.variables.set("r_ext", 0.0);
            this.variables.set("v_wind", 0.0);
        } else if (system === "UniverseGuide") {
            this.variables.set("M", 1 * M_sun);
            this.variables.set("r", 1.496e11);                  // AU
            this.variables.set("z", 0.0);
            this.variables.set("t_default", 4.35e17);           // Hubble time
            this.variables.set("SFR", 0.0);
            this.variables.set("M0", this.variables.get("M"));
            this.variables.set("M_visible", this.variables.get("M"));
            this.variables.set("M_ext", 0.0);
            this.variables.set("r_ext", 0.0);
            this.variables.set("v_wind", 0.0);
        }

        this.variables.set("rho_fluid", 1e-20);  // Default, override if needed
        this.variables.set("V", 1.0 / this.variables.get("rho_fluid"));
        this.variables.set("M_DM", 0.85 * this.variables.get("M"));  // Default fraction
        this.variables.set("M_visible", 0.15 * this.variables.get("M"));  // Adjust if needed
    }

    // Dynamic variable operations
    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
        } else {
            console.error(`Variable '${name}' not found. Adding with value ${value}`);
            this.variables.set(name, value);
        }

        if (name === "M") {
            this.variables.set("M0", value);
            this.variables.set("M_DM", 0.85 * value);
            this.variables.set("M_visible", 0.15 * value);
        } else if (name === "rho_fluid") {
            this.variables.set("V", 1.0 / value);
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

    // Compute H(t, z) in s^-1
    computeHtz(z) {
        const Hz_kms = this.variables.get("H0") * Math.sqrt(
            this.variables.get("Omega_m") * Math.pow(1.0 + z, 3) +
            this.variables.get("Omega_Lambda")
        );
        return (Hz_kms * 1e3) / this.variables.get("Mpc_to_m");
    }

    // Compute F_env(t): System-specific environmental term
    computeF_env(t) {
        let f_env = 1.0;  // Base
        if (this.current_system === "MagnetarSGR1745") {
            const M_mag = 1e40;  // J (est magnetic energy)
            const D_t = Math.exp(-t / (1e3 * this.variables.get("year_to_s")));  // Decay
            f_env += M_mag / (this.variables.get("M") * this.variables.get("c") * this.variables.get("c")) + D_t;
        } else if (this.current_system === "SagittariusA") {
            const omega_dot = 1e-3;  // rad/s (est spin)
            f_env += (Math.pow(this.variables.get("G") * this.variables.get("M"), 2) /
                     (Math.pow(this.variables.get("c"), 4) * this.variables.get("r"))) * Math.pow(omega_dot, 2);
        } else if (this.current_system === "TapestryStarbirth" || this.current_system === "Westerlund2") {
            f_env += this.variables.get("rho_fluid") * Math.pow(this.variables.get("v_wind"), 2);
        } else if (this.current_system === "PillarsCreation") {
            const E_t = 1.0 - Math.exp(-t / (2e6 * this.variables.get("year_to_s")));  // Erosion
            f_env += this.variables.get("rho_fluid") * Math.pow(this.variables.get("v_wind"), 2) * E_t;
        } else if (this.current_system === "RingsRelativity") {
            const L_t = 1.0 + 0.1 * Math.sin(2 * this.variables.get("pi") * t / this.variables.get("t_Hubble"));  // Lensing variation
            f_env += L_t;
        } else if (this.current_system === "UniverseGuide") {
            f_env += 0.0;  // Minimal
        }
        return f_env;
    }

    // Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral_psi_total * (2 pi / t_Hubble)
    computeQuantumTerm(t_Hubble_val) {
        const sqrt_unc = Math.sqrt(this.variables.get("Delta_x_Delta_p"));
        const integral_val = this.variables.get("integral_psi_total");
        return (this.variables.get("hbar") / sqrt_unc) * integral_val * (2 * this.variables.get("pi") / t_Hubble_val);
    }

    // Fluid term: rho_fluid * V * g_base
    computeFluidTerm(g_base) {
        return this.variables.get("rho_fluid") * this.variables.get("V") * g_base;
    }

    // Ug sum: Ug1 = G M / r^2, Ug2=0, Ug3' = G M_ext / r_ext^2, Ug4 = Ug1 * f_sc
    computeUgSum(r) {
        const G = this.variables.get("G");
        const M = this.variables.get("M");
        const Ug1 = (G * M) / (r * r);
        this.variables.set("Ug1", Ug1);
        this.variables.set("Ug2", 0.0);
        const Ug3_prime = (this.variables.get("M_ext") > 0) ?
            (G * this.variables.get("M_ext")) / (this.variables.get("r_ext") * this.variables.get("r_ext")) : 0.0;
        this.variables.set("Ug3", Ug3_prime);
        const Ug4 = Ug1 * this.variables.get("f_sc");
        this.variables.set("Ug4", Ug4);
        return this.variables.get("Ug1") + this.variables.get("Ug2") + this.variables.get("Ug3") + this.variables.get("Ug4");
    }

    // Star formation factor: (SFR * t_yr) / M0
    computeMsfFactor(t) {
        if (this.variables.get("SFR") === 0.0) return 0.0;
        const t_yr = t / this.variables.get("year_to_s");
        return (this.variables.get("SFR") * t_yr) / this.variables.get("M0");
    }

    // DM pert term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3)
    computeDMPertTerm(r) {
        const pert = this.variables.get("delta_rho_over_rho") + 3 * this.variables.get("G") * this.variables.get("M") / Math.pow(r, 3);
        return (this.variables.get("M_visible") + this.variables.get("M_DM")) * pert;
    }

    // Full compressed computation
    computeG(t) {
        this.variables.set("t", t);
        const z = this.variables.get("z");
        const Hz = this.computeHtz(z);
        const expansion = 1.0 + Hz * t;
        const sc_correction = 1.0 - (this.variables.get("B") / this.variables.get("B_crit"));
        const f_env = this.computeF_env(t);
        const msf_factor = this.computeMsfFactor(t);
        const m_factor = 1.0 + msf_factor;
        const r = this.variables.get("r");

        // Base gravity with expansion, SC, F_env, M(t)
        const g_base = (this.variables.get("G") * this.variables.get("M") * m_factor / (r * r)) * expansion * sc_correction * f_env;

        // Ug sum
        const ug_sum = this.computeUgSum(r);

        // Cosmological
        const lambda_term = this.variables.get("Lambda") * (this.variables.get("c") * this.variables.get("c")) / 3.0;

        // Quantum (psi_total)
        const quantum_term = this.computeQuantumTerm(this.variables.get("t_Hubble"));

        // Fluid
        const fluid_term = this.computeFluidTerm(g_base);

        // DM pert
        const dm_pert_term = this.computeDMPertTerm(r);

        // Total: Sum all
        return g_base + ug_sum + lambda_term + quantum_term + fluid_term + dm_pert_term;
    }

    // Get equation text (descriptive)
    getEquationText() {
        return "g_UQFF(r, t) = (G * M(t) / r^2) * (1 + H(t, z)) * (1 - B(t) / B_crit) * (1 + F_env(t)) + (Ug1 + Ug2 + Ug3' + Ug4) + (Lambda * c^2 / 3) + " +
               "(hbar / sqrt(Delta_x * Delta_p)) * ∫(ψ_total H ψ_total dV) * (2π / t_Hubble) + ρ_fluid * V * g + (M_visible + M_DM) * (δρ/ρ + 3 G M / r^3)\n" +
               "Where H(t, z) = H_0 * sqrt(Ω_m (1+z)^3 + Ω_Λ); M(t) = M * (1 + M_sf(t)); M_sf(t) = (SFR * t_yr) / M0;\n" +
               "F_env(t) = system-specific (e.g., ρ v_wind^2 for Starbirth, E(t) for Pillars); Ug3' = G M_ext / r_ext^2;\n" +
               "ψ_total = combined waves (magnetic + standing + quantum).\n" +
               "Special Terms:\n" +
               "- Compression: Unified H(t,z), F_env(t) modular, Ug3' generalized, ψ_total consolidated.\n" +
               "- Adaptations: Magnetar (M_BH, decay); SgrA* (GW spin); Starbirth/Westerlund2 (winds); Pillars (erosion); Rings (lensing); UniverseGuide (solar).\n" +
               "Solutions: Varies by system/t; e.g., Magnetar t=1kyr ~1e12 m/s² (B_crit dominant).\n" +
               "From UQFF Cycle 2: Streamlines 7 systems, reduces redundancy.";
    }

    // Print variables
    printVariables() {
        console.log(`Current Variables for ${this.current_system}:`);
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }

    // Install method for integration
    install_uqff_module() {
        console.log("MultiCompressedUQFFModule installed for compressed UQFF calculations");
        console.log("Supports 7 astrophysical systems: MagnetarSGR1745, SagittariusA, TapestryStarbirth, Westerlund2, PillarsCreation, RingsRelativity, UniverseGuide");
        console.log("Includes modular F_env(t) environmental terms and unified H(t,z) cosmology");
    }
}

module.exports = { MultiCompressedUQFFModule };