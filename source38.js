// Source38.js - Compressed Resonance UQFF Module (JavaScript Implementation)
// Based on CompressedResonanceUQFFModule.cpp
// Implements UQFF compressed and resonance terms with dynamic variable management

class CompressedResonanceUQFFModule {
    constructor() {
        // Dynamic variables using Map (equivalent to std::map)
        this.variables = new Map();

        // Initialize all variables with UQFF-specific values for compressed/resonance
        // Base constants (UQFF universal)
        this.variables.set("c", 3e8);                           // m/s
        this.variables.set("pi", 3.141592653589793);            // pi
        this.variables.set("E_vac", 7.09e-36);                  // J/m^3 (plasmotic vacuum)
        this.variables.set("hbar", 1.0546e-34);                 // J s
        this.variables.set("f_TRZ", 0.1);                       // Time-reversal

        // Compressed parameters (streamlined DPM, THz, vac_diff, super)
        this.variables.set("f_DPM", 1e12);                      // Hz
        this.variables.set("f_THz", 1e12);                      // Hz
        this.variables.set("f_vac_diff", 0.143);                // Hz
        this.variables.set("f_super", 1.411e16);                // Hz
        this.variables.set("I", 1e21);                          // A
        this.variables.set("A_vort", 3.142e8);                  // m^2
        this.variables.set("omega_1", 1e-3);                    // rad/s
        this.variables.set("omega_2", -1e-3);                   // rad/s
        this.variables.set("v_exp", 1e3);                       // m/s
        this.variables.set("E_0", 6.381e-36);                   // J/m^3
        this.variables.set("V_sys", 4.189e12);                  // m^3

        // Resonance parameters (aether, U_g4i, osc, quantum, fluid, exp)
        this.variables.set("f_aether", 1e4);                    // Hz
        this.variables.set("f_react", 1e10);                    // Hz (U_g4i)
        this.variables.set("f_quantum", 1.445e-17);             // Hz
        this.variables.set("f_fluid", 1.269e-14);               // Hz
        this.variables.set("f_exp", 1.373e-8);                  // Hz
        this.variables.set("f_osc", 4.57e14);                   // Hz
        this.variables.set("k", 1e20);                          // m^-1
        this.variables.set("omega_osc", 1e15);                  // rad/s
        this.variables.set("x", 0.0);                           // m
        this.variables.set("A", 1e-10);                         // Amplitude
        this.variables.set("rho_fluid", 1e-21);                 // kg/m^3
        this.variables.set("V", 1e3);                           // m^3
        this.variables.set("delta_rho", 0.1 * this.variables.get("rho_fluid"));
        this.variables.set("rho", this.variables.get("rho_fluid"));

        // Superconductive integrated
        this.variables.set("B_crit", 1e11);                     // T
        this.variables.set("f_sc", 1.0);                        // Factor

        // Quantum
        this.variables.set("Delta_x", 1e-10);                   // m
        this.variables.set("Delta_p", this.variables.get("hbar") / this.variables.get("Delta_x"));
        this.variables.set("integral_psi", 1.0);
    }

    // Dynamic variable operations
    updateVariable(name, value) {
        this.variables.set(name, value);
        // Recompute dependents
        if (name === "Delta_x") {
            this.variables.set("Delta_p", this.variables.get("hbar") / value);
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

    // Compute Compressed Term: Sum streamlined DPM + THz + vac_diff + super
    computeCompressedTerm() {
        const F_DPM = this.variables.get("I") * this.variables.get("A_vort") *
                     (this.variables.get("omega_1") - this.variables.get("omega_2"));
        const a_DPM = (F_DPM * this.variables.get("f_DPM") * this.variables.get("E_vac")) /
                     (this.variables.get("c") * this.variables.get("V_sys"));
        const a_THz = (this.variables.get("f_THz") * this.variables.get("E_vac") *
                      this.variables.get("v_exp") * a_DPM) /
                     ((this.variables.get("E_vac") / 10) * this.variables.get("c"));
        const a_vac_diff = (this.variables.get("E_0") * this.variables.get("f_vac_diff") *
                           this.variables.get("V_sys") * a_DPM) / this.variables.get("hbar");
        const a_super = (this.variables.get("hbar") * this.variables.get("f_super") *
                        this.variables.get("f_DPM") * a_DPM) /
                       (this.variables.get("E_vac") * this.variables.get("c"));
        return a_DPM + a_THz + a_vac_diff + a_super;
    }

    // Compute Resonance Term: Sum aether + U_g4i + osc + quantum + fluid + exp
    computeResonanceTerm(t) {
        const F_DPM = this.variables.get("I") * this.variables.get("A_vort") *
                     (this.variables.get("omega_1") - this.variables.get("omega_2"));
        const a_DPM = (F_DPM * this.variables.get("f_DPM") * this.variables.get("E_vac")) /
                     (this.variables.get("c") * this.variables.get("V_sys"));

        const a_aether = this.variables.get("f_aether") * 1e-8 * this.variables.get("f_DPM") *
                        (1 + this.variables.get("f_TRZ")) * a_DPM;

        const Ug1_proxy = 1.0;
        const a_u_g4i = this.variables.get("f_sc") * Ug1_proxy * this.variables.get("f_react") *
                       a_DPM / (this.variables.get("E_vac") * this.variables.get("c"));

        const cos_term = 2 * this.variables.get("A") *
                        Math.cos(this.variables.get("k") * this.variables.get("x")) *
                        Math.cos(this.variables.get("omega_osc") * t);

        // Complex exponential: A * exp(i * (k*x - omega*t))
        const phase = this.variables.get("k") * this.variables.get("x") - this.variables.get("omega_osc") * t;
        const real_exp = this.variables.get("A") * Math.cos(phase);
        const exp_factor = (2 * this.variables.get("pi") / 13.8);
        const a_osc = cos_term + exp_factor * real_exp;

        const a_quantum = (this.variables.get("f_quantum") * this.variables.get("E_vac") * a_DPM) /
                         ((this.variables.get("E_vac") / 10) * this.variables.get("c"));
        const a_fluid = (this.variables.get("f_fluid") * this.variables.get("E_vac") * this.variables.get("V")) /
                       ((this.variables.get("E_vac") / 10) * this.variables.get("c"));
        const a_exp = (this.variables.get("f_exp") * this.variables.get("E_vac") * a_DPM) /
                     ((this.variables.get("E_vac") / 10) * this.variables.get("c"));

        return a_aether + a_u_g4i + a_osc + a_quantum + a_fluid + a_exp;
    }

    // Compute SC Integrated: (1 - B / B_crit) * f_sc
    computeSCIntegrated(B) {
        return (1.0 - (B / this.variables.get("B_crit"))) * this.variables.get("f_sc");
    }

    // Full Compressed + Resonance with SC: (compressed + resonance) * SC * (1 + f_TRZ)
    computeCompressedResTerm(t, B) {
        const comp = this.computeCompressedTerm();
        const res = this.computeResonanceTerm(t);
        const sc_int = this.computeSCIntegrated(B);
        const tr_factor = 1.0 + this.variables.get("f_TRZ");
        return (comp + res) * sc_int * tr_factor;
    }

    // Get equation text (descriptive)
    getEquationText() {
        return "Compressed Terms: a_comp = a_DPM + a_THz + a_vac_diff + a_super\n" +
               "Resonance Terms: a_res = a_aether + U_g4i + a_osc + a_quantum + a_fluid + a_exp\n" +
               "Full: g_comp_res = (a_comp + a_res) * SC_int * (1 + f_TRZ)\n" +
               "Where SC_int = (1 - B / B_crit) * f_sc\n" +
               "Special Terms: UQFF compressed/resonance via plasmotic vacuum; no SM; for systems 10-16.\n" +
               "Solutions: Example g_comp_res ~1e-40 m/s² (micro-scale).\n" +
               "Adaptations: Scaled frequencies for nebulae/SMBH/starbirth.";
    }

    // Print all current variables (for debugging/updates)
    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }
}

module.exports = CompressedResonanceUQFFModule;