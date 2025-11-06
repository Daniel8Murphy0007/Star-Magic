// CrabResonanceUQFFModule.js
// JavaScript implementation of the Master Universal Gravity Equation (UQFF Resonance) for Crab Nebula Evolution.
// This module implements comprehensive resonance physics with dynamic variable management.
// All variables are stored in a Map for dynamic addition/subtraction/update.
// Includes all resonance terms: DPM resonance, THz pipeline resonance, Aether-mediated resonance,
// U_g4i reactive resonance, quantum resonance, fluid resonance, oscillatory resonance, cosmic expansion resonance.
// Superconductive correction integrated.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Resonance terms use real part of exp; frequencies from pulsar spin/wind.
// Crab params: M=4.6 Msun, r0=5.2e16 m, v_exp=1.5e6 m/s, f_DPM=1e12 Hz (pulsar-aligned).

class CrabResonanceUQFFModule {
    constructor() {
        // Initialize Map for dynamic variables
        this.variables = new Map();

        // Base constants (UQFF universal)
        this.variables.set("c", 3e8);                           // m/s
        this.variables.set("pi", Math.PI);                     // pi
        this.variables.set("E_vac", 7.09e-36);                 // J/m^3 (plasmotic vacuum energy density)
        this.variables.set("hbar", 1.0546e-34);                // J s
        this.variables.set("f_TRZ", 0.1);                      // Time-reversal correction

        // Crab Nebula parameters
        const M_sun_val = 1.989e30;                            // kg
        this.variables.set("M_sun", M_sun_val);
        this.variables.set("M", 4.6 * M_sun_val);              // Total mass kg
        this.variables.set("r0", 5.2e16);                      // m (initial radius)
        this.variables.set("v_exp", 1.5e6);                    // m/s (expansion velocity)

        // Resonance parameters (pulsar-driven)
        this.variables.set("f_DPM", 1e12);                     // Hz (DPM, aligned with 30 Hz pulsar scaled)
        this.variables.set("f_THz", 1e12);                     // Hz (THz hole)
        this.variables.set("f_aether", 1e4);                   // Hz (Aether-mediated)
        this.variables.set("f_react", 1e10);                   // Hz (U_g4i reactive)
        this.variables.set("f_quantum", 1.445e-17);            // Hz (quantum wave)
        this.variables.set("f_fluid", 1.269e-14);              // Hz (filament fluid)
        this.variables.set("f_exp", 1.373e-8);                 // Hz (expansion)
        this.variables.set("f_osc", 30.2 * 60);                // Hz (pulsar 30.2 Hz * 60 for res scale)
        this.variables.set("I", 1e21);                         // A (current proxy from wind)
        this.variables.set("A_vort", 3.142e8);                 // m^2 (vortical area proxy)
        this.variables.set("omega_1", 1e-3);                   // rad/s
        this.variables.set("omega_2", -1e-3);                  // rad/s
        this.variables.set("E_0", 6.381e-36);                  // J/m^3
        this.variables.set("f_vac_diff", 0.143);               // Hz
        this.variables.set("V_sys", 4.189e12);                 // m^3 (proxy)

        // Superconductive resonance integrated
        this.variables.set("B_crit", 1e11);                    // T
        this.variables.set("f_sc", 1.0);                       // Factor

        // Oscillatory/resonant
        this.variables.set("k", 1e20);                         // m^-1
        this.variables.set("omega_osc", 1e15);                 // rad/s (synchrotron scale)
        this.variables.set("x", 0.0);                          // m
        this.variables.set("A", 1e-10);                        // Amplitude

        // Fluid/DM proxies
        this.variables.set("rho_fluid", 1e-21);                // kg/m^3 (filaments)
        this.variables.set("V", 1e3);                          // m^3
        this.variables.set("delta_rho", 0.1 * this.variables.get("rho_fluid"));
        this.variables.set("rho", this.variables.get("rho_fluid"));

        // Quantum
        this.variables.set("Delta_x", 1e-10);                  // m
        this.variables.set("Delta_p", this.variables.get("hbar") / this.variables.get("Delta_x"));
        this.variables.set("integral_psi", 1.0);
    }

    // Dynamic variable operations
    updateVariable(name, value) {
        this.variables.set(name, value);
        // Auto-update dependent variables
        if (name === "Delta_x") {
            this.variables.set("Delta_p", this.variables.get("hbar") / value);
        }
    }

    addToVariable(name, delta) {
        const current = this.variables.get(name) || 0;
        this.variables.set(name, current + delta);
        // Auto-update dependent variables
        if (name === "Delta_x") {
            this.variables.set("Delta_p", this.variables.get("hbar") / this.variables.get(name));
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    getVariable(name) {
        return this.variables.get(name) || 0;
    }

    // Compute DPM Resonance Term: a_DPM_res = (F_DPM * f_DPM * E_vac) / (c * V_sys)
    computeDPMResTerm() {
        const F_DPM = this.variables.get("I") * this.variables.get("A_vort") *
                     (this.variables.get("omega_1") - this.variables.get("omega_2"));
        const t = this.variables.get("t") || 0;
        const r_t = this.variables.get("r0") + this.variables.get("v_exp") * t;  // r(t) proxy
        const V_sys_t = (4.0 / 3.0) * this.variables.get("pi") * Math.pow(r_t, 3);  // Updated volume
        return (F_DPM * this.variables.get("f_DPM") * this.variables.get("E_vac")) /
               (this.variables.get("c") * V_sys_t);
    }

    // Compute THz Resonance Term: a_THz_res = (f_THz * E_vac * v_exp * a_DPM_res) / (E_vac / 10 * c)
    computeTHzResTerm() {
        const a_DPM_res = this.computeDPMResTerm();
        return (this.variables.get("f_THz") * this.variables.get("E_vac") * this.variables.get("v_exp") * a_DPM_res) /
               ((this.variables.get("E_vac") / 10.0) * this.variables.get("c"));
    }

    // Compute Aether Resonance Term: a_aether_res = f_aether * 1e-8 * f_DPM * (1 + f_TRZ) * a_DPM_res
    computeAetherResTerm() {
        const a_DPM_res = this.computeDPMResTerm();
        return this.variables.get("f_aether") * 1e-8 * this.variables.get("f_DPM") *
               (1 + this.variables.get("f_TRZ")) * a_DPM_res;
    }

    // Compute U_g4i Reactive Resonance Term: U_g4i_res = f_sc * Ug1_proxy * f_react * a_DPM_res / (E_vac * c)
    computeU_g4iResTerm() {
        const Ug1_proxy = 1.0;  // Normalized
        const a_DPM_res = this.computeDPMResTerm();
        return this.variables.get("f_sc") * Ug1_proxy * this.variables.get("f_react") * a_DPM_res /
               (this.variables.get("E_vac") * this.variables.get("c"));
    }

    // Compute Quantum Resonance Term: a_quantum_res = (f_quantum * E_vac * a_DPM_res) / ((E_vac / 10) * c)
    computeQuantumResTerm() {
        const a_DPM_res = this.computeDPMResTerm();
        return (this.variables.get("f_quantum") * this.variables.get("E_vac") * a_DPM_res) /
               ((this.variables.get("E_vac") / 10.0) * this.variables.get("c"));
    }

    // Compute Fluid Resonance Term: a_fluid_res = (f_fluid * E_vac * V * a_DPM_res) / ((E_vac / 10) * c)
    computeFluidResTerm() {
        const a_DPM_res = this.computeDPMResTerm();
        return (this.variables.get("f_fluid") * this.variables.get("E_vac") * this.variables.get("V") * a_DPM_res) /
               ((this.variables.get("E_vac") / 10.0) * this.variables.get("c"));
    }

    // Compute Oscillatory Resonance Term: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
    computeOscResTerm(t) {
        const cos_term = 2 * this.variables.get("A") *
                        Math.cos(this.variables.get("k") * this.variables.get("x")) *
                        Math.cos(this.variables.get("omega_osc") * t);
        // Complex exponential: exp(i * phase) = cos(phase) + i*sin(phase), take real part
        const phase = this.variables.get("k") * this.variables.get("x") - this.variables.get("omega_osc") * t;
        const real_exp = this.variables.get("A") * Math.cos(phase);
        const exp_factor = (2 * this.variables.get("pi") / 13.8);
        return cos_term + exp_factor * real_exp;
    }

    // Compute Expansion Resonance Term: a_exp_res = (f_exp * E_vac * a_DPM_res) / ((E_vac / 10) * c)
    computeExpResTerm() {
        const a_DPM_res = this.computeDPMResTerm();
        return (this.variables.get("f_exp") * this.variables.get("E_vac") * a_DPM_res) /
               ((this.variables.get("E_vac") / 10.0) * this.variables.get("c"));
    }

    // Compute SC Resonance Integrated: (1 - B / B_crit) * f_sc
    computeSCResIntegrated(B) {
        return (1.0 - (B / this.variables.get("B_crit"))) * this.variables.get("f_sc");
    }

    // Full g_UQFF Resonance: Sum resonance terms * SC * (1 + f_TRZ)
    computeG(t, B) {
        this.variables.set("t", t);
        const a_DPM_res = this.computeDPMResTerm();
        const a_THz_res = this.computeTHzResTerm();
        const a_aether_res = this.computeAetherResTerm();
        const a_u_g4i_res = this.computeU_g4iResTerm();
        const a_quantum_res = this.computeQuantumResTerm();
        const a_fluid_res = this.computeFluidResTerm();
        const a_osc_res = this.computeOscResTerm(t);
        const a_exp_res = this.computeExpResTerm();
        const sc_int = this.computeSCResIntegrated(B);
        const tr_factor = 1.0 + this.variables.get("f_TRZ");
        const res_sum = a_DPM_res + a_THz_res + a_aether_res + a_u_g4i_res +
                       a_quantum_res + a_fluid_res + a_osc_res + a_exp_res;
        return res_sum * sc_int * tr_factor;
    }

    // Get equation text (descriptive)
    getEquationText() {
        return "g_Crab_Res(t, B) = [a_DPM_res + a_THz_res + a_aether_res + U_g4i_res + a_quantum_res + a_fluid_res + a_osc_res + a_exp_res] * SC_int * (1 + f_TRZ)\n" +
               "Where:\n" +
               "- a_DPM_res = (F_DPM * f_DPM * E_vac) / (c * V_sys(t)); F_DPM = I * A * (ω1 - ω2); V_sys(t) = 4/3 π r(t)^3, r(t)=r0 + v_exp t\n" +
               "- a_THz_res = (f_THz * E_vac * v_exp * a_DPM_res) / (E_vac/10 * c)\n" +
               "- a_aether_res = f_aether * 1e-8 * f_DPM * (1 + f_TRZ) * a_DPM_res\n" +
               "- U_g4i_res = f_sc * Ug1_proxy * f_react * a_DPM_res / (E_vac * c)\n" +
               "- a_quantum_res = (f_quantum * E_vac * a_DPM_res) / (E_vac/10 * c)\n" +
               "- a_fluid_res = (f_fluid * E_vac * V * a_DPM_res) / (E_vac/10 * c)\n" +
               "- a_osc_res = 2 A cos(k x) cos(ω t) + (2π / 13.8) A Re[exp(i (k x - ω t))]\n" +
               "- a_exp_res = (f_exp * E_vac * a_DPM_res) / (E_vac/10 * c)\n" +
               "- SC_int = (1 - B / B_crit) * f_sc\n" +
               "Special Terms: UQFF resonance via plasmotic vacuum; Aether replaces dark energy; no SM terms; pulsar-driven f_osc.\n" +
               "Solutions: At t=971 yr, B=1e-8 T, g ≈ 1e-40 m/s² (resonance micro-scale, wind proxy).\n" +
               "Adaptations: Resonance focus for Crab wisps/shocks per Hubble/Chandra.";
    }

    // Print variables (for debugging)
    printVariables() {
        console.log("Current Variables:");
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value}`);
        }
    }
}

module.exports = CrabResonanceUQFFModule;