// Source13.js - SGR 1745-2900 Magnetar Module (JavaScript Implementation)
// Based on source13.cpp - Master Universal Gravity Equation (MUGE) for SGR 1745-2900
// Includes ALL terms: base gravity, cosmic expansion, BH influence, UQFF Ug, Lambda,
// EM, GW, quantum uncertainty, fluid dynamics, oscillatory waves, DM/density perturbations,
// magnetic energy, and decay energy

class MagnetarSGR1745_2900 {
    constructor() {
        this.initializeDefaults();
    }

    initializeDefaults() {
        // Core UQFF parameters (matching C++ implementation)
        this.G = 6.6743e-11;
        this.M = 1.4 * 1.989e30; // 1.4 M_sun in kg
        this.r = 1e4; // 10 km radius
        this.Hz = 2.269e-18; // H(z) in s^-1
        this.B0 = 2e10; // Initial B field
        this.B = this.B0; // Current B field (static model)
        this.tau_B = 4000 * 3.15576e7; // B decay timescale (not used)
        this.B_crit = 1e11; // Critical B field
        this.Lambda = 1.1e-52; // Cosmological constant
        this.c_light = 3e8; // Speed of light
        this.q_charge = 1.602e-19; // Proton charge
        this.v_surf = 1e6; // Surface velocity
        this.f_sc = 1 - (this.B / this.B_crit); // Superconductive factor
        this.rho_vac_UA = 7.09e-36; // UA vacuum density
        this.rho_vac_SCm = 7.09e-37; // SCm vacuum density
        this.P_init = 3.76; // Pulse period in seconds
        this.tau_Omega = 10000 * 3.15576e7; // Omega decay timescale
        this.scale_EM = 1e-12; // EM scaling factor
        this.proton_mass = 1.673e-27; // Proton mass
        this.M_BH = 4e6 * 1.989e30; // Sgr A* mass
        this.r_BH = 2.83e16; // Distance to Sgr A*
        this.mu0 = 4 * Math.PI * 1e-7; // Vacuum permeability
        this.L0_W = 5e28; // Initial luminosity (W)
        this.tau_decay = 3.5 * 365.25 * 24 * 3600; // 3.5 years in seconds

        // Full term parameters
        this.hbar = 1.0546e-34; // Reduced Planck's constant
        this.t_Hubble = 13.8e9 * 3.15576e7; // Hubble time in seconds
        this.t_Hubble_gyr = 13.8; // Hubble time in Gyr
        this.delta_x = 1e-10; // Position uncertainty
        this.delta_p = this.hbar / this.delta_x; // Momentum uncertainty
        this.integral_psi = 1.0; // Wavefunction integral
        this.rho_fluid = 1e17; // Fluid density
        this.A_osc = 1e10; // Oscillatory amplitude
        this.k_osc = 1.0 / this.r; // Wave number
        this.omega_osc = 2 * Math.PI / this.P_init; // Angular frequency
        this.x_pos = this.r; // Position for oscillation
        this.M_DM_factor = 0.1; // Dark matter mass fraction
        this.delta_rho_over_rho = 1e-5; // Density perturbation

        this.updateCache();
    }

    updateCache() {
        this.ug1_base = (this.G * this.M) / (this.r * this.r);
        this.f_sc = 1 - (this.B / this.B_crit);
    }

    setVariable(varName, newValue) {
        if (this.hasOwnProperty(varName)) {
            this[varName] = newValue;
            if (varName === 'B0') this.B = newValue;
            this.updateCache();
            return true;
        }
        return false;
    }

    getVariable(varName) {
        return this.hasOwnProperty(varName) ? this[varName] : 0.0;
    }

    addToVariable(varName, delta) {
        return this.setVariable(varName, this.getVariable(varName) + delta);
    }

    subtractFromVariable(varName, delta) {
        return this.addToVariable(varName, -delta);
    }

    B_t(t) {
        return this.B; // Static B field for this model
    }

    Omega_t(t) {
        return (2 * Math.PI / this.P_init) * Math.exp(-t / this.tau_Omega);
    }

    dOmega_dt(t) {
        const omega0 = 2 * Math.PI / this.P_init;
        return omega0 * (-1.0 / this.tau_Omega) * Math.exp(-t / this.tau_Omega);
    }

    compute_Ug() {
        const Ug1 = this.ug1_base;
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const Ug4 = Ug1 * this.f_sc;
        return Ug1 + Ug2 + Ug3 + Ug4;
    }

    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    compute_M_mag() {
        const V = this.compute_V();
        return (this.B_t(0) * this.B_t(0) / (2 * this.mu0)) * V;
    }

    compute_cumulative_D(t) {
        const exp_term = Math.exp(-t / this.tau_decay);
        return this.L0_W * this.tau_decay * (1 - exp_term);
    }

    compute_g_Magnetar(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        const Bt = this.B_t(t);
        const dOdt = this.dOmega_dt(t);

        // f_sc update
        const current_f_sc = 1 - (Bt / this.B_crit);

        // Term 1: Base + H(z) + B corrections
        const corr_H = 1 + this.Hz * t;
        const corr_B = current_f_sc;
        const term1 = this.ug1_base * corr_H * corr_B;

        // BH term
        const term_BH = (this.G * this.M_BH) / (this.r_BH * this.r_BH);

        // Term 2: UQFF Ug
        const term2 = this.compute_Ug();

        // Term 3: Lambda
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled EM (v x B magnitude)
        const cross_vB = this.v_surf * Bt;
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const term4 = em_base * this.scale_EM;

        // Term 5: GW
        const gw_prefactor = (this.G * this.M * this.M) / (Math.pow(this.c_light, 4) * this.r);
        const term5 = gw_prefactor * (dOdt * dOdt);

        // Quantum uncertainty term
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Fluid term (effective acceleration)
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * this.ug1_base) / this.M;

        // Oscillatory terms (real parts)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term (converted to acceleration)
        const M_dm = this.M * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * this.M / (this.r * this.r * this.r);
        const term_dm_force_like = (this.M + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / this.M;

        // Magnetic energy term (effective g)
        const M_mag = this.compute_M_mag();
        const term_mag = M_mag / (this.M * this.r);

        // Decay term (cumulative energy effective g)
        const cum_D = this.compute_cumulative_D(t);
        const term_decay = cum_D / (this.M * this.r);

        // Total g_Magnetar (all terms summed)
        return term1 + term_BH + term2 + term3 + term4 + term5 + term_q + term_fluid + term_osc + term_DM + term_mag + term_decay;
    }

    printParameters() {
        console.log("SGR 1745-2900 Parameters:");
        console.log(`G: ${this.G}, M: ${this.M}, r: ${this.r}`);
        console.log(`Hz: ${this.Hz}, B: ${this.B}, M_BH: ${this.M_BH}, r_BH: ${this.r_BH}`);
        console.log(`L0_W: ${this.L0_W}, tau_decay: ${this.tau_decay}`);
        console.log(`f_sc: ${this.f_sc}, rho_fluid: ${this.rho_fluid}, M_DM_factor: ${this.M_DM_factor}`);
        console.log(`A_osc: ${this.A_osc}, delta_rho_over_rho: ${this.delta_rho_over_rho}`);
        const M_mag = this.compute_M_mag();
        console.log(`M_mag (J): ${M_mag}, ug1_base: ${this.ug1_base}`);
    }

    exampleAtOneYear() {
        const t_example = 1.0 * 365.25 * 24 * 3600;
        return this.compute_g_Magnetar(t_example);
    }

    getParameters() {
        return {
            name: "SGR 1745-2900",
            G: this.G,
            M: this.M,
            r: this.r,
            Hz: this.Hz,
            B: this.B,
            M_BH: this.M_BH,
            r_BH: this.r_BH,
            L0_W: this.L0_W,
            tau_decay: this.tau_decay,
            f_sc: this.f_sc,
            ug1_base: this.ug1_base
        };
    }
}

module.exports = { MagnetarSGR1745_2900 };