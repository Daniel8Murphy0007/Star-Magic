// Source14.js - SGR 0501+4516 Magnetar Module (JavaScript Implementation)
// Based on source14.cpp - Master Universal Gravity Equation (MUGE) for SGR 0501+4516
// Time-Reversal Magnetar with enhanced framework including f_TRZ factor
// Includes ALL terms: base gravity, cosmic expansion, magnetic decay, UQFF Ug with f_TRZ,
// Lambda, scaled EM with UA correction, GW, quantum uncertainty, fluid dynamics,
// oscillatory waves, and DM/density perturbations

class MagnetarSGR0501_4516 {
    constructor() {
        this.initializeDefaults();
    }

    initializeDefaults() {
        // Core UQFF parameters (matching C++ implementation)
        this.G = 6.6743e-11;
        this.M = 1.4 * 1.989e30; // 1.4 M_sun in kg
        this.r = 2e4; // 20 km radius (larger than SGR 1745-2900)
        this.H0 = 2.184e-18; // Different Hubble constant
        this.B0 = 1e10; // Initial B field (weaker than SGR 1745-2900)
        this.tau_B = 4000 * 3.15576e7; // B decay timescale (4000 years)
        this.B_crit = 1e11; // Critical B field
        this.Lambda = 1.1e-52; // Cosmological constant
        this.c_light = 3e8; // Speed of light
        this.q_charge = 1.602e-19; // Proton charge
        this.v_surf = 1e6; // Surface velocity
        this.f_TRZ = 0.1; // Time-reversal factor (unique to SGR 0501+4516)
        this.rho_vac_UA = 7.09e-36; // UA vacuum density
        this.rho_vac_SCm = 7.09e-37; // SCm vacuum density
        this.P_init = 5.0; // Initial rotation period (slower than SGR 1745-2900)
        this.tau_Omega = 10000 * 3.15576e7; // Omega decay timescale
        this.scale_EM = 1e-12; // EM scaling factor
        this.proton_mass = 1.673e-27; // Proton mass

        // Full term parameters
        this.hbar = 1.0546e-34; // Reduced Planck's constant
        this.t_Hubble = 13.8e9 * 3.15576e7; // Hubble time in seconds
        this.delta_x = 1e-10; // Position uncertainty
        this.delta_p = this.hbar / this.delta_x; // Momentum uncertainty
        this.integral_psi = 1.0; // Wavefunction integral
        this.rho_fluid = 1e17; // Fluid density
        this.A_osc = 1e10; // Oscillatory amplitude
        this.k_osc = 1.0 / this.r; // Wave number
        this.omega_osc = 2 * Math.PI / this.P_init; // Angular frequency
        this.x_pos = this.r; // Position for oscillation
        this.t_Hubble_gyr = 13.8; // Hubble time in Gyr
        this.M_DM_factor = 0.1; // Dark matter mass fraction
        this.delta_rho_over_rho = 1e-5; // Density perturbation

        this.updateCache();
    }

    updateCache() {
        this.ug1_base = (this.G * this.M) / (this.r * this.r);
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
        return this.B0 * Math.exp(-t / this.tau_B);
    }

    Omega_t(t) {
        return (2 * Math.PI / this.P_init) * Math.exp(-t / this.tau_Omega);
    }

    dOmega_dt(t) {
        const omega0 = 2 * Math.PI / this.P_init;
        return omega0 * (-1.0 / this.tau_Omega) * Math.exp(-t / this.tau_Omega);
    }

    compute_Ug(Bt) {
        const Ug1 = this.ug1_base;
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const Ug4 = Ug1 * (1 - Bt / this.B_crit);
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ);
    }

    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    compute_g_Magnetar(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        const Bt = this.B_t(t);
        const dOdt = this.dOmega_dt(t);

        // Term 1: Base + H0 + B corrections
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - Bt / this.B_crit;
        const term1 = this.ug1_base * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ (time-reversal enhancement)
        const term2 = this.compute_Ug(Bt);

        // Term 3: Lambda
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled EM with UA correction
        const cross_vB = this.v_surf * Bt;
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        const term4 = (em_base * corr_UA) * this.scale_EM;

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

        // Total g_Magnetar (all terms summed)
        return term1 + term2 + term3 + term4 + term5 + term_q + term_fluid + term_osc + term_DM;
    }

    printParameters() {
        console.log("SGR 0501+4516 Parameters:");
        console.log(`G: ${this.G}, M: ${this.M}, r: ${this.r}`);
        console.log(`H0: ${this.H0}, B0: ${this.B0}, tau_B: ${this.tau_B}`);
        console.log(`f_TRZ: ${this.f_TRZ}, rho_fluid: ${this.rho_fluid}, M_DM_factor: ${this.M_DM_factor}`);
        console.log(`A_osc: ${this.A_osc}, delta_rho_over_rho: ${this.delta_rho_over_rho}`);
        console.log(`ug1_base: ${this.ug1_base}`);
    }

    exampleAt5000Years() {
        const t_example = 5000 * 3.15576e7;
        return this.compute_g_Magnetar(t_example);
    }

    getParameters() {
        return {
            name: "SGR 0501+4516",
            G: this.G,
            M: this.M,
            r: this.r,
            H0: this.H0,
            B0: this.B0,
            tau_B: this.tau_B,
            f_TRZ: this.f_TRZ,
            ug1_base: this.ug1_base
        };
    }
}

module.exports = { MagnetarSGR0501_4516 };