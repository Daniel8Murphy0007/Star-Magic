// Source15.js - Sagittarius A* SMBH Module (JavaScript Implementation)
// Based on source15.cpp - Master Universal Gravity Equation (MUGE) for Sgr A*
// Includes ALL terms: base gravity with mass growth M(t), cosmic expansion (H_0), magnetic decay,
// UQFF Ug components with f_TRZ, Lambda, quantum uncertainty, EM (with B(t)), fluid dynamics,
// oscillatory waves, DM/density perturbations with precession sin(30°), and GW term
// Maintains full dynamics: mass accretion, magnetic field decay, spin evolution

class SMBHSgrAStar {
    constructor() {
        this.initializeDefaults();
    }

    initializeDefaults() {
        // Core UQFF parameters (matching C++ implementation)
        this.G = 6.6743e-11;
        this.M_initial = 4.3e6 * 1.989e30; // 4.3 million solar masses in kg
        this.r = 1.27e10; // Schwarzschild radius (m)
        this.H0 = 2.184e-18; // Hubble constant (s^-1)
        this.B0_G = 1e4; // Initial magnetic field (Gauss)
        this.tau_B = 1e6 * 3.15576e7; // B decay timescale (1 million years)
        this.B_crit = 1e11; // Critical B field (T)
        this.Lambda = 1.1e-52; // Cosmological constant
        this.c_light = 3e8; // Speed of light
        this.q_charge = 1.602e-19; // Proton charge
        this.v_surf = 1e6; // Surface velocity (arbitrary for BH)
        this.f_TRZ = 0.1; // Time-reversal factor
        this.M_dot_0 = 0.01; // Initial mass accretion rate factor
        this.tau_acc = 9e9 * 3.15576e7; // Accretion timescale (9 billion years)
        this.spin_factor = 0.3; // Spin factor
        this.tau_Omega = 9e9 * 3.15576e7; // Omega decay timescale

        // Full term parameters
        this.hbar = 1.0546e-34; // Reduced Planck's constant
        this.t_Hubble = 13.8e9 * 3.15576e7; // Hubble time in seconds
        this.t_Hubble_gyr = 13.8; // Hubble time in Gyr
        this.delta_x = 1e-10; // Position uncertainty
        this.delta_p = this.hbar / this.delta_x; // Momentum uncertainty
        this.integral_psi = 1.0; // Wavefunction integral
        this.rho_fluid = 1e17; // Fluid density (accretion disk)
        this.A_osc = 1e6; // Oscillatory amplitude (scaled for BH)
        this.k_osc = 1.0 / this.r; // Wave number
        this.omega_osc = 2 * Math.PI / (this.r / this.c_light); // Orbital-like frequency
        this.x_pos = this.r; // Position for oscillation
        this.M_DM_factor = 0.1; // Dark matter mass fraction
        this.delta_rho_over_rho = 1e-5; // Density perturbation
        this.precession_angle_deg = 30.0; // Precession angle (degrees)

        this.updateCache();
    }

    updateCache() {
        this.ug1_base = (this.G * this.M_initial) / (this.r * this.r);
    }

    setVariable(varName, newValue) {
        if (this.hasOwnProperty(varName)) {
            this[varName] = newValue;
            if (varName === 'B0_G') this.B = newValue;
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

    M_t(t) {
        const M_dot = this.M_dot_0 * Math.exp(-t / this.tau_acc);
        return this.M_initial * (1 + M_dot);
    }

    B_t(t) {
        const B_G = this.B0_G * Math.exp(-t / this.tau_B);
        return B_G * 1e-4; // Convert Gauss to Tesla
    }

    Omega_t(t) {
        const omega0 = this.spin_factor * this.c_light / this.r;
        return omega0 * Math.exp(-t / this.tau_Omega);
    }

    dOmega_dt(t) {
        const omega0 = this.spin_factor * this.c_light / this.r;
        return omega0 * (-1.0 / this.tau_Omega) * Math.exp(-t / this.tau_Omega);
    }

    compute_Ug(Mt, Bt) {
        const Ug1 = (this.G * Mt) / (this.r * this.r);
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const corr_B = 1 - Bt / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ);
    }

    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    compute_g_SgrA(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        const Mt = this.M_t(t);
        const Bt = this.B_t(t);
        const dOdt = this.dOmega_dt(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);

        // Term 1: Base + H0 + B corrections
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - Bt / this.B_crit;
        const term1 = ug1_t * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ
        const term2 = this.compute_Ug(Mt, Bt);

        // Term 3: Lambda
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: EM (v x B, no scaling or UA here)
        const cross_vB = this.v_surf * Bt;
        const em_base = this.q_charge * cross_vB / 1.673e-27; // Acceleration
        const term4 = em_base;

        // Term 5: GW
        const gw_prefactor = (this.G * Mt * Mt) / (Math.pow(this.c_light, 4) * this.r);
        const term5 = gw_prefactor * (dOdt * dOdt);

        // Quantum uncertainty term
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Fluid term (effective acceleration)
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * ug1_t) / Mt;

        // Oscillatory terms (real parts)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term with precession (converted to acceleration)
        const M_dm = Mt * this.M_DM_factor;
        const sin_prec = Math.sin(this.precession_angle_deg * Math.PI / 180.0);
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2 * sin_prec);
        const term_DM = term_dm_force_like / Mt;

        // Total g_SgrA (all terms summed)
        return term1 + term2 + term3 + term4 + term5 + term_q + term_fluid + term_osc + term_DM;
    }

    printParameters() {
        console.log("Sgr A* Parameters:");
        console.log(`G: ${this.G}, M_initial: ${this.M_initial}, r: ${this.r}`);
        console.log(`H0: ${this.H0}, B0_G: ${this.B0_G}, tau_B: ${this.tau_B}`);
        console.log(`f_TRZ: ${this.f_TRZ}, M_dot_0: ${this.M_dot_0}, tau_acc: ${this.tau_acc}`);
        console.log(`rho_fluid: ${this.rho_fluid}, M_DM_factor: ${this.M_DM_factor}`);
        console.log(`A_osc: ${this.A_osc}, precession_angle_deg: ${this.precession_angle_deg}`);
        console.log(`ug1_base: ${this.ug1_base}`);
    }

    exampleAt4_5Gyr() {
        const t_example = 4.5e9 * 3.15576e7;
        return this.compute_g_SgrA(t_example);
    }

    getParameters() {
        return {
            name: "Sagittarius A*",
            G: this.G,
            M_initial: this.M_initial,
            r: this.r,
            H0: this.H0,
            B0_G: this.B0_G,
            tau_B: this.tau_B,
            f_TRZ: this.f_TRZ,
            M_dot_0: this.M_dot_0,
            tau_acc: this.tau_acc,
            spin_factor: this.spin_factor,
            ug1_base: this.ug1_base
        };
    }
}

module.exports = { SMBHSgrAStar };