/**
 * PillarsOfCreation.js - JavaScript implementation of Pillars of Creation (Eagle Nebula)
 * Converted from source18.cpp, maintaining all time-dependent dynamics
 */

class PillarsOfCreation {
    constructor() {
        this.initializeDefaults();
    }

    initializeDefaults() {
        // Core parameters (converted from C++)
        this.G = 6.6743e-11;              // Gravitational constant
        this.M_initial = 10100.0 * 1.989e30;  // Initial mass (10,100 Msun)
        this.r = 4.731e16;                // Radius (5 ly)
        this.H0 = 2.184e-18;              // Hubble constant
        this.B = 1e-6;                    // Static magnetic field
        this.B_crit = 1e11;               // Critical B field
        this.Lambda = 1.1e-52;            // Cosmological constant
        this.c_light = 3e8;               // Speed of light
        this.q_charge = 1.602e-19;        // Charge (proton)
        this.gas_v = 1e5;                 // Gas velocity for EM
        this.f_TRZ = 0.1;                 // Time-reversal factor
        this.M_dot_factor = 1e4 / 10100.0; // Star formation factor
        this.tau_SF = 1e6 * 3.156e7;      // Star formation timescale
        this.E_0 = 0.1;                   // Initial erosion factor
        this.tau_erosion = 1e6 * 3.156e7; // Erosion timescale
        this.rho_wind = 1e-21;            // Wind density
        this.v_wind = 2e6;                // Wind velocity
        this.rho_fluid = 1e-21;           // Fluid density
        this.rho_vac_UA = 7.09e-36;       // UA vacuum density
        this.rho_vac_SCm = 7.09e-37;      // SCm vacuum density
        this.scale_EM = 1e-12;            // EM scaling factor
        this.proton_mass = 1.673e-27;     // Proton mass

        // Full terms parameters
        this.hbar = 1.0546e-34;           // Reduced Planck's constant
        this.t_Hubble = 13.8e9 * 3.156e7; // Hubble time
        this.t_Hubble_gyr = 13.8;         // Hubble time in Gyr
        this.delta_x = 1e-10;             // Position uncertainty
        this.delta_p = this.hbar / this.delta_x; // Momentum uncertainty
        this.integral_psi = 1.0;          // Wavefunction integral
        this.A_osc = 1e-10;               // Oscillatory amplitude
        this.k_osc = 1.0 / this.r;        // Wave number
        this.omega_osc = 2 * Math.PI / (this.r / this.c_light); // Angular frequency
        this.x_pos = this.r;              // Position for oscillation
        this.M_DM_factor = 0.1;           // Dark matter mass fraction
        this.delta_rho_over_rho = 1e-5;   // Density perturbation fraction

        this.updateCache();
    }

    updateCache() {
        this.ug1_base = (this.G * this.M_initial) / (this.r * this.r);
    }

    // Time-dependent mass growth M(t)
    M_t(t) {
        const M_dot = this.M_dot_factor * Math.exp(-t / this.tau_SF);
        return this.M_initial * (1 + M_dot);
    }

    // Time-dependent erosion E(t)
    E_t(t) {
        return this.E_0 * Math.exp(-t / this.tau_erosion);
    }

    // Universal gravity components Ug
    compute_Ug(Mt) {
        const Ug1 = (this.G * Mt) / (this.r * this.r);
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const corr_B = 1 - this.B / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ);
    }

    // Volume computation for fluid dynamics
    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    // Main MUGE computation with ALL terms
    compute_g_Pillars(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        const Mt = this.M_t(t);
        const Et = this.E_t(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);

        // Term 1: Base gravity + H0 + B + E corrections
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - this.B / this.B_crit;
        const corr_E = 1 - Et;
        const term1 = ug1_t * corr_H * corr_B * corr_E;

        // Term 2: UQFF Ug with f_TRZ
        const term2 = this.compute_Ug(Mt);

        // Term 3: Cosmological constant Lambda
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled EM with UA vacuum corrections
        const cross_vB = this.gas_v * this.B; // Magnitude, assuming perpendicular
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        const term4 = (em_base * corr_UA) * this.scale_EM;

        // Quantum uncertainty term (Heisenberg principle)
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Fluid dynamics term
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * ug1_t) / Mt;

        // Oscillatory waves (real exponential terms)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // Dark matter and density perturbations
        const M_dm = Mt * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / Mt;

        // Stellar wind feedback (pressure/density acceleration)
        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        const term_wind = wind_pressure / this.rho_fluid;

        // Total MUGE result (all terms summed)
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind;
    }

    // Parameter display for debugging
    printParameters() {
        console.log("Pillars of Creation Parameters:");
        console.log(`M_initial: ${this.M_initial.toExponential(3)} kg (${(this.M_initial / 1.989e30).toFixed(0)} Msun)`);
        console.log(`r: ${this.r.toExponential(3)} m (${(this.r / 9.461e15).toFixed(1)} ly)`);
        console.log(`B: ${this.B.toExponential(3)} T`);
        console.log(`M_dot_factor: ${this.M_dot_factor.toExponential(3)}`);
        console.log(`tau_SF: ${(this.tau_SF / (1e6 * 3.156e7)).toFixed(1)} Myr`);
        console.log(`E_0: ${this.E_0}`);
        console.log(`tau_erosion: ${(this.tau_erosion / (1e6 * 3.156e7)).toFixed(1)} Myr`);
        console.log(`rho_wind: ${this.rho_wind.toExponential(3)} kg/m³`);
        console.log(`v_wind: ${(this.v_wind / 1e6).toFixed(0)} Mm/s`);
    }
}

module.exports = PillarsOfCreation;