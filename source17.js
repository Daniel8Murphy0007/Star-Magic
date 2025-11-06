/**
 * Westerlund 2 Super Star Cluster UQFF Module
 * Converted from source17.cpp
 * Maintains all time-dependent dynamics and MUGE physics
 */

class Westerlund2 {
    constructor() {
        this.initializeDefaults();
    }

    initializeDefaults() {
        // Core parameters from UQFF document
        this.G = 6.6743e-11;                    // Gravitational constant
        this.M_sun = 1.989e30;                  // Solar mass in kg
        this.M_initial_sun = 30000.0;          // Initial mass in solar masses
        this.M_initial = this.M_initial_sun * this.M_sun;  // Initial mass in kg
        this.r = 9.461e16;                      // Radius (10 ly in meters)
        this.H0 = 2.184e-18;                    // Hubble constant (s^-1)
        this.B = 1e-5;                          // Static magnetic field (T)
        this.B_crit = 1e11;                     // Critical B field (T)
        this.Lambda = 1.1e-52;                  // Cosmological constant
        this.c_light = 3e8;                     // Speed of light
        this.q_charge = 1.602e-19;              // Proton charge
        this.gas_v = 1e5;                       // Gas velocity for EM (m/s)
        this.f_TRZ = 0.1;                       // Time-reversal factor
        this.M_dot_factor = 1e5 / this.M_initial_sun;  // Star formation factor
        this.tau_SF = 2e6 * 3.156e7;            // Star formation timescale (2 Myr)
        this.rho_wind = 1e-20;                  // Wind density (kg/m^3)
        this.v_wind = 2e6;                      // Wind velocity (m/s)
        this.rho_fluid = 1e-20;                 // Fluid density (kg/m^3)
        this.rho_vac_UA = 7.09e-36;             // UA vacuum density
        this.rho_vac_SCm = 7.09e-37;            // SCm vacuum density
        this.scale_EM = 1e-12;                  // EM scaling factor
        this.proton_mass = 1.673e-27;           // Proton mass

        // Full MUGE terms
        this.hbar = 1.0546e-34;                 // Reduced Planck's constant
        this.t_Hubble = 13.8e9 * 3.156e7;       // Hubble time (s)
        this.t_Hubble_gyr = 13.8;               // Hubble time in Gyr
        this.delta_x = 1e-10;                   // Position uncertainty (m)
        this.delta_p = this.hbar / this.delta_x; // Momentum uncertainty
        this.integral_psi = 1.0;                // Wavefunction integral
        this.A_osc = 1e-9;                      // Oscillatory amplitude
        this.k_osc = 1.0 / this.r;              // Wave number
        this.omega_osc = 2 * Math.PI / (this.r / this.c_light); // Angular frequency
        this.x_pos = this.r;                    // Position for oscillation
        this.M_DM_factor = 0.1;                 // Dark matter mass fraction
        this.delta_rho_over_rho = 1e-5;         // Density perturbation

        this.updateCache();
    }

    updateCache() {
        this.ug1_base = (this.G * this.M_initial) / (this.r * this.r);
    }

    // Dynamic variable management
    setVariable(varName, newValue) {
        if (this.hasOwnProperty(varName)) {
            this[varName] = newValue;
            if (varName === 'delta_x') {
                this.delta_p = this.hbar / this.delta_x;
            }
            this.updateCache();
            return true;
        }
        return false;
    }

    addToVariable(varName, delta) {
        const current = this[varName] || 0;
        return this.setVariable(varName, current + delta);
    }

    subtractFromVariable(varName, delta) {
        return this.addToVariable(varName, -delta);
    }

    // Time-dependent mass growth M(t)
    M_t(t) {
        const M_dot = this.M_dot_factor * Math.exp(-t / this.tau_SF);
        return this.M_initial * (1 + M_dot);
    }

    // UQFF Ug components
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
    compute_g_Westerlund2(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative");
            return 0.0;
        }

        const Mt = this.M_t(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);

        // Term 1: Base gravity + H0 + B corrections
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - this.B / this.B_crit;
        const term1 = ug1_t * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ
        const term2 = this.compute_Ug(Mt);

        // Term 3: Cosmological constant Lambda
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled EM with UA corrections
        const cross_vB = this.gas_v * this.B;
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        const term4 = (em_base * corr_UA) * this.scale_EM;

        // Quantum uncertainty term
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Fluid dynamics term
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * ug1_t) / Mt;

        // Oscillatory waves (real parts)
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

        // Total MUGE result
        const g_Westerlund2 = term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind;

        return {
            g_Westerlund2: g_Westerlund2,
            components: {
                term1: term1,      // Base + corrections
                term2: term2,      // UQFF Ug
                term3: term3,      // Lambda
                term4: term4,      // EM
                term_q: term_q,    // Quantum
                term_fluid: term_fluid,  // Fluid
                term_osc: term_osc,      // Oscillatory
                term_DM: term_DM,        // DM/density
                term_wind: term_wind     // Wind feedback
            },
            diagnostics: {
                Mt: Mt,
                M_solar: Mt / this.M_sun,
                corr_H: corr_H,
                corr_B: corr_B,
                corr_UA: corr_UA,
                massGrowthFactor: Mt / this.M_initial,
                starFormationTimescale: this.tau_SF / (3.156e7 * 1e6), // Myr
                windPressure: wind_pressure,
                time_years: t / (365.25 * 24 * 3600)
            }
        };
    }

    // Example computation for testing
    exampleAt1Myr() {
        const t_example = 1e6 * 3.156e7; // 1 Myr in seconds
        return this.compute_g_Westerlund2(t_example);
    }
}

module.exports = Westerlund2;