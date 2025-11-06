/**
 * ================================================================================================
 * JavaScript Module: source16.js - Tapestry of Blazing Starbirth (NGC 2014 & NGC 2020)
 *
 * Description: Converted from source16.cpp for Star-Magic UQFF Framework
 *              Implements Master Universal Gravity Equation (MUGE) for star-forming region evolution
 *              in the Large Magellanic Cloud. Maintains ALL time-dependent dynamics including:
 *              - M(t) star formation mass growth
 *              - Cosmic expansion (H_0)
 *              - Magnetic field corrections
 *              - UQFF Ug components with f_TRZ
 *              - Lambda cosmological constant
 *              - Quantum uncertainty terms
 *              - Scaled electromagnetic interactions
 *              - Fluid dynamics
 *              - Oscillatory waves
 *              - Dark matter/density perturbations
 *              - Stellar wind feedback
 *
 * Integration: Used in index.js for cross-system UQFF analysis
 *              Instantiate: const starbirth = new StarbirthTapestry()
 *              Compute: g = starbirth.compute_g_Starbirth(t)
 *
 * Author: Converted from C++ source16.cpp for JavaScript UQFF framework
 * Date: November 5, 2025
 * ================================================================================================
 */

class StarbirthTapestry {
    constructor() {
        // Initialize with UQFF default values
        this.initializeDefaults();
    }

    initializeDefaults() {
        // Core gravitational parameters
        this.G = 6.6743e-11; // Gravitational constant

        // System parameters (NGC 2014 & NGC 2020)
        const M_sun = 1.989e30;
        const M_initial_sun = 240.0; // Initial mass in solar masses
        this.M_initial = M_initial_sun * M_sun; // Convert to kg

        const ly_to_m = 9.461e15;
        this.r = 10.0 * ly_to_m; // 10 light years radius

        // Cosmological parameters
        this.H0 = 2.184e-18; // Hubble constant (s^-1)
        this.Lambda = 1.1e-52; // Cosmological constant

        // Magnetic field parameters
        this.B = 1e-6; // Static magnetic field (T)
        this.B_crit = 1e11; // Critical B field (T)

        // Fundamental constants
        this.c_light = 3e8; // Speed of light (m/s)
        this.q_charge = 1.602e-19; // Proton charge (C)
        this.proton_mass = 1.673e-27; // Proton mass (kg)
        this.hbar = 1.0546e-34; // Reduced Planck's constant

        // Time-reversal and scaling factors
        this.f_TRZ = 0.1; // Time-reversal factor
        this.scale_EM = 1e-12; // EM scaling factor

        // Star formation dynamics
        const gas_mass_sun = 10000.0;
        this.M_dot_factor = gas_mass_sun / M_initial_sun; // Star formation factor
        this.tau_SF = 5e6 * 3.156e7; // Star formation timescale (5 Myr in seconds)

        // Gas and wind parameters
        this.gas_v = 1e5; // Gas velocity for EM (m/s)
        this.rho_wind = 1e-21; // Wind density (kg/m^3)
        this.v_wind = 2e6; // Wind velocity (m/s)
        this.rho_fluid = 1e-21; // Fluid density (kg/m^3)

        // Vacuum energy densities
        this.rho_vac_UA = 7.09e-36; // UA vacuum density (J/m^3)
        this.rho_vac_SCm = 7.09e-37; // SCm vacuum density (J/m^3)

        // Quantum uncertainty parameters
        this.delta_x = 1e-10; // Position uncertainty (m)
        this.delta_p = this.hbar / this.delta_x; // Momentum uncertainty
        this.integral_psi = 1.0; // Wavefunction integral approximation

        // Oscillatory wave parameters
        this.A_osc = 1e-10; // Oscillatory amplitude (m/s^2)
        this.k_osc = 1.0 / this.r; // Wave number (1/m)
        this.omega_osc = 2 * Math.PI / (this.r / this.c_light); // Angular frequency (rad/s)
        this.x_pos = this.r; // Position for oscillation (m)

        // Hubble time parameters
        this.t_Hubble = 13.8e9 * 3.156e7; // Hubble time in seconds
        this.t_Hubble_gyr = 13.8; // Hubble time in Gyr

        // Dark matter parameters
        this.M_DM_factor = 0.1; // Dark matter mass fraction
        this.delta_rho_over_rho = 1e-5; // Density perturbation fraction

        // Update cached values
        this.updateCache();
    }

    updateCache() {
        this.ug1_base = (this.G * this.M_initial) / (this.r * this.r);
    }

    // Universal setter for any variable
    setVariable(varName, newValue) {
        if (this.hasOwnProperty(varName)) {
            this[varName] = newValue;
            this.updateCache();
            return true;
        }
        console.error(`Error: Unknown variable '${varName}'`);
        return false;
    }

    // Addition method for variables
    addToVariable(varName, delta) {
        const currentValue = this.getVariable(varName);
        return this.setVariable(varName, currentValue + delta);
    }

    // Subtraction method for variables
    subtractFromVariable(varName, delta) {
        return this.addToVariable(varName, -delta);
    }

    // Getter for any variable
    getVariable(varName) {
        if (this.hasOwnProperty(varName)) {
            return this[varName];
        }
        console.error(`Error: Unknown variable '${varName}'`);
        return 0.0;
    }

    // Time-dependent mass growth M(t)
    M_t(t) {
        const M_dot = this.M_dot_factor * Math.exp(-t / this.tau_SF);
        return this.M_initial * (1 + M_dot);
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

    // Main MUGE computation - Master Universal Gravity Equation
    compute_g_Starbirth(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
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

        // Term 4: Scaled EM with UA vacuum correction
        const cross_vB = this.gas_v * this.B; // Magnitude, assuming perpendicular
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        const term4 = (em_base * corr_UA) * this.scale_EM;

        // Quantum uncertainty term
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Fluid dynamics term (effective acceleration)
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * ug1_t) / Mt;

        // Oscillatory wave terms (real parts)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // Dark matter and density perturbation term
        const M_dm = Mt * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / Mt;

        // Stellar wind feedback term (pressure/density for acceleration)
        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        const term_wind = wind_pressure / this.rho_fluid;

        // Total g_Starbirth (sum of all terms)
        const g_total = term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind;

        return g_total;
    }

    // Debug output method
    printParameters() {
        console.log("Tapestry of Blazing Starbirth Parameters:");
        console.log(`G: ${this.G}, M_initial: ${this.M_initial}, r: ${this.r}`);
        console.log(`H0: ${this.H0}, B: ${this.B}, B_crit: ${this.B_crit}`);
        console.log(`f_TRZ: ${this.f_TRZ}, M_dot_factor: ${this.M_dot_factor}, tau_SF: ${this.tau_SF}`);
        console.log(`rho_fluid: ${this.rho_fluid}, rho_wind: ${this.rho_wind}, v_wind: ${this.v_wind}`);
        console.log(`gas_v: ${this.gas_v}, M_DM_factor: ${this.M_DM_factor}`);
        console.log(`A_osc: ${this.A_osc}, delta_rho_over_rho: ${this.delta_rho_over_rho}`);
        console.log(`ug1_base: ${this.ug1_base}`);
    }

    // Example computation at t=2.5 Myr
    exampleAt2_5Myr() {
        const t_example = 2.5e6 * 3.156e7; // 2.5 Myr in seconds
        return this.compute_g_Starbirth(t_example);
    }
}

module.exports = StarbirthTapestry;