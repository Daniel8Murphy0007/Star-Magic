// Source19.js - Rings of Relativity (GAL-CLUS-022058s Einstein Ring) Module
// JavaScript implementation of the Master Universal Gravity Equation (MUGE)
// Converted from source19.cpp with full physics fidelity

class RingsOfRelativity {
    constructor(params = {}) {
        // Initialize with default UQFF values or provided parameters
        this.G = params.G || 6.6743e-11; // Gravitational constant
        this.M_sun = params.M_sun || CONSTANTS.SOLAR_MASS;
        this.M = params.mass || 1e14 * this.M_sun; // Lensing mass (kg) - 1e14 Msun galaxy cluster
        this.r = params.radius || 3.086e20; // Einstein radius (m) ~10 kpc
        this.z_lens = params.z_lens || 0.5; // Lens redshift

        // Hubble parameter calculation
        const Hz_kms = 70 * Math.sqrt(0.3 * Math.pow(1 + this.z_lens, 3) + 0.7); // km/s/Mpc
        this.Hz = (Hz_kms * 1000) / 3.086e19; // Convert to s^-1

        this.B = params.magneticField || 1e-5; // Static magnetic field (T)
        this.B_crit = params.B_crit || 1e11; // Critical B field (T)
        this.Lambda = params.Lambda || 1.1e-52; // Cosmological constant
        this.c_light = params.c_light || 3e8; // Speed of light
        this.q_charge = params.q_charge || 1.602e-19; // Charge (proton)
        this.gas_v = params.gas_v || 1e5; // Gas velocity for EM (m/s)
        this.f_TRZ = params.f_TRZ || 0.1; // Time-reversal factor
        this.L_factor = params.L_factor || 0.67; // Lensing factor (D_LS / D_S)

        // Vacuum densities
        this.rho_vac_UA = params.rho_vac_UA || 7.09e-36; // UA vacuum density (J/m^3)
        this.rho_vac_SCm = params.rho_vac_SCm || 7.09e-37; // SCm vacuum density (J/m^3)
        this.scale_EM = params.scale_EM || 1e-12; // EM scaling factor
        this.proton_mass = params.proton_mass || 1.673e-27; // Proton mass

        // Full MUGE terms
        this.hbar = params.hbar || 1.0546e-34; // Reduced Planck's constant
        this.t_Hubble = params.t_Hubble || 13.8e9 * 3.156e7; // Hubble time (s)
        this.t_Hubble_gyr = params.t_Hubble_gyr || 13.8; // Hubble time in Gyr
        this.delta_x = params.delta_x || 1e-10; // Position uncertainty (m)
        this.delta_p = params.delta_p || this.hbar / this.delta_x; // Momentum uncertainty
        this.integral_psi = params.integral_psi || 1.0; // Wavefunction integral approximation
        this.rho_fluid = params.rho_fluid || 1e-21; // Fluid density (kg/m^3)
        this.A_osc = params.A_osc || 1e-12; // Oscillatory amplitude (m/s^2)
        this.k_osc = params.k_osc || 1.0 / this.r; // Wave number (1/m)
        this.omega_osc = params.omega_osc || 2 * Math.PI / (this.r / this.c_light); // Angular frequency
        this.x_pos = params.x_pos || this.r; // Position for oscillation (m)
        this.M_DM_factor = params.M_DM_factor || 0.1; // Dark matter mass fraction
        this.delta_rho_over_rho = params.delta_rho_over_rho || 1e-5; // Density perturbation
        this.rho_wind = params.rho_wind || 1e-21; // Wind density (kg/m^3)
        this.v_wind = params.v_wind || 2e6; // Wind velocity (m/s)

        // Update cached values
        this.updateCache();
    }

    // Update cached computations for efficiency
    updateCache() {
        this.ug1_base = (this.G * this.M) / (this.r * this.r);
        this.L_t = ((this.G * this.M) / (Math.pow(this.c_light, 2) * this.r)) * this.L_factor;
    }

    // Universal setter for any variable
    setVariable(varName, newValue) {
        if (this.hasOwnProperty(varName)) {
            this[varName] = newValue;
            this.updateCache();
            return true;
        }
        return false;
    }

    // Addition method for variables
    addToVariable(varName, delta) {
        if (this.hasOwnProperty(varName)) {
            this[varName] += delta;
            this.updateCache();
            return true;
        }
        return false;
    }

    // Subtraction method for variables
    subtractFromVariable(varName, delta) {
        return this.addToVariable(varName, -delta);
    }

    // Getter for any variable
    getVariable(varName) {
        return this[varName] || 0;
    }

    // Ug terms computation (Universal Gravity components)
    compute_Ug(Mt) {
        const Ug1 = this.ug1_base;
        const Ug2 = 0.0; // No time derivative term
        const Ug3 = 0.0; // No moon term
        const corr_B = 1 - this.B / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ);
    }

    // Volume computation for fluid dynamics
    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    // Main MUGE computation with ALL terms
    compute_g_Rings(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        // Term 1: Base gravity + Hubble + Magnetic + Lensing corrections
        const corr_H = 1 + this.Hz * t;
        const corr_B = 1 - this.B / this.B_crit;
        const corr_L = 1 + this.L_t;
        const term1 = this.ug1_base * corr_H * corr_B * corr_L;

        // Term 2: UQFF Universal Gravity with time-reversal factor
        const term2 = this.compute_Ug(0); // No mass variation for lensing

        // Term 3: Cosmological constant (Lambda)
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled electromagnetic with UA enhancement
        const cross_vB = this.gas_v * this.B; // Magnitude assuming perpendicular
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        const term4 = (em_base * corr_UA) * this.scale_EM;

        // Quantum uncertainty term
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Fluid dynamics term (effective acceleration)
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * this.ug1_base) / this.M;

        // Oscillatory terms (real parts of complex exponentials)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // Dark matter and density perturbation term
        const M_dm = this.M * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * this.M / (this.r * this.r * this.r);
        const term_dm_force_like = (this.M + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / this.M;

        // Stellar wind feedback term (pressure/density for acceleration)
        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        const term_wind = wind_pressure / this.rho_fluid;

        // Total g_Rings (sum of all MUGE components)
        const g_Rings = term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind;

        // Return comprehensive results
        return {
            g_Rings: g_Rings,
            components: {
                term1: term1,      // Base + Hubble + B + Lensing
                term2: term2,      // Universal Gravity (Ug)
                term3: term3,      // Dark Energy (Lambda)
                term4: term4,      // Electromagnetic + UA
                term_q: term_q,    // Quantum Uncertainty
                term_fluid: term_fluid, // Fluid Dynamics
                term_osc: term_osc,     // Oscillatory Waves
                term_DM: term_DM,       // Dark Matter/Density
                term_wind: term_wind    // Stellar Wind Feedback
            },
            diagnostics: {
                hubbleCorrection: corr_H,
                magneticCorrection: corr_B,
                lensingCorrection: corr_L,
                lensingFactor: this.L_factor,
                einsteinRadius: this.r,
                redshift: this.z_lens,
                lensMass: this.M,
                ug1_base: this.ug1_base,
                L_t: this.L_t
            }
        };
    }

    // Print parameters for debugging
    printParameters() {
        console.log("Rings of Relativity Parameters:");
        console.log(`Mass: ${(this.M / this.M_sun).toExponential(2)} Msun (Galaxy Cluster)`);
        console.log(`Einstein Radius: ${(this.r / 1000).toFixed(0)} km (~10 kpc)`);
        console.log(`Redshift z: ${this.z_lens}`);
        console.log(`Hubble Parameter Hz: ${this.Hz.toExponential(2)} s^-1`);
        console.log(`Lensing Factor L_factor: ${this.L_factor}`);
        console.log(`Lensing Amplification L_t: ${this.L_t.toExponential(2)}`);
        console.log(`Cluster Gas Density: ${this.rho_fluid.toExponential(2)} kg/m³`);
        console.log(`Gas Velocity: ${(this.gas_v / 1e5).toFixed(1)} × 10^5 m/s`);
        console.log(`Magnetic Field: ${this.B.toExponential(2)} T`);
        console.log(`Time-Reversal Factor f_TRZ: ${this.f_TRZ}`);
    }
}

module.exports = { RingsOfRelativity };