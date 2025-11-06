// Source20.js - Galaxy NGC 2525 Module (Complete MUGE Implementation)
// Enhanced: November 05, 2025 - Full MUGE physics with self-expanding capabilities

class PhysicsTerm {
    constructor() {
        this.dynamicParameters = new Map();
        this.metadata = new Map();
        this.enableLogging = false;
        this.learningRate = 0.001;
    }

    compute(t, params) {
        // Base implementation - override in derived classes
        return 0.0;
    }

    getName() {
        return "BasePhysicsTerm";
    }

    getDescription() {
        return "Base physics term";
    }

    validate(params) {
        return true;
    }

    setDynamicParameter(key, value) {
        this.dynamicParameters.set(key, value);
    }

    getDynamicParameter(key, defaultValue = 0.0) {
        return this.dynamicParameters.get(key) || defaultValue;
    }
}

class DynamicVacuumTerm extends PhysicsTerm {
    constructor(amplitude = 1e-10, frequency = 1e-15) {
        super();
        this.amplitude = amplitude;
        this.frequency = frequency;
    }

    compute(t, params) {
        const rho_vac_UA = params.get('rho_vac_UA') || 7.09e-36;
        return this.amplitude * rho_vac_UA * Math.sin(this.frequency * t);
    }

    getName() {
        return "DynamicVacuum";
    }

    getDescription() {
        return "Time-varying vacuum energy";
    }
}

class QuantumCouplingTerm extends PhysicsTerm {
    constructor(coupling_strength = 1e-40) {
        super();
        this.coupling_strength = coupling_strength;
    }

    compute(t, params) {
        const hbar = params.get('hbar') || 1.0546e-34;
        const M = params.get('M') || 1.989e30;
        const r = params.get('r') || 1e4;
        return this.coupling_strength * (hbar * hbar) / (M * r * r) * Math.cos(t / 1e6);
    }

    getName() {
        return "QuantumCoupling";
    }

    getDescription() {
        return "Non-local quantum effects";
    }
}

class GalaxyNGC2525 {
    constructor() {
        // Self-expanding framework
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;

        this.metadata.set('enhanced', 'true');
        this.metadata.set('version', '2.0-Enhanced');
        this.metadata.set('galaxy', 'NGC 2525');
        this.metadata.set('morphology', 'Barred Spiral');

        this.initializeDefaults();
    }

    initializeDefaults() {
        // Core UQFF parameters for NGC 2525
        this.G = 6.6743e-11; // Gravitational constant
        const M_sun = 1.989e30;
        this.M = (1e10 + 2.25e7) * M_sun; // Total galaxy mass (kg)
        this.r = 2.836e20; // Galaxy radius (m)
        this.z_gal = 0.016; // Galaxy redshift
        const Hz_kms = 70 * Math.sqrt(0.3 * Math.pow(1 + this.z_gal, 3) + 0.7); // km/s/Mpc
        this.Hz = (Hz_kms * 1000) / 3.086e19; // Hubble parameter (s^-1)
        this.B = 1e-5; // Static magnetic field (T)
        this.B_crit = 1e11; // Critical B field (T)
        this.Lambda = 1.1e-52; // Cosmological constant
        this.c_light = 3e8; // Speed of light
        this.q_charge = 1.602e-19; // Charge (proton)
        this.gas_v = 1e5; // Gas velocity for EM (m/s)
        this.f_TRZ = 0.1; // Time-reversal factor
        this.M_BH = 2.25e7 * M_sun; // Black hole mass (kg)
        this.r_BH = 1.496e11; // Black hole influence radius (m)
        this.M_SN0 = 1.4 * M_sun; // Initial SN mass (kg)
        this.tau_SN = 1 * 3.156e7; // SN decay timescale (s)
        this.rho_vac_UA = 7.09e-36; // UA vacuum density (J/m^3)
        this.rho_vac_SCm = 7.09e-37; // SCm vacuum density (J/m^3)
        this.scale_EM = 1e-12; // EM scaling factor
        this.proton_mass = 1.673e-27; // Proton mass for EM acceleration

        // Full MUGE terms parameters
        this.hbar = 1.0546e-34; // Reduced Planck's constant
        this.t_Hubble = 13.8e9 * 3.156e7; // Hubble time (s)
        this.t_Hubble_gyr = 13.8; // Hubble time in Gyr
        this.delta_x = 1e-10; // Position uncertainty (m)
        this.delta_p = this.hbar / this.delta_x; // Momentum uncertainty (kg m/s)
        this.integral_psi = 1.0; // Wavefunction integral approximation
        this.rho_fluid = 1e-21; // Fluid density (kg/m^3)
        this.A_osc = 1e-10; // Oscillatory amplitude (m/s^2)
        this.k_osc = 1.0 / this.r; // Wave number (1/m)
        this.omega_osc = 2 * Math.PI / (this.r / this.c_light); // Angular frequency (rad/s)
        this.x_pos = this.r; // Position for oscillation (m)
        this.M_DM_factor = 0.1; // Dark matter mass fraction
        this.delta_rho_over_rho = 1e-5; // Density perturbation fraction

        this.updateCache();
    }

    updateCache() {
        this.ug1_base = (this.G * this.M) / (this.r * this.r);
        this.g_BH = (this.G * this.M_BH) / (this.r_BH * this.r_BH);
    }

    setVariable(varName, newValue) {
        const validVars = [
            'G', 'M', 'r', 'Hz', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'M_BH', 'r_BH', 'M_SN0', 'tau_SN', 'rho_vac_UA',
            'rho_vac_SCm', 'scale_EM', 'proton_mass', 'z_gal', 'hbar', 't_Hubble',
            't_Hubble_gyr', 'delta_x', 'delta_p', 'integral_psi', 'rho_fluid',
            'A_osc', 'k_osc', 'omega_osc', 'x_pos', 'M_DM_factor', 'delta_rho_over_rho'
        ];

        if (validVars.includes(varName)) {
            this[varName] = newValue;
            this.updateCache();
            return true;
        } else {
            console.error(`Error: Unknown variable '${varName}'.`);
            return false;
        }
    }

    addToVariable(varName, delta) {
        const currentValue = this.getVariable(varName);
        return this.setVariable(varName, currentValue + delta);
    }

    subtractFromVariable(varName, delta) {
        return this.addToVariable(varName, -delta);
    }

    getVariable(varName) {
        const validVars = [
            'G', 'M', 'r', 'Hz', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'M_BH', 'r_BH', 'M_SN0', 'tau_SN', 'rho_vac_UA',
            'rho_vac_SCm', 'scale_EM', 'proton_mass', 'z_gal', 'hbar', 't_Hubble',
            't_Hubble_gyr', 'delta_x', 'delta_p', 'integral_psi', 'rho_fluid',
            'A_osc', 'k_osc', 'omega_osc', 'x_pos', 'M_DM_factor', 'delta_rho_over_rho'
        ];

        if (validVars.includes(varName)) {
            return this[varName];
        } else {
            console.error(`Error: Unknown variable '${varName}'.`);
            return 0.0;
        }
    }

    M_SN_t(t) {
        return this.M_SN0 * Math.exp(-t / this.tau_SN);
    }

    compute_Ug() {
        const Ug1 = this.ug1_base;
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const corr_B = 1 - this.B / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ);
    }

    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    compute_g_NGC2525(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        const MSNt = this.M_SN_t(t);

        // Term 1: Base + Hz + B corrections
        const corr_H = 1 + this.Hz * t;
        const corr_B = 1 - this.B / this.B_crit;
        const term1 = this.ug1_base * corr_H * corr_B;

        // BH term
        const term_BH = this.g_BH;

        // Term 2: UQFF Ug with f_TRZ
        const term2 = this.compute_Ug();

        // Term 3: Lambda
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled EM with UA
        const cross_vB = this.gas_v * this.B; // Magnitude, assuming perpendicular
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        const term4 = (em_base * corr_UA) * this.scale_EM;

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

        // SN mass loss term (negative acceleration)
        const term_SN = -(this.G * MSNt) / (this.r * this.r);

        // Total g_NGC2525 (all terms summed)
        const total_g = term1 + term_BH + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_SN;

        // Add dynamic terms if enabled
        if (this.enableDynamicTerms) {
            for (const term of this.dynamicTerms) {
                total_g += term.compute(t, this.getAllParameters());
            }
        }

        return total_g;
    }

    getAllParameters() {
        const params = new Map();
        params.set('G', this.G);
        params.set('M', this.M);
        params.set('r', this.r);
        params.set('Hz', this.Hz);
        params.set('B', this.B);
        params.set('B_crit', this.B_crit);
        params.set('Lambda', this.Lambda);
        params.set('c_light', this.c_light);
        params.set('q_charge', this.q_charge);
        params.set('gas_v', this.gas_v);
        params.set('f_TRZ', this.f_TRZ);
        params.set('M_BH', this.M_BH);
        params.set('r_BH', this.r_BH);
        params.set('M_SN0', this.M_SN0);
        params.set('tau_SN', this.tau_SN);
        params.set('rho_vac_UA', this.rho_vac_UA);
        params.set('rho_vac_SCm', this.rho_vac_SCm);
        params.set('scale_EM', this.scale_EM);
        params.set('proton_mass', this.proton_mass);
        params.set('z_gal', this.z_gal);
        params.set('hbar', this.hbar);
        params.set('t_Hubble', this.t_Hubble);
        params.set('t_Hubble_gyr', this.t_Hubble_gyr);
        params.set('delta_x', this.delta_x);
        params.set('delta_p', this.delta_p);
        params.set('integral_psi', this.integral_psi);
        params.set('rho_fluid', this.rho_fluid);
        params.set('A_osc', this.A_osc);
        params.set('k_osc', this.k_osc);
        params.set('omega_osc', this.omega_osc);
        params.set('x_pos', this.x_pos);
        params.set('M_DM_factor', this.M_DM_factor);
        params.set('delta_rho_over_rho', this.delta_rho_over_rho);

        // Add dynamic parameters
        for (const [key, value] of this.dynamicParameters) {
            params.set(key, value);
        }

        return params;
    }

    registerDynamicTerm(term) {
        if (term instanceof PhysicsTerm) {
            this.dynamicTerms.push(term);
            if (this.enableLogging) {
                console.log(`Registered dynamic term: ${term.getName()}`);
            }
        }
    }

    setDynamicParameter(key, value) {
        this.dynamicParameters.set(key, value);
    }

    getDynamicParameter(key, defaultValue = 0.0) {
        return this.dynamicParameters.get(key) || defaultValue;
    }

    exportState(filename) {
        const state = {
            parameters: Object.fromEntries(this.getAllParameters()),
            dynamicTerms: this.dynamicTerms.map(term => ({
                name: term.getName(),
                description: term.getDescription(),
                parameters: Object.fromEntries(term.dynamicParameters)
            })),
            metadata: Object.fromEntries(this.metadata),
            enableDynamicTerms: this.enableDynamicTerms,
            enableLogging: this.enableLogging,
            learningRate: this.learningRate
        };

        // In Node.js, you would write to file here
        if (this.enableLogging) {
            console.log(`State exported to ${filename}`);
        }
        return state;
    }

    setEnableLogging(enable) {
        this.enableLogging = enable;
    }

    setLearningRate(rate) {
        this.learningRate = rate;
    }

    printParameters() {
        console.log("NGC 2525 Parameters:");
        console.log(`G: ${this.G}, M: ${this.M}, r: ${this.r}`);
        console.log(`Hz: ${this.Hz}, B: ${this.B}, B_crit: ${this.B_crit}`);
        console.log(`f_TRZ: ${this.f_TRZ}, M_BH: ${this.M_BH}, r_BH: ${this.r_BH}`);
        console.log(`M_SN0: ${this.M_SN0}, tau_SN: ${this.tau_SN}`);
        console.log(`rho_fluid: ${this.rho_fluid}, M_DM_factor: ${this.M_DM_factor}`);
        console.log(`A_osc: ${this.A_osc}, delta_rho_over_rho: ${this.delta_rho_over_rho}`);
        console.log(`ug1_base: ${this.ug1_base}, g_BH: ${this.g_BH}`);
    }

    exampleAt7Years() {
        const t_example = 7 * 3.156e7;
        return this.compute_g_NGC2525(t_example);
    }
}

module.exports = { GalaxyNGC2525 };