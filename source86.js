// source86.js - Source86 UQFF Module
// JavaScript implementation maintaining all dynamics from source86.cpp
// Models multiple astronomical systems with compressed and resonance-based UQFF models
// Supports systems: Magnetar SGR 1745-2900, Sagittarius A*, Tapestry of Blazing Starbirth, Westerlund 2, Pillars of Creation, Rings of Relativity, Students Guide to the Universe
// No SM gravity/magnetics - pure UQFF frequency domain calculations

class PhysicsTerm {
    constructor() {
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;
    }

    compute(t, params) {
        throw new Error('compute method must be implemented by subclass');
    }

    getName() {
        throw new Error('getName method must be implemented by subclass');
    }

    getDescription() {
        throw new Error('getDescription method must be implemented by subclass');
    }

    validate(params) {
        return true;
    }

    registerDynamicTerm(term) {
        if (this.enableDynamicTerms) {
            this.dynamicTerms.push(term);
            if (this.enableLogging) {
                console.log(`Registered dynamic term: ${term.getName()}`);
            }
        }
    }

    setDynamicParameter(key, value) {
        this.dynamicParameters.set(key, value);
    }

    getDynamicParameter(key) {
        return this.dynamicParameters.get(key);
    }

    setEnableLogging(enable) {
        this.enableLogging = enable;
    }

    setLearningRate(rate) {
        this.learningRate = rate;
    }

    exportState(filename) {
        const state = {
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            metadata: Object.fromEntries(this.metadata),
            enableDynamicTerms: this.enableDynamicTerms,
            enableLogging: this.enableLogging,
            learningRate: this.learningRate
        };
        // In Node.js environment, this would write to file
        if (typeof require !== 'undefined') {
            const fs = require('fs');
            fs.writeFileSync(filename, JSON.stringify(state, null, 2));
        }
        return state;
    }
}

class DynamicVacuumTerm extends PhysicsTerm {
    constructor(amplitude = 1e-10, frequency = 1e-15) {
        super();
        this.amplitude = amplitude;
        this.frequency = frequency;
        this.metadata.set('type', 'vacuum_energy');
        this.metadata.set('description', 'Time-varying vacuum energy density');
    }

    compute(t, params) {
        const rho_vac = params.get('rho_vac_UA') || 7.09e-36;
        return this.amplitude * rho_vac * Math.sin(this.frequency * t);
    }

    getName() {
        return 'DynamicVacuum';
    }

    getDescription() {
        return 'Time-varying vacuum energy with sinusoidal modulation';
    }
}

class QuantumCouplingTerm extends PhysicsTerm {
    constructor(coupling_strength = 1e-40) {
        super();
        this.coupling_strength = coupling_strength;
        this.metadata.set('type', 'quantum_coupling');
        this.metadata.set('description', 'Non-local quantum coupling effects');
    }

    compute(t, params) {
        const hbar = params.get('hbar') || 1.0546e-34;
        const M = params.get('M') || 1.989e30;
        const r = params.get('r') || 1e4;
        return this.coupling_strength * (hbar * hbar) / (M * r * r) * Math.cos(t / 1e6);
    }

    getName() {
        return 'QuantumCoupling';
    }

    getDescription() {
        return 'Non-local quantum effects with oscillatory coupling';
    }
}

// System type enumeration
const SystemType = {
    MAGNETAR_SGR_1745_2900: 0,
    SAGITTARIUS_A: 1,
    TAPESTRY_BLAZING_STARBIRTH: 2,
    WESTERLUND_2: 3,
    PILLARS_CREATION: 4,
    RINGS_RELATIVITY: 5,
    STUDENTS_GUIDE_UNIVERSE: 6
};

class Source86UQFFModule {
    constructor(systemType = SystemType.MAGNETAR_SGR_1745_2900) {
        // Initialize dynamic framework
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;
        this.metadata.set('enhanced', 'true');
        this.metadata.set('version', '2.0-Enhanced');
        this.metadata.set('system', 'Multi-System MUGE');

        // Core variables (maintaining all dynamics from C++)
        this.variables = new Map();

        // Universal constants
        this.variables.set('G', 6.6743e-11);                    // m^3 kg^-1 s^-2
        this.variables.set('c', 3e8);                           // m/s
        this.variables.set('hbar', 1.0546e-34);                 // J s
        this.variables.set('Lambda', 1.1e-52);                  // m^-2
        this.variables.set('q', 1.602e-19);                     // C
        this.variables.set('pi', Math.PI);
        this.variables.set('t_Hubble', 4.35e17);                // s
        this.variables.set('H0', 2.269e-18);                    // s^-1 (70 km/s/Mpc)
        this.variables.set('Omega_m', 0.3);
        this.variables.set('Omega_Lambda', 0.7);
        this.variables.set('year_to_s', 3.156e7);
        this.variables.set('M_sun', 1.989e30);                  // kg

        // Quantum defaults
        this.variables.set('Delta_x', 1e-10);
        this.variables.set('Delta_p', this.variables.get('hbar') / this.variables.get('Delta_x'));
        this.variables.set('integral_psi', 2.176e-18);          // J, normalized

        // Resonance defaults
        this.variables.set('Evac_neb', 7.09e-36);               // J/m^3
        this.variables.set('Evac_ISM', 7.09e-37);               // J/m^3
        this.variables.set('Delta_Evac', 6.381e-36);            // J/m^3
        this.variables.set('v_exp', 1e3);                       // m/s
        this.variables.set('f_THz', 1e12);                      // Hz, placeholder
        this.variables.set('f_DPM', 1e9);                       // Hz
        this.variables.set('FDPM', 6.284e29);                   // A m^2
        this.variables.set('F_super', 6.287e-19);               // dimensionless
        this.variables.set('UA_SCm', 10.0);                     // scaling
        this.variables.set('omega_i', 1e-8);                    // rad/s
        this.variables.set('k4', 1.0);
        this.variables.set('f_react', 1e10);                    // Hz
        this.variables.set('E_react', 1e-20);                   // J
        this.variables.set('f_quantum', 1.445e-17);             // Hz
        this.variables.set('f_Aether', 1.576e-35);              // Hz
        this.variables.set('f_fluid', 1.269e-14);               // Hz
        this.variables.set('f_osc', 4.57e14);                   // Hz
        this.variables.set('f_exp', 1e-18);                     // Hz
        this.variables.set('f_TRZ', 0.1);                       // dimensionless

        // Fluid/DM defaults
        this.variables.set('rho_fluid', 1e-20);                 // kg/m^3
        this.variables.set('V', 1e3);                           // m^3
        this.variables.set('g_local', 9.8);                     // m/s^2
        this.variables.set('DM_fraction', 0.85);
        this.variables.set('delta_rho_over_rho', 1e-5);

        // Ug defaults (negligible except Ug3' where applicable)
        this.variables.set('Ug1', 0.0);
        this.variables.set('Ug2', 0.0);
        this.variables.set('Ug3_prime', 0.0);
        this.variables.set('Ug4', 0.0);

        // Initialize with default dynamic terms
        this.registerDynamicTerm(new DynamicVacuumTerm());
        this.registerDynamicTerm(new QuantumCouplingTerm());

        // Set initial system
        this.setSystem(systemType);
    }

    // Set system and update params
    setSystem(systemType) {
        this.currentSystem = systemType;
        switch (systemType) {
            case SystemType.MAGNETAR_SGR_1745_2900:
                this.variables.set('M', 1.5 * this.variables.get('M_sun'));
                this.variables.set('r', 1e4);                       // m
                this.variables.set('z', 0.0009);
                this.variables.set('B', 1e10);                      // T
                this.variables.set('B_crit', 1e11);                 // T
                this.variables.set('r_BH', 2.84e15);                // m to Sgr A*
                this.variables.set('M_BH', 4.1e6 * this.variables.get('M_sun'));
                this.variables.set('t', 3.799e10);                  // s
                this.variables.set('rho_fluid', 1e-15);
                this.variables.set('V', 4.189e12);
                this.variables.set('g_local', 10.0);
                this.variables.set('M_DM', 0.0);
                this.variables.set('M_visible', this.variables.get('M'));
                this.variables.set('Ug3_prime', (this.variables.get('G') * this.variables.get('M_BH')) / (this.variables.get('r_BH') * this.variables.get('r_BH')));
                this.variables.set('F_env', 0.0);                   // Mmag + D(t) negligible
                this.variables.set('v_wind', 0.0);                  // No wind
                break;
            case SystemType.SAGITTARIUS_A:
                this.variables.set('M', 4.1e6 * this.variables.get('M_sun'));
                this.variables.set('r', 1.18e10);                   // m (event horizon approx)
                this.variables.set('z', 0.00034);
                this.variables.set('B', 1e-5);
                this.variables.set('B_crit', 1e11);
                this.variables.set('t', 1e6 * this.variables.get('year_to_s'));
                this.variables.set('rho_fluid', 1e-20);
                this.variables.set('V', 1e3);
                this.variables.set('g_local', 1e-6);
                this.variables.set('M_DM', 0.85 * this.variables.get('M'));
                this.variables.set('M_visible', 0.15 * this.variables.get('M'));
                this.variables.set('spin_adjust', Math.sin(30.0 * this.variables.get('pi') / 180.0));  // sin(30)
                this.variables.set('dOmega_dt', 1e-3);              // rad/s, placeholder for GW
                this.variables.set('F_env', 0.0);
                this.variables.set('v_wind', 8e3);
                break;
            case SystemType.TAPESTRY_BLAZING_STARBIRTH:
                this.variables.set('M', 2000 * this.variables.get('M_sun'));
                this.variables.set('r', 1.18e17);
                this.variables.set('z', 0.00034);
                this.variables.set('B', 1e-5);
                this.variables.set('B_crit', 1e11);
                this.variables.set('t', 1e6 * this.variables.get('year_to_s'));
                this.variables.set('rho_fluid', 1e-20);
                this.variables.set('V', 1e3);
                this.variables.set('g_local', 1e-12);
                this.variables.set('M_DM', 0.85 * this.variables.get('M'));
                this.variables.set('M_visible', 0.15 * this.variables.get('M'));
                this.variables.set('rho', this.variables.get('rho_fluid'));
                this.variables.set('F_env', 0.0);
                this.variables.set('v_wind', 8e3);
                break;
            case SystemType.WESTERLUND_2:
                this.variables.set('M', 3000 * this.variables.get('M_sun'));
                this.variables.set('r', 2e17);
                this.variables.set('z', 0.001);
                this.variables.set('B', 1e-5);
                this.variables.set('B_crit', 1e11);
                this.variables.set('t', 2e6 * this.variables.get('year_to_s'));
                this.variables.set('rho_fluid', 1e-20);
                this.variables.set('V', 1e3);
                this.variables.set('g_local', 1e-12);
                this.variables.set('M_DM', 0.85 * this.variables.get('M'));
                this.variables.set('M_visible', 0.15 * this.variables.get('M'));
                this.variables.set('rho', this.variables.get('rho_fluid'));
                this.variables.set('F_env', 0.0);
                this.variables.set('v_wind', 1e4);
                break;
            case SystemType.PILLARS_CREATION:
                this.variables.set('M', 800 * this.variables.get('M_sun'));
                this.variables.set('r', 1e17);
                this.variables.set('z', 0.002);
                this.variables.set('B', 1e-6);
                this.variables.set('B_crit', 1e11);
                this.variables.set('t', 1e6 * this.variables.get('year_to_s'));
                this.variables.set('rho_fluid', 1e-19);
                this.variables.set('V', 1e4);
                this.variables.set('g_local', 1e-11);
                this.variables.set('M_DM', 0.85 * this.variables.get('M'));
                this.variables.set('M_visible', 0.15 * this.variables.get('M'));
                this.variables.set('E_t', 0.1);                     // Erosion term
                this.variables.set('rho', this.variables.get('rho_fluid'));
                this.variables.set('F_env', 0.0);
                this.variables.set('v_wind', 8e3);
                break;
            case SystemType.RINGS_RELATIVITY:
                this.variables.set('M', 1e6 * this.variables.get('M_sun'));
                this.variables.set('r', 1e16);
                this.variables.set('z', 0.01);
                this.variables.set('B', 1e-4);
                this.variables.set('B_crit', 1e11);
                this.variables.set('t', 1e7 * this.variables.get('year_to_s'));
                this.variables.set('rho_fluid', 1e-21);
                this.variables.set('V', 1e5);
                this.variables.set('g_local', 1e-10);
                this.variables.set('M_DM', 0.85 * this.variables.get('M'));
                this.variables.set('M_visible', 0.15 * this.variables.get('M'));
                this.variables.set('L_t', 0.05);                    // Lensing term
                this.variables.set('F_env', 0.0);
                this.variables.set('v_wind', 5e3);
                break;
            case SystemType.STUDENTS_GUIDE_UNIVERSE:
                this.variables.set('M', 1 * this.variables.get('M_sun'));
                this.variables.set('r', 1e11);                      // AU scale
                this.variables.set('z', 0.0);
                this.variables.set('B', 1e-5);
                this.variables.set('B_crit', 1e11);
                this.variables.set('t', 1e9 * this.variables.get('year_to_s'));
                this.variables.set('rho_fluid', 1e-25);
                this.variables.set('V', 1e12);
                this.variables.set('g_local', 1e-11);
                this.variables.set('M_DM', 0.27 * this.variables.get('M'));
                this.variables.set('M_visible', 0.73 * this.variables.get('M'));
                this.variables.set('F_env', 0.0);
                this.variables.set('v_wind', 0.0);
                break;
        }

        // Update dependent variables
        this.variables.set('Delta_p', this.variables.get('hbar') / this.variables.get('Delta_x'));
        this.variables.set('M_visible', (1.0 - this.variables.get('DM_fraction')) * this.variables.get('M'));
        this.variables.set('M_DM', this.variables.get('DM_fraction') * this.variables.get('M'));

        // System-specific updates
        if (this.currentSystem === SystemType.MAGNETAR_SGR_1745_2900 || this.currentSystem === SystemType.SAGITTARIUS_A) {
            if (this.variables.has('M_BH') && this.variables.has('r_BH')) {
                this.variables.set('Ug3_prime', (this.variables.get('G') * this.variables.get('M_BH')) / (this.variables.get('r_BH') * this.variables.get('r_BH')));
            }
        }
    }

    // Dynamic variable operations
    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
        } else {
            if (this.enableLogging) {
                console.log(`Variable '${name}' not found. Adding.`);
            }
            this.variables.set(name, value);
        }

        // Update dependent variables
        if (name === 'Delta_x') {
            this.variables.set('Delta_p', this.variables.get('hbar') / value);
        } else if (name === 'M') {
            this.variables.set('M_visible', (1.0 - this.variables.get('DM_fraction')) * value);
            this.variables.set('M_DM', this.variables.get('DM_fraction') * value);
        } else if (name === 'DM_fraction') {
            this.variables.set('M_visible', (1.0 - value) * this.variables.get('M'));
            this.variables.set('M_DM', value * this.variables.get('M'));
        }

        // System-specific updates
        if ((this.currentSystem === SystemType.MAGNETAR_SGR_1745_2900 || this.currentSystem === SystemType.SAGITTARIUS_A) &&
            this.variables.has('M_BH') && this.variables.has('r_BH')) {
            this.variables.set('Ug3_prime', (this.variables.get('G') * this.variables.get('M_BH')) / (this.variables.get('r_BH') * this.variables.get('r_BH')));
        }
    }

    addToVariable(name, delta) {
        const current = this.variables.get(name) || 0;
        this.variables.set(name, current + delta);
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    // Compute H(t,z)
    computeH(t, z) {
        const Hz = this.variables.get('H0') * Math.sqrt(
            this.variables.get('Omega_m') * Math.pow(1 + z, 3) +
            this.variables.get('Omega_Lambda')
        );
        return Hz * t;
    }

    // Ug sum (Ug3' for external, others 0)
    computeUgSum() {
        return this.variables.get('Ug1') + this.variables.get('Ug2') +
               this.variables.get('Ug3_prime') + this.variables.get('Ug4');
    }

    // Lambda term
    computeLambdaTerm() {
        return (this.variables.get('Lambda') * this.variables.get('c') * this.variables.get('c')) / 3.0;
    }

    // Quantum term
    computeQuantumTerm() {
        const unc = Math.sqrt(this.variables.get('Delta_x') * this.variables.get('Delta_p'));
        const integral_val = this.variables.get('integral_psi');
        return (this.variables.get('hbar') / unc) * integral_val * (2 * this.variables.get('pi') / this.variables.get('t_Hubble'));
    }

    // Fluid term
    computeFluidTerm(g_base) {
        return this.variables.get('rho_fluid') * this.variables.get('V') * g_base;
    }

    // DM term
    computeDMTerm() {
        const pert = this.variables.get('delta_rho_over_rho');
        const curv = 3 * this.variables.get('G') * this.variables.get('M') / (this.variables.get('r') * this.variables.get('r') * this.variables.get('r'));
        return (this.variables.get('M_visible') + this.variables.get('M_DM')) * (pert + curv);
    }

    // Resonant term (cos + Re[exp])
    computeResonantTerm(t) {
        const A = this.variables.get('A') || 1e-10;  // Assume added if needed, default 1e-10
        const k = this.variables.get('k') || 1e20;  // 1e20
        const omega = this.variables.get('omega') || 1e15;  // 1e15
        const x = 0.0;
        const cos_term = 2 * A * Math.cos(k * x) * Math.cos(omega * t);
        // Complex exponential: A * exp(i(k*x - ω*t))
        const real_exp = A * Math.cos(k * x - omega * t);
        const exp_factor = (2 * this.variables.get('pi') / 13.8);
        return cos_term + exp_factor * real_exp;
    }

    // EM term q (v x B) magnitude
    computeEMTerm() {
        const v = this.variables.get('v_wind');
        const B = this.variables.get('B');
        const scale_macro = this.variables.get('scale_macro') || 1e-12;
        return (this.variables.get('q') * v * B) / 1.673e-27 * scale_macro;  // Scaled
    }

    // System-specific term (e.g., wind, erosion, lensing)
    computeSystemSpecificTerm(t) {
        let term = 0.0;
        switch (this.currentSystem) {
            case SystemType.SAGITTARIUS_A:
                term += (this.variables.get('G') * this.variables.get('M') * this.variables.get('M') /
                        (this.variables.get('c') * this.variables.get('c') * this.variables.get('c') * this.variables.get('c') * this.variables.get('r'))) *
                        Math.pow(this.variables.get('dOmega_dt'), 2);
                term *= this.variables.get('spin_adjust');
                break;
            case SystemType.TAPESTRY_BLAZING_STARBIRTH:
            case SystemType.WESTERLUND_2:
                term += this.variables.get('rho') * Math.pow(this.variables.get('v_wind'), 2);
                break;
            case SystemType.PILLARS_CREATION:
                term += this.variables.get('rho') * Math.pow(this.variables.get('v_wind'), 2) * (1 - this.variables.get('E_t'));
                break;
            case SystemType.RINGS_RELATIVITY:
                term += this.variables.get('rho_fluid') * this.variables.get('V') * this.variables.get('g_local') * (1 + this.variables.get('L_t'));
                break;
            case SystemType.STUDENTS_GUIDE_UNIVERSE:
                term = 0.0;  // Simplified
                break;
            default:
                term += this.variables.get('rho_fluid') * Math.pow(this.variables.get('v_wind'), 2);  // Default wind
        }
        return term;
    }

    // Compressed MUGE
    computeG_compressed(t) {
        this.variables.set('t', t);
        const Hz_t = this.computeH(t, this.variables.get('z'));
        const expansion = 1.0 + Hz_t;
        const sc_correction = 1.0 - (this.variables.get('B') / this.variables.get('B_crit'));
        const env_factor = 1.0 + this.variables.get('F_env');
        const g_base = (this.variables.get('G') * this.variables.get('M') / (this.variables.get('r') * this.variables.get('r'))) *
                       expansion * sc_correction * env_factor;

        const ug_sum = this.computeUgSum();
        const lambda_term = this.computeLambdaTerm();
        const quantum_term = this.computeQuantumTerm();
        const em_term = this.computeEMTerm();
        const fluid_term = this.computeFluidTerm(g_base);
        const resonant_term = this.computeResonantTerm(t);
        const dm_term = this.computeDMTerm();
        const sys_term = this.computeSystemSpecificTerm(t);

        // Add dynamic terms
        let total = g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + sys_term;
        for (const term of this.dynamicTerms) {
            total += term.compute(t, this.variables);
        }

        return total;
    }

    // Resonance helpers
    computeADPM() {
        return this.variables.get('c') * this.variables.get('V') * this.variables.get('FDPM') *
               this.variables.get('f_DPM') * this.variables.get('Evac_neb');
    }

    computeATHz() {
        return (this.variables.get('Evac_ISM') / this.variables.get('c')) * this.variables.get('f_THz') *
               this.variables.get('Evac_neb') * this.variables.get('v_exp') * this.computeADPM();
    }

    computeAvacDiff() {
        return (this.variables.get('Evac_neb') / (this.variables.get('c') * this.variables.get('c'))) *
               this.variables.get('Delta_Evac') * Math.pow(this.variables.get('v_exp'), 2) * this.computeADPM();
    }

    computeASuperFreq() {
        return (this.variables.get('Evac_neb') / this.variables.get('c')) * this.variables.get('F_super') *
               this.variables.get('f_THz') * this.computeADPM();
    }

    computeAAetherRes() {
        return this.variables.get('UA_SCm') * this.variables.get('omega_i') * this.variables.get('f_THz') *
               this.computeADPM() * (1 + this.variables.get('f_TRZ'));
    }

    computeUg4i() {
        return this.variables.get('k4') * this.variables.get('E_react') * this.variables.get('f_react') *
               this.computeADPM() / (this.variables.get('Evac_neb') * this.variables.get('c'));
    }

    computeAQuantumFreq() {
        return (this.variables.get('Evac_ISM') / this.variables.get('c')) * this.variables.get('f_quantum') *
               this.variables.get('Evac_neb') * this.computeADPM();
    }

    computeAAetherFreq() {
        return (this.variables.get('Evac_ISM') / this.variables.get('c')) * this.variables.get('f_Aether') *
               this.variables.get('Evac_neb') * this.computeADPM();
    }

    computeAFluidFreq() {
        return (this.variables.get('Evac_ISM') / this.variables.get('c')) * this.variables.get('f_fluid') *
               this.variables.get('Evac_neb') * this.variables.get('V');
    }

    computeOscTerm(t) {
        const A = 1e-10;  // Default
        const omega = this.variables.get('f_osc') * 2 * this.variables.get('pi');
        return 2 * A * Math.cos(omega * t);  // Simplified osc
    }

    computeAExpFreq() {
        return (this.variables.get('Evac_ISM') / this.variables.get('c')) * this.variables.get('f_exp') *
               this.variables.get('Evac_neb') * this.computeADPM();
    }

    computeFTRZ() {
        return this.variables.get('f_TRZ');
    }

    // Resonance MUGE
    computeG_resonance(t) {
        const aDPM = this.computeADPM();
        const aTHz = this.computeATHz();
        const aVacDiff = this.computeAvacDiff();
        const aSuperFreq = this.computeASuperFreq();
        const aAetherRes = this.computeAAetherRes();
        const ug4i = this.computeUg4i();
        const aQuantumFreq = this.computeAQuantumFreq();
        const aAetherFreq = this.computeAAetherFreq();
        const aFluidFreq = this.computeAFluidFreq();
        const oscTerm = this.computeOscTerm(t);
        const aExpFreq = this.computeAExpFreq();
        const fTRZ = this.computeFTRZ();

        return aDPM + aTHz + aVacDiff + aSuperFreq + aAetherRes + ug4i + aQuantumFreq +
               aAetherFreq + aFluidFreq + oscTerm + aExpFreq + fTRZ;
    }

    // Dynamic term management
    registerDynamicTerm(term) {
        if (this.enableDynamicTerms) {
            this.dynamicTerms.push(term);
            if (this.enableLogging) {
                console.log(`Registered dynamic term: ${term.getName()}`);
            }
        }
    }

    setDynamicParameter(key, value) {
        this.dynamicParameters.set(key, value);
    }

    getDynamicParameter(key) {
        return this.dynamicParameters.get(key);
    }

    setEnableLogging(enable) {
        this.enableLogging = enable;
        for (const term of this.dynamicTerms) {
            term.setEnableLogging(enable);
        }
    }

    setLearningRate(rate) {
        this.learningRate = rate;
        for (const term of this.dynamicTerms) {
            term.setLearningRate(rate);
        }
    }

    exportState(filename) {
        const state = {
            variables: Object.fromEntries(this.variables),
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            metadata: Object.fromEntries(this.metadata),
            enableDynamicTerms: this.enableDynamicTerms,
            enableLogging: this.enableLogging,
            learningRate: this.learningRate,
            currentSystem: this.currentSystem,
            dynamicTerms: this.dynamicTerms.map(term => ({
                name: term.getName(),
                description: term.getDescription()
            }))
        };

        // In Node.js environment, write to file
        if (typeof require !== 'undefined') {
            const fs = require('fs');
            fs.writeFileSync(filename, JSON.stringify(state, null, 2));
        }

        return state;
    }

    // Equation descriptions
    getEquationText_compressed() {
        let sys_name;
        switch (this.currentSystem) {
            case SystemType.MAGNETAR_SGR_1745_2900: sys_name = "Magnetar SGR 1745-2900"; break;
            case SystemType.SAGITTARIUS_A: sys_name = "Sagittarius A*"; break;
            case SystemType.TAPESTRY_BLAZING_STARBIRTH: sys_name = "Tapestry of Blazing Starbirth"; break;
            case SystemType.WESTERLUND_2: sys_name = "Westerlund 2"; break;
            case SystemType.PILLARS_CREATION: sys_name = "Pillars of Creation"; break;
            case SystemType.RINGS_RELATIVITY: sys_name = "Rings of Relativity"; break;
            case SystemType.STUDENTS_GUIDE_UNIVERSE: sys_name = "Students Guide to the Universe"; break;
            default: sys_name = "Generic";
        }
        return `Compressed MUGE for ${sys_name}:
g(r,t) = (G M(t)/r^2) (1 + H(t,z)) (1 - B/B_crit) (1 + F_env(t)) + (Ug1 + Ug2 + Ug3' + Ug4) + (Lambda c^2 / 3) + (hbar / sqrt(Delta_x Delta_p)) ∫(ψ H ψ dV) (2π / t_Hubble) + q (v × B) + ρ_fluid V g + 2 A cos(k x) cos(ω t) + (2π/13.8) A Re[exp(i (k x - ω t))] + (M_vis + M_DM) (δρ/ρ + 3 G M / r^3) + SysTerm
SysTerm: e.g., for Magnetar: G M_BH / r_BH^2; for Sgr A*: (G M² / c⁴ r) (dΩ/dt)² sin(30); for Starbirth: ρ v_wind²
Variables: As in map; Approximations: Ug1=Ug2=Ug4=0, integral normalized=1.0 scaled.`;
    }

    getEquationText_resonance() {
        let sys_name;
        switch (this.currentSystem) {
            case SystemType.MAGNETAR_SGR_1745_2900: sys_name = "Magnetar SGR 1745-2900"; break;
            case SystemType.SAGITTARIUS_A: sys_name = "Sagittarius A*"; break;
            case SystemType.TAPESTRY_BLAZING_STARBIRTH: sys_name = "Tapestry of Blazing Starbirth"; break;
            case SystemType.WESTERLUND_2: sys_name = "Westerlund 2"; break;
            case SystemType.PILLARS_CREATION: sys_name = "Pillars of Creation"; break;
            case SystemType.RINGS_RELATIVITY: sys_name = "Rings of Relativity"; break;
            case SystemType.STUDENTS_GUIDE_UNIVERSE: sys_name = "Students Guide to the Universe"; break;
            default: sys_name = "Generic";
        }
        return `Resonance MUGE for ${sys_name}:
g(r,t) = a_DPM + a_THz + a_vac_diff + a_super_freq + a_aether_res + Ug4i + a_quantum_freq + a_Aether_freq + a_fluid_freq + Osc_term + a_exp_freq + f_TRZ
Where a_DPM = c V_sys F_DPM f_DPM E_vac,neb; a_THz = (E_vac,ISM / c) f_THz E_vac,neb v_exp a_DPM; etc.
Variables: Resonance freqs tuned (f_THz=1e12 Hz, etc.); Osc_term ≈ 2 A cos(ω t); f_TRZ=0.1.
Integration: Sum yields effective g ~1e-11 m/s² for nebulae, dominated by fluid/resonant terms.`;
    }

    // Print all variables
    printVariables() {
        console.log('Source86 Variables for system', this.currentSystem, ':');
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }
}

// Export the module
module.exports = { Source86UQFFModule, SystemType };