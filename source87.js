// source87.js - Source87 UQFF Module
// JavaScript implementation maintaining all dynamics from source87.cpp
// Resonance-based superconductive MUGE (UQFF) for multiple astronomical systems
// Uses frequency-driven dynamics via plasmotic vacuum energy and resonances
// Excludes SM gravity/magnetics - pure UQFF frequency domain calculations
// Supports 12 systems: Magnetar SGR 1745-2900, Sagittarius A*, Tapestry of Blazing Starbirth, Westerlund 2, Pillars of Creation, Rings of Relativity, Students Guide to the Universe, NGC 2525, NGC 3603, Bubble Nebula, Antennae Galaxies, Horsehead Nebula

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
        return 'Non-local quantum effects with time-dependent coupling';
    }
}

// System Type Enum
const SystemType = {
    MAGNETAR_SGR_1745_2900: 0,
    SAGITTARIUS_A: 1,
    TAPESTRY_BLAZING_STARBIRTH: 2,
    WESTERLUND_2: 3,
    PILLARS_CREATION: 4,
    RINGS_RELATIVITY: 5,
    STUDENTS_GUIDE_UNIVERSE: 6,
    NGC_2525: 7,
    NGC_3603: 8,
    BUBBLE_NEBULA: 9,
    ANTENNAE_GALAXIES: 10,
    HORSEHEAD_NEBULA: 11
};

class Source87UQFFModule {
    constructor(systemType = SystemType.MAGNETAR_SGR_1745_2900) {
        this.variables = new Map();
        this.currentSystem = systemType;

        // Self-expanding framework members
        this.dynamicParameters = new Map();
        this.dynamicTerms = [];
        this.metadata = new Map();
        this.enableDynamicTerms = true;
        this.enableLogging = false;
        this.learningRate = 0.001;

        this.initializeConstants();
        this.setSystem(systemType);
    }

    initializeConstants() {
        // Universal constants
        this.variables.set('c', 3e8);                           // m/s
        this.variables.set('pi', Math.PI);
        this.variables.set('H0', 2.269e-18);                   // s^-1
        this.variables.set('Omega_m', 0.3);
        this.variables.set('Omega_Lambda', 0.7);
        this.variables.set('G', 6.6743e-11);                   // For f_fluid tuning
        this.variables.set('M_sun', 1.989e30);                 // kg
        this.variables.set('year_to_s', 3.156e7);

        // Resonance parameters
        this.variables.set('Evac_neb', 7.09e-36);              // J/m^3
        this.variables.set('Evac_ISM', 7.09e-37);              // J/m^3
        this.variables.set('Delta_Evac', 6.381e-36);           // J/m^3
        this.variables.set('v_exp', 1e3);                      // m/s default
        this.variables.set('f_DPM', 1e12);                     // Hz
        this.variables.set('f_THz', 1e12);                     // Hz
        this.variables.set('f_quantum', 1.445e-17);            // Hz
        this.variables.set('f_Aether', 1.576e-35);             // Hz
        this.variables.set('f_fluid', 1e-14);                  // Hz default
        this.variables.set('f_react', 1e10);                   // Hz
        this.variables.set('f_osc', 4.57e14);                  // Hz
        this.variables.set('F_super', 6.287e-19);              // dimensionless
        this.variables.set('UA_SCm', 10.0);                    // scaling
        this.variables.set('omega_i', 1e-8);                   // rad/s
        this.variables.set('k4', 1.0);                         // dimensionless
        this.variables.set('E_react_base', 1e46);              // J
        this.variables.set('decay_rate', 5e-4);                // s^-1
        this.variables.set('f_TRZ', 0.1);                      // dimensionless

        // Vortices defaults
        this.variables.set('I', 1e21);                         // A
        this.variables.set('A_vort', 1e8);                     // m^2 default
        this.variables.set('omega1', 1e-3);                    // rad/s
        this.variables.set('omega2', -1e-3);                   // rad/s
    }

    setSystem(systemType) {
        this.currentSystem = systemType;
        const M_sun = this.variables.get('M_sun');

        switch (systemType) {
            case SystemType.MAGNETAR_SGR_1745_2900:
                this.variables.set('M', 1.5 * M_sun);
                this.variables.set('r', 1e4);
                this.variables.set('z', 0.0009);
                this.variables.set('t', 3.799e10);
                this.variables.set('I', 1e21);
                this.variables.set('A_vort', this.variables.get('pi') * this.variables.get('r') * this.variables.get('r'));
                this.variables.set('omega1', 1e-3);
                this.variables.set('omega2', -1e-3);
                this.variables.set('v_exp', 1e3);
                this.variables.set('Vsys', 4.189e12);
                this.variables.set('f_fluid', 1.269e-14);
                break;
            case SystemType.SAGITTARIUS_A:
                this.variables.set('M', 4.1e6 * M_sun);
                this.variables.set('r', 1.18e10);
                this.variables.set('z', 0.00034);
                this.variables.set('t', 1e6 * this.variables.get('year_to_s'));
                this.variables.set('I', 1e22);
                this.variables.set('A_vort', this.variables.get('pi') * this.variables.get('r') * this.variables.get('r'));
                this.variables.set('omega1', 1e-4);
                this.variables.set('omega2', -1e-4);
                this.variables.set('v_exp', 5e3);
                this.variables.set('Vsys', 4.0 / 3.0 * this.variables.get('pi') * Math.pow(this.variables.get('r'), 3));
                this.variables.set('f_fluid', 1e-12);
                break;
            case SystemType.TAPESTRY_BLAZING_STARBIRTH:
                this.variables.set('M', 2000 * M_sun);
                this.variables.set('r', 1.18e17);
                this.variables.set('z', 0.00034);
                this.variables.set('t', 1e6 * this.variables.get('year_to_s'));
                this.variables.set('I', 1e23);
                this.variables.set('A_vort', this.variables.get('pi') * this.variables.get('r') * this.variables.get('r'));
                this.variables.set('omega1', 1e-5);
                this.variables.set('omega2', -1e-5);
                this.variables.set('v_exp', 1e4);
                this.variables.set('Vsys', 4.0 / 3.0 * this.variables.get('pi') * Math.pow(this.variables.get('r'), 3));
                this.variables.set('f_fluid', 1e-14);
                break;
            case SystemType.WESTERLUND_2:
                this.variables.set('M', 3000 * M_sun);
                this.variables.set('r', 2e17);
                this.variables.set('z', 0.001);
                this.variables.set('t', 2e6 * this.variables.get('year_to_s'));
                this.variables.set('I', 1e23);
                this.variables.set('A_vort', this.variables.get('pi') * this.variables.get('r') * this.variables.get('r'));
                this.variables.set('omega1', 1e-5);
                this.variables.set('omega2', -1e-5);
                this.variables.set('v_exp', 1e4);
                this.variables.set('Vsys', 4.0 / 3.0 * this.variables.get('pi') * Math.pow(this.variables.get('r'), 3));
                this.variables.set('f_fluid', 1e-14);
                break;
            case SystemType.PILLARS_CREATION:
                this.variables.set('M', 800 * M_sun);
                this.variables.set('r', 1e17);
                this.variables.set('z', 0.002);
                this.variables.set('t', 1e6 * this.variables.get('year_to_s'));
                this.variables.set('I', 1e22);
                this.variables.set('A_vort', this.variables.get('pi') * this.variables.get('r') * this.variables.get('r'));
                this.variables.set('omega1', 1e-5);
                this.variables.set('omega2', -1e-5);
                this.variables.set('v_exp', 8e3);
                this.variables.set('Vsys', 4.0 / 3.0 * this.variables.get('pi') * Math.pow(this.variables.get('r'), 3));
                this.variables.set('f_fluid', 1e-13);
                break;
            case SystemType.RINGS_RELATIVITY:
                this.variables.set('M', 1e6 * M_sun);
                this.variables.set('r', 1e16);
                this.variables.set('z', 0.01);
                this.variables.set('t', 1e7 * this.variables.get('year_to_s'));
                this.variables.set('I', 1e23);
                this.variables.set('A_vort', this.variables.get('pi') * this.variables.get('r') * this.variables.get('r'));
                this.variables.set('omega1', 1e-6);
                this.variables.set('omega2', -1e-6);
                this.variables.set('v_exp', 5e3);
                this.variables.set('Vsys', 4.0 / 3.0 * this.variables.get('pi') * Math.pow(this.variables.get('r'), 3));
                this.variables.set('f_fluid', 1e-12);
                break;
            case SystemType.STUDENTS_GUIDE_UNIVERSE:
                this.variables.set('M', 1 * M_sun);
                this.variables.set('r', 1e11);
                this.variables.set('z', 0.0);
                this.variables.set('t', 1e9 * this.variables.get('year_to_s'));
                this.variables.set('I', 1e20);
                this.variables.set('A_vort', this.variables.get('pi') * this.variables.get('r') * this.variables.get('r'));
                this.variables.set('omega1', 1e-2);
                this.variables.set('omega2', -1e-2);
                this.variables.set('v_exp', 1e2);
                this.variables.set('Vsys', 4.0 / 3.0 * this.variables.get('pi') * Math.pow(this.variables.get('r'), 3));
                this.variables.set('f_fluid', 1e-10);
                break;
            case SystemType.NGC_2525:
                this.variables.set('M', 1e10 * M_sun);
                this.variables.set('r', 1e20);
                this.variables.set('z', 0.001);
                this.variables.set('t', 1e9 * this.variables.get('year_to_s'));
                this.variables.set('I', 1e24);
                this.variables.set('A_vort', this.variables.get('pi') * this.variables.get('r') * this.variables.get('r'));
                this.variables.set('omega1', 1e-6);
                this.variables.set('omega2', -1e-6);
                this.variables.set('v_exp', 1e5);
                this.variables.set('Vsys', 1.543e64);
                this.variables.set('f_fluid', 8.457e-4);
                break;
            case SystemType.NGC_3603:
                this.variables.set('M', 2000 * M_sun);
                this.variables.set('r', 1.18e17);
                this.variables.set('z', 0.00034);
                this.variables.set('t', 1e6 * this.variables.get('year_to_s'));
                this.variables.set('I', 1e23);
                this.variables.set('A_vort', this.variables.get('pi') * this.variables.get('r') * this.variables.get('r'));
                this.variables.set('omega1', 1e-5);
                this.variables.set('omega2', -1e-5);
                this.variables.set('v_exp', 1e4);
                this.variables.set('Vsys', 4.0 / 3.0 * this.variables.get('pi') * Math.pow(this.variables.get('r'), 3));
                this.variables.set('f_fluid', 1e-14);
                break;
            case SystemType.BUBBLE_NEBULA:
                this.variables.set('M', 100 * M_sun);
                this.variables.set('r', 4.73e16);
                this.variables.set('z', 0.0);
                this.variables.set('t', 1e5 * this.variables.get('year_to_s'));
                this.variables.set('I', 1e21);
                this.variables.set('A_vort', this.variables.get('pi') * this.variables.get('r') * this.variables.get('r'));
                this.variables.set('omega1', 1e-3);
                this.variables.set('omega2', -1e-3);
                this.variables.set('v_exp', 5e4);
                this.variables.set('Vsys', 4.0 / 3.0 * this.variables.get('pi') * Math.pow(this.variables.get('r'), 3));
                this.variables.set('f_fluid', 8.457e-14);
                break;
            case SystemType.ANTENNAE_GALAXIES:
                this.variables.set('M', 5e10 * M_sun);
                this.variables.set('r', 4.629e21);
                this.variables.set('z', 0.001);
                this.variables.set('t', 13.8e9 * this.variables.get('year_to_s'));
                this.variables.set('I', 1e24);
                this.variables.set('A_vort', this.variables.get('pi') * this.variables.get('r') * this.variables.get('r'));
                this.variables.set('omega1', 1e-6);
                this.variables.set('omega2', -1e-6);
                this.variables.set('v_exp', 2e5);
                this.variables.set('Vsys', 4.0 / 3.0 * this.variables.get('pi') * Math.pow(this.variables.get('r'), 3));
                this.variables.set('f_fluid', 4.228e-4);
                break;
            case SystemType.HORSEHEAD_NEBULA:
                this.variables.set('M', 100 * M_sun);
                this.variables.set('r', 9.46e15);
                this.variables.set('z', 0.0);
                this.variables.set('t', 1e6 * this.variables.get('year_to_s'));
                this.variables.set('I', 1e21);
                this.variables.set('A_vort', this.variables.get('pi') * this.variables.get('r') * this.variables.get('r'));
                this.variables.set('omega1', 1e-3);
                this.variables.set('omega2', -1e-3);
                this.variables.set('v_exp', 2e3);
                this.variables.set('Vsys', 4.0 / 3.0 * this.variables.get('pi') * Math.pow(this.variables.get('r'), 3));
                this.variables.set('f_fluid', 1e-13);
                break;
        }

        // Compute derived values
        this.variables.set('FDPM', this.computeFDPM());
        if (!this.variables.has('Vsys')) {
            this.variables.set('Vsys', this.computeVsys());
        }
    }

    // Dynamic variable operations
    updateVariable(name, value) {
        this.variables.set(name, value);
        if (name === 'r') {
            this.variables.set('A_vort', this.variables.get('pi') * value * value);
            this.variables.set('FDPM', this.computeFDPM());
            this.variables.set('Vsys', this.computeVsys());
        } else if (['I', 'omega1', 'omega2'].includes(name)) {
            this.variables.set('FDPM', this.computeFDPM());
        }
    }

    addToVariable(name, delta) {
        this.updateVariable(name, this.variables.get(name) + delta);
    }

    subtractFromVariable(name, delta) {
        this.updateVariable(name, this.variables.get(name) - delta);
    }

    // Helper methods
    computeHz(z) {
        const H0 = this.variables.get('H0');
        const Omega_m = this.variables.get('Omega_m');
        const Omega_Lambda = this.variables.get('Omega_Lambda');
        return H0 * Math.sqrt(Omega_m * Math.pow(1 + z, 3) + Omega_Lambda);
    }

    computeFDPM() {
        const I = this.variables.get('I');
        const A_vort = this.variables.get('A_vort');
        const omega1 = this.variables.get('omega1');
        const omega2 = this.variables.get('omega2');
        return I * A_vort * Math.abs(omega1 - omega2);
    }

    computeVsys() {
        const pi = this.variables.get('pi');
        const r = this.variables.get('r');
        return 4.0 / 3.0 * pi * Math.pow(r, 3);
    }

    computeEreact(t) {
        const E_react_base = this.variables.get('E_react_base');
        const decay_rate = this.variables.get('decay_rate');
        return E_react_base * Math.exp(-decay_rate * t);
    }

    computeFexp(t) {
        const z = this.variables.get('z');
        const Hz_val = this.computeHz(z);
        const Ht = Hz_val * t;
        const pi = this.variables.get('pi');
        return Ht / (2 * pi);
    }

    // Resonance term computations
    computeADPM() {
        const fdpm = this.variables.get('FDPM');
        const f_dpm = this.variables.get('f_DPM');
        const evac_neb = this.variables.get('Evac_neb');
        const c = this.variables.get('c');
        const vsys = this.variables.get('Vsys');
        return (fdpm * f_dpm * evac_neb) / (c * vsys);
    }

    computeATHz() {
        const a_dpm = this.computeADPM();
        const f_thz = this.variables.get('f_THz');
        const evac_neb = this.variables.get('Evac_neb');
        const v_exp = this.variables.get('v_exp');
        const evac_ism = this.variables.get('Evac_ISM');
        const c = this.variables.get('c');
        return (evac_ism / c) * f_thz * evac_neb * v_exp * a_dpm;
    }

    computeAvacDiff() {
        const a_dpm = this.computeADPM();
        const delta_evac = this.variables.get('Delta_Evac');
        const v_exp = this.variables.get('v_exp');
        const evac_neb = this.variables.get('Evac_neb');
        const c = this.variables.get('c');
        return (evac_neb / (c * c)) * delta_evac * (v_exp * v_exp) * a_dpm;
    }

    computeASuperFreq() {
        const a_dpm = this.computeADPM();
        const f_super = this.variables.get('F_super');
        const f_thz = this.variables.get('f_THz');
        const evac_neb = this.variables.get('Evac_neb');
        const c = this.variables.get('c');
        return (evac_neb / c) * f_super * f_thz * a_dpm;
    }

    computeAAetherRes() {
        const a_dpm = this.computeADPM();
        const ua_scm = this.variables.get('UA_SCm');
        const omega_i = this.variables.get('omega_i');
        const f_thz = this.variables.get('f_THz');
        const f_trz = this.variables.get('f_TRZ');
        return ua_scm * omega_i * f_thz * a_dpm * (1.0 + f_trz);
    }

    computeUg4i(t) {
        const a_dpm = this.computeADPM();
        const k4 = this.variables.get('k4');
        const e_react = this.computeEreact(t);
        const f_react = this.variables.get('f_react');
        const evac_neb = this.variables.get('Evac_neb');
        const c = this.variables.get('c');
        return (evac_neb / c) * k4 * e_react * f_react * a_dpm;
    }

    computeAQuantumFreq() {
        const a_dpm = this.computeADPM();
        const f_quantum = this.variables.get('f_quantum');
        const evac_neb = this.variables.get('Evac_neb');
        const evac_ism = this.variables.get('Evac_ISM');
        const c = this.variables.get('c');
        return (evac_ism / c) * f_quantum * evac_neb * a_dpm;
    }

    computeAAetherFreq() {
        const a_dpm = this.computeADPM();
        const f_aether = this.variables.get('f_Aether');
        const evac_neb = this.variables.get('Evac_neb');
        const evac_ism = this.variables.get('Evac_ISM');
        const c = this.variables.get('c');
        return (evac_ism / c) * f_aether * evac_neb * a_dpm;
    }

    computeAFluidFreq() {
        const f_fluid = this.variables.get('f_fluid');
        const evac_neb = this.variables.get('Evac_neb');
        const vsys = this.variables.get('Vsys');
        const evac_ism = this.variables.get('Evac_ISM');
        const c = this.variables.get('c');
        return (evac_ism / c) * f_fluid * evac_neb * vsys;
    }

    computeOscTerm(t) {
        return 0.0; // Approximation = 0
    }

    computeAExpFreq(t) {
        const f_exp = this.computeFexp(t);
        const a_dpm = this.computeADPM();
        const evac_neb = this.variables.get('Evac_neb');
        const evac_ism = this.variables.get('Evac_ISM');
        const c = this.variables.get('c');
        return (evac_ism / c) * f_exp * evac_neb * a_dpm;
    }

    // Core computation: Resonance MUGE g(r, t)
    computeG_resonance(t) {
        this.variables.set('t', t);
        const a_dpm = this.computeADPM();
        const a_thz = this.computeATHz();
        const a_vac_diff = this.computeAvacDiff();
        const a_super_freq = this.computeASuperFreq();
        const a_aether_res = this.computeAAetherRes();
        const ug4i = this.computeUg4i(t);
        const a_quantum_freq = this.computeAQuantumFreq();
        const a_aether_freq = this.computeAAetherFreq();
        const a_fluid_freq = this.computeAFluidFreq();
        const osc_term = this.computeOscTerm(t);
        const a_exp_freq = this.computeAExpFreq(t);
        const f_trz = this.variables.get('f_TRZ');

        return a_dpm + a_thz + a_vac_diff + a_super_freq + a_aether_res + ug4i +
               a_quantum_freq + a_aether_freq + a_fluid_freq + osc_term + a_exp_freq + f_trz;
    }

    // Output descriptive text
    getEquationText() {
        let sys_name;
        switch (this.currentSystem) {
            case SystemType.MAGNETAR_SGR_1745_2900: sys_name = "Magnetar SGR 1745-2900"; break;
            case SystemType.SAGITTARIUS_A: sys_name = "Sagittarius A*"; break;
            case SystemType.TAPESTRY_BLAZING_STARBIRTH: sys_name = "Tapestry of Blazing Starbirth"; break;
            case SystemType.WESTERLUND_2: sys_name = "Westerlund 2"; break;
            case SystemType.PILLARS_CREATION: sys_name = "Pillars of Creation"; break;
            case SystemType.RINGS_RELATIVITY: sys_name = "Rings of Relativity"; break;
            case SystemType.STUDENTS_GUIDE_UNIVERSE: sys_name = "Student's Guide to the Universe"; break;
            case SystemType.NGC_2525: sys_name = "NGC 2525"; break;
            case SystemType.NGC_3603: sys_name = "NGC 3603"; break;
            case SystemType.BUBBLE_NEBULA: sys_name = "Bubble Nebula"; break;
            case SystemType.ANTENNAE_GALAXIES: sys_name = "Antennae Galaxies"; break;
            case SystemType.HORSEHEAD_NEBULA: sys_name = "Horsehead Nebula"; break;
            default: sys_name = "Generic System";
        }
        return `Resonance Superconductive MUGE for ${sys_name}:
g(r,t) = a_DPM + a_THz + a_vac_diff + a_super_freq + a_aether_res + Ug4i + a_quantum_freq + a_Aether_freq + a_fluid_freq + Osc_term + a_exp_freq + f_TRZ
a_DPM = (F_DPM f_DPM E_vac,neb) / (c V_sys); a_THz = (E_vac,ISM / c) f_THz E_vac,neb v_exp a_DPM;
a_vac_diff = (E_vac,neb / c^2) ΔE_vac v_exp^2 a_DPM; a_super_freq = (E_vac,neb / c) F_super f_THz a_DPM;
a_aether_res = [U_A' : SC_m] ω_i f_THz a_DPM (1 + f_TRZ); Ug4i = (E_vac,neb / c) k4 E_react(t) f_react a_DPM;
a_quantum_freq = (E_vac,ISM / c) f_quantum E_vac,neb a_DPM; similar for a_Aether_freq, a_exp_freq;
a_fluid_freq = (E_vac,ISM / c) f_fluid E_vac,neb V_sys; Osc_term ≈ 0; f_exp = H(z) t / (2π).
Aether models expansion; tuned params yield g ~10^{-9} to 10^{35} m/s² (fluid dominant for large systems).`;
    }

    // Print all current variables
    printVariables() {
        console.log(`Resonance Variables for System ${this.currentSystem}:`);
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential(6)}`);
        }
    }

    // Self-expanding framework methods
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
            variables: Object.fromEntries(this.variables),
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            metadata: Object.fromEntries(this.metadata),
            currentSystem: this.currentSystem,
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

module.exports = {
    Source87UQFFModule,
    SystemType,
    PhysicsTerm,
    DynamicVacuumTerm,
    QuantumCouplingTerm
};