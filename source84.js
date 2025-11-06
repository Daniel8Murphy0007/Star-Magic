// LENRCalibUQFFModule.js
// JavaScript implementation of LENR neutron production calibration module
// Preserves all C++ dynamics: neutron production rate η, Um magnetism, non-local terms,
// scenario calibration (hydride/wires/corona), and dynamic variable management

class LENRCalibUQFFModule {
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

        // Core variables (preserved from C++)
        this.variables = new Map();

        // Universal constants
        this.variables.set('pi', Math.PI);
        this.variables.set('year_to_s', 3.156e7);               // s/yr
        this.variables.set('r', 1e-10);                         // m (default)
        this.variables.set('S_S_q', 1.0);                       // Non-local base

        // UQFF params
        this.variables.set('rho_vac_SCm', 7.09e-37);            // J/m³
        this.variables.set('rho_vac_UA', 7.09e-36);             // J/m³
        this.variables.set('rho_vac_UA_prime', 1e-23);          // For UA':SCm
        this.variables.set('gamma', 0.00005);                   // day^-1
        this.variables.set('t_n', 0.0);                         // days
        this.variables.set('P_scm', 1.0);                       // Polarization
        this.variables.set('E_react_0', 1e46);                  // Initial
        this.variables.set('omega_c', 2 * this.variables.get('pi') / 3.96e8);  // rad/s
        this.variables.set('f_heaviside', 0.01);
        this.variables.set('f_quasi', 0.01);

        // Calib defaults (overridden by scenario)
        this.variables.set('k_eta', 1e13);                      // cm^-2/s
        this.variables.set('t', 1.0 * this.variables.get('year_to_s'));  // 1 yr s
        this.variables.set('n', 1);                             // State

        this.current_scenario = 'hydride';
    }

    // Set scenario: Load calibrated params
    setScenario(scen_name) {
        this.current_scenario = scen_name;
        if (scen_name === 'hydride') {
            this.variables.set('k_eta', 1e13);  // cm^-2/s
            this.variables.set('E_target', 2e11);  // V/m
        } else if (scen_name === 'wires') {
            this.variables.set('k_eta', 1e8);
            this.variables.set('E_target', 28.8e11);
        } else if (scen_name === 'corona') {
            this.variables.set('k_eta', 7e-3);
            this.variables.set('E_target', 1.2e-3);
        }
    }

    // Dynamic variable operations
    updateVariable(name, value) {
        this.variables.set(name, value);
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
        }
    }

    subtractFromVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) - delta);
        }
    }

    // Mu_j computation
    computeMuJ(t) {
        const omega_c = this.variables.get('omega_c');
        return (1e3 + 0.4 * Math.sin(omega_c * t)) * 3.38e20;
    }

    // E_react computation
    computeEReact(t) {
        const E_react_0 = this.variables.get('E_react_0');
        const year_to_s = this.variables.get('year_to_s');
        return E_react_0 * Math.exp(-0.0005 * t / year_to_s);
    }

    // Um computation (universal magnetism)
    computeUm(t, r, n) {
        const mu_j = this.computeMuJ(t);
        const gamma = this.variables.get('gamma');
        const pi = this.variables.get('pi');
        const t_n = this.variables.get('t_n');
        const P_scm = this.variables.get('P_scm');
        const f_heaviside = this.variables.get('f_heaviside');
        const f_quasi = this.variables.get('f_quasi');

        const term1 = mu_j / r;
        const term2 = 1.0 - Math.exp(-gamma * (t / 86400) * Math.cos(pi * t_n));
        const factor = P_scm * this.computeEReact(t) * (1.0 + 1e13 * f_heaviside) * (1.0 + f_quasi);

        return term1 * term2 * factor;
    }

    // Electric field computation
    computeElectricField(um_val, rho_vac_val, r_val) {
        return um_val / (rho_vac_val * r_val);
    }

    // Delta_n computation
    computeDeltaN(n) {
        const pi = this.variables.get('pi');
        return Math.pow(2 * pi, n / 6.0);
    }

    // Rho_vac UA':SCm computation
    computeRhoVacUAScm(n, t) {
        const non_local = this.computeNonLocalExp(n, t);
        const rho_vac_UA_prime = this.variables.get('rho_vac_UA_prime');
        return rho_vac_UA_prime * Math.pow(0.1, n) * non_local;
    }

    // Non-local exponential computation
    computeNonLocalExp(n, t) {
        const pi = this.variables.get('pi');
        const year_to_s = this.variables.get('year_to_s');
        const S_S_q = this.variables.get('S_S_q');

        const exp_inner = Math.exp(-pi - t / year_to_s);
        const base = Math.pow(S_S_q, n) * Math.pow(2, 6);
        return Math.exp(-base * exp_inner);
    }

    // Eta computation (neutron production rate)
    computeEtaInternal(um_val, rho_vac_val, n, t) {
        const non_local = this.computeNonLocalExp(n, t);
        const k_eta = this.variables.get('k_eta');
        return k_eta * non_local * (um_val / rho_vac_val);
    }

    // Core computeEta (main interface)
    computeEta(t, n = 1) {
        this.variables.set('t', t);
        this.variables.set('n', n);
        const r = this.variables.get('r');
        const um = this.computeUm(t, r, n);
        const rho_vac_ua = this.variables.get('rho_vac_UA');
        return this.computeEtaInternal(um, rho_vac_ua, n, t);
    }

    // Unified field force computation (gravity equivalent for LENR)
    computeG(t = 0) {
        const eta = this.computeEta(t);
        const rho_vac = this.variables.get('rho_vac_UA');
        // Convert neutron production rate to unified field force
        return eta * rho_vac * 1e-20; // Scaled for consistency
    }

    // Environmental forces computation
    computeEnvironmentalForces(t = 0) {
        const um = this.computeUm(t, this.variables.get('r'), this.variables.get('n'));
        const rho_vac = this.variables.get('rho_vac_UA');
        const electric_field = this.computeElectricField(um, rho_vac, this.variables.get('r'));

        return {
            um_magnetic: um,
            electric_field: electric_field,
            vacuum_energy: rho_vac,
            neutron_rate: this.computeEta(t)
        };
    }

    // Equation text
    getEquationText() {
        return "η(t, n) = k_η * exp(-[S S_q]^n 2^6 e^(-π - t/yr)) * U_m / ρ_vac,[UA]\n" +
               "U_m(t,r,n) = μ [μ_j / r * (1 - e^{-γ t cos(π t_n)}) * ρ^j ] * P_scm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)\n" +
               "μ_j(t) = (10^3 + 0.4 sin(ω_c t)) * 3.38e20; E_react(t) = 10^46 e^{-0.0005 t/yr}\n" +
               "Δ_n = (2π)^{n/6}; ρ_vac,[UA']:SCm](n,t) = 10^{-23} (0.1)^n exp(-[S S_q]^n 2^6 e^(-π - t/yr))\n" +
               "E = U_m / (ρ_vac,[UA] r); Insights: Calib k_η for 100% accuracy; hydride η=1e13 cm^{-2}/s, E=2e11 V/m.\n" +
               "Adaptations: Pramana 2008; Scenarios: hydride/wires/corona. Solutions: η ~1e13 cm^{-2}/s (non-local dominant).";
    }

    // Print all current variables (for debugging)
    printVariables() {
        console.log(`LENR Calib Scenario: ${this.current_scenario}`);
        console.log('Variables:');
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential(6)}`);
        }
    }

    // Export state for cross-module communication
    exportState(filename = 'lenr_calib_state.txt') {
        const state = {
            scenario: this.current_scenario,
            variables: Object.fromEntries(this.variables),
            metadata: Object.fromEntries(this.metadata)
        };

        // In Node.js environment, this would write to file
        // For now, return the state object
        return state;
    }

    // Register dynamic term
    registerDynamicTerm(term) {
        this.dynamicTerms.push(term);
    }

    // Set dynamic parameter
    setDynamicParameter(name, value) {
        this.dynamicParameters.set(name, value);
    }

    // Get dynamic parameter
    getDynamicParameter(name) {
        return this.dynamicParameters.get(name);
    }

    // Set learning rate for optimization
    setLearningRate(rate) {
        this.learningRate = rate;
    }

    // Enable/disable logging
    setEnableLogging(enable) {
        this.enableLogging = enable;
    }
}

// Export for Node.js
if (typeof module !== 'undefined' && module.exports) {
    module.exports = LENRCalibUQFFModule;
}

// Example usage:
// const mod = new LENRCalibUQFFModule();
// mod.setScenario('hydride');
// const eta = mod.computeEta(3.156e7); // 1 year in seconds
// console.log('Neutron rate:', eta, 'cm^-2/s');
// mod.updateVariable('k_eta', 1.1e13);
// mod.printVariables();