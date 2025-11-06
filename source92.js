class Source92UQFFModule {
    constructor() {
        this.variables = new Map();

        // Universal constants for buoyancy coupling
        this.variables.set('beta', 0.6);                        // β_i uniform (unitless)
        this.variables.set('Omega_g', 7.3e-16);                 // rad/s (galactic spin)
        this.variables.set('M_bh', 8.15e36);                    // kg (black hole mass)
        this.variables.set('d_g', 2.55e20);                     // m (galactic distance)
        this.variables.set('E_react', 1.0);                     // Reactive energy (normalized)
        this.variables.set('epsilon_sw', 0.001);                // Swirl factor
        this.variables.set('rho_vac_sw', 8e-21);                // J/m³
        this.variables.set('U_UA', 1.0);                        // Universal Aether factor
        this.variables.set('t_n', 0.0);                         // Time node (s)
        this.variables.set('pi', Math.PI);

        // U_gi defaults (example from doc for Ug1; others placeholder)
        this.variables.set('U_g1', 1.39e26);                    // J/m³
        this.variables.set('U_g2', 1e25);                       // Placeholder J/m³
        this.variables.set('U_g3', 1e24);                       // Placeholder J/m³
        this.variables.set('U_g4', 1e23);                       // Placeholder J/m³
    }

    // Dynamic variable operations
    setDynamicParameter(name, value) {
        this.variables.set(name, value);
    }

    addToVariable(name, delta) {
        const current = this.variables.get(name) || 0;
        this.variables.set(name, current + delta);
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    // Compute β_i (uniform 0.6)
    computeBeta(i) {
        return this.variables.get('beta');  // Uniform for all i
    }

    // Compute U_bi for specific i
    computeU_bi(i) {
        const ug_key = 'U_g' + i;
        const U_gi = this.variables.get(ug_key) || 1e26; // Default if not found
        const beta_i = this.computeBeta(i);
        const M_bh_over_d_g = this.variables.get('M_bh') / this.variables.get('d_g');
        const swirl_factor = 1.0 + this.variables.get('epsilon_sw') * this.variables.get('rho_vac_sw');
        const cos_term = Math.cos(this.variables.get('pi') * this.variables.get('t_n'));
        return -beta_i * U_gi * this.variables.get('Omega_g') * M_bh_over_d_g *
               this.variables.get('E_react') * swirl_factor * this.variables.get('U_UA') * cos_term;
    }

    // Compute all U_bi
    computeAllU_bi() {
        const u_bi = [];
        for (let i = 1; i <= 4; i++) {
            u_bi.push(this.computeU_bi(i));
        }
        return u_bi;
    }

    // Contribution to F_U (sum of -β_i terms)
    computeF_U_contribution() {
        const all_u_bi = this.computeAllU_bi();
        return all_u_bi.reduce((sum, val) => sum + val, 0);
    }

    // Get equation text
    getEquationText() {
        return 'U_bi = -β_i * U_gi * Ω_g * (M_bh / d_g) * E_react * (1 + ε_sw * ρ_vac,sw) * U_UA * cos(π t_n)\n' +
               'Where β_i = 0.6 (unitless, uniform for i=1-4: Ug1-Ug4);\n' +
               'Opposes gravity: 60% scaling of gravitational term.\n' +
               'In F_U: Σ [k_i U_gi - β_i U_gi Ω_g (M_bh/d_g) E_react] + other terms.\n' +
               'Role: Stabilizes systems (e.g., molecular clouds, nebulae); counteracts Ug collapse.\n' +
               'Example Ug1: U_b1 ≈ -1.94e27 J/m³ (at t_n=0, Sun params).\n' +
               'UQFF: Uniform buoyancy across scales; tunable for refinements.';
    }

    // Print all current variables
    printVariables() {
        console.log('Current Buoyancy Variables:');
        for (const [key, value] of this.variables) {
            console.log(`   ${key} = ${value.toExponential(3)}`);
        }
    }

    // Print U_bi
    printU_bi() {
        const all_u_bi = this.computeAllU_bi();
        console.log('Universal Buoyancy Terms U_bi (J/m³):');
        for (let i = 1; i <= 4; i++) {
            console.log(`   U_b${i} = ${all_u_bi[i-1].toExponential(3)}`);
        }
        console.log(`F_U Buoyancy Contribution (sum): ${this.computeF_U_contribution().toExponential(3)}`);
    }
}

module.exports = { Source92UQFFModule };