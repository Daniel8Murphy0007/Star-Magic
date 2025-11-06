// Source32.js - Crab Nebula UQFF Module (JavaScript Implementation from source32.cpp)

class CrabUQFFModule {
    constructor() {
        // Initialize all variables with Crab Nebula defaults (matching C++ source32.cpp)
        this.variables = new Map();

        // Base constants (universal)
        this.variables.set('G', 6.6743e-11);                    // m^3 kg^-1 s^-2
        this.variables.set('c', 3e8);                           // m/s
        this.variables.set('hbar', 1.0546e-34);                 // J s
        this.variables.set('Lambda', 1.1e-52);                  // m^-2 (cosmological constant)
        this.variables.set('q', 1.602e-19);                     // C (electron charge)
        this.variables.set('pi', 3.141592653589793);            // pi
        this.variables.set('t_Hubble', 13.8e9 * 3.156e7);       // s (13.8 Gyr)

        // Crab Nebula parameters
        const M_sun_val = 1.989e30;                             // kg
        this.variables.set('M_sun', M_sun_val);
        this.variables.set('M', 4.6 * M_sun_val);               // Total mass kg
        this.variables.set('M_visible', this.variables.get('M')); // Visible mass (ejecta + pulsar)
        this.variables.set('M_DM', 0.0);                        // No significant DM
        this.variables.set('r0', 5.2e16);                       // m (initial radius)
        this.variables.set('v_exp', 1.5e6);                     // m/s (expansion velocity)

        // Hubble/cosmology
        this.variables.set('H0', 70.0);                         // km/s/Mpc
        this.variables.set('Mpc_to_m', 3.086e22);               // m/Mpc
        this.variables.set('z', 0.0015);                        // Redshift
        this.variables.set('Omega_m', 0.3);
        this.variables.set('Omega_Lambda', 0.7);
        this.variables.set('t', 971 * 3.156e7);                 // Default t=971 years s (since 1054 AD)

        // Nebula dynamics
        this.variables.set('rho_fluid', 1e-21);                 // kg/m^3 (filament density)
        this.variables.set('V', 1e3);                           // m^3 (arbitrary volume scale)
        this.variables.set('v_shock', 1.5e6);                   // m/s (shock velocity)
        this.variables.set('P_pulsar', 5e31);                   // W (pulsar luminosity)
        this.variables.set('delta_rho', 0.1 * this.variables.get('rho_fluid')); // Perturbation
        this.variables.set('rho', this.variables.get('rho_fluid')); // Mean density

        // EM/magnetic/superconductivity
        this.variables.set('B', 1e-8);                          // T (nebula average magnetic field)
        this.variables.set('B_crit', 1e11);                     // T (10^15 G ≈ 1e11 T)
        this.variables.set('m_e', 9.11e-31);                    // kg (electron mass)

        // Quantum terms
        this.variables.set('Delta_x', 1e-10);                   // m (position uncertainty, atomic scale)
        this.variables.set('Delta_p', this.variables.get('hbar') / this.variables.get('Delta_x')); // Momentum uncertainty (Heisenberg)
        this.variables.set('integral_psi', 1.0);                // Normalized <psi|H|psi> dV ≈ E_ground (simplified to 1 for unitless)

        // Resonant/oscillatory terms
        this.variables.set('A', 1e-10);                         // Amplitude (arbitrary small)
        this.variables.set('k', 1e20);                          // m^-1 (wave number, short wavelength)
        this.variables.set('omega', 1e15);                      // rad/s (high freq, e.g., synchrotron)
        this.variables.set('x', 0.0);                           // m (position, central)

        // Ug subterms (computed dynamically, but init placeholders)
        this.variables.set('Ug1', 0.0);  // Will be G M / r^2
        this.variables.set('Ug2', 0.0);  // d^2 Phi / dt^2 ≈ 0 (negligible)
        this.variables.set('Ug3', 0.0);  // G M_moon / r_moon^2 ≈ 0 (no moon)
        this.variables.set('Ug4', 0.0);  // Ug1 * f_sc, f_sc=1

        // Scale factors (from streamlining)
        this.variables.set('scale_macro', 1e-12);               // For macro effects
        this.variables.set('f_TRZ', 0.1);                       // Time-reversal factor
        this.variables.set('f_sc', 1.0);                        // Superconductive factor
    }

    // Dynamic variable operations (maintaining C++ interface)
    updateVariable(name, value) {
        if (this.variables.has(name)) {
            this.variables.set(name, value);
            // Recompute dependent vars if needed
            if (name === 'Delta_x') {
                this.variables.set('Delta_p', this.variables.get('hbar') / value);
            } else if (name === 'M') {
                this.variables.set('M_visible', value); // For nebula
                this.variables.set('M_DM', 0.0);
            }
        } else {
            console.log(`Variable '${name}' not found. Adding with value ${value}`);
            this.variables.set(name, value);
        }
    }

    getVariable(name) {
        if (this.variables.has(name)) {
            return this.variables.get(name);
        } else {
            console.log(`Variable '${name}' not found.`);
            return undefined;
        }
    }

    addToVariable(name, delta) {
        if (this.variables.has(name)) {
            this.variables.set(name, this.variables.get(name) + delta);
        } else {
            console.log(`Variable '${name}' not found. Adding with delta ${delta}`);
            this.variables.set(name, delta);
        }
    }

    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }

    // Compute H(z) in s^-1
    computeHz() {
        const Hz_kms = this.variables.get('H0') * Math.sqrt(
            this.variables.get('Omega_m') * Math.pow(1.0 + this.variables.get('z'), 3) +
            this.variables.get('Omega_Lambda')
        );
        return (Hz_kms * 1e3) / this.variables.get('Mpc_to_m');
    }

    // Compute Ug sum: Ug1 = G M / r^2, Ug4 = Ug1 * f_sc, others 0
    computeUgSum() {
        const r = this.variables.get('r0') + this.variables.get('v_exp') * this.variables.get('t');  // Use current r(t)
        const Ug1 = (this.variables.get('G') * this.variables.get('M')) / (r * r);
        this.variables.set('Ug1', Ug1);
        this.variables.set('Ug4', Ug1 * this.variables.get('f_sc'));
        return this.variables.get('Ug1') + this.variables.get('Ug2') + this.variables.get('Ug3') + this.variables.get('Ug4');
    }

    // Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
    computeQuantumTerm(t_Hubble_val) {
        const unc = Math.sqrt(this.variables.get('Delta_x') * this.variables.get('Delta_p'));
        const integral_val = this.variables.get('integral_psi');
        return (this.variables.get('hbar') / unc) * integral_val * (2 * this.variables.get('pi') / t_Hubble_val);
    }

    // Fluid term: rho_fluid * V * g (g approx base g_grav)
    computeFluidTerm(g_base) {
        return this.variables.get('rho_fluid') * this.variables.get('V') * g_base;
    }

    // Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
    computeResonantTerm(t) {
        const cos_term = 2 * this.variables.get('A') * Math.cos(this.variables.get('k') * this.variables.get('x')) * Math.cos(this.variables.get('omega') * t);
        const exp_term_real = this.variables.get('A') * Math.cos(this.variables.get('k') * this.variables.get('x') - this.variables.get('omega') * t);
        const exp_factor = (2 * this.variables.get('pi') / 13.8);
        return cos_term + exp_factor * exp_term_real;
    }

    // DM term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3)
    computeDMTerm() {
        const r = this.variables.get('r0') + this.variables.get('v_exp') * this.variables.get('t');  // Use current r(t)
        const pert = this.variables.get('delta_rho') / this.variables.get('rho');
        const curv = 3 * this.variables.get('G') * this.variables.get('M') / (r * r * r);
        return (this.variables.get('M_visible') + this.variables.get('M_DM')) * (pert + curv);
    }

    // Wind term: (P_pulsar / (4 pi r^2)) * (1 + v_shock / c) / rho_fluid * scale_macro
    computeWindTerm(r) {
        const pressure = (this.variables.get('P_pulsar') / (4 * this.variables.get('pi') * r * r)) * (1.0 + this.variables.get('v_shock') / this.variables.get('c'));
        return (pressure / this.variables.get('rho_fluid')) * this.variables.get('scale_macro');
    }

    // Magnetic term: (q * v_shock * B) / m_e * scale_macro
    computeMagTerm() {
        const force = this.variables.get('q') * this.variables.get('v_shock') * this.variables.get('B');
        return (force / this.variables.get('m_e')) * this.variables.get('scale_macro');
    }

    // Full computation: g_UQFF(r, t) with all Crab terms
    computeG(t) {
        this.variables.set('t', t);  // Update t
        const r = this.variables.get('r0') + this.variables.get('v_exp') * t;  // r(t)
        const Hz = this.computeHz();
        const expansion = 1.0 + Hz * t;
        const sc_correction = 1.0 - (this.variables.get('B') / this.variables.get('B_crit'));
        const tr_factor = 1.0 + this.variables.get('f_TRZ');

        // Base gravity with expansion, SC, TR
        const g_base = (this.variables.get('G') * this.variables.get('M') / (r * r)) * expansion * sc_correction * tr_factor;

        // Ug sum
        const ug_sum = this.computeUgSum();

        // Cosmological
        const lambda_term = this.variables.get('Lambda') * (this.variables.get('c') * this.variables.get('c')) / 3.0;

        // Quantum
        const quantum_term = this.computeQuantumTerm(this.variables.get('t_Hubble'));

        // EM Lorentz (magnitude v_shock B)
        const em_base = this.variables.get('q') * this.variables.get('v_shock') * this.variables.get('B') / 1.673e-27;
        const em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * this.variables.get('scale_macro');

        // Fluid (uses g_base approx)
        const fluid_term = this.computeFluidTerm(g_base);

        // Resonant
        const resonant_term = this.computeResonantTerm(t);

        // DM
        const dm_term = this.computeDMTerm();

        // Wind
        const wind_term = this.computeWindTerm(r);

        // Mag
        const mag_term = this.computeMagTerm();

        // Total: Sum all
        return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + wind_term + mag_term;
    }

    // Get equation text (descriptive)
    getEquationText() {
        return "g_Crab(r, t) = (G * M / r(t)^2) * (1 + H(z) * t) * (1 - B / B_crit) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 / 3) + " +
            "(hbar / sqrt(Delta_x * Delta_p)) * ∫(ψ* H ψ dV) * (2π / t_Hubble) + q (v × B) + ρ_fluid * V * g + " +
            "2 A cos(k x) cos(ω t) + (2π / 13.8) A exp(i (k x - ω t)) + (M_visible + M_DM) * (δρ/ρ + 3 G M / r^3) + a_wind + M_mag\n" +
            "Where r(t) = r0 + v_exp * t; a_wind = [P_pulsar / (4π r^2) * (1 + v_shock / c)] / ρ * 1e-12; M_mag = (q v_shock B) / m_e * 1e-12\n" +
            "Special Terms:\n" +
            "- Quantum: Heisenberg uncertainty with normalized wavefunction integral (ground state approx) for particle quantum effects.\n" +
            "- Fluid: Nebular filament density-volume-gravity coupling.\n" +
            "- Resonant: Oscillatory Aether-mediated waves (real part of complex exp) for wisp dynamics.\n" +
            "- DM: Visible mass (ejecta + pulsar) with density perturbations and curvature term (M_DM=0).\n" +
            "- Superconductivity: (1 - B/B_crit) for quantum field effects near pulsar.\n" +
            "- Pulsar Wind: a_wind from relativistic wind pressure, dominant outward force.\n" +
            "- Magnetic: M_mag from Lorentz force on electrons in nebula fields.\n" +
            "Solutions: Numerical evaluation at t=971 yr yields ~1.481e6 m/s² (a_wind dominant; g_grav ~2e-13; micro terms ~1e-10 to 1e-3).\n" +
            "Adaptations for Crab: Pulsar-driven remnant with r(t); z=0.0015; v_shock=1.5e6 m/s boosts wind/mag.";
    }

    // Print variables
    printVariables() {
        console.log('Current Crab Variables:');
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }
}

module.exports = { CrabUQFFModule };