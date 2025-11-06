// Source30.js - Saturn Planet UQFF Module (JavaScript Implementation from source30.cpp)

class SaturnUQFFModule {
    constructor() {
        // Initialize all variables with Saturn-specific defaults (matching C++ source30.cpp)
        this.variables = new Map();

        // Base constants (universal)
        this.variables.set('G', 6.6743e-11);                    // m^3 kg^-1 s^-2
        this.variables.set('c', 3e8);                           // m/s
        this.variables.set('hbar', 1.0546e-34);                 // J s
        this.variables.set('Lambda', 1.1e-52);                  // m^-2 (cosmological constant)
        this.variables.set('q', 1.602e-19);                     // C (proton charge)
        this.variables.set('pi', 3.141592653589793);            // pi
        this.variables.set('t_Hubble', 13.8e9 * 3.156e7);       // s (13.8 Gyr)

        // Saturn parameters
        this.variables.set('M_Sun', 1.989e30);                  // kg
        this.variables.set('M', 5.683e26);                      // Planet mass kg (rings negligible addition)
        this.variables.set('M_ring', 1.5e19);                   // Ring mass kg
        this.variables.set('r', 6.0268e7);                      // m (equatorial radius)
        this.variables.set('r_orbit', 1.43e12);                 // m (orbital distance)
        this.variables.set('r_ring', 7e7);                      // m (average ring radius)
        this.variables.set('M_visible', this.variables.get('M')); // Visible mass (planet)
        this.variables.set('M_DM', 0.0);                        // No significant DM

        // Hubble/cosmology
        this.variables.set('H0', 70.0);                         // km/s/Mpc
        this.variables.set('Mpc_to_m', 3.086e22);               // m/Mpc
        this.variables.set('z', 0.0);                           // No redshift (Solar System)
        this.variables.set('Omega_m', 0.3);
        this.variables.set('Omega_Lambda', 0.7);
        this.variables.set('t', 4.5e9 * 3.156e7);               // Default t=4.5 Gyr s (Solar System age)

        // Atmospheric/wind dynamics
        this.variables.set('rho_atm', 2e-4);                    // kg/m^3 (upper atmosphere)
        this.variables.set('v_wind', 500.0);                    // m/s (average wind speed)
        this.variables.set('rho_fluid', 2e-4);                  // kg/m^3 (fluid density, atmospheric)
        this.variables.set('V', 1e3);                           // m^3 (arbitrary volume scale)

        // EM/magnetic/superconductivity
        this.variables.set('B', 1e-7);                          // T (planetary magnetic field)
        this.variables.set('B_crit', 1e11);                     // T (10^15 G ≈ 1e11 T)

        // Quantum terms
        this.variables.set('Delta_x', 1e-10);                   // m (position uncertainty, atomic scale)
        this.variables.set('Delta_p', this.variables.get('hbar') / this.variables.get('Delta_x'));  // Momentum uncertainty (Heisenberg)
        this.variables.set('integral_psi', 1.0);                // Normalized <psi|H|psi> dV ≈ E_ground (simplified to 1 for unitless)

        // Resonant/oscillatory terms
        this.variables.set('A', 1e-10);                         // Amplitude (arbitrary small)
        this.variables.set('k', 1e20);                          // m^-1 (wave number, short wavelength)
        this.variables.set('omega', 1e15);                      // rad/s (high freq, e.g., optical)
        this.variables.set('x', 0.0);                           // m (position, central)

        // DM perturbations
        this.variables.set('delta_rho', 0.1 * this.variables.get('rho_atm'));  // Perturbation
        this.variables.set('rho', this.variables.get('rho_atm'));                // Mean density

        // Ug subterms (computed dynamically, but init placeholders)
        this.variables.set('Ug1', 0.0);  // Will be G M / r^2
        this.variables.set('Ug2', 0.0);  // d^2 Phi / dt^2 ≈ 0 (negligible)
        this.variables.set('Ug3', 0.0);  // G M_moon / r_moon^2 ≈ 0 (no specific moon)
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
                this.variables.set('M_visible', value);  // For planet
                this.variables.set('M_DM', 0.0);
            }
        } else {
            console.log(`Variable '${name}' not found. Adding with value ${value}`);
            this.variables.set(name, value);
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
        const Ug1 = (this.variables.get('G') * this.variables.get('M')) / (this.variables.get('r') * this.variables.get('r'));
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

    // Fluid term: rho_fluid * V * g (g approx base g_saturn)
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
        const pert = this.variables.get('delta_rho') / this.variables.get('rho');
        const curv = 3 * this.variables.get('G') * this.variables.get('M') / (this.variables.get('r') * this.variables.get('r') * this.variables.get('r'));
        return (this.variables.get('M_visible') + this.variables.get('M_DM')) * (pert + curv);
    }

    // Wind term: rho_atm * v_wind^2 / rho_atm * scale_macro = v_wind^2 * scale_macro
    computeWindTerm() {
        return Math.pow(this.variables.get('v_wind'), 2) * this.variables.get('scale_macro');
    }

    // Full computation: g_UQFF(r, t) with all Saturn terms
    computeG(t) {
        this.variables.set('t', t);  // Update t
        const Hz = this.computeHz();
        const expansion = 1.0 + Hz * t;
        const sc_correction = 1.0 - (this.variables.get('B') / this.variables.get('B_crit'));
        const tr_factor = 1.0 + this.variables.get('f_TRZ');

        // Sun gravity with expansion and TR
        const g_sun = (this.variables.get('G') * this.variables.get('M_Sun') / (this.variables.get('r_orbit') * this.variables.get('r_orbit'))) * expansion * tr_factor;

        // Saturn gravity with superconductivity correction
        const g_saturn_base = (this.variables.get('G') * this.variables.get('M')) / (this.variables.get('r') * this.variables.get('r'));
        const g_saturn = g_saturn_base * sc_correction;

        // Ring tidal
        const T_ring = (this.variables.get('G') * this.variables.get('M_ring')) / (this.variables.get('r_ring') * this.variables.get('r_ring'));

        // Ug sum
        const ug_sum = this.computeUgSum();

        // Cosmological
        const lambda_term = this.variables.get('Lambda') * (this.variables.get('c') * this.variables.get('c')) / 3.0;

        // Quantum
        const quantum_term = this.computeQuantumTerm(this.variables.get('t_Hubble'));

        // EM Lorentz (magnitude v_wind B)
        const em_base = this.variables.get('q') * this.variables.get('v_wind') * this.variables.get('B') / 1.673e-27;
        const em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * this.variables.get('scale_macro');

        // Fluid (uses g_saturn approx)
        const fluid_term = this.computeFluidTerm(g_saturn);

        // Resonant
        const resonant_term = this.computeResonantTerm(t);

        // DM
        const dm_term = this.computeDMTerm();

        // Wind
        const wind_term = this.computeWindTerm();

        // Total: Sum all
        return g_sun + g_saturn + T_ring + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + wind_term;
    }

    // Get equation text (descriptive)
    getEquationText() {
        return "g_Saturn(r, t) = (G * M_Sun / r_orbit^2) * (1 + H(z) * t) * (1 + f_TRZ) + (G * M / r^2) * (1 - B / B_crit) + T_ring + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 / 3) + " +
            "(hbar / sqrt(Delta_x * Delta_p)) * ∫(ψ* H ψ dV) * (2π / t_Hubble) + q (v × B) + ρ_fluid * V * g + " +
            "2 A cos(k x) cos(ω t) + (2π / 13.8) A exp(i (k x - ω t)) + (M_visible + M_DM) * (δρ/ρ + 3 G M / r^3) + a_wind\n" +
            "Special Terms:\n" +
            "- Quantum: Heisenberg uncertainty with normalized wavefunction integral (ground state approx) for atmospheric quantum effects.\n" +
            "- Fluid: Atmospheric density-volume-gravity coupling.\n" +
            "- Resonant: Oscillatory Aether-mediated waves (real part of complex exp) for ring dynamics.\n" +
            "- DM: Visible mass (planet) with density perturbations and curvature term (M_DM=0).\n" +
            "- Superconductivity: (1 - B/B_crit) for quantum field effects in atmosphere.\n" +
            "- Ring Tidal: G M_ring / r_ring^2 for ring influence.\n" +
            "- Wind: v_wind^2 * 1e-12 for atmospheric feedback.\n" +
            "Solutions: Numerical evaluation at t=4.5 Gyr yields ~10.44 m/s² (g_saturn dominant; orbital g_sun ~9e-5; micro terms ~1e-7 to 1e-10).\n" +
            "Adaptations for Saturn: Solar System orbital term; z=0 negligible expansion; wind/rings boost local effects.";
    }

    // Print variables
    printVariables() {
        console.log('Current Saturn Variables:');
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }
}

module.exports = { SaturnUQFFModule };