// Source29.js - Sombrero Galaxy M104 UQFF Module (JavaScript Implementation from source29.cpp)

class SombreroUQFFModule {
    constructor() {
        // Initialize all variables with Sombrero-specific defaults (matching C++ source29.cpp)
        this.variables = new Map();

        // Base constants (universal)
        this.variables.set('G', 6.6743e-11);                    // m^3 kg^-1 s^-2
        this.variables.set('c', 3e8);                           // m/s
        this.variables.set('hbar', 1.0546e-34);                 // J s
        this.variables.set('Lambda', 1.1e-52);                  // m^-2 (cosmological constant)
        this.variables.set('q', 1.602e-19);                     // C (proton charge)
        this.variables.set('pi', 3.141592653589793);            // pi
        this.variables.set('t_Hubble', 13.8e9 * 3.156e7);       // s (13.8 Gyr)

        // Sombrero galaxy parameters
        const M_sun_val = 1.989e30;                             // kg
        this.variables.set('M_sun', M_sun_val);
        this.variables.set('M', 1e11 * M_sun_val);              // Total mass kg (incl. DM)
        this.variables.set('M_visible', 0.8 * this.variables.get('M'));  // Visible mass fraction (est. for bulge/arms)
        this.variables.set('M_DM', 0.2 * this.variables.get('M'));       // Dark matter mass (halo dominant but lower fraction)
        this.variables.set('r', 2.36e20);                       // m (half diameter ~25k ly)
        this.variables.set('M_BH', 1e9 * M_sun_val);            // SMBH kg
        this.variables.set('r_BH', 1e15);                       // m (core scale)

        // Hubble/cosmology
        this.variables.set('H0', 70.0);                         // km/s/Mpc
        this.variables.set('Mpc_to_m', 3.086e22);               // m/Mpc
        this.variables.set('z', 0.0063);                        // Redshift
        this.variables.set('Omega_m', 0.3);
        this.variables.set('Omega_Lambda', 0.7);
        this.variables.set('t', 10e9 * 3.156e7);                // Default t=10 Gyr s

        // Dust/fluid dynamics
        this.variables.set('rho_dust', 1e-20);                  // kg/m^3
        this.variables.set('v_orbit', 2e5);                     // m/s
        this.variables.set('rho_mass', 1e-21);                  // kg/m^3 (ISM)
        this.variables.set('rho_fluid', 1e-21);                 // kg/m^3 (fluid density, dust lane-like)
        this.variables.set('V', 1e3);                           // m^3 (arbitrary volume scale)

        // EM/magnetic/superconductivity
        this.variables.set('B', 1e-5);                          // T (galactic field)
        this.variables.set('B_crit', 1e15 * 1e-4);              // T (10^15 G = 1e11 T, but doc 10^15 G ≈1e11 T)

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
        this.variables.set('delta_rho', 0.1 * this.variables.get('rho_mass'));  // Perturbation
        this.variables.set('rho', this.variables.get('rho_mass'));               // Mean density

        // Ug subterms (computed dynamically, but init placeholders)
        this.variables.set('Ug1', 0.0);  // Will be G M / r^2
        this.variables.set('Ug2', 0.0);  // d^2 Phi / dt^2 ≈ 0 (negligible)
        this.variables.set('Ug3', 0.0);  // G M_moon / r_moon^2 ≈ 0 (no moon)
        this.variables.set('Ug4', 0.0);  // Ug1 * f_sc, f_sc=1

        // Scale factors (from streamlining)
        this.variables.set('scale_macro', 1e-12);               // For macro effects
        this.variables.set('f_TRZ', 0.1);                       // Time-reversal factor (from May09)
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
                this.variables.set('M_visible', 0.8 * value);
                this.variables.set('M_DM', 0.2 * value);
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

    // Fluid term: rho_fluid * V * g (g approx base grav)
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

    // Dust term: rho_dust * v_orbit^2 / rho_mass * scale_macro (as a_dust)
    computeDustTerm() {
        const force_dust = this.variables.get('rho_dust') * (this.variables.get('v_orbit') * this.variables.get('v_orbit'));
        return (force_dust / this.variables.get('rho_mass')) * this.variables.get('scale_macro');
    }

    // Full computation: g_UQFF(r, t) = ... all terms, incl. superconductivity (1 - B/B_crit) on base
    computeG(t) {
        this.variables.set('t', t);  // Update t
        const Hz = this.computeHz();
        const expansion = 1.0 + Hz * t;
        const sc_correction = 1.0 - (this.variables.get('B') / this.variables.get('B_crit'));
        const tr_factor = 1.0 + this.variables.get('f_TRZ');

        // Base gravity with expansion, SC, TR
        const g_base = ((this.variables.get('G') * this.variables.get('M') / (this.variables.get('r') * this.variables.get('r'))) * expansion * sc_correction) * tr_factor;

        // BH term
        const g_BH = (this.variables.get('G') * this.variables.get('M_BH')) / (this.variables.get('r_BH') * this.variables.get('r_BH'));

        // Ug sum
        const ug_sum = this.computeUgSum();

        // Cosmological
        const lambda_term = this.variables.get('Lambda') * (this.variables.get('c') * this.variables.get('c')) / 3.0;

        // Quantum
        const quantum_term = this.computeQuantumTerm(this.variables.get('t_Hubble'));

        // EM Lorentz (magnitude v B)
        const em_base = this.variables.get('q') * this.variables.get('v_orbit') * this.variables.get('B') / 1.673e-27;
        const em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * this.variables.get('scale_macro');

        // Fluid (uses g_base approx)
        const fluid_term = this.computeFluidTerm(g_base);

        // Resonant
        const resonant_term = this.computeResonantTerm(t);

        // DM
        const dm_term = this.computeDMTerm();

        // Dust
        const dust_term = this.computeDustTerm();

        // Total: Sum all
        return g_base + g_BH + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + dust_term;
    }

    // Get equation text (descriptive)
    getEquationText() {
        return "g_Sombrero(r, t) = (G * M / r^2) * (1 + H(z) * t) * (1 - B / B_crit) * (1 + f_TRZ) + (G * M_BH / r_BH^2) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 / 3) + " +
            "(hbar / sqrt(Delta_x * Delta_p)) * ∫(ψ* H ψ dV) * (2π / t_Hubble) + q (v × B) + ρ_fluid * V * g + " +
            "2 A cos(k x) cos(ω t) + (2π / 13.8) A exp(i (k x - ω t)) + (M_visible + M_DM) * (δρ/ρ + 3 G M / r^3) + D_dust\n" +
            "Special Terms:\n" +
            "- Quantum: Heisenberg uncertainty with normalized wavefunction integral (ground state approx).\n" +
            "- Fluid: ISM-like density-volume-gravity coupling in dust lane.\n" +
            "- Resonant: Oscillatory Aether-mediated waves (real part of complex exp) for globular cluster dynamics.\n" +
            "- DM: Visible+dark mass with density perturbations and curvature term for halo.\n" +
            "- Superconductivity: (1 - B/B_crit) for quantum field effects.\n" +
            "Solutions: Numerical evaluation at t=10 Gyr yields ~0.535 m/s² (dust/BH dominant; full sum includes micro terms ~1e-10 to 1e-3).\n" +
            "Adaptations for Sombrero: Virgo Cluster z=0.0063; prominent dust lane boosts D_dust; SMBH=1e9 Msun shapes bulge.";
    }

    // Print variables
    printVariables() {
        console.log('Current Sombrero Variables:');
        for (const [key, value] of this.variables) {
            console.log(`${key} = ${value.toExponential()}`);
        }
    }
}

module.exports = { SombreroUQFFModule };