// source50.js - UQFF Module for Compressed and Resonance Equations
// Implements dynamic variable management for multiple astrophysical systems

class SystemData {
    constructor(name, description, M, r, z, t, V, F_env, v_exp, I, A, omega1, omega2, M_sun = 0, r_orbit = 0) {
        this.name = name;
        this.description = description;
        this.M = M;
        this.r = r;
        this.z = z;
        this.t = t;
        this.V = V;
        this.F_env = F_env;
        this.v_exp = v_exp;
        this.I = I;
        this.A = A;
        this.omega1 = omega1;
        this.omega2 = omega2;
        this.M_sun = M_sun;
        this.r_orbit = r_orbit;

        // Dynamic variables map (equivalent to std::map<std::string, double>)
        this.vars = new Map();

        // Initialize with base values
        this.vars.set("M", M);
        this.vars.set("r", r);
        this.vars.set("z", z);
        this.vars.set("t", t);
        this.vars.set("V", V);
        this.vars.set("F_env", F_env);
        this.vars.set("v_exp", v_exp);
        this.vars.set("I", I);
        this.vars.set("A", A);
        this.vars.set("omega1", omega1);
        this.vars.set("omega2", omega2);
        this.vars.set("M_sun", M_sun);
        this.vars.set("r_orbit", r_orbit);

        // Computed variables
        this.vars.set("H_t_z", 0);
        this.vars.set("one_plus_H_t", 0);
        this.vars.set("B_adjust", 0);
        this.vars.set("one_plus_F_env", 0);
        this.vars.set("Lambda_c2_3", 0);
        this.vars.set("hbar_over_sqrt_delta", 0);
        this.vars.set("quantum_term", 0);
        this.vars.set("rho_V_g", 0);
        this.vars.set("three_G_M_over_r3", 0);
        this.vars.set("density_pert", 0);
        this.vars.set("M_vis_DM_pert", 0);
        this.vars.set("H_z", 0);
    }
}

class UQFFModule {
    constructor() {
        this.systems = new Map(); // Map<string, SystemData>
        this.installed = false;

        // Physical constants
        this.G = 6.67430e-11;
        this.c = 299792458;
        this.hbar = 1.0545718e-34;
        this.H0 = 2.2e-18; // Hubble constant in s^-1
        this.Lambda = 1.1056e-52; // Cosmological constant
        this.pi = Math.PI;
        this.g_earth = 9.80665;
        this.B_t = 1e-9; // Current magnetic field
        this.B_crit = 1e-8; // Critical magnetic field
        this.Delta_x_Delta_p = 1e-68; // Uncertainty principle
        this.integral_psi = 1e-10; // Wave function integral
        this.t_Hubble = 4.35e17; // Hubble time
        this.rho_fluid = 1e3; // Fluid density
        this.delta_rho_over_rho = 1e-5; // Density perturbation
        this.M_DM_default = 1e40; // Default DM mass

        // Resonance constants
        this.E_vac_neb = 1e-9; // Vacuum energy in nebula
        this.E_vac_ISM = 1e-10; // Vacuum energy in ISM
        this.Delta_E_vac = 1e-8; // Vacuum energy differential
        this.F_super = 1e12; // Superconductor frequency factor
        this.k_4 = 1e-40; // Aether constant
        this.omega_i = 1e15; // Initial frequency
        this.UA_SC_m = 1e-20; // Superconductive correction
        this.f_quantum = 1e14; // Quantum frequency
        this.f_Aether = 1e13; // Aether frequency
        this.f_TRZ = 1e-15; // TRZ frequency
    }

    // Compute volume if not provided (4/3 pi r^3)
    compute_volume(r) {
        return (4.0 / 3.0) * this.pi * r * r * r;
    }

    // Install all systems with defaults
    install_uqff_module() {
        // Clear existing
        this.systems.clear();

        // Hubble Sees Galaxies Galore
        let V_gal = this.compute_volume(1.543e21);
        this.systems.set("Hubble Sees Galaxies Galore", new SystemData(
            "Hubble Sees Galaxies Galore", "Hubble Deep Field observations, capturing thousands of galaxies.",
            1.989e41, 1.543e21, 1.0, 4.35e17, V_gal, 0.0, 1e5, 1e24, 7.487e42, 1e-6, -1e-6
        ));

        // The Stellar Forge
        let V_stellar = this.compute_volume(9.46e16);
        this.systems.set("The Stellar Forge", new SystemData(
            "The Stellar Forge", "Star-forming region in Large Magellanic Cloud (30 Doradus Nebula).",
            1.989e34, 9.46e16, 0.00005, 6.312e13, V_stellar, 0.0, 1e4, 1e22, 8.508e35, 1e-2, -1e-2
        ));

        // Hubble Mosaic of the Majestic Sombrero Galaxy
        let V_sombrero = this.compute_volume(4.73e20);
        this.systems.set("Hubble Mosaic of the Majestic Sombrero Galaxy", new SystemData(
            "Hubble Mosaic of the Majestic Sombrero Galaxy", "Sombrero Galaxy (M104), peculiar galaxy with dust lane.",
            1.591e42, 4.73e20, 0.002, 4.35e17, V_sombrero, 0.0, 2e5, 1e24, 7.487e42, 1e-6, -1e-6
        ));

        // Saturn
        let V_saturn = this.compute_volume(6.027e7);
        this.systems.set("Saturn", new SystemData(
            "Saturn", "Hubble observations of Saturn, rings and atmosphere.",
            5.68e26, 6.027e7, 0.0, 4.35e17, V_saturn, 0.0, 5e3, 1e20, 7.032e22, 1e-4, -1e-4, 1.989e30, 1.36e12
        ));

        // New Stars Shed Light on the Past
        this.systems.set("New Stars Shed Light on the Past", new SystemData(
            "New Stars Shed Light on the Past", "Star-forming region in Small Magellanic Cloud (N90).",
            1.989e34, 9.46e16, 0.00006, 6.312e13, V_stellar, 0.0, 1e4, 1e22, 8.508e35, 1e-2, -1e-2
        ));

        // The Crab Nebula
        let V_crab = this.compute_volume(5.203e16);
        this.systems.set("The Crab Nebula", new SystemData(
            "The Crab Nebula", "Supernova remnant formed in 1054 CE.",
            9.945e30, 5.203e16, 0.00002, 3.064e10, V_crab, 0.0, 1.34e6, 1e22, 8.508e35, 1e-2, -1e-2
        ));

        // Student's Guide to the Universe
        let V_guide = this.compute_volume(1.496e11);
        this.systems.set("Student's Guide to the Universe", new SystemData(
            "Student's Guide to the Universe", "General framework using solar mass and AU-scale.",
            1.989e30, 1.496e11, 0.0, 4.35e17, V_guide, 0.0, 3e4, 1e20, 7.032e22, 1e-4, -1e-4
        ));

        // Additional systems
        let V_lagoon = this.compute_volume(5.203e17);
        this.systems.set("The Lagoon Nebula", new SystemData(
            "The Lagoon Nebula", "Emission nebula with star formation.",
            1.989e34, 5e16, 0.0001, 6.312e13, 5.913e53, 0.0, 1e4, 1e22, 8.508e35, 1e-2, -1e-2
        ));

        let V_spirals = this.compute_volume(1.543e21);
        this.systems.set("Spirals and Supernovae", new SystemData(
            "Spirals and Supernovae", "Galactic spirals and supernova dynamics.",
            1.989e41, 1.543e21, 0.002, 4.35e17, V_spirals, 0.0, 2e5, 1e24, 7.487e42, 1e-6, -1e-6
        ));

        let V_ngc = this.compute_volume(1.514e16);
        this.systems.set("NGC 6302 (Butterfly Nebula)", new SystemData(
            "NGC 6302 (Butterfly Nebula)", "Planetary nebula with bipolar outflows.",
            1.989e30, 1.514e16, 0.00001, 3.156e11, 1.458e48, 0.0, 2e4, 1e21, 7.207e32, 1e-3, -1e-3
        ));

        let V_orion = this.compute_volume(1.135e17);
        this.systems.set("Orion Nebula", new SystemData(
            "Orion Nebula", "Stellar nursery near Earth.",
            3.978e33, 1.135e17, 0.00004, 3.156e13, 6.132e51, 0.0, 1e4, 1e22, 4.047e34, 1e-2, -1e-2
        ));

        // Update volumes if not computed
        for (let [name, sys] of this.systems) {
            if (sys.vars.get("V") === 0) {
                sys.vars.set("V", this.compute_volume(sys.vars.get("r")));
            }
        }

        this.installed = true;
        console.log("UQFF Module Installed: All systems initialized with defaults.");
    }

    // Update a variable (additive or set)
    updateVariable(system_name, var_name, value, is_add = false) {
        let sys = this.systems.get(system_name);
        if (sys) {
            if (is_add) {
                sys.vars.set(var_name, sys.vars.get(var_name) + value);
            } else {
                sys.vars.set(var_name, value);
            }

            // Propagate updates to dependent vars
            if (var_name === "z") {
                let zz = sys.vars.get("z");
                sys.vars.set("H_t_z", this.H0 * (0.3 * Math.pow(1 + zz, 3) + 0.7));
                sys.vars.set("one_plus_H_t", 1 + sys.vars.get("H_t_z") * sys.vars.get("t"));
            }
        }
    }

    // Add to variable
    addToVariable(system_name, var_name, value) {
        this.updateVariable(system_name, var_name, value, true);
    }

    // Subtract from variable
    subtractFromVariable(system_name, var_name, delta) {
        this.updateVariable(system_name, var_name, -delta, true);
    }

    // Get variable value
    getVariable(system_name, var_name) {
        let sys = this.systems.get(system_name);
        return sys ? sys.vars.get(var_name) : undefined;
    }

    // Print system text/description
    printSystemText(system_name) {
        let sys = this.systems.get(system_name);
        if (sys) {
            console.log("System:", sys.name);
            console.log("Description:", sys.description);
            console.log("Key Variables:");
            for (let [key, value] of sys.vars) {
                console.log(`  ${key} = ${value.toExponential(3)}`);
            }
        } else {
            console.log("System not found:", system_name);
        }
    }

    // Compute Compressed MUGE
    computeCompressedMUGE(system_name, updates = new Map()) {
        let sys = this.systems.get(system_name);
        if (!sys) return 0.0;

        // Create local vars copy
        let local_vars = new Map(sys.vars);

        // Apply updates
        for (let [key, value] of updates) {
            local_vars.set(key, value);
        }

        // Update from struct if not in map
        if (!local_vars.has("M")) local_vars.set("M", sys.M);
        if (!local_vars.has("r")) local_vars.set("r", sys.r);
        if (!local_vars.has("z")) local_vars.set("z", sys.z);
        if (!local_vars.has("t")) local_vars.set("t", sys.t);
        if (!local_vars.has("V")) local_vars.set("V", sys.V);
        if (!local_vars.has("F_env")) local_vars.set("F_env", sys.F_env);
        if (!local_vars.has("M_sun")) local_vars.set("M_sun", sys.M_sun);
        if (!local_vars.has("r_orbit")) local_vars.set("r_orbit", sys.r_orbit);

        // Recompute dependents
        let zz = local_vars.get("z");
        let tt = local_vars.get("t");
        local_vars.set("H_t_z", this.H0 * (0.3 * Math.pow(1 + zz, 3) + 0.7));
        local_vars.set("one_plus_H_t", 1 + local_vars.get("H_t_z") * tt);
        local_vars.set("B_adjust", 1 - this.B_t / this.B_crit);
        local_vars.set("one_plus_F_env", 1 + local_vars.get("F_env"));
        local_vars.set("Lambda_c2_3", this.Lambda * this.c * this.c / 3);
        local_vars.set("hbar_over_sqrt_delta", this.hbar / Math.sqrt(this.Delta_x_Delta_p));
        local_vars.set("quantum_term", local_vars.get("hbar_over_sqrt_delta") * this.integral_psi * (2 * this.pi / this.t_Hubble));
        local_vars.set("rho_V_g", this.rho_fluid * local_vars.get("V") * this.g_earth);

        let MM = local_vars.get("M");
        let rr = local_vars.get("r");
        local_vars.set("three_G_M_over_r3", 3 * this.G * MM / (rr * rr * rr));
        local_vars.set("density_pert", this.delta_rho_over_rho + local_vars.get("three_G_M_over_r3"));
        local_vars.set("M_vis_DM_pert", (MM + this.M_DM_default) * local_vars.get("density_pert"));

        // Gravity base term
        let grav_base = (this.G * MM / (rr * rr)) * local_vars.get("one_plus_H_t") *
                       local_vars.get("B_adjust") * local_vars.get("one_plus_F_env");

        if (sys.M_sun > 0) { // Planetary: add orbital
            let orb_grav = (this.G * local_vars.get("M_sun") / (local_vars.get("r_orbit") * local_vars.get("r_orbit"))) *
                          local_vars.get("one_plus_H_t");
            grav_base += orb_grav;
        }

        // Gravity modes (0 as per doc)
        let U_g_sum = 0.0;

        // Full sum
        let muge = grav_base + U_g_sum + local_vars.get("Lambda_c2_3") + local_vars.get("quantum_term") +
                  local_vars.get("rho_V_g") + local_vars.get("M_vis_DM_pert");

        console.log(`Compressed MUGE for ${system_name}: ${muge.toExponential(3)} m/s^2`);
        console.log(`  Breakdown: grav_base=${grav_base.toExponential(3)}, U_g_sum=${U_g_sum.toExponential(3)}, Lambda=${local_vars.get("Lambda_c2_3").toExponential(3)}, quantum=${local_vars.get("quantum_term").toExponential(3)}, fluid=${local_vars.get("rho_V_g").toExponential(3)}, pert=${local_vars.get("M_vis_DM_pert").toExponential(3)}`);

        return muge;
    }

    // Compute Resonance MUGE
    computeResonanceMUGE(system_name, updates = new Map()) {
        let sys = this.systems.get(system_name);
        if (!sys) return 0.0;

        // Create local vars copy
        let local_vars = new Map(sys.vars);

        // Apply updates
        for (let [key, value] of updates) {
            local_vars.set(key, value);
        }

        // Update basics
        if (!local_vars.has("M")) local_vars.set("M", sys.M);
        if (!local_vars.has("r")) local_vars.set("r", sys.r);
        if (!local_vars.has("V")) local_vars.set("V", sys.V);
        if (!local_vars.has("v_exp")) local_vars.set("v_exp", sys.v_exp);
        if (!local_vars.has("I")) local_vars.set("I", sys.I);
        if (!local_vars.has("A")) local_vars.set("A", sys.A);
        if (!local_vars.has("omega1")) local_vars.set("omega1", sys.omega1);
        if (!local_vars.has("omega2")) local_vars.set("omega2", sys.omega2);
        if (!local_vars.has("z")) local_vars.set("z", sys.z);

        let zz = local_vars.get("z");
        if (zz === 0) {
            local_vars.set("H_z", this.H0);
        } else {
            local_vars.set("H_z", this.H0 * (0.3 * Math.pow(1 + zz, 3) + 0.7));
        }

        let II = local_vars.get("I");
        let AA = local_vars.get("A");
        let delta_omega = local_vars.get("omega1") - local_vars.get("omega2");
        let F_DPM = II * AA * delta_omega;
        local_vars.set("f_DPM", 1e12); // Hz fixed
        let a_DPM = F_DPM * local_vars.get("f_DPM") * this.E_vac_neb / (this.c * local_vars.get("V"));

        // THz Hole Resonance
        let a_THz = local_vars.get("f_DPM") * this.E_vac_neb * local_vars.get("v_exp") * a_DPM / (this.E_vac_ISM * this.c);

        // Plasmotic Vacuum Energy Density Differential
        let a_vac_diff = this.Delta_E_vac * Math.pow(local_vars.get("v_exp"), 2) * a_DPM / (this.E_vac_neb * this.c * this.c);

        // Superconductor Frequency Interaction
        let a_super_freq = this.F_super * local_vars.get("f_DPM") * a_DPM / (this.E_vac_neb * this.c);

        // Aether-Mediated Resonance
        let a_aether_res = this.k_4 * this.omega_i * local_vars.get("f_DPM") * a_DPM * (1 + this.UA_SC_m * 0.1);

        // Reactive Dynamics U_g4i (0 as per doc)
        let U_g4i = 0.0;

        // Quantum Wave Dynamics
        let a_quantum_freq = this.f_quantum * this.E_vac_neb * a_DPM / (this.E_vac_ISM * this.c);

        // Aether Effect
        let a_Aether_freq = this.f_Aether * this.E_vac_neb * a_DPM / (this.E_vac_ISM * this.c);

        // Fluid Dynamics
        let f_fluid = (this.G * local_vars.get("M") / (local_vars.get("r") * local_vars.get("r"))) / (2 * this.pi);
        let a_fluid_freq = f_fluid * this.E_vac_neb * local_vars.get("V") / (this.E_vac_ISM * this.c);

        // Oscillatory Components (0)
        let Osc_term = 0.0;

        // Cosmic Expansion
        let f_exp = local_vars.get("H_z") * local_vars.get("t") / (2 * this.pi);
        let a_exp_freq = f_exp * this.E_vac_neb * a_DPM / (this.E_vac_ISM * this.c);

        // Final sum
        let muge = a_DPM + a_THz + a_vac_diff + a_super_freq + a_aether_res + U_g4i +
                  a_quantum_freq + a_Aether_freq + a_fluid_freq + Osc_term + a_exp_freq + this.f_TRZ;

        console.log(`Resonance MUGE for ${system_name}: ${muge.toExponential(3)} m/s^2`);
        console.log(`  Breakdown: a_DPM=${a_DPM.toExponential(3)}, a_THz=${a_THz.toExponential(3)}, a_vac_diff=${a_vac_diff.toExponential(3)}, a_super_freq=${a_super_freq.toExponential(3)}, a_aether_res=${a_aether_res.toExponential(3)}, U_g4i=${U_g4i.toExponential(3)}, a_quantum=${a_quantum_freq.toExponential(3)}, a_Aether=${a_Aether_freq.toExponential(3)}, a_fluid=${a_fluid_freq.toExponential(3)}, a_exp=${a_exp_freq.toExponential(3)}, f_TRZ=${this.f_TRZ.toExponential(3)}`);

        return muge;
    }

    // Example usage
    exampleUsage() {
        this.install_uqff_module();
        this.computeCompressedMUGE("Hubble Sees Galaxies Galore");
        this.computeResonanceMUGE("Hubble Sees Galaxies Galore");
        this.printSystemText("Hubble Sees Galaxies Galore");
        this.addToVariable("Hubble Sees Galaxies Galore", "M", 1e40);
        this.computeCompressedMUGE("Hubble Sees Galaxies Galore");
    }
}

// Export the module
module.exports = {
    UQFFModule: UQFFModule,
    SystemData: SystemData
};