// Source22.js - SN 1987A Supernova Module (JavaScript Stub)

class SN1987A {
    constructor() {
        this.name = "SN 1987A";
        this.progenitor_mass = 20 * 1.989e30; // Solar masses
        this.ejecta_mass = 15 * 1.989e30; // Solar masses
        this.explosion_energy = 1e44; // Joules
        this.distance = 1.5e5 * 3.086e16; // Light years to meters
        this.age = (2025 - 1987) * 365.25 * 24 * 3600; // seconds since explosion
    }

    compute_g_Supernova(time) {
        const G = 6.67430e-11;

        // Ejecta gravity
        const g_ejecta = G * this.ejecta_mass / (1e17 * 1e17);

        // Shock wave dynamics
        const shock_term = this.explosion_energy / (4 * Math.PI * 1e17 * 1e17);

        // Neutron star remnant
        const ns_mass = 1.4 * 1.989e30;
        const g_ns = G * ns_mass / (1e4 * 1e4);

        return g_ejecta + shock_term + g_ns;
    }

    getParameters() {
        return {
            name: this.name,
            progenitor_mass: this.progenitor_mass,
            ejecta_mass: this.ejecta_mass,
            explosion_energy: this.explosion_energy,
            distance: this.distance,
            age: this.age
        };
    }
}

module.exports = { SN1987A };