// Source23.js - IC 1396 Elephant Trunk Nebula Module (JavaScript Stub)

class IC1396ElephantTrunk {
    constructor() {
        this.name = "IC 1396 Elephant Trunk";
        this.length = 2e17; // meters
        this.density = 100; // Particles/cm³
        this.temperature = 50; // K
        this.ionizing_flux = 1e8; // Photons/m²/s
        this.magnetic_field = 1e-4; // T
    }

    compute_g_ElephantTrunk(time) {
        const G = 6.67430e-11;

        // Gas pillar gravity
        const mass = this.density * this.length * this.length * this.length * 1.67e-27;
        const g_gas = G * mass / (this.length * this.length);

        // Photoionization effects
        const ionization_term = this.ionizing_flux * 1e-20;

        // Magnetic support
        const magnetic_term = this.magnetic_field * this.magnetic_field / (4 * Math.PI * 1e-7);

        return g_gas + ionization_term + magnetic_term;
    }

    getParameters() {
        return {
            name: this.name,
            length: this.length,
            density: this.density,
            temperature: this.temperature,
            ionizing_flux: this.ionizing_flux,
            magnetic_field: this.magnetic_field
        };
    }
}

module.exports = { IC1396ElephantTrunk };