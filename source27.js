// Source27.js - RS Puppis Cepheid Variable Module (JavaScript Stub)

class RSPuppis {
    constructor() {
        this.name = "RS Puppis";
        this.mass = 9 * 1.989e30; // Solar masses
        this.radius = 200 * 6.96e8; // Solar radii
        this.period = 41.4 * 24 * 3600; // days to seconds
        this.pulsation_amplitude = 0.1; // Magnitude
    }

    compute_g_VariableStar(time) {
        const G = 6.67430e-11;
        const g_base = G * this.mass / (this.radius * this.radius);
        return g_base + Math.sin(2 * Math.PI * time / this.period) * 1e-8;
    }

    getParameters() {
        return {
            name: this.name,
            mass: this.mass,
            radius: this.radius,
            period: this.period,
            pulsation_amplitude: this.pulsation_amplitude
        };
    }
}

module.exports = { RSPuppis };