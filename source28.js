// Source28.js - NGC 602 Young Cluster SMC Module (JavaScript Stub)

class NGC602 {
    constructor() {
        this.name = "NGC 602";
        this.stars_count = 1000;
        this.age = 5e6; // years
        this.metallicity = 0.2; // Solar units
        this.distance = 2e5 * 3.086e16; // Light years
    }

    compute_g_YoungCluster(time) {
        const G = 6.67430e-11;
        const mass = this.stars_count * 0.5 * 1.989e30; // Average 0.5 M_sun
        return G * mass / (1e17 * 1e17) + 1e-10;
    }

    getParameters() {
        return {
            name: this.name,
            stars_count: this.stars_count,
            age: this.age,
            metallicity: this.metallicity,
            distance: this.distance
        };
    }
}

module.exports = { NGC602 };