// source5.js - JavaScript port of Source5.cpp with 2.0-Enhanced Framework
// Self-Expanding UQFF Module with PhysicsTerm Plugin Architecture

// ============================================================================
// UQFF 2.0-Enhanced Framework: PhysicsTerm Classes
// ============================================================================

class PhysicsTerm {
    compute(params) {
        throw new Error("Must implement compute() method");
    }
    
    description() {
        throw new Error("Must implement description() method");
    }
    
    version() {
        return "1.0";
    }
}

class DarkMatterHaloTerm extends PhysicsTerm {
    constructor(M_halo, r_scale) {
        super();
        this.M_halo = M_halo;      // Halo mass (kg)
        this.r_scale = r_scale;    // Scale radius (m)
    }
    
    compute(params) {
        const r = params.r || 0;
        if (r === 0 || this.r_scale === 0) return 0.0;
        
        // NFW profile contribution
        const G = 6.67430e-11;
        const x = r / this.r_scale;
        const rho_0 = this.M_halo / (4.0 * Math.PI * Math.pow(this.r_scale, 3) * (Math.log(2.0) - 0.5));
        return G * this.M_halo * Math.log(1 + x) / (r * x);
    }
    
    description() {
        return "Dark matter halo contribution (NFW profile)";
    }
}

class VacuumEnergyTerm extends PhysicsTerm {
    constructor(E_vac_scale, lambda) {
        super();
        this.E_vac_scale = E_vac_scale;  // Vacuum energy scale
        this.lambda = lambda;             // Coupling strength
    }
    
    compute(params) {
        const t = params.t || 0;
        // Time-varying vacuum energy contribution
        return this.lambda * this.E_vac_scale * (1.0 + 0.1 * Math.sin(1e-10 * t));
    }
    
    description() {
        return "Vacuum energy fluctuation term";
    }
}

// ============================================================================
// Data Structures
// ============================================================================

class CelestialBody {
    constructor(name, Ms, Rs, Rb, Ts_surface, omega_s, Bs_avg, SCm_density, QUA, Pcore, PSCm, omega_c) {
        this.name = name;
        this.Ms = Ms;                    // Mass (kg)
        this.Rs = Rs;                    // Radius (m)
        this.Rb = Rb;                    // Bubble radius (m)
        this.Ts_surface = Ts_surface;    // Surface temperature (K)
        this.omega_s = omega_s;          // Rotation rate (rad/s)
        this.Bs_avg = Bs_avg;            // Average surface magnetic field (T)
        this.SCm_density = SCm_density;  // SCm density (kg/m^3)
        this.QUA = QUA;                  // Trapped Universal Aether charge (C)
        this.Pcore = Pcore;              // Core penetration factor
        this.PSCm = PSCm;                // SCm penetration factor
        this.omega_c = omega_c;          // Cycle frequency (rad/s)
    }
}

class ResonanceParams {
    constructor() {
        this.fDPM = 1e12;
        this.fTHz = 1e12;
        this.Evac_neb = 7.09e-36;
        this.Evac_ISM = 7.09e-37;
        this.Delta_Evac = 6.381e-36;
        this.Fsuper = 6.287e-19;
        this.UA_SCM = 10;
        this.omega_i = 1e-8;
        this.k4_res = 1.0;
        this.freact = 1e10;
        this.fquantum = 1.445e-17;
        this.fAether = 1.576e-35;
        this.fosc = 4.57e14;
        this.fTRZ = 0.1;
        this.c_res = 3e8;
    }
}

class MUGESystem {
    constructor(name, I, A, omega1, omega2, Vsys, vexp, t, z, ffluid, M, r, B, Bcrit, rho_fluid, g_local, M_DM, delta_rho_rho) {
        this.name = name;
        this.I = I;
        this.A = A;
        this.omega1 = omega1;
        this.omega2 = omega2;
        this.Vsys = Vsys;
        this.vexp = vexp;
        this.t = t;
        this.z = z;
        this.ffluid = ffluid;
        this.M = M;
        this.r = r;
        this.B = B;
        this.Bcrit = Bcrit;
        this.rho_fluid = rho_fluid;
        this.g_local = g_local;
        this.M_DM = M_DM;
        this.delta_rho_rho = delta_rho_rho;
    }
}

class FluidSolver {
    constructor(N = 32) {
        this.N = N;
        const size = (N + 2) * (N + 2);
        this.u = new Array(size).fill(0);
        this.v = new Array(size).fill(0);
        this.u_prev = new Array(size).fill(0);
        this.v_prev = new Array(size).fill(0);
        this.dens = new Array(size).fill(0);
        this.dens_prev = new Array(size).fill(0);
        this.dt_ns = 0.1;
        this.visc = 0.0001;
    }
    
    IX(i, j) {
        return i + (this.N + 2) * j;
    }
    
    add_source(x, s) {
        for (let i = 0; i < x.length; i++) {
            x[i] += this.dt_ns * s[i];
        }
    }
    
    set_bnd(b, x) {
        for (let i = 1; i <= this.N; i++) {
            x[this.IX(0, i)] = (b === 1) ? -x[this.IX(1, i)] : x[this.IX(1, i)];
            x[this.IX(this.N + 1, i)] = (b === 1) ? -x[this.IX(this.N, i)] : x[this.IX(this.N, i)];
            x[this.IX(i, 0)] = (b === 2) ? -x[this.IX(i, 1)] : x[this.IX(i, 1)];
            x[this.IX(i, this.N + 1)] = (b === 2) ? -x[this.IX(i, this.N)] : x[this.IX(i, this.N)];
        }
        x[this.IX(0, 0)] = 0.5 * (x[this.IX(1, 0)] + x[this.IX(0, 1)]);
        x[this.IX(0, this.N + 1)] = 0.5 * (x[this.IX(1, this.N + 1)] + x[this.IX(0, this.N)]);
        x[this.IX(this.N + 1, 0)] = 0.5 * (x[this.IX(this.N, 0)] + x[this.IX(this.N + 1, 1)]);
        x[this.IX(this.N + 1, this.N + 1)] = 0.5 * (x[this.IX(this.N, this.N + 1)] + x[this.IX(this.N + 1, this.N)]);
    }
    
    diffuse(b, x, x0, diff) {
        const a = this.dt_ns * diff * this.N * this.N;
        for (let k = 0; k < 20; k++) {
            for (let i = 1; i <= this.N; i++) {
                for (let j = 1; j <= this.N; j++) {
                    x[this.IX(i, j)] = (x0[this.IX(i, j)] + a * (x[this.IX(i - 1, j)] + x[this.IX(i + 1, j)] +
                        x[this.IX(i, j - 1)] + x[this.IX(i, j + 1)])) / (1 + 4 * a);
                }
            }
            this.set_bnd(b, x);
        }
    }
    
    advect(b, d, d0) {
        for (let i = 1; i <= this.N; i++) {
            for (let j = 1; j <= this.N; j++) {
                let x = i - this.dt_ns * this.N * this.u[this.IX(i, j)];
                let y = j - this.dt_ns * this.N * this.v[this.IX(i, j)];
                if (x < 0.5) x = 0.5; if (x > this.N + 0.5) x = this.N + 0.5;
                if (y < 0.5) y = 0.5; if (y > this.N + 0.5) y = this.N + 0.5;
                const i0 = Math.floor(x), i1 = i0 + 1;
                const j0 = Math.floor(y), j1 = j0 + 1;
                const s1 = x - i0, s0 = 1 - s1;
                const t1 = y - j0, t0 = 1 - t1;
                d[this.IX(i, j)] = s0 * (t0 * d0[this.IX(i0, j0)] + t1 * d0[this.IX(i0, j1)]) +
                    s1 * (t0 * d0[this.IX(i1, j0)] + t1 * d0[this.IX(i1, j1)]);
            }
        }
        this.set_bnd(b, d);
    }
    
    project(u, v, p, div) {
        const h = 1.0 / this.N;
        for (let i = 1; i <= this.N; i++) {
            for (let j = 1; j <= this.N; j++) {
                div[this.IX(i, j)] = -0.5 * h * (u[this.IX(i + 1, j)] - u[this.IX(i - 1, j)] + v[this.IX(i, j + 1)] - v[this.IX(i, j - 1)]);
                p[this.IX(i, j)] = 0;
            }
        }
        this.set_bnd(0, div); this.set_bnd(0, p);
        for (let k = 0; k < 20; k++) {
            for (let i = 1; i <= this.N; i++) {
                for (let j = 1; j <= this.N; j++) {
                    p[this.IX(i, j)] = (div[this.IX(i, j)] + p[this.IX(i - 1, j)] + p[this.IX(i + 1, j)] +
                        p[this.IX(i, j - 1)] + p[this.IX(i, j + 1)]) / 4;
                }
            }
            this.set_bnd(0, p);
        }
        for (let i = 1; i <= this.N; i++) {
            for (let j = 1; j <= this.N; j++) {
                u[this.IX(i, j)] -= 0.5 * (p[this.IX(i + 1, j)] - p[this.IX(i - 1, j)]) / h;
                v[this.IX(i, j)] -= 0.5 * (p[this.IX(i, j + 1)] - p[this.IX(i, j - 1)]) / h;
            }
        }
        this.set_bnd(1, u); this.set_bnd(2, v);
    }
    
    step(uqff_g = 0.0) {
        for (let i = 1; i <= this.N; i++) {
            for (let j = 1; j <= this.N; j++) {
                this.v[this.IX(i, j)] += this.dt_ns * uqff_g;
            }
        }
        this.diffuse(1, this.u_prev, this.u, this.visc);
        this.diffuse(2, this.v_prev, this.v, this.visc);
        this.project(this.u_prev, this.v_prev, this.u, this.v);
        this.advect(1, this.u, this.u_prev);
        this.advect(2, this.v, this.v_prev);
        this.project(this.u, this.v, this.u_prev, this.v_prev);
    }
    
    add_jet_force(force) {
        for (let i = Math.floor(this.N / 4); i <= Math.floor(3 * this.N / 4); i++) {
            this.v[this.IX(i, Math.floor(this.N / 2))] += force;
        }
    }
    
    print_velocity_field() {
        console.log("Velocity field (magnitude):");
        for (let j = this.N; j >= 1; j--) {
            let line = "";
            for (let i = 1; i <= this.N; i++) {
                const mag = Math.sqrt(this.u[this.IX(i, j)] ** 2 + this.v[this.IX(i, j)] ** 2);
                line += (mag > 1.0) ? '#' : (mag > 0.5) ? '+' : (mag > 0.1) ? '.' : ' ';
            }
            console.log(line);
        }
    }
}

// ============================================================================
// Constants
// ============================================================================

const PI = Math.PI;
const c = 3.0e8;
const G = 6.67430e-11;

let v_SCm = 0.99 * c;
let rho_A = 1e-23;
let rho_sw = 8e-21;
let v_sw = 5e5;
let QA = 1e-10;
let Qs = 0.0;
let kappa = 0.0005;
let alpha = 0.001;
let gamma = 0.00005;
let delta_sw = 0.01;
let epsilon_sw = 0.001;
let delta_def = 0.01;
let HSCm = 1.0;
let UUA = 1.0;
let eta = 1e-22;
let k1 = 1.5, k2 = 1.2, k3 = 1.8, k4 = 2.0;
let beta_i = 0.6;
let rho_v = 6e-27;
let C_concentration = 1.0;
let f_feedback = 0.1;
const num_strings = 1e9;
let Ts00 = 1.27e3 + 1.11e7;

// ============================================================================
// Helper Functions
// ============================================================================

function step_function(r, Rb) {
    return (r > Rb) ? 1.0 : 0.0;
}

function compute_Ereact(t, rho_SCm, v_SCm, rho_A, kappa) {
    if (rho_A === 0.0) throw new Error("Division by zero in rho_A");
    return (rho_SCm * v_SCm * v_SCm / rho_A) * Math.exp(-kappa * t);
}

function compute_mu_s(t, Bs, omega_c, Rs, SCm_contrib = 1e3) {
    const Bs_t = Bs + 0.4 * Math.sin(omega_c * t) + SCm_contrib;
    return Bs_t * Math.pow(Rs, 3);
}

function compute_grad_Ms_r(Ms, Rs) {
    if (Rs === 0.0) throw new Error("Division by zero in Rs");
    return G * Ms / (Rs * Rs);
}

function compute_Bj(t, omega_c, SCm_contrib = 1e3) {
    return 1e-3 + 0.4 * Math.sin(omega_c * t) + SCm_contrib;
}

function compute_omega_s_t(t, omega_s, omega_c) {
    return omega_s - 0.4e-6 * Math.sin(omega_c * t);
}

function compute_mu_j(t, omega_c, Rs, SCm_contrib = 1e3) {
    const Bj = compute_Bj(t, omega_c, SCm_contrib);
    return Bj * Math.pow(Rs, 3);
}

// ============================================================================
// Main UQFF Functions
// ============================================================================

function compute_Ug1(body, r, t, tn, alpha, delta_def, k1) {
    if (r <= 0.0) throw new Error("Invalid r value");
    const mu_s = compute_mu_s(t, body.Bs_avg, body.omega_c, body.Rs);
    const grad_Ms_r = compute_grad_Ms_r(body.Ms, body.Rs);
    const defect = 1.0 + delta_def * Math.sin(0.001 * t);
    return k1 * mu_s * grad_Ms_r * Math.exp(-alpha * t) * Math.cos(PI * tn) * defect;
}

function compute_Ug2(body, r, t, tn, k2, QA, delta_sw, v_sw, HSCm, rho_A, kappa) {
    if (r === 0.0) throw new Error("Division by zero in r");
    const Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A, kappa);
    const S = step_function(r, body.Rb);
    const wind_mod = 1.0 + delta_sw * v_sw;
    return k2 * (QA + body.QUA) * body.Ms / (r * r) * S * wind_mod * HSCm * Ereact;
}

function compute_Ug3(body, r, t, tn, theta, rho_A, kappa, k3) {
    const Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A, kappa);
    const omega_s_t = compute_omega_s_t(t, body.omega_s, body.omega_c);
    const Bj = compute_Bj(t, body.omega_c);
    return k3 * Bj * Math.cos(omega_s_t * t * PI) * body.Pcore * Ereact;
}

function compute_Ug4(t, tn, rho_v, C_concentration, Mbh, dg, alpha, f_feedback, k4) {
    const decay = Math.exp(-alpha * t);
    const cycle = Math.cos(PI * tn);
    return k4 * rho_v * C_concentration * Mbh / dg * decay * cycle * (1 + f_feedback);
}

function compute_Um(body, t, tn, rj, gamma, rho_A, kappa, num_strings, phi_hat = 1.0) {
    if (rj === 0.0) throw new Error("Division by zero in rj");
    const Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A, kappa);
    const mu_j = compute_mu_j(t, body.omega_c, body.Rs);
    const decay = 1.0 - Math.exp(-gamma * t * Math.cos(PI * tn));
    const single = mu_j / rj * decay * phi_hat;
    return single * num_strings * body.PSCm * Ereact;
}

// ============================================================================
// MUGE Functions
// ============================================================================

function compute_compressed_base(sys) {
    if (sys.r === 0.0) throw new Error("Division by zero in r");
    return G * sys.M / (sys.r * sys.r);
}

function compute_compressed_expansion(sys, H0 = 2.269e-18) {
    const H_tz = H0 * sys.t;
    return 1 + H_tz;
}

function compute_compressed_super_adj(sys) {
    if (sys.Bcrit === 0.0) throw new Error("Division by zero in Bcrit");
    return 1 - sys.B / sys.Bcrit;
}

function compute_compressed_env() {
    return 1.0;
}

function compute_compressed_Ug_sum() {
    return 0.0;
}

function compute_compressed_cosm(Lambda = 1.1e-52) {
    return Lambda * c * c / 3.0;
}

function compute_compressed_quantum(hbar = 1.0546e-34, Delta_x_p = 1e-68, integral_psi = 2.176e-18, tHubble = 4.35e17) {
    if (Delta_x_p === 0.0) throw new Error("Division by zero in Delta_x_p");
    return (hbar / Delta_x_p) * integral_psi * (2 * PI / tHubble);
}

function compute_compressed_fluid(sys) {
    return sys.rho_fluid * sys.Vsys * sys.g_local;
}

function compute_compressed_perturbation(sys) {
    if (sys.r === 0.0) throw new Error("Division by zero in r^3");
    return (sys.M + sys.M_DM) * (sys.delta_rho_rho + 3 * G * sys.M / (sys.r * sys.r * sys.r));
}

function compute_compressed_MUGE(sys) {
    const base = compute_compressed_base(sys);
    const expansion = compute_compressed_expansion(sys);
    const super_adj = compute_compressed_super_adj(sys);
    const env = compute_compressed_env();
    const adjusted_base = base * expansion * super_adj * env;
    
    const Ug_sum = compute_compressed_Ug_sum();
    const cosm = compute_compressed_cosm();
    const quantum = compute_compressed_quantum();
    const fluid = compute_compressed_fluid(sys);
    const perturbation = compute_compressed_perturbation(sys);
    
    return adjusted_base + Ug_sum + cosm + quantum + fluid + perturbation;
}

function compute_aDPM(sys, res) {
    const FDPM = sys.I * sys.A * (sys.omega1 - sys.omega2);
    return FDPM * res.fDPM * res.Evac_neb * res.c_res * sys.Vsys;
}

function compute_aTHz(aDPM, sys, res) {
    return res.fTHz * res.Evac_neb * sys.vexp * aDPM / res.Evac_ISM / res.c_res;
}

function compute_avac_diff(aDPM, sys, res) {
    return res.Delta_Evac * sys.vexp * sys.vexp * aDPM / res.Evac_neb / (res.c_res * res.c_res);
}

function compute_asuper_freq(aDPM, res) {
    return res.Fsuper * res.fTHz * aDPM / res.Evac_neb / res.c_res;
}

function compute_aaether_res(aDPM, res) {
    return res.UA_SCM * res.omega_i * res.fTHz * aDPM * (1 + res.fTRZ);
}

function compute_Ug4i(aDPM, sys, res) {
    const Ereact = 1046 * Math.exp(-0.0005 * sys.t);
    return res.k4_res * Ereact * res.freact * aDPM / res.Evac_neb * res.c_res;
}

function compute_aquantum_freq(aDPM, res) {
    return res.fquantum * res.Evac_neb * aDPM / res.Evac_ISM / res.c_res;
}

function compute_aAether_freq(aDPM, res) {
    return res.fAether * res.Evac_neb * aDPM / res.Evac_ISM / res.c_res;
}

function compute_afluid_freq(sys, res) {
    return sys.ffluid * res.Evac_neb * sys.Vsys / res.Evac_ISM / res.c_res;
}

function compute_Osc_term() {
    return 0.0;
}

function compute_aexp_freq(aDPM, sys, res, H_z = 2.270e-18) {
    const fexp = 2 * PI * H_z * sys.t;
    return fexp * res.Evac_neb * aDPM / res.Evac_ISM / res.c_res;
}

function compute_fTRZ(res) {
    return res.fTRZ;
}

function compute_a_wormhole(r, b = 1.0, f_worm = 1.0, Evac_neb = 7.09e-36) {
    return f_worm * Evac_neb * (1.0 / (b * b + r * r));
}

function compute_resonance_MUGE(sys, res) {
    const aDPM = compute_aDPM(sys, res);
    const aTHz = compute_aTHz(aDPM, sys, res);
    const avac_diff = compute_avac_diff(aDPM, sys, res);
    const asuper_freq = compute_asuper_freq(aDPM, res);
    const aaether_res = compute_aaether_res(aDPM, res);
    const Ug4i = compute_Ug4i(aDPM, sys, res);
    const aquantum_freq = compute_aquantum_freq(aDPM, res);
    const aAether_freq = compute_aAether_freq(aDPM, res);
    const afluid_freq = compute_afluid_freq(sys, res);
    const Osc_term = compute_Osc_term();
    const aexp_freq = compute_aexp_freq(aDPM, sys, res);
    const fTRZ = compute_fTRZ(res);
    const a_worm = compute_a_wormhole(sys.r);
    
    return aDPM + aTHz + avac_diff + asuper_freq + aaether_res + Ug4i + aquantum_freq + aAether_freq + afluid_freq + Osc_term + aexp_freq + fTRZ + a_worm;
}

function compute_Ubi(Ugi, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn) {
    const wind_mod = 1.0 + epsilon_sw * rho_sw;
    return -beta_i * Ugi * Omega_g * Mbh / dg * wind_mod * UUA * Math.cos(PI * tn);
}

function compute_A_mu_nu(tn, eta, Ts00) {
    const mod = eta * Ts00 * Math.cos(PI * tn);
    return [
        [1.0 + mod, 0.0 + mod, 0.0 + mod, 0.0 + mod],
        [0.0 + mod, -1.0 + mod, 0.0 + mod, 0.0 + mod],
        [0.0 + mod, 0.0 + mod, -1.0 + mod, 0.0 + mod],
        [0.0 + mod, 0.0 + mod, 0.0 + mod, -1.0 + mod]
    ];
}

function compute_FU(body, r, t, tn, theta, Omega_g = 7.3e-16, Mbh = 8.15e36, dg = 2.55e20) {
    const Ug1 = compute_Ug1(body, r, t, tn, alpha, delta_def, k1);
    const Ug2 = compute_Ug2(body, r, t, tn, k2, QA, delta_sw, v_sw, HSCm, rho_A, kappa);
    const Ug3 = compute_Ug3(body, r, t, tn, theta, rho_A, kappa, k3);
    const Ug4 = compute_Ug4(t, tn, rho_v, C_concentration, Mbh, dg, alpha, f_feedback, k4);
    const sum_Ugi = Ug1 + Ug2 + Ug3 + Ug4;
    
    const Ubi1 = compute_Ubi(Ug1, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
    const Ubi2 = compute_Ubi(Ug2, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
    const Ubi3 = compute_Ubi(Ug3, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
    const Ubi4 = compute_Ubi(Ug4, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
    const sum_Ubi = Ubi1 + Ubi2 + Ubi3 + Ubi4;
    
    const Um = compute_Um(body, t, tn, body.Rb, gamma, rho_A, kappa, num_strings);
    
    const A = compute_A_mu_nu(tn, eta, Ts00);
    const A_scalar = A[0][0] + A[1][1] + A[2][2] + A[3][3];
    
    return sum_Ugi + sum_Ubi + Um + A_scalar;
}

function simulate_quasar_jet(initial_velocity) {
    const solver = new FluidSolver();
    solver.add_jet_force(initial_velocity / 10.0);
    
    const res = new ResonanceParams();
    const sagA = new MUGESystem("Sagittarius A*", 1e23, 2.813e30, 1e-5, -1e-5, 3.552e45, 5e6, 3.786e14, 0.0009, 3.465e-8, 8.155e36, 1e12, 1e-5, 1e-4, 1e-20, 1e-5, 1e37, 1e-3);
    const uqff_g = compute_resonance_MUGE(sagA, res);
    
    console.log(`Simulating quasar jet with Navier-Stokes (10 steps) using UQFF g=${uqff_g}...`);
    for (let step = 0; step < 10; step++) {
        solver.step(uqff_g / 1e30);
    }
    solver.print_velocity_field();
}

// ============================================================================
// UQFFModule5JS: Self-Expanding Module
// ============================================================================

class UQFFModule5JS {
    constructor() {
        this.dynamic_terms = [];
        this.dynamic_parameters = new Map();
        this.metadata = new Map([
            ["version", "2.0-Enhanced"],
            ["created", "2025-11-08"],
            ["framework", "Self-Expanding UQFF"]
        ]);
        this.learning_rate = 0.01;
        this.logging_enabled = false;
    }
    
    log(message) {
        if (this.logging_enabled) {
            console.log(`[UQFFModule5JS] ${message}`);
        }
    }
    
    registerDynamicTerm(term) {
        this.log(`Registering dynamic term: ${term.description()}`);
        this.dynamic_terms.push(term);
    }
    
    setDynamicParameter(name, value) {
        this.log(`Setting parameter: ${name} = ${value}`);
        this.dynamic_parameters.set(name, value);
    }
    
    getDynamicParameter(name, default_val = 0.0) {
        return this.dynamic_parameters.has(name) ? this.dynamic_parameters.get(name) : default_val;
    }
    
    computeDynamicContributions(params) {
        let sum = 0.0;
        for (const term of this.dynamic_terms) {
            sum += term.compute(params);
        }
        return sum;
    }
    
    exportState(filename) {
        const fs = require('fs');
        let output = "# UQFFModule5JS State Export\n";
        output += `# Version: ${this.metadata.get("version")}\n`;
        output += `# Created: ${this.metadata.get("created")}\n\n`;
        
        output += "[Parameters]\n";
        for (const [key, value] of this.dynamic_parameters) {
            output += `${key} = ${value}\n`;
        }
        
        output += "\n[Terms]\n";
        for (let i = 0; i < this.dynamic_terms.length; i++) {
            output += `Term_${i} = ${this.dynamic_terms[i].description()}\n`;
        }
        
        output += "\n[Metadata]\n";
        for (const [key, value] of this.metadata) {
            output += `${key} = ${value}\n`;
        }
        
        fs.writeFileSync(filename, output);
        this.log(`State exported to: ${filename}`);
    }
    
    setLearningRate(rate) {
        this.learning_rate = rate;
        this.log(`Learning rate set to: ${rate}`);
    }
    
    setEnableLogging(enable) {
        this.logging_enabled = enable;
    }
    
    printInfo() {
        console.log("=== UQFFModule5JS Info ===");
        console.log(`Version: ${this.metadata.get("version")}`);
        console.log(`Dynamic Terms: ${this.dynamic_terms.length}`);
        console.log(`Dynamic Parameters: ${this.dynamic_parameters.size}`);
        console.log(`Learning Rate: ${this.learning_rate}`);
        console.log(`Logging: ${this.logging_enabled ? "Enabled" : "Disabled"}`);
    }
    
    compute_Ug1_enhanced(body, r, t, tn, alpha, delta_def, k1) {
        const original = compute_Ug1(body, r, t, tn, alpha, delta_def, k1);
        const params = {r, t, tn};
        const dynamic = this.computeDynamicContributions(params);
        return original + dynamic;
    }
    
    compute_Ug2_enhanced(body, r, t, tn, k2, QA, delta_sw, v_sw, HSCm, rho_A, kappa) {
        const original = compute_Ug2(body, r, t, tn, k2, QA, delta_sw, v_sw, HSCm, rho_A, kappa);
        const params = {r, t, tn};
        const dynamic = this.computeDynamicContributions(params);
        return original + dynamic;
    }
    
    compute_Ug3_enhanced(body, r, t, tn, theta, rho_A, kappa, k3) {
        const original = compute_Ug3(body, r, t, tn, theta, rho_A, kappa, k3);
        const params = {r, t, tn, theta};
        const dynamic = this.computeDynamicContributions(params);
        return original + dynamic;
    }
    
    compute_MUGE_enhanced(sys, res, use_compressed = true) {
        const original = use_compressed ? compute_compressed_MUGE(sys) : compute_resonance_MUGE(sys, res);
        const params = {t: sys.t, r: sys.r, M: sys.M};
        const dynamic = this.computeDynamicContributions(params);
        return original + dynamic;
    }
}

// ============================================================================
// Exports
// ============================================================================

module.exports = {
    PhysicsTerm,
    DarkMatterHaloTerm,
    VacuumEnergyTerm,
    UQFFModule5JS,
    CelestialBody,
    ResonanceParams,
    MUGESystem,
    FluidSolver,
    compute_Ug1,
    compute_Ug2,
    compute_Ug3,
    compute_Ug4,
    compute_Um,
    compute_compressed_base,
    compute_compressed_expansion,
    compute_compressed_super_adj,
    compute_compressed_env,
    compute_compressed_Ug_sum,
    compute_compressed_cosm,
    compute_compressed_quantum,
    compute_compressed_fluid,
    compute_compressed_perturbation,
    compute_compressed_MUGE,
    compute_aDPM,
    compute_aTHz,
    compute_avac_diff,
    compute_asuper_freq,
    compute_aaether_res,
    compute_Ug4i,
    compute_aquantum_freq,
    compute_aAether_freq,
    compute_afluid_freq,
    compute_Osc_term,
    compute_aexp_freq,
    compute_fTRZ,
    compute_a_wormhole,
    compute_resonance_MUGE,
    compute_Ubi,
    compute_A_mu_nu,
    compute_FU,
    simulate_quasar_jet
};
