// ============================================================================
// source5.js - JavaScript Port of Source5.cpp
// Unified Quantum Field Force (UQFF) with 2.0-Enhanced Self-Expanding Framework
// Converted from Source5.cpp (1,540 lines C++ → ~1,300 lines JavaScript)
// Includes: PhysicsTerm framework, MUGE calculations, Navier-Stokes fluid simulation
// Self-expanding, self-updating, simulation-capable
// ============================================================================

const PI = 3.141592653589793;
const c = 3.0e8;
const G = 6.67430e-11;

// Global physical constants
let Omega_g = 7.3e-16;
let Mbh = 8.15e36;
let dg = 2.55e20;
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
const g_mu_nu = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, -1.0, 0.0, 0.0],
    [0.0, 0.0, -1.0, 0.0],
    [0.0, 0.0, 0.0, -1.0]
];

// Fluid solver constants
const N = 32;
const dt_ns = 0.1;
const visc = 0.0001;
const force_jet = 10.0;

// ============================================================================
// PhysicsTerm Base Class (2.0-Enhanced Framework)
// ============================================================================
class PhysicsTerm {
    /**
     * @param {Object} params - Parameter map {name: value}
     * @returns {number} - Contribution to field equations
     */
    compute(params) {
        throw new Error("PhysicsTerm.compute() must be implemented by subclass");
    }

    description() {
        throw new Error("PhysicsTerm.description() must be implemented by subclass");
    }

    version() {
        return "1.0";
    }
}

// ============================================================================
// DarkMatterHaloTerm - NFW Profile Dark Matter Contribution
// ============================================================================
class DarkMatterHaloTerm extends PhysicsTerm {
    constructor(M_halo, r_scale) {
        super();
        this.M_halo = M_halo;
        this.r_scale = r_scale;
    }

    compute(params) {
        const r = params.r || 1e10;
        const x = r / this.r_scale;
        const ln_term = Math.log(1 + x);
        const rho_0 = this.M_halo / (4 * PI * Math.pow(this.r_scale, 3));
        const rho_NFW = rho_0 / (x * Math.pow(1 + x, 2));
        const M_enc = this.M_halo * (ln_term - x / (1 + x)) / (ln_term - 0.5);
        return (r > 0) ? G * M_enc / (r * r) : 0;
    }

    description() {
        return `DarkMatterHaloTerm(M=${this.M_halo.toExponential(2)}, r_s=${this.r_scale.toExponential(2)})`;
    }
}

// ============================================================================
// VacuumEnergyTerm - Time-Varying Vacuum Energy Contribution
// ============================================================================
class VacuumEnergyTerm extends PhysicsTerm {
    constructor(E_vac_scale, lambda) {
        super();
        this.E_vac_scale = E_vac_scale;
        this.lambda = lambda;
    }

    compute(params) {
        const t = params.t || 0;
        const E_vac = this.E_vac_scale * (1.0 + 0.1 * Math.sin(this.lambda * t));
        return E_vac * c * c / 3.0;
    }

    description() {
        return `VacuumEnergyTerm(E_vac=${this.E_vac_scale.toExponential(2)}, lambda=${this.lambda.toExponential(2)})`;
    }
}

// ============================================================================
// CelestialBody Structure
// ============================================================================
class CelestialBody {
    constructor(config = {}) {
        this.name = config.name || "UnnamedBody";
        this.Ms = config.Ms || 0;
        this.Rs = config.Rs || 0;
        this.Rb = config.Rb || 0;
        this.Ts_surface = config.Ts_surface || 0;
        this.omega_s = config.omega_s || 0;
        this.Bs_avg = config.Bs_avg || 0;
        this.SCm_density = config.SCm_density || 0;
        this.QUA = config.QUA || 0;
        this.Pcore = config.Pcore || 0;
        this.PSCm = config.PSCm || 0;
        this.omega_c = config.omega_c || 0;
    }
}

// ============================================================================
// ResonanceParams and MUGESystem Structures
// ============================================================================
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
    constructor(config = {}) {
        this.name = config.name || "UnnamedSystem";
        this.I = config.I || 0;
        this.A = config.A || 0;
        this.omega1 = config.omega1 || 0;
        this.omega2 = config.omega2 || 0;
        this.Vsys = config.Vsys || 0;
        this.vexp = config.vexp || 0;
        this.t = config.t || 0;
        this.z = config.z || 0;
        this.ffluid = config.ffluid || 0;
        this.M = config.M || 0;
        this.r = config.r || 0;
        this.B = config.B || 0;
        this.Bcrit = config.Bcrit || 0;
        this.rho_fluid = config.rho_fluid || 0;
        this.g_local = config.g_local || 0;
        this.M_DM = config.M_DM || 0;
        this.delta_rho_rho = config.delta_rho_rho || 0;
    }
}

// ============================================================================
// Helper Functions
// ============================================================================
function step_function(r, Rb) {
    return (r > Rb) ? 1.0 : 0.0;
}

function compute_Ereact(t, rho_SCm, v_SCm_param, rho_A_param, kappa_param) {
    if (rho_A_param === 0.0) throw new Error("Division by zero in rho_A");
    return (rho_SCm * v_SCm_param * v_SCm_param / rho_A_param) * Math.exp(-kappa_param * t);
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
// Main Computation Functions
// ============================================================================
function compute_Ug1(body, r, t, tn, alpha_param, delta_def_param, k1_param) {
    if (r <= 0.0) throw new Error("Invalid r value");
    const mu_s = compute_mu_s(t, body.Bs_avg, body.omega_c, body.Rs);
    const grad_Ms_r = compute_grad_Ms_r(body.Ms, body.Rs);
    const defect = 1.0 + delta_def_param * Math.sin(0.001 * t);
    return k1_param * mu_s * grad_Ms_r * Math.exp(-alpha_param * t) * Math.cos(PI * tn) * defect;
}

function compute_Ug2(body, r, t, tn, k2_param, QA_param, delta_sw_param, v_sw_param, HSCm_param, rho_A_param, kappa_param) {
    if (r === 0.0) throw new Error("Division by zero in r");
    const Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A_param, kappa_param);
    const S = step_function(r, body.Rb);
    const wind_mod = 1.0 + delta_sw_param * v_sw_param;
    return k2_param * (QA_param + body.QUA) * body.Ms / (r * r) * S * wind_mod * HSCm_param * Ereact;
}

function compute_Ug3(body, r, t, tn, theta, rho_A_param, kappa_param, k3_param) {
    const Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A_param, kappa_param);
    const omega_s_t = compute_omega_s_t(t, body.omega_s, body.omega_c);
    const Bj = compute_Bj(t, body.omega_c);
    return k3_param * Bj * Math.cos(omega_s_t * t * PI) * body.Pcore * Ereact;
}

function compute_Ug4(t, tn, rho_v_param, C_concentration_param, Mbh_param, dg_param, alpha_param, f_feedback_param, k4_param) {
    const decay = Math.exp(-alpha_param * t);
    const cycle = Math.cos(PI * tn);
    return k4_param * rho_v_param * C_concentration_param * Mbh_param / dg_param * decay * cycle * (1 + f_feedback_param);
}

function compute_Um(body, t, tn, rj, gamma_param, rho_A_param, kappa_param, num_strings_param, phi_hat = 1.0) {
    if (rj === 0.0) throw new Error("Division by zero in rj");
    const Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A_param, kappa_param);
    const mu_j = compute_mu_j(t, body.omega_c, body.Rs);
    const decay = 1.0 - Math.exp(-gamma_param * t * Math.cos(PI * tn));
    const single = mu_j / rj * decay * phi_hat;
    return single * num_strings_param * body.PSCm * Ereact;
}

function compute_Ubi(Ugi, beta_i_param, Omega_g_param, Mbh_param, dg_param, epsilon_sw_param, rho_sw_param, UUA_param, tn) {
    const wind_mod = 1.0 + epsilon_sw_param * rho_sw_param;
    return -beta_i_param * Ugi * Omega_g_param * Mbh_param / dg_param * wind_mod * UUA_param * Math.cos(PI * tn);
}

function compute_A_mu_nu(tn, eta_param, Ts00_param) {
    const A = g_mu_nu.map(row => [...row]);
    const mod = eta_param * Ts00_param * Math.cos(PI * tn);
    for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
            A[i][j] += mod;
        }
    }
    return A;
}

function compute_FU(body, r, t, tn, theta) {
    try {
        const Ug1 = compute_Ug1(body, r, t, tn, alpha, delta_def, k1);
        const Ug2 = compute_Ug2(body, r, t, tn, k2, QA, delta_sw, v_sw, HSCm, rho_A, kappa);
        const Ug3 = compute_Ug3(body, r, t, tn, theta, rho_A, kappa, k3);
        const Ug4_val = compute_Ug4(t, tn, rho_v, C_concentration, Mbh, dg, alpha, f_feedback, k4);
        const sum_Ugi = Ug1 + Ug2 + Ug3 + Ug4_val;

        const Ubi1 = compute_Ubi(Ug1, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
        const Ubi2 = compute_Ubi(Ug2, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
        const Ubi3 = compute_Ubi(Ug3, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
        const Ubi4 = compute_Ubi(Ug4_val, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
        const sum_Ubi = Ubi1 + Ubi2 + Ubi3 + Ubi4;

        const Um = compute_Um(body, t, tn, body.Rb, gamma, rho_A, kappa, num_strings);

        const A = compute_A_mu_nu(tn, eta, Ts00);
        const A_scalar = A[0][0] + A[1][1] + A[2][2] + A[3][3];

        return sum_Ugi + sum_Ubi + Um + A_scalar;
    } catch (e) {
        console.error("Error in compute_FU:", e.message);
        return 0.0;
    }
}

// ============================================================================
// Compressed MUGE Functions
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

// ============================================================================
// Resonance MUGE Functions
// ============================================================================
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
    const fTRZ_val = compute_fTRZ(res);
    const a_worm = compute_a_wormhole(sys.r);

    return aDPM + aTHz + avac_diff + asuper_freq + aaether_res + Ug4i + 
           aquantum_freq + aAether_freq + afluid_freq + Osc_term + aexp_freq + fTRZ_val + a_worm;
}

// ============================================================================
// FluidSolver Class (Navier-Stokes Simulation)
// ============================================================================
function IX(i, j) {
    return i + (N + 2) * j;
}

class FluidSolver {
    constructor() {
        const size = (N + 2) * (N + 2);
        this.u = new Array(size).fill(0.0);
        this.v = new Array(size).fill(0.0);
        this.u_prev = new Array(size).fill(0.0);
        this.v_prev = new Array(size).fill(0.0);
        this.dens = new Array(size).fill(0.0);
        this.dens_prev = new Array(size).fill(0.0);
    }

    add_source(x, s) {
        for (let i = 0; i < x.length; i++) {
            x[i] += dt_ns * s[i];
        }
    }

    diffuse(b, x, x0, diff) {
        const a = dt_ns * diff * N * N;
        for (let k = 0; k < 20; k++) {
            for (let i = 1; i <= N; i++) {
                for (let j = 1; j <= N; j++) {
                    x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                                       x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
                                  (1 + 4 * a);
                }
            }
            this.set_bnd(b, x);
        }
    }

    advect(b, d, d0) {
        for (let i = 1; i <= N; i++) {
            for (let j = 1; j <= N; j++) {
                let x = i - dt_ns * N * this.u[IX(i, j)];
                let y = j - dt_ns * N * this.v[IX(i, j)];
                if (x < 0.5) x = 0.5;
                if (x > N + 0.5) x = N + 0.5;
                if (y < 0.5) y = 0.5;
                if (y > N + 0.5) y = N + 0.5;
                const i0 = Math.floor(x);
                const i1 = i0 + 1;
                const j0 = Math.floor(y);
                const j1 = j0 + 1;
                const s1 = x - i0;
                const s0 = 1 - s1;
                const t1 = y - j0;
                const t0 = 1 - t1;
                d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                              s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
            }
        }
        this.set_bnd(b, d);
    }

    project(u, v, p, div) {
        const h = 1.0 / N;
        for (let i = 1; i <= N; i++) {
            for (let j = 1; j <= N; j++) {
                div[IX(i, j)] = -0.5 * h * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + 
                                            v[IX(i, j + 1)] - v[IX(i, j - 1)]);
                p[IX(i, j)] = 0;
            }
        }
        this.set_bnd(0, div);
        this.set_bnd(0, p);
        for (let k = 0; k < 20; k++) {
            for (let i = 1; i <= N; i++) {
                for (let j = 1; j <= N; j++) {
                    p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] +
                                   p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4;
                }
            }
            this.set_bnd(0, p);
        }
        for (let i = 1; i <= N; i++) {
            for (let j = 1; j <= N; j++) {
                u[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
                v[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
            }
        }
        this.set_bnd(1, u);
        this.set_bnd(2, v);
    }

    set_bnd(b, x) {
        for (let i = 1; i <= N; i++) {
            x[IX(0, i)] = (b === 1) ? -x[IX(1, i)] : x[IX(1, i)];
            x[IX(N + 1, i)] = (b === 1) ? -x[IX(N, i)] : x[IX(N, i)];
            x[IX(i, 0)] = (b === 2) ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, N + 1)] = (b === 2) ? -x[IX(i, N)] : x[IX(i, N)];
        }
        x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
        x[IX(0, N + 1)] = 0.5 * (x[IX(1, N + 1)] + x[IX(0, N)]);
        x[IX(N + 1, 0)] = 0.5 * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
        x[IX(N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
    }

    step(uqff_g = 0.0) {
        for (let i = 1; i <= N; i++) {
            for (let j = 1; j <= N; j++) {
                this.v[IX(i, j)] += dt_ns * uqff_g;
            }
        }

        this.diffuse(1, this.u_prev, this.u, visc);
        this.diffuse(2, this.v_prev, this.v, visc);
        this.project(this.u_prev, this.v_prev, this.u, this.v);
        this.advect(1, this.u, this.u_prev);
        this.advect(2, this.v, this.v_prev);
        this.project(this.u, this.v, this.u_prev, this.v_prev);
    }

    add_jet_force(force) {
        for (let i = Math.floor(N / 4); i <= Math.floor(3 * N / 4); i++) {
            this.v[IX(i, Math.floor(N / 2))] += force;
        }
    }

    print_velocity_field() {
        console.log("Velocity field (magnitude):");
        for (let j = N; j >= 1; j--) {
            let row = "";
            for (let i = 1; i <= N; i++) {
                const mag = Math.sqrt(this.u[IX(i, j)] * this.u[IX(i, j)] + 
                                      this.v[IX(i, j)] * this.v[IX(i, j)]);
                const sym = (mag > 1.0) ? '#' : (mag > 0.5) ? '+' : (mag > 0.1) ? '.' : ' ';
                row += sym;
            }
            console.log(row);
        }
    }

    getVelocityMagnitude(i, j) {
        return Math.sqrt(this.u[IX(i, j)] * this.u[IX(i, j)] + 
                        this.v[IX(i, j)] * this.v[IX(i, j)]);
    }
}

// ============================================================================
// UQFFModule5JS - Self-Expanding Unified Field Module
// ============================================================================
class UQFFModule5JS {
    constructor() {
        this.dynamic_terms = [];
        this.dynamic_parameters = new Map();
        this.metadata = {
            version: "2.0-Enhanced",
            created: "2025-11-08",
            framework: "Self-Expanding UQFF",
            source: "source5.js converted from Source5.cpp"
        };
        this.learning_rate = 0.01;
        this.logging_enabled = false;
    }

    log(message) {
        if (this.logging_enabled) {
            console.log(`[UQFFModule5JS] ${message}`);
        }
    }

    // ========================================================================
    // Self-Expanding Framework Methods
    // ========================================================================
    
    /**
     * Install a new physics term at runtime
     * @param {PhysicsTerm} term - Physics term instance
     */
    installPhysicsTerm(term) {
        if (!(term instanceof PhysicsTerm)) {
            throw new Error("Term must extend PhysicsTerm base class");
        }
        this.log(`Installing dynamic term: ${term.description()}`);
        this.dynamic_terms.push(term);
    }

    /**
     * Set a dynamic parameter
     * @param {string} name - Parameter name
     * @param {number} value - Parameter value
     */
    setDynamicParameter(name, value) {
        this.log(`Setting parameter: ${name} = ${value}`);
        this.dynamic_parameters.set(name, value);
    }

    /**
     * Get a dynamic parameter
     * @param {string} name - Parameter name
     * @param {number} defaultVal - Default value if not found
     * @returns {number}
     */
    getDynamicParameter(name, defaultVal = 0.0) {
        return this.dynamic_parameters.has(name) ? this.dynamic_parameters.get(name) : defaultVal;
    }

    /**
     * Compute all dynamic term contributions
     * @param {Object} params - Parameter map
     * @returns {number}
     */
    computeDynamicContributions(params) {
        let sum = 0.0;
        for (const term of this.dynamic_terms) {
            sum += term.compute(params);
        }
        return sum;
    }

    /**
     * Export state for cross-module communication
     * @returns {Object}
     */
    exportState() {
        const state = {
            version: this.metadata.version,
            created: this.metadata.created,
            framework: this.metadata.framework,
            parameters: Object.fromEntries(this.dynamic_parameters),
            terms: this.dynamic_terms.map(t => t.description()),
            learning_rate: this.learning_rate
        };
        this.log("State exported");
        return state;
    }

    /**
     * Import state from another module
     * @param {Object} state - State object
     */
    importState(state) {
        if (state.parameters) {
            for (const [key, value] of Object.entries(state.parameters)) {
                this.dynamic_parameters.set(key, value);
            }
        }
        if (state.learning_rate) {
            this.learning_rate = state.learning_rate;
        }
        this.log("State imported");
    }

    /**
     * Set learning rate for optimization
     * @param {number} rate - Learning rate
     */
    setLearningRate(rate) {
        this.learning_rate = rate;
        this.log(`Learning rate set to: ${rate}`);
    }

    /**
     * Enable/disable logging
     * @param {boolean} enable - Enable flag
     */
    setEnableLogging(enable) {
        this.logging_enabled = enable;
    }

    /**
     * Get module information
     * @returns {Object}
     */
    getInfo() {
        return {
            version: this.metadata.version,
            dynamicTerms: this.dynamic_terms.length,
            dynamicParameters: this.dynamic_parameters.size,
            learningRate: this.learning_rate,
            logging: this.logging_enabled
        };
    }

    /**
     * Print module information
     */
    printInfo() {
        console.log("=== UQFFModule5JS Info ===");
        console.log(`Version: ${this.metadata.version}`);
        console.log(`Dynamic Terms: ${this.dynamic_terms.length}`);
        console.log(`Dynamic Parameters: ${this.dynamic_parameters.size}`);
        console.log(`Learning Rate: ${this.learning_rate}`);
        console.log(`Logging: ${this.logging_enabled ? 'Enabled' : 'Disabled'}`);
    }

    /**
     * Tune a parameter using gradient-based optimization
     * @param {string} paramName - Parameter to tune
     * @param {Function} lossFunction - Loss function
     * @param {number} iterations - Number of iterations
     */
    tuneParameter(paramName, lossFunction, iterations = 100) {
        let value = this.getDynamicParameter(paramName, 1.0);
        const epsilon = 1e-8;

        for (let i = 0; i < iterations; i++) {
            const loss = lossFunction(value);
            const lossPlus = lossFunction(value + epsilon);
            const gradient = (lossPlus - loss) / epsilon;
            value -= this.learning_rate * gradient;
            this.setDynamicParameter(paramName, value);
        }

        this.log(`Parameter ${paramName} tuned to ${value}`);
        return value;
    }

    /**
     * Optimize all dynamic parameters
     * @param {Function} lossFunction - Loss function
     * @param {number} iterations - Number of iterations
     */
    optimize(lossFunction, iterations = 100) {
        this.log("Starting optimization...");
        for (const [paramName] of this.dynamic_parameters) {
            this.tuneParameter(paramName, lossFunction, iterations);
        }
        this.log("Optimization complete");
    }

    // ========================================================================
    // Enhanced Compute Functions (with dynamic contributions)
    // ========================================================================
    
    compute_Ug1_enhanced(body, r, t, tn, alpha_param, delta_def_param, k1_param) {
        const original = compute_Ug1(body, r, t, tn, alpha_param, delta_def_param, k1_param);
        const params = {r, t, tn};
        const dynamic = this.computeDynamicContributions(params);
        return original + dynamic;
    }

    compute_Ug2_enhanced(body, r, t, tn, k2_param, QA_param, delta_sw_param, v_sw_param, HSCm_param, rho_A_param, kappa_param) {
        const original = compute_Ug2(body, r, t, tn, k2_param, QA_param, delta_sw_param, v_sw_param, HSCm_param, rho_A_param, kappa_param);
        const params = {r, t, tn};
        const dynamic = this.computeDynamicContributions(params);
        return original + dynamic;
    }

    compute_Ug3_enhanced(body, r, t, tn, theta, rho_A_param, kappa_param, k3_param) {
        const original = compute_Ug3(body, r, t, tn, theta, rho_A_param, kappa_param, k3_param);
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

    compute_FU_enhanced(body, r, t, tn, theta) {
        const original = compute_FU(body, r, t, tn, theta);
        const params = {r, t, tn, theta};
        const dynamic = this.computeDynamicContributions(params);
        return original + dynamic;
    }

    // ========================================================================
    // Simulation Methods
    // ========================================================================
    
    /**
     * Simulate quasar jet using Navier-Stokes fluid dynamics
     * @param {number} initial_velocity - Initial jet velocity
     * @param {number} steps - Number of simulation steps
     * @returns {FluidSolver}
     */
    simulate_quasar_jet(initial_velocity, steps = 10) {
        try {
            const solver = new FluidSolver();
            solver.add_jet_force(initial_velocity / 10.0);

            const res = new ResonanceParams();
            const sagA = new MUGESystem({
                name: "Sagittarius A*",
                I: 1e23,
                A: 2.813e30,
                omega1: 1e-5,
                omega2: -1e-5,
                Vsys: 3.552e45,
                vexp: 5e6,
                t: 3.786e14,
                ffluid: 3.465e-8,
                r: 1e12,
                M: 8.155e36,
                B: 1e-5,
                Bcrit: 1e-4,
                rho_fluid: 1e-20,
                g_local: 1e-5,
                M_DM: 1e37,
                delta_rho_rho: 1e-3
            });

            const uqff_g = compute_resonance_MUGE(sagA, res);
            this.log(`Simulating quasar jet (${steps} steps) with UQFF g=${uqff_g}...`);

            for (let step = 0; step < steps; step++) {
                solver.step(uqff_g / 1e30);
            }

            return solver;
        } catch (e) {
            console.error("Error in simulate_quasar_jet:", e.message);
            return null;
        }
    }

    /**
     * Simulate celestial body evolution over time
     * @param {CelestialBody} body - Celestial body
     * @param {number} t_start - Start time
     * @param {number} t_end - End time
     * @param {number} dt - Time step
     * @returns {Array} Time series data
     */
    simulate_celestial_evolution(body, t_start, t_end, dt) {
        const results = [];
        const r = body.Rb;
        const theta = 0.0;

        for (let t = t_start; t < t_end; t += dt) {
            const tn = t;
            const FU = this.compute_FU_enhanced(body, r, t, tn, theta);
            const Ug1 = this.compute_Ug1_enhanced(body, r, t, tn, alpha, delta_def, k1);
            const Ug2 = this.compute_Ug2_enhanced(body, r, t, tn, k2, QA, delta_sw, v_sw, HSCm, rho_A, kappa);
            const Ug3 = this.compute_Ug3_enhanced(body, r, t, tn, theta, rho_A, kappa, k3);

            results.push({
                t,
                FU,
                Ug1,
                Ug2,
                Ug3,
                body_name: body.name
            });
        }

        this.log(`Simulated ${results.length} time steps for ${body.name}`);
        return results;
    }

    /**
     * Simulate MUGE system evolution
     * @param {MUGESystem} sys - MUGE system
     * @param {ResonanceParams} res - Resonance parameters
     * @param {number} t_start - Start time
     * @param {number} t_end - End time
     * @param {number} dt - Time step
     * @returns {Array} Time series data
     */
    simulate_MUGE_evolution(sys, res, t_start, t_end, dt) {
        const results = [];
        const original_t = sys.t;

        for (let t = t_start; t < t_end; t += dt) {
            sys.t = t;
            const compressed_g = compute_compressed_MUGE(sys);
            const resonance_g = compute_resonance_MUGE(sys, res);
            const enhanced_g = this.compute_MUGE_enhanced(sys, res, true);

            results.push({
                t,
                compressed_g,
                resonance_g,
                enhanced_g,
                system_name: sys.name
            });
        }

        sys.t = original_t; // Restore original time
        this.log(`Simulated ${results.length} time steps for ${sys.name}`);
        return results;
    }

    /**
     * Run comprehensive field analysis for a celestial body
     * @param {CelestialBody} body - Celestial body
     * @param {number} r - Radial distance
     * @param {number} t - Time
     * @returns {Object} Analysis results
     */
    analyze_celestial_body(body, r, t) {
        const tn = t;
        const theta = 0.0;

        const FU = this.compute_FU_enhanced(body, r, t, tn, theta);
        const Ug1 = this.compute_Ug1_enhanced(body, r, t, tn, alpha, delta_def, k1);
        const Ug2 = this.compute_Ug2_enhanced(body, r, t, tn, k2, QA, delta_sw, v_sw, HSCm, rho_A, kappa);
        const Ug3 = this.compute_Ug3_enhanced(body, r, t, tn, theta, rho_A, kappa, k3);
        const Ug4_val = compute_Ug4(t, tn, rho_v, C_concentration, Mbh, dg, alpha, f_feedback, k4);
        const Um = compute_Um(body, t, tn, body.Rb, gamma, rho_A, kappa, num_strings);

        return {
            body_name: body.name,
            t,
            r,
            FU,
            components: {
                Ug1,
                Ug2,
                Ug3,
                Ug4: Ug4_val,
                Um
            },
            total_field_strength: FU,
            dynamic_contributions: this.computeDynamicContributions({r, t, tn, theta})
        };
    }

    /**
     * Run comprehensive MUGE system analysis
     * @param {MUGESystem} sys - MUGE system
     * @param {ResonanceParams} res - Resonance parameters
     * @returns {Object} Analysis results
     */
    analyze_MUGE_system(sys, res) {
        const compressed_g = compute_compressed_MUGE(sys);
        const resonance_g = compute_resonance_MUGE(sys, res);
        const enhanced_g = this.compute_MUGE_enhanced(sys, res, true);

        const aDPM = compute_aDPM(sys, res);
        const components = {
            aDPM,
            aTHz: compute_aTHz(aDPM, sys, res),
            avac_diff: compute_avac_diff(aDPM, sys, res),
            asuper_freq: compute_asuper_freq(aDPM, res),
            aaether_res: compute_aaether_res(aDPM, res),
            Ug4i: compute_Ug4i(aDPM, sys, res),
            aquantum_freq: compute_aquantum_freq(aDPM, res),
            aAether_freq: compute_aAether_freq(aDPM, res),
            afluid_freq: compute_afluid_freq(sys, res),
            Osc_term: compute_Osc_term(),
            aexp_freq: compute_aexp_freq(aDPM, sys, res),
            fTRZ: compute_fTRZ(res),
            a_wormhole: compute_a_wormhole(sys.r)
        };

        return {
            system_name: sys.name,
            compressed_g,
            resonance_g,
            enhanced_g,
            components,
            dynamic_contributions: this.computeDynamicContributions({t: sys.t, r: sys.r, M: sys.M})
        };
    }
}

// ============================================================================
// Default Celestial Bodies
// ============================================================================
function createDefaultBodies() {
    return [
        new CelestialBody({
            name: "Sun",
            Ms: 1.989e30,
            Rs: 6.96e8,
            Rb: 1.496e13,
            Ts_surface: 5778.0,
            omega_s: 2.5e-6,
            Bs_avg: 1e-4,
            SCm_density: 1e15,
            QUA: 1e-11,
            Pcore: 1.0,
            PSCm: 1.0,
            omega_c: 2 * PI / (11.0 * 365.25 * 24 * 3600)
        }),
        new CelestialBody({
            name: "Earth",
            Ms: 5.972e24,
            Rs: 6.371e6,
            Rb: 1e7,
            Ts_surface: 288.0,
            omega_s: 7.292e-5,
            Bs_avg: 3e-5,
            SCm_density: 1e12,
            QUA: 1e-12,
            Pcore: 1e-3,
            PSCm: 1e-3,
            omega_c: 2 * PI / (1.0 * 365.25 * 24 * 3600)
        }),
        new CelestialBody({
            name: "Jupiter",
            Ms: 1.898e27,
            Rs: 6.9911e7,
            Rb: 1e8,
            Ts_surface: 165.0,
            omega_s: 1.76e-4,
            Bs_avg: 4e-4,
            SCm_density: 1e13,
            QUA: 1e-11,
            Pcore: 1e-3,
            PSCm: 1e-3,
            omega_c: 2 * PI / (11.86 * 365.25 * 24 * 3600)
        }),
        new CelestialBody({
            name: "Neptune",
            Ms: 1.024e26,
            Rs: 2.4622e7,
            Rb: 5e7,
            Ts_surface: 72.0,
            omega_s: 1.08e-4,
            Bs_avg: 1e-4,
            SCm_density: 1e11,
            QUA: 1e-13,
            Pcore: 1e-3,
            PSCm: 1e-3,
            omega_c: 2 * PI / (164.8 * 365.25 * 24 * 3600)
        })
    ];
}

// ============================================================================
// Default MUGE Systems
// ============================================================================
function createDefaultMUGESystems() {
    return [
        new MUGESystem({
            name: "Magnetar SGR 1745-2900",
            I: 1e21,
            A: 3.142e8,
            omega1: 1e-3,
            omega2: -1e-3,
            Vsys: 4.189e12,
            vexp: 1e3,
            t: 3.799e10,
            z: 0.0009,
            ffluid: 1.269e-14,
            M: 2.984e30,
            r: 1e4,
            B: 1e10,
            Bcrit: 1e11,
            rho_fluid: 1e-15,
            g_local: 10.0,
            M_DM: 0.0,
            delta_rho_rho: 1e-5
        }),
        new MUGESystem({
            name: "Sagittarius A*",
            I: 1e23,
            A: 2.813e30,
            omega1: 1e-5,
            omega2: -1e-5,
            Vsys: 3.552e45,
            vexp: 5e6,
            t: 3.786e14,
            z: 0.0009,
            ffluid: 3.465e-8,
            M: 8.155e36,
            r: 1e12,
            B: 1e-5,
            Bcrit: 1e-4,
            rho_fluid: 1e-20,
            g_local: 1e-5,
            M_DM: 1e37,
            delta_rho_rho: 1e-3
        }),
        new MUGESystem({
            name: "Tapestry of Blazing Starbirth",
            I: 1e22,
            A: 1e35,
            omega1: 1e-4,
            omega2: -1e-4,
            Vsys: 1e53,
            vexp: 1e4,
            t: 3.156e13,
            z: 0.0,
            ffluid: 1e-12,
            M: 1.989e35,
            r: 3.086e17,
            B: 1e-4,
            Bcrit: 1e-3,
            rho_fluid: 1e-21,
            g_local: 1e-8,
            M_DM: 1e35,
            delta_rho_rho: 1e-4
        })
    ];
}

// ============================================================================
// Exports
// ============================================================================
module.exports = {
    // Classes
    PhysicsTerm,
    DarkMatterHaloTerm,
    VacuumEnergyTerm,
    CelestialBody,
    ResonanceParams,
    MUGESystem,
    FluidSolver,
    UQFFModule5JS,

    // Computation Functions
    compute_Ug1,
    compute_Ug2,
    compute_Ug3,
    compute_Ug4,
    compute_Um,
    compute_Ubi,
    compute_FU,
    compute_compressed_MUGE,
    compute_resonance_MUGE,
    compute_aDPM,
    compute_aTHz,
    compute_avac_diff,
    compute_asuper_freq,
    compute_aaether_res,
    compute_Ug4i,
    compute_aquantum_freq,
    compute_aAether_freq,
    compute_afluid_freq,
    compute_a_wormhole,

    // Factory Functions
    createDefaultBodies,
    createDefaultMUGESystems,

    // Constants
    PI,
    c,
    G,
    N
};
