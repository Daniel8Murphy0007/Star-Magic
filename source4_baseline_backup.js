// source4.js - Unified Field Theory Implementation (Star Magic - The Quest for Unity)
// JavaScript conversion of source4.cpp with complete physics fidelity
// Implements: Unified Field Equation (FU), MUGE (compressed & resonance), Navier-Stokes fluid solver
// Copyright ©2025 Daniel T. Murphy, daniel.murphy00@gmail.com – All Rights Reserved

// ========== PHYSICAL CONSTANTS ==========
const PI = Math.PI;
const c = 3.0e8;              // Speed of light (m/s)
const G = 6.67430e-11;        // Gravitational constant (m³/kg·s²)

// Galactic parameters
const Omega_g = 7.3e-16;      // Galactic spin rate (rad/s)
const M_bh = 8.15e36;         // Black hole mass (kg)
const d_g = 8.178e3;          // Galactic center distance (pc)

// Speculative shared parameters (from compressed UQFF, resonance equations)
const v_SCm = 0.99 * c;       // [SCm] relativistic velocity (m/s)
const rho_A = 1e-23;          // Aether density (kg/m³)
const rho_sw = 1e-21;         // Solar wind density (kg/m³)
const v_sw = 4e5;             // Solar wind velocity (m/s)
const Q_A = 1.602e-19;        // Aether charge unit (C)
const Q_s = 1.602e-19;        // Star charge unit (C)
const alpha_t = 1e-9;         // Exponential decay rate (s⁻¹)
const beta_i = 0.1;           // Buoyancy coupling
const kappa_t = 5e-4;         // Reactor decay rate (day⁻¹)
const n_t = 1.0;              // Temporal cycle count
const eta_A = 0.01;           // Aether modulation factor
const C_vac = 1.0;            // Vacuum coupling constant
const f_feedback = 1.0;       // Feedback factor
const H_SCM = 1.0;            // [SCm] Heaviside step result
const U_UA = 1.0;             // Universal Aether contribution
const S_wind = 1.0;           // Solar wind modulation
const num_magnetic_strings = 10; // String count

// Vacuum energy density
const rho_v = 6e-27;          // kg/m³

// Stress-energy tensor component
const Ts00 = 1.27e3 + 1.11e7; // Example value (Pa)

// Coupling constants
const k1 = 1.0, k2 = 1.0, k3 = 1.0, k4 = 1.0;

// Background Aether metric (4x4 identity-like tensor)
const g_mu_nu = [
  [1, 0, 0, 0],
  [0, -1, 0, 0],
  [0, 0, -1, 0],
  [0, 0, 0, -1]
];

// ========== CELESTIAL BODY STRUCTURE ==========
class CelestialBody {
  constructor(name, Ms, Rs, Rb, Ts_surface, omega_s, Bs_avg, SCm_density, QUA, Pcore, PSCm, omega_c) {
    this.name = name;
    this.Ms = Ms;                 // Mass (kg)
    this.Rs = Rs;                 // Radius (m)
    this.Rb = Rb;                 // Bubble radius (m)
    this.Ts_surface = Ts_surface; // Surface temperature (K)
    this.omega_s = omega_s;       // Rotation rate (rad/s)
    this.Bs_avg = Bs_avg;         // Average magnetic field (T)
    this.SCm_density = SCm_density; // [SCm] density (kg/m³)
    this.QUA = QUA;               // Unified aether charge (C)
    this.Pcore = Pcore;           // Core penetration factor
    this.PSCm = PSCm;             // [SCm] penetration factor
    this.omega_c = omega_c;       // Cycle frequency (rad/s)
  }
}

// ========== UTILITY FUNCTIONS ==========
function step_function(r, Rb) {
  return (r > Rb) ? 1.0 : 0.0;
}

function compute_Ereact(t, rho_SCm, v_scm_local, rho_A_local, kappa) {
  return (rho_SCm * v_scm_local * v_scm_local / rho_A_local) * Math.exp(-kappa * t);
}

function compute_mu_s(t, Bs_t, Rs) {
  return Bs_t * Math.pow(Rs, 3);
}

function compute_grad_Ms_r(Ms, Rs) {
  return G * Ms / (Rs * Rs);
}

function compute_Bj(t, omega_c, SCm_density) {
  return 1e-3 + 0.4 * Math.sin(omega_c * t) + SCm_density;
}

function compute_omega_s_t(t, omega_s, decay_factor = 1e-9) {
  return omega_s * Math.exp(-decay_factor * t);
}

function compute_mu_j(t, Bj, omega_c) {
  return Bj * Math.sin(omega_c * t);
}

// ========== MAIN PHYSICS FUNCTIONS ==========

// Universal Gravity Component 1: Dipole + defect modulation
function compute_Ug1(t, body, r, defect_factor = 1.0) {
  const Bs_t = body.Bs_avg;
  const mu_s = compute_mu_s(t, Bs_t, body.Rs);
  const grad_Ms_r = compute_grad_Ms_r(body.Ms, body.Rs);
  const decay = Math.exp(-alpha_t * t);
  const temporal_mod = Math.cos(PI * t * n_t);
  return k1 * mu_s * grad_Ms_r * decay * temporal_mod * defect_factor;
}

// Universal Gravity Component 2: Charge + reactor + step function
function compute_Ug2(t, body, r) {
  const step = step_function(r, body.Rb);
  const Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A, kappa_t);
  return k2 * body.QUA * body.Ms / (r * r) * S_wind * H_SCM * Ereact * step;
}

// Universal Gravity Component 3: Magnetic strings + core penetration
function compute_Ug3(t, body, r) {
  const Bj = compute_Bj(t, body.omega_c, body.SCm_density);
  const omega_s_t = compute_omega_s_t(t, body.omega_s);
  const Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A, kappa_t);
  return k3 * Bj * Math.cos(omega_s_t * PI) * body.Pcore * Ereact;
}

// Universal Gravity Component 4: Vacuum energy + black hole + feedback
function compute_Ug4(t, body) {
  const decay = Math.exp(-alpha_t * t);
  const temporal_mod = Math.cos(PI * t * n_t);
  return k4 * rho_v * C_vac * M_bh / d_g * decay * temporal_mod * f_feedback;
}

// Universal Buoyancy: Negative coupling to each Ug component
function compute_Ubi(Ugi, t) {
  const temporal_mod = Math.cos(PI * t * n_t);
  return -beta_i * Ugi * Omega_g * M_bh / d_g * S_wind * U_UA * temporal_mod;
}

// Universal Magnetism: String-based, optimized
function compute_Um(t, body, r_j) {
  const mu_j = compute_mu_j(t, body.Bs_avg, body.omega_c);
  const decay = Math.exp(-alpha_t * t);
  const Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A, kappa_t);
  return (mu_j / r_j) * decay * num_magnetic_strings * body.PSCm * Ereact;
}

// Aether Tensor: Modulated by stress-energy
function compute_A_mu_nu(t) {
  const temporal_mod = Math.cos(PI * t * n_t);
  const result = [];
  for (let i = 0; i < 4; i++) {
    result[i] = [];
    for (let j = 0; j < 4; j++) {
      result[i][j] = g_mu_nu[i][j] + eta_A * Ts00 * temporal_mod;
    }
  }
  return result;
}

// Master Unified Field Equation
function compute_FU(t, body, r, r_j) {
  const Ug1 = compute_Ug1(t, body, r);
  const Ug2 = compute_Ug2(t, body, r);
  const Ug3 = compute_Ug3(t, body, r);
  const Ug4 = compute_Ug4(t, body);
  
  const Ubi1 = compute_Ubi(Ug1, t);
  const Ubi2 = compute_Ubi(Ug2, t);
  const Ubi3 = compute_Ubi(Ug3, t);
  const Ubi4 = compute_Ubi(Ug4, t);
  
  const Um = compute_Um(t, body, r_j);
  
  const A_tensor = compute_A_mu_nu(t);
  const trace_A = A_tensor[0][0] + A_tensor[1][1] + A_tensor[2][2] + A_tensor[3][3];
  
  return Ug1 + Ug2 + Ug3 + Ug4 + Ubi1 + Ubi2 + Ubi3 + Ubi4 + Um + trace_A;
}

// ========== NAVIER-STOKES FLUID SOLVER ==========
class FluidSolver {
  constructor(N = 32, dt = 0.1, visc = 0.0001, force_jet = 10.0) {
    this.N = N;
    this.dt = dt;
    this.visc = visc;
    this.force_jet = force_jet;
    
    const size = (N + 2) * (N + 2);
    this.u = new Array(size).fill(0);      // Velocity x
    this.v = new Array(size).fill(0);      // Velocity y
    this.dens = new Array(size).fill(0);   // Density
    this.u_prev = new Array(size).fill(0);
    this.v_prev = new Array(size).fill(0);
    this.dens_prev = new Array(size).fill(0);
  }
  
  IX(i, j) {
    return i + (this.N + 2) * j;
  }
  
  add_source(x, s) {
    for (let i = 0; i < x.length; i++) {
      x[i] += this.dt * s[i];
    }
  }
  
  set_bnd(b, x) {
    const N = this.N;
    for (let i = 1; i <= N; i++) {
      x[this.IX(0, i)] = (b === 1) ? -x[this.IX(1, i)] : x[this.IX(1, i)];
      x[this.IX(N + 1, i)] = (b === 1) ? -x[this.IX(N, i)] : x[this.IX(N, i)];
      x[this.IX(i, 0)] = (b === 2) ? -x[this.IX(i, 1)] : x[this.IX(i, 1)];
      x[this.IX(i, N + 1)] = (b === 2) ? -x[this.IX(i, N)] : x[this.IX(i, N)];
    }
    x[this.IX(0, 0)] = 0.5 * (x[this.IX(1, 0)] + x[this.IX(0, 1)]);
    x[this.IX(0, N + 1)] = 0.5 * (x[this.IX(1, N + 1)] + x[this.IX(0, N)]);
    x[this.IX(N + 1, 0)] = 0.5 * (x[this.IX(N, 0)] + x[this.IX(N + 1, 1)]);
    x[this.IX(N + 1, N + 1)] = 0.5 * (x[this.IX(N, N + 1)] + x[this.IX(N + 1, N)]);
  }
  
  diffuse(b, x, x0, diff) {
    const a = this.dt * diff * this.N * this.N;
    for (let k = 0; k < 20; k++) {
      for (let i = 1; i <= this.N; i++) {
        for (let j = 1; j <= this.N; j++) {
          x[this.IX(i, j)] = (x0[this.IX(i, j)] + a * (
            x[this.IX(i - 1, j)] + x[this.IX(i + 1, j)] +
            x[this.IX(i, j - 1)] + x[this.IX(i, j + 1)]
          )) / (1 + 4 * a);
        }
      }
      this.set_bnd(b, x);
    }
  }
  
  advect(b, d, d0) {
    const N = this.N;
    for (let i = 1; i <= N; i++) {
      for (let j = 1; j <= N; j++) {
        let x = i - this.dt * N * this.u[this.IX(i, j)];
        let y = j - this.dt * N * this.v[this.IX(i, j)];
        
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
        
        d[this.IX(i, j)] = s0 * (t0 * d0[this.IX(i0, j0)] + t1 * d0[this.IX(i0, j1)]) +
                           s1 * (t0 * d0[this.IX(i1, j0)] + t1 * d0[this.IX(i1, j1)]);
      }
    }
    this.set_bnd(b, d);
  }
  
  project(u, v, p, div) {
    const N = this.N;
    const h = 1.0 / N;
    
    for (let i = 1; i <= N; i++) {
      for (let j = 1; j <= N; j++) {
        div[this.IX(i, j)] = -0.5 * h * (
          u[this.IX(i + 1, j)] - u[this.IX(i - 1, j)] +
          v[this.IX(i, j + 1)] - v[this.IX(i, j - 1)]
        );
        p[this.IX(i, j)] = 0;
      }
    }
    
    this.set_bnd(0, div);
    this.set_bnd(0, p);
    
    for (let k = 0; k < 20; k++) {
      for (let i = 1; i <= N; i++) {
        for (let j = 1; j <= N; j++) {
          p[this.IX(i, j)] = (div[this.IX(i, j)] +
            p[this.IX(i - 1, j)] + p[this.IX(i + 1, j)] +
            p[this.IX(i, j - 1)] + p[this.IX(i, j + 1)]) / 4;
        }
      }
      this.set_bnd(0, p);
    }
    
    for (let i = 1; i <= N; i++) {
      for (let j = 1; j <= N; j++) {
        u[this.IX(i, j)] -= 0.5 * (p[this.IX(i + 1, j)] - p[this.IX(i - 1, j)]) / h;
        v[this.IX(i, j)] -= 0.5 * (p[this.IX(i, j + 1)] - p[this.IX(i, j - 1)]) / h;
      }
    }
    
    this.set_bnd(1, u);
    this.set_bnd(2, v);
  }
  
  step(uqff_g = 0.0) {
    // Add UQFF gravity-like force as body force in v (vertical direction)
    for (let i = 1; i <= this.N; i++) {
      for (let j = 1; j <= this.N; j++) {
        this.v[this.IX(i, j)] += this.dt * uqff_g;
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
    // Add force in center as a jet (simulating [SCm] expulsion)
    for (let i = Math.floor(this.N / 4); i <= Math.floor(3 * this.N / 4); i++) {
      this.v[this.IX(i, Math.floor(this.N / 2))] += force;
    }
  }
  
  print_velocity_field() {
    console.log("Velocity field (magnitude):");
    for (let j = this.N; j >= 1; j--) {
      let row = "";
      for (let i = 1; i <= this.N; i++) {
        const mag = Math.sqrt(
          this.u[this.IX(i, j)] * this.u[this.IX(i, j)] +
          this.v[this.IX(i, j)] * this.v[this.IX(i, j)]
        );
        const sym = (mag > 1.0) ? '#' : (mag > 0.5) ? '+' : (mag > 0.1) ? '.' : ' ';
        row += sym;
      }
      console.log(row);
    }
  }
}

// ========== RESONANCE PARAMETERS ==========
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

// ========== MUGE SYSTEM ==========
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

// ========== COMPRESSED MUGE FUNCTIONS ==========
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

// ========== RESONANCE MUGE FUNCTIONS ==========
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
  
  return aDPM + aTHz + avac_diff + asuper_freq + aaether_res + Ug4i + 
         aquantum_freq + aAether_freq + afluid_freq + Osc_term + aexp_freq + fTRZ + a_worm;
}

// ========== SIMULATION FUNCTIONS ==========
function simulate_quasar_jet(initial_velocity) {
  const solver = new FluidSolver();
  solver.add_jet_force(initial_velocity / 10.0);
  
  // Integrate UQFF: Use Sgr A* as example
  const res = new ResonanceParams();
  const sagA = new MUGESystem(
    "Sagittarius A*",
    1e23, 2.813e30, 1e-5, -1e-5, 3.552e45, 5e6, 3.786e14, 0.0009, 3.465e-8,
    8.155e36, 1e12, 1e-5, 1e-4, 1e-20, 1e-5, 1e37, 1e-3
  );
  const uqff_g = compute_resonance_MUGE(sagA, res);
  
  console.log(`Simulating quasar jet with Navier-Stokes (10 steps) using UQFF g=${uqff_g}...`);
  for (let step = 0; step < 10; step++) {
    solver.step(uqff_g / 1e30); // Scale g to avoid numerical blowup
  }
  solver.print_velocity_field();
}

// ========== MAIN EXECUTION ==========
function main() {
  console.log("========================================");
  console.log("Star Magic - Unified Field Theory (JS)");
  console.log("Copyright ©2025 Daniel T. Murphy");
  console.log("========================================\n");
  
  // Define celestial body (Sun)
  const sun = new CelestialBody(
    "Sun",
    1.989e30,                    // Ms
    6.96e8,                      // Rs
    1.496e13,                    // Rb
    5778.0,                      // Ts_surface
    2.5e-6,                      // omega_s
    1e-4,                        // Bs_avg
    1e15,                        // SCm_density
    1e-11,                       // QUA
    1.0,                         // Pcore
    1.0,                         // PSCm
    2 * PI / (11.0 * 365.25 * 24 * 3600) // omega_c (11-year cycle)
  );
  
  // Test unified field computation
  const t = 0;
  const r = 1.496e11;  // 1 AU
  const r_j = 1e10;
  
  const FU = compute_FU(t, sun, r, r_j);
  console.log(`Unified Field (FU) at t=${t}s, r=${r}m: ${FU.toExponential(3)}\n`);
  
  // Simulate quasar jet using Navier-Stokes
  simulate_quasar_jet(v_SCm);
  
  // Define MUGE systems
  const res_params = new ResonanceParams();
  
  const sgr1745 = new MUGESystem(
    "Magnetar SGR 1745-2900",
    1e21, 3.142e8, 1e-3, -1e-3, 4.189e12, 1e3, 3.799e10, 0.0009, 1.269e-14,
    2.984e30, 1e4, 1e10, 1e11, 1e-15, 10.0, 0.0, 1e-5
  );
  
  const sagA = new MUGESystem(
    "Sagittarius A*",
    1e23, 2.813e30, 1e-5, -1e-5, 3.552e45, 5e6, 3.786e14, 0.0009, 3.465e-8,
    8.155e36, 1e12, 1e-5, 1e-4, 1e-20, 1e-5, 1e37, 1e-3
  );
  
  const tapestry = new MUGESystem(
    "Tapestry of Blazing Starbirth",
    1e22, 1e35, 1e-4, -1e-4, 1e53, 1e4, 3.156e13, 0.0, 1e-12,
    1.989e35, 3.086e17, 1e-4, 1e-3, 1e-21, 1e-8, 1e35, 1e-4
  );
  
  const westerlund = new MUGESystem(
    "Westerlund 2",
    1e22, 1e35, 1e-4, -1e-4, 1e53, 1e4, 3.156e13, 0.0, 1e-12,
    1.989e35, 3.086e17, 1e-4, 1e-3, 1e-21, 1e-8, 1e35, 1e-4
  );
  
  const pillars = new MUGESystem(
    "Pillars of Creation",
    1e21, 2.813e32, 1e-3, -1e-3, 3.552e48, 2e3, 3.156e13, 0.0, 8.457e-14,
    1.989e32, 9.46e15, 1e-4, 1e-3, 1e-21, 1e-8, 0.0, 1e-5
  );
  
  const rings = new MUGESystem(
    "Rings of Relativity",
    1e22, 1e35, 1e-4, -1e-4, 1e54, 1e5, 3.156e14, 0.01, 1e-9,
    1.989e36, 3.086e17, 1e-5, 1e-4, 1e-20, 1e-5, 1e36, 1e-3
  );
  
  const student_guide = new MUGESystem(
    "Student's Guide to the Universe",
    1e24, 1e52, 1e-6, -1e-6, 1e80, 3e8, 4.35e17, 0.0, 1e-18,
    1e53, 1e26, 1e-10, 1e-9, 1e-30, 1e-10, 1e53, 1e-6
  );
  
  const muge_systems = [sgr1745, sagA, tapestry, westerlund, pillars, rings, student_guide];
  
  console.log("\n========================================");
  console.log("MUGE Calculations (Compressed & Resonance)");
  console.log("========================================\n");
  
  for (const sys of muge_systems) {
    const compressed_g = compute_compressed_MUGE(sys);
    const resonance_g = compute_resonance_MUGE(sys, res_params);
    console.log(`${sys.name}:`);
    console.log(`  Compressed MUGE g: ${compressed_g.toExponential(3)} m/s²`);
    console.log(`  Resonance MUGE g:  ${resonance_g.toExponential(3)} m/s²\n`);
  }
  
  console.log("========================================");
  console.log("All computations completed successfully!");
  console.log("========================================");
}

// ========== EXPORTS ==========
module.exports = {
  // Classes
  CelestialBody,
  FluidSolver,
  ResonanceParams,
  MUGESystem,
  
  // Unified Field Functions
  compute_Ug1,
  compute_Ug2,
  compute_Ug3,
  compute_Ug4,
  compute_Ubi,
  compute_Um,
  compute_A_mu_nu,
  compute_FU,
  
  // MUGE Functions
  compute_compressed_MUGE,
  compute_resonance_MUGE,
  
  // Simulation
  simulate_quasar_jet,
  
  // Main
  main
};

// Run main if executed directly
if (require.main === module) {
  main();
}

// Watermark: ©2025 Daniel T. Murphy, daniel.murphy00@gmail.com – All Rights Reserved
