// Source161.js - UQFFBuoyancyModule161
// Multi-system: J1610+1811, PLCK G287.0+32.9, PSZ2 G181.06+48.47, ASKAP J1832-0911, Sonification Collection
// SIMULATION-READY: Quasar/cluster systems with parallel computing optimization
// Copyright - Daniel T. Murphy, Nov 2025

const { addEnhancedDynamics } = require('./enhanced_dynamics.js');

function complexAdd(a, b) { return {re: a.re + b.re, im: a.im + b.im}; }

class UQFFBuoyancyModule161 {
  constructor(systemName = "J1610") {
    this.systemConfigs = {
      "J1610": { M: 7.96e40, r: 9.26e20, L_X: 1e42, B0: 1e-9, rho_gas: 1e-21, T_val: 1e7, omega0: 1e-14, x2: -1.35e172 },
      "PLCK_G287": { M: 1.59e45, r: 9.26e24, L_X: 1e45, B0: 1e-9, rho_gas: 1e-23, T_val: 1e8, omega0: 1e-15, x2: -2.27e172 },
      "PSZ2_G181": { M: 1.39e45, r: 7.72e24, L_X: 1e44, B0: 1e-9, rho_gas: 1e-23, T_val: 8e7, omega0: 1e-15, x2: -2.27e172 },
      "ASKAP_J1832": { M: 3.978e40, r: 4.63e20, L_X: 1e41, B0: 1e-9, rho_gas: 1e-21, T_val: 1e7, omega0: 1e-14, x2: -1.35e172 },
      "Sonification": { M: 1.989e41, r: 1.54e21, L_X: 1e40, B0: 1e-9, rho_gas: 1e-21, T_val: 1e7, omega0: 1e-14, x2: -1.35e172 }
    };
    
    this.currentSystem = systemName;
    this.DPM_momentum = 0.93; this.k_LENR = 1e-10; this.E_cm = 2.18e-6;
    this.G = 6.67430e-11; this.c = 299792458.0; this.k_B = 1.380649e-23; this.m_n = 1.674927498e-27;
    
    this.setSystem(systemName);
    this.dynamicTerms = []; this.dynamicParameters = new Map();
    this.metadata = new Map([["enhanced", true], ["version", "2.0.0"], ["simulation_ready", true]]);
    this.enableLogging = false;
  }
  
  setSystem(systemName) {
    const config = this.systemConfigs[systemName];
    if (!config) throw new Error(`Unknown system: ${systemName}`);
    this.currentSystem = systemName;
    this.variables = new Map(Object.entries(config).map(([k,v]) => [k, {re:v, im:0}]));
    this.variables.set("DPM_momentum", {re: this.DPM_momentum, im: 0});
    this.variables.set("k_LENR", {re: this.k_LENR, im: 0});
    this.variables.set("t", {re: 0, im: 0});
    this.variables.set("V", {re: 1e6, im: 0});
    this.variables.set("k_DE", {re: 1e-14, im: 0});
  }
  
  computeF(t) {
    const M = this.variables.get("M").re, r = this.variables.get("r").re;
    const V = this.variables.get("V").re, L_X = this.variables.get("L_X").re;
    const omega = this.variables.get("omega0").re, x2 = this.variables.get("x2").re;
    const DPM = this.variables.get("DPM_momentum").re;
    
    const F_grav = this.G * M * M / (r * r);
    const F_dpm_re = DPM * F_grav * Math.cos(omega * t) * x2;
    const F_vel = 0.5 * M * V * V / r;
    const F_xray = L_X / (this.c * 4 * Math.PI * r * r);
    
    let F_total = {re: F_grav + F_vel + F_xray + F_dpm_re, im: 0};
    for (const term of this.dynamicTerms) F_total = complexAdd(F_total, term.compute(this.variables));
    return F_total;
  }
  
  computeBuoyancy(t) {
    const rho = this.variables.get("rho_gas").re, T = this.variables.get("T_val").re, r = this.variables.get("r").re;
    return {re: rho * this.k_B * T / this.m_n * (4/3) * Math.PI * r**3, im: 0};
  }
  
  getAvailableSystems() { return Object.keys(this.systemConfigs); }
  clone() { return new UQFFBuoyancyModule161(this.currentSystem); }
  updateVariable(n, v) { if (this.variables.has(n)) this.variables.set(n, {re:v, im:0}); }
  setEnableLogging(e) { this.enableLogging = e; }
  registerDynamicTerm(t) { this.dynamicTerms.push(t); }
}

const domainExpansion = {
  expandQuasarScale(qF, lF) {
    this.variables.get("M").re *= qF;
    this.variables.get("L_X").re *= lF;
  },
  expandClusterScale(cF, icmF) {
    this.variables.get("M").re *= cF;
    this.variables.get("rho_gas").re *= icmF;
    this.variables.get("T_val").re *= Math.pow(cF/icmF, 1/3);
  },
  expandMultiSystemScale(sysName, mF, rF) {
    this.setSystem(sysName);
    this.variables.get("M").re *= mF;
    this.variables.get("r").re *= rF;
    this.variables.get("rho_gas").re *= mF / (rF**3);
  },
  optimizeForMetric(m) {
    const p = {"standard":{DPM_momentum:0.93,k_LENR:1e-10},"simulation_fast":{DPM_momentum:1.0,k_LENR:1e-10},"quasar_mode":{omega0:1e-14,x2:-1.35e172},"cluster_mode":{omega0:1e-15,x2:-2.27e172}}[m];
    if (p) Object.entries(p).forEach(([k,v]) => this.variables.has(k) && (this.variables.get(k).re = v));
  }
};

addEnhancedDynamics(UQFFBuoyancyModule161, "Multi_System_161", domainExpansion);
module.exports = UQFFBuoyancyModule161;
