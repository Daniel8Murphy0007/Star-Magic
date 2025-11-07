// Source157.js - UQFFBuoyancyModule157
// Multi-system buoyancy: M104, NGC 4839, Chandra/Webb, NGC 346, NGC 1672
// SIMULATION-READY: Thread-safe, parallel computing, fast computation paths
// Copyright - Daniel T. Murphy, Nov 2025

const { addEnhancedDynamics } = require('./enhanced_dynamics.js');

function complexAdd(a, b) { return {re: a.re + b.re, im: a.im + b.im}; }
function complexScale(a, s) { return {re: a.re*s, im: a.im*s}; }

class UQFFBuoyancyModule157 {
  constructor(systemName = "M104") {
    this.systemConfigs = {
      "M104": { M: 1.59e41, r: 1.23e21, L_X: 1e39, B0: 1e-9, rho_gas: 1e-21, T_val: 1e7, omega0: 1e-14, x2: -1.35e172 },
      "NGC4839": { M: 7.96e44, r: 3.09e22, L_X: 1e44, B0: 1e-9, rho_gas: 1e-23, T_val: 1e8, omega0: 1e-15, x2: -2.27e172 },
      "Chandra_Webb": { M: 3.978e40, r: 6.17e21, L_X: 1e38, B0: 1e-9, rho_gas: 1e-22, T_val: 5e7, omega0: 1e-15, x2: -1.35e172 },
      "NGC346": { M: 1.989e35, r: 6.17e18, L_X: 1e33, B0: 1e-8, rho_gas: 1e-19, T_val: 1e4, omega0: 1e-10, x2: -3.40e172 },
      "NGC1672": { M: 5.97e40, r: 1.54e21, L_X: 1e38, B0: 1e-9, rho_gas: 1e-21, T_val: 1e7, omega0: 1e-14, x2: -1.35e172 }
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
    this.variables.set("n_e", {re: 1e3, im: 0});
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
    const P = rho * this.k_B * T / this.m_n;
    return {re: P * (4/3) * Math.PI * r**3, im: 0};
  }
  
  getAvailableSystems() { return Object.keys(this.systemConfigs); }
  clone() { return new UQFFBuoyancyModule157(this.currentSystem); }
  updateVariable(n, v) { if (this.variables.has(n)) this.variables.set(n, {re:v, im:0}); }
  setEnableLogging(e) { this.enableLogging = e; }
  registerDynamicTerm(t) { this.dynamicTerms.push(t); }
}

const domainExpansion = {
  expandMultiSystemScale(sysName, mF, rF) {
    this.setSystem(sysName);
    const M = this.variables.get("M").re, r = this.variables.get("r").re;
    this.variables.get("M").re = M * mF;
    this.variables.get("r").re = r * rF;
    this.variables.get("rho_gas").re *= mF / (rF**3);
  },
  expandBuoyancyScale(pF, tF) {
    this.variables.get("T_val").re *= tF;
    this.variables.get("rho_gas").re *= pF / tF;
  },
  expandDynamicsScale(dF, lF) {
    this.variables.get("DPM_momentum").re *= dF;
    this.variables.get("k_LENR").re *= lF;
  },
  optimizeForMetric(m) {
    const p = {"standard":{DPM_momentum:0.93,k_LENR:1e-10},"simulation_fast":{DPM_momentum:1.0,k_LENR:1e-10}}[m];
    if (p) Object.entries(p).forEach(([k,v]) => this.variables.has(k) && (this.variables.get(k).re = v));
  }
};

addEnhancedDynamics(UQFFBuoyancyModule157, "Multi_System_157", domainExpansion);
module.exports = UQFFBuoyancyModule157;
