// Source156.js - UQFFBuoyancyCNBModule
// Multi-system UQFF buoyancy with Cosmic Neutrino Background (CNB) integration
// Enhanced with full 25-method self-expansion framework
// Systems: J1610+1811, PLCK G287.0+32.9, PSZ2 G181.06+48.47, ASKAP J1832-0911, Sonification Collection, Centaurus A
// Copyright - Daniel T. Murphy, Nov 2025

const { addEnhancedDynamics } = require('./enhanced_dynamics.js');

// Complex number helpers
function complexAdd(a, b) { return {re: a.re + b.re, im: a.im + b.im}; }
function complexSub(a, b) { return {re: a.re - b.re, im: a.im - b.im}; }
function complexMul(a, b) { return {re: a.re*b.re - a.im*b.im, im: a.re*b.im + a.im*b.re}; }
function complexScale(a, s) { return {re: a.re*s, im: a.im*s}; }

class UQFFBuoyancyCNBModule {
  constructor(systemName = "J1610") {
    this.systemConfigs = {
      "J1610": { M: 2.785e30, r: 3.09e15, L_X: 1e30, B0: 1e-9, rho_gas: 1e-20, omega0: 1e-12, x2: -1.35e172 },
      "PLCK_G287": { M: 1.989e44, r: 3.09e22, L_X: 1e44, B0: 1e-9, rho_gas: 1e-24, omega0: 1e-15, x2: -1.35e172 },
      "PSZ2_G181": { M: 1.989e44, r: 3.09e22, L_X: 1e44, B0: 1e-9, rho_gas: 1e-24, omega0: 1e-15, x2: -1.35e172 },
      "ASKAP_J1832": { M: 2.785e30, r: 4.63e16, L_X: 1e30, B0: 1e-8, rho_gas: 1e-21, omega0: 1e-12, x2: -3.40e172 },
      "Sonification": { M: 1.989e31, r: 6.17e16, L_X: 1e31, B0: 1e-8, rho_gas: 1e-21, omega0: 1e-12, x2: -3.40e172 },
      "CentaurusA": { M: 1.094e38, r: 6.17e17, L_X: 1e37, B0: 1e-9, rho_gas: 1e-22, omega0: 1e-14, x2: -1.35e172 }
    };
    
    this.currentSystem = systemName;
    const config = this.systemConfigs[systemName] || this.systemConfigs["J1610"];
    Object.assign(this, config);
    
    this.T_val = 1e7; this.V = 1e5; this.n_e = 1e6;
    this.k_DE = 1e-16; this.k_act = 1e-14;
    this.DPM_momentum = 0.93; this.DPM_gravity = 1.0; this.DPM_stability = 0.01;
    this.k_LENR = 1e-10; this.beta_i = 0.6;
    this.F_CNB = 9.07e-42; // Cosmic Neutrino Background force
    
    this.G = 6.67430e-11; this.c = 299792458.0; this.hbar = 1.054571817e-34;
    this.m_n = 1.674927498e-27; this.k_B = 1.380649e-23;
    
    this.variables = new Map();
    this.initializeVariables();
    this.dynamicTerms = []; this.dynamicParameters = new Map();
    
    this.metadata = new Map();
    this.metadata.set("enhanced", true);
    this.metadata.set("version", "2.0.0");
    this.metadata.set("object_type", "multi_system_buoyancy_cnb");
    this.metadata.set("system_name", systemName);
    this.enableLogging = false;
  }
  
  initializeVariables() {
    ["M", "r", "L_X", "B0", "omega0", "rho_gas", "T_val", "V", "n_e", "k_DE", "k_act",
     "DPM_momentum", "DPM_gravity", "DPM_stability", "k_LENR", "x2", "beta_i", "F_CNB"].forEach(k => {
      this.variables.set(k, {re: this[k], im: 0});
    });
    this.variables.set("t", {re: 0, im: 0});
  }
  
  updateVariable(name, value) { if (this.variables.has(name)) this.variables.set(name, {re: value, im: 0}); }
  addToVariable(name, delta) {
    if (this.variables.has(name)) {
      const c = this.variables.get(name);
      this.variables.set(name, {re: c.re + delta, im: c.im});
    }
  }
  
  setSystem(systemName) {
    if (this.systemConfigs[systemName]) {
      this.currentSystem = systemName;
      Object.entries(this.systemConfigs[systemName]).forEach(([k, v]) => {
        this[k] = v;
        if (this.variables.has(k)) this.variables.get(k).re = v;
      });
      this.metadata.set("system_name", systemName);
    }
  }
  
  computeDPM_resonance(t) {
    const M = this.variables.get("M").re, r = this.variables.get("r").re;
    const omega = this.variables.get("omega0").re, DPM = this.variables.get("DPM_momentum").re;
    const F_grav = this.G * M * M / (r * r), x2 = this.variables.get("x2").re;
    return {re: DPM * F_grav * Math.cos(omega * t) * x2, im: 0};
  }
  
  computeLENRTerm() {
    const k = this.variables.get("k_LENR").re, omega_L = 2 * Math.PI * 1.25e12;
    const omega0 = this.variables.get("omega0").re;
    return {re: k * Math.pow(omega_L / omega0, 2), im: 0};
  }
  
  computeF(t) {
    this.variables.set("t", {re: t, im: 0});
    const M = this.variables.get("M").re, r = this.variables.get("r").re;
    const L_X = this.variables.get("L_X").re, F_CNB = this.variables.get("F_CNB").re;
    
    const F_grav = this.G * M * M / (r * r);
    const F_dpm = this.computeDPM_resonance(t);
    const F_lenr = this.computeLENRTerm();
    const F_xray = L_X / (this.c * 4 * Math.PI * r * r);
    
    let F_total = {re: F_grav + F_xray + F_CNB, im: 0};
    F_total = complexAdd(F_total, F_dpm);
    F_total = complexAdd(F_total, F_lenr);
    
    this.dynamicTerms.forEach(term => { F_total = complexAdd(F_total, term.compute(this.variables)); });
    return F_total;
  }
  
  computeBuoyancy(t) {
    const rho = this.variables.get("rho_gas").re, T = this.variables.get("T_val").re;
    const r = this.variables.get("r").re, beta = this.variables.get("beta_i").re;
    const pressure = rho * this.k_B * T / this.m_n;
    const volume = (4/3) * Math.PI * r * r * r;
    return {re: beta * pressure * volume, im: 0};
  }
  
  setEnableLogging(enable) { this.enableLogging = enable; }
  registerDynamicTerm(term) { this.dynamicTerms.push(term); }
  setDynamicParameter(name, value) { this.dynamicParameters.set(name, value); }
  getDynamicParameter(name) { return this.dynamicParameters.get(name); }
  
  // Clone for parallel processing
  clone() {
    const cloned = new UQFFBuoyancyCNBModule(this.currentSystem);
    cloned.variables = new Map(this.variables);
    cloned.dynamicParameters = new Map(this.dynamicParameters);
    cloned.metadata = new Map(this.metadata);
    cloned.enableLogging = this.enableLogging;
    return cloned;
  }
}

const domainExpansion = {
  expandSystemScale(massFactor, radiusFactor) {
    const M = this.variables.get("M").re, r = this.variables.get("r").re, rho = this.variables.get("rho_gas").re;
    this.variables.get("M").re = M * massFactor;
    this.variables.get("r").re = r * radiusFactor;
    this.variables.get("rho_gas").re = rho * massFactor / (radiusFactor ** 3);
    if (this.enableLogging) console.log(`Expanded system: M×${massFactor}, r×${radiusFactor}`);
  },
  expandBuoyancyScale(buoyancyFactor, CNBFactor) {
    const beta = this.variables.get("beta_i").re, F_CNB = this.variables.get("F_CNB").re;
    this.variables.get("beta_i").re = beta * buoyancyFactor;
    this.variables.get("F_CNB").re = F_CNB * CNBFactor;
    if (this.enableLogging) console.log(`Expanded buoyancy: β×${buoyancyFactor}, CNB×${CNBFactor}`);
  },
  expandMultiSystemScale(factor) {
    Object.keys(this.systemConfigs).forEach(s => {
      this.systemConfigs[s].M *= factor;
      this.systemConfigs[s].r *= Math.sqrt(factor);
    });
    if (this.enableLogging) console.log(`Scaled all systems: ×${factor}`);
  },
  optimizeForMetric(metricName) {
    const presets = {
      "standard": { DPM_momentum: 0.93, k_LENR: 1e-10, F_CNB: 9.07e-42 },
      "high_energy": { DPM_momentum: 1.5, k_LENR: 5e-10, F_CNB: 5e-41 },
      "quiescent": { DPM_momentum: 0.5, k_LENR: 1e-11, F_CNB: 1e-42 }
    };
    if (presets[metricName]) {
      Object.entries(presets[metricName]).forEach(([k, v]) => {
        if (this.variables.has(k)) this.variables.get(k).re = v;
      });
      if (this.enableLogging) console.log(`Optimized: ${metricName}`);
    }
  }
};

addEnhancedDynamics(UQFFBuoyancyCNBModule, "Multi_System_CNB", domainExpansion);
module.exports = UQFFBuoyancyCNBModule;
