// source147.js - NGC2207UQFFModule
// Interacting galaxy system NGC 2207 and IC 2163
// Enhanced with full 25-method self-expansion framework
// M=3.978e40 kg, r=4.40e20 m, omega0=1e-12, x2=-1.35e172
// Copyright - Daniel T. Murphy, Nov 2025

const { addEnhancedDynamics } = require('./enhanced_dynamics.js');

// Complex number helpers
function complexAdd(a, b) { return {re: a.re + b.re, im: a.im + b.im}; }
function complexSub(a, b) { return {re: a.re - b.re, im: a.im - b.im}; }
function complexMul(a, b) { return {re: a.re*b.re - a.im*b.im, im: a.re*b.im + a.im*b.re}; }
function complexDiv(a, b) {
  const denom = b.re*b.re + b.im*b.im;
  return {re: (a.re*b.re + a.im*b.im)/denom, im: (a.im*b.re - a.re*b.im)/denom};
}
function complexPow(base, exponent) {
  const r = Math.sqrt(base.re*base.re + base.im*base.im);
  const theta = Math.atan2(base.im, base.re);
  const newR = Math.pow(r, exponent);
  const newTheta = theta * exponent;
  return {re: newR * Math.cos(newTheta), im: newR * Math.sin(newTheta)};
}
function complexScale(a, s) { return {re: a.re*s, im: a.im*s}; }
function complexNeg(a) { return {re: -a.re, im: -a.im}; }

class NGC2207UQFFModule {
  constructor() {
    // Core physics parameters from C++ source147.cpp
    this.M = 3.978e40;           // Total mass (kg) - interacting galaxy system
    this.r = 4.40e20;            // Separation distance (m)
    this.rho_gas = 1e-21;        // Gas density (kg/m^3)
    this.L_X = 1e41;             // X-ray luminosity (W)
    this.V = 1e5;                // Tidal velocity (m/s) - 100 km/s
    this.B0 = 1e-10;             // Magnetic field (T)
    this.n_e = 1e6;              // Electron density (m^-3)
    this.omega0 = 1e-12;         // Angular frequency (rad/s) - galactic
    this.T_val = 1e7;            // Temperature (K) - hot gas
    this.k_DE = 1e-16;           // Dark energy coupling (J/m^3)
    this.k_act = 1e-14;          // Active galaxy coupling
    this.x2 = -1.35e172;         // UQFF coupling constant
    
    // UQFF parameters
    this.DPM_momentum = 0.25;
    this.k_LENR = 1e-10;
    this.E_cm = 2.18e-6;         // 13.6 TeV
    this.epsilon = 0.01;
    
    // Constants
    this.G = 6.67430e-11;
    this.c = 299792458.0;
    this.hbar = 1.054571817e-34;
    this.m_n = 1.674927498e-27;
    this.m_e = 9.1093837015e-31;
    this.k_B = 1.380649e-23;
    this.mu0 = 1.25663706212e-6;
    
    // Initialize variable storage (for enhanced dynamics)
    this.variables = new Map();
    this.variables.set("M", {re: this.M, im: 0});
    this.variables.set("r", {re: this.r, im: 0});
    this.variables.set("rho_gas", {re: this.rho_gas, im: 0});
    this.variables.set("L_X", {re: this.L_X, im: 0});
    this.variables.set("V", {re: this.V, im: 0});
    this.variables.set("B0", {re: this.B0, im: 0});
    this.variables.set("n_e", {re: this.n_e, im: 0});
    this.variables.set("omega0", {re: this.omega0, im: 0});
    this.variables.set("T_val", {re: this.T_val, im: 0});
    this.variables.set("k_DE", {re: this.k_DE, im: 0});
    this.variables.set("k_act", {re: this.k_act, im: 0});
    this.variables.set("DPM_momentum", {re: this.DPM_momentum, im: 0});
    this.variables.set("k_LENR", {re: this.k_LENR, im: 0});
    this.variables.set("t", {re: 0, im: 0});
    
    // Dynamic terms and parameters
    this.dynamicTerms = [];
    this.dynamicParameters = new Map();
    
    // Metadata
    this.metadata = new Map();
    this.metadata.set("enhanced", true);
    this.metadata.set("version", "2.0.0");
    this.metadata.set("object_type", "interacting_galaxies");
    this.metadata.set("system_name", "NGC2207_IC2163");
    
    // Logging
    this.enableLogging = false;
  }
  
  // ========== Core Physics Methods ==========
  
  updateVariable(name, value) {
    if (this.variables.has(name)) {
      this.variables.set(name, {re: value, im: 0});
    }
  }
  
  addToVariable(name, delta) {
    if (this.variables.has(name)) {
      const current = this.variables.get(name);
      this.variables.set(name, {re: current.re + delta, im: current.im});
    }
  }
  
  computeDPM_resonance(t) {
    const M_val = this.variables.get("M").re;
    const r_val = this.variables.get("r").re;
    const omega = this.variables.get("omega0").re;
    
    const F_grav = this.G * M_val * M_val / (r_val * r_val);
    const osc = Math.cos(omega * t);
    return {re: this.DPM_momentum * F_grav * osc * this.x2, im: 0};
  }
  
  computeLENRTerm(t) {
    const T = this.variables.get("T_val").re;
    const rho = this.variables.get("rho_gas").re;
    
    const thermal = this.k_B * T;
    const density_factor = rho / this.m_n;
    return {re: this.k_LENR * thermal * density_factor, im: 0};
  }
  
  computeIntegrand(t) {
    const F_dpm = this.computeDPM_resonance(t);
    const F_lenr = this.computeLENRTerm(t);
    const sum = complexAdd(F_dpm, F_lenr);
    
    const omega = this.variables.get("omega0").re;
    const phase = {re: Math.cos(omega * t), im: Math.sin(omega * t)};
    return complexMul(sum, phase);
  }
  
  computeF(t) {
    this.variables.set("t", {re: t, im: 0});
    
    const M_val = this.variables.get("M").re;
    const r_val = this.variables.get("r").re;
    const L_X_val = this.variables.get("L_X").re;
    const V_val = this.variables.get("V").re;
    const k_DE_val = this.variables.get("k_DE").re;
    const k_act_val = this.variables.get("k_act").re;
    
    // Gravitational force
    const F_grav = this.G * M_val * M_val / (r_val * r_val);
    
    // DPM resonance
    const F_dpm = this.computeDPM_resonance(t);
    
    // LENR contribution
    const F_lenr = this.computeLENRTerm(t);
    
    // X-ray luminosity pressure
    const F_xray = L_X_val / (this.c * 4.0 * Math.PI * r_val * r_val);
    
    // Tidal interaction
    const F_tidal = 0.5 * M_val * V_val * V_val / r_val;
    
    // Dark energy term
    const F_DE = k_DE_val * r_val * r_val;
    
    // Active nucleus contribution
    const F_active = k_act_val * Math.sqrt(L_X_val) / r_val;
    
    // Total force
    let F_total = {re: F_grav + F_xray + F_tidal + F_DE + F_active, im: 0};
    F_total = complexAdd(F_total, F_dpm);
    F_total = complexAdd(F_total, F_lenr);
    
    // Add dynamic terms
    this.dynamicTerms.forEach(term => {
      F_total = complexAdd(F_total, term.compute(this.variables));
    });
    
    return F_total;
  }
  
  computeBuoyancy(t) {
    const rho = this.variables.get("rho_gas").re;
    const T = this.variables.get("T_val").re;
    const r_val = this.variables.get("r").re;
    
    const pressure = rho * this.k_B * T / this.m_n;
    const volume = (4.0/3.0) * Math.PI * r_val * r_val * r_val;
    return {re: pressure * volume, im: 0};
  }
  
  computeSuperconductive(t) {
    const B = this.variables.get("B0").re;
    const n = this.variables.get("n_e").re;
    
    const lambda_London = Math.sqrt(this.m_e / (this.mu0 * n * Math.pow(1.602e-19, 2)));
    const penetration = Math.exp(-this.r / lambda_London);
    return {re: B * B * penetration / (2.0 * this.mu0), im: 0};
  }
  
  setEnableLogging(enable) {
    this.enableLogging = enable;
  }
  
  registerDynamicTerm(term) {
    this.dynamicTerms.push(term);
  }
  
  setDynamicParameter(name, value) {
    this.dynamicParameters.set(name, value);
  }
  
  getDynamicParameter(name) {
    return this.dynamicParameters.get(name);
  }
  
  // Clone for parallel processing
  clone() {
    const cloned = new NGC2207UQFFModule();
    cloned.variables = new Map(this.variables);
    cloned.dynamicParameters = new Map(this.dynamicParameters);
    cloned.metadata = new Map(this.metadata);
    cloned.enableDynamicTerms = this.enableDynamicTerms;
    cloned.enableLogging = this.enableLogging;
    cloned.learningRate = this.learningRate;
    return cloned;
  }
}

// ========== Domain-Specific Expansion Methods ==========
const domainExpansion = {
  expandGalaxyScale(massFactor, radiusFactor) {
    const M_val = this.variables.get("M").re;
    const r_val = this.variables.get("r").re;
    const rho_val = this.variables.get("rho_gas").re;
    const L_X_val = this.variables.get("L_X").re;
    
    this.variables.get("M").re = M_val * massFactor;
    this.variables.get("r").re = r_val * radiusFactor;
    
    // Scaling relationships
    this.variables.get("rho_gas").re = rho_val * massFactor / (radiusFactor * radiusFactor * radiusFactor);
    this.variables.get("L_X").re = L_X_val * massFactor * massFactor;
    
    if (this.enableLogging) {
      console.log(`Expanded galaxy scale: M×${massFactor}, r×${radiusFactor}`);
    }
  },
  
  expandForceScale(dpmFactor, lenrFactor) {
    this.variables.get("DPM_momentum").re *= dpmFactor;
    this.variables.get("k_LENR").re *= lenrFactor;
    
    if (this.enableLogging) {
      console.log(`Expanded force scale: DPM×${dpmFactor}, LENR×${lenrFactor}`);
    }
  },
  
  expandInteractionScale(tidalFactor, luminosityFactor) {
    const V_val = this.variables.get("V").re;
    const L_X_val = this.variables.get("L_X").re;
    const k_DE_val = this.variables.get("k_DE").re;
    const k_act_val = this.variables.get("k_act").re;
    const B0_val = this.variables.get("B0").re;
    
    this.variables.get("V").re = V_val * tidalFactor;
    this.variables.get("L_X").re = L_X_val * luminosityFactor;
    this.variables.get("k_DE").re = k_DE_val * tidalFactor;
    this.variables.get("k_act").re = k_act_val * Math.sqrt(luminosityFactor);
    this.variables.get("B0").re = B0_val * Math.sqrt(luminosityFactor);
    
    if (this.enableLogging) {
      console.log(`Expanded interaction scale: tidal×${tidalFactor}, L×${luminosityFactor}`);
    }
  },
  
  optimizeForMetric(metricName) {
    const presets = {
      "standard_ngc2207": {
        DPM_momentum: 0.25,
        k_LENR: 1e-10,
        k_DE: 1e-16,
        k_act: 1e-14
      },
      "close_approach": {
        DPM_momentum: 0.50,
        k_LENR: 2e-10,
        k_DE: 5e-16,
        k_act: 5e-14,
        V: 2e5
      },
      "starburst_phase": {
        DPM_momentum: 0.75,
        k_LENR: 5e-10,
        k_DE: 1e-15,
        k_act: 1e-13,
        L_X: 5e41,
        rho_gas: 5e-21
      },
      "tidal_arms": {
        DPM_momentum: 0.40,
        k_LENR: 3e-10,
        k_DE: 3e-16,
        k_act: 3e-14,
        V: 1.5e5
      },
      "quiescent": {
        DPM_momentum: 0.10,
        k_LENR: 5e-11,
        k_DE: 5e-17,
        k_act: 5e-15,
        V: 5e4,
        L_X: 5e40
      }
    };
    
    if (presets[metricName]) {
      Object.entries(presets[metricName]).forEach(([key, value]) => {
        if (this.variables.has(key)) {
          this.variables.get(key).re = value;
        }
      });
      if (this.enableLogging) {
        console.log(`Optimized for metric: ${metricName}`);
      }
    }
  }
};

// Apply enhanced dynamics
addEnhancedDynamics(NGC2207UQFFModule, "NGC2207_IC2163", domainExpansion);

module.exports = NGC2207UQFFModule;
