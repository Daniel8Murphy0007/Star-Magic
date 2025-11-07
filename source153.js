// source153.js - Abell2256UQFFModule153
// Abell 2256 galaxy cluster (alternative implementation with full buoyancy)
// Enhanced with full 25-method self-expansion framework
// M=1.23e45 kg, r=3.93e22 m, omega0=1e-15, x2=-1.35e172, E_cm=2.18e-6 J (13.6 TeV)
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

class Abell2256UQFFModule153 {
  constructor() {
    // Core physics parameters from C++ source153.cpp
    this.M = 1.23e45;            // Total cluster mass (kg)
    this.r = 3.93e22;            // Cluster radius (m)
    this.L_X = 3.7e37;           // X-ray luminosity (W)
    this.B0 = 1e-9;              // Magnetic field (T)
    this.omega0 = 1e-15;         // Angular frequency (rad/s)
    this.rho_gas = 5e-24;        // ICM density (kg/m^3)
    this.V = 1e-3;               // Particle velocity (m/s)
    this.T_val = 1e8;            // Temperature (K) - very hot ICM
    this.n_e = 1e3;              // Electron density (m^-3)
    this.k_DE = 1e-30;           // Dark energy coupling
    this.k_act = 1e-6;           // Active coupling
    this.x2 = -1.35e172;         // UQFF coupling constant
    
    // UQFF parameters
    this.DPM_momentum = 0.93;
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
    
    // Initialize variable storage
    this.variables = new Map();
    this.variables.set("M", {re: this.M, im: 0});
    this.variables.set("r", {re: this.r, im: 0});
    this.variables.set("L_X", {re: this.L_X, im: 0});
    this.variables.set("B0", {re: this.B0, im: 0});
    this.variables.set("omega0", {re: this.omega0, im: 0});
    this.variables.set("rho_gas", {re: this.rho_gas, im: 0});
    this.variables.set("V", {re: this.V, im: 0});
    this.variables.set("T_val", {re: this.T_val, im: 0});
    this.variables.set("n_e", {re: this.n_e, im: 0});
    this.variables.set("k_DE", {re: this.k_DE, im: 0});
    this.variables.set("k_act", {re: this.k_act, im: 0});
    this.variables.set("DPM_momentum", {re: this.DPM_momentum, im: 0});
    this.variables.set("k_LENR", {re: this.k_LENR, im: 0});
    this.variables.set("t", {re: 0, im: 0});
    
    this.dynamicTerms = [];
    this.dynamicParameters = new Map();
    
    this.metadata = new Map();
    this.metadata.set("enhanced", true);
    this.metadata.set("version", "2.0.0");
    this.metadata.set("object_type", "galaxy_cluster");
    this.metadata.set("system_name", "Abell_2256_v2");
    
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
    const V_val = this.variables.get("V").re;
    const L_X_val = this.variables.get("L_X").re;
    const k_DE_val = this.variables.get("k_DE").re;
    
    // Gravitational force
    const F_grav = this.G * M_val * M_val / (r_val * r_val);
    
    // DPM resonance
    const F_dpm = this.computeDPM_resonance(t);
    
    // LENR contribution
    const F_lenr = this.computeLENRTerm(t);
    
    // ICM pressure
    const F_icm = 0.5 * M_val * V_val * V_val / r_val;
    
    // X-ray radiation pressure
    const F_xray = L_X_val / (this.c * 4.0 * Math.PI * r_val * r_val);
    
    // Dark energy term
    const F_DE = k_DE_val * r_val * r_val;
    
    // Total force
    let F_total = {re: F_grav + F_icm + F_xray + F_DE, im: 0};
    F_total = complexAdd(F_total, F_dpm);
    F_total = complexAdd(F_total, F_lenr);
    
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
}

// ========== Domain-Specific Expansion Methods ==========
const domainExpansion = {
  expandClusterScale(massFactor, radiusFactor) {
    const M_val = this.variables.get("M").re;
    const r_val = this.variables.get("r").re;
    const rho_val = this.variables.get("rho_gas").re;
    const L_X_val = this.variables.get("L_X").re;
    
    this.variables.get("M").re = M_val * massFactor;
    this.variables.get("r").re = r_val * radiusFactor;
    this.variables.get("rho_gas").re = rho_val * massFactor / (radiusFactor ** 3);
    this.variables.get("L_X").re = L_X_val * massFactor * massFactor;
    
    if (this.enableLogging) {
      console.log(`Expanded cluster: M×${massFactor}, r×${radiusFactor}`);
    }
  },
  
  expandICMScale(temperatureFactor, densityFactor) {
    const T_val = this.variables.get("T_val").re;
    const rho_val = this.variables.get("rho_gas").re;
    const V_val = this.variables.get("V").re;
    const L_X_val = this.variables.get("L_X").re;
    
    this.variables.get("T_val").re = T_val * temperatureFactor;
    this.variables.get("rho_gas").re = rho_val * densityFactor;
    this.variables.get("V").re = V_val * Math.sqrt(temperatureFactor);
    this.variables.get("L_X").re = L_X_val * densityFactor * densityFactor * 
                                     Math.sqrt(temperatureFactor);
    
    if (this.enableLogging) {
      console.log(`Expanded ICM: T×${temperatureFactor}, ρ×${densityFactor}`);
    }
  },
  
  expandMergerScale(activityFactor, magneticFactor) {
    const k_act_val = this.variables.get("k_act").re;
    const B0_val = this.variables.get("B0").re;
    const k_DE_val = this.variables.get("k_DE").re;
    
    this.variables.get("k_act").re = k_act_val * activityFactor;
    this.variables.get("B0").re = B0_val * magneticFactor;
    this.variables.get("k_DE").re = k_DE_val * activityFactor;
    
    if (this.enableLogging) {
      console.log(`Expanded merger: activity×${activityFactor}, B×${magneticFactor}`);
    }
  },
  
  optimizeForMetric(metricName) {
    const presets = {
      "standard": {
        DPM_momentum: 0.93,
        k_LENR: 1e-10,
        k_DE: 1e-30
      },
      "merger": {
        DPM_momentum: 1.3,
        k_LENR: 5e-10,
        k_DE: 5e-30,
        k_act: 5e-6,
        T_val: 2e8,
        V: 1e-2
      },
      "relaxed": {
        DPM_momentum: 0.6,
        k_LENR: 2e-11,
        k_DE: 2e-31,
        k_act: 1e-7,
        T_val: 5e7
      },
      "cool_core": {
        DPM_momentum: 0.75,
        k_LENR: 5e-11,
        T_val: 3e7,
        rho_gas: 1e-23,
        L_X: 1e38
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

addEnhancedDynamics(Abell2256UQFFModule153, "Abell_2256_v2", domainExpansion);

module.exports = Abell2256UQFFModule153;
