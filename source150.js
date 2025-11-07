// source150.js - SPTCLJ2215UQFFModule
// SPT-CL J2215-3537 massive galaxy cluster
// Enhanced with full 25-method self-expansion framework
// M=1.46e45 kg, r=3.09e22 m, x2=-2.27e172
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

class SPTCLJ2215UQFFModule {
  constructor() {
    // Core physics parameters from C++ source150.cpp
    this.M = 1.46e45;            // Total cluster mass (kg) - one of most massive known
    this.r = 3.09e22;            // Cluster radius (m) - ~1 Mpc
    this.rho_gas = 1e-23;        // ICM density (kg/m^3)
    this.L_X = 1e45;             // X-ray luminosity (W) - very luminous ICM
    this.V = 1e6;                // Velocity dispersion (m/s) - 1000 km/s
    this.B0 = 1e-9;              // Magnetic field (T)
    this.n_e = 1e3;              // Electron density (m^-3)
    this.omega0 = 1e-15;         // Angular frequency (rad/s) - cluster scale
    this.T_val = 1e8;            // Temperature (K) - very hot ICM
    this.k_DE = 1e-14;           // Dark energy coupling (J/m^3)
    this.k_act = 1e-12;          // Active coupling
    this.x2 = -2.27e172;         // UQFF coupling constant (massive clusters)
    
    // UQFF parameters
    this.DPM_momentum = 0.35;
    this.k_LENR = 1e-10;
    this.E_cm = 2.18e-6;
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
    
    this.dynamicTerms = [];
    this.dynamicParameters = new Map();
    
    this.metadata = new Map();
    this.metadata.set("enhanced", true);
    this.metadata.set("version", "2.0.0");
    this.metadata.set("object_type", "galaxy_cluster");
    this.metadata.set("system_name", "SPT_CL_J2215_3537");
    
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
    
    // Dark energy term (significant at cluster scales)
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
    
    // Scaling relationships for clusters
    this.variables.get("rho_gas").re = rho_val * massFactor / (radiusFactor * radiusFactor * radiusFactor);
    
    // L_X scales as M^2 * T^0.5 for clusters
    this.variables.get("L_X").re = L_X_val * massFactor * massFactor;
    
    if (this.enableLogging) {
      console.log(`Expanded cluster scale: M×${massFactor}, r×${radiusFactor}`);
    }
  },
  
  expandICMScale(temperatureFactor, densityFactor) {
    const T_val = this.variables.get("T_val").re;
    const rho_val = this.variables.get("rho_gas").re;
    const L_X_val = this.variables.get("L_X").re;
    const V_val = this.variables.get("V").re;
    
    this.variables.get("T_val").re = T_val * temperatureFactor;
    this.variables.get("rho_gas").re = rho_val * densityFactor;
    
    // Velocity dispersion scales with sqrt(T)
    this.variables.get("V").re = V_val * Math.sqrt(temperatureFactor);
    
    // X-ray luminosity scales with rho^2 * sqrt(T)
    this.variables.get("L_X").re = L_X_val * densityFactor * densityFactor * Math.sqrt(temperatureFactor);
    
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
        DPM_momentum: 0.35,
        k_LENR: 1e-10,
        k_DE: 1e-14
      },
      "merger": {
        DPM_momentum: 0.70,
        k_LENR: 5e-10,
        k_DE: 5e-14,
        k_act: 5e-12,
        T_val: 2e8,
        V: 2e6
      },
      "relaxed": {
        DPM_momentum: 0.15,
        k_LENR: 2e-11,
        k_DE: 2e-15,
        k_act: 1e-13,
        T_val: 5e7
      },
      "cool_core": {
        DPM_momentum: 0.25,
        k_LENR: 5e-11,
        T_val: 3e7,
        rho_gas: 1e-22,
        L_X: 5e44
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

addEnhancedDynamics(SPTCLJ2215UQFFModule, "SPT_CL_J2215_3537", domainExpansion);

module.exports = SPTCLJ2215UQFFModule;
