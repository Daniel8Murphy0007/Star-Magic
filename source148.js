// source148.js - RAquariiUQFFModule
// R Aquarii symbiotic star system with stellar wind outflows
// Enhanced with full 25-method self-expansion framework
// M=3.978e30 kg, r=2.18e15 m, V=1e5 m/s (100 km/s), x2=-3.40e172
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

class RAquariiUQFFModule {
  constructor() {
    // Core physics parameters from C++ source148.cpp
    this.M = 3.978e30;           // Total mass (kg) - symbiotic binary system
    this.r = 2.18e15;            // Separation distance (m) - ~14.5 AU
    this.V = 1e5;                // Wind velocity (m/s) - 100 km/s
    this.rho_gas = 1e-18;        // Gas density (kg/m^3)
    this.L_X = 1e30;             // X-ray luminosity (W)
    this.B0 = 1e-8;              // Magnetic field (T)
    this.n_e = 1e9;              // Electron density (m^-3)
    this.omega0 = 1e-8;          // Angular frequency (rad/s) - binary orbital
    this.T_val = 1e5;            // Temperature (K) - stellar wind
    this.k_DE = 1e-20;           // Dark energy coupling (J/m^3)
    this.k_act = 1e-18;          // Active coupling
    this.x2 = -3.40e172;         // UQFF coupling constant (symbiotic systems)
    
    // UQFF parameters
    this.DPM_momentum = 0.20;
    this.k_LENR = 1e-11;
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
    this.variables.set("V", {re: this.V, im: 0});
    this.variables.set("rho_gas", {re: this.rho_gas, im: 0});
    this.variables.set("L_X", {re: this.L_X, im: 0});
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
    this.metadata.set("object_type", "symbiotic_star");
    this.metadata.set("system_name", "R_Aquarii");
    
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
    const rho_val = this.variables.get("rho_gas").re;
    const L_X_val = this.variables.get("L_X").re;
    
    // Gravitational force
    const F_grav = this.G * M_val * M_val / (r_val * r_val);
    
    // DPM resonance
    const F_dpm = this.computeDPM_resonance(t);
    
    // LENR contribution
    const F_lenr = this.computeLENRTerm(t);
    
    // Wind ram pressure
    const F_wind = rho_val * V_val * V_val * 4.0 * Math.PI * r_val * r_val;
    
    // X-ray pressure
    const F_xray = L_X_val / (this.c * 4.0 * Math.PI * r_val * r_val);
    
    // Total force
    let F_total = {re: F_grav + F_wind + F_xray, im: 0};
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
  expandStellarScale(massFactor, luminosityFactor) {
    const M_val = this.variables.get("M").re;
    const L_X_val = this.variables.get("L_X").re;
    
    this.variables.get("M").re = M_val * massFactor;
    this.variables.get("L_X").re = L_X_val * luminosityFactor;
    
    if (this.enableLogging) {
      console.log(`Expanded stellar scale: M×${massFactor}, L×${luminosityFactor}`);
    }
  },
  
  expandWindScale(velocityFactor, densityFactor) {
    const V_val = this.variables.get("V").re;
    const rho_val = this.variables.get("rho_gas").re;
    
    this.variables.get("V").re = V_val * velocityFactor;
    this.variables.get("rho_gas").re = rho_val * densityFactor;
    
    if (this.enableLogging) {
      console.log(`Expanded wind scale: V×${velocityFactor}, ρ×${densityFactor}`);
    }
  },
  
  expandSymbioticScale(orbitFactor, activityFactor) {
    const r_val = this.variables.get("r").re;
    const omega = this.variables.get("omega0").re;
    const k_act_val = this.variables.get("k_act").re;
    
    this.variables.get("r").re = r_val * orbitFactor;
    this.variables.get("omega0").re = omega / Math.sqrt(orbitFactor * orbitFactor * orbitFactor);
    this.variables.get("k_act").re = k_act_val * activityFactor;
    
    if (this.enableLogging) {
      console.log(`Expanded symbiotic scale: orbit×${orbitFactor}, activity×${activityFactor}`);
    }
  },
  
  optimizeForMetric(metricName) {
    const presets = {
      "standard": {
        DPM_momentum: 0.20,
        k_LENR: 1e-11,
        V: 1e5
      },
      "outburst": {
        DPM_momentum: 0.50,
        k_LENR: 5e-11,
        V: 2e5,
        L_X: 5e30,
        rho_gas: 5e-18
      },
      "quiescent": {
        DPM_momentum: 0.05,
        k_LENR: 2e-12,
        V: 5e4,
        L_X: 1e29
      },
      "jet_active": {
        DPM_momentum: 0.75,
        k_LENR: 1e-10,
        V: 3e5,
        k_act: 1e-17
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

addEnhancedDynamics(RAquariiUQFFModule, "R_Aquarii", domainExpansion);

module.exports = RAquariiUQFFModule;
