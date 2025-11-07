// source152.js - VelaPulsarUQFFModule
// Vela Pulsar (PSR B0833-45) - supernova remnant neutron star
// Enhanced with full 25-method self-expansion framework
// M=2.8e30 kg (~1.4 M☉), r=1.7e17 m, omega0=1e-12, x2=-3.40e172
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

class VelaPulsarUQFFModule {
  constructor() {
    // Core physics parameters from C++ source152.cpp
    this.M = 2.8e30;             // Neutron star mass (kg) - ~1.4 solar masses
    this.r = 1.7e17;             // Pulsar radius scale (m) - ~10 km typical NS
    this.L_X = 1e27;             // X-ray luminosity (W)
    this.B0 = 3e-8;              // Surface magnetic field (T) - 30 μT
    this.omega0 = 1e-12;         // Angular frequency (rad/s) - pulsar rotation
    this.rho_gas = 1e-23;        // Surrounding gas density (kg/m^3)
    this.V = 1e-3;               // Particle velocity (m/s)
    this.T_val = 1e6;            // Temperature (K) - X-ray emitting gas
    this.n_e = 1e10;             // Electron density (m^-3)
    this.k_DE = 1e-30;           // Dark energy coupling
    this.k_act = 1e-6;           // Active coupling
    this.x2 = -3.40e172;         // UQFF coupling constant (pulsars)
    
    // UQFF parameters
    this.DPM_momentum = 0.93;
    this.k_LENR = 1e-10;
    this.E_cm = 3.0264e-8;       // 189 GeV
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
    this.metadata.set("object_type", "pulsar");
    this.metadata.set("system_name", "Vela_Pulsar");
    
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
    const B0_val = this.variables.get("B0").re;
    const omega_val = this.variables.get("omega0").re;
    const L_X_val = this.variables.get("L_X").re;
    
    // Gravitational force
    const F_grav = this.G * M_val * M_val / (r_val * r_val);
    
    // DPM resonance
    const F_dpm = this.computeDPM_resonance(t);
    
    // LENR contribution
    const F_lenr = this.computeLENRTerm(t);
    
    // Magnetic dipole radiation
    const F_mag = (this.mu0 / (6.0 * Math.PI * this.c)) * 
                   B0_val * B0_val * r_val ** 6 * omega_val ** 4;
    
    // X-ray luminosity pressure
    const F_xray = L_X_val / (this.c * 4.0 * Math.PI * r_val * r_val);
    
    // Total force
    let F_total = {re: F_grav + F_mag + F_xray, im: 0};
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
  
  // Clone for parallel processing
  clone() {
    const cloned = new VelaPulsarUQFFModule();
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
  expandPulsarScale(massFactor, radiusFactor) {
    const M_val = this.variables.get("M").re;
    const r_val = this.variables.get("r").re;
    const rho_val = this.variables.get("rho_gas").re;
    
    this.variables.get("M").re = M_val * massFactor;
    this.variables.get("r").re = r_val * radiusFactor;
    this.variables.get("rho_gas").re = rho_val / (radiusFactor ** 3);
    
    if (this.enableLogging) {
      console.log(`Expanded pulsar: M×${massFactor}, r×${radiusFactor}`);
    }
  },
  
  expandMagneticScale(fieldFactor, rotationFactor) {
    const B0_val = this.variables.get("B0").re;
    const omega_val = this.variables.get("omega0").re;
    const L_X_val = this.variables.get("L_X").re;
    
    this.variables.get("B0").re = B0_val * fieldFactor;
    this.variables.get("omega0").re = omega_val * rotationFactor;
    
    // Magnetic dipole radiation scales as B^2 * omega^4
    this.variables.get("L_X").re = L_X_val * fieldFactor * fieldFactor * 
                                     rotationFactor ** 4;
    
    if (this.enableLogging) {
      console.log(`Expanded magnetic: B×${fieldFactor}, ω×${rotationFactor}`);
    }
  },
  
  expandWindScale(velocityFactor, densityFactor) {
    const V_val = this.variables.get("V").re;
    const rho_val = this.variables.get("rho_gas").re;
    
    this.variables.get("V").re = V_val * velocityFactor;
    this.variables.get("rho_gas").re = rho_val * densityFactor;
    
    if (this.enableLogging) {
      console.log(`Expanded wind: V×${velocityFactor}, ρ×${densityFactor}`);
    }
  },
  
  optimizeForMetric(metricName) {
    const presets = {
      "standard": {
        DPM_momentum: 0.93,
        k_LENR: 1e-10,
        B0: 3e-8
      },
      "young_pulsar": {
        DPM_momentum: 1.2,
        k_LENR: 5e-10,
        B0: 1e-7,
        omega0: 1e-11,
        L_X: 1e28
      },
      "millisecond": {
        DPM_momentum: 0.5,
        k_LENR: 1e-11,
        B0: 1e-9,
        omega0: 1e-10,
        L_X: 1e26
      },
      "magnetar": {
        DPM_momentum: 1.5,
        k_LENR: 1e-9,
        B0: 1e-6,
        L_X: 1e30
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

addEnhancedDynamics(VelaPulsarUQFFModule, "Vela_Pulsar", domainExpansion);

module.exports = VelaPulsarUQFFModule;
