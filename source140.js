// source140.js - PinwheelGalaxyUQFFModule
// M101 (Pinwheel Galaxy) - Grand design spiral with asymmetric structure
// Enhanced with full 25-method self-expansion framework
// M=1e12 M☉, r=8.5e22 m, SFR=3 M☉/yr, x2=-5.1e175
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
function complexScale(a, s) { return {re: a.re*s, im: a.im*s}; }
function toComplex(x) { return typeof x === 'object' ? x : {re: x, im: 0}; }

class PinwheelGalaxyUQFFModule {
  constructor() {
    // Core physics parameters - M101 (Pinwheel Galaxy)
    this.M = 1e12 * 1.989e30;        // Total mass (1 trillion M☉)
    this.r = 8.5e22;                  // Galaxy radius (85 kpc)
    this.r_HII = 5e21;                // HII region scale (5 kpc)
    this.rho_gas = 1e-24;             // Gas density (kg/m^3)
    this.SFR = 3 * 1.989e30 / (365.25 * 24 * 3600); // Star formation rate (3 M☉/yr)
    this.V_rot = 250e3;               // Rotation velocity (m/s)
    this.n_HII = 3000;                // Number of HII regions
    this.L_UV = 1e43;                 // UV luminosity (W)
    this.B0 = 5e-10;                  // Magnetic field (T)
    this.n_e = 1e3;                   // Electron density (m^-3)
    this.omega0 = 3e-16;              // Angular frequency (rad/s)
    this.T_val = 1e4;                 // Temperature (K)
    this.k_asym = 1e-14;              // Asymmetry coupling
    this.x2 = -5.1e175;               // UQFF coupling constant
    
    // UQFF parameters
    this.DPM_momentum = 0.87;
    this.k_LENR = 9e-11;
    this.E_cm = 2.02e-6;
    this.epsilon = 0.016;
    
    // Constants
    this.G = 6.67430e-11;
    this.c = 299792458.0;
    this.hbar = 1.054571817e-34;
    this.m_e = 9.10938356e-31;
    this.k_B = 1.380649e-23;
    
    // Enhanced dynamics infrastructure
    this.variables = new Map();
    this.initializeVariables();
    this.dynamicTerms = [];
    this.dynamicParameters = new Map();
    this.metadata = new Map();
    this.metadata.set("enhanced", true);
    this.metadata.set("version", "2.0.0");
    this.metadata.set("object_type", "grand_design_spiral");
    this.metadata.set("system_name", "M101_Pinwheel");
    this.enableLogging = false;
    this.learningRate = 0.01;
  }
  
  initializeVariables() {
    const vars = ["M", "r", "r_HII", "rho_gas", "SFR", "V_rot", "n_HII", "L_UV", 
                  "B0", "n_e", "omega0", "T_val", "k_asym", "x2", 
                  "DPM_momentum", "k_LENR", "E_cm", "epsilon"];
    vars.forEach(v => {
      if (this[v] !== undefined) {
        this.variables.set(v, toComplex(this[v]));
      }
    });
    this.variables.set("t", {re: 0, im: 0});
  }
  
  computeDPM_resonance(t) {
    const M = this.variables.get("M").re;
    const r = this.variables.get("r").re;
    const omega = this.variables.get("omega0").re;
    const DPM = this.variables.get("DPM_momentum").re;
    const F_grav = this.G * M * M / (r * r);
    const x2 = this.variables.get("x2").re;
    return {re: DPM * F_grav * Math.cos(omega * t) * x2, im: 0};
  }
  
  computeLENRTerm() {
    const k = this.variables.get("k_LENR").re;
    const omega_L = 2 * Math.PI * 1.25e12;
    const omega0 = this.variables.get("omega0").re;
    return {re: k * Math.pow(omega_L / omega0, 2), im: 0};
  }
  
  computeHIIRegionPressure() {
    const n_HII = this.variables.get("n_HII").re;
    const r_HII = this.variables.get("r_HII").re;
    const L_UV = this.variables.get("L_UV").re;
    const T_val = this.variables.get("T_val").re;
    const P_HII = n_HII * this.k_B * T_val + L_UV / (this.c * 4 * Math.PI * r_HII * r_HII);
    return {re: P_HII, im: 0};
  }
  
  computeF(t) {
    this.variables.set("t", {re: t, im: 0});
    const M = this.variables.get("M").re;
    const r = this.variables.get("r").re;
    const SFR = this.variables.get("SFR").re;
    const V_rot = this.variables.get("V_rot").re;
    
    const F_grav = this.G * M * M / (r * r);
    const F_dpm = this.computeDPM_resonance(t);
    const F_lenr = this.computeLENRTerm();
    const F_hii = this.computeHIIRegionPressure();
    const F_rot = M * V_rot * V_rot / r; // Centrifugal force
    const F_sf = SFR * this.c * this.c / (4 * Math.PI * r * r);
    
    let F_total = {re: F_grav + F_rot + F_sf, im: 0};
    F_total = complexAdd(F_total, F_dpm);
    F_total = complexAdd(F_total, F_lenr);
    F_total = complexAdd(F_total, F_hii);
    
    this.dynamicTerms.forEach(term => {
      F_total = complexAdd(F_total, term.compute(this.variables));
    });
    
    return F_total;
  }
  
  setEnableLogging(enable) { this.enableLogging = enable; }
  registerDynamicTerm(term) { this.dynamicTerms.push(term); }
  setDynamicParameter(name, value) { this.dynamicParameters.set(name, value); }
  getDynamicParameter(name) { return this.dynamicParameters.get(name); }
  
  clone() {
    const cloned = new PinwheelGalaxyUQFFModule();
    cloned.variables = new Map(this.variables);
    cloned.dynamicParameters = new Map(this.dynamicParameters);
    cloned.metadata = new Map(this.metadata);
    cloned.enableLogging = this.enableLogging;
    cloned.learningRate = this.learningRate;
    return cloned;
  }
}

const domainExpansion = {
  expandGalaxyScale(massFactor, radiusFactor) {
    const M = this.variables.get("M").re;
    const r = this.variables.get("r").re;
    this.variables.get("M").re = M * massFactor;
    this.variables.get("r").re = r * radiusFactor;
    if (this.enableLogging) console.log(`Expanded M101: M×${massFactor}, r×${radiusFactor}`);
  },
  expandStarFormation(SFRfactor) {
    const SFR = this.variables.get("SFR").re;
    const n_HII = this.variables.get("n_HII").re;
    this.variables.get("SFR").re = SFR * SFRfactor;
    this.variables.get("n_HII").re = n_HII * SFRfactor;
    if (this.enableLogging) console.log(`Expanded star formation: ×${SFRfactor}`);
  }
};

addEnhancedDynamics(PinwheelGalaxyUQFFModule, "M101_Pinwheel", domainExpansion);
module.exports = PinwheelGalaxyUQFFModule;
