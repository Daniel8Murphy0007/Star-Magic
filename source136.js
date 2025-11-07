// source136.js - NGC1316FornaxAUQFFModule
// NGC 1316 (Fornax A) - Giant elliptical galaxy with merger remnants
// Enhanced with full 25-method self-expansion framework
// M=1.5e11 M☉, r=8e22 m, merger age ~3 Gyr, x2=-2.8e174
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

class NGC1316FornaxAUQFFModule {
  constructor() {
    // Core physics parameters - NGC 1316 (Fornax A)
    this.M = 1.5e11 * 1.989e30;     // Total mass (150 billion M☉)
    this.M_shell = 2e10 * 1.989e30;  // Shell mass from merger
    this.r = 8e22;                   // Galaxy radius (80 kpc)
    this.r_shell = 2e22;             // Shell radius (20 kpc)
    this.rho_gas = 5e-25;            // Gas density (kg/m^3)
    this.T_merger = 3e9 * 365.25 * 24 * 3600; // 3 Gyr in seconds
    this.V_infall = 500e3;           // Infall velocity (m/s)
    this.L_radio = 1e39;             // Radio luminosity (W)
    this.B0 = 5e-10;                 // Magnetic field (T)
    this.n_e = 1e3;                  // Electron density (m^-3)
    this.omega0 = 5e-17;             // Angular frequency (rad/s)
    this.T_val = 5e6;                // Temperature (K)
    this.k_merger = 1e-14;           // Merger coupling constant
    this.x2 = -2.8e174;              // UQFF coupling constant
    
    // UQFF parameters
    this.DPM_momentum = 0.88;
    this.k_LENR = 5e-11;
    this.E_cm = 1.95e-6;
    this.epsilon = 0.015;
    
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
    this.metadata.set("object_type", "merger_elliptical_galaxy");
    this.metadata.set("system_name", "NGC1316_Fornax_A");
    this.enableLogging = false;
    this.learningRate = 0.01;
  }
  
  initializeVariables() {
    const vars = ["M", "M_shell", "r", "r_shell", "rho_gas", "T_merger", "V_infall", 
                  "L_radio", "B0", "n_e", "omega0", "T_val", "k_merger", "x2", 
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
  
  computeMergerShells(t) {
    const M_shell = this.variables.get("M_shell").re;
    const r_shell = this.variables.get("r_shell").re;
    const V_infall = this.variables.get("V_infall").re;
    const k_merger = this.variables.get("k_merger").re;
    const T_merger = this.variables.get("T_merger").re;
    const age_fraction = t / T_merger;
    return {re: k_merger * M_shell * V_infall * V_infall / (r_shell * r_shell) * (1 - age_fraction), im: 0};
  }
  
  computeF(t) {
    this.variables.set("t", {re: t, im: 0});
    const M = this.variables.get("M").re;
    const r = this.variables.get("r").re;
    const L_radio = this.variables.get("L_radio").re;
    
    const F_grav = this.G * M * M / (r * r);
    const F_dpm = this.computeDPM_resonance(t);
    const F_lenr = this.computeLENRTerm();
    const F_merger = this.computeMergerShells(t);
    const F_radio = L_radio / (this.c * 4 * Math.PI * r * r);
    
    let F_total = {re: F_grav + F_radio, im: 0};
    F_total = complexAdd(F_total, F_dpm);
    F_total = complexAdd(F_total, F_lenr);
    F_total = complexAdd(F_total, F_merger);
    
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
    const cloned = new NGC1316FornaxAUQFFModule();
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
    if (this.enableLogging) console.log(`Expanded NGC1316: M×${massFactor}, r×${radiusFactor}`);
  },
  expandMergerDynamics(shellMassFactor) {
    const M_shell = this.variables.get("M_shell").re;
    this.variables.get("M_shell").re = M_shell * shellMassFactor;
    if (this.enableLogging) console.log(`Expanded merger shells: M×${shellMassFactor}`);
  }
};

addEnhancedDynamics(NGC1316FornaxAUQFFModule, "NGC1316_Fornax_A", domainExpansion);
module.exports = NGC1316FornaxAUQFFModule;
