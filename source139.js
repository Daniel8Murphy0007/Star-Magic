// source139.js - WhirlpoolGalaxyUQFFModule
// M51 (Whirlpool Galaxy) - Grand design spiral interacting with NGC 5195
// Enhanced with full 25-method self-expansion framework
// M=1.6e11 M☉, r=4e22 m, tidal interaction, x2=-4.5e174
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

class WhirlpoolGalaxyUQFFModule {
  constructor() {
    // Core physics parameters - M51 (Whirlpool Galaxy)
    this.M = 1.6e11 * 1.989e30;      // M51 mass (160 billion M☉)
    this.M_companion = 8e10 * 1.989e30; // NGC 5195 mass (80 billion M☉)
    this.r = 4e22;                    // M51 radius (40 kpc)
    this.d_sep = 1e22;                // Separation distance (10 kpc)
    this.rho_gas = 2e-24;             // Gas density (kg/m^3)
    this.SFR = 3 * 1.989e30 / (365.25 * 24 * 3600); // Star formation rate (3 M☉/yr)
    this.V_rot = 200e3;               // Rotation velocity (m/s)
    this.B0 = 1e-9;                   // Magnetic field (T)
    this.n_e = 5e3;                   // Electron density (m^-3)
    this.omega0 = 5e-16;              // Angular frequency (rad/s)
    this.T_val = 1e4;                 // Temperature (K)
    this.k_tidal = 1e-13;             // Tidal coupling constant
    this.x2 = -4.5e174;               // UQFF coupling constant
    
    // UQFF parameters
    this.DPM_momentum = 0.89;
    this.k_LENR = 7e-11;
    this.E_cm = 1.98e-6;
    this.epsilon = 0.014;
    
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
    this.metadata.set("object_type", "interacting_spiral_galaxy");
    this.metadata.set("system_name", "M51_Whirlpool");
    this.enableLogging = false;
    this.learningRate = 0.01;
  }
  
  initializeVariables() {
    const vars = ["M", "M_companion", "r", "d_sep", "rho_gas", "SFR", "V_rot", 
                  "B0", "n_e", "omega0", "T_val", "k_tidal", "x2", 
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
  
  computeTidalForce(t) {
    const M = this.variables.get("M").re;
    const M_comp = this.variables.get("M_companion").re;
    const d_sep = this.variables.get("d_sep").re;
    const k_tidal = this.variables.get("k_tidal").re;
    const omega = this.variables.get("omega0").re;
    const F_tidal = 2 * this.G * M * M_comp / Math.pow(d_sep, 3) * this.r;
    return {re: k_tidal * F_tidal * (1 + 0.1 * Math.cos(omega * t)), im: 0};
  }
  
  computeF(t) {
    this.variables.set("t", {re: t, im: 0});
    const M = this.variables.get("M").re;
    const r = this.variables.get("r").re;
    const SFR = this.variables.get("SFR").re;
    
    const F_grav = this.G * M * M / (r * r);
    const F_dpm = this.computeDPM_resonance(t);
    const F_lenr = this.computeLENRTerm();
    const F_tidal = this.computeTidalForce(t);
    const F_sf = SFR * this.c * this.c / (4 * Math.PI * r * r); // Star formation energy
    
    let F_total = {re: F_grav + F_sf, im: 0};
    F_total = complexAdd(F_total, F_dpm);
    F_total = complexAdd(F_total, F_lenr);
    F_total = complexAdd(F_total, F_tidal);
    
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
    const cloned = new WhirlpoolGalaxyUQFFModule();
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
    if (this.enableLogging) console.log(`Expanded M51: M×${massFactor}, r×${radiusFactor}`);
  },
  expandTidalInteraction(separationFactor) {
    const d_sep = this.variables.get("d_sep").re;
    this.variables.get("d_sep").re = d_sep * separationFactor;
    if (this.enableLogging) console.log(`Expanded tidal separation: ×${separationFactor}`);
  }
};

addEnhancedDynamics(WhirlpoolGalaxyUQFFModule, "M51_Whirlpool", domainExpansion);
module.exports = WhirlpoolGalaxyUQFFModule;
