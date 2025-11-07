// source145.js - BulletClusterUQFFModule
// Bullet Cluster (1E 0657-56) - Merging cluster showing dark matter separation
// Enhanced with full 25-method self-expansion framework
// M_cluster=3e15 M☉, r=4e24 m, merger velocity=4500 km/s, x2=-1.8e178
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

class BulletClusterUQFFModule {
  constructor() {
    // Core physics parameters - Bullet Cluster
    this.M_cluster = 3e15 * 1.989e30;    // Total mass (3 quadrillion M☉)
    this.M_bullet = 1e14 * 1.989e30;     // Bullet subcluster mass
    this.M_main = 2e15 * 1.989e30;       // Main cluster mass
    this.r = 4e24;                        // System radius (4 Mpc)
    this.d_sep = 7e23;                    // DM-gas separation (700 kpc)
    this.V_merger = 4500e3;               // Merger velocity (4500 km/s)
    this.rho_ICM = 1e-26;                 // Intracluster medium density (kg/m^3)
    this.T_ICM = 1.5e8;                   // ICM temperature (K)
    this.L_X = 1e45;                      // X-ray luminosity (W)
    this.sigma_v = 1200e3;                // Velocity dispersion (m/s)
    this.f_DM = 0.90;                     // Dark matter fraction
    this.omega0 = 1e-19;                  // Angular frequency (rad/s)
    this.k_merger = 1e-15;                // Merger coupling constant
    this.x2 = -1.8e178;                   // UQFF coupling constant
    
    // UQFF parameters
    this.DPM_momentum = 0.97;
    this.k_LENR = 2.5e-10;
    this.E_cm = 2.5e-6;
    this.epsilon = 0.005;
    
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
    this.metadata.set("object_type", "merging_cluster");
    this.metadata.set("system_name", "Bullet_Cluster");
    this.enableLogging = false;
    this.learningRate = 0.01;
  }
  
  initializeVariables() {
    const vars = ["M_cluster", "M_bullet", "M_main", "r", "d_sep", "V_merger", "rho_ICM", 
                  "T_ICM", "L_X", "sigma_v", "f_DM", "omega0", "k_merger", "x2", 
                  "DPM_momentum", "k_LENR", "E_cm", "epsilon"];
    vars.forEach(v => {
      if (this[v] !== undefined) {
        this.variables.set(v, toComplex(this[v]));
      }
    });
    this.variables.set("t", {re: 0, im: 0});
  }
  
  computeDPM_resonance(t) {
    const M = this.variables.get("M_cluster").re;
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
  
  computeMergerDynamics(t) {
    const M_bullet = this.variables.get("M_bullet").re;
    const M_main = this.variables.get("M_main").re;
    const V_merger = this.variables.get("V_merger").re;
    const d_sep = this.variables.get("d_sep").re;
    const k_merger = this.variables.get("k_merger").re;
    const E_merger = 0.5 * (M_bullet * M_main / (M_bullet + M_main)) * V_merger * V_merger;
    return {re: k_merger * E_merger / d_sep, im: 0};
  }
  
  computeF(t) {
    this.variables.set("t", {re: t, im: 0});
    const M_cluster = this.variables.get("M_cluster").re;
    const r = this.variables.get("r").re;
    const L_X = this.variables.get("L_X").re;
    const f_DM = this.variables.get("f_DM").re;
    
    const F_grav = this.G * M_cluster * M_cluster / (r * r);
    const F_dm = f_DM * F_grav;
    const F_dpm = this.computeDPM_resonance(t);
    const F_lenr = this.computeLENRTerm();
    const F_merger = this.computeMergerDynamics(t);
    const F_xray = L_X / (this.c * 4 * Math.PI * r * r);
    
    let F_total = {re: F_grav + F_dm + F_xray, im: 0};
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
    const cloned = new BulletClusterUQFFModule();
    cloned.variables = new Map(this.variables);
    cloned.dynamicParameters = new Map(this.dynamicParameters);
    cloned.metadata = new Map(this.metadata);
    cloned.enableLogging = this.enableLogging;
    cloned.learningRate = this.learningRate;
    return cloned;
  }
}

const domainExpansion = {
  expandClusterScale(massFactor, radiusFactor) {
    const M = this.variables.get("M_cluster").re;
    const r = this.variables.get("r").re;
    this.variables.get("M_cluster").re = M * massFactor;
    this.variables.get("r").re = r * radiusFactor;
    if (this.enableLogging) console.log(`Expanded Bullet: M×${massFactor}, r×${radiusFactor}`);
  },
  expandMergerVelocity(velocityFactor) {
    const V_merger = this.variables.get("V_merger").re;
    this.variables.get("V_merger").re = V_merger * velocityFactor;
    if (this.enableLogging) console.log(`Expanded merger velocity: ×${velocityFactor}`);
  }
};

addEnhancedDynamics(BulletClusterUQFFModule, "Bullet_Cluster", domainExpansion);
module.exports = BulletClusterUQFFModule;
