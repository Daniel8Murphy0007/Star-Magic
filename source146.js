// source146.js - PerseusClusterUQFFModule
// Perseus Galaxy Cluster - Massive cluster with central AGN
// Enhanced with full 25-method self-expansion framework
// M_cluster=1.3e15 M☉, r=2e24 m, central AGN outbursts, x2=-1.1e178
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

class PerseusClusterUQFFModule {
  constructor() {
    // Core physics parameters - Perseus Cluster
    this.M_cluster = 1.3e15 * 1.989e30;  // Total mass (1.3 quadrillion M☉)
    this.M_NGC1275 = 1e12 * 1.989e30;    // NGC 1275 (Perseus A) mass
    this.M_BH = 3.4e8 * 1.989e30;        // Central SMBH mass
    this.r = 2e24;                        // Cluster radius (2 Mpc)
    this.r_core = 1.5e23;                 // Cool core radius (150 kpc)
    this.rho_ICM = 3e-26;                 // Intracluster medium density (kg/m^3)
    this.T_ICM = 6e7;                     // ICM temperature (K)
    this.L_X = 1e45;                      // X-ray luminosity (W)
    this.P_AGN = 1e43;                    // AGN outburst power (W)
    this.t_outburst = 1e7 * 365.25 * 24 * 3600; // Outburst timescale (10 Myr)
    this.sigma_v = 1350e3;                // Velocity dispersion (m/s)
    this.B0 = 2e-9;                       // Magnetic field (T)
    this.omega0 = 2e-19;                  // Angular frequency (rad/s)
    this.k_AGN = 1e-14;                   // AGN coupling constant
    this.x2 = -1.1e178;                   // UQFF coupling constant
    
    // UQFF parameters
    this.DPM_momentum = 0.98;
    this.k_LENR = 1.8e-10;
    this.E_cm = 2.45e-6;
    this.epsilon = 0.006;
    
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
    this.metadata.set("object_type", "cool_core_cluster");
    this.metadata.set("system_name", "Perseus_Cluster");
    this.enableLogging = false;
    this.learningRate = 0.01;
  }
  
  initializeVariables() {
    const vars = ["M_cluster", "M_NGC1275", "M_BH", "r", "r_core", "rho_ICM", "T_ICM", 
                  "L_X", "P_AGN", "t_outburst", "sigma_v", "B0", "omega0", "k_AGN", "x2", 
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
  
  computeAGNFeedback(t) {
    const P_AGN = this.variables.get("P_AGN").re;
    const t_outburst = this.variables.get("t_outburst").re;
    const r_core = this.variables.get("r_core").re;
    const k_AGN = this.variables.get("k_AGN").re;
    const omega = 2 * Math.PI / t_outburst;
    const duty_cycle = 0.5 * (1 + Math.cos(omega * t)); // Periodic outbursts
    return {re: k_AGN * P_AGN * duty_cycle / (4 * Math.PI * r_core * r_core), im: 0};
  }
  
  computeF(t) {
    this.variables.set("t", {re: t, im: 0});
    const M_cluster = this.variables.get("M_cluster").re;
    const r = this.variables.get("r").re;
    const L_X = this.variables.get("L_X").re;
    
    const F_grav = this.G * M_cluster * M_cluster / (r * r);
    const F_dpm = this.computeDPM_resonance(t);
    const F_lenr = this.computeLENRTerm();
    const F_agn = this.computeAGNFeedback(t);
    const F_xray = L_X / (this.c * 4 * Math.PI * r * r);
    
    let F_total = {re: F_grav + F_xray, im: 0};
    F_total = complexAdd(F_total, F_dpm);
    F_total = complexAdd(F_total, F_lenr);
    F_total = complexAdd(F_total, F_agn);
    
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
    const cloned = new PerseusClusterUQFFModule();
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
    if (this.enableLogging) console.log(`Expanded Perseus: M×${massFactor}, r×${radiusFactor}`);
  },
  expandAGNActivity(powerFactor) {
    const P_AGN = this.variables.get("P_AGN").re;
    this.variables.get("P_AGN").re = P_AGN * powerFactor;
    if (this.enableLogging) console.log(`Expanded AGN power: ×${powerFactor}`);
  }
};

addEnhancedDynamics(PerseusClusterUQFFModule, "Perseus_Cluster", domainExpansion);
module.exports = PerseusClusterUQFFModule;
