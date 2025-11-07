// source141.js - TriangulumGalaxyUQFFModule
// M33 (Triangulum Galaxy) - Small spiral in Local Group
// Enhanced with full 25-method self-expansion framework
// M=5e10 M☉, r=3e22 m, satellite of Andromeda, x2=-2.1e174
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

class TriangulumGalaxyUQFFModule {
  constructor() {
    // Core physics parameters - M33 (Triangulum Galaxy)
    this.M = 5e10 * 1.989e30;        // Total mass (50 billion M☉)
    this.M_M31 = 1.5e12 * 1.989e30;  // Andromeda mass (M31 host)
    this.r = 3e22;                    // Galaxy radius (30 kpc)
    this.d_M31 = 2e23;                // Distance to M31 (200 kpc)
    this.rho_gas = 5e-25;             // Gas density (kg/m^3)
    this.SFR = 0.45 * 1.989e30 / (365.25 * 24 * 3600); // Star formation rate (0.45 M☉/yr)
    this.V_rot = 100e3;               // Rotation velocity (m/s)
    this.L_Halpha = 1e40;             // H-alpha luminosity (W)
    this.B0 = 3e-10;                  // Magnetic field (T)
    this.n_e = 500;                   // Electron density (m^-3)
    this.omega0 = 1e-15;              // Angular frequency (rad/s)
    this.T_val = 8000;                // Temperature (K)
    this.k_LG = 1e-15;                // Local Group coupling
    this.x2 = -2.1e174;               // UQFF coupling constant
    
    // UQFF parameters
    this.DPM_momentum = 0.85;
    this.k_LENR = 4e-11;
    this.E_cm = 1.88e-6;
    this.epsilon = 0.018;
    
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
    this.metadata.set("object_type", "local_group_spiral");
    this.metadata.set("system_name", "M33_Triangulum");
    this.enableLogging = false;
    this.learningRate = 0.01;
  }
  
  initializeVariables() {
    const vars = ["M", "M_M31", "r", "d_M31", "rho_gas", "SFR", "V_rot", "L_Halpha", 
                  "B0", "n_e", "omega0", "T_val", "k_LG", "x2", 
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
  
  computeLocalGroupTidal(t) {
    const M = this.variables.get("M").re;
    const M_M31 = this.variables.get("M_M31").re;
    const d_M31 = this.variables.get("d_M31").re;
    const k_LG = this.variables.get("k_LG").re;
    const omega = this.variables.get("omega0").re;
    const F_tidal = this.G * M * M_M31 / (d_M31 * d_M31);
    return {re: k_LG * F_tidal * (1 + 0.05 * Math.cos(omega * t)), im: 0};
  }
  
  computeF(t) {
    this.variables.set("t", {re: t, im: 0});
    const M = this.variables.get("M").re;
    const r = this.variables.get("r").re;
    const SFR = this.variables.get("SFR").re;
    const L_Halpha = this.variables.get("L_Halpha").re;
    
    const F_grav = this.G * M * M / (r * r);
    const F_dpm = this.computeDPM_resonance(t);
    const F_lenr = this.computeLENRTerm();
    const F_lg = this.computeLocalGroupTidal(t);
    const F_sf = SFR * this.c * this.c / (4 * Math.PI * r * r);
    const F_halpha = L_Halpha / (this.c * 4 * Math.PI * r * r);
    
    let F_total = {re: F_grav + F_sf + F_halpha, im: 0};
    F_total = complexAdd(F_total, F_dpm);
    F_total = complexAdd(F_total, F_lenr);
    F_total = complexAdd(F_total, F_lg);
    
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
    const cloned = new TriangulumGalaxyUQFFModule();
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
    if (this.enableLogging) console.log(`Expanded M33: M×${massFactor}, r×${radiusFactor}`);
  },
  expandLocalGroupInteraction(distanceFactor) {
    const d_M31 = this.variables.get("d_M31").re;
    this.variables.get("d_M31").re = d_M31 * distanceFactor;
    if (this.enableLogging) console.log(`Expanded LG distance: ×${distanceFactor}`);
  }
};

addEnhancedDynamics(TriangulumGalaxyUQFFModule, "M33_Triangulum", domainExpansion);
module.exports = TriangulumGalaxyUQFFModule;
