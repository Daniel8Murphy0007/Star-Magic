// source155.js - UQFFBuoyancyModule
// Multi-system UQFF buoyancy module supporting multiple astronomical objects
// Enhanced with full 25-method self-expansion framework
// Systems: CrabNebula, ElGordo, and custom configurations
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

class UQFFBuoyancyModule {
  constructor(systemName = "CrabNebula") {
    // System configurations
    this.systemConfigs = {
      "CrabNebula": {
        M: 1e31,
        r: 4.63e16,
        L_X: 1e31,
        B0: 1e-8,
        rho_gas: 1e-21,
        t: 9.504e6,
        omega0: 1e-12,
        x2: -3.40e172
      },
      "ElGordo": {
        M: 4.97e45,
        r: 3.09e22,
        L_X: 2e38,
        B0: 1e-10,
        rho_gas: 1e-24,
        t: 2.21e16,
        omega0: 1e-15,
        x2: -2.27e172
      }
    };
    
    // Set current system
    this.currentSystem = systemName;
    const config = this.systemConfigs[systemName] || this.systemConfigs["CrabNebula"];
    
    // Apply configuration
    this.M = config.M;
    this.r = config.r;
    this.L_X = config.L_X;
    this.B0 = config.B0;
    this.rho_gas = config.rho_gas;
    this.omega0 = config.omega0;
    this.x2 = config.x2;
    
    // Common parameters
    this.T_val = 1e7;
    this.V = 1e5;
    this.n_e = 1e6;
    this.k_DE = 1e-16;
    this.k_act = 1e-14;
    
    // UQFF parameters
    this.DPM_momentum = 0.93;
    this.DPM_gravity = 1.0;
    this.DPM_stability = 0.01;
    this.k_LENR = 1e-10;
    this.k_neutron = 1e10;
    this.sigma_n = 1e-4;
    this.k_rel = 1e-10;
    this.E_cm = 3.0264e-8;
    this.E_cm_astro = 1.24e24;
    
    // Buoyancy parameters
    this.beta_i = 0.6;
    this.V_infl_UA = 1e-6;
    this.rho_vac_UA = 7.09e-36;
    this.rho_vac_A = 1e-30;
    this.a_universal = 1e12;
    
    // Superconductive parameters
    this.lambda_i = 1.0;
    this.rho_vac_SCm = 7.09e-37;
    this.omega_s = 2.5e-6;
    this.f_TRZ = 0.1;
    this.t_scale = 1e16;
    
    // Constants
    this.G = 6.67430e-11;
    this.c = 299792458.0;
    this.hbar = 1.054571817e-34;
    this.m_n = 1.674927498e-27;
    this.m_e = 9.1093837015e-31;
    this.k_B = 1.380649e-23;
    this.mu0 = 1.25663706212e-6;
    this.q = 1.602176634e-19;
    
    // Initialize variable storage
    this.variables = new Map();
    this.initializeVariables();
    
    this.dynamicTerms = [];
    this.dynamicParameters = new Map();
    
    this.metadata = new Map();
    this.metadata.set("enhanced", true);
    this.metadata.set("version", "2.0.0");
    this.metadata.set("object_type", "multi_system_buoyancy");
    this.metadata.set("system_name", systemName);
    
    this.enableLogging = false;
  }
  
  initializeVariables() {
    this.variables.set("M", {re: this.M, im: 0});
    this.variables.set("r", {re: this.r, im: 0});
    this.variables.set("L_X", {re: this.L_X, im: 0});
    this.variables.set("B0", {re: this.B0, im: 0});
    this.variables.set("omega0", {re: this.omega0, im: 0});
    this.variables.set("rho_gas", {re: this.rho_gas, im: 0});
    this.variables.set("T_val", {re: this.T_val, im: 0});
    this.variables.set("V", {re: this.V, im: 0});
    this.variables.set("n_e", {re: this.n_e, im: 0});
    this.variables.set("k_DE", {re: this.k_DE, im: 0});
    this.variables.set("k_act", {re: this.k_act, im: 0});
    this.variables.set("DPM_momentum", {re: this.DPM_momentum, im: 0});
    this.variables.set("DPM_gravity", {re: this.DPM_gravity, im: 0});
    this.variables.set("DPM_stability", {re: this.DPM_stability, im: 0});
    this.variables.set("k_LENR", {re: this.k_LENR, im: 0});
    this.variables.set("x2", {re: this.x2, im: 0});
    this.variables.set("beta_i", {re: this.beta_i, im: 0});
    this.variables.set("t", {re: 0, im: 0});
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
  
  setSystem(systemName) {
    if (this.systemConfigs[systemName]) {
      this.currentSystem = systemName;
      const config = this.systemConfigs[systemName];
      
      Object.entries(config).forEach(([key, value]) => {
        if (this.variables.has(key)) {
          this.variables.get(key).re = value;
        }
        this[key] = value;
      });
      
      this.metadata.set("system_name", systemName);
      
      if (this.enableLogging) {
        console.log(`Switched to system: ${systemName}`);
      }
    }
  }
  
  addSystemConfig(name, config) {
    this.systemConfigs[name] = config;
    if (this.enableLogging) {
      console.log(`Added system configuration: ${name}`);
    }
  }
  
  computeDPM_resonance(t) {
    const M_val = this.variables.get("M").re;
    const r_val = this.variables.get("r").re;
    const omega = this.variables.get("omega0").re;
    const DPM_mom = this.variables.get("DPM_momentum").re;
    
    const F_grav = this.G * M_val * M_val / (r_val * r_val);
    const osc = Math.cos(omega * t);
    const x2_val = this.variables.get("x2").re;
    return {re: DPM_mom * F_grav * osc * x2_val, im: 0};
  }
  
  computeLENRTerm(t) {
    const k_LENR = this.variables.get("k_LENR").re;
    const omega_LENR = 2 * Math.PI * 1.25e12;
    const omega0 = this.variables.get("omega0").re;
    
    return {re: k_LENR * Math.pow(omega_LENR / omega0, 2), im: 0};
  }
  
  computeIntegrand(t) {
    const M_val = this.variables.get("M").re;
    const r_val = this.variables.get("r").re;
    const L_X_val = this.variables.get("L_X").re;
    const B0_val = this.variables.get("B0").re;
    const V_val = this.variables.get("V").re;
    const DPM_mom = this.variables.get("DPM_momentum").re;
    const DPM_grav = this.variables.get("DPM_gravity").re;
    const DPM_stab = this.variables.get("DPM_stability").re;
    const k_DE = this.variables.get("k_DE").re;
    const k_act = this.variables.get("k_act").re;
    
    const theta = Math.PI / 4;
    const phi = Math.PI / 4;
    const omega_act = 2 * Math.PI * 300;
    
    // All UQFF terms
    const F0 = 1.83e71;
    const term_base = -F0;
    const term_mom = (this.m_e * this.c * this.c / (r_val * r_val)) * DPM_mom * Math.cos(theta);
    const term_grav = (this.G * M_val / (r_val * r_val)) * DPM_grav;
    const term_vac = this.rho_vac_UA * DPM_stab;
    const term_LENR = this.computeLENRTerm(t).re;
    const term_act = k_act * Math.cos(omega_act * t + phi);
    const term_DE = k_DE * L_X_val;
    const term_res = 2 * this.q * B0_val * V_val * Math.sin(theta) * this.computeDPM_resonance(t).re;
    const term_neut = this.k_neutron * this.sigma_n;
    const term_rel = this.k_rel * Math.pow(this.E_cm_astro / this.E_cm, 2);
    
    const total = term_base + term_mom + term_grav + term_vac + term_LENR + 
                  term_act + term_DE + term_res + term_neut + term_rel;
    
    return {re: total, im: 0};
  }
  
  computeFBi(t) {
    const integ = this.computeIntegrand(t);
    const x2_val = this.variables.get("x2").re;
    return complexScale(integ, x2_val);
  }
  
  computeF(t) {
    this.variables.set("t", {re: t, im: 0});
    
    let F_total = this.computeFBi(t);
    
    // Add dynamic terms
    this.dynamicTerms.forEach(term => {
      F_total = complexAdd(F_total, term.compute(this.variables));
    });
    
    return F_total;
  }
  
  computeBuoyancy(t) {
    const rho = this.variables.get("rho_gas").re;
    const T = this.variables.get("T_val").re;
    const r_val = this.variables.get("r").re;
    const beta_i = this.variables.get("beta_i").re;
    
    const pressure = rho * this.k_B * T / this.m_n;
    const volume = (4.0/3.0) * Math.PI * r_val * r_val * r_val;
    const buoyancy = beta_i * pressure * volume;
    
    return {re: buoyancy, im: 0};
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
  expandSystemScale(massFactor, radiusFactor) {
    const M_val = this.variables.get("M").re;
    const r_val = this.variables.get("r").re;
    const rho_val = this.variables.get("rho_gas").re;
    
    this.variables.get("M").re = M_val * massFactor;
    this.variables.get("r").re = r_val * radiusFactor;
    this.variables.get("rho_gas").re = rho_val * massFactor / (radiusFactor ** 3);
    
    if (this.enableLogging) {
      console.log(`Expanded system: M×${massFactor}, r×${radiusFactor}`);
    }
  },
  
  expandBuoyancyScale(buoyancyFactor, densityFactor) {
    const beta_i = this.variables.get("beta_i").re;
    const rho_val = this.variables.get("rho_gas").re;
    
    this.variables.get("beta_i").re = beta_i * buoyancyFactor;
    this.variables.get("rho_gas").re = rho_val * densityFactor;
    
    if (this.enableLogging) {
      console.log(`Expanded buoyancy: β×${buoyancyFactor}, ρ×${densityFactor}`);
    }
  },
  
  expandMultiSystemScale(systemScaleFactor) {
    // Scale all systems uniformly
    Object.keys(this.systemConfigs).forEach(sysName => {
      const config = this.systemConfigs[sysName];
      config.M *= systemScaleFactor;
      config.r *= Math.sqrt(systemScaleFactor);
      config.rho_gas /= Math.sqrt(systemScaleFactor);
    });
    
    if (this.enableLogging) {
      console.log(`Scaled all systems by factor: ${systemScaleFactor}`);
    }
  },
  
  optimizeForMetric(metricName) {
    const presets = {
      "standard": {
        DPM_momentum: 0.93,
        DPM_gravity: 1.0,
        DPM_stability: 0.01,
        k_LENR: 1e-10
      },
      "high_energy": {
        DPM_momentum: 1.5,
        DPM_gravity: 1.2,
        DPM_stability: 0.05,
        k_LENR: 5e-10
      },
      "quiescent": {
        DPM_momentum: 0.5,
        DPM_gravity: 0.8,
        DPM_stability: 0.005,
        k_LENR: 1e-11
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

addEnhancedDynamics(UQFFBuoyancyModule, "Multi_System_Buoyancy", domainExpansion);

module.exports = UQFFBuoyancyModule;
