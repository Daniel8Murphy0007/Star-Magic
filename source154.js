// source154.js - HydrogenResonanceUQFFModule
// Hydrogen Resonance for Periodic Table of Elements (PToE) with Surface Magnetic Field
// Enhanced with full 25-method self-expansion framework
// Multi-element support (Z=1-118), surface magnetic field coupling
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

class HydrogenResonanceUQFFModule {
  constructor() {
    // Surface Magnetic Field parameters (Sun baseline)
    this.B_s_min = 1e-4;         // T (quiet Sun)
    this.B_s_max = 0.4;          // T (sunspot)
    this.B_ref = 0.4;            // T (reference max)
    this.k_3 = 1.8;              // Coupling constant
    this.omega_s = 2.5e-6;       // rad/s (solar cycle)
    this.P_core = 1.0;           // Core pressure factor
    this.E_react = 1e46;         // J (reaction energy)
    
    // PToE parameters (default Z=1, A=1 - Hydrogen)
    this.Z = 1;                  // Atomic number
    this.A = 1;                  // Mass number
    this.E_bind = 0;             // Binding energy (MeV) - H has none
    this.f_res = 1.0;            // Resonance factor
    this.f_dp = 0.1;             // Deep pairing factor
    this.phi_dp = Math.PI/4;     // Pairing phase
    
    // UQFF parameters
    this.DPM_momentum = 0.5;
    this.k_LENR = 1e-12;
    this.x2 = -1.35e172;
    
    // Constants
    this.G = 6.67430e-11;
    this.c = 299792458.0;
    this.hbar = 1.054571817e-34;
    this.m_n = 1.674927498e-27;
    this.m_e = 9.1093837015e-31;
    this.k_B = 1.380649e-23;
    this.mu0 = 1.25663706212e-6;
    this.m_p = 1.672621898e-27;
    
    // Initialize variable storage
    this.variables = new Map();
    this.variables.set("B_s_min", {re: this.B_s_min, im: 0});
    this.variables.set("B_s_max", {re: this.B_s_max, im: 0});
    this.variables.set("B_ref", {re: this.B_ref, im: 0});
    this.variables.set("k_3", {re: this.k_3, im: 0});
    this.variables.set("omega_s", {re: this.omega_s, im: 0});
    this.variables.set("P_core", {re: this.P_core, im: 0});
    this.variables.set("E_react", {re: this.E_react, im: 0});
    this.variables.set("Z", {re: this.Z, im: 0});
    this.variables.set("A", {re: this.A, im: 0});
    this.variables.set("E_bind", {re: this.E_bind, im: 0});
    this.variables.set("f_res", {re: this.f_res, im: 0});
    this.variables.set("f_dp", {re: this.f_dp, im: 0});
    this.variables.set("phi_dp", {re: this.phi_dp, im: 0});
    this.variables.set("DPM_momentum", {re: this.DPM_momentum, im: 0});
    this.variables.set("k_LENR", {re: this.k_LENR, im: 0});
    this.variables.set("t", {re: 0, im: 0});
    this.variables.set("B_s", {re: this.B_ref, im: 0});
    
    this.dynamicTerms = [];
    this.dynamicParameters = new Map();
    
    this.metadata = new Map();
    this.metadata.set("enhanced", true);
    this.metadata.set("version", "2.0.0");
    this.metadata.set("object_type", "hydrogen_resonance");
    this.metadata.set("system_name", "PToE_H_Resonance");
    
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
  
  setElement(Z, A, E_bind_MeV) {
    this.Z = Z;
    this.A = A;
    this.E_bind = E_bind_MeV * 1.60218e-13; // Convert MeV to J
    this.variables.set("Z", {re: Z, im: 0});
    this.variables.set("A", {re: A, im: 0});
    this.variables.set("E_bind", {re: this.E_bind, im: 0});
    
    // Adjust resonance factors based on nuclear properties
    const N = A - Z;
    this.f_res = 1.0 + 0.1 * Math.log(A);
    
    // Pairing term (even-odd effect)
    if (Z % 2 === 0 && N % 2 === 0) {
      this.f_dp = 0.2; // Even-even (more stable)
    } else if (Z % 2 === 1 && N % 2 === 1) {
      this.f_dp = 0.05; // Odd-odd (less stable)
    } else {
      this.f_dp = 0.1; // Even-odd or odd-even
    }
    
    this.variables.set("f_res", {re: this.f_res, im: 0});
    this.variables.set("f_dp", {re: this.f_dp, im: 0});
  }
  
  computeB_j(t, B_s) {
    const B_ref = this.variables.get("B_ref").re;
    const omega_s = this.variables.get("omega_s").re;
    const base_b = B_ref + 0.4 * Math.sin(omega_s * t);
    return base_b * (B_s / B_ref);
  }
  
  computeU_g3(t, B_s) {
    const k_3 = this.variables.get("k_3").re;
    const omega_s = this.variables.get("omega_s").re;
    const P_core = this.variables.get("P_core").re;
    const E_react = this.variables.get("E_react").re;
    const B_j = this.computeB_j(t, B_s);
    
    return k_3 * B_j * Math.cos(omega_s * t * Math.PI) * P_core * E_react;
  }
  
  computeA_res(Z, A) {
    const alpha = 1.0 / 137.036;
    const base = Math.pow(Z, 2) * alpha / A;
    return {re: base, im: 0};
  }
  
  computeF_res(E_bind, A) {
    if (A === 0 || E_bind === 0) return {re: 1.0, im: 0};
    const binding_per_nucleon = E_bind / A;
    const factor = 1.0 + 0.01 * binding_per_nucleon / (1.60218e-13); // Normalize to MeV
    return {re: factor, im: 0};
  }
  
  computeU_dp(A1, A2) {
    const f_dp = this.variables.get("f_dp").re;
    const phi_dp = this.variables.get("phi_dp").re;
    const coupling = f_dp * Math.sqrt(A1 * A2) * Math.cos(phi_dp);
    return {re: coupling * 1e-13, im: 0}; // Scale to nuclear energies
  }
  
  computeHRes(Z, A, t) {
    this.variables.set("t", {re: t, im: 0});
    const E_bind = this.variables.get("E_bind").re;
    
    const A_res = this.computeA_res(Z, A);
    const F_res = this.computeF_res(E_bind, A);
    const U_dp = this.computeU_dp(A, A);
    
    // Combine resonance terms
    let H_res = complexMul(A_res, F_res);
    H_res = complexAdd(H_res, U_dp);
    
    // Add LENR contribution
    const k_LENR = this.variables.get("k_LENR").re;
    const LENR_term = {re: k_LENR * Z * A, im: 0};
    H_res = complexAdd(H_res, LENR_term);
    
    // Surface magnetic field contribution
    const B_s = this.variables.get("B_s").re;
    const U_g3 = this.computeU_g3(t, B_s);
    const mag_contrib = {re: U_g3 * 1e-46, im: 0}; // Scale down for combination
    H_res = complexAdd(H_res, mag_contrib);
    
    // Add dynamic terms
    this.dynamicTerms.forEach(term => {
      H_res = complexAdd(H_res, term.compute(this.variables));
    });
    
    return H_res;
  }
  
  computeF(t) {
    const Z = this.variables.get("Z").re;
    const A = this.variables.get("A").re;
    return this.computeHRes(Z, A, t);
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
  expandElementScale(massFactor, bindingFactor) {
    const A_val = this.variables.get("A").re;
    const E_bind_val = this.variables.get("E_bind").re;
    
    this.variables.get("A").re = Math.round(A_val * massFactor);
    this.variables.get("E_bind").re = E_bind_val * bindingFactor;
    
    if (this.enableLogging) {
      console.log(`Expanded element: A×${massFactor}, E_bind×${bindingFactor}`);
    }
  },
  
  expandMagneticScale(fieldFactor, cycleFactor) {
    const B_ref_val = this.variables.get("B_ref").re;
    const omega_s_val = this.variables.get("omega_s").re;
    
    this.variables.get("B_ref").re = B_ref_val * fieldFactor;
    this.variables.get("B_s_min").re *= fieldFactor;
    this.variables.get("B_s_max").re *= fieldFactor;
    this.variables.get("omega_s").re = omega_s_val * cycleFactor;
    
    if (this.enableLogging) {
      console.log(`Expanded magnetic: B×${fieldFactor}, ω×${cycleFactor}`);
    }
  },
  
  expandResonanceScale(resonanceFactor, pairingFactor) {
    const f_res_val = this.variables.get("f_res").re;
    const f_dp_val = this.variables.get("f_dp").re;
    
    this.variables.get("f_res").re = f_res_val * resonanceFactor;
    this.variables.get("f_dp").re = f_dp_val * pairingFactor;
    
    if (this.enableLogging) {
      console.log(`Expanded resonance: f_res×${resonanceFactor}, f_dp×${pairingFactor}`);
    }
  },
  
  optimizeForMetric(metricName) {
    const presets = {
      "hydrogen": {
        Z: 1,
        A: 1,
        E_bind: 0,
        f_res: 1.0,
        f_dp: 0.0
      },
      "helium": {
        Z: 2,
        A: 4,
        E_bind: 28.3 * 1.60218e-13,
        f_res: 1.1,
        f_dp: 0.2
      },
      "iron": {
        Z: 26,
        A: 56,
        E_bind: 492.3 * 1.60218e-13,
        f_res: 1.4,
        f_dp: 0.2
      },
      "uranium": {
        Z: 92,
        A: 238,
        E_bind: 1801.7 * 1.60218e-13,
        f_res: 1.8,
        f_dp: 0.2
      },
      "quiet_sun": {
        B_ref: 0.4,
        B_s: 1e-4,
        omega_s: 2.5e-6
      },
      "sunspot": {
        B_ref: 0.4,
        B_s: 0.4,
        omega_s: 2.5e-6
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

addEnhancedDynamics(HydrogenResonanceUQFFModule, "H_Resonance_PToE", domainExpansion);

module.exports = HydrogenResonanceUQFFModule;
