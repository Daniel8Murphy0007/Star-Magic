// source132.js - Butterfly Nebula (NGC 6302) ENHANCED + SIMULATION-READY
// UQFF Force with integral calculations, LENR/resonance/buoyancy dynamics
// Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

const { addEnhancedDynamics } = require('./enhanced_dynamics.js');

class ButterflyNebulaUQFFModule {
  constructor() {
    const M_sun = 1.989e30;
    this.variables = new Map([
      ["G", 6.6743e-11], ["c", 2.998e8], ["m_e", 9.109e-31], ["hbar", 1.0546e-34],
      ["mu_B", 9.274e-24], ["e", 1.602e-19], ["M_sun", M_sun], ["q", 1.602e-19], ["pi", Math.PI],
      ["M", 0.64 * M_sun], ["r", 3.22e19], ["x1", 0.0], ["x2", 3.22e19], ["level", 13.0],
      ["F0", 1.0], ["theta", 0.0], ["DPM_momentum", 1.0], ["DPM_gravity", 1.0], ["DPM_stability", 0.01],
      ["rho_vac_UA", 7.09e-36], ["k_LENR", 1.0], ["omega_LENR", 7.85e12], ["omega_0", 1e-12],
      ["k_act", 1.0], ["omega_act", 1.0], ["k_DE", 1.0], ["L_x", 1.0], ["B_0", 1.0], ["V", 1.0],
      ["g", 9.8], ["k_neutron", 1e10], ["sigma_n", 1e-4], ["k_rel", 1.0], ["E_cm", 1.0],
      ["E_cm_eff", 1.0], ["F_Sweet_vac", 7.09e-39], ["F_Kozima", 7.85e30], ["t", 0.0]
    ]);
    
    this.dynamicParameters = new Map(); this.dynamicTerms = [];
    this.metadata = new Map([["enhanced", true], ["version", "2.0.0"], ["simulation_ready", true]]);
    this.enableDynamicTerms = true; this.enableLogging = false; this.learningRate = 0.001;
  }
  
  updateVariable(n, v) { this.variables.set(n, v); }
  addToVariable(n, d) { this.variables.has(n) ? this.variables.set(n, this.variables.get(n) + d) : this.variables.set(n, d); }
  subtractFromVariable(n, d) { this.addToVariable(n, -d); }
  
  computeDPM_momentum_term(r) {
    const m_e_c2 = this.variables.get("m_e") * Math.pow(this.variables.get("c"), 2);
    return (m_e_c2 / (r * r)) * this.variables.get("DPM_momentum") * Math.cos(this.variables.get("theta"));
  }
  
  computeDPM_gravity_term(r) {
    return (this.variables.get("G") * this.variables.get("M") / (r * r)) * this.variables.get("DPM_gravity");
  }
  
  computeDPM_stability_term() { return this.variables.get("rho_vac_UA") * this.variables.get("DPM_stability"); }
  
  computeLENR_term(r, t) {
    return this.variables.get("k_LENR") * Math.cos(this.variables.get("omega_LENR") * t) * Math.exp(-r / this.variables.get("L_x"));
  }
  
  computeActivation_term(r, t) {
    return this.variables.get("k_act") * Math.cos(this.variables.get("omega_act") * t) * Math.exp(-r / (2 * this.variables.get("L_x")));
  }
  
  computeDE_term(r) { return this.variables.get("k_DE") * r; }
  
  computeEM_term(r, t) {
    const B = this.variables.get("B_0") * Math.exp(-r / this.variables.get("L_x"));
    return this.variables.get("mu_B") * B * Math.cos(this.variables.get("omega_0") * t);
  }
  
  computeNeutron_term(r) {
    return this.variables.get("k_neutron") * Math.exp(-this.variables.get("sigma_n") * r);
  }
  
  computeRelativistic_term(r) {
    return this.variables.get("k_rel") * (this.variables.get("E_cm_eff") / this.variables.get("E_cm")) * Math.exp(-r / this.variables.get("L_x"));
  }
  
  computeIntegrand(r, t) {
    return -this.variables.get("F0") + this.computeDPM_momentum_term(r) + this.computeDPM_gravity_term(r) + 
           this.computeDPM_stability_term() + this.computeLENR_term(r, t) + this.computeActivation_term(r, t) + 
           this.computeDE_term(r) + this.computeEM_term(r, t) + this.computeNeutron_term(r) + 
           this.computeRelativistic_term(r) + this.variables.get("F_Sweet_vac") + this.variables.get("F_Kozima");
  }
  
  computeForceIntegral(t, steps = 1000) {
    const x1 = this.variables.get("x1"), x2 = this.variables.get("x2");
    const dx = (x2 - x1) / steps;
    let sum = 0;
    for (let i = 0; i < steps; i++) {
      const r = x1 + (i + 0.5) * dx;
      sum += this.computeIntegrand(r, t) * dx;
    }
    return sum;
  }
  
  clone() { const m = new ButterflyNebulaUQFFModule(); m.variables = new Map(this.variables); return m; }
  setEnableLogging(e) { this.enableLogging = e; }
  registerDynamicTerm(t) { this.dynamicTerms.push(t); }
}

const domainExpansion = {
  expandNebulaScale(mF, rF) {
    this.variables.set("M", this.variables.get("M") * mF);
    this.variables.set("r", this.variables.get("r") * rF);
    this.variables.set("x2", this.variables.get("x2") * rF);
  },
  expandExpansionScale(vF, lF) {
    this.variables.set("V", this.variables.get("V") * vF);
    this.variables.set("L_x", this.variables.get("L_x") * lF);
  },
  expandLENRScale(kF, omegaF) {
    this.variables.set("k_LENR", this.variables.get("k_LENR") * kF);
    this.variables.set("omega_LENR", this.variables.get("omega_LENR") * omegaF);
  },
  optimizeForMetric(m) {
    const p = {"standard":{DPM_momentum:1.0,k_LENR:1.0},"simulation_fast":{DPM_momentum:1.0,k_LENR:1.0},"high_expansion":{V:10.0}}[m];
    if (p) Object.entries(p).forEach(([k,v]) => this.variables.has(k) && this.variables.set(k, v));
  }
};

addEnhancedDynamics(ButterflyNebulaUQFFModule, "Butterfly_Nebula", domainExpansion);
module.exports = ButterflyNebulaUQFFModule;
