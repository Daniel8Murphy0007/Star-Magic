// source131.js - SCM Velocity Module (ENHANCED + SIMULATION-READY)
// [SCm] Velocity (v_SCm) in UQFF Framework with self-expansion + parallel optimization
// Computes v_SCm = 1e8 m/s (~c/3); E_react = ρ_vac,[SCm] v_SCm² / ρ_vac,A * exp(-κ t)
// Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

const { addEnhancedDynamics } = require('./enhanced_dynamics.js');

class ScmVelocityModule {
  constructor() {
    // ========== CORE PARAMETERS (Original UQFF - Preserved Exactly) ==========
    this.variables = new Map([
      ["v_sc m", 1e8],                      // m/s (~c/3)
      ["rho_vac_SCm", 7.09e-37],           // J/m³
      ["rho_vac_A", 1e-23],                // J/m³
      ["kappa_day", 0.0005],               // day⁻¹
      ["day_to_s", 86400.0],               // s/day
      ["t_day", 0.0],                      // days
      ["mu_over_rj", 2.26e10],             // T m²
      ["P_SCm", 1.0],                      // Normalized
      ["heaviside_f", 1e11 + 1.0],         // 1 + 10^13 * 0.01
      ["quasi_f", 1.01],                   // 1 + 0.01
      ["one_minus_exp", 0.0]               // At t=0
    ]);
    
    // Derived values
    const v_scm = this.variables.get("v_sc m");
    const rho_SCm = this.variables.get("rho_vac_SCm");
    const rho_A = this.variables.get("rho_vac_A");
    this.variables.set("E_react_base", rho_SCm * Math.pow(v_scm, 2) / rho_A);
    this.variables.set("kappa_s", this.variables.get("kappa_day") / this.variables.get("day_to_s"));
    
    // ========== SELF-EXPANDING FRAMEWORK ==========
    this.dynamicParameters = new Map();
    this.dynamicTerms = [];
    this.metadata = new Map([["enhanced", true], ["version", "2.0.0"], ["simulation_ready", true]]);
    this.enableDynamicTerms = true;
    this.enableLogging = false;
    this.learningRate = 0.001;
  }
  
  // Update variable with automatic recalculation
  updateVariable(name, value) {
    if (this.variables.has(name)) {
      this.variables.set(name, value);
      if (name === "v_sc m" || name === "rho_vac_SCm" || name === "rho_vac_A") {
        const v_scm = this.variables.get("v_sc m");
        const rho_SCm = this.variables.get("rho_vac_SCm");
        const rho_A = this.variables.get("rho_vac_A");
        this.variables.set("E_react_base", rho_SCm * Math.pow(v_scm, 2) / rho_A);
      } else if (name === "kappa_day") {
        this.variables.set("kappa_s", value / this.variables.get("day_to_s"));
      }
    } else {
      this.enableLogging && console.log(`Adding variable '${name}' = ${value}`);
      this.variables.set(name, value);
    }
  }
  
  addToVariable(name, delta) {
    if (this.variables.has(name)) {
      this.variables.set(name, this.variables.get(name) + delta);
      if (name === "v_sc m" || name === "rho_vac_SCm" || name === "rho_vac_A") {
        const v_scm = this.variables.get("v_sc m"), rho_SCm = this.variables.get("rho_vac_SCm"), rho_A = this.variables.get("rho_vac_A");
        this.variables.set("E_react_base", rho_SCm * Math.pow(v_scm, 2) / rho_A);
      } else if (name === "kappa_day") {
        this.variables.set("kappa_s", this.variables.get("kappa_day") / this.variables.get("day_to_s"));
      }
    } else {
      this.variables.set(name, delta);
    }
  }
  
  subtractFromVariable(name, delta) { this.addToVariable(name, -delta); }
  
  // Core physics computations (EXACT from original)
  computeV_scm() { return this.variables.get("v_sc m"); }
  
  computeE_react(t_day) {
    this.variables.set("t_day", t_day);
    const arg = -this.variables.get("kappa_day") * t_day;
    return this.variables.get("E_react_base") * Math.exp(arg);
  }
  
  computeUmExample(t_day) {
    const e_react = this.computeE_react(t_day);
    const one_minus_exp = this.variables.get("one_minus_exp");
    const phi_hat = 1.0;
    const p_scm = this.variables.get("P_SCm");
    const heaviside_f = this.variables.get("heaviside_f");
    const quasi_f = this.variables.get("quasi_f");
    return (this.variables.get("mu_over_rj") * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
  }
  
  getEquationText() {
    return "E_react = [ρ_vac,[SCm] v_SCm² / ρ_vac,A] * exp(-κ t) (t days);\n" +
           "v_SCm = 1e8 m/s (~c/3, [SCm] propagation speed);\n" +
           "Scales reactivity in U_m, U_bi, U_i, U_gi via E_react.\n" +
           "Example t=0: E_react=1e46 J; t=2000 days: ~3.68e45 J (~36.8%).\n" +
           "U_m (t=0): ≈2.28e65 J/m³; t=2000: ≈8.39e64 J/m³.\n" +
           "Role: [SCm] dynamic speed for relativistic effects; jets/energy transfer.\n" +
           "UQFF: Subluminal propagation; [SCm]-[UA] reactions in nebulae/formation.";
  }
  
  printVariables() {
    console.log("Current Variables:");
    for (const [key, value] of this.variables) {
      console.log(`${key} = ${typeof value === 'number' ? value.toExponential() : value}`);
    }
  }
  
  printVelocityEffects(t_day = 2000.0) {
    const v = this.computeV_scm();
    const e_react = this.computeE_react(t_day);
    const um_ex = this.computeUmExample(t_day);
    const fraction = e_react / this.variables.get("E_react_base");
    console.log(`[SCm] Velocity Effects at t=${t_day} days:`);
    console.log(`v_SCm = ${v.toExponential()} m/s`);
    console.log(`E_react = ${e_react.toExponential()} J (${fraction.toFixed(4)} of initial)`);
    console.log(`U_m example = ${um_ex.toExponential()} J/m³`);
  }
  
  // Simulation-ready methods
  clone() {
    const m = new ScmVelocityModule();
    m.variables = new Map(this.variables);
    m.dynamicParameters = new Map(this.dynamicParameters);
    return m;
  }
  
  setEnableLogging(e) { this.enableLogging = e; }
  registerDynamicTerm(t) { this.dynamicTerms.push(t); }
}

const domainExpansion = {
  expandVelocityScale(vF, densityF) {
    this.variables.set("v_sc m", this.variables.get("v_sc m") * vF);
    this.variables.set("rho_vac_SCm", this.variables.get("rho_vac_SCm") * densityF);
    const v_scm = this.variables.get("v_sc m"), rho_SCm = this.variables.get("rho_vac_SCm"), rho_A = this.variables.get("rho_vac_A");
    this.variables.set("E_react_base", rho_SCm * Math.pow(v_scm, 2) / rho_A);
  },
  expandVacuumScale(rhoF, decayF) {
    this.variables.set("rho_vac_A", this.variables.get("rho_vac_A") * rhoF);
    this.variables.set("kappa_day", this.variables.get("kappa_day") * decayF);
    this.variables.set("kappa_s", this.variables.get("kappa_day") / this.variables.get("day_to_s"));
  },
  expandReactionScale(reactF) {
    this.variables.set("P_SCm", this.variables.get("P_SCm") * reactF);
    this.variables.set("heaviside_f", this.variables.get("heaviside_f") * reactF);
  },
  optimizeForMetric(m) {
    const p = {"standard":{v_scm:1e8,kappa_day:0.0005},"simulation_fast":{v_scm:1e8,kappa_day:0.0005},"high_velocity":{v_scm:2e8},"slow_decay":{kappa_day:0.0001}}[m];
    if (p) Object.entries(p).forEach(([k,v]) => this.variables.has(k) && this.updateVariable(k, v));
  }
};

addEnhancedDynamics(ScmVelocityModule, "SCM_Velocity", domainExpansion);
module.exports = ScmVelocityModule;
