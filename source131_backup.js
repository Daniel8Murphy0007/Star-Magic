// ScmVelocityModule.js
// JavaScript implementation of the [SCm] Velocity (v_SCm) in the Universal Quantum Field Superconductive Framework (UQFF).
// Computes v_SCm = 1e8 m/s (~c/3); scales in E_react = ρ_vac,[SCm] v_SCm² / ρ_vac,A * exp(-κ t) for U_m, U_bi, etc.
// Variables: v_sc m=1e8 m/s, ρ_vac,[SCm]=7.09e-37 J/m³, ρ_vac,A=1e-23 J/m³, κ=0.0005 day⁻¹
// Role: [SCm] dynamic speed for relativistic effects; jets/energy transfer in UQFF framework
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

class ScmVelocityModule {
  constructor() {
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    this.variables = new Map();
    
    // Universal constants
    this.variables.set("v_sc m", 1e8);                     // m/s (~c/3)
    this.variables.set("rho_vac_SCm", 7.09e-37);           // J/m³
    this.variables.set("rho_vac_A", 1e-23);                // J/m³
    this.variables.set("kappa_day", 0.0005);               // day⁻¹
    this.variables.set("day_to_s", 86400.0);               // s/day
    this.variables.set("t_day", 0.0);                      // days
    
    // Derived - calculate E_react_base
    const v_scm = this.variables.get("v_sc m");
    const rho_SCm = this.variables.get("rho_vac_SCm");
    const rho_A = this.variables.get("rho_vac_A");
    this.variables.set("E_react_base", rho_SCm * Math.pow(v_scm, 2) / rho_A);
    
    this.variables.set("mu_over_rj", 2.26e10);             // T m² (example)
    this.variables.set("P_SCm", 1.0);                      // Normalized
    this.variables.set("heaviside_f", 1e11 + 1.0);         // 1 + 10^13 * 0.01
    this.variables.set("quasi_f", 1.01);                   // 1 + 0.01
    this.variables.set("one_minus_exp", 0.0);              // At t=0
    
    // Derived kappa_s
    this.variables.set("kappa_s", this.variables.get("kappa_day") / this.variables.get("day_to_s"));
    
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    this.dynamicParameters = new Map();
    this.dynamicTerms = [];
    this.metadata = new Map();
    this.enableDynamicTerms = true;
    this.enableLogging = false;
    this.learningRate = 0.001;
    this.metadata.set("enhanced", "true");
    this.metadata.set("version", "2.0-Enhanced");
  }
  
  // Update variable with automatic recalculation of derived variables
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
      console.error(`Variable '${name}' not found. Adding with value ${value}`);
      this.variables.set(name, value);
    }
  }
  
  // Add delta to variable
  addToVariable(name, delta) {
    if (this.variables.has(name)) {
      this.variables.set(name, this.variables.get(name) + delta);
      if (name === "v_sc m" || name === "rho_vac_SCm" || name === "rho_vac_A") {
        const v_scm = this.variables.get("v_sc m");
        const rho_SCm = this.variables.get("rho_vac_SCm");
        const rho_A = this.variables.get("rho_vac_A");
        this.variables.set("E_react_base", rho_SCm * Math.pow(v_scm, 2) / rho_A);
      } else if (name === "kappa_day") {
        this.variables.set("kappa_s", this.variables.get("kappa_day") / this.variables.get("day_to_s"));
      }
    } else {
      console.error(`Variable '${name}' not found. Adding with delta ${delta}`);
      this.variables.set(name, delta);
    }
  }
  
  // Subtract delta from variable
  subtractFromVariable(name, delta) {
    this.addToVariable(name, -delta);
  }
  
  // Compute v_SCm (m/s)
  computeV_scm() {
    return this.variables.get("v_sc m");
  }
  
  // Compute E_react = E_react_base * exp(-κ t)
  computeE_react(t_day) {
    this.variables.set("t_day", t_day);
    const arg = -this.variables.get("kappa_day") * t_day;
    return this.variables.get("E_react_base") * Math.exp(arg);
  }
  
  // Simplified U_m example with E_react
  computeUmExample(t_day) {
    const e_react = this.computeE_react(t_day);
    const one_minus_exp = this.variables.get("one_minus_exp");
    const phi_hat = 1.0;
    const p_scm = this.variables.get("P_SCm");
    const heaviside_f = this.variables.get("heaviside_f");
    const quasi_f = this.variables.get("quasi_f");
    return (this.variables.get("mu_over_rj") * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
  }
  
  // Get equation text
  getEquationText() {
    return "E_react = [ρ_vac,[SCm] v_SCm² / ρ_vac,A] * exp(-κ t) (t days);\n" +
           "v_SCm = 1e8 m/s (~c/3, [SCm] propagation speed);\n" +
           "Scales reactivity in U_m, U_bi, U_i, U_gi via E_react.\n" +
           "Example t=0: E_react=1e46 J; t=2000 days: ~3.68e45 J (~36.8%).\n" +
           "U_m (t=0): ≈2.28e65 J/m³; t=2000: ≈8.39e64 J/m³.\n" +
           "Role: [SCm] dynamic speed for relativistic effects; jets/energy transfer.\n" +
           "UQFF: Subluminal propagation; [SCm]-[UA] reactions in nebulae/formation.";
  }
  
  // Print all current variables
  printVariables() {
    console.log("Current Variables:");
    for (const [key, value] of this.variables) {
      console.log(`${key} = ${value.toExponential()}`);
    }
  }
  
  // Print velocity effects at specified time
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
}

module.exports = ScmVelocityModule;
