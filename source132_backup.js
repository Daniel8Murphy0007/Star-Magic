// ButterflyNebulaUQFFModule.js
// JavaScript implementation of the UQFF Force for NGC 6302 (Butterfly Nebula).
// Computes F_U_Bi_i,enhanced as integral from x1 to x2 of [-F0 + DPM terms + LENR + activation + DE + EM + neutron + rel + Sweet + Kozima].
// Variables: M=0.64 M_sun, r=3.22e19 m, level=13; ~ -2.09e212 N at t=0
// Integrates LENR/resonance/buoyancy for nebula force; [SCm]/[UA] dynamics
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

class ButterflyNebulaUQFFModule {
  constructor() {
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    this.variables = new Map();
    
    // Universal constants
    this.variables.set("G", 6.6743e-11);                    // m^3 kg^-1 s^-2
    this.variables.set("c", 2.998e8);                       // m/s
    this.variables.set("m_e", 9.109e-31);                   // kg
    this.variables.set("hbar", 1.0546e-34);                 // J s
    this.variables.set("mu_B", 9.274e-24);                  // J/T
    this.variables.set("e", 1.602e-19);                     // C
    this.variables.set("M_sun", 1.989e30);                  // kg
    this.variables.set("q", 1.602e-19);                     // C
    this.variables.set("pi", 3.141592653589793);
    
    // Nebula-specific params
    this.variables.set("M", 0.64 * this.variables.get("M_sun"));  // kg (central star)
    this.variables.set("r", 3.22e19);                       // m (distance)
    this.variables.set("x1", 0.0);                          // m (integral lower)
    this.variables.set("x2", 3.22e19);                      // m (upper)
    this.variables.set("level", 13.0);                      // Quantum level
    this.variables.set("F0", 1.0);                          // Base force (normalized)
    this.variables.set("theta", 0.0);                       // rad (angle)
    this.variables.set("DPM_momentum", 1.0);                // Normalized
    this.variables.set("DPM_gravity", 1.0);                 // Normalized
    this.variables.set("DPM_stability", 0.01);              // Normalized
    this.variables.set("rho_vac_UA", 7.09e-36);             // J/m³
    this.variables.set("k_LENR", 1.0);                      // Coupling
    this.variables.set("omega_LENR", 7.85e12);              // Hz
    this.variables.set("omega_0", 1e-12);                   // Hz (reference)
    this.variables.set("k_act", 1.0);                       // Activation coupling
    this.variables.set("omega_act", 1.0);                   // rad/s
    this.variables.set("k_DE", 1.0);                        // DE coupling
    this.variables.set("L_x", 1.0);                         // Length scale
    this.variables.set("B_0", 1.0);                         // T
    this.variables.set("V", 1.0);                           // m/s
    this.variables.set("g", 9.8);                           // m/s²
    this.variables.set("k_neutron", 1e10);                  // Neutron coupling
    this.variables.set("sigma_n", 1e-4);                    // Barn
    this.variables.set("k_rel", 1.0);                       // Rel coupling
    this.variables.set("E_cm", 1.0);                        // eV
    this.variables.set("E_cm_eff", 1.0);                    // Enhanced eV
    this.variables.set("F_Sweet_vac", 7.09e-39);            // N (negligible)
    this.variables.set("F_Kozima", 7.85e30);                // N
    this.variables.set("t", 0.0);                           // s
    
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
  
  // Update variable
  updateVariable(name, value) {
    this.variables.set(name, value);
  }
  
  // Add delta to variable
  addToVariable(name, delta) {
    if (this.variables.has(name)) {
      this.variables.set(name, this.variables.get(name) + delta);
    } else {
      console.error(`Variable '${name}' not found. Adding with delta ${delta}`);
      this.variables.set(name, delta);
    }
  }
  
  // Subtract delta from variable
  subtractFromVariable(name, delta) {
    this.addToVariable(name, -delta);
  }
  
  // DPM momentum term
  computeDPM_momentum_term(r) {
    const m_e_c2 = this.variables.get("m_e") * Math.pow(this.variables.get("c"), 2);
    return (m_e_c2 / (r * r)) * this.variables.get("DPM_momentum") * Math.cos(this.variables.get("theta"));
  }
  
  // DPM gravity term
  computeDPM_gravity_term(r) {
    return (this.variables.get("G") * this.variables.get("M") / (r * r)) * this.variables.get("DPM_gravity");
  }
  
  // DPM stability term
  computeDPM_stability_term() {
    return this.variables.get("rho_vac_UA") * this.variables.get("DPM_stability");
  }
  
  // LENR term
  computeLENR_term() {
    const ratio = Math.pow(this.variables.get("omega_LENR") / this.variables.get("omega_0"), 2);
    return this.variables.get("k_LENR") * ratio;
  }
  
  // Activation term
  computeActivation_term(t) {
    return this.variables.get("k_act") * Math.cos(this.variables.get("omega_act") * t);
  }
  
  // DE term
  computeDE_term(L_x) {
    return this.variables.get("k_DE") * L_x;
  }
  
  // EM term
  computeEM_term() {
    const q_v_B = 2 * this.variables.get("q") * this.variables.get("B_0") * this.variables.get("V") * Math.sin(this.variables.get("theta"));
    const g_mu_B = this.variables.get("g") * this.variables.get("mu_B") * this.variables.get("B_0") / (this.variables.get("hbar") * this.variables.get("omega_0"));
    return q_v_B * g_mu_B;
  }
  
  // Neutron term
  computeNeutron_term() {
    return this.variables.get("k_neutron") * this.variables.get("sigma_n");
  }
  
  // Relativistic term
  computeRel_term(E_cm_eff) {
    const ratio = Math.pow(E_cm_eff / this.variables.get("E_cm"), 2);
    return this.variables.get("k_rel") * ratio;
  }
  
  // Sweet vacuum term
  computeSweet_vac_term() {
    return this.variables.get("F_Sweet_vac");
  }
  
  // Kozima term
  computeKozima_term() {
    return this.variables.get("F_Kozima");
  }
  
  // Full integrand at position x and time t
  computeIntegrand(x, t) {
    return -this.variables.get("F0") + 
           this.computeDPM_momentum_term(x) + 
           this.computeDPM_gravity_term(x) + 
           this.computeDPM_stability_term() +
           this.computeLENR_term() + 
           this.computeActivation_term(t) + 
           this.computeDE_term(this.variables.get("L_x")) + 
           this.computeEM_term() +
           this.computeNeutron_term() + 
           this.computeRel_term(this.variables.get("E_cm_eff")) + 
           this.computeSweet_vac_term() + 
           this.computeKozima_term();
  }
  
  // Numerical integral using trapezoidal rule
  computeIntegral(x1, x2, t, n_points = 1000) {
    const dx = (x2 - x1) / n_points;
    let integral = 0.0;
    for (let i = 0; i <= n_points; i++) {
      const x = x1 + i * dx;
      const weight = (i === 0 || i === n_points) ? 0.5 : 1.0;
      integral += weight * this.computeIntegrand(x, t);
    }
    return integral * dx;
  }
  
  // Main F_U_Bi_i,enhanced computation
  computeF_U_Bi(x1, x2, t) {
    return this.computeIntegral(x1, x2, t);
  }
  
  // Get equation text
  getEquationText() {
    return "F_U_Bi_i,enhanced = ∫_{x1}^{x2} [-F0 + (m_e c^2 / r^2) DPM_mom cosθ + (G M / r^2) DPM_grav + ρ_[UA] DPM_stab + k_LENR (ω_LENR/ω_0)^2 + k_act cos(ω_act t) + k_DE L_x + 2 q B_0 V sinθ (g μ_B B_0 / ℏ ω_0) + k_neutron σ_n + k_rel (E_cm,eff / E_cm)^2 + F_Sweet,vac + F_Kozima] dx\n" +
           "NGC 6302: M=0.64 M_sun, r=3.22e19 m, level=13; ~ -2.09e212 N (repulsive stabilization).\n" +
           "Sweet: ρ_[UA] DPM_stab V ≈7.09e-39 N (negligible); Kozima: k_n σ_n (ω_LENR/ω_0) ≈7.85e30 N.\n" +
           "UQFF: Integrates LENR/resonance/buoyancy for nebula force; [SCm]/[UA] dynamics.";
  }
  
  // Print all current variables
  printVariables() {
    console.log("NGC 6302 Variables:");
    for (const [key, value] of this.variables) {
      console.log(`${key} = ${value.toExponential()}`);
    }
  }
}

module.exports = ButterflyNebulaUQFFModule;
