// Wolfram-Enhanced Physics Terms from source82.cpp
// Generated: November 27, 2025
// Source: SMBHUQFFModule - SMBH M-Ïƒ Relation UQFF Integration (372 lines)
// Total Classes: 15 (SMBH-specific terms from DynamicVacuum, QuantumCoupling, and MUGE framework)

#include <cmath>
#include <string>
#include <map>
#include <memory>
#include <vector>

// ============================================================================
// BASE PHYSICS TERM INTERFACE
// ============================================================================

class PhysicsTerm {
public:
    virtual ~PhysicsTerm() {}
    virtual double compute(double t, const std::map<std::string, double>& params) const = 0;
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
    virtual bool validate(const std::map<std::string, double>& params) const { return true; }
};

// ============================================================================
// SOURCE82-SPECIFIC CLASSES (SMBH M-Ïƒ RELATION TERMS)
// ============================================================================

class SMBHDynamicVacuumTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // SMBH-specific vacuum energy oscillation
        // E_vac(t) = amplitude * rho_vac * sin(frequency * t)
        auto it_amp = params.find("amplitude");
        auto it_rho = params.find("rho_vac_UA");
        auto it_freq = params.find("frequency");
        
        double amplitude = (it_amp != params.end()) ? it_amp->second : 1e-10;
        double rho_vac = (it_rho != params.end()) ? it_rho->second : 7.09e-36; // kg/mÂ³
        double frequency = (it_freq != params.end()) ? it_freq->second : 1e-15; // Hz
        
        return amplitude * rho_vac * std::sin(frequency * t);
    }
    
    std::string getName() const override { return "SMBH_DynamicVacuum"; }
    
    std::string getDescription() const override {
        return "SMBH time-varying vacuum energy: E_vac = amp * rho_vac * sin(freq*t)";
    }
    
    bool validate(const std::map<std::string, double>& params) const override {
        return params.count("rho_vac_UA") && params.at("rho_vac_UA") > 0;
    }
};

class SMBHQuantumCouplingTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Quantum coupling correction for SMBH
        // F_q = coupling_strength * (M_bh / M_sun)^(1/3) * exp(-r/lambda_Q)
        auto it_coupling = params.find("coupling_strength");
        auto it_M_bh = params.find("M_bh");
        auto it_M_sun = params.find("M_sun");
        auto it_r = params.find("r");
        auto it_lambda = params.find("lambda_Q");
        
        double coupling = (it_coupling != params.end()) ? it_coupling->second : 1.23e-40;
        double M_bh = (it_M_bh != params.end()) ? it_M_bh->second : 1e11; // Solar masses
        double M_sun = (it_M_sun != params.end()) ? it_M_sun->second : 1.989e30; // kg
        double r = (it_r != params.end()) ? it_r->second : 1e3; // parsecs
        double lambda_Q = (it_lambda != params.end()) ? it_lambda->second : 1e-10; // m
        
        double mass_ratio = std::pow(M_bh / M_sun, 1.0/3.0);
        return coupling * mass_ratio * std::exp(-r / lambda_Q);
    }
    
    std::string getName() const override { return "SMBH_QuantumCoupling"; }
    
    std::string getDescription() const override {
        return "SMBH quantum coupling: F_q = g_coupling * (M_bh/M_sun)^(1/3) * exp(-r/lambda)";
    }
    
    bool validate(const std::map<std::string, double>& params) const override {
        return params.count("M_bh") && params.at("M_bh") > 0;
    }
};

class SMBHMSigmaRelationTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // M-Ïƒ relation: M_bh = M_0 * (sigma / sigma_0)^alpha
        // Canonical: M_0 = 1.3e8 M_sun, sigma_0 = 200 km/s, alpha = 4.24
        auto it_M0 = params.find("M_0");
        auto it_sigma = params.find("sigma");
        auto it_sigma0 = params.find("sigma_0");
        auto it_alpha = params.find("alpha_M_sigma");
        
        double M_0 = (it_M0 != params.end()) ? it_M0->second : 1.3e8; // M_sun
        double sigma = (it_sigma != params.end()) ? it_sigma->second : 300.0; // km/s
        double sigma_0 = (it_sigma0 != params.end()) ? it_sigma0->second : 200.0; // km/s
        double alpha = (it_alpha != params.end()) ? it_alpha->second : 4.24;
        
        return M_0 * std::pow(sigma / sigma_0, alpha);
    }
    
    std::string getName() const override { return "SMBH_MSigmaRelation"; }
    
    std::string getDescription() const override {
        return "M-Ïƒ relation: M_bh = M_0 * (sigma/sigma_0)^4.24";
    }
    
    bool validate(const std::map<std::string, double>& params) const override {
        return params.count("sigma") && params.at("sigma") > 0;
    }
};

class SMBHBulgeGravityTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Bulge gravitational potential influence on SMBH
        // Î¦_bulge = -G * M_bulge / R_bulge
        auto it_G = params.find("G");
        auto it_M_bulge = params.find("M_bulge");
        auto it_R_bulge = params.find("R_bulge");
        
        double G = (it_G != params.end()) ? it_G->second : 6.674e-11; // mÂ³/kg/sÂ²
        double M_bulge = (it_M_bulge != params.end()) ? it_M_bulge->second : 1e11 * 1.989e30; // kg
        double R_bulge = (it_R_bulge != params.end()) ? it_R_bulge->second : 1e3 * 3.086e16; // m
        
        return -G * M_bulge / R_bulge;
    }
    
    std::string getName() const override { return "SMBH_BulgeGravity"; }
    
    std::string getDescription() const override {
        return "Bulge gravitational potential: Î¦ = -G*M_bulge/R_bulge";
    }
    
    bool validate(const std::map<std::string, double>& params) const override {
        return params.count("R_bulge") && params.at("R_bulge") > 0;
    }
};

class SMBHUg1Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Ug1: Primary unified gravity component (Newtonian base)
        // Ug1 = G * M_bh / r^2
        auto it_G = params.find("G");
        auto it_M_bh = params.find("M_bh");
        auto it_r = params.find("r");
        
        double G = (it_G != params.end()) ? it_G->second : 6.674e-11;
        double M_bh = (it_M_bh != params.end()) ? it_M_bh->second : 1e11 * 1.989e30; // kg
        double r = (it_r != params.end()) ? it_r->second : 1e3 * 3.086e16; // m
        
        return G * M_bh / (r * r);
    }
    
    std::string getName() const override { return "SMBH_Ug1"; }
    
    std::string getDescription() const override {
        return "SMBH Ug1 (Newtonian): G*M_bh/r^2";
    }
};

class SMBHUg2Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Ug2: General relativistic correction
        // Ug2 = Ug1 * (3 * G * M_bh) / (r * c^2)
        auto it_G = params.find("G");
        auto it_M_bh = params.find("M_bh");
        auto it_r = params.find("r");
        auto it_c = params.find("c");
        
        double G = (it_G != params.end()) ? it_G->second : 6.674e-11;
        double M_bh = (it_M_bh != params.end()) ? it_M_bh->second : 1e11 * 1.989e30;
        double r = (it_r != params.end()) ? it_r->second : 1e3 * 3.086e16;
        double c = (it_c != params.end()) ? it_c->second : 2.998e8; // m/s
        
        double Ug1 = G * M_bh / (r * r);
        return Ug1 * (3.0 * G * M_bh) / (r * c * c);
    }
    
    std::string getName() const override { return "SMBH_Ug2"; }
    
    std::string getDescription() const override {
        return "SMBH Ug2 (GR correction): Ug1 * 3GM/(rc^2)";
    }
};

class SMBHUg3Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Ug3: Quantum gravity correction (monopole-like)
        // Ug3 = k_monopole * exp(-r/lambda_Q)
        auto it_k = params.find("k_monopole");
        auto it_r = params.find("r");
        auto it_lambda = params.find("lambda_Q");
        
        double k = (it_k != params.end()) ? it_k->second : 1e-50;
        double r = (it_r != params.end()) ? it_r->second : 1e3 * 3.086e16;
        double lambda_Q = (it_lambda != params.end()) ? it_lambda->second : 1e-10;
        
        return k * std::exp(-r / lambda_Q);
    }
    
    std::string getName() const override { return "SMBH_Ug3"; }
    
    std::string getDescription() const override {
        return "SMBH Ug3 (quantum gravity): k_monopole * exp(-r/lambda_Q)";
    }
};

class SMBHUg4Term : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Ug4: Vacuum energy contribution
        // Ug4 = rho_vac * c^2 * (r^3) / (3 * M_bh)
        auto it_rho = params.find("rho_vac_UA");
        auto it_c = params.find("c");
        auto it_r = params.find("r");
        auto it_M_bh = params.find("M_bh");
        
        double rho_vac = (it_rho != params.end()) ? it_rho->second : 7.09e-36;
        double c = (it_c != params.end()) ? it_c->second : 2.998e8;
        double r = (it_r != params.end()) ? it_r->second : 1e3 * 3.086e16;
        double M_bh = (it_M_bh != params.end()) ? it_M_bh->second : 1e11 * 1.989e30;
        
        return rho_vac * c * c * (r * r * r) / (3.0 * M_bh);
    }
    
    std::string getName() const override { return "SMBH_Ug4"; }
    
    std::string getDescription() const override {
        return "SMBH Ug4 (vacuum energy): rho_vac*c^2*r^3/(3*M_bh)";
    }
};

class SMBHReactorEfficiencyTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Reactor energy efficiency (exponential decay model)
        // E_react(t) = E_0 * exp(-t / tau_decay)
        auto it_E0 = params.find("E_react_0");
        auto it_tau = params.find("tau_decay");
        
        double E_0 = (it_E0 != params.end()) ? it_E0->second : 1e52; // J
        double tau = (it_tau != params.end()) ? it_tau->second : 4.543e9 * 3.156e7; // yr -> s
        
        return E_0 * std::exp(-t / tau);
    }
    
    std::string getName() const override { return "SMBH_ReactorEfficiency"; }
    
    std::string getDescription() const override {
        return "Reactor efficiency decay: E_react = E_0 * exp(-t/tau)";
    }
};

class SMBHPseudoMonopoleTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Pseudo-monopole shift (state-dependent)
        // delta_n(state) for states 1-26
        auto it_state = params.find("state_n");
        auto it_shift = params.find("delta_shift");
        
        double state_n = (it_state != params.end()) ? it_state->second : 1.0;
        double shift = (it_shift != params.end()) ? it_shift->second : 1e-30;
        
        // Simplified: linear in state number
        return shift * state_n;
    }
    
    std::string getName() const override { return "SMBH_PseudoMonopole"; }
    
    std::string getDescription() const override {
        return "Pseudo-monopole state shift: delta_n(state)";
    }
};

class SMBHRedshiftCorrectionTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Cosmological redshift correction for observed mass
        // M_obs = M_intrinsic * (1 + z)
        auto it_M = params.find("M_bh");
        auto it_z = params.find("redshift");
        
        double M_bh = (it_M != params.end()) ? it_M->second : 1e11;
        double z = (it_z != params.end()) ? it_z->second : 0.0; // z=0-6
        
        return M_bh * (1.0 + z);
    }
    
    std::string getName() const override { return "SMBH_RedshiftCorrection"; }
    
    std::string getDescription() const override {
        return "Redshift mass correction: M_obs = M_intrinsic * (1+z)";
    }
};

class SMBHUiTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Ui: Internal energy term
        // Ui = 3/2 * k_B * T_virial * (M_bh / m_p)
        auto it_k = params.find("k_B");
        auto it_T = params.find("T_virial");
        auto it_M_bh = params.find("M_bh");
        auto it_mp = params.find("m_proton");
        
        double k_B = (it_k != params.end()) ? it_k->second : 1.381e-23; // J/K
        double T = (it_T != params.end()) ? it_T->second : 1e7; // K
        double M_bh = (it_M_bh != params.end()) ? it_M_bh->second : 1e11 * 1.989e30; // kg
        double m_p = (it_mp != params.end()) ? it_mp->second : 1.673e-27; // kg
        
        return 1.5 * k_B * T * (M_bh / m_p);
    }
    
    std::string getName() const override { return "SMBH_Ui"; }
    
    std::string getDescription() const override {
        return "SMBH internal energy: Ui = 3/2 * k_B * T * (M_bh/m_p)";
    }
};

class SMBHUmTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Um: Magnetic energy term (26-layer compressed gravity)
        // Um = B^2 / (2 * mu_0) * V_bulge
        auto it_B = params.find("B_field");
        auto it_mu0 = params.find("mu_0");
        auto it_R = params.find("R_bulge");
        
        double B = (it_B != params.end()) ? it_B->second : 1e-6; // Tesla
        double mu_0 = (it_mu0 != params.end()) ? it_mu0->second : 4.0 * M_PI * 1e-7; // H/m
        double R = (it_R != params.end()) ? it_R->second : 1e3 * 3.086e16; // m
        
        double V_bulge = (4.0 / 3.0) * M_PI * R * R * R;
        return (B * B) / (2.0 * mu_0) * V_bulge;
    }
    
    std::string getName() const override { return "SMBH_Um"; }
    
    std::string getDescription() const override {
        return "SMBH magnetic energy: Um = B^2/(2*mu_0) * V_bulge";
    }
};

class SMBHOmegaSGalacticTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // omega_s: Galactic-scale angular velocity
        // omega_s = v_circ / R_bulge
        auto it_v = params.find("v_circ");
        auto it_R = params.find("R_bulge");
        
        double v_circ = (it_v != params.end()) ? it_v->second : 220e3; // m/s (220 km/s)
        double R = (it_R != params.end()) ? it_R->second : 1e3 * 3.086e16; // m
        
        return v_circ / R;
    }
    
    std::string getName() const override { return "SMBH_OmegaS"; }
    
    std::string getDescription() const override {
        return "Galactic angular velocity: omega_s = v_circ / R_bulge";
    }
};

class SMBHCosmicTimeTerm : public PhysicsTerm {
public:
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Cosmic time approximation (age of universe at redshift z)
        // t_cosmic â‰ˆ t_0 / (1 + z)^(3/2) for matter-dominated
        auto it_t0 = params.find("t_universe");
        auto it_z = params.find("redshift");
        
        double t_0 = (it_t0 != params.end()) ? it_t0->second : 4.543e9 * 3.156e7; // 4.543 Gyr in seconds
        double z = (it_z != params.end()) ? it_z->second : 0.0;
        
        return t_0 / std::pow(1.0 + z, 1.5);
    }
    
    std::string getName() const override { return "SMBH_CosmicTime"; }
    
    std::string getDescription() const override {
        return "Cosmic time: t_cosmic = t_0 / (1+z)^(3/2)";
    }
};

// ============================================================================
// WOLFRAM EXPORT UTILITIES
// ============================================================================

std::string exportToWolfram_Source82() {
    std::string wolfram_code;
    
    wolfram_code += "(* Wolfram Language Export - source82_wolfram.cpp *)\n";
    wolfram_code += "(* SMBH M-Ïƒ Relation UQFF Integration - 15 Physics Terms *)\n";
    wolfram_code += "(* Generated: November 27, 2025 *)\n\n";
    
    wolfram_code += "source82Sector = {\n";
    wolfram_code += "  (* Dynamic Vacuum Energy *)\n";
    wolfram_code += "  SMBHDynamicVacuum[t_, amp_, rhoVac_, freq_] := amp * rhoVac * Sin[freq * t],\n\n";
    
    wolfram_code += "  (* Quantum Coupling *)\n";
    wolfram_code += "  SMBHQuantumCoupling[Mbh_, Msun_, r_, lambda_, g_] := g * (Mbh/Msun)^(1/3) * Exp[-r/lambda],\n\n";
    
    wolfram_code += "  (* M-Ïƒ Relation *)\n";
    wolfram_code += "  SMBHMSigmaRelation[sigma_, M0_, sigma0_, alpha_] := M0 * (sigma/sigma0)^alpha,\n\n";
    
    wolfram_code += "  (* Bulge Gravity Potential *)\n";
    wolfram_code += "  SMBHBulgeGravity[G_, Mbulge_, Rbulge_] := -G * Mbulge / Rbulge,\n\n";
    
    wolfram_code += "  (* Unified Gravity Components *)\n";
    wolfram_code += "  SMBHUg1[G_, Mbh_, r_] := G * Mbh / r^2,\n";
    wolfram_code += "  SMBHUg2[G_, Mbh_, r_, c_] := (G * Mbh / r^2) * (3 * G * Mbh) / (r * c^2),\n";
    wolfram_code += "  SMBHUg3[k_, r_, lambda_] := k * Exp[-r/lambda],\n";
    wolfram_code += "  SMBHUg4[rhoVac_, c_, r_, Mbh_] := rhoVac * c^2 * r^3 / (3 * Mbh),\n\n";
    
    wolfram_code += "  (* Reactor Efficiency *)\n";
    wolfram_code += "  SMBHReactorEfficiency[t_, E0_, tau_] := E0 * Exp[-t/tau],\n\n";
    
    wolfram_code += "  (* Pseudo-Monopole Shift *)\n";
    wolfram_code += "  SMBHPseudoMonopole[stateN_, shift_] := shift * stateN,\n\n";
    
    wolfram_code += "  (* Redshift Correction *)\n";
    wolfram_code += "  SMBHRedshiftCorrection[Mbh_, z_] := Mbh * (1 + z),\n\n";
    
    wolfram_code += "  (* Internal Energy *)\n";
    wolfram_code += "  SMBHUi[kB_, T_, Mbh_, mp_] := 3/2 * kB * T * (Mbh/mp),\n\n";
    
    wolfram_code += "  (* Magnetic Energy *)\n";
    wolfram_code += "  SMBHUm[B_, mu0_, Rbulge_] := (B^2 / (2 * mu0)) * (4/3 * Pi * Rbulge^3),\n\n";
    
    wolfram_code += "  (* Galactic Angular Velocity *)\n";
    wolfram_code += "  SMBHOmegaS[vCirc_, Rbulge_] := vCirc / Rbulge,\n\n";
    
    wolfram_code += "  (* Cosmic Time *)\n";
    wolfram_code += "  SMBHCosmicTime[t0_, z_] := t0 / (1 + z)^(3/2)\n";
    wolfram_code += "};\n\n";
    
    wolfram_code += "(* Master SMBH UQFF Equation *)\n";
    wolfram_code += "MasterSMBHUQFF[params_] := Module[{},\n";
    wolfram_code += "  (* Sum all components *)\n";
    wolfram_code += "  Total[{\n";
    wolfram_code += "    SMBHUg1[params[\"G\"], params[\"Mbh\"], params[\"r\"]],\n";
    wolfram_code += "    SMBHUg2[params[\"G\"], params[\"Mbh\"], params[\"r\"], params[\"c\"]],\n";
    wolfram_code += "    SMBHUg3[params[\"k_monopole\"], params[\"r\"], params[\"lambda_Q\"]],\n";
    wolfram_code += "    SMBHUg4[params[\"rho_vac\"], params[\"c\"], params[\"r\"], params[\"Mbh\"]],\n";
    wolfram_code += "    SMBHUi[params[\"kB\"], params[\"T\"], params[\"Mbh\"], params[\"mp\"]],\n";
    wolfram_code += "    SMBHUm[params[\"B\"], params[\"mu0\"], params[\"Rbulge\"]]\n";
    wolfram_code += "  }]\n";
    wolfram_code += "];\n\n";
    
    wolfram_code += "(* End source82_wolfram.cpp export *)\n";
    
    return wolfram_code;
}

// End of source82_wolfram.cpp

// ===== VIRGO CLUSTER COSMOLOGICAL TERMS =====

class VirgoClusterMassTerm
{
private:
    double G;              // Gravitational constant (m┬│/kg┬╖s┬▓)
    double M_cluster;      // Total cluster mass (kg, ~1.2e15 M_sun)
    double R_virial;       // Virial radius (m, ~2.2 Mpc = 6.79e22 m)
    double z;              // Redshift (0.0036 for Virgo center)

public:
    VirgoClusterMassTerm(double mass = 1.2e15 * M_SUN, double r_vir = 2.2 * MPC_TO_M, double redshift = 0.0036)
        : G(G_CONST), M_cluster(mass), R_virial(r_vir), z(redshift) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // Enclosed mass profile: M(<r) = M_cluster * (r/R_virial)^3 * (1 + (R_virial/r))^(-2)
        // Simplified NFW-like profile for cluster
        double x = r / R_virial;
        double M_enclosed = M_cluster * (x * x * x) / (1.0 + x) / (1.0 + x);
        
        // Return gravitational acceleration: a = G * M_enclosed / r┬▓
        return (G * M_enclosed) / (r * r);
    }

    std::string toWolfram() const
    {
        return "VirgoClusterMass[r_, Mcluster_: 1.2*^15 * 1.989*^30, Rvirial_: 2.2 * 3.086*^22, G_: 6.6743*^-11] := "
               "Module[{x, Menclosed}, "
               "x = r / Rvirial; "
               "Menclosed = Mcluster * x^3 / (1 + x)^2; "
               "G * Menclosed / r^2]";
    }

    std::string getSignature() const { return "VirgoClusterMassTerm(r, params)"; }
    std::string getCategory() const { return "gravity"; }
    std::string getName() const { return "VirgoClusterMass"; }
    std::string getDescription() const { return "Virgo Cluster total mass gravitational acceleration: a = G┬╖M(<r)/r┬▓ with M_cluster~1.2e15 M_sun"; }
};

// ========================================
// CLASS 821: VirgoClusterIntraclusterMediumTerm
// Category: gas_dynamics
// Physics: Intracluster medium (ICM) - hot X-ray emitting gas at T~2-4 keV
// ICM contains ~15% of cluster baryonic mass, creates X-ray emission and SZ effect
// ========================================
class VirgoClusterIntraclusterMediumTerm
{
private:
    double T_ICM;          // ICM temperature (keV, ~2.3 keV for Virgo)
    double n_e0;           // Central electron density (mΓü╗┬│, ~3e3)
    double r_c;            // Core radius (m, ~40 kpc = 1.23e21 m)
    double beta;           // Beta model parameter (~0.5 for Virgo)

public:
    VirgoClusterIntraclusterMediumTerm(double temp = 2.3, double ne_central = 3e3, double core_r = 40 * KPC_TO_M, double beta_param = 0.5)
        : T_ICM(temp), n_e0(ne_central), r_c(core_r), beta(beta_param) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // Beta model: n_e(r) = n_e0 * (1 + (r/r_c)┬▓)^(-3╬▓/2)
        double n_e = n_e0 * std::pow(1.0 + (r * r) / (r_c * r_c), -1.5 * beta);
        
        // Pressure: P = n_e * k_B * T (convert T from keV to Kelvin: T_K = T_keV * 1.16e7)
        double T_K = T_ICM * 1.16e7;
        double P_ICM = n_e * K_BOLTZ * T_K;
        
        return P_ICM;
    }

    std::string toWolfram() const
    {
        return "VirgoICM[r_, T_: 2.3, ne0_: 3000, rc_: 40 * 3.086*^19, beta_: 0.5] := "
               "Module[{ne, TK, kB}, "
               "kB = 1.381*^-23; "
               "ne = ne0 * (1 + (r/rc)^2)^(-1.5 * beta); "
               "TK = T * 1.16*^7; "
               "ne * kB * TK]";
    }

    std::string getSignature() const { return "VirgoClusterIntraclusterMediumTerm(r, params)"; }
    std::string getCategory() const { return "gas_dynamics"; }
    std::string getName() const { return "VirgoClusterICM"; }
    std::string getDescription() const { return "ICM beta-model pressure: P = n_e(r)┬╖k_B┬╖T with T~2.3 keV, ╬▓~0.5"; }
};

// ========================================
// CLASS 822: VirgoClusterGravitationalPotentialTerm
// Category: gravity
// Physics: Cluster gravitational potential well ╬ª(r) = -G┬╖M(<r)/r
// Determines escape velocity and binding energy of member galaxies
// ========================================
class VirgoClusterGravitationalPotentialTerm
{
private:
    double G;              // Gravitational constant
    double M_cluster;      // Total cluster mass
    double R_virial;       // Virial radius
    double c_NFW;          // NFW concentration parameter (~4 for clusters)

public:
    VirgoClusterGravitationalPotentialTerm(double mass = 1.2e15 * M_SUN, double r_vir = 2.2 * MPC_TO_M, double concentration = 4.0)
        : G(G_CONST), M_cluster(mass), R_virial(r_vir), c_NFW(concentration) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // NFW potential: ╬ª(r) = -G┬╖M_virial ┬╖ ln(1 + c┬╖x) / (x ┬╖ f(c))
        // where x = r/R_virial, f(c) = ln(1+c) - c/(1+c)
        double x = r / R_virial;
        double f_c = std::log(1.0 + c_NFW) - c_NFW / (1.0 + c_NFW);
        
        // Avoid division by zero at r=0
        if (x < 1e-10) x = 1e-10;
        
        double Phi = -G * M_cluster * std::log(1.0 + c_NFW * x) / (x * f_c * R_virial);
        
        return Phi;
    }

    std::string toWolfram() const
    {
        return "VirgoClusterPotential[r_, Mcluster_: 1.2*^15 * 1.989*^30, Rvirial_: 2.2 * 3.086*^22, c_: 4] := "
               "Module[{x, fc, G}, "
               "G = 6.6743*^-11; "
               "x = Max[r / Rvirial, 10^-10]; "
               "fc = Log[1 + c] - c/(1 + c); "
               "-G * Mcluster * Log[1 + c*x] / (x * fc * Rvirial)]";
    }

    std::string getSignature() const { return "VirgoClusterGravitationalPotentialTerm(r, params)"; }
    std::string getCategory() const { return "gravity"; }
    std::string getName() const { return "VirgoClusterPotential"; }
    std::string getDescription() const { return "NFW gravitational potential: ╬ª(r) = -G┬╖M┬╖ln(1+c┬╖x)/(x┬╖f(c)┬╖R_vir) with c~4"; }
};

// ========================================
// CLASS 823: VirgoClusterDarkMatterTerm
// Category: dark_matter
// Physics: Dark matter halo profile (NFW) - ~85% of cluster mass
// ╧ü_DM(r) = ╧ü_s / ((r/r_s)(1 + r/r_s)┬▓)
// ========================================
class VirgoClusterDarkMatterTerm
{
private:
    double rho_s;          // Scale density (kg/m┬│, ~3e-23 for Virgo)
    double r_s;            // Scale radius (m, ~550 kpc = 1.70e22 m)
    double M_DM;           // Total dark matter mass (~1e15 M_sun)

public:
    VirgoClusterDarkMatterTerm(double scale_density = 3e-23, double scale_radius = 550 * KPC_TO_M, double dm_mass = 1e15 * M_SUN)
        : rho_s(scale_density), r_s(scale_radius), M_DM(dm_mass) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // NFW density profile: ╧ü(r) = ╧ü_s / ((r/r_s) * (1 + r/r_s)┬▓)
        double x = r / r_s;
        if (x < 1e-10) x = 1e-10;  // Avoid singularity at r=0
        
        double rho_DM = rho_s / (x * (1.0 + x) * (1.0 + x));
        
        return rho_DM;
    }

    std::string toWolfram() const
    {
        return "VirgoClusterDarkMatter[r_, rhoS_: 3*^-23, rS_: 550 * 3.086*^19] := "
               "Module[{x}, "
               "x = Max[r/rS, 10^-10]; "
               "rhoS / (x * (1 + x)^2)]";
    }

    std::string getSignature() const { return "VirgoClusterDarkMatterTerm(r, params)"; }
    std::string getCategory() const { return "dark_matter"; }
    std::string getName() const { return "VirgoClusterDarkMatter"; }
    std::string getDescription() const { return "NFW dark matter profile: ╧ü(r) = ╧ü_s/((r/r_s)┬╖(1+r/r_s)┬▓) with r_s~550 kpc"; }
};

// ========================================
// CLASS 824: VirgoClusterM87JetTerm
// Category: agn
// Physics: M87 (Virgo A) central AGN relativistic jet energy injection
// Jet power ~1e44 erg/s = 1e37 W, affects ICM heating and cluster evolution
// ========================================
class VirgoClusterM87JetTerm
{
private:
    double L_jet;          // Jet luminosity (W, ~1e37 for M87)
    double theta_jet;      // Jet opening angle (rad, ~0.1)
    double v_jet;          // Jet velocity (m/s, ~0.99c)
    double r_jet_base;     // Jet base radius (m, ~100 pc = 3.086e18 m)

public:
    VirgoClusterM87JetTerm(double luminosity = 1e37, double angle = 0.1, double velocity = 0.99 * C_LIGHT, double base_r = 100 * 3.086e16)
        : L_jet(luminosity), theta_jet(angle), v_jet(velocity), r_jet_base(base_r) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // Jet energy density: u_jet(r) = L_jet / (4╧Ç┬╖r┬▓┬╖v_jet┬╖╬⌐_jet)
        // where ╬⌐_jet = 2╧Ç(1 - cos(╬╕_jet)) is the solid angle
        double Omega_jet = 2.0 * M_PI * (1.0 - std::cos(theta_jet));
        
        // Energy density in jet cone
        double u_jet = L_jet / (4.0 * M_PI * r * r * v_jet * (Omega_jet / (4.0 * M_PI)));
        
        // Lorentz factor for relativistic correction
        double gamma = 1.0 / std::sqrt(1.0 - (v_jet * v_jet) / (C_LIGHT * C_LIGHT));
        
        return u_jet * gamma;
    }

    std::string toWolfram() const
    {
        return "VirgoM87Jet[r_, Ljet_: 10^37, theta_: 0.1, vjet_: 0.99 * 2.998*^8] := "
               "Module[{OmegaJet, uJet, gamma, c}, "
               "c = 2.998*^8; "
               "OmegaJet = 2 * Pi * (1 - Cos[theta]); "
               "uJet = Ljet / (4 * Pi * r^2 * vjet * (OmegaJet/(4*Pi))); "
               "gamma = 1/Sqrt[1 - (vjet/c)^2]; "
               "uJet * gamma]";
    }

    std::string getSignature() const { return "VirgoClusterM87JetTerm(r, params)"; }
    std::string getCategory() const { return "agn"; }
    std::string getName() const { return "VirgoM87Jet"; }
    std::string getDescription() const { return "M87 AGN jet energy density: u_jet = L_jet┬╖╬│/(4╧Ç┬╖r┬▓┬╖v_jet┬╖╬⌐) with L~10┬│Γü╖ W"; }
};

// ========================================
// CLASS 825: VirgoClusterTidalStrippingTerm
// Category: dynamics
// Physics: Tidal stripping of infalling galaxies by cluster potential
// Tidal radius r_t = r┬╖(M_gal/(3┬╖M_cluster(<r)))^(1/3)
// ========================================
class VirgoClusterTidalStrippingTerm
{
private:
    double M_cluster;      // Cluster mass (kg)
    double R_virial;       // Virial radius (m)
    double G;              // Gravitational constant

public:
    VirgoClusterTidalStrippingTerm(double mass = 1.2e15 * M_SUN, double r_vir = 2.2 * MPC_TO_M)
        : M_cluster(mass), R_virial(r_vir), G(G_CONST) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // Get galaxy mass from params, default to 1e11 M_sun (typical spiral)
        double M_gal = (params.count("M_gal") ? params.at("M_gal") : 1e11 * M_SUN);
        
        // Enclosed cluster mass at radius r (NFW-like approximation)
        double x = r / R_virial;
        double M_enclosed = M_cluster * (x * x * x) / (1.0 + x) / (1.0 + x);
        
        // Tidal radius: r_t = r * (M_gal / (3 * M_enclosed))^(1/3)
        double r_tidal = r * std::pow(M_gal / (3.0 * M_enclosed), 1.0 / 3.0);
        
        // Tidal acceleration at galaxy edge
        double a_tidal = 2.0 * G * M_enclosed * r_tidal / (r * r * r);
        
        return a_tidal;
    }

    std::string toWolfram() const
    {
        return "VirgoTidalStripping[r_, Mgal_, Mcluster_: 1.2*^15 * 1.989*^30, Rvirial_: 2.2 * 3.086*^22] := "
               "Module[{x, Menclosed, rtidal, G}, "
               "G = 6.6743*^-11; "
               "x = r / Rvirial; "
               "Menclosed = Mcluster * x^3 / (1 + x)^2; "
               "rtidal = r * (Mgal / (3 * Menclosed))^(1/3); "
               "2 * G * Menclosed * rtidal / r^3]";
    }

    std::string getSignature() const { return "VirgoClusterTidalStrippingTerm(r, params)"; }
    std::string getCategory() const { return "dynamics"; }
    std::string getName() const { return "VirgoTidalStripping"; }
    std::string getDescription() const { return "Tidal stripping acceleration: a_tidal = 2┬╖G┬╖M(<r)┬╖r_t/r┬│ with r_t from Jacobi radius"; }
};

// ========================================
// CLASS 826: VirgoClusterVirialTerm
// Category: thermodynamics
// Physics: Virial theorem equilibrium: 2K + U = 0, ╧â_v┬▓ = G┬╖M_vir/(3┬╖R_vir)
// Velocity dispersion ╧â_v ~ 700 km/s for Virgo
// ========================================
class VirgoClusterVirialTerm
{
private:
    double M_virial;       // Virial mass (kg)
    double R_virial;       // Virial radius (m)
    double sigma_v;        // Velocity dispersion (m/s, ~700 km/s)
    double G;              // Gravitational constant

public:
    VirgoClusterVirialTerm(double mass = 1.2e15 * M_SUN, double r_vir = 2.2 * MPC_TO_M, double sigma = 700e3)
        : M_virial(mass), R_virial(r_vir), sigma_v(sigma), G(G_CONST) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // r parameter included for API consistency (unused for virial equilibrium)
        (void)r;
        
        // Virial relation: ╧â┬▓ = G┬╖M_vir / (3┬╖R_vir)
        // Returns virial equilibrium parameter (should be ~1 if in equilibrium)
        double sigma_virial = std::sqrt(G * M_virial / (3.0 * R_virial));
        
        // Ratio of observed to virial velocity dispersion
        double virial_ratio = sigma_v / sigma_virial;
        
        return virial_ratio;
    }

    // Helper method (extends base interface for specific virial mass calculations)
    double computeVirialMass() const
    {
        // M_vir = 3┬╖╧â┬▓┬╖R_vir / G
        return 3.0 * sigma_v * sigma_v * R_virial / G;
    }

    std::string toWolfram() const
    {
        return "VirgoClusterVirial[Mvir_: 1.2*^15 * 1.989*^30, Rvir_: 2.2 * 3.086*^22, sigmaV_: 700000] := "
               "Module[{sigmaVirial, G}, "
               "G = 6.6743*^-11; "
               "sigmaVirial = Sqrt[G * Mvir / (3 * Rvir)]; "
               "sigmaV / sigmaVirial]";
    }

    std::string getSignature() const { return "VirgoClusterVirialTerm(r, params)"; }
    std::string getCategory() const { return "thermodynamics"; }
    std::string getName() const { return "VirgoClusterVirial"; }
    std::string getDescription() const { return "Virial equilibrium ratio: ╧â_obs/╧â_vir where ╧â_vir┬▓ = G┬╖M/(3┬╖R), ╧â~700 km/s"; }
};

// ========================================
// CLASS 827: VirgoClusterXRayLuminosityTerm
// Category: radiation
// Physics: X-ray emission from hot ICM: L_X Γê¥ n_e┬▓┬╖╬¢(T)┬╖V
// Virgo X-ray luminosity ~10^43 erg/s = 10^36 W
// ========================================
class VirgoClusterXRayLuminosityTerm
{
private:
    double n_e0;           // Central electron density (mΓü╗┬│)
    double T_ICM;          // ICM temperature (keV)
    double r_c;            // Core radius (m)
    double beta;           // Beta model parameter

public:
    VirgoClusterXRayLuminosityTerm(double ne_central = 3e3, double temp = 2.3, double core_r = 40 * KPC_TO_M, double beta_param = 0.5)
        : n_e0(ne_central), T_ICM(temp), r_c(core_r), beta(beta_param) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // Electron density profile (beta model)
        double n_e = n_e0 * std::pow(1.0 + (r * r) / (r_c * r_c), -1.5 * beta);
        
        // Cooling function ╬¢(T) ~ 3e-23 * sqrt(T_keV) [W┬╖m┬│] (simplified)
        double Lambda_T = 3e-23 * std::sqrt(T_ICM);
        
        // X-ray emissivity: ╬╡_X = n_e┬▓ * ╬¢(T) [W/m┬│]
        double epsilon_X = n_e * n_e * Lambda_T;
        
        return epsilon_X;
    }

    double computeTotalLuminosity(double R_max) const
    {
        // Integrate emissivity over volume (simplified spherical)
        // L_X = Γê½ ╬╡_X(r) 4╧Çr┬▓ dr
        // Using analytical beta-model integration approximation
        double L_X = 4.0 * M_PI * n_e0 * n_e0 * 3e-23 * std::sqrt(T_ICM) * 
                     std::pow(r_c, 3) * M_PI / (6.0 * beta);
        return L_X;
    }

    std::string toWolfram() const
    {
        return "VirgoXRay[r_, ne0_: 3000, T_: 2.3, rc_: 40 * 3.086*^19, beta_: 0.5] := "
               "Module[{ne, LambdaT}, "
               "ne = ne0 * (1 + (r/rc)^2)^(-1.5 * beta); "
               "LambdaT = 3*^-23 * Sqrt[T]; "
               "ne^2 * LambdaT]";
    }

    std::string getSignature() const { return "VirgoClusterXRayLuminosityTerm(r, params)"; }
    std::string getCategory() const { return "radiation"; }
    std::string getName() const { return "VirgoXRay"; }
    std::string getDescription() const { return "X-ray emissivity: ╬╡_X = n_e┬▓┬╖╬¢(T) with ╬¢(T)~3e-23┬╖ΓêÜT W┬╖m┬│"; }
};

// ========================================
// CLASS 828: VirgoClusterVelocityDispersionTerm
// Category: kinematics
// Physics: Galaxy velocity dispersion profile ╧â(r)
// Central ╧â ~ 700 km/s, decreases with radius
// ========================================
class VirgoClusterVelocityDispersionTerm
{
private:
    double sigma_0;        // Central velocity dispersion (m/s, ~700 km/s)
    double r_sigma;        // Scale radius for dispersion (m, ~500 kpc)
    double G;              // Gravitational constant
    double M_cluster;      // Cluster mass

public:
    VirgoClusterVelocityDispersionTerm(double sigma_central = 700e3, double scale_r = 500 * KPC_TO_M, double mass = 1.2e15 * M_SUN)
        : sigma_0(sigma_central), r_sigma(scale_r), G(G_CONST), M_cluster(mass) {}

    double compute(double r, const std::map<std::string, double>& params) const
    {
        // Velocity dispersion profile: ╧â(r) = ╧â_0 / sqrt(1 + (r/r_╧â)┬▓)
        double sigma_r = sigma_0 / std::sqrt(1.0 + (r * r) / (r_sigma * r_sigma));
        
        return sigma_r;
    }

    double computeDynamicalMass(double r) const
    {
        // M(<r) = r┬╖╧â┬▓(r) / G (simplified spherical Jeans equation)
        double sigma_r = sigma_0 / std::sqrt(1.0 + (r * r) / (r_sigma * r_sigma));
        return r * sigma_r * sigma_r / G;
    }

    std::string toWolfram() const
    {
        return "VirgoVelocityDispersion[r_, sigma0_: 700000, rSigma_: 500 * 3.086*^19] := "
               "sigma0 / Sqrt[1 + (r/rSigma)^2]";
    }

    std::string getSignature() const { return "VirgoClusterVelocityDispersionTerm(r, params)"; }
    std::string getCategory() const { return "kinematics"; }
    std::string getName() const { return "VirgoVelocityDispersion"; }
    std::string getDescription() const { return "Velocity dispersion profile: ╧â(r) = ╧â_0/ΓêÜ(1+(r/r_╧â)┬▓) with ╧â_0~700 km/s"; }
};

// ========================================
// CLASS 829: SMBHMSigmaRelationTerm
// Category: scaling_relations
// Physics: M-╧â relation from source82 SMBH module
// M_BH = 1.9e8 * (╧â/200 km/s)^4.38 M_sun (McConnell & Ma 2013)
// Links SMBH mass to host galaxy velocity dispersion
// ========================================
