#ifndef CELESTIAL_BODY_H
#define CELESTIAL_BODY_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <functional>

// ============================================================================
// Structure Definitions
// ============================================================================

struct ResonanceParams
{
    double fDPM = 1e12;
    double fTHz = 1e12;
    double Evac_neb = 7.09e-36;
    double Evac_ISM = 7.09e-37;
    double Delta_Evac = 6.381e-36;
    double Fsuper = 6.287e-19;
    double UA_SCM = 10;
    double omega_i = 1e-8;
    double k4_res = 1.0;
    double freact = 1e10;
    double fquantum = 1.445e-17;
    double fAether = 1.576e-35;
    double fosc = 4.57e14;
    double fTRZ = 0.1;
    double c_res = 3e8;
};

struct MUGESystem
{
    std::string name;
    double I;
    double A;
    double omega1;
    double omega2;
    double Vsys;
    double vexp;
    double t;
    double z;
    double ffluid;
    double M;
    double r;
    double B;
    double Bcrit;
    double rho_fluid;
    double g_local;
    double M_DM;
    double delta_rho_rho;
};

struct CelestialBody
{
    std::string name;
    double Ms;          // Mass (kg)
    double Rs;          // Radius (m)
    double Rb;          // Bubble radius (e.g., heliosphere or magnetosphere, m)
    double Ts_surface;  // Surface temperature (K)
    double omega_s;     // Rotation rate (rad/s)
    double Bs_avg;      // Average surface magnetic field (T)
    double SCm_density; // SCm density (kg/m^3)
    double QUA;         // Trapped Universal Aether charge (C)
    double Pcore;       // Planetary core penetration factor
    double PSCm;        // SCm penetration factor
    double omega_c;     // Cycle frequency (rad/s)
};

// ============================================================================
// UQFF 2.0-Enhanced Framework: Self-Expanding Physics Module
// ============================================================================

// Base class for extensible physics terms
class PhysicsTerm
{
public:
    virtual ~PhysicsTerm() = default;
    virtual double compute(const std::map<std::string, double> &params) const = 0;
    virtual std::string description() const = 0;
    virtual std::string version() const { return "1.0"; }
};

// Dark matter halo contribution term
class DarkMatterHaloTerm : public PhysicsTerm
{
private:
    double M_halo;  // Halo mass (kg)
    double r_scale; // Scale radius (m)

public:
    DarkMatterHaloTerm(double mass, double scale);
    double compute(const std::map<std::string, double> &params) const override;
    std::string description() const override;
};

// Vacuum energy fluctuation term
class VacuumEnergyTerm : public PhysicsTerm
{
private:
    double E_vac_scale; // Vacuum energy scale
    double lambda;      // Coupling strength

public:
    VacuumEnergyTerm(double e_scale, double coupling);
    double compute(const std::map<std::string, double> &params) const override;
    std::string description() const override;
};

// ============================================================================
// UQFFModule5: Self-Expanding Unified Field Module
// ============================================================================
class UQFFModule5
{
private:
    std::vector<std::unique_ptr<PhysicsTerm>> dynamic_terms;
    std::map<std::string, double> dynamic_parameters;
    std::map<std::string, std::string> metadata;
    double learning_rate;
    bool logging_enabled;
    void log(const std::string &message) const;

public:
    UQFFModule5();
    void registerDynamicTerm(std::unique_ptr<PhysicsTerm> term);
    void setDynamicParameter(const std::string &name, double value);
    double getDynamicParameter(const std::string &name, double default_val = 0.0) const;
    double computeDynamicContributions(const std::map<std::string, double> &params) const;
    void exportState(const std::string &filename) const;
    void setLearningRate(double rate);
    void setEnableLogging(bool enable);
    void printInfo() const;

    // Enhanced compute functions
    double compute_Ug1_enhanced(const CelestialBody &body, double r, double t, double tn,
                                double alpha, double delta_def, double k1);
    double compute_Ug2_enhanced(const CelestialBody &body, double r, double t, double tn,
                                double k2, double QA, double delta_sw, double v_sw,
                                double HSCm, double rho_A, double kappa);
    double compute_Ug3_enhanced(const CelestialBody &body, double r, double t, double tn,
                                double theta, double rho_A, double kappa, double k3);
    double compute_MUGE_enhanced(const MUGESystem &sys, const ResonanceParams &res, bool use_compressed = true);
};

// ============================================================================
// Function Declarations
// ============================================================================

// CelestialBody functions
double compute_Ug1(const CelestialBody &body, double r, double t, double tn, double alpha, double delta_def, double k1);
double compute_Ug2(const CelestialBody &body, double r, double t, double tn, double k2, double QA, double delta_sw, double v_sw, double HSCm, double rho_A, double kappa);
double compute_Ug3(const CelestialBody &body, double r, double t, double tn, double theta, double rho_A, double kappa, double k3);
double compute_Um(const CelestialBody &body, double t, double tn, double rj, double gamma, double rho_A, double kappa, double num_strings, double phi_hat = 1.0);
void output_json_params(const CelestialBody &body);
std::vector<CelestialBody> load_bodies(const std::string &filename);

// MUGE functions
double compute_compressed_MUGE(const MUGESystem &sys);
double compute_resonance_MUGE(const MUGESystem &sys, const ResonanceParams &res);
std::vector<MUGESystem> load_muge_systems(const std::string &filename);

// Helper functions
double step_function(double r, double Rb);
double compute_Ereact(double t, double rho_SCm, double v_SCm, double rho_A, double kappa);
double compute_mu_s(double t, double Bs, double omega_c, double Rs, double SCm_contrib = 1e3);
double compute_grad_Ms_r(double Ms, double Rs);
double compute_Bj(double t, double omega_c, double SCm_contrib = 1e3);
double compute_omega_s_t(double t, double omega_s, double omega_c);
double compute_mu_j(double t, double omega_c, double Rs, double SCm_contrib = 1e3);

#endif // CELESTIAL_BODY_H
