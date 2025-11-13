#ifndef MUGE_H
#define MUGE_H

#include <string>
#include <vector>

// Already defined in CelestialBody.h but included here for standalone use
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

// Compressed MUGE term functions
double compute_compressed_base(const MUGESystem &sys);
double compute_compressed_expansion(const MUGESystem &sys, double H0 = 2.269e-18);
double compute_compressed_super_adj(const MUGESystem &sys);
double compute_compressed_env();
double compute_compressed_Ug_sum();
double compute_compressed_cosm(double Lambda = 1.1e-52);
double compute_compressed_quantum(double hbar = 1.0546e-34, double Delta_x_p = 1e-68, double integral_psi = 2.176e-18, double tHubble = 4.35e17);
double compute_compressed_fluid(const MUGESystem &sys);
double compute_compressed_perturbation(const MUGESystem &sys);

// Resonance MUGE term functions
double compute_aDPM(const MUGESystem &sys, const ResonanceParams &res);
double compute_aTHz(double aDPM, const MUGESystem &sys, const ResonanceParams &res);
double compute_avac_diff(double aDPM, const MUGESystem &sys, const ResonanceParams &res);
double compute_asuper_freq(double aDPM, const ResonanceParams &res);
double compute_aaether_res(double aDPM, const ResonanceParams &res);
double compute_Ug4i(double aDPM, const MUGESystem &sys, const ResonanceParams &res);
double compute_aquantum_freq(double aDPM, const ResonanceParams &res);
double compute_aAether_freq(double aDPM, const ResonanceParams &res);
double compute_afluid_freq(const MUGESystem &sys, const ResonanceParams &res);
double compute_Osc_term();
double compute_aexp_freq(double aDPM, const MUGESystem &sys, const ResonanceParams &res, double H_z = 2.270e-18);
double compute_fTRZ(const ResonanceParams &res);
double compute_a_wormhole(double r, double b = 1.0, double f_worm = 1.0, double Evac_neb = 7.09e-36);

// Full MUGE functions
double compute_compressed_MUGE(const MUGESystem &sys);
double compute_resonance_MUGE(const MUGESystem &sys, const ResonanceParams &res);

std::vector<MUGESystem> load_muge_systems(const std::string &filename);

#endif // MUGE_H
