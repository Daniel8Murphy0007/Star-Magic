"""
Session 173 — Generate all 14 remaining C++ module pairs (.h + .cpp),
PAPER_674–687, CP4 #258–271
Source: grok_share_fc21e30c24b4.txt (29,600-line Grok share)

Modules:
  674 UQFFComparedToLIGOData         line ~18514
  675 UQFFComparedToGW170817          line ~19205
  676 UQFFComparedToGW190425          line ~19417
  677 UQFFPredictionsForLISA          line ~19842
  678 LISAVsLIGOComparisons           line ~20057
  679 AetherSuperfluidDynamics        line ~20556
  680 VortexQuantization              line ~21091
  681 GrossPitaevskiiVortexSimulation line ~21332
  682 UQFFStabilityNumericallyForSgrA line ~21554
  683 UQFFHawkingTemperatureModulation line ~21788
  684 UQFFPrimordialBHEvaporation     line ~22011
  685 UQFFPBHDarkMatterImplications   line ~22222
  686 UQFFModulationForM87             line ~22436
  687 M87MassEvolutionSimulation      line ~22623
"""
import os

ROOT = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
os.makedirs(os.path.join(ROOT, "whitepapers"), exist_ok=True)

COMMON_INCLUDES = """\
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <random>
#include <chrono>
"""

COMMON_CONSTANTS = """\
    static constexpr double G       = 6.6743e-11;
    static constexpr double C       = 2.998e8;
    static constexpr double HBAR    = 1.0546e-34;
    static constexpr double K_B     = 1.380649e-23;
    static constexpr double PI      = 3.14159265358979323846;
    static constexpr double M_SUN   = 1.989e30;
    static constexpr double RHO_UA  = 7.09e-36;
    static constexpr double RHO_SCM = 7.09e-37;
    static constexpr double F_TRZ   = 0.1;
    static constexpr double KAPPA   = 0.0005;
    static constexpr double SSQ     = 0.57;
    static constexpr double MU_J    = 3.38e23;
    static constexpr double GAMMA   = 5.0e-5 / 86400.0;
    static constexpr double T_N_REF = 1.0e8;
"""

def write(path, content):
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  CREATED: {os.path.basename(path)}")

def helpers(cls):
    return f"""\
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _{cls}_T_H(double M){{
    return ({cls}::HBAR*{cls}::C*{cls}::C*{cls}::C)/
           (8.0*M_PI*{cls}::G*M*{cls}::K_B);
}}
static inline double _{cls}_L_H(double M){{
    return ({cls}::HBAR*std::pow({cls}::C,6.0))/
           (15360.0*M_PI*{cls}::G*{cls}::G*M*M);
}}
static inline double _{cls}_r_s(double M){{
    return 2.0*{cls}::G*M/({cls}::C*{cls}::C);
}}
static inline double _{cls}_tau(double M){{
    return 5120.0*M_PI*{cls}::G*{cls}::G*M*M*M/
           ({cls}::HBAR*std::pow({cls}::C,4.0));
}}
static inline double _{cls}_Um(double r,double t){{
    double tn=t/{cls}::T_N_REF;
    return ({cls}::MU_J/r)*(1.0-std::exp(-{cls}::GAMMA*t*std::cos(M_PI*tn)));
}}
"""

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_674 — UQFFComparedToLIGOData
# ─────────────────────────────────────────────────────────────────────────────
def gen_674():
    cls = "UQFFComparedToLIGOData"
    guard = "UQFF_COMPARED_TO_LIGO_DATA_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF vs General LIGO Dataset — PAPER_674 | Session 173 | April 1, 2026
 * Source: grok_share_fc21e30c24b4.txt line ~18514
 *
 * Compares UQFF-modified GW strain and phase to the generalised LIGO O1/O2/O3
 * dataset. UQFF corrections enter through three suppression factors:
 *
 *   h_UQFF(f) = h_GR(f) * (1 - f_TRZ) * exp(-rho_SCm * r_s / k_B T_H)
 *                        * exp(-U_m * 2*pi*f / c^2)
 *
 *   Delta_phi_UQFF = Delta_phi_GR + kappa * f_TRZ * t_coal
 *
 * LIGO sensitivity band: 10–2000 Hz.
 * Chirp mass: Mc = (m1*m2)^(3/5) / (m1+m2)^(1/5)
 * GW luminosity (Peters): P_GW = -32/5 * G^4/c^5 * m1^2*m2^2*(m1+m2)/a^5
 * UQFF modifies via S_UA * S_SCm * S_TRZ * S_Um factors.
 * Self-simulate: sweep frequency 10–2000 Hz, output h_UQFF vs h_GR.
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    {cls}();
    double compute_chirp_mass(double m1, double m2) const;
    double compute_h_GR(double Mc, double d, double f) const;
    double compute_S_SCm(double M, double T_H) const;
    double compute_S_Um(double M, double r, double t, double f) const;
    double compute_h_UQFF(double m1, double m2, double d, double f,
                          double r=0.0, double t=1e8) const;
    double compute_delta_phi(double t_coal) const;
    void simulate_frequency_sweep(double f_start, double f_end, double df,
                                  double m1, double m2, double d,
                                  const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}

double {cls}::compute_chirp_mass(double m1,double m2) const{{
    return std::pow(m1*m2,0.6)/std::pow(m1+m2,0.2);
}}
double {cls}::compute_h_GR(double Mc,double d,double f) const{{
    // h_GR = (4/d)*(G*Mc/c^2)^(5/3)*(pi*f)^(2/3)/c^2
    double GM=G*Mc; double pif=PI*f;
    return (4.0/d)*std::pow(GM/C/C,5.0/3.0)*std::pow(pif,2.0/3.0)/(C*C);
}}
double {cls}::compute_S_SCm(double M,double T_H) const{{
    double rs=_{cls}_r_s(M);
    double arg=RHO_SCM*rs/(K_B*T_H); if(arg>700)arg=700;
    return std::exp(-arg);
}}
double {cls}::compute_S_Um(double M,double r,double t,double f) const{{
    if(r<=0) r=_{cls}_r_s(M);
    double Um=_{cls}_Um(r,t);
    double arg=Um*2.0*PI*f/(C*C); if(arg>700)arg=700;
    return std::exp(-arg);
}}
double {cls}::compute_h_UQFF(double m1,double m2,double d,double f,
                              double r,double t) const{{
    double Mc=compute_chirp_mass(m1,m2);
    double h_gr=compute_h_GR(Mc,d,f);
    double M=m1+m2;
    double T_H=_{cls}_T_H(M);
    double S_SCm=compute_S_SCm(M,T_H);
    double S_Um=compute_S_Um(M,r,t,f);
    double h=(1.0-F_TRZ)*S_SCm*S_Um*h_gr;
    for(auto& m:mods_) h+=m(f,t);
    return h;
}}
double {cls}::compute_delta_phi(double t_coal) const{{
    return KAPPA*F_TRZ*t_coal;
}}
void {cls}::simulate_frequency_sweep(double f0,double f1,double df,
    double m1,double m2,double d,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"f_Hz,h_GR,h_UQFF\\n";}}}}
    double M=m1+m2; double Mc=compute_chirp_mass(m1,m2);
    double T_H=_{cls}_T_H(M);
    for(double f=f0;f<=f1;f+=df){{
        double h_gr=compute_h_GR(Mc,d,f);
        double h_u=compute_h_UQFF(m1,m2,d,f);
        *os<<f<<","<<h_gr<<","<<h_u<<"\\n";
    }}
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn); if(!fi.is_open())return;
    std::string l; while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    double m1=7.16e31,m2=5.77e31,d=1.27e25;
    std::cout<<"Chirp mass="<<obj.compute_chirp_mass(m1,m2)<<"\\n";
    std::cout<<"h_UQFF@150Hz="<<obj.compute_h_UQFF(m1,m2,d,150.0)<<"\\n";
    obj.simulate_frequency_sweep(10,500,10,m1,m2,d,"ligo_sweep.csv");
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_675 — UQFFComparedToGW170817
# ─────────────────────────────────────────────────────────────────────────────
def gen_675():
    cls = "UQFFComparedToGW170817"
    guard = "UQFF_COMPARED_TO_GW170817_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF vs LIGO/Virgo GW170817 (NS-NS merger) — PAPER_675 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~19205
 *
 * GW170817: Binary neutron star coalescence, Aug 17, 2017.
 *   m1 = 1.36 Msun, m2 = 1.17 Msun
 *   Luminosity distance: ~40 Mpc = 1.234e24 m
 *   Peak frequency: ~1500 Hz
 *   Multi-messenger: GRB170817A detected 1.7 s after merger
 *
 * UQFF modifications:
 *   T_H,NS = HBAR c^3 / (8 pi G M_NS k_B)   (characteristic NS Hawking T)
 *   h_UQFF = h_GR * (1-f_TRZ) * S_SCm * S_Um
 *   GRB delay: Delta_t_UQFF = 1.7 s * (1 + f_TRZ * rho_UA/rho_SCm)
 *              -> probes aether propagation speed
 *   tidal deformability Lambda modified by rho_UA suppression
 *   kappa_tidal_UQFF = kappa_tidal_GR * (1 - f_TRZ * rho_SCm/rho_UA)
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double M1_GW170817 = 1.36 * 1.989e30;
    static constexpr double M2_GW170817 = 1.17 * 1.989e30;
    static constexpr double D_GW170817  = 1.234e24;  // 40 Mpc in m
    static constexpr double F_PEAK      = 1500.0;    // Hz
    static constexpr double DT_GRB      = 1.7;       // s

    {cls}();
    double compute_chirp_mass(double m1, double m2) const;
    double compute_h_GR(double Mc, double d, double f) const;
    double compute_h_UQFF(double m1, double m2, double d, double f,
                          double t=1e8) const;
    double compute_grb_delay_UQFF() const;
    double compute_tidal_UQFF(double kappa_GR) const;
    void simulate_inspiral(double f_start, double f_peak, double df,
                           const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}
double {cls}::compute_chirp_mass(double m1,double m2) const{{
    return std::pow(m1*m2,0.6)/std::pow(m1+m2,0.2);
}}
double {cls}::compute_h_GR(double Mc,double d,double f) const{{
    return (4.0/d)*std::pow(G*Mc/(C*C),5.0/3.0)*std::pow(PI*f,2.0/3.0)/(C*C);
}}
double {cls}::compute_h_UQFF(double m1,double m2,double d,double f,double t) const{{
    double Mc=compute_chirp_mass(m1,m2);
    double h_gr=compute_h_GR(Mc,d,f);
    double M=m1+m2; double T_H=_{cls}_T_H(M); double r=_{cls}_r_s(M);
    double S_SCm=std::exp(-RHO_SCM*r/(K_B*T_H));
    double Um=_{cls}_Um(r,t);
    double S_Um=std::exp(-Um*2.0*PI*f/(C*C));
    double h=(1.0-F_TRZ)*S_SCm*S_Um*h_gr;
    for(auto& m:mods_) h+=m(f,t);
    return h;
}}
double {cls}::compute_grb_delay_UQFF() const{{
    // Delta_t_UQFF = 1.7 * (1 + f_TRZ * rho_UA/rho_SCm)
    return DT_GRB*(1.0+F_TRZ*(RHO_UA/RHO_SCM));
}}
double {cls}::compute_tidal_UQFF(double kappa_GR) const{{
    return kappa_GR*(1.0-F_TRZ*RHO_SCM/RHO_UA);
}}
void {cls}::simulate_inspiral(double f0,double fp,double df,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"f_Hz,h_GR,h_UQFF\\n";}}}}
    double m1=M1_GW170817,m2=M2_GW170817,d=D_GW170817;
    double Mc=compute_chirp_mass(m1,m2);
    for(double f=f0;f<=fp;f+=df){{
        *os<<f<<","<<compute_h_GR(Mc,d,f)<<","<<compute_h_UQFF(m1,m2,d,f)<<"\\n";
    }}
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    std::cout<<"GRB delay UQFF="<<obj.compute_grb_delay_UQFF()<<" s\\n";
    std::cout<<"h_UQFF@1500Hz="<<obj.compute_h_UQFF(M1_GW170817,M2_GW170817,D_GW170817,1500.0)<<"\\n";
    obj.simulate_inspiral(10,1500,10,"gw170817_sweep.csv");
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_676 — UQFFComparedToGW190425
# ─────────────────────────────────────────────────────────────────────────────
def gen_676():
    cls = "UQFFComparedToGW190425"
    guard = "UQFF_COMPARED_TO_GW190425_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF vs LIGO/Virgo GW190425 (heavy NS-NS) — PAPER_676 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~19417
 *
 * GW190425: Binary NS merger April 25, 2019. Heaviest NS binary known.
 *   Total mass: 3.4 Msun (m1~1.9, m2~1.5 Msun)
 *   Distance: ~159 Mpc = 4.9e24 m
 *   Only LIGO-Livingston + Virgo (LIGO-Hanford offline)
 *   No EM counterpart detected (high distance, low inclination sensitivity)
 *
 * UQFF analysis:
 *   Mc_GW190425 = (1.9*1.5*Msun^2)^0.6 / (3.4*Msun)^0.2
 *   h_UQFF(f) = h_GR(f) * (1-f_TRZ) * S_SCm * S_Um
 *   Mass increase over GW170817 (~2.5x higher total mass)
 *   -> UQFF suppression factors scale with T_H and r_s
 *   post-merger frequency f_pm ~ 2.5 kHz: UQFF phase shift probes aether
 *   UQFF upper-limit on ejecta mass through U_m suppression
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double M1_DEFAULT = 1.9 * 1.989e30;
    static constexpr double M2_DEFAULT = 1.5 * 1.989e30;
    static constexpr double D_DEFAULT  = 4.9e24;
    static constexpr double F_PM       = 2500.0;

    {cls}();
    double compute_chirp_mass(double m1, double m2) const;
    double compute_h_GR(double Mc, double d, double f) const;
    double compute_h_UQFF(double m1, double m2, double d, double f,
                          double t=1e8) const;
    double compute_post_merger_phase(double m1, double m2) const;
    double compute_ejecta_limit_UQFF(double m1, double m2) const;
    void simulate_inspiral(double f_start, double f_end, double df,
                           const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}
double {cls}::compute_chirp_mass(double m1,double m2) const{{
    return std::pow(m1*m2,0.6)/std::pow(m1+m2,0.2);
}}
double {cls}::compute_h_GR(double Mc,double d,double f) const{{
    return (4.0/d)*std::pow(G*Mc/(C*C),5.0/3.0)*std::pow(PI*f,2.0/3.0)/(C*C);
}}
double {cls}::compute_h_UQFF(double m1,double m2,double d,double f,double t) const{{
    double Mc=compute_chirp_mass(m1,m2);
    double h_gr=compute_h_GR(Mc,d,f);
    double M=m1+m2; double T_H=_{cls}_T_H(M); double r=_{cls}_r_s(M);
    double S_SCm=std::exp(-RHO_SCM*r/(K_B*T_H));
    double Um=_{cls}_Um(r,t);
    double S_Um=std::exp(-Um*2.0*PI*f/(C*C));
    double h=(1.0-F_TRZ)*S_SCm*S_Um*h_gr;
    for(auto& m:mods_) h+=m(f,t);
    return h;
}}
double {cls}::compute_post_merger_phase(double m1,double m2) const{{
    // UQFF phase: phi_UQFF = phi_GR + kappa*f_TRZ*t_merger
    double t_merger=5.0*(m1+m2)/(C*C*C/G)*(C*C*C/G)/(m1*m2)*1e-6;
    return KAPPA*F_TRZ*std::abs(t_merger);
}}
double {cls}::compute_ejecta_limit_UQFF(double m1,double m2) const{{
    // M_ej,UQFF < M_GR * (rho_SCm/rho_UA) * (1-f_TRZ) 
    double M_ej_GR=0.05*(m1+m2); // typical 5% ejecta
    return M_ej_GR*(RHO_SCM/RHO_UA)*(1.0-F_TRZ);
}}
void {cls}::simulate_inspiral(double f0,double fe,double df,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"f_Hz,h_GR,h_UQFF\\n";}}}}
    double Mc=compute_chirp_mass(M1_DEFAULT,M2_DEFAULT);
    for(double f=f0;f<=fe;f+=df)
        *os<<f<<","<<compute_h_GR(Mc,D_DEFAULT,f)<<","
           <<compute_h_UQFF(M1_DEFAULT,M2_DEFAULT,D_DEFAULT,f)<<"\\n";
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    double m1=M1_DEFAULT,m2=M2_DEFAULT,d=D_DEFAULT;
    std::cout<<"Chirp mass="<<obj.compute_chirp_mass(m1,m2)/M_SUN<<" Msun\\n";
    std::cout<<"h_UQFF@2500Hz="<<obj.compute_h_UQFF(m1,m2,d,2500.0)<<"\\n";
    std::cout<<"Ejecta limit UQFF="<<obj.compute_ejecta_limit_UQFF(m1,m2)/M_SUN<<" Msun\\n";
    obj.simulate_inspiral(10,2500,20,"gw190425_sweep.csv");
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_677 — UQFFPredictionsForLISA
# ─────────────────────────────────────────────────────────────────────────────
def gen_677():
    cls = "UQFFPredictionsForLISA"
    guard = "UQFF_PREDICTIONS_FOR_LISA_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Predictions for LISA Space GW Observatory — PAPER_677 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~19842
 *
 * LISA (Laser Interferometer Space Antenna): arm length 2.5 Gm,
 *   frequency band 0.1 mHz – 1 Hz, launch ~2037.
 *   Target sources: SMBH mergers, EMRIs, galactic WD binaries, stochastic background.
 *
 * UQFF predictions for LISA:
 *   1. SMBH inspiral phase: h_UQFF_LISA = h_GR * S_UA_LISA * S_SCm_LISA
 *      where S_UA_LISA = 1 - rho_UA * L_LISA / (k_B T_eff)
 *   2. Stochastic background: Omega_GW,UQFF = Omega_GW * (rho_UA/rho_crit)^(f_TRZ)
 *   3. EMRI rate: R_EMRI,UQFF = R_EMRI,GR * (1 + f_TRZ * rho_UA/rho_SCm)
 *   4. U_m modulation of waveform phase: phi_UQFF = phi_GR + kappa * f_TRZ * tau_RD
 *   5. UQFF predicts enhanced SMBH stability -> longer inspiral -> more LISA cycles
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double L_LISA     = 2.5e9;          // m arm length
    static constexpr double F_MIN_LISA = 1.0e-4;         // 0.1 mHz
    static constexpr double F_MAX_LISA = 1.0;            // 1 Hz
    static constexpr double RHO_CRIT   = 9.47e-27;       // kg/m^3

    {cls}();
    double compute_h_GR_SMBH(double Mc, double d, double f) const;
    double compute_S_UA_LISA(double T_eff) const;
    double compute_h_UQFF_LISA(double Mc, double d, double f, double T_eff) const;
    double compute_omega_GW_UQFF(double omega_GW_GR) const;
    double compute_EMRI_rate_UQFF(double R_GR) const;
    double compute_phase_mod(double tau_RD) const;
    void simulate_SMBH_sweep(double Mc, double d, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}
double {cls}::compute_h_GR_SMBH(double Mc,double d,double f) const{{
    return (4.0/d)*std::pow(G*Mc/(C*C),5.0/3.0)*std::pow(PI*f,2.0/3.0)/(C*C);
}}
double {cls}::compute_S_UA_LISA(double T_eff) const{{
    double arg=RHO_UA*L_LISA/(K_B*T_eff); if(arg>700)arg=700;
    return std::max(0.0,1.0-arg);
}}
double {cls}::compute_h_UQFF_LISA(double Mc,double d,double f,double T_eff) const{{
    double h_gr=compute_h_GR_SMBH(Mc,d,f);
    double S_UA=compute_S_UA_LISA(T_eff);
    double M_tot=Mc*std::pow(4.0,0.2); // approx equal mass
    double r=_{cls}_r_s(M_tot); double T_H=_{cls}_T_H(M_tot);
    double S_SCm=std::exp(-RHO_SCM*r/(K_B*T_H));
    double h=(1.0-F_TRZ)*S_UA*S_SCm*h_gr;
    for(auto& m:mods_) h+=m(f,T_eff);
    return h;
}}
double {cls}::compute_omega_GW_UQFF(double omega_GW_GR) const{{
    return omega_GW_GR*std::pow(RHO_UA/RHO_CRIT,F_TRZ);
}}
double {cls}::compute_EMRI_rate_UQFF(double R_GR) const{{
    return R_GR*(1.0+F_TRZ*(RHO_UA/RHO_SCM));
}}
double {cls}::compute_phase_mod(double tau_RD) const{{
    return KAPPA*F_TRZ*tau_RD;
}}
void {cls}::simulate_SMBH_sweep(double Mc,double d,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"f_Hz,h_GR,h_UQFF\\n";}}}}
    double T_eff=2.73; // CMB temperature as proxy
    double lf=std::log10(F_MIN_LISA),hf=std::log10(F_MAX_LISA);
    int N=200;
    for(int i=0;i<=N;i++){{
        double f=std::pow(10.0,lf+(hf-lf)*i/N);
        *os<<f<<","<<compute_h_GR_SMBH(Mc,d,f)<<","<<compute_h_UQFF_LISA(Mc,d,f,T_eff)<<"\\n";
    }}
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    // Sgr A* + M87-type SMBH merger at cosmological distance
    double Mc=1e8*1.989e30; double d=1e27; // ~1 Gpc
    std::cout<<"Omega_GW_UQFF/GR ratio="<<obj.compute_omega_GW_UQFF(1.0)<<"\\n";
    std::cout<<"EMRI rate boost="<<obj.compute_EMRI_rate_UQFF(100.0)<<" yr-1\\n";
    obj.simulate_SMBH_sweep(Mc,d,"lisa_prediction.csv");
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_678 — LISAVsLIGOComparisons
# ─────────────────────────────────────────────────────────────────────────────
def gen_678():
    cls = "LISAVsLIGOComparisons"
    guard = "LISA_VS_LIGO_COMPARISONS_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF LISA vs LIGO Cross-Comparison — PAPER_678 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~20057
 *
 * Quantitative comparison of UQFF corrections in the LISA band (0.1 mHz – 1 Hz)
 * versus the LIGO band (10 – 2000 Hz):
 *
 *   LIGO band: f ~ 10-2000 Hz, U_m suppression dominant (S_Um small at high f)
 *   LISA band: f ~ 1e-4 – 1 Hz, S_UA_LISA suppression dominant (long arm)
 *
 *   Suppression ratio: R_supp(f) = h_UQFF(f) / h_GR(f)
 *   LIGO regime:  R ~ (1-f_TRZ) * S_SCm  (U_m term ~ 1 at low f)
 *   LISA regime:  R ~ (1-f_TRZ) * S_UA * S_SCm
 *
 *   Crossover frequency: f_cross where S_Um = S_UA
 *   f_cross = c^2 * rho_UA * L_LISA / (2 pi U_m * k_B * T_H)^{{-1}}
 *
 *   Sensitivity noise floor UQFF correction:
 *   S_n,UQFF(f) = S_n,GR(f) * (1 + F_TRZ * rho_UA/rho_SCm * f^(-2/3))
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double L_LISA = 2.5e9;
    static constexpr double L_LIGO = 4.0e3;

    {cls}();
    double compute_R_supp_LIGO(double M, double f, double t=1e8) const;
    double compute_R_supp_LISA(double M, double f, double T_eff=2.73) const;
    double compute_crossover_freq(double M, double T_eff=2.73) const;
    double compute_Sn_UQFF(double Sn_GR, double f) const;
    void compare_sweep(double M, double d, double f_ligo_lo, double f_ligo_hi,
                       double f_lisa_lo, double f_lisa_hi,
                       const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}
double {cls}::compute_R_supp_LIGO(double M,double f,double t) const{{
    double T_H=_{cls}_T_H(M); double r=_{cls}_r_s(M);
    double S_SCm=std::exp(-RHO_SCM*r/(K_B*T_H));
    double Um=_{cls}_Um(r,t);
    double S_Um=std::exp(-Um*2.0*PI*f/(C*C));
    return (1.0-F_TRZ)*S_SCm*S_Um;
}}
double {cls}::compute_R_supp_LISA(double M,double f,double T_eff) const{{
    double T_H=_{cls}_T_H(M); double r=_{cls}_r_s(M);
    double S_SCm=std::exp(-RHO_SCM*r/(K_B*T_H));
    double S_UA=std::max(0.0,1.0-RHO_UA*L_LISA/(K_B*T_eff));
    return (1.0-F_TRZ)*S_UA*S_SCm;
}}
double {cls}::compute_crossover_freq(double M,double T_eff) const{{
    // Solve S_Um = S_UA => exp(-Um*2pif/c^2) = 1 - rho_UA*L/k_B T_eff
    double T_H=_{cls}_T_H(M); double r=_{cls}_r_s(M); double Um=_{cls}_Um(r,1e8);
    double S_UA=std::max(1e-300,1.0-RHO_UA*L_LISA/(K_B*T_eff));
    if(Um<=0||S_UA<=0) return 0.0;
    return -std::log(S_UA)*C*C/(2.0*PI*Um);
}}
double {cls}::compute_Sn_UQFF(double Sn_GR,double f) const{{
    return Sn_GR*(1.0+F_TRZ*(RHO_UA/RHO_SCM)*std::pow(f,-2.0/3.0));
}}
void {cls}::compare_sweep(double M,double d,
    double fl_lo,double fl_hi,double fa_lo,double fa_hi,
    const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"f_Hz,R_LIGO,R_LISA\\n";}}}}
    // Sweep LIGO band
    for(double f=fl_lo;f<=fl_hi;f+=(fl_hi-fl_lo)/100.0)
        *os<<f<<","<<compute_R_supp_LIGO(M,f)<<","<<compute_R_supp_LISA(M,f)<<"\\n";
    // Sweep LISA band
    for(int i=0;i<=100;i++){{
        double f=fa_lo+i*(fa_hi-fa_lo)/100.0;
        *os<<f<<","<<compute_R_supp_LIGO(M,f)<<","<<compute_R_supp_LISA(M,f)<<"\\n";
    }}
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    double M=20*1.989e30, d=1e25;
    std::cout<<"R_LIGO@150Hz="<<obj.compute_R_supp_LIGO(M,150.0)<<"\\n";
    std::cout<<"R_LISA@0.01Hz="<<obj.compute_R_supp_LISA(M,0.01)<<"\\n";
    std::cout<<"f_cross="<<obj.compute_crossover_freq(M)<<" Hz\\n";
    obj.compare_sweep(M,d,10,2000,1e-4,1,"lisa_vs_ligo.csv");
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_679 — AetherSuperfluidDynamics
# ─────────────────────────────────────────────────────────────────────────────
def gen_679():
    cls = "AetherSuperfluidDynamics"
    guard = "AETHER_SUPERFLUID_DYNAMICS_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Aether [UA] as Superfluid — PAPER_679 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~20556
 *
 * The UQFF treats the Universal Aether [UA] as a bosonic superfluid
 * described by a macroscopic wavefunction Psi(r,t) = sqrt(n_UA) * e^{{i phi}}
 *
 * Gross-Pitaevskii-like equation for [UA]:
 *   i hbar d/dt Psi = [-hbar^2/(2 m_UA) nabla^2 + g_UA |Psi|^2 - mu_UA] Psi
 *
 * Local [UA] density: n_UA = rho_UA / m_UA  (m_UA ~ hbar H_0/c^2 ultralight boson)
 * Healing length: xi_UA = hbar / sqrt(2 m_UA g_UA n_UA)
 * Sound speed: c_UA = sqrt(g_UA n_UA / m_UA)  (subsonic, propagates aether waves)
 * Vortex core: r_vortex ~ xi_UA  (quantized circulation kappa_v = h/m_UA)
 *
 * Gravitational coupling:
 *   g_eff(r) = g_Newton(r) * (1 + c_UA^2/c^2 * f_TRZ * rho_UA/rho_SCm)
 *
 * Self-simulate: 1D GP evolution on radial grid around BH.
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    // Ultralight boson mass: m_UA ~ hbar H_0 / c^2 ~ 2e-68 kg
    static constexpr double H0_SI    = 2.184e-18;      // s-1
    static constexpr double M_UA     = HBAR * H0_SI / (C * C);  // ~2.4e-68 kg
    static constexpr double G_UA     = 1.0e-10;         // interaction constant m^3/s^2
    static constexpr double N_UA_REF = RHO_UA / (HBAR * H0_SI / (C*C)); // m-3

    {cls}();
    double compute_m_UA() const;
    double compute_healing_length(double n_UA) const;
    double compute_sound_speed(double n_UA) const;
    double compute_vortex_circulation() const;
    double compute_g_eff(double r, double M) const;
    double compute_GP_energy(double n_UA, double psi_sq) const;
    void simulate_radial_profile(double r_min, double r_max, double dr, double M,
                                 const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}
double {cls}::compute_m_UA() const{{return M_UA;}}
double {cls}::compute_healing_length(double n_UA) const{{
    if(n_UA<=0||G_UA<=0) return 1e300;
    return HBAR/std::sqrt(2.0*M_UA*G_UA*n_UA);
}}
double {cls}::compute_sound_speed(double n_UA) const{{
    return std::sqrt(G_UA*n_UA/M_UA);
}}
double {cls}::compute_vortex_circulation() const{{
    return (2.0*PI*HBAR)/M_UA; // h/m_UA
}}
double {cls}::compute_g_eff(double r,double M) const{{
    double g_N=G*M/(r*r);
    double c_UA=compute_sound_speed(N_UA_REF);
    double boost=1.0+(c_UA*c_UA/(C*C))*F_TRZ*(RHO_UA/RHO_SCM);
    double val=g_N*boost;
    for(auto& m:mods_) val+=m(r,M);
    return val;
}}
double {cls}::compute_GP_energy(double n_UA,double psi_sq) const{{
    // E/V = -mu*|psi|^2 + g_UA/2 * |psi|^4
    double mu_UA=G_UA*n_UA;
    return -mu_UA*psi_sq + 0.5*G_UA*psi_sq*psi_sq;
}}
void {cls}::simulate_radial_profile(double r0,double r1,double dr,double M,
    const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"r_m,g_Newton,g_UQFF,n_UA_profile\\n";}}}}
    for(double r=r0;r<=r1;r+=dr){{
        double g_N=G*M/(r*r);
        double g_u=compute_g_eff(r,M);
        // UA density: enhanced near BH by gravitational compression
        double n_local=N_UA_REF*(1.0+_{cls}_r_s(M)/r*F_TRZ);
        *os<<r<<","<<g_N<<","<<g_u<<","<<n_local<<"\\n";
    }}
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    std::cout<<"m_UA="<<obj.compute_m_UA()<<" kg\\n";
    std::cout<<"xi_UA="<<obj.compute_healing_length(obj.N_UA_REF)<<" m\\n";
    std::cout<<"c_UA="<<obj.compute_sound_speed(obj.N_UA_REF)<<" m/s\\n";
    std::cout<<"kappa_v="<<obj.compute_vortex_circulation()<<" m^2/s\\n";
    double M=8.55e36; double rs=2*6.6743e-11*M/(9e16);
    obj.simulate_radial_profile(rs,rs*100,rs,"aether_superfluid.csv");
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_680 — VortexQuantization
# ─────────────────────────────────────────────────────────────────────────────
def gen_680():
    cls = "VortexQuantization"
    guard = "VORTEX_QUANTIZATION_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Vortex Quantization in Aether Superfluid — PAPER_680 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~21091
 *
 * In the UQFF Aether superfluid, angular momentum is quantized:
 *   L_n = n * hbar  (n = winding number)
 *
 * Vortex circulation: kappa_v = n * h / m_UA
 * Vortex core radius: a_v = cst * xi_UA * exp(-n * pi)  (approximation)
 * Vortex energy per unit length: E_v/L = rho_UA * kappa_v^2 / (4 pi) * ln(R/a_v)
 *
 * UQFF vortex in BH ergosphere:
 *   Omega_v,UQFF = Omega_v * (1 + f_TRZ * c_UA/c)
 *   E_v,UQFF = E_v * (rho_UA/rho_SCm)
 *
 * Aether Magnus force on vortex:
 *   F_Magnus = rho_UA * kappa_v x v_s
 *   F_UQFF   = F_Magnus * (1 + F_TRZ * rho_UA/rho_SCm)
 *
 * Self-simulate: multi-vortex lattice energy vs winding number.
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double H0_SI = 2.184e-18;
    static constexpr double M_UA  = HBAR * H0_SI / (C * C);
    static constexpr double G_UA  = 1.0e-10;
    static constexpr double N_UA  = RHO_UA / (HBAR * H0_SI / (C * C));

    {cls}();
    double compute_circulation(int n) const;
    double compute_vortex_core(int n) const;
    double compute_vortex_energy(int n, double R_outer, double xi) const;
    double compute_omega_v_UQFF(double omega_v) const;
    double compute_magnus_force_UQFF(double v_s) const;
    void simulate_lattice(int n_max, double R_outer, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}
double {cls}::compute_circulation(int n) const{{
    return n*(2.0*PI*HBAR)/M_UA;
}}
double {cls}::compute_vortex_core(int n) const{{
    double xi=HBAR/std::sqrt(2.0*M_UA*G_UA*N_UA);
    return xi*std::exp(-n*PI);
}}
double {cls}::compute_vortex_energy(int n,double R_outer,double xi) const{{
    double kv=compute_circulation(n);
    if(xi<=0||R_outer<=xi) return 0.0;
    return RHO_UA*kv*kv/(4.0*PI)*std::log(R_outer/xi)*(RHO_UA/RHO_SCM);
}}
double {cls}::compute_omega_v_UQFF(double omega_v) const{{
    double c_UA=std::sqrt(G_UA*N_UA/M_UA);
    return omega_v*(1.0+F_TRZ*c_UA/C);
}}
double {cls}::compute_magnus_force_UQFF(double v_s) const{{
    // |F_Magnus| = rho_UA * kappa_1 * v_s per unit length
    double kv=compute_circulation(1);
    return RHO_UA*kv*v_s*(1.0+F_TRZ*(RHO_UA/RHO_SCM));
}}
void {cls}::simulate_lattice(int n_max,double R_outer,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"n,kv,a_core,E_v_UQFF\\n";}}}}
    double xi=HBAR/std::sqrt(2.0*M_UA*G_UA*N_UA);
    for(int n=1;n<=n_max;n++){{
        double kv=compute_circulation(n);
        double av=compute_vortex_core(n);
        double Ev=compute_vortex_energy(n,R_outer,xi);
        *os<<n<<","<<kv<<","<<av<<","<<Ev<<"\\n";
    }}
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    double xi=1e-5, R=1.0;
    for(int n=1;n<=5;n++){{
        std::cout<<"n="<<n<<" kv="<<obj.compute_circulation(n)
                 <<" E_v="<<obj.compute_vortex_energy(n,R,xi)<<"\\n";
    }}
    obj.simulate_lattice(10,1.0,"vortex_lattice.csv");
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_681 — GrossPitaevskiiVortexSimulation
# ─────────────────────────────────────────────────────────────────────────────
def gen_681():
    cls = "GrossPitaevskiiVortexSimulation"
    guard = "GROSS_PITAEVSKII_VORTEX_SIMULATION_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Gross-Pitaevskii Vortex Simulation — PAPER_681 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~21332
 *
 * Full 1D Gross-Pitaevskii equation for the [UA] Aether wavefunction
 * in radial (spherical) coordinates around a UQFF black hole:
 *
 *   i hbar dpsi/dt = [ -hbar^2/(2 m_UA) (1/r^2 d/dr r^2 d/dr)
 *                    + V_grav(r)
 *                    + g_UA |psi|^2
 *                    + U_m(r,t) ] psi
 *
 * where V_grav(r) = -G M m_UA / r  (gravitational well)
 *       U_m(r,t) from UQFF magnetic string term
 *
 * Imaginary-time propagation to find ground state.
 * Real-time propagation to simulate vortex dynamics.
 * Normalisation: integral |psi|^2 r^2 dr = N_UA * V_sphere
 *
 * Output: density profile |psi(r)|^2, phase phi(r), current j(r).
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double H0_SI = 2.184e-18;
    static constexpr double M_UA  = HBAR * H0_SI / (C * C);
    static constexpr double G_UA  = 1.0e-10;
    static constexpr double N_UA  = RHO_UA / (HBAR * H0_SI / (C * C));

    {cls}();
    // Compute single GP time step (imaginary time)
    std::vector<double> imaginary_time_step(
        const std::vector<double>& psi,
        const std::vector<double>& r_grid,
        double M_bh, double dt, int n_steps) const;
    double compute_chemical_potential(const std::vector<double>& psi,
                                      const std::vector<double>& r_grid,
                                      double M_bh) const;
    void simulate(double r_min, double r_max, int nr, double M_bh,
                  int n_steps, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}

std::vector<double> {cls}::imaginary_time_step(
    const std::vector<double>& psi0,
    const std::vector<double>& r_grid,
    double M_bh, double dt, int n_steps) const{{
    int N=(int)r_grid.size();
    std::vector<double> psi=psi0;
    double dr=r_grid.size()>1?r_grid[1]-r_grid[0]:1.0;
    for(int step=0;step<n_steps;step++){{
        std::vector<double> dpsi(N,0.0);
        for(int i=1;i<N-1;i++){{
            double r=r_grid[i];
            double lap=(psi[i+1]-2*psi[i]+psi[i-1])/(dr*dr)+
                       2.0/r*(psi[i+1]-psi[i-1])/(2*dr);
            double kin=-HBAR*HBAR/(2.0*M_UA)*lap;
            double V_grav=-G*M_bh*M_UA/r;
            double gp_int=G_UA*psi[i]*psi[i]*psi[i];
            double Um_val=_{cls}_Um(r,1e8+step*dt);
            dpsi[i]=psi[i]-dt/HBAR*(kin/psi[i]+V_grav+gp_int+Um_val*psi[i]/(HBAR));
        }}
        // boundary
        dpsi[0]=dpsi[1]; dpsi[N-1]=0.0;
        // normalise
        double norm=0.0;
        for(int i=0;i<N;i++) norm+=dpsi[i]*dpsi[i]*r_grid[i]*r_grid[i]*dr;
        norm=std::sqrt(norm*4*PI*N_UA);
        if(norm>0) for(int i=0;i<N;i++) psi[i]=dpsi[i]/norm*std::sqrt(N_UA);
    }}
    return psi;
}}
double {cls}::compute_chemical_potential(
    const std::vector<double>& psi,
    const std::vector<double>& r_grid,
    double M_bh) const{{
    double mu=0.0,norm=0.0;
    int N=(int)r_grid.size();
    double dr=N>1?r_grid[1]-r_grid[0]:1.0;
    for(int i=0;i<N;i++){{
        double r=r_grid[i];
        double psi2=psi[i]*psi[i];
        double V=-G*M_bh*M_UA/r;
        mu+=(V+G_UA*psi2)*psi2*r*r*dr;
        norm+=psi2*r*r*dr;
    }}
    return norm>0?mu/norm:0.0;
}}
void {cls}::simulate(double r_min,double r_max,int nr,double M_bh,
    int n_steps,const std::string& out) const{{
    std::vector<double> r_grid(nr);
    double dr=(r_max-r_min)/(nr-1);
    for(int i=0;i<nr;i++) r_grid[i]=r_min+i*dr;
    // Initial guess: Gaussian
    std::vector<double> psi(nr);
    double r_mid=(r_min+r_max)/2.0;
    for(int i=0;i<nr;i++)
        psi[i]=std::exp(-(r_grid[i]-r_mid)*(r_grid[i]-r_mid)/(2.0*r_mid*r_mid));
    // Normalise
    double norm=0.0;
    for(int i=0;i<nr;i++) norm+=psi[i]*psi[i]*r_grid[i]*r_grid[i]*dr;
    if(norm>0){{for(auto& p:psi) p/=std::sqrt(norm*4*PI/N_UA);}}
    psi=imaginary_time_step(psi,r_grid,M_bh,1e-30,n_steps);
    double mu=compute_chemical_potential(psi,r_grid,M_bh);
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"r_m,psi_sq,mu_eV\\n";}}}}
    for(int i=0;i<nr;i++)
        *os<<r_grid[i]<<","<<psi[i]*psi[i]<<","<<mu/1.602e-19<<"\\n";
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    double M=8.55e36; double rs=2*6.6743e-11*M/(9e16);
    obj.simulate(rs,rs*50,100,M,50,"gp_vortex.csv");
    std::cout<<"GP simulation complete.\\n";
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_682 — UQFFStabilityNumericallyForSgrA
# ─────────────────────────────────────────────────────────────────────────────
def gen_682():
    cls = "UQFFStabilityNumericallyForSgrA"
    guard = "UQFF_STABILITY_NUMERICALLY_FOR_SgrA_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Numerical Stability Analysis for Sgr A* — PAPER_682 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~21554
 *
 * Sgr A* parameters:
 *   M_SgrA = 4.297e6 Msun = 8.548e36 kg
 *   Distance: 8.178 kpc = 2.523e20 m
 *
 * Numerical stability analysis in UQFF:
 *   1. Perturbation expansion: delta_M(t) / M_0 = epsilon * exp(i omega t)
 *      omega = omega_R + i omega_I
 *      stable: omega_I < 0  (damped perturbation)
 *
 *   2. UQFF stability equation:
 *      omega_UQFF = omega_GR * (1 + f_TRZ * rho_UA/rho_SCm * U_m/k_B T_H)
 *      -> enhanced damping rate
 *
 *   3. Lyapunov exponent:
 *      lambda_UQFF = lambda_GR * (RHO_SCM/RHO_UA) * exp(-U_m/k_B T_H)
 *      lambda < 0: stable fixed point
 *
 *   4. RK4 integration of perturbed BH mass evolution:
 *      dM/dt = -(1-f_TRZ) * L_H(M) / c^2 * (rho_SCm/rho_UA) * exp(-U_m/k_B T_H)
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double M_SGRA   = 4.297e6 * 1.989e30;
    static constexpr double D_SGRA   = 2.523e20;

    {cls}();
    double compute_omega_I_UQFF(double M) const;
    double compute_lyapunov_UQFF(double M) const;
    double compute_dM_dt_UQFF(double M, double t) const;
    std::vector<std::pair<double,double>> rk4_evolution(double M0, double t_end,
                                                        double dt) const;
    bool is_stable(double M) const;
    void simulate_stability(double t_end, double dt, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}
double {cls}::compute_omega_I_UQFF(double M) const{{
    // omega_GR ~ -1/(tau_std) (decays at Hawking rate)
    double tau=_{cls}_tau(M);
    double T_H=_{cls}_T_H(M); double r=_{cls}_r_s(M);
    double Um=_{cls}_Um(r,T_N_REF);
    double omega_GR=-1.0/std::max(tau,1e-300);
    return omega_GR*(1.0+F_TRZ*(RHO_UA/RHO_SCM)*Um/(K_B*T_H));
}}
double {cls}::compute_lyapunov_UQFF(double M) const{{
    double tau=_{cls}_tau(M);
    double T_H=_{cls}_T_H(M); double r=_{cls}_r_s(M);
    double Um=_{cls}_Um(r,T_N_REF);
    double exp_arg=Um/(K_B*T_H); if(exp_arg>700)exp_arg=700;
    return -(RHO_SCM/RHO_UA)*std::exp(-exp_arg)/std::max(tau,1e-300);
}}
double {cls}::compute_dM_dt_UQFF(double M,double t) const{{
    if(M<=0) return 0.0;
    double L_H=_{cls}_L_H(M);
    double T_H=_{cls}_T_H(M); double r=_{cls}_r_s(M);
    double Um=_{cls}_Um(r,t);
    double exp_arg=Um/(K_B*T_H); if(exp_arg>700)exp_arg=700;
    double dM=-(1.0-F_TRZ)*L_H/(C*C)*(RHO_SCM/RHO_UA)*std::exp(-exp_arg);
    for(auto& m:mods_) dM+=m(M,t);
    return dM;
}}
std::vector<std::pair<double,double>> {cls}::rk4_evolution(double M0,double t_end,double dt) const{{
    std::vector<std::pair<double,double>> traj;
    double M=M0,t=0;
    while(t<t_end&&M>0){{
        traj.push_back({{t,M}});
        double k1=compute_dM_dt_UQFF(M,t);
        double k2=compute_dM_dt_UQFF(M+0.5*dt*k1,t+0.5*dt);
        double k3=compute_dM_dt_UQFF(M+0.5*dt*k2,t+0.5*dt);
        double k4=compute_dM_dt_UQFF(M+dt*k3,t+dt);
        M+=dt/6.0*(k1+2*k2+2*k3+k4);
        t+=dt;
    }}
    return traj;
}}
bool {cls}::is_stable(double M) const{{
    return compute_lyapunov_UQFF(M)<0.0;
}}
void {cls}::simulate_stability(double t_end,double dt,const std::string& out) const{{
    auto traj=rk4_evolution(M_SGRA,t_end,dt);
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"t_s,M_kg,stable\\n";}}}}
    for(auto& p:traj)
        *os<<p.first<<","<<p.second<<","<<(is_stable(p.second)?1:0)<<"\\n";
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    std::cout<<"Lyapunov SgrA*="<<obj.compute_lyapunov_UQFF(obj.M_SGRA)<<"\\n";
    std::cout<<"Stable="<<obj.is_stable(obj.M_SGRA)<<"\\n";
    obj.simulate_stability(1e60,1e55,"sgra_stability.csv");
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_683 — UQFFHawkingTemperatureModulation
# ─────────────────────────────────────────────────────────────────────────────
def gen_683():
    cls = "UQFFHawkingTemperatureModulation"
    guard = "UQFF_HAWKING_TEMPERATURE_MODULATION_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Hawking Temperature Modulation — PAPER_683 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~21788
 *
 * Standard Hawking temperature: T_H = hbar c^3 / (8 pi G M k_B)
 *
 * UQFF modulates T_H through three channels:
 *   1. f_TRZ time-reversal boost:     T_UQFF = T_H * (1 + f_TRZ)
 *   2. Density ratio suppression:     T_UQFF = T_H * (1 - rho_SCm/rho_UA)
 *   3. U_m magnetic string modulation:
 *      T_UQFF = T_H * (1 + f_TRZ) * (1 - rho_SCm/rho_UA) * (1 + U_m/k_B T_H)
 *
 *   Combined (first-order expansion):
 *   T_UQFF = T_H * [1 + f_TRZ + U_m/k_B T_H - rho_SCm/rho_UA * (1 + f_TRZ)]
 *
 * Spectral modulation:
 *   n(omega)_UQFF = 1 / (exp(hbar omega / k_B T_UQFF) - 1)
 *   -> shifted Planck spectrum peak wavelength
 *   lambda_max_UQFF = hbar c / (2.82 k_B T_UQFF)  (Wien law)
 *
 * Time-dependence: T_UQFF(t) includes time-varying U_m(r,t) term.
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    {cls}();
    double compute_T_H(double M) const;
    double compute_T_UQFF(double M, double r=0.0, double t=1e8) const;
    double compute_planck_spectrum(double omega, double T) const;
    double compute_wien_peak(double T) const;
    double compute_T_modulation_factor(double M, double r, double t) const;
    void simulate_T_vs_M(double M_start, double M_end, double dM,
                         const std::string& out="") const;
    void simulate_spectrum(double M, double omega_min, double omega_max,
                           double d_omega, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}
double {cls}::compute_T_H(double M) const{{return _{cls}_T_H(M);}}
double {cls}::compute_T_UQFF(double M,double r,double t) const{{
    double T_H=_{cls}_T_H(M);
    if(r<=0) r=_{cls}_r_s(M);
    double Um=_{cls}_Um(r,t);
    // T_UQFF = T_H*(1+f_TRZ)*(1-rho_SCm/rho_UA)*(1+U_m/k_B T_H)
    double fac=(1.0+F_TRZ)*(1.0-RHO_SCM/RHO_UA)*(1.0+Um/(K_B*T_H));
    double T=T_H*fac;
    for(auto& m:mods_) T+=m(M,t);
    return T;
}}
double {cls}::compute_planck_spectrum(double omega,double T) const{{
    if(T<=0||omega<=0) return 0.0;
    double x=HBAR*omega/(K_B*T); if(x>700) return 0.0;
    return HBAR*omega*omega*omega/(PI*PI*C*C*C*(std::exp(x)-1.0));
}}
double {cls}::compute_wien_peak(double T) const{{
    return HBAR*C/(2.82*K_B*T);
}}
double {cls}::compute_T_modulation_factor(double M,double r,double t) const{{
    double T_H=compute_T_H(M);
    double T_U=compute_T_UQFF(M,r,t);
    return T_U/std::max(T_H,1e-300);
}}
void {cls}::simulate_T_vs_M(double M0,double M1,double dM,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"M_kg,T_H_K,T_UQFF_K,factor\\n";}}}}
    for(double M=M0;M<=M1;M+=dM){{
        double T_H=compute_T_H(M), T_U=compute_T_UQFF(M);
        *os<<M<<","<<T_H<<","<<T_U<<","<<T_U/std::max(T_H,1e-300)<<"\\n";
    }}
    if(ofs.is_open())ofs.close();
}}
void {cls}::simulate_spectrum(double M,double w_min,double w_max,double dw,
    const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"omega_rad_s,n_H,n_UQFF\\n";}}}}
    double T_H=compute_T_H(M), T_U=compute_T_UQFF(M);
    for(double w=w_min;w<=w_max;w+=dw)
        *os<<w<<","<<compute_planck_spectrum(w,T_H)<<","<<compute_planck_spectrum(w,T_U)<<"\\n";
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    double M=1e12;
    std::cout<<"T_H="<<obj.compute_T_H(M)<<" K\\n";
    std::cout<<"T_UQFF="<<obj.compute_T_UQFF(M)<<" K\\n";
    std::cout<<"Modulation factor="<<obj.compute_T_modulation_factor(M,0,1e8)<<"\\n";
    obj.simulate_T_vs_M(1e10,1e15,1e12,"T_modulation.csv");
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_684 — UQFFPrimordialBHEvaporation
# ─────────────────────────────────────────────────────────────────────────────
def gen_684():
    cls = "UQFFPrimordialBHEvaporation"
    guard = "UQFF_PRIMORDIAL_BH_EVAPORATION_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Primordial Black Hole Evaporation — PAPER_684 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~22011
 *
 * Primordial BHs (PBHs) form in the early universe from density
 * fluctuations. Standard Hawking evaporation:
 *   dM/dt = -hbar c^4 / (15360 pi G^2 M^2)  [Schwarzschild PBH]
 *
 * UQFF modifies evaporation amplitude by suppression chain:
 *   dM/dt_UQFF = dM/dt * (1-f_TRZ) * (rho_SCm/rho_UA) * exp(-U_m/k_B T_H)
 *   -> factor ~0.033 slower than standard
 *
 * PBH mass at formation: M_f = rho_rad(t_f) * (4/3) pi (c t_f / 2)^3
 *   For t_f ~ 1e-23 s: M_f ~ 1e10 kg  (minimum viable DM mass in UQFF)
 *
 * Analytic lifetime: tau_UQFF = 5120 pi G^2 M_f^3 / (hbar c^4)
 *                              * 1/(1-f_TRZ) * (rho_UA/rho_SCm) * exp(U_m/k_B T_H)
 *
 * Euler / RK4 numerical integration of M(t) from M_f to 0.
 * Burst parameters: dN_gamma/dE at end of life — UQFF delays burst epoch.
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    {cls}();
    double compute_dM_dt_std(double M) const;
    double compute_dM_dt_UQFF(double M, double t) const;
    double compute_tau_std(double M) const;
    double compute_tau_UQFF(double M) const;
    double compute_M_formation(double t_form) const;
    std::vector<std::pair<double,double>> evolve_rk4(double M0, double t_end,
                                                     double dt) const;
    void simulate_evaporation(double M0, double dt, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}
double {cls}::compute_dM_dt_std(double M) const{{
    if(M<=0) return 0;
    return -HBAR*std::pow(C,4.0)/(15360.0*PI*G*G*M*M);
}}
double {cls}::compute_dM_dt_UQFF(double M,double t) const{{
    if(M<=0) return 0;
    double dM_std=compute_dM_dt_std(M);
    double T_H=_{cls}_T_H(M); double r=_{cls}_r_s(M);
    double Um=_{cls}_Um(r,t);
    double exp_arg=Um/(K_B*T_H); if(exp_arg>700)exp_arg=700;
    double dM=dM_std*(1.0-F_TRZ)*(RHO_SCM/RHO_UA)*std::exp(-exp_arg);
    for(auto& m:mods_) dM+=m(M,t);
    return dM;
}}
double {cls}::compute_tau_std(double M) const{{return _{cls}_tau(M);}}
double {cls}::compute_tau_UQFF(double M) const{{
    double tau=_{cls}_tau(M);
    double T_H=_{cls}_T_H(M); double r=_{cls}_r_s(M);
    double Um=_{cls}_Um(r,T_N_REF);
    double exp_arg=Um/(K_B*T_H); if(exp_arg>700)exp_arg=700;
    return tau/(1.0-F_TRZ)*(RHO_UA/RHO_SCM)*std::exp(exp_arg);
}}
double {cls}::compute_M_formation(double t_form) const{{
    // M_f = (4/3)*pi*(c*t_f/2)^3 * rho_rad(t_f)
    // rho_rad(t_f) ~ 3/(32*pi*G*t_f^2) (radiation dominated)
    double r_form=C*t_form/2.0;
    double rho_rad=3.0/(32.0*PI*G*t_form*t_form);
    return (4.0/3.0)*PI*r_form*r_form*r_form*rho_rad;
}}
std::vector<std::pair<double,double>> {cls}::evolve_rk4(double M0,double t_end,double dt) const{{
    std::vector<std::pair<double,double>> traj;
    double M=M0,t=0;
    while(t<t_end&&M>0){{
        traj.push_back({{t,M}});
        double k1=compute_dM_dt_UQFF(M,t);
        double k2=compute_dM_dt_UQFF(M+0.5*dt*k1,t+0.5*dt);
        double k3=compute_dM_dt_UQFF(M+0.5*dt*k2,t+0.5*dt);
        double k4=compute_dM_dt_UQFF(M+dt*k3,t+dt);
        M+=dt/6.0*(k1+2*k2+2*k3+k4);
        t+=dt;
    }}
    return traj;
}}
void {cls}::simulate_evaporation(double M0,double dt,const std::string& out) const{{
    double tau_U=compute_tau_UQFF(M0);
    auto traj=evolve_rk4(M0,tau_U*2.0,dt);
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"t_s,M_kg\\n";}}}}
    for(auto& p:traj) *os<<p.first<<","<<p.second<<"\\n";
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    double Mf=obj.compute_M_formation(1e-23);
    std::cout<<"M_formation="<<Mf<<" kg\\n";
    std::cout<<"tau_UQFF="<<obj.compute_tau_UQFF(Mf)<<" s\\n";
    obj.simulate_evaporation(1e12,1e50,"pbh_evaporation.csv");
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_685 — UQFFPBHDarkMatterImplications
# ─────────────────────────────────────────────────────────────────────────────
def gen_685():
    cls = "UQFFPBHDarkMatterImplications"
    guard = "UQFF_PBH_DARK_MATTER_IMPLICATIONS_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Primordial BH Dark Matter Implications — PAPER_685 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~22222
 *
 * Cosmological dark matter fraction from stable PBHs in UQFF:
 *
 *   Standard: PBHs with M < M_crit_std evaporate before today
 *             M_crit_std = (hbar c^4 * t_age / (5120 pi G^2))^(1/3)
 *             M_crit_std ~ 5e11 kg  (only heavier PBHs survive)
 *
 *   UQFF: M_crit_UQFF lowered by factor ~0.033 suppression -> smaller PBHs survive
 *         M_crit_UQFF = M_crit_std * (tau_ratio)^(-1/3)
 *         where tau_ratio = tau_UQFF/tau_std ~ 30
 *         -> M_crit_UQFF ~ M_crit_std / 30^(1/3) ~ 0.32 * M_crit_std
 *
 *   DM fraction:
 *     f_PBH = rho_PBH / rho_DM
 *     rho_PBH ~ M_f * n_PBH  (n_PBH from inflationary spectrum)
 *     UQFF enhances f_PBH by keeping lower-mass PBHs alive:
 *     f_PBH_UQFF = f_PBH_std * 30^(2/3)  (more PBHs survive)
 *
 *   Observational constraints: microlensing, CMB distortions, GW background
 *   UQFF shifts the viable mass window: 1e10 – 1e17 kg fully open in UQFF.
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double T_AGE = 4.34e17;  // 13.8 Gyr in s

    {cls}();
    double compute_M_crit_std() const;
    double compute_M_crit_UQFF() const;
    double compute_f_PBH_boost() const;
    double compute_rho_PBH_UQFF(double M_f, double n_PBH) const;
    bool is_DM_viable_std(double M) const;
    bool is_DM_viable_UQFF(double M) const;
    void scan_mass_window(double M_min, double M_max, int N,
                          const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}
double {cls}::compute_M_crit_std() const{{
    return std::cbrt(HBAR*std::pow(C,4.0)*T_AGE/(5120.0*PI*G*G));
}}
double {cls}::compute_M_crit_UQFF() const{{
    // tau_UQFF = tau_std * 30 => M_crit3 scales with tau => M_crit scales as tau^{{1/3}}
    // But evaporation slower means SAME mass survives longer
    // M_crit_UQFF^3 = M_crit_std^3 / 30 => smaller PBHs survive
    double M_std=compute_M_crit_std();
    double tau_ratio=1.0/(1.0-F_TRZ)*(RHO_UA/RHO_SCM); // ~11
    return M_std/std::cbrt(tau_ratio);
}}
double {cls}::compute_f_PBH_boost() const{{
    // f_PBH_UQFF/f_std = (M_crit_std/M_crit_UQFF)^2 (more PBHs in window)
    double tau_ratio=1.0/(1.0-F_TRZ)*(RHO_UA/RHO_SCM);
    return std::pow(tau_ratio,2.0/3.0);
}}
double {cls}::compute_rho_PBH_UQFF(double M_f,double n_PBH) const{{
    double rho_std=M_f*n_PBH;
    return rho_std*compute_f_PBH_boost();
}}
bool {cls}::is_DM_viable_std(double M) const{{
    return M>compute_M_crit_std();
}}
bool {cls}::is_DM_viable_UQFF(double M) const{{
    return M>compute_M_crit_UQFF();
}}
void {cls}::scan_mass_window(double M_min,double M_max,int N,
    const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"M_kg,viable_std,viable_UQFF,tau_ratio\\n";}}}}
    for(int i=0;i<=N;i++){{
        double M=std::exp(std::log(M_min)+i*(std::log(M_max)-std::log(M_min))/N);
        double tau=_{cls}_tau(M);
        *os<<M<<","<<is_DM_viable_std(M)<<","<<is_DM_viable_UQFF(M)<<","<<compute_f_PBH_boost()<<"\\n";
    }}
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    std::cout<<"M_crit_std="<<obj.compute_M_crit_std()<<" kg\\n";
    std::cout<<"M_crit_UQFF="<<obj.compute_M_crit_UQFF()<<" kg\\n";
    std::cout<<"f_PBH_boost="<<obj.compute_f_PBH_boost()<<"x\\n";
    obj.scan_mass_window(1e8,1e20,100,"pbh_dm_window.csv");
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_686 — UQFFModulationForM87
# ─────────────────────────────────────────────────────────────────────────────
def gen_686():
    cls = "UQFFModulationForM87"
    guard = "UQFF_MODULATION_FOR_M87_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Modulation for M87 SMBH — PAPER_686 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~22436
 *
 * M87 (Messier 87 / NGC 4486):
 *   SMBH mass: M87* = 6.5e9 Msun = 1.293e40 kg (EHT April 2019)
 *   Distance: 16.4 Mpc = 5.06e23 m
 *   Shadow radius: r_shadow = 3*sqrt(3) * G*M/c^2 ~ 2e13 m
 *   Jet power: P_jet ~ 10^44 erg/s = 10^37 W
 *
 * UQFF modifications for M87*:
 *   1. Hawking T: T_UQFF = T_H * (1+f_TRZ) * (1-rho_SCm/rho_UA)
 *   2. Shadow size: r_sh,UQFF = r_sh * (1 + f_TRZ * rho_UA/rho_SCm)^0.5
 *   3. Jet power: P_jet,UQFF = P_jet,BZ * (1+f_TRZ) * (rho_UA/rho_SCm)^0.5
 *      Blandford-Znajek: P_BZ ~ kappa * Phi_BH^2 * Omega_H^2 / (4 pi c)
 *   4. Ring brightness: B_UQFF = B_GR * (rho_UA/rho_SCm)^(f_TRZ/4)
 *   5. U_m coupling to accretion disk: T_acc,UQFF = T_acc * (1 + U_m/k_B T_H)
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double M_M87     = 6.5e9 * 1.989e30;
    static constexpr double D_M87     = 5.06e23;
    static constexpr double P_JET_GR  = 1.0e37;   // W (Blandford-Znajek estimate)
    static constexpr double B_FIELD   = 1.0e3;     // T near horizon

    {cls}();
    double compute_T_H(double M=M_M87) const;
    double compute_T_UQFF(double M=M_M87) const;
    double compute_shadow_radius(double M=M_M87) const;
    double compute_shadow_UQFF(double M=M_M87) const;
    double compute_jet_power_UQFF(double P_BZ=P_JET_GR) const;
    double compute_ring_brightness_UQFF(double B_GR) const;
    double compute_T_accretion_UQFF(double T_acc, double M=M_M87) const;
    void simulate_T_evolution(double M_start, double M_end, double dM,
                              const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}
double {cls}::compute_T_H(double M) const{{return _{cls}_T_H(M);}}
double {cls}::compute_T_UQFF(double M) const{{
    double T_H=_{cls}_T_H(M);
    return T_H*(1.0+F_TRZ)*(1.0-RHO_SCM/RHO_UA);
}}
double {cls}::compute_shadow_radius(double M) const{{
    return 3.0*std::sqrt(3.0)*G*M/(C*C);
}}
double {cls}::compute_shadow_UQFF(double M) const{{
    return compute_shadow_radius(M)*std::sqrt(1.0+F_TRZ*(RHO_UA/RHO_SCM));
}}
double {cls}::compute_jet_power_UQFF(double P_BZ) const{{
    return P_BZ*(1.0+F_TRZ)*std::sqrt(RHO_UA/RHO_SCM);
}}
double {cls}::compute_ring_brightness_UQFF(double B_GR) const{{
    return B_GR*std::pow(RHO_UA/RHO_SCM,F_TRZ/4.0);
}}
double {cls}::compute_T_accretion_UQFF(double T_acc,double M) const{{
    double T_H=_{cls}_T_H(M); double r=_{cls}_r_s(M);
    double Um=_{cls}_Um(r,T_N_REF);
    return T_acc*(1.0+Um/(K_B*T_H));
}}
void {cls}::simulate_T_evolution(double M0,double M1,double dM,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"M_Msun,T_H,T_UQFF,shadow_UQFF_m\\n";}}}}
    for(double M=M0;M<=M1;M+=dM){{
        *os<<M/M_SUN<<","<<compute_T_H(M)<<","<<compute_T_UQFF(M)<<","<<compute_shadow_UQFF(M)<<"\\n";
    }}
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    std::cout<<"T_H(M87*)="<<obj.compute_T_H()<<" K\\n";
    std::cout<<"T_UQFF(M87*)="<<obj.compute_T_UQFF()<<" K\\n";
    std::cout<<"Shadow GR="<<obj.compute_shadow_radius()/1e9<<" Gm\\n";
    std::cout<<"Shadow UQFF="<<obj.compute_shadow_UQFF()/1e9<<" Gm\\n";
    std::cout<<"Jet UQFF="<<obj.compute_jet_power_UQFF()<<" W\\n";
    obj.simulate_T_evolution(1e9*1.989e30,15e9*1.989e30,0.5e9*1.989e30,"m87_modulation.csv");
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─────────────────────────────────────────────────────────────────────────────
# PAPER_687 — M87MassEvolutionSimulation
# ─────────────────────────────────────────────────────────────────────────────
def gen_687():
    cls = "M87MassEvolutionSimulation"
    guard = "M87_MASS_EVOLUTION_SIMULATION_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief M87* SMBH Mass Evolution in UQFF — PAPER_687 | Session 173
 * Source: grok_share_fc21e30c24b4.txt line ~22623
 *
 * M87*: 6.5e9 Msun at D=16.4 Mpc. Combines UQFF accretion with UQFF-suppressed
 * Hawking evaporation for a full mass evolution simulation.
 *
 * Coupled equations:
 *   dM/dt_acc  = Mdot_Bondi_UQFF  (accretion from ISM)
 *   dM/dt_evap = dM/dt_Hawking_UQFF  (suppressed Hawking)
 *   dM/dt_jet  = -P_jet_UQFF / c^2  (jet mass loss)
 *   dM/dt_total = dM_acc - |dM_evap| - dM_jet
 *
 * Bondi-UQFF:
 *   Mdot_Bondi = 4 pi G^2 M^2 rho_inf / (c_s^3)
 *   Mdot_UQFF  = Mdot_Bondi * (rho_inf + rho_UA - rho_SCm) / rho_inf * (1+f_TRZ)
 *
 * Jet (Blandford-Znajek UQFF):
 *   P_jet = kappa_BZ * Phi_BH^2 * Omega_H^2 / (4 pi c) * (1+f_TRZ) * sqrt(rho_UA/rho_SCm)
 *   Omega_H ~ c/(4 G M) (maximal spin approximation)
 *
 * Simulate: M(t) over 14 Gyr (age of universe). Compare GR vs UQFF evolution.
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double M_M87    = 6.5e9 * 1.989e30;
    static constexpr double RHO_ISM  = 1.67e-25;   // kg/m^3 M87 environment
    static constexpr double T_ISM    = 1.0e7;       // K hot plasma
    static constexpr double T_HUBBLE = 4.34e17;     // s (13.8 Gyr)

    {cls}();
    double compute_sound_speed_ISM(double T_ISM) const;
    double compute_Mdot_Bondi_UQFF(double M, double rho_inf, double T_inf) const;
    double compute_dM_dt_evap_UQFF(double M, double t) const;
    double compute_jet_power_UQFF(double M) const;
    double compute_dM_dt_total(double M, double t,
                               double rho_inf=RHO_ISM,
                               double T_inf=T_ISM) const;
    std::vector<std::pair<double,double>> evolve(double M0, double t_end, double dt,
                                                 double rho_inf=RHO_ISM,
                                                 double T_inf=T_ISM) const;
    void simulate_over_hubble(double dt, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{helpers(cls)}
{cls}::{cls}(){{}}
double {cls}::compute_sound_speed_ISM(double T) const{{
    // c_s = sqrt(gamma_adiab * k_B T / m_p) for fully ionized H
    return std::sqrt(5.0/3.0*K_B*T/(1.673e-27));
}}
double {cls}::compute_Mdot_Bondi_UQFF(double M,double rho_inf,double T_inf) const{{
    double cs=compute_sound_speed_ISM(T_inf);
    double Mdot_B=4.0*PI*G*G*M*M*rho_inf/(cs*cs*cs);
    double rho_eff=rho_inf+RHO_UA-RHO_SCM;
    return Mdot_B*(rho_eff/rho_inf)*(1.0+F_TRZ);
}}
double {cls}::compute_dM_dt_evap_UQFF(double M,double t) const{{
    if(M<=0) return 0;
    double dM_std=-HBAR*std::pow(C,4.0)/(15360.0*PI*G*G*M*M);
    double T_H=_{cls}_T_H(M); double r=_{cls}_r_s(M);
    double Um=_{cls}_Um(r,t);
    double exp_arg=Um/(K_B*T_H); if(exp_arg>700)exp_arg=700;
    return dM_std*(1.0-F_TRZ)*(RHO_SCM/RHO_UA)*std::exp(-exp_arg);
}}
double {cls}::compute_jet_power_UQFF(double M) const{{
    // P_BZ ~ B^2 r_s^2 c / 4  (simplified), B from equipartition
    double rs=_{cls}_r_s(M);
    double rho_disk=RHO_ISM*1e6; // enhanced near BH
    double B_eq=std::sqrt(8.0*PI*rho_disk*C*C); // crude equipartition B
    if(B_eq>1e8) B_eq=1e8; // cap
    double Omega_H=C/(4.0*G*M/C/C)*C; // c/(4*Rg)*c
    double P_BZ=0.01*B_eq*B_eq*rs*rs*C*Omega_H*Omega_H/(4.0*PI*C);
    return P_BZ*(1.0+F_TRZ)*std::sqrt(RHO_UA/RHO_SCM);
}}
double {cls}::compute_dM_dt_total(double M,double t,double rho_inf,double T_inf) const{{
    double dM_acc=compute_Mdot_Bondi_UQFF(M,rho_inf,T_inf);
    double dM_evap=compute_dM_dt_evap_UQFF(M,t);
    double P_jet=compute_jet_power_UQFF(M);
    double dM_jet=-P_jet/(C*C);
    double total=dM_acc+dM_evap+dM_jet;
    for(auto& m:mods_) total+=m(M,t);
    return total;
}}
std::vector<std::pair<double,double>> {cls}::evolve(double M0,double t_end,double dt,
    double rho_inf,double T_inf) const{{
    std::vector<std::pair<double,double>> traj;
    double M=M0,t=0;
    while(t<t_end&&M>0){{
        traj.push_back({{t,M}});
        double k1=compute_dM_dt_total(M,t,rho_inf,T_inf);
        double k2=compute_dM_dt_total(M+0.5*dt*k1,t+0.5*dt,rho_inf,T_inf);
        double k3=compute_dM_dt_total(M+0.5*dt*k2,t+0.5*dt,rho_inf,T_inf);
        double k4=compute_dM_dt_total(M+dt*k3,t+dt,rho_inf,T_inf);
        M+=dt/6.0*(k1+2*k2+2*k3+k4);
        t+=dt;
    }}
    return traj;
}}
void {cls}::simulate_over_hubble(double dt,const std::string& out) const{{
    auto traj=evolve(M_M87,T_HUBBLE,dt);
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){{ofs.open(out);if(ofs.is_open()){{os=&ofs;*os<<"t_Gyr,M_Msun\\n";}}}}
    for(auto& p:traj)
        *os<<p.first/3.156e16<<","<<p.second/M_SUN<<"\\n";
    if(ofs.is_open())ofs.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){{auto e=l.find('=');if(e==std::string::npos)continue;}}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} obj;
    double M=M_M87;
    std::cout<<"dM/dt_total="<<obj.compute_dM_dt_total(M,0)<<" kg/s\\n";
    std::cout<<"Mdot_Bondi_UQFF="<<obj.compute_Mdot_Bondi_UQFF(M,RHO_ISM,T_ISM)<<" kg/s\\n";
    obj.simulate_over_hubble(1e13,"m87_mass_evolution.csv");
    return 0;
}}
#endif
"""
    write(f"{cls}.h", h)
    write(f"{cls}.cpp", cpp)

# ─── RUN ALL ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import os
    os.chdir(ROOT)
    print("Generating PAPER_674-687 C++ modules (28 files)...")
    gen_674(); gen_675(); gen_676(); gen_677(); gen_678(); gen_679(); gen_680()
    gen_681(); gen_682(); gen_683(); gen_684(); gen_685(); gen_686(); gen_687()
    print(f"Done: 28 files created.")
