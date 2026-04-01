"""
Session 172 — Generate all 14 remaining C++ module pairs (.h + .cpp),
14 CP4 append scripts, and 14 whitepaper .md files.
Papers: PAPER_660–673, CP4 #244–257
"""
import os, textwrap

ROOT = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
WP   = os.path.join(ROOT, "whitepapers")
os.makedirs(WP, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
#  COMMON HEADER (shared constants, boilerplate)
# ─────────────────────────────────────────────────────────────────────────────
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
    static constexpr double G    = 6.6743e-11;
    static constexpr double C    = 2.998e8;
    static constexpr double HBAR = 1.0546e-34;
    static constexpr double K_B  = 1.380649e-23;
    static constexpr double PI   = 3.14159265358979323846;
    static constexpr double RHO_UA  = 7.09e-36;   // J/m3 [UA] vacuum density
    static constexpr double RHO_SCM = 7.09e-37;   // J/m3 [SCm] vacuum density
    static constexpr double F_TRZ   = 0.1;         // time-reversal factor
    static constexpr double KAPPA   = 0.0005;      // kappa day-1
    static constexpr double SSQ     = 0.57;        // [SSq]
    static constexpr double MU_J    = 3.38e23;     // J/T magnetic string j=1
    static constexpr double GAMMA   = 5.0e-5/86400.0; // s-1
    static constexpr double T_N_REF = 1.0e8;       // s normalisation
    static constexpr double M_SUN   = 1.989e30;    // kg
"""

def write(path, content):
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  CREATED: {os.path.basename(path)}")

def common_impl_head(classname, guard):
    return f"""\
// UQFF helper: {classname}.cpp utils
static double _T_H(double M){{
    return ({classname}::HBAR*{classname}::C*{classname}::C*{classname}::C)/
           (8.0*M_PI*{classname}::G*M*{classname}::K_B);
}}
static double _L_H(double M){{
    return ({classname}::HBAR*std::pow({classname}::C,6.0))/
           (15360.0*M_PI*{classname}::G*{classname}::G*M*M);
}}
static double _r_s(double M){{
    return 2.0*{classname}::G*M/({classname}::C*{classname}::C);
}}
static double _tau_std(double M){{
    return (5120.0*M_PI*{classname}::G*{classname}::G*M*M*M)/
           ({classname}::HBAR*std::pow({classname}::C,4.0));
}}
static double _U_m(double r, double t){{
    double tn = t/{classname}::T_N_REF;
    return ({classname}::MU_J/r)*(1.0-std::exp(-{classname}::GAMMA*t*std::cos(M_PI*tn)));
}}
"""

# ─────────────────────────────────────────────────────────────────────────────
#  660 — WhiteHoleRadiationUQFF
# ─────────────────────────────────────────────────────────────────────────────
def gen_660():
    cls = "WhiteHoleRadiationUQFF"
    guard = "WHITE_HOLE_RADIATION_UQFF_H"
    h = f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF White Hole Radiation — PAPER_660 | Session 172 | April 2, 2026
 *
 * Source: grok_share_fc21e30c24b4.txt — WhiteHoleRadiation class (line 26585)
 *
 * 6-step derivation (time-reversed Hawking in UQFF):
 *  Step 1: L_WH ~ -L_H  (time-reversed expulsion, magnitude = L_H)
 *  Step 2: UQFF inversion: r_s,UQFF = r_s(1-rho_SCm/rho_UA)
 *  Step 3: f_TRZ boost: L_WH' = L_H * (1+f_TRZ)
 *  Step 4: Aether amplification: L_WH'' = L_WH' * (rho_UA/rho_SCm)  ~x10
 *  Step 5: U_m channeling: L_WH,UQFF = L_WH'' * exp(U_m / k_B T_H)
 *  Step 6: Full formula:
 *    L_WH,UQFF = (hbar c^6)/(15360 pi G^2 M^2)*(1+f_TRZ)*(rho_UA/rho_SCm)*exp(U_m/k_B T_H)
 * Numerical (Sgr A*): L_WH,UQFF ~ 3e-3 W  (vs L_H ~1e-29 W)
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    {cls}();
    double compute_L_H(double M) const;
    double compute_T_H(double M) const;
    double compute_r_s_UQFF(double M) const;
    double compute_L_WH_prime(double L_H) const;
    double compute_L_WH_double_prime(double L_WH_prime) const;
    double compute_L_WH_UQFF(double M, double r=0.0, double t=0.0) const;
    void simulate_over_M(double M_start, double M_end, double dM, const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
    mutable std::mt19937 rng_{{42u}};
}};
#endif // {guard}
"""
    cpp = f"""#include "{cls}.h"
{cls}::{cls}(){{}}

double {cls}::compute_L_H(double M) const {{ return _L_H(M); }}
double {cls}::compute_T_H(double M) const {{ return _T_H(M); }}
double {cls}::compute_r_s_UQFF(double M) const {{ return _r_s(M)*(1.0-RHO_SCM/RHO_UA); }}
double {cls}::compute_L_WH_prime(double L_H) const {{ return L_H*(1.0+F_TRZ); }}
double {cls}::compute_L_WH_double_prime(double L_WH_prime) const {{ return L_WH_prime*(RHO_UA/RHO_SCM); }}
double {cls}::compute_L_WH_UQFF(double M, double r, double t) const {{
    double L_H  = compute_L_H(M);
    double T_H  = compute_T_H(M);
    if(r<=0.0) r = _r_s(M);
    double Um   = _U_m(r, t>0?t:T_N_REF);
    double exp_arg = Um/(K_B*T_H); if(exp_arg>700) exp_arg=700;
    double L = L_H*(1.0+F_TRZ)*(RHO_UA/RHO_SCM)*std::exp(exp_arg);
    for(auto& m:mods_) L+=m(M);
    return L;
}}
void {cls}::simulate_over_M(double M0,double M1,double dM,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"M_kg,L_WH_UQFF_W\\n";}}}}
    for(double M=M0;M<=M1;M+=dM) *os<<M<<","<<compute_L_WH_UQFF(M)<<"\\n";
    if(f.is_open())f.close();
}}
void {cls}::add_mod(std::function<double(double)> mod){{mods_.push_back(std::move(mod));}}
void {cls}::update_from_file(const std::string& fn){{
    std::ifstream fi(fn); if(!fi.is_open())return; std::string l;
    while(std::getline(fi,l)){{ auto e=l.find('='); if(e==std::string::npos)continue; }}
}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} whr;
    double M=8.55e36; // Sgr A*
    std::cout<<"L_H="<<whr.compute_L_H(M)<<"\\n";
    std::cout<<"T_H="<<whr.compute_T_H(M)<<"\\n";
    std::cout<<"L_WH_UQFF="<<whr.compute_L_WH_UQFF(M)<<"\\n";
    whr.simulate_over_M(1e30,1e40,1e38,"wh_radiation_sweep.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  661 — UQFFPBHDarkMatter
# ─────────────────────────────────────────────────────────────────────────────
def gen_661():
    cls="UQFFPBHDarkMatter"
    guard="UQFF_PBH_DARK_MATTER_H"
    h=f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Primordial Black Hole Dark Matter — PAPER_661 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFPBHDarkMatterImplications (line 22222)
 * Audit: SESSION_172_AUDIT_HELPER.md — PBH Dark Matter section
 *
 * Standard Hawking evaporation: tau_std = 5120 pi G² M³ / (hbar c⁴)
 * PBH with M < 1e12 kg: tau < age of universe → evaporate → gamma rays → constrain DM
 * M ~ 1e10-1e15 g: constrained window
 *
 * UQFF elevates lifetime by ~11x (same factor as LQC critical density):
 *   Step 1: tau_std (Hawking)
 *   Step 2: tau' = tau_std / (1 - f_TRZ)          x1.11  (negentropic)
 *   Step 3: tau'' = tau' * (rho_UA/rho_SCm)        x10    (aether suppression)
 *   Step 4: tau_UQFF = tau'' * exp(U_m / k_B T_H)  x2.7   (string anchor)
 *   Total factor: ~30x (PBHs in "evaporation window" become viable DM)
 *
 *   f_PBH_UQFF = Omega_PBH / Omega_DM  (enhanced by tau factor)
 *   UQFF reopens DM window: M ~ 1e10-1e15 g now stable
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    {cls}();
    double compute_tau_standard(double M) const;
    double compute_tau_prime(double tau_std) const;
    double compute_tau_double_prime(double tau_prime) const;
    double compute_tau_UQFF(double M) const;
    double compute_f_PBH_enhancement(double M) const;
    double compute_T_H(double M) const;
    bool is_DM_candidate(double M) const;
    void simulate_lifetime_sweep(double M_start, double M_end, double dM, const std::string& out="") const;
    void add_modifier(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
}};
#endif // {guard}
"""
    cpp=f"""#include "{cls}.h"
{cls}::{cls}(){{}}
double {cls}::compute_tau_standard(double M) const {{ return _tau_std(M); }}
double {cls}::compute_T_H(double M) const {{ return _T_H(M); }}
double {cls}::compute_tau_prime(double tau_std) const {{ return tau_std/(1.0-F_TRZ); }}
double {cls}::compute_tau_double_prime(double tau_prime) const {{ return tau_prime*(RHO_UA/RHO_SCM); }}
double {cls}::compute_tau_UQFF(double M) const {{
    double tau=compute_tau_double_prime(compute_tau_prime(compute_tau_standard(M)));
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=Um/(K_B*T_H); if(ex>700)ex=700;
    tau*=std::exp(ex);
    for(auto& m:mods_) tau*=m(M);
    return tau;
}}
double {cls}::compute_f_PBH_enhancement(double M) const {{
    // Ratio tau_UQFF/tau_std ~ elevation factor applied to f_PBH
    double tau_s=compute_tau_standard(M);
    if(tau_s<1e-300) return 1.0;
    return compute_tau_UQFF(M)/tau_s;
}}
bool {cls}::is_DM_candidate(double M) const {{
    // Standard: only M > ~1e12 kg survive; UQFF: window M > ~1e11 kg
    return compute_tau_UQFF(M) > 4.35e17; // older than Hubble time
}}
void {cls}::simulate_lifetime_sweep(double M0,double M1,double dM,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"M_kg,tau_std_s,tau_UQFF_s,DM_candidate\\n";}}}}
    for(double M=M0;M<=M1;M+=dM)
        *os<<M<<","<<compute_tau_standard(M)<<","<<compute_tau_UQFF(M)<<","<<(is_DM_candidate(M)?1:0)<<"\\n";
    if(f.is_open())f.close();
}}
void {cls}::add_modifier(std::function<double(double)> mod){{mods_.push_back(std::move(mod));}}
void {cls}::update_from_file(const std::string& fn){{std::ifstream fi(fn);}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} pbh;
    double M=1e12; // Small PBH 1e12 kg
    std::cout<<"tau_std="<<pbh.compute_tau_standard(M)<<" s\\n";
    std::cout<<"tau_UQFF="<<pbh.compute_tau_UQFF(M)<<" s\\n";
    std::cout<<"DM candidate: "<<(pbh.is_DM_candidate(M)?"yes":"no")<<"\\n";
    pbh.simulate_lifetime_sweep(1e9,1e15,1e12,"pbh_dm_sweep.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  662 — UQFFHawkingDerivation
# ─────────────────────────────────────────────────────────────────────────────
def gen_662():
    cls="UQFFHawkingDerivation"
    guard="UQFF_HAWKING_DERIVATION_H"
    h=f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Hawking Radiation Derivation — PAPER_662 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFHawkingRadiationDerivation (line 25551)
 *
 * Step-by-step UQFF modification of Hawking radiation:
 *  Standard: T_H = hbar c³ / (8π G M k_B)
 *            L_H = hbar c⁶ / (15360 π G² M²)
 *            dM/dt = -L_H / c²
 *  UQFF Mods:
 *    T_UQFF = T_H * (1 + f_TRZ) * (1 - rho_SCm/rho_UA)
 *    L_UQFF = L_H * exp(-U_m / k_B T_H)     [string damping]
 *    dM/dt_UQFF = -L_UQFF / c²              [suppressed evaporation]
 *
 * Explains non-evaporating BHs; also generates τ_UQFF estimates.
 * Virtual pair mechanism: [UA] vacuum modulates pair production at horizon.
 * [SCm] Type-II B_crit~1e11 T modulates thermal state.
 * f_TRZ negentropic correction reverses pair annihilations.
 * U_m damps Boltzmann factor → suppressed emission.
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    {cls}();
    double compute_T_standard(double M) const;
    double compute_T_UQFF(double M) const;
    double compute_L_standard(double M) const;
    double compute_L_UQFF(double M) const;
    double compute_dM_dt_standard(double M) const;
    double compute_dM_dt_UQFF(double M) const;
    void simulate_evaporation(double M0, double t0, double t1, double dt,
                              const std::string& out="") const;
    void add_term(std::function<double(double,double)> term);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> terms_;
}};
#endif // {guard}
"""
    cpp=f"""#include "{cls}.h"
{cls}::{cls}(){{}}
double {cls}::compute_T_standard(double M) const {{ return _T_H(M); }}
double {cls}::compute_T_UQFF(double M) const {{
    return _T_H(M)*(1.0+F_TRZ)*(1.0-RHO_SCM/RHO_UA);
}}
double {cls}::compute_L_standard(double M) const {{ return _L_H(M); }}
double {cls}::compute_L_UQFF(double M) const {{
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=-Um/(K_B*T_H); if(ex<-700)ex=-700;
    double L=_L_H(M)*std::exp(ex);
    for(auto& t:terms_) L+=t(M,T_H);
    return L;
}}
double {cls}::compute_dM_dt_standard(double M) const {{ return -_L_H(M)/(C*C); }}
double {cls}::compute_dM_dt_UQFF(double M) const {{ return -compute_L_UQFF(M)/(C*C); }}
void {cls}::simulate_evaporation(double M0,double t0,double t1,double dt,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"t_s,M_kg,T_UQFF_K,L_UQFF_W\\n";}}}}
    double M=M0;
    for(double t=t0;t<=t1&&M>0;t+=dt){{
        double L=compute_L_UQFF(M);
        *os<<t<<","<<M<<","<<compute_T_UQFF(M)<<","<<L<<"\\n";
        M+=compute_dM_dt_UQFF(M)*dt;
        if(M<0)M=0;
    }}
    if(f.is_open())f.close();
}}
void {cls}::add_term(std::function<double(double,double)> t){{terms_.push_back(std::move(t));}}
void {cls}::update_from_file(const std::string& fn){{std::ifstream fi(fn);}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} hd;
    double M=2e30;
    std::cout<<"T_std="<<hd.compute_T_standard(M)<<" K\\n";
    std::cout<<"T_UQFF="<<hd.compute_T_UQFF(M)<<" K\\n";
    std::cout<<"L_std="<<hd.compute_L_standard(M)<<" W\\n";
    std::cout<<"L_UQFF="<<hd.compute_L_UQFF(M)<<" W\\n";
    std::cout<<"dM/dt_UQFF="<<hd.compute_dM_dt_UQFF(M)<<" kg/s\\n";
    hd.simulate_evaporation(2e30,0.0,1e67,1e57,"hawking_evap_uqff.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  663 — UQFFBlackHoleInversion
# ─────────────────────────────────────────────────────────────────────────────
def gen_663():
    cls="UQFFBlackHoleInversion"
    guard="UQFF_BLACK_HOLE_INVERSION_H"
    h=f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Black Hole Inversion Probability — PAPER_663 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFBlackHoleInversion (line 24810)
 *
 * UQFF BH interior inversion via [UA]/[SCm] gradient flip.
 * Derives full stochastic Theta_inv criterion (log-normal distribution).
 *
 * Core equations:
 *   r_s,UQFF = r_s * (1 - delta_rho)   delta_rho = rho_SCm/rho_UA ~ 0.1
 *   E_inv,UQFF = G M² / r_s,UQFF
 *   P_inv = f_TRZ * exp(-E_inv / k_B T_H)
 *   Phi_inv = (1/delta_rho) * (G M / c) * (1 + f_TRZ)   [gradient drives flux]
 *   S_Um = exp(U_m / k_B T_H)
 *   Theta_inv = P_inv * Phi_inv * S_Um
 *   If Theta_inv > 1 → inversion occurs
 *
 * Stochastic distribution:
 *   delta_rho, f_TRZ, U_m sampled from Gaussians  → Theta_inv log-normal
 *   P_invert = Prob(Theta_inv > 1)  (Monte-Carlo or analytic CDF)
 *
 * Numerical (Sgr A*): P_invert ~ 0.95 with Gaussian variability
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    {cls}(unsigned int seed=42u);
    double compute_r_s_UQFF(double M) const;
    double compute_E_inv(double M) const;
    double compute_P_inv(double M) const;
    double compute_Phi_inv(double M) const;
    double compute_Theta_inv(double M) const;
    double compute_P_invert_MC(double M, int n=5000, double sigma=0.05) const;
    void simulate_P_invert_sweep(double M0, double M1, double dM, const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    mutable std::mt19937 rng_;
    std::vector<std::function<double(double)>> mods_;
}};
#endif // {guard}
"""
    cpp=f"""#include "{cls}.h"
{cls}::{cls}(unsigned int seed):rng_(seed){{}}
double {cls}::compute_r_s_UQFF(double M) const {{ return _r_s(M)*(1.0-RHO_SCM/RHO_UA); }}
double {cls}::compute_E_inv(double M) const {{
    double rs=compute_r_s_UQFF(M); return G*M*M/rs;
}}
double {cls}::compute_P_inv(double M) const {{
    double T_H=_T_H(M); double E=compute_E_inv(M);
    double ex=-E/(K_B*T_H); if(ex<-700)ex=-700;
    return F_TRZ*std::exp(ex);
}}
double {cls}::compute_Phi_inv(double M) const {{
    double delta_rho=RHO_SCM/RHO_UA;
    return (1.0/delta_rho)*(G*M/C)*(1.0+F_TRZ);
}}
double {cls}::compute_Theta_inv(double M) const {{
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=Um/(K_B*T_H); if(ex>700)ex=700;
    double S_Um=std::exp(ex);
    double Theta=compute_P_inv(M)*compute_Phi_inv(M)*S_Um;
    for(auto& m:mods_) Theta+=m(M);
    return Theta;
}}
double {cls}::compute_P_invert_MC(double M,int n,double sigma) const{{
    std::normal_distribution<double> nd(0.0,1.0);
    int cnt=0;
    for(int i=0;i<n;i++){{
        // Perturb vacuum densities
        double pUA =RHO_UA *std::exp(sigma*nd(rng_));
        double pSCM=RHO_SCM*std::exp(sigma*nd(rng_));
        double T_H=HBAR*C*C*C/(8.0*PI*G*M*K_B);
        double rs_u=_r_s(M)*(1.0-pSCM/pUA);
        double E=G*M*M/rs_u;
        double P=0.1*std::exp(std::max(-E/(K_B*T_H),-700.0));
        double Ph=(pUA/pSCM)*(G*M/C)*(1.0+0.1);
        double r=_r_s(M); double t=T_N_REF;
        double tn=t/T_N_REF;
        double Um=(MU_J/r)*(1.0-std::exp(-GAMMA*t*std::cos(PI*tn)));
        double ex2=Um/(K_B*T_H); if(ex2>700)ex2=700;
        if(P*Ph*std::exp(ex2)>1.0) cnt++;
    }}
    return (double)cnt/n;
}}
void {cls}::simulate_P_invert_sweep(double M0,double M1,double dM,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"M_kg,Theta_inv,P_invert\\n";}}}}
    for(double M=M0;M<=M1;M+=dM)
        *os<<M<<","<<compute_Theta_inv(M)<<","<<compute_P_invert_MC(M,200,0.05)<<"\\n";
    if(f.is_open())f.close();
}}
void {cls}::add_mod(std::function<double(double)> mod){{mods_.push_back(std::move(mod));}}
void {cls}::update_from_file(const std::string& fn){{std::ifstream fi(fn);}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} inv;
    double M=8.55e36;
    std::cout<<"Theta_inv="<<inv.compute_Theta_inv(M)<<"\\n";
    std::cout<<"P_invert_MC="<<inv.compute_P_invert_MC(M,5000,0.05)<<"\\n";
    inv.simulate_P_invert_sweep(1e30,1e40,1e38,"bh_inversion_sweep.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  664 — WhiteHoleStabilityUQFF
# ─────────────────────────────────────────────────────────────────────────────
def gen_664():
    cls="WhiteHoleStabilityUQFF"
    guard="WHITE_HOLE_STABILITY_UQFF_H"
    h=f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF White Hole Stability — PAPER_664 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — WhiteHoleStabilityInUQFF (line 27097)
 *
 * 4-proof derivation of white hole stability in UQFF:
 *  WH "explodes": L_WH ~ L_H (reversed), dM/dt ~ L_WH/c² > 0
 *
 *  Proof 1 — f_TRZ Negentropic Confinement:
 *    L' = L_WH * (1 - f_TRZ)           [~10% reduction]
 *    tau' = tau_std / (1 - f_TRZ)       [x1.11]
 *
 *  Proof 2 — Aether/SCm Density Gradient:
 *    L'' = L' * |1 - rho_UA/rho_SCm|^{-1} = L' * 0.1   [x10 confinement]
 *    tau'' = tau' * |1 - rho_UA/rho_SCm|                [x10 longer]
 *
 *  Proof 3 — U_m Magnetic String Anchoring:
 *    L_UQFF = L'' * exp(-U_m / k_B |T_WH|)             [exponential suppression]
 *    tau_UQFF = tau'' * exp(U_m / k_B |T_WH|)          [x2.7 at U_m/k_B T_H=1]
 *
 *  Proof 4 — Combined:
 *    tau_UQFF = tau_std / (1-f_TRZ) * |1-rho_UA/rho_SCm| * exp(U_m/k_B|T_WH|)
 *    Factor: 1.11 * 10 * 2.718 ~ 30.2   (Sgr A*: effectively eternal)
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    {cls}();
    double compute_L_WH(double M) const;
    double compute_T_WH(double M) const;
    double compute_L_prime(double L_WH) const;
    double compute_L_double_prime(double L_prime) const;
    double compute_L_UQFF(double L_double_prime, double T_WH) const;
    double compute_tau_standard(double M) const;
    double compute_tau_UQFF(double M) const;
    double compute_stability_factor(double M) const;
    void simulate_over_mass(double M0, double M1, double dM, const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
}};
#endif // {guard}
"""
    cpp=f"""#include "{cls}.h"
{cls}::{cls}(){{}}
double {cls}::compute_L_WH(double M) const {{ return _L_H(M); }}
double {cls}::compute_T_WH(double M) const {{ return -_T_H(M); }} // negative conceptual
double {cls}::compute_L_prime(double L) const {{ return L*(1.0-F_TRZ); }}
double {cls}::compute_L_double_prime(double L2) const {{
    double ratio=std::abs(1.0-RHO_UA/RHO_SCM); // ~9
    if(ratio<1e-10) return L2;
    return L2/ratio;  // L' * 0.1 confinement
}}
double {cls}::compute_L_UQFF(double Ldp, double T_WH) const {{
    double absT=std::abs(T_WH);
    double r=1e10; // representative r
    double Um=_U_m(r,T_N_REF);
    double ex=-Um/(K_B*absT); if(ex<-700)ex=-700;
    return Ldp*std::exp(ex);
}}
double {cls}::compute_tau_standard(double M) const {{
    return _r_s(M)/C;  // instability timescale
}}
double {cls}::compute_tau_UQFF(double M) const {{
    double tau=compute_tau_standard(M);
    double T_WH=std::abs(compute_T_WH(M));
    // Proof 1
    tau/=(1.0-F_TRZ);
    // Proof 2
    double grad=std::abs(1.0-RHO_UA/RHO_SCM);
    tau*=grad;
    // Proof 3
    double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=Um/(K_B*T_WH); if(ex>700)ex=700;
    tau*=std::exp(ex);
    for(auto& m:mods_) tau*=m(M);
    return tau;
}}
double {cls}::compute_stability_factor(double M) const {{
    double tau_s=compute_tau_standard(M);
    if(tau_s<1e-300) return 1.0;
    return compute_tau_UQFF(M)/tau_s;
}}
void {cls}::simulate_over_mass(double M0,double M1,double dM,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"M_kg,tau_std_s,tau_UQFF_s,factor\\n";}}}}
    for(double M=M0;M<=M1;M+=dM)
        *os<<M<<","<<compute_tau_standard(M)<<","<<compute_tau_UQFF(M)<<","<<compute_stability_factor(M)<<"\\n";
    if(f.is_open())f.close();
}}
void {cls}::add_mod(std::function<double(double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{std::ifstream fi(fn);}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} whs;
    double M=8.55e36;
    std::cout<<"tau_std="<<whs.compute_tau_standard(M)<<" s\\n";
    std::cout<<"tau_UQFF="<<whs.compute_tau_UQFF(M)<<" s\\n";
    std::cout<<"stability_factor="<<whs.compute_stability_factor(M)<<"\\n";
    whs.simulate_over_mass(1e30,1e40,1e38,"wh_stability_sweep.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  665 — UQFFSuppressionEquationsHawking
# ─────────────────────────────────────────────────────────────────────────────
def gen_665():
    cls="UQFFSuppressionEquationsHawking"
    guard="UQFF_SUPPRESSION_EQUATIONS_HAWKING_H"
    h=f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Hawking Radiation Suppression Equations — PAPER_665 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFSuppressionEquations (line 23548)
 *
 * Derives suppression factors on Hawking radiation from all UQFF mechanisms:
 *   S_1 = (1 + f_TRZ)         negentropic modulation of T_H
 *   S_2 = (1 - rho_SCm/rho_UA) aether-superconductive damping
 *   S_3 = exp(-U_m / k_B T_H)  magnetic string exponential suppression
 *   S_total = S_1 * S_2 * S_3  (total multiplicative suppression on L)
 *
 *   L_UQFF = L_H * S_2 * S_3  (S_1 affects T, not directly L)
 *   T_UQFF = T_H * S_1 * S_2
 *   dT_H/dM = -T_H / M        standard
 *   dT_UQFF/dM = -T_UQFF / M * (S_1 * S_2)
 *
 * Sensitivity sweep: S over rho_UA/rho_SCm ratio from 2 to 20
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    {cls}();
    double compute_S1() const;
    double compute_S2() const;
    double compute_S3(double M) const;
    double compute_S_total(double M) const;
    double compute_T_UQFF(double M) const;
    double compute_L_UQFF(double M) const;
    double compute_dT_dM_standard(double M) const;
    double compute_dT_dM_UQFF(double M) const;
    void sensitivity_sweep_rho_ratio(double ratio_min, double ratio_max, double dR,
                                     const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
}};
#endif // {guard}
"""
    cpp=f"""#include "{cls}.h"
{cls}::{cls}(){{}}
double {cls}::compute_S1() const {{ return 1.0+F_TRZ; }}
double {cls}::compute_S2() const {{ return 1.0-RHO_SCM/RHO_UA; }}
double {cls}::compute_S3(double M) const {{
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=-Um/(K_B*T_H); if(ex<-700)ex=-700;
    return std::exp(ex);
}}
double {cls}::compute_S_total(double M) const {{ return compute_S1()*compute_S2()*compute_S3(M); }}
double {cls}::compute_T_UQFF(double M) const {{ return _T_H(M)*compute_S1()*compute_S2(); }}
double {cls}::compute_L_UQFF(double M) const {{
    double L=_L_H(M)*compute_S2()*compute_S3(M);
    for(auto& m:mods_) L*=m(M);
    return L;
}}
double {cls}::compute_dT_dM_standard(double M) const {{ return -_T_H(M)/M; }}
double {cls}::compute_dT_dM_UQFF(double M) const {{
    return compute_dT_dM_standard(M)*compute_S1()*compute_S2();
}}
void {cls}::sensitivity_sweep_rho_ratio(double r0,double r1,double dR,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"rho_ratio,S2,S_total\\n";}}}}
    double M=8.55e36;
    for(double R=r0;R<=r1;R+=dR){{
        double s2=1.0-1.0/R;
        double t_h=_T_H(M); double r2=_r_s(M);
        double Um=_U_m(r2,T_N_REF);
        double ex=-Um/(K_B*t_h); if(ex<-700)ex=-700;
        double s3=std::exp(ex);
        *os<<R<<","<<s2<<","<<compute_S1()*s2*s3<<"\\n";
    }}
    if(f.is_open())f.close();
}}
void {cls}::add_mod(std::function<double(double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{std::ifstream fi(fn);}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} se;
    double M=8.55e36;
    std::cout<<"S1="<<se.compute_S1()<<"\\n";
    std::cout<<"S2="<<se.compute_S2()<<"\\n";
    std::cout<<"S3="<<se.compute_S3(M)<<"\\n";
    std::cout<<"S_total="<<se.compute_S_total(M)<<"\\n";
    std::cout<<"T_UQFF="<<se.compute_T_UQFF(M)<<"\\n";
    se.sensitivity_sweep_rho_ratio(2.0,20.0,1.0,"suppression_sweep.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  666 — UQFFGWSuppression
# ─────────────────────────────────────────────────────────────────────────────
def gen_666():
    cls="UQFFGWSuppression"
    guard="UQFF_GW_SUPPRESSION_H"
    h=f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Gravitational Wave Power Suppression — PAPER_666 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFGWSuppression (line 23812)
 *
 * Standard GW power (Peters formula for circular orbit):
 *   P_GW = 32/5 * G^4/c^5 * m1^2 * m2^2 * (m1+m2) / r^5
 *
 * UQFF suppression mechanisms:
 *   S_UA  = 1 - (rho_UA / rho_critical)    [aether absorption]
 *   S_SCm = exp(-rho_SCm * r_s / (k_B T_H)) [superconductive damping]
 *   S_TRZ = (1 - f_TRZ)                    [negentropic reversal]
 *   S_Um  = exp(-U_m / (omega_GW * c))      [string impedance; omega_GW = c/r_s]
 *
 *   P_GW_UQFF = P_GW * S_UA * S_SCm * S_TRZ * S_Um
 *
 * Strain suppression: h_UQFF = h_standard * sqrt(P_GW_UQFF / P_GW)
 * Validation: GW150914 — compare UQFF strain vs LIGO observed
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double RHO_CRIT = 9.47e-27; // kg/m3 cosmological critical density
    {cls}();
    double compute_P_GW_standard(double m1, double m2, double r) const;
    double compute_S_UA() const;
    double compute_S_SCm(double M) const;
    double compute_S_TRZ() const;
    double compute_S_Um(double M) const;
    double compute_P_GW_UQFF(double m1, double m2, double r) const;
    double compute_h_suppression(double m1, double m2, double r) const;
    void simulate_P_sweep(double r_min, double r_max, double dr,
                          double m1, double m2, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp=f"""#include "{cls}.h"
{cls}::{cls}(){{}}
double {cls}::compute_P_GW_standard(double m1,double m2,double r) const{{
    return (32.0/5.0)*std::pow(G,4)/std::pow(C,5)*(m1*m1)*(m2*m2)*(m1+m2)/std::pow(r,5);
}}
double {cls}::compute_S_UA() const{{
    return 1.0-RHO_UA/RHO_CRIT;
}}
double {cls}::compute_S_SCm(double M) const{{
    double T_H=_T_H(M); double rs=_r_s(M);
    double ex=-RHO_SCM*rs/(K_B*T_H); if(ex<-700)ex=-700;
    return std::exp(ex);
}}
double {cls}::compute_S_TRZ() const{{ return 1.0-F_TRZ; }}
double {cls}::compute_S_Um(double M) const{{
    double rs=_r_s(M); double omega_GW=C/rs;
    double Um=_U_m(rs,T_N_REF);
    double ex=-Um/(omega_GW*C); if(ex<-700)ex=-700;
    return std::exp(ex);
}}
double {cls}::compute_P_GW_UQFF(double m1,double m2,double r) const{{
    double M=m1+m2;
    double P=compute_P_GW_standard(m1,m2,r);
    P*=compute_S_UA()*compute_S_SCm(M)*compute_S_TRZ()*compute_S_Um(M);
    for(auto& m:mods_) P*=m(m1+m2,r);
    return P;
}}
double {cls}::compute_h_suppression(double m1,double m2,double r) const{{
    double P_std=compute_P_GW_standard(m1,m2,r);
    double P_uqff=compute_P_GW_UQFF(m1,m2,r);
    if(P_std<1e-300) return 1.0;
    return std::sqrt(P_uqff/P_std);
}}
void {cls}::simulate_P_sweep(double r0,double r1,double dr,double m1,double m2,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"r_m,P_std_W,P_UQFF_W,h_supp\\n";}}}}
    for(double r=r0;r<=r1;r+=dr)
        *os<<r<<","<<compute_P_GW_standard(m1,m2,r)<<","<<compute_P_GW_UQFF(m1,m2,r)<<","<<compute_h_suppression(m1,m2,r)<<"\\n";
    if(f.is_open())f.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{std::ifstream fi(fn);}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} gws;
    // GW150914: m1=36 Msun, m2=29 Msun, r~350 Rsun
    double m1=36*1.989e30, m2=29*1.989e30, r=350*6.96e8;
    std::cout<<"P_GW_std="<<gws.compute_P_GW_standard(m1,m2,r)<<" W\\n";
    std::cout<<"P_GW_UQFF="<<gws.compute_P_GW_UQFF(m1,m2,r)<<" W\\n";
    std::cout<<"h_suppression="<<gws.compute_h_suppression(m1,m2,r)<<"\\n";
    gws.simulate_P_sweep(1e8,1e12,1e8,m1,m2,"gw_suppression_sweep.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  667 — UQFFBlackHoleStabilityProofs
# ─────────────────────────────────────────────────────────────────────────────
def gen_667():
    cls="UQFFBlackHoleStabilityProofs"
    guard="UQFF_BLACK_HOLE_STABILITY_PROOFS_H"
    h=f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF BH Stability Mathematical Proofs — PAPER_667 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFBlackHoleStabilityProofs (line 23305)
 *
 * 4 proofs that UQFF black holes are more stable than standard GR:
 *
 * Proof 1 — Negentropic Confinement (f_TRZ):
 *   L_UQFF' = L_H * (1 - f_TRZ)    tau' = tau / (1 - f_TRZ) x1.11
 *
 * Proof 2 — Aether-SCm Gradient Barrier:
 *   E_barrier = k_B T_H * (rho_SCm / rho_UA)   (inverted: SCm/UA confinement)
 *   tau'' = tau' * (rho_UA / rho_SCm)           x10
 *
 * Proof 3 — U_m String Anchoring:
 *   tau_UQFF = tau'' * exp(U_m / k_B T_H)       x2.718
 *
 * Proof 4 — Combined:
 *   tau_UQFF = tau / (1 - f_TRZ) * (rho_UA/rho_SCm) * exp(U_m/k_B T_H)
 *   Factor = 1.11 * 10 * 2.718 ~ 30  (Sgr A*: eternal stability)
 *
 * These proofs apply to BH stability (not WH); dual of WhiteHoleStabilityUQFF.
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    {cls}();
    double compute_tau_Hawking(double M) const;          // Standard Hawking τ
    double compute_E_barrier(double M) const;
    double compute_tau_proof1(double M) const;
    double compute_tau_proof2(double M) const;
    double compute_tau_proof3(double M) const;
    double compute_tau_UQFF(double M) const;
    double compute_stability_factor(double M) const;
    void prove_stability(double M) const;                 // Print proof chain
    void simulate_stability_sweep(double M0, double M1, double dM, const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
}};
#endif // {guard}
"""
    cpp=f"""#include "{cls}.h"
{cls}::{cls}(){{}}
double {cls}::compute_tau_Hawking(double M) const{{return _tau_std(M);}}
double {cls}::compute_E_barrier(double M) const{{return K_B*_T_H(M)*(RHO_SCM/RHO_UA);}}
double {cls}::compute_tau_proof1(double M) const{{return _tau_std(M)/(1.0-F_TRZ);}}
double {cls}::compute_tau_proof2(double M) const{{return compute_tau_proof1(M)*(RHO_UA/RHO_SCM);}}
double {cls}::compute_tau_proof3(double M) const{{
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=Um/(K_B*T_H); if(ex>700)ex=700;
    return compute_tau_proof2(M)*std::exp(ex);
}}
double {cls}::compute_tau_UQFF(double M) const{{
    double tau=compute_tau_proof3(M);
    for(auto& m:mods_) tau*=m(M);
    return tau;
}}
double {cls}::compute_stability_factor(double M) const{{
    double ts=_tau_std(M); if(ts<1e-300) return 1.0;
    return compute_tau_UQFF(M)/ts;
}}
void {cls}::prove_stability(double M) const{{
    std::cout<<"=== UQFF BH Stability Proofs for M="<<M<<" kg ===\\n"
             <<"  tau_Hawking = "<<compute_tau_Hawking(M)<<" s\\n"
             <<"  Proof 1 (f_TRZ): tau' = "<<compute_tau_proof1(M)<<" s (x"<<1.0/(1.0-F_TRZ)<<")\\n"
             <<"  Proof 2 (rho):   tau'' = "<<compute_tau_proof2(M)<<" s (x"<<RHO_UA/RHO_SCM<<")\\n"
             <<"  Proof 3 (U_m):   tau_UQFF = "<<compute_tau_proof3(M)<<" s\\n"
             <<"  Proof 4 factor = "<<compute_stability_factor(M)<<"\\n";
}}
void {cls}::simulate_stability_sweep(double M0,double M1,double dM,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"M_kg,tau_Hawking_s,tau_UQFF_s,factor\\n";}}}}
    for(double M=M0;M<=M1;M+=dM)
        *os<<M<<","<<compute_tau_Hawking(M)<<","<<compute_tau_UQFF(M)<<","<<compute_stability_factor(M)<<"\\n";
    if(f.is_open())f.close();
}}
void {cls}::add_mod(std::function<double(double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{std::ifstream fi(fn);}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} sp;
    double M=8.55e36;
    sp.prove_stability(M);
    sp.simulate_stability_sweep(1e20,1e40,1e38,"bh_stability_proofs.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  668 — UQFFStabilityPrimordialBH
# ─────────────────────────────────────────────────────────────────────────────
def gen_668():
    cls="UQFFStabilityPrimordialBH"
    guard="UQFF_STABILITY_PRIMORDIAL_BH_H"
    h=f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Primordial BH Stability Analysis — PAPER_668 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFStabilityForPrimordialBHs (line 25798)
 *
 * PBH stability analysis in UQFF context:
 *  Standard: τ_std = 5120 π G² M³ / (ħ c⁴)
 *  M < 1e12 kg → τ < universe age → evaporate
 *  M ~ 1e10-1e15 g constrained by gamma-ray observations
 *
 * UQFF step-by-step:
 *   1. τ_std (Hawking)
 *   2. τ' = τ / (1 - f_TRZ)      x1.11  negentropic
 *   3. τ'' = τ' * (rho_UA/rho_SCm)  x10
 *   4. τ_UQFF = τ'' * exp(U_m/k_B T_H) x2.7
 *   Total: ~30x
 *
 * Numerical check M=1e12 kg:
 *   τ_std ~ 1e10 yr; τ_UQFF ~ 3e11 yr > universe age → UQFF promotes DM candidate
 *
 * Mass categories: stable (τ_UQFF > t_H), marginal (0.1 t_H < τ < t_H), evaporating
 * t_H = 4.35e17 s (Hubble time)
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double T_HUBBLE = 4.35e17; // s
    {cls}();
    double compute_tau_std(double M) const;
    double compute_tau_UQFF(double M) const;
    std::string classify(double M) const;    // stable/marginal/evaporating
    double pbh_dm_window_min_mass_uqff() const;  // minimum DM-viable mass
    void simulate_mass_stability_map(double M0, double M1, int n_pts, const std::string& out="") const;
    void add_modifier(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
}};
#endif // {guard}
"""
    cpp=f"""#include "{cls}.h"
{cls}::{cls}(){{}}
double {cls}::compute_tau_std(double M) const{{return _tau_std(M);}}
double {cls}::compute_tau_UQFF(double M) const{{
    double tau=_tau_std(M)/(1.0-F_TRZ);
    tau*=(RHO_UA/RHO_SCM);
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=Um/(K_B*T_H); if(ex>700)ex=700;
    tau*=std::exp(ex);
    for(auto& m:mods_) tau*=m(M);
    return tau;
}}
std::string {cls}::classify(double M) const{{
    double tu=compute_tau_UQFF(M);
    if(tu>=T_HUBBLE) return "stable_DM";
    if(tu>=0.1*T_HUBBLE) return "marginal";
    return "evaporating";
}}
double {cls}::pbh_dm_window_min_mass_uqff() const{{
    // Binary search for minimum M such that tau_UQFF >= T_HUBBLE
    double lo=1e6, hi=1e25;
    for(int i=0;i<80;i++){{
        double mid=(lo+hi)/2.0;
        if(compute_tau_UQFF(mid)>=T_HUBBLE) hi=mid; else lo=mid;
    }}
    return (lo+hi)/2.0;
}}
void {cls}::simulate_mass_stability_map(double M0,double M1,int n,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"M_kg,tau_std_s,tau_UQFF_s,class\\n";}}}}
    double dM=(M1-M0)/std::max(n-1,1);
    for(int i=0;i<n;i++){{
        double M=M0+i*dM;
        *os<<M<<","<<compute_tau_std(M)<<","<<compute_tau_UQFF(M)<<","<<classify(M)<<"\\n";
    }}
    if(f.is_open())f.close();
}}
void {cls}::add_modifier(std::function<double(double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{std::ifstream fi(fn);}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} pbhs;
    double M=1e12;
    std::cout<<"tau_std="<<pbhs.compute_tau_std(M)<<" s\\n";
    std::cout<<"tau_UQFF="<<pbhs.compute_tau_UQFF(M)<<" s\\n";
    std::cout<<"class="<<pbhs.classify(M)<<"\\n";
    std::cout<<"Min DM mass UQFF="<<pbhs.pbh_dm_window_min_mass_uqff()<<" kg\\n";
    pbhs.simulate_mass_stability_map(1e8,1e18,100,"pbh_stability_map.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  669 — UQFFComparedToGW150914
# ─────────────────────────────────────────────────────────────────────────────
def gen_669():
    cls="UQFFComparedToGW150914"
    guard="UQFF_COMPARED_TO_GW150914_H"
    h=f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF Waveform vs LIGO GW150914 — PAPER_669 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFComparedToGW150914 (line 25030)
 *
 * GW150914: First direct GW detection (Sep 14, 2015, LIGO)
 *   m1 = 36 M☉, m2 = 29 M☉, M_final = 62 M☉, radiated = 3 M☉ c²
 *   Peak frequency: f_peak ~ 150 Hz, strain: h ~ 1e-21
 *   Distance: d ~ 410 Mpc
 *
 * UQFF modifications to strain:
 *   h_UQFF(t) = h_GR(t) * (1 - f_TRZ) * S_SCm * exp(-U_m * omega / c²)
 *   omega = 2 pi f_GW  (instantaneous angular frequency)
 *   S_SCm = exp(-rho_SCm * lambda_GW)   lambda_GW = c/f_GW
 *   Phase shift: dphi_UQFF = dphi_GR + kappa * f_TRZ * t
 *
 * Chirp mass: M_chirp = (m1*m2)^(3/5) / (m1+m2)^(1/5)
 * Inspiral frequency evolution: df/dt = 96/5 * pi^(8/3) * (G M_chirp)^(5/3) / c^5 * f^(11/3)
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    // GW150914 canonical parameters
    static constexpr double M1_GW150914 = 36.0*1.989e30;
    static constexpr double M2_GW150914 = 29.0*1.989e30;
    static constexpr double D_GW150914  = 410.0*3.086e22; // 410 Mpc in m
    static constexpr double F_PEAK      = 150.0;          // Hz
    {cls}();
    double compute_chirp_mass(double m1, double m2) const;
    double compute_h_GR(double Mc, double d, double f) const;
    double compute_h_UQFF(double Mc, double d, double f) const;
    double compute_df_dt(double Mc, double f) const;
    double compute_dphi_UQFF(double f, double t) const;
    double compute_S_SCm(double f) const;
    void simulate_inspiral(double f_start, double f_end, double df,
                           double m1, double m2, double d, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp=f"""#include "{cls}.h"
{cls}::{cls}(){{}}
double {cls}::compute_chirp_mass(double m1,double m2) const{{
    return std::pow(m1*m2,3.0/5.0)/std::pow(m1+m2,1.0/5.0);
}}
double {cls}::compute_h_GR(double Mc,double d,double f) const{{
    // h ~ 4/d * (G Mc/c²)^(5/3) * (pi f)^(2/3) / c  (leading order)
    double gmc=G*Mc/(C*C);
    return (4.0/d)*std::pow(gmc,5.0/3.0)*std::pow(PI*f,2.0/3.0)/C;
}}
double {cls}::compute_S_SCm(double f) const{{
    double lam=C/f;
    double ex=-RHO_SCM*lam; if(ex<-700)ex=-700;
    return std::exp(ex);
}}
double {cls}::compute_h_UQFF(double Mc,double d,double f) const{{
    double h=compute_h_GR(Mc,d,f);
    double omega=2.0*PI*f;
    double Um=MU_J/(_r_s(2*Mc));
    double ex=-Um*omega/(C*C); if(ex<-700)ex=-700;
    h*=(1.0-F_TRZ)*compute_S_SCm(f)*std::exp(ex);
    for(auto& m:mods_) h*=m(Mc,f);
    return h;
}}
double {cls}::compute_df_dt(double Mc,double f) const{{
    return (96.0/5.0)*std::pow(PI,8.0/3.0)*std::pow(G*Mc,5.0/3.0)/std::pow(C,5.0)*std::pow(f,11.0/3.0);
}}
double {cls}::compute_dphi_UQFF(double f,double t) const{{
    return 2.0*PI*f*t + KAPPA*F_TRZ*t;
}}
void {cls}::simulate_inspiral(double f0,double f1,double df,double m1,double m2,double d,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"f_Hz,h_GR,h_UQFF\\n";}}}}
    double Mc=compute_chirp_mass(m1,m2);
    for(double freq=f0;freq<=f1;freq+=df)
        *os<<freq<<","<<compute_h_GR(Mc,d,freq)<<","<<compute_h_UQFF(Mc,d,freq)<<"\\n";
    if(f.is_open())f.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{std::ifstream fi(fn);}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} gw;
    double Mc=gw.compute_chirp_mass(M1_GW150914,M2_GW150914);
    std::cout<<"Chirp mass="<<Mc/1.989e30<<" Msun\\n";
    std::cout<<"h_GR(150Hz)="<<gw.compute_h_GR(Mc,D_GW150914,150.0)<<"\\n";
    std::cout<<"h_UQFF(150Hz)="<<gw.compute_h_UQFF(Mc,D_GW150914,150.0)<<"\\n";
    gw.simulate_inspiral(10.0,300.0,1.0,M1_GW150914,M2_GW150914,D_GW150914,"gw150914_uqff.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  670 — UQFFBlackHoleAccretionModel
# ─────────────────────────────────────────────────────────────────────────────
def gen_670():
    cls="UQFFBlackHoleAccretionModel"
    guard="UQFF_BLACK_HOLE_ACCRETION_MODEL_H"
    h=f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF BH Accretion Model — PAPER_670 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFBlackHoleAccretionModel (line 25268)
 *
 * Bondi accretion in UQFF:
 *   Mdot_Bondi = 4 pi lambda_B (G M)² rho_inf / cs³
 *   lambda_B = 1/4 (polytropic index gamma_ad=5/3)
 *   cs = sound speed = sqrt(gamma_ad k_B T_inf / m_p)  [proton mass m_p]
 *
 * UQFF modifications:
 *   rho_eff = rho_inf + rho_UA - rho_SCm  (aether mass contribution)
 *   S_TRZ   = (1 + f_TRZ)                 (negentropic inflow boost)
 *   S_Um    = 1 - exp(-U_m / k_B T_inf)   (magnetic string impedance)
 *   Mdot_UQFF = Mdot_Bondi * (rho_eff/rho_inf) * S_TRZ * S_Um
 *
 * Eddington accretion:
 *   L_Edd = 4 pi G M m_p c / sigma_T   [sigma_T Thompson cross-section]
 *   Mdot_Edd = L_Edd / (eta c²)         [eta=0.1 efficiency]
 *   Mdot_UQFF vs. Mdot_Edd ratio
 *
 * Simulation: M(t) evolution from dM/dt = Mdot_UQFF - Mdot_Hawking
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double M_P = 1.6726e-27;   // kg proton mass
    static constexpr double SIGMA_T = 6.652e-29; // m2 Thompson
    static constexpr double ETA = 0.1;           // radiative efficiency
    {cls}();
    double compute_cs(double T_inf) const;
    double compute_Mdot_Bondi(double M, double rho_inf, double T_inf) const;
    double compute_Mdot_UQFF(double M, double rho_inf, double T_inf) const;
    double compute_L_Edd(double M) const;
    double compute_Mdot_Edd(double M) const;
    double compute_Eddington_ratio(double M, double rho_inf, double T_inf) const;
    void simulate_M_evolution(double M0, double t0, double t1, double dt,
                              double rho_inf, double T_inf, const std::string& out="") const;
    void add_mod(std::function<double(double,double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double,double)>> mods_;
}};
#endif // {guard}
"""
    cpp=f"""#include "{cls}.h"
{cls}::{cls}(){{}}
double {cls}::compute_cs(double T_inf) const{{
    return std::sqrt(5.0/3.0*K_B*T_inf/M_P);
}}
double {cls}::compute_Mdot_Bondi(double M,double rho_inf,double T_inf) const{{
    double cs=compute_cs(T_inf);
    return 4.0*PI*0.25*G*G*M*M*rho_inf/(cs*cs*cs);
}}
double {cls}::compute_Mdot_UQFF(double M,double rho_inf,double T_inf) const{{
    double rho_eff=rho_inf+RHO_UA-RHO_SCM;
    double s_trz=1.0+F_TRZ;
    double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=-Um/(K_B*T_inf); if(ex<-700)ex=-700;
    double s_um=1.0-std::exp(ex);
    double Mdot=compute_Mdot_Bondi(M,rho_inf,T_inf)*(rho_eff/rho_inf)*s_trz*s_um;
    for(auto& m:mods_) Mdot*=m(M,T_inf);
    return Mdot;
}}
double {cls}::compute_L_Edd(double M) const{{
    return 4.0*PI*G*M*M_P*C/SIGMA_T;
}}
double {cls}::compute_Mdot_Edd(double M) const{{
    return compute_L_Edd(M)/(ETA*C*C);
}}
double {cls}::compute_Eddington_ratio(double M,double rho_inf,double T_inf) const{{
    double Md_Edd=compute_Mdot_Edd(M);
    if(Md_Edd<1e-300) return 0.0;
    return compute_Mdot_UQFF(M,rho_inf,T_inf)/Md_Edd;
}}
void {cls}::simulate_M_evolution(double M0,double t0,double t1,double dt,double rho_inf,double T_inf,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"t_s,M_kg,Mdot_UQFF,Mdot_Edd\\n";}}}}
    double M=M0;
    for(double t=t0;t<=t1;t+=dt){{
        double Md=compute_Mdot_UQFF(M,rho_inf,T_inf);
        double Ml=-_L_H(M)/(C*C);
        *os<<t<<","<<M<<","<<Md<<","<<compute_Mdot_Edd(M)<<"\\n";
        M+=(Md+Ml)*dt;
        if(M<0) M=0;
    }}
    if(f.is_open())f.close();
}}
void {cls}::add_mod(std::function<double(double,double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{std::ifstream fi(fn);}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} acc;
    double M=8.55e36, rho=1e-19, T=1e7; // typical SMBH environment
    std::cout<<"Mdot_Bondi="<<acc.compute_Mdot_Bondi(M,rho,T)<<" kg/s\\n";
    std::cout<<"Mdot_UQFF="<<acc.compute_Mdot_UQFF(M,rho,T)<<" kg/s\\n";
    std::cout<<"Edd_ratio="<<acc.compute_Eddington_ratio(M,rho,T)<<"\\n";
    acc.simulate_M_evolution(8.55e36,0,3.15e13,3.15e10,rho,T,"accretion_evolution.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  671 — UQFFDMDtDerivation
# ─────────────────────────────────────────────────────────────────────────────
def gen_671():
    cls="UQFFDMDtDerivation"
    guard="UQFF_DM_DT_DERIVATION_H"
    h=f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF dM/dt Full Derivation — PAPER_671 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFDMDtFromUQFFSteps (line 26116)
 *
 * Full step-by-step derivation of dM/dt from UQFF mechanics:
 *
 * Standard Hawking:
 *   dM/dt = -L_H / c²  = -hbar c⁴ / (15360 pi G² M² c²)
 *
 * UQFF Step 1 — f_TRZ negentropic reversal:
 *   (dM/dt)' = dM/dt * (1 - f_TRZ)         [suppresses evaporation ~10%]
 *
 * UQFF Step 2 — Aether density quench:
 *   (dM/dt)'' = (dM/dt)' * (rho_SCm / rho_UA)  [x0.1 further suppression]
 *
 * UQFF Step 3 — U_m string anchor:
 *   (dM/dt)_UQFF = (dM/dt)'' * exp(-U_m / k_B T_H)  [exponential suppression]
 *
 * Combined:
 *   (dM/dt)_UQFF = dM/dt * (1-f_TRZ) * (rho_SCm/rho_UA) * exp(-U_m/k_B T_H)
 *   Suppression factor: 0.9 * 0.1 * e^(-1) ~ 0.033  (30x slower evaporation)
 *
 * Time-integrated:
 *   M(t) via Euler solver; also analytic approximation for near-constant suppressor:
 *   M(t) ≈ (M0³ - 3 * |dM/dt|_UQFF / c² * t * c²)^(1/3)
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    {cls}();
    double compute_dM_dt_standard(double M) const;
    double compute_dM_dt_step1(double M) const;
    double compute_dM_dt_step2(double M) const;
    double compute_dM_dt_UQFF(double M) const;
    double compute_suppression_factor(double M) const;
    double compute_M_at_t(double M0, double t) const;  // analytic approx
    void simulate_evaporation(double M0, double t0, double t1, double dt, const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
}};
#endif // {guard}
"""
    cpp=f"""#include "{cls}.h"
{cls}::{cls}(){{}}
double {cls}::compute_dM_dt_standard(double M) const{{return -_L_H(M)/(C*C);}}
double {cls}::compute_dM_dt_step1(double M) const{{return compute_dM_dt_standard(M)*(1.0-F_TRZ);}}
double {cls}::compute_dM_dt_step2(double M) const{{return compute_dM_dt_step1(M)*(RHO_SCM/RHO_UA);}}
double {cls}::compute_dM_dt_UQFF(double M) const{{
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=-Um/(K_B*T_H); if(ex<-700)ex=-700;
    double dMdt=compute_dM_dt_step2(M)*std::exp(ex);
    for(auto& m:mods_) dMdt*=m(M);
    return dMdt;
}}
double {cls}::compute_suppression_factor(double M) const{{
    double std_rate=compute_dM_dt_standard(M);
    if(std::abs(std_rate)<1e-300) return 1.0;
    return compute_dM_dt_UQFF(M)/std_rate;
}}
double {cls}::compute_M_at_t(double M0,double t) const{{
    // Analytic: dM/dt ~ -A/M² → M(t) = (M0³ - 3 A t)^(1/3)
    // A = |dM/dt| * M² evaluated at M0 (constant suppression approx)
    double rate=std::abs(compute_dM_dt_UQFF(M0));
    double A=rate*M0*M0;
    double val=M0*M0*M0-3.0*A*t;
    if(val<0) return 0.0;
    return std::cbrt(val);
}}
void {cls}::simulate_evaporation(double M0,double t0,double t1,double dt,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"t_s,M_kg,dM_dt_UQFF,suppression\\n";}}}}
    double M=M0;
    for(double t=t0;t<=t1&&M>0;t+=dt){{
        double dMdt=compute_dM_dt_UQFF(M);
        *os<<t<<","<<M<<","<<dMdt<<","<<compute_suppression_factor(M)<<"\\n";
        M+=dMdt*dt; if(M<0)M=0;
    }}
    if(f.is_open())f.close();
}}
void {cls}::add_mod(std::function<double(double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{std::ifstream fi(fn);}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} dm;
    double M=2e30;
    std::cout<<"dM/dt_std="<<dm.compute_dM_dt_standard(M)<<" kg/s\\n";
    std::cout<<"dM/dt_UQFF="<<dm.compute_dM_dt_UQFF(M)<<" kg/s\\n";
    std::cout<<"suppression="<<dm.compute_suppression_factor(M)<<"\\n";
    std::cout<<"M at 1e67 s (analytic)="<<dm.compute_M_at_t(2e30,1e67)<<" kg\\n";
    dm.simulate_evaporation(2e30,0,1e67,1e57,"dmdt_evaporation.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  672 — UQFFEvaporationTimescale
# ─────────────────────────────────────────────────────────────────────────────
def gen_672():
    cls="UQFFEvaporationTimescale"
    guard="UQFF_EVAPORATION_TIMESCALE_H"
    h=f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF BH Evaporation Timescale — PAPER_672 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — UQFFEvaporationTimescale (line 26358)
 *
 * Standard evaporation timescale:
 *   tau_Hawking = 5120 pi G² M³ / (hbar c⁴)
 *   For M = M☉: tau ~ 2.1e74 years (much > universe age)
 *   For M = 1e11 kg:  tau ~ 3.4e-7 s  (instant)
 *   Boundary mass M_evap (tau = t_Hubble): M_evap ~ 1e11 kg
 *
 * UQFF timescale:
 *   tau_UQFF = tau_Hawking * (1/(1-f_TRZ)) * (rho_UA/rho_SCm) * exp(U_m/k_B T_H)
 *   Factor ~ 30 (same as stability proofs)
 *   UQFF boundary mass: M_evap,UQFF = M_evap * (1/factor)^(1/3) ~ M_evap / 3.1
 *
 * Universe-age crossings:
 *   Standard: M_cross,std  = (hbar c⁴ t_H / 5120 pi G²)^(1/3)
 *   UQFF:     M_cross,UQFF = M_cross,std / (factor)^(1/3)
 *
 * Sensitivity to U_m: tau_UQFF / tau_std vs U_m/k_B T_H
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double T_HUBBLE = 4.35e17; // s
    {cls}();
    double compute_tau_standard(double M) const;
    double compute_tau_UQFF(double M) const;
    double compute_factor(double M) const;
    double compute_M_cross_standard() const;
    double compute_M_cross_UQFF() const;
    void sensitivity_Um_sweep(double Um_min, double Um_max, double dUm,
                              double M, const std::string& out="") const;
    void simulate_timescale_sweep(double M0, double M1, int n_pts, const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
}};
#endif // {guard}
"""
    cpp=f"""#include "{cls}.h"
{cls}::{cls}(){{}}
double {cls}::compute_tau_standard(double M) const{{return _tau_std(M);}}
double {cls}::compute_tau_UQFF(double M) const{{
    double tau=_tau_std(M)/(1.0-F_TRZ)*(RHO_UA/RHO_SCM);
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=Um/(K_B*T_H); if(ex>700)ex=700;
    tau*=std::exp(ex);
    for(auto& m:mods_) tau*=m(M);
    return tau;
}}
double {cls}::compute_factor(double M) const{{
    double ts=compute_tau_standard(M);
    if(ts<1e-300) return 1.0;
    return compute_tau_UQFF(M)/ts;
}}
double {cls}::compute_M_cross_standard() const{{
    return std::cbrt(HBAR*std::pow(C,4)*T_HUBBLE/(5120.0*PI*G*G));
}}
double {cls}::compute_M_cross_UQFF() const{{
    // tau_UQFF ~ 30 * tau_std; tau_std(M) = 5120 pi G² M³ / (hbar c^4)
    // tau_UQFF = t_H → M^3 = t_H * hbar c^4 / (5120 pi G² * 30)
    double factor=1.0/(1.0-F_TRZ)*(RHO_UA/RHO_SCM)*std::exp(1.0); // ~30 approx
    return std::cbrt(HBAR*std::pow(C,4)*T_HUBBLE/(5120.0*PI*G*G*factor));
}}
void {cls}::sensitivity_Um_sweep(double Um0,double Um1,double dUm,double M,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"Um_ratio,tau_UQFF_s,factor\\n";}}}}
    double T_H=_T_H(M);
    for(double Um_r=Um0;Um_r<=Um1;Um_r+=dUm){{
        double Um_val=Um_r*K_B*T_H;
        double tau=_tau_std(M)/(1.0-F_TRZ)*(RHO_UA/RHO_SCM)*std::exp(std::min(Um_r,700.0));
        *os<<Um_r<<","<<tau<<","<<tau/_tau_std(M)<<"\\n";
    }}
    if(f.is_open())f.close();
}}
void {cls}::simulate_timescale_sweep(double M0,double M1,int n,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"M_kg,tau_std_s,tau_UQFF_s,factor\\n";}}}}
    double dM=(M1-M0)/std::max(n-1,1);
    for(int i=0;i<n;i++){{
        double M=M0+i*dM;
        *os<<M<<","<<compute_tau_standard(M)<<","<<compute_tau_UQFF(M)<<","<<compute_factor(M)<<"\\n";
    }}
    if(f.is_open())f.close();
}}
void {cls}::add_mod(std::function<double(double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{std::ifstream fi(fn);}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} evap;
    std::cout<<"M_cross_std="<<evap.compute_M_cross_standard()<<" kg\\n";
    std::cout<<"M_cross_UQFF="<<evap.compute_M_cross_UQFF()<<" kg\\n";
    double M=1e12;
    std::cout<<"tau_std(1e12kg)="<<evap.compute_tau_standard(M)<<" s\\n";
    std::cout<<"tau_UQFF(1e12kg)="<<evap.compute_tau_UQFF(M)<<" s\\n";
    evap.sensitivity_Um_sweep(0.0,5.0,0.1,M,"evap_Um_sensitivity.csv");
    evap.simulate_timescale_sweep(1e8,1e20,100,"evap_timescale_sweep.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  673 — UQFFAdvancementsAndTHzHoles
# ─────────────────────────────────────────────────────────────────────────────
def gen_673():
    cls="UQFFAdvancementsAndTHzHoles"
    guard="UQFF_ADVANCEMENTS_AND_THZ_HOLES_H"
    h=f"""#ifndef {guard}
#define {guard}
{COMMON_INCLUDES}
/**
 * @class {cls}
 * @brief UQFF THz Holes + Red Dwarf Reactor Meta-Module — PAPER_673 | Session 172
 *
 * Source: grok_share_fc21e30c24b4.txt — LearningAndFrameworkAdvancement (line 28923)
 *
 * Meta-module synthesising UQFF advancements discovered in grok_share_fc21e30c24b4.txt:
 *
 * 1. THz Hole Analogy:
 *    THz "holes" = quasi-particles in superconductors that mimic BH horizons at THz freq.
 *    f_THz = k_B T_c / (2 pi hbar)  [critical temperature THz frequency]
 *    T_c ~ 100 K → f_THz ~ 2 THz
 *    UQFF analogy: rho_SCm / (m_e c²) ~ electron pair density at THz scale
 *    L_THz_UQFF = L_H * (f_THz / f_Hawking)^4 * (rho_UA / rho_SCm)
 *
 * 2. Red Dwarf Reactor UQFF:
 *    Red dwarf core: T_core ~ 1e7 K, rho_core ~ 1e5 kg/m³
 *    UQFF fusion suppression: Gamma_UQFF = sigma_fus * n² * (1 - rho_SCm/rho_UA)
 *    Lifetime extension: tau_RD_UQFF = tau_std * (rho_UA/rho_SCm) * (1+f_TRZ)
 *    tau_std ~ 1e13 yr; tau_RD_UQFF ~ 1.1e14 yr
 *
 * 3. Framework Advancement Score (FAS):
 *    FAS = N_papers * (1 + f_TRZ) * sqrt(rho_UA/rho_SCm)
 *    Tracks learning advancement across UQFF sessions
 *
 * 4. Self-Consistent Cycle:
 *    KB7 → BH_Bounce → BH_WH_Transition → WH_Stability → THz_Holes → back to KB7
 */
class {cls} {{
public:
{COMMON_CONSTANTS}
    static constexpr double M_E = 9.109e-31;      // kg electron mass
    static constexpr double T_C_SC = 100.0;        // K superconductor critical temp
    {cls}();
    // THz hole section
    double compute_f_THz(double T_c=T_C_SC) const;
    double compute_f_Hawking_M(double M) const;
    double compute_L_THz_UQFF(double M, double T_c=T_C_SC) const;
    // Red dwarf section
    double compute_tau_RD_standard() const;
    double compute_tau_RD_UQFF() const;
    double compute_fusion_suppression(double rho_core=1e5, double T_core=1e7) const;
    // Framework tracking
    double compute_FAS(int n_papers=673) const;
    // Cycle cross-reference
    void print_self_consistent_cycle() const;
    // Sweep
    void simulate_THz_sweep(double M0, double M1, double dM, double T_c, const std::string& out="") const;
    void add_mod(std::function<double(double)> mod);
    void update_from_file(const std::string& f);
private:
    std::vector<std::function<double(double)>> mods_;
}};
#endif // {guard}
"""
    cpp=f"""#include "{cls}.h"
{cls}::{cls}(){{}}
double {cls}::compute_f_THz(double T_c) const{{
    return K_B*T_c/(2.0*PI*HBAR);
}}
double {cls}::compute_f_Hawking_M(double M) const{{
    // f_H = k_B T_H / (2 pi hbar)
    return K_B*_T_H(M)/(2.0*PI*HBAR);
}}
double {cls}::compute_L_THz_UQFF(double M,double T_c) const{{
    double f_THz_val=compute_f_THz(T_c);
    double f_H=compute_f_Hawking_M(M);
    if(f_H<1e-300) f_H=1e-300;
    double ratio_f4=std::pow(f_THz_val/f_H,4.0); if(ratio_f4>1e200) ratio_f4=1e200;
    double L=_L_H(M)*ratio_f4*(RHO_UA/RHO_SCM);
    for(auto& m:mods_) L*=m(M);
    return L;
}}
double {cls}::compute_tau_RD_standard() const{{
    // tau_RD_std ~ 1e13 yr in seconds
    return 1e13*3.156e7;
}}
double {cls}::compute_tau_RD_UQFF() const{{
    return compute_tau_RD_standard()*(RHO_UA/RHO_SCM)*(1.0+F_TRZ);
}}
double {cls}::compute_fusion_suppression(double rho_core, double T_core) const{{
    // Gamma_UQFF = sigma_fus * n^2 * (1 - rho_SCm/rho_UA)
    // sigma_fus ~ 1e-34 m^2 at T_core=1e7 K (pp chain)
    double sigma_fus=1e-34;
    double n=rho_core/M_E; // rough number density
    return sigma_fus*n*n*(1.0-RHO_SCM/RHO_UA);
}}
double {cls}::compute_FAS(int n_papers) const{{
    return n_papers*(1.0+F_TRZ)*std::sqrt(RHO_UA/RHO_SCM);
}}
void {cls}::print_self_consistent_cycle() const{{
    std::cout<<"UQFF Self-Consistent Cycle (PAPER_657-673):\\n"
             <<"  KB7 (657) -> BH_Bounce (658) -> BH_WH_Trans (659) -> WH_Rad (660)\\n"
             <<"  -> PBH_DM (661) -> Hawking (662) -> BH_Inv (663) -> WH_Stab (664)\\n"
             <<"  -> Suppression (665) -> GW_Supp (666) -> BH_Stab (667) -> PBH_Stab (668)\\n"
             <<"  -> GW150914 (669) -> Accretion (670) -> dM/dt (671) -> tau_evap (672)\\n"
             <<"  -> THz_Holes (673) -> back to KB7\\n"
             <<"  FAS = "<<compute_FAS(673)<<"\\n";
}}
void {cls}::simulate_THz_sweep(double M0,double M1,double dM,double T_c,const std::string& out) const{{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){{f.open(out);if(f.is_open()){{os=&f;*os<<"M_kg,f_Hawking_Hz,L_THz_UQFF_W\\n";}}}}
    for(double M=M0;M<=M1;M+=dM)
        *os<<M<<","<<compute_f_Hawking_M(M)<<","<<compute_L_THz_UQFF(M,T_c)<<"\\n";
    if(f.is_open())f.close();
}}
void {cls}::add_mod(std::function<double(double)> m){{mods_.push_back(std::move(m));}}
void {cls}::update_from_file(const std::string& fn){{std::ifstream fi(fn);}}
#ifdef STANDALONE_{cls.upper()}
int main(){{
    {cls} thz;
    std::cout<<"f_THz(100K)="<<thz.compute_f_THz()<<" Hz\\n";
    std::cout<<"tau_RD_std="<<thz.compute_tau_RD_standard()<<"s\\n";
    std::cout<<"tau_RD_UQFF="<<thz.compute_tau_RD_UQFF()<<"s\\n";
    std::cout<<"FAS="<<thz.compute_FAS(673)<<"\\n";
    thz.print_self_consistent_cycle();
    thz.simulate_THz_sweep(1e15,1e42,1e40,100.0,"thz_holes_sweep.csv");
    return 0;
}}
#endif
"""
    write(os.path.join(ROOT,f"{cls}.h"), h+common_impl_head(cls,guard))
    write(os.path.join(ROOT,f"{cls}.cpp"), cpp)

# ─────────────────────────────────────────────────────────────────────────────
#  RUN ALL GENERATORS
# ─────────────────────────────────────────────────────────────────────────────
os.chdir(ROOT)
print("Generating C++ modules...")
gen_660(); gen_661(); gen_662(); gen_663(); gen_664()
gen_665(); gen_666(); gen_667(); gen_668(); gen_669()
gen_670(); gen_671(); gen_672(); gen_673()
print("Done: 28 C++ files created")
