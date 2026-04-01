#include "UQFFPBHDarkMatterImplications.h"
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _UQFFPBHDarkMatterImplications_T_H(double M){
    return (UQFFPBHDarkMatterImplications::HBAR*UQFFPBHDarkMatterImplications::C*UQFFPBHDarkMatterImplications::C*UQFFPBHDarkMatterImplications::C)/
           (8.0*M_PI*UQFFPBHDarkMatterImplications::G*M*UQFFPBHDarkMatterImplications::K_B);
}
static inline double _UQFFPBHDarkMatterImplications_L_H(double M){
    return (UQFFPBHDarkMatterImplications::HBAR*std::pow(UQFFPBHDarkMatterImplications::C,6.0))/
           (15360.0*M_PI*UQFFPBHDarkMatterImplications::G*UQFFPBHDarkMatterImplications::G*M*M);
}
static inline double _UQFFPBHDarkMatterImplications_r_s(double M){
    return 2.0*UQFFPBHDarkMatterImplications::G*M/(UQFFPBHDarkMatterImplications::C*UQFFPBHDarkMatterImplications::C);
}
static inline double _UQFFPBHDarkMatterImplications_tau(double M){
    return 5120.0*M_PI*UQFFPBHDarkMatterImplications::G*UQFFPBHDarkMatterImplications::G*M*M*M/
           (UQFFPBHDarkMatterImplications::HBAR*std::pow(UQFFPBHDarkMatterImplications::C,4.0));
}
static inline double _UQFFPBHDarkMatterImplications_Um(double r,double t){
    double tn=t/UQFFPBHDarkMatterImplications::T_N_REF;
    return (UQFFPBHDarkMatterImplications::MU_J/r)*(1.0-std::exp(-UQFFPBHDarkMatterImplications::GAMMA*t*std::cos(M_PI*tn)));
}

UQFFPBHDarkMatterImplications::UQFFPBHDarkMatterImplications(){}
double UQFFPBHDarkMatterImplications::compute_M_crit_std() const{
    return std::cbrt(HBAR*std::pow(C,4.0)*T_AGE/(5120.0*PI*G*G));
}
double UQFFPBHDarkMatterImplications::compute_M_crit_UQFF() const{
    // tau_UQFF = tau_std * 30 => M_crit3 scales with tau => M_crit scales as tau^{1/3}
    // But evaporation slower means SAME mass survives longer
    // M_crit_UQFF^3 = M_crit_std^3 / 30 => smaller PBHs survive
    double M_std=compute_M_crit_std();
    double tau_ratio=1.0/(1.0-F_TRZ)*(RHO_UA/RHO_SCM); // ~11
    return M_std/std::cbrt(tau_ratio);
}
double UQFFPBHDarkMatterImplications::compute_f_PBH_boost() const{
    // f_PBH_UQFF/f_std = (M_crit_std/M_crit_UQFF)^2 (more PBHs in window)
    double tau_ratio=1.0/(1.0-F_TRZ)*(RHO_UA/RHO_SCM);
    return std::pow(tau_ratio,2.0/3.0);
}
double UQFFPBHDarkMatterImplications::compute_rho_PBH_UQFF(double M_f,double n_PBH) const{
    double rho_std=M_f*n_PBH;
    return rho_std*compute_f_PBH_boost();
}
bool UQFFPBHDarkMatterImplications::is_DM_viable_std(double M) const{
    return M>compute_M_crit_std();
}
bool UQFFPBHDarkMatterImplications::is_DM_viable_UQFF(double M) const{
    return M>compute_M_crit_UQFF();
}
void UQFFPBHDarkMatterImplications::scan_mass_window(double M_min,double M_max,int N,
    const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"M_kg,viable_std,viable_UQFF,tau_ratio\n";}}
    for(int i=0;i<=N;i++){
        double M=std::exp(std::log(M_min)+i*(std::log(M_max)-std::log(M_min))/N);
        double tau=_UQFFPBHDarkMatterImplications_tau(M);
        *os<<M<<","<<is_DM_viable_std(M)<<","<<is_DM_viable_UQFF(M)<<","<<compute_f_PBH_boost()<<"\n";
    }
    if(ofs.is_open())ofs.close();
}
void UQFFPBHDarkMatterImplications::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void UQFFPBHDarkMatterImplications::update_from_file(const std::string& fn){
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_UQFFPBHDARKMATTERIMPLICATIONS
int main(){
    UQFFPBHDarkMatterImplications obj;
    std::cout<<"M_crit_std="<<obj.compute_M_crit_std()<<" kg\n";
    std::cout<<"M_crit_UQFF="<<obj.compute_M_crit_UQFF()<<" kg\n";
    std::cout<<"f_PBH_boost="<<obj.compute_f_PBH_boost()<<"x\n";
    obj.scan_mass_window(1e8,1e20,100,"pbh_dm_window.csv");
    return 0;
}
#endif
