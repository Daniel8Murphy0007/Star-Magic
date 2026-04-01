#include "UQFFHawkingTemperatureModulation.h"
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _UQFFHawkingTemperatureModulation_T_H(double M){
    return (UQFFHawkingTemperatureModulation::HBAR*UQFFHawkingTemperatureModulation::C*UQFFHawkingTemperatureModulation::C*UQFFHawkingTemperatureModulation::C)/
           (8.0*M_PI*UQFFHawkingTemperatureModulation::G*M*UQFFHawkingTemperatureModulation::K_B);
}
static inline double _UQFFHawkingTemperatureModulation_L_H(double M){
    return (UQFFHawkingTemperatureModulation::HBAR*std::pow(UQFFHawkingTemperatureModulation::C,6.0))/
           (15360.0*M_PI*UQFFHawkingTemperatureModulation::G*UQFFHawkingTemperatureModulation::G*M*M);
}
static inline double _UQFFHawkingTemperatureModulation_r_s(double M){
    return 2.0*UQFFHawkingTemperatureModulation::G*M/(UQFFHawkingTemperatureModulation::C*UQFFHawkingTemperatureModulation::C);
}
static inline double _UQFFHawkingTemperatureModulation_tau(double M){
    return 5120.0*M_PI*UQFFHawkingTemperatureModulation::G*UQFFHawkingTemperatureModulation::G*M*M*M/
           (UQFFHawkingTemperatureModulation::HBAR*std::pow(UQFFHawkingTemperatureModulation::C,4.0));
}
static inline double _UQFFHawkingTemperatureModulation_Um(double r,double t){
    double tn=t/UQFFHawkingTemperatureModulation::T_N_REF;
    return (UQFFHawkingTemperatureModulation::MU_J/r)*(1.0-std::exp(-UQFFHawkingTemperatureModulation::GAMMA*t*std::cos(M_PI*tn)));
}

UQFFHawkingTemperatureModulation::UQFFHawkingTemperatureModulation(){}
double UQFFHawkingTemperatureModulation::compute_T_H(double M) const{return _UQFFHawkingTemperatureModulation_T_H(M);}
double UQFFHawkingTemperatureModulation::compute_T_UQFF(double M,double r,double t) const{
    double T_H=_UQFFHawkingTemperatureModulation_T_H(M);
    if(r<=0) r=_UQFFHawkingTemperatureModulation_r_s(M);
    double Um=_UQFFHawkingTemperatureModulation_Um(r,t);
    // T_UQFF = T_H*(1+f_TRZ)*(1-rho_SCm/rho_UA)*(1+U_m/k_B T_H)
    double fac=(1.0+F_TRZ)*(1.0-RHO_SCM/RHO_UA)*(1.0+Um/(K_B*T_H));
    double T=T_H*fac;
    for(auto& m:mods_) T+=m(M,t);
    return T;
}
double UQFFHawkingTemperatureModulation::compute_planck_spectrum(double omega,double T) const{
    if(T<=0||omega<=0) return 0.0;
    double x=HBAR*omega/(K_B*T); if(x>700) return 0.0;
    return HBAR*omega*omega*omega/(PI*PI*C*C*C*(std::exp(x)-1.0));
}
double UQFFHawkingTemperatureModulation::compute_wien_peak(double T) const{
    return HBAR*C/(2.82*K_B*T);
}
double UQFFHawkingTemperatureModulation::compute_T_modulation_factor(double M,double r,double t) const{
    double T_H=compute_T_H(M);
    double T_U=compute_T_UQFF(M,r,t);
    return T_U/std::max(T_H,1e-300);
}
void UQFFHawkingTemperatureModulation::simulate_T_vs_M(double M0,double M1,double dM,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"M_kg,T_H_K,T_UQFF_K,factor\n";}}
    for(double M=M0;M<=M1;M+=dM){
        double T_H=compute_T_H(M), T_U=compute_T_UQFF(M);
        *os<<M<<","<<T_H<<","<<T_U<<","<<T_U/std::max(T_H,1e-300)<<"\n";
    }
    if(ofs.is_open())ofs.close();
}
void UQFFHawkingTemperatureModulation::simulate_spectrum(double M,double w_min,double w_max,double dw,
    const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"omega_rad_s,n_H,n_UQFF\n";}}
    double T_H=compute_T_H(M), T_U=compute_T_UQFF(M);
    for(double w=w_min;w<=w_max;w+=dw)
        *os<<w<<","<<compute_planck_spectrum(w,T_H)<<","<<compute_planck_spectrum(w,T_U)<<"\n";
    if(ofs.is_open())ofs.close();
}
void UQFFHawkingTemperatureModulation::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void UQFFHawkingTemperatureModulation::update_from_file(const std::string& fn){
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_UQFFHAWKINGTEMPERATUREMODULATION
int main(){
    UQFFHawkingTemperatureModulation obj;
    double M=1e12;
    std::cout<<"T_H="<<obj.compute_T_H(M)<<" K\n";
    std::cout<<"T_UQFF="<<obj.compute_T_UQFF(M)<<" K\n";
    std::cout<<"Modulation factor="<<obj.compute_T_modulation_factor(M,0,1e8)<<"\n";
    obj.simulate_T_vs_M(1e10,1e15,1e12,"T_modulation.csv");
    return 0;
}
#endif
