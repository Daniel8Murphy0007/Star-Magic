#include "UQFFComparedToGW190425.h"
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _UQFFComparedToGW190425_T_H(double M){
    return (UQFFComparedToGW190425::HBAR*UQFFComparedToGW190425::C*UQFFComparedToGW190425::C*UQFFComparedToGW190425::C)/
           (8.0*M_PI*UQFFComparedToGW190425::G*M*UQFFComparedToGW190425::K_B);
}
static inline double _UQFFComparedToGW190425_L_H(double M){
    return (UQFFComparedToGW190425::HBAR*std::pow(UQFFComparedToGW190425::C,6.0))/
           (15360.0*M_PI*UQFFComparedToGW190425::G*UQFFComparedToGW190425::G*M*M);
}
static inline double _UQFFComparedToGW190425_r_s(double M){
    return 2.0*UQFFComparedToGW190425::G*M/(UQFFComparedToGW190425::C*UQFFComparedToGW190425::C);
}
static inline double _UQFFComparedToGW190425_tau(double M){
    return 5120.0*M_PI*UQFFComparedToGW190425::G*UQFFComparedToGW190425::G*M*M*M/
           (UQFFComparedToGW190425::HBAR*std::pow(UQFFComparedToGW190425::C,4.0));
}
static inline double _UQFFComparedToGW190425_Um(double r,double t){
    double tn=t/UQFFComparedToGW190425::T_N_REF;
    return (UQFFComparedToGW190425::MU_J/r)*(1.0-std::exp(-UQFFComparedToGW190425::GAMMA*t*std::cos(M_PI*tn)));
}

UQFFComparedToGW190425::UQFFComparedToGW190425(){}
double UQFFComparedToGW190425::compute_chirp_mass(double m1,double m2) const{
    return std::pow(m1*m2,0.6)/std::pow(m1+m2,0.2);
}
double UQFFComparedToGW190425::compute_h_GR(double Mc,double d,double f) const{
    return (4.0/d)*std::pow(G*Mc/(C*C),5.0/3.0)*std::pow(PI*f,2.0/3.0)/(C*C);
}
double UQFFComparedToGW190425::compute_h_UQFF(double m1,double m2,double d,double f,double t) const{
    double Mc=compute_chirp_mass(m1,m2);
    double h_gr=compute_h_GR(Mc,d,f);
    double M=m1+m2; double T_H=_UQFFComparedToGW190425_T_H(M); double r=_UQFFComparedToGW190425_r_s(M);
    double S_SCm=std::exp(-RHO_SCM*r/(K_B*T_H));
    double Um=_UQFFComparedToGW190425_Um(r,t);
    double S_Um=std::exp(-Um*2.0*PI*f/(C*C));
    double h=(1.0-F_TRZ)*S_SCm*S_Um*h_gr;
    for(auto& m:mods_) h+=m(f,t);
    return h;
}
double UQFFComparedToGW190425::compute_post_merger_phase(double m1,double m2) const{
    // UQFF phase: phi_UQFF = phi_GR + kappa*f_TRZ*t_merger
    double t_merger=5.0*(m1+m2)/(C*C*C/G)*(C*C*C/G)/(m1*m2)*1e-6;
    return KAPPA*F_TRZ*std::abs(t_merger);
}
double UQFFComparedToGW190425::compute_ejecta_limit_UQFF(double m1,double m2) const{
    // M_ej,UQFF < M_GR * (rho_SCm/rho_UA) * (1-f_TRZ) 
    double M_ej_GR=0.05*(m1+m2); // typical 5% ejecta
    return M_ej_GR*(RHO_SCM/RHO_UA)*(1.0-F_TRZ);
}
void UQFFComparedToGW190425::simulate_inspiral(double f0,double fe,double df,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"f_Hz,h_GR,h_UQFF\n";}}
    double Mc=compute_chirp_mass(M1_DEFAULT,M2_DEFAULT);
    for(double f=f0;f<=fe;f+=df)
        *os<<f<<","<<compute_h_GR(Mc,D_DEFAULT,f)<<","
           <<compute_h_UQFF(M1_DEFAULT,M2_DEFAULT,D_DEFAULT,f)<<"\n";
    if(ofs.is_open())ofs.close();
}
void UQFFComparedToGW190425::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void UQFFComparedToGW190425::update_from_file(const std::string& fn){
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_UQFFCOMPAREDTOGW190425
int main(){
    UQFFComparedToGW190425 obj;
    double m1=M1_DEFAULT,m2=M2_DEFAULT,d=D_DEFAULT;
    std::cout<<"Chirp mass="<<obj.compute_chirp_mass(m1,m2)/M_SUN<<" Msun\n";
    std::cout<<"h_UQFF@2500Hz="<<obj.compute_h_UQFF(m1,m2,d,2500.0)<<"\n";
    std::cout<<"Ejecta limit UQFF="<<obj.compute_ejecta_limit_UQFF(m1,m2)/M_SUN<<" Msun\n";
    obj.simulate_inspiral(10,2500,20,"gw190425_sweep.csv");
    return 0;
}
#endif
