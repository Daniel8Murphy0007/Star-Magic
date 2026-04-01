#include "UQFFComparedToGW170817.h"
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _UQFFComparedToGW170817_T_H(double M){
    return (UQFFComparedToGW170817::HBAR*UQFFComparedToGW170817::C*UQFFComparedToGW170817::C*UQFFComparedToGW170817::C)/
           (8.0*M_PI*UQFFComparedToGW170817::G*M*UQFFComparedToGW170817::K_B);
}
static inline double _UQFFComparedToGW170817_L_H(double M){
    return (UQFFComparedToGW170817::HBAR*std::pow(UQFFComparedToGW170817::C,6.0))/
           (15360.0*M_PI*UQFFComparedToGW170817::G*UQFFComparedToGW170817::G*M*M);
}
static inline double _UQFFComparedToGW170817_r_s(double M){
    return 2.0*UQFFComparedToGW170817::G*M/(UQFFComparedToGW170817::C*UQFFComparedToGW170817::C);
}
static inline double _UQFFComparedToGW170817_tau(double M){
    return 5120.0*M_PI*UQFFComparedToGW170817::G*UQFFComparedToGW170817::G*M*M*M/
           (UQFFComparedToGW170817::HBAR*std::pow(UQFFComparedToGW170817::C,4.0));
}
static inline double _UQFFComparedToGW170817_Um(double r,double t){
    double tn=t/UQFFComparedToGW170817::T_N_REF;
    return (UQFFComparedToGW170817::MU_J/r)*(1.0-std::exp(-UQFFComparedToGW170817::GAMMA*t*std::cos(M_PI*tn)));
}

UQFFComparedToGW170817::UQFFComparedToGW170817(){}
double UQFFComparedToGW170817::compute_chirp_mass(double m1,double m2) const{
    return std::pow(m1*m2,0.6)/std::pow(m1+m2,0.2);
}
double UQFFComparedToGW170817::compute_h_GR(double Mc,double d,double f) const{
    return (4.0/d)*std::pow(G*Mc/(C*C),5.0/3.0)*std::pow(PI*f,2.0/3.0)/(C*C);
}
double UQFFComparedToGW170817::compute_h_UQFF(double m1,double m2,double d,double f,double t) const{
    double Mc=compute_chirp_mass(m1,m2);
    double h_gr=compute_h_GR(Mc,d,f);
    double M=m1+m2; double T_H=_UQFFComparedToGW170817_T_H(M); double r=_UQFFComparedToGW170817_r_s(M);
    double S_SCm=std::exp(-RHO_SCM*r/(K_B*T_H));
    double Um=_UQFFComparedToGW170817_Um(r,t);
    double S_Um=std::exp(-Um*2.0*PI*f/(C*C));
    double h=(1.0-F_TRZ)*S_SCm*S_Um*h_gr;
    for(auto& m:mods_) h+=m(f,t);
    return h;
}
double UQFFComparedToGW170817::compute_grb_delay_UQFF() const{
    // Delta_t_UQFF = 1.7 * (1 + f_TRZ * rho_UA/rho_SCm)
    return DT_GRB*(1.0+F_TRZ*(RHO_UA/RHO_SCM));
}
double UQFFComparedToGW170817::compute_tidal_UQFF(double kappa_GR) const{
    return kappa_GR*(1.0-F_TRZ*RHO_SCM/RHO_UA);
}
void UQFFComparedToGW170817::simulate_inspiral(double f0,double fp,double df,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"f_Hz,h_GR,h_UQFF\n";}}
    double m1=M1_GW170817,m2=M2_GW170817,d=D_GW170817;
    double Mc=compute_chirp_mass(m1,m2);
    for(double f=f0;f<=fp;f+=df){
        *os<<f<<","<<compute_h_GR(Mc,d,f)<<","<<compute_h_UQFF(m1,m2,d,f)<<"\n";
    }
    if(ofs.is_open())ofs.close();
}
void UQFFComparedToGW170817::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void UQFFComparedToGW170817::update_from_file(const std::string& fn){
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_UQFFCOMPAREDTOGW170817
int main(){
    UQFFComparedToGW170817 obj;
    std::cout<<"GRB delay UQFF="<<obj.compute_grb_delay_UQFF()<<" s\n";
    std::cout<<"h_UQFF@1500Hz="<<obj.compute_h_UQFF(M1_GW170817,M2_GW170817,D_GW170817,1500.0)<<"\n";
    obj.simulate_inspiral(10,1500,10,"gw170817_sweep.csv");
    return 0;
}
#endif
