#include "UQFFComparedToLIGOData.h"
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _UQFFComparedToLIGOData_T_H(double M){
    return (UQFFComparedToLIGOData::HBAR*UQFFComparedToLIGOData::C*UQFFComparedToLIGOData::C*UQFFComparedToLIGOData::C)/
           (8.0*M_PI*UQFFComparedToLIGOData::G*M*UQFFComparedToLIGOData::K_B);
}
static inline double _UQFFComparedToLIGOData_L_H(double M){
    return (UQFFComparedToLIGOData::HBAR*std::pow(UQFFComparedToLIGOData::C,6.0))/
           (15360.0*M_PI*UQFFComparedToLIGOData::G*UQFFComparedToLIGOData::G*M*M);
}
static inline double _UQFFComparedToLIGOData_r_s(double M){
    return 2.0*UQFFComparedToLIGOData::G*M/(UQFFComparedToLIGOData::C*UQFFComparedToLIGOData::C);
}
static inline double _UQFFComparedToLIGOData_tau(double M){
    return 5120.0*M_PI*UQFFComparedToLIGOData::G*UQFFComparedToLIGOData::G*M*M*M/
           (UQFFComparedToLIGOData::HBAR*std::pow(UQFFComparedToLIGOData::C,4.0));
}
static inline double _UQFFComparedToLIGOData_Um(double r,double t){
    double tn=t/UQFFComparedToLIGOData::T_N_REF;
    return (UQFFComparedToLIGOData::MU_J/r)*(1.0-std::exp(-UQFFComparedToLIGOData::GAMMA*t*std::cos(M_PI*tn)));
}

UQFFComparedToLIGOData::UQFFComparedToLIGOData(){}

double UQFFComparedToLIGOData::compute_chirp_mass(double m1,double m2) const{
    return std::pow(m1*m2,0.6)/std::pow(m1+m2,0.2);
}
double UQFFComparedToLIGOData::compute_h_GR(double Mc,double d,double f) const{
    // h_GR = (4/d)*(G*Mc/c^2)^(5/3)*(pi*f)^(2/3)/c^2
    double GM=G*Mc; double pif=PI*f;
    return (4.0/d)*std::pow(GM/C/C,5.0/3.0)*std::pow(pif,2.0/3.0)/(C*C);
}
double UQFFComparedToLIGOData::compute_S_SCm(double M,double T_H) const{
    double rs=_UQFFComparedToLIGOData_r_s(M);
    double arg=RHO_SCM*rs/(K_B*T_H); if(arg>700)arg=700;
    return std::exp(-arg);
}
double UQFFComparedToLIGOData::compute_S_Um(double M,double r,double t,double f) const{
    if(r<=0) r=_UQFFComparedToLIGOData_r_s(M);
    double Um=_UQFFComparedToLIGOData_Um(r,t);
    double arg=Um*2.0*PI*f/(C*C); if(arg>700)arg=700;
    return std::exp(-arg);
}
double UQFFComparedToLIGOData::compute_h_UQFF(double m1,double m2,double d,double f,
                              double r,double t) const{
    double Mc=compute_chirp_mass(m1,m2);
    double h_gr=compute_h_GR(Mc,d,f);
    double M=m1+m2;
    double T_H=_UQFFComparedToLIGOData_T_H(M);
    double S_SCm=compute_S_SCm(M,T_H);
    double S_Um=compute_S_Um(M,r,t,f);
    double h=(1.0-F_TRZ)*S_SCm*S_Um*h_gr;
    for(auto& m:mods_) h+=m(f,t);
    return h;
}
double UQFFComparedToLIGOData::compute_delta_phi(double t_coal) const{
    return KAPPA*F_TRZ*t_coal;
}
void UQFFComparedToLIGOData::simulate_frequency_sweep(double f0,double f1,double df,
    double m1,double m2,double d,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"f_Hz,h_GR,h_UQFF\n";}}
    double M=m1+m2; double Mc=compute_chirp_mass(m1,m2);
    double T_H=_UQFFComparedToLIGOData_T_H(M);
    for(double f=f0;f<=f1;f+=df){
        double h_gr=compute_h_GR(Mc,d,f);
        double h_u=compute_h_UQFF(m1,m2,d,f);
        *os<<f<<","<<h_gr<<","<<h_u<<"\n";
    }
    if(ofs.is_open())ofs.close();
}
void UQFFComparedToLIGOData::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void UQFFComparedToLIGOData::update_from_file(const std::string& fn){
    std::ifstream fi(fn); if(!fi.is_open())return;
    std::string l; while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_UQFFCOMPAREDTOLIGODATA
int main(){
    UQFFComparedToLIGOData obj;
    double m1=7.16e31,m2=5.77e31,d=1.27e25;
    std::cout<<"Chirp mass="<<obj.compute_chirp_mass(m1,m2)<<"\n";
    std::cout<<"h_UQFF@150Hz="<<obj.compute_h_UQFF(m1,m2,d,150.0)<<"\n";
    obj.simulate_frequency_sweep(10,500,10,m1,m2,d,"ligo_sweep.csv");
    return 0;
}
#endif
