#include "UQFFModulationForM87.h"
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _UQFFModulationForM87_T_H(double M){
    return (UQFFModulationForM87::HBAR*UQFFModulationForM87::C*UQFFModulationForM87::C*UQFFModulationForM87::C)/
           (8.0*M_PI*UQFFModulationForM87::G*M*UQFFModulationForM87::K_B);
}
static inline double _UQFFModulationForM87_L_H(double M){
    return (UQFFModulationForM87::HBAR*std::pow(UQFFModulationForM87::C,6.0))/
           (15360.0*M_PI*UQFFModulationForM87::G*UQFFModulationForM87::G*M*M);
}
static inline double _UQFFModulationForM87_r_s(double M){
    return 2.0*UQFFModulationForM87::G*M/(UQFFModulationForM87::C*UQFFModulationForM87::C);
}
static inline double _UQFFModulationForM87_tau(double M){
    return 5120.0*M_PI*UQFFModulationForM87::G*UQFFModulationForM87::G*M*M*M/
           (UQFFModulationForM87::HBAR*std::pow(UQFFModulationForM87::C,4.0));
}
static inline double _UQFFModulationForM87_Um(double r,double t){
    double tn=t/UQFFModulationForM87::T_N_REF;
    return (UQFFModulationForM87::MU_J/r)*(1.0-std::exp(-UQFFModulationForM87::GAMMA*t*std::cos(M_PI*tn)));
}

UQFFModulationForM87::UQFFModulationForM87(){}
double UQFFModulationForM87::compute_T_H(double M) const{return _UQFFModulationForM87_T_H(M);}
double UQFFModulationForM87::compute_T_UQFF(double M) const{
    double T_H=_UQFFModulationForM87_T_H(M);
    return T_H*(1.0+F_TRZ)*(1.0-RHO_SCM/RHO_UA);
}
double UQFFModulationForM87::compute_shadow_radius(double M) const{
    return 3.0*std::sqrt(3.0)*G*M/(C*C);
}
double UQFFModulationForM87::compute_shadow_UQFF(double M) const{
    return compute_shadow_radius(M)*std::sqrt(1.0+F_TRZ*(RHO_UA/RHO_SCM));
}
double UQFFModulationForM87::compute_jet_power_UQFF(double P_BZ) const{
    return P_BZ*(1.0+F_TRZ)*std::sqrt(RHO_UA/RHO_SCM);
}
double UQFFModulationForM87::compute_ring_brightness_UQFF(double B_GR) const{
    return B_GR*std::pow(RHO_UA/RHO_SCM,F_TRZ/4.0);
}
double UQFFModulationForM87::compute_T_accretion_UQFF(double T_acc,double M) const{
    double T_H=_UQFFModulationForM87_T_H(M); double r=_UQFFModulationForM87_r_s(M);
    double Um=_UQFFModulationForM87_Um(r,T_N_REF);
    return T_acc*(1.0+Um/(K_B*T_H));
}
void UQFFModulationForM87::simulate_T_evolution(double M0,double M1,double dM,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"M_Msun,T_H,T_UQFF,shadow_UQFF_m\n";}}
    for(double M=M0;M<=M1;M+=dM){
        *os<<M/M_SUN<<","<<compute_T_H(M)<<","<<compute_T_UQFF(M)<<","<<compute_shadow_UQFF(M)<<"\n";
    }
    if(ofs.is_open())ofs.close();
}
void UQFFModulationForM87::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void UQFFModulationForM87::update_from_file(const std::string& fn){
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_UQFFMODULATIONFORM87
int main(){
    UQFFModulationForM87 obj;
    std::cout<<"T_H(M87*)="<<obj.compute_T_H()<<" K\n";
    std::cout<<"T_UQFF(M87*)="<<obj.compute_T_UQFF()<<" K\n";
    std::cout<<"Shadow GR="<<obj.compute_shadow_radius()/1e9<<" Gm\n";
    std::cout<<"Shadow UQFF="<<obj.compute_shadow_UQFF()/1e9<<" Gm\n";
    std::cout<<"Jet UQFF="<<obj.compute_jet_power_UQFF()<<" W\n";
    obj.simulate_T_evolution(1e9*1.989e30,15e9*1.989e30,0.5e9*1.989e30,"m87_modulation.csv");
    return 0;
}
#endif
