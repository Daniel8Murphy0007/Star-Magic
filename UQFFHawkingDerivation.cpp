#include "UQFFHawkingDerivation.h"
UQFFHawkingDerivation::UQFFHawkingDerivation(){}
double UQFFHawkingDerivation::compute_T_standard(double M) const { return _T_H(M); }
double UQFFHawkingDerivation::compute_T_UQFF(double M) const {
    return _T_H(M)*(1.0+F_TRZ)*(1.0-RHO_SCM/RHO_UA);
}
double UQFFHawkingDerivation::compute_L_standard(double M) const { return _L_H(M); }
double UQFFHawkingDerivation::compute_L_UQFF(double M) const {
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=-Um/(K_B*T_H); if(ex<-700)ex=-700;
    double L=_L_H(M)*std::exp(ex);
    for(auto& t:terms_) L+=t(M,T_H);
    return L;
}
double UQFFHawkingDerivation::compute_dM_dt_standard(double M) const { return -_L_H(M)/(C*C); }
double UQFFHawkingDerivation::compute_dM_dt_UQFF(double M) const { return -compute_L_UQFF(M)/(C*C); }
void UQFFHawkingDerivation::simulate_evaporation(double M0,double t0,double t1,double dt,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"t_s,M_kg,T_UQFF_K,L_UQFF_W\n";}}
    double M=M0;
    for(double t=t0;t<=t1&&M>0;t+=dt){
        double L=compute_L_UQFF(M);
        *os<<t<<","<<M<<","<<compute_T_UQFF(M)<<","<<L<<"\n";
        M+=compute_dM_dt_UQFF(M)*dt;
        if(M<0)M=0;
    }
    if(f.is_open())f.close();
}
void UQFFHawkingDerivation::add_term(std::function<double(double,double)> t){terms_.push_back(std::move(t));}
void UQFFHawkingDerivation::update_from_file(const std::string& fn){std::ifstream fi(fn);}
#ifdef STANDALONE_UQFFHAWKINGDERIVATION
int main(){
    UQFFHawkingDerivation hd;
    double M=2e30;
    std::cout<<"T_std="<<hd.compute_T_standard(M)<<" K\n";
    std::cout<<"T_UQFF="<<hd.compute_T_UQFF(M)<<" K\n";
    std::cout<<"L_std="<<hd.compute_L_standard(M)<<" W\n";
    std::cout<<"L_UQFF="<<hd.compute_L_UQFF(M)<<" W\n";
    std::cout<<"dM/dt_UQFF="<<hd.compute_dM_dt_UQFF(M)<<" kg/s\n";
    hd.simulate_evaporation(2e30,0.0,1e67,1e57,"hawking_evap_uqff.csv");
    return 0;
}
#endif
