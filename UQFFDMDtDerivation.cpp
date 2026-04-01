#include "UQFFDMDtDerivation.h"
UQFFDMDtDerivation::UQFFDMDtDerivation(){}
double UQFFDMDtDerivation::compute_dM_dt_standard(double M) const{return -_L_H(M)/(C*C);}
double UQFFDMDtDerivation::compute_dM_dt_step1(double M) const{return compute_dM_dt_standard(M)*(1.0-F_TRZ);}
double UQFFDMDtDerivation::compute_dM_dt_step2(double M) const{return compute_dM_dt_step1(M)*(RHO_SCM/RHO_UA);}
double UQFFDMDtDerivation::compute_dM_dt_UQFF(double M) const{
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=-Um/(K_B*T_H); if(ex<-700)ex=-700;
    double dMdt=compute_dM_dt_step2(M)*std::exp(ex);
    for(auto& m:mods_) dMdt*=m(M);
    return dMdt;
}
double UQFFDMDtDerivation::compute_suppression_factor(double M) const{
    double std_rate=compute_dM_dt_standard(M);
    if(std::abs(std_rate)<1e-300) return 1.0;
    return compute_dM_dt_UQFF(M)/std_rate;
}
double UQFFDMDtDerivation::compute_M_at_t(double M0,double t) const{
    // Analytic: dM/dt ~ -A/M² → M(t) = (M0³ - 3 A t)^(1/3)
    // A = |dM/dt| * M² evaluated at M0 (constant suppression approx)
    double rate=std::abs(compute_dM_dt_UQFF(M0));
    double A=rate*M0*M0;
    double val=M0*M0*M0-3.0*A*t;
    if(val<0) return 0.0;
    return std::cbrt(val);
}
void UQFFDMDtDerivation::simulate_evaporation(double M0,double t0,double t1,double dt,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"t_s,M_kg,dM_dt_UQFF,suppression\n";}}
    double M=M0;
    for(double t=t0;t<=t1&&M>0;t+=dt){
        double dMdt=compute_dM_dt_UQFF(M);
        *os<<t<<","<<M<<","<<dMdt<<","<<compute_suppression_factor(M)<<"\n";
        M+=dMdt*dt; if(M<0)M=0;
    }
    if(f.is_open())f.close();
}
void UQFFDMDtDerivation::add_mod(std::function<double(double)> m){mods_.push_back(std::move(m));}
void UQFFDMDtDerivation::update_from_file(const std::string& fn){std::ifstream fi(fn);}
#ifdef STANDALONE_UQFFDMDTDERIVATION
int main(){
    UQFFDMDtDerivation dm;
    double M=2e30;
    std::cout<<"dM/dt_std="<<dm.compute_dM_dt_standard(M)<<" kg/s\n";
    std::cout<<"dM/dt_UQFF="<<dm.compute_dM_dt_UQFF(M)<<" kg/s\n";
    std::cout<<"suppression="<<dm.compute_suppression_factor(M)<<"\n";
    std::cout<<"M at 1e67 s (analytic)="<<dm.compute_M_at_t(2e30,1e67)<<" kg\n";
    dm.simulate_evaporation(2e30,0,1e67,1e57,"dmdt_evaporation.csv");
    return 0;
}
#endif
