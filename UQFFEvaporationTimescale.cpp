#include "UQFFEvaporationTimescale.h"
UQFFEvaporationTimescale::UQFFEvaporationTimescale(){}
double UQFFEvaporationTimescale::compute_tau_standard(double M) const{return _tau_std(M);}
double UQFFEvaporationTimescale::compute_tau_UQFF(double M) const{
    double tau=_tau_std(M)/(1.0-F_TRZ)*(RHO_UA/RHO_SCM);
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=Um/(K_B*T_H); if(ex>700)ex=700;
    tau*=std::exp(ex);
    for(auto& m:mods_) tau*=m(M);
    return tau;
}
double UQFFEvaporationTimescale::compute_factor(double M) const{
    double ts=compute_tau_standard(M);
    if(ts<1e-300) return 1.0;
    return compute_tau_UQFF(M)/ts;
}
double UQFFEvaporationTimescale::compute_M_cross_standard() const{
    return std::cbrt(HBAR*std::pow(C,4)*T_HUBBLE/(5120.0*PI*G*G));
}
double UQFFEvaporationTimescale::compute_M_cross_UQFF() const{
    // tau_UQFF ~ 30 * tau_std; tau_std(M) = 5120 pi G² M³ / (hbar c^4)
    // tau_UQFF = t_H → M^3 = t_H * hbar c^4 / (5120 pi G² * 30)
    double factor=1.0/(1.0-F_TRZ)*(RHO_UA/RHO_SCM)*std::exp(1.0); // ~30 approx
    return std::cbrt(HBAR*std::pow(C,4)*T_HUBBLE/(5120.0*PI*G*G*factor));
}
void UQFFEvaporationTimescale::sensitivity_Um_sweep(double Um0,double Um1,double dUm,double M,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"Um_ratio,tau_UQFF_s,factor\n";}}
    double T_H=_T_H(M);
    for(double Um_r=Um0;Um_r<=Um1;Um_r+=dUm){
        double Um_val=Um_r*K_B*T_H;
        double tau=_tau_std(M)/(1.0-F_TRZ)*(RHO_UA/RHO_SCM)*std::exp(std::min(Um_r,700.0));
        *os<<Um_r<<","<<tau<<","<<tau/_tau_std(M)<<"\n";
    }
    if(f.is_open())f.close();
}
void UQFFEvaporationTimescale::simulate_timescale_sweep(double M0,double M1,int n,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"M_kg,tau_std_s,tau_UQFF_s,factor\n";}}
    double dM=(M1-M0)/std::max(n-1,1);
    for(int i=0;i<n;i++){
        double M=M0+i*dM;
        *os<<M<<","<<compute_tau_standard(M)<<","<<compute_tau_UQFF(M)<<","<<compute_factor(M)<<"\n";
    }
    if(f.is_open())f.close();
}
void UQFFEvaporationTimescale::add_mod(std::function<double(double)> m){mods_.push_back(std::move(m));}
void UQFFEvaporationTimescale::update_from_file(const std::string& fn){std::ifstream fi(fn);}
#ifdef STANDALONE_UQFFEVAPORATIONTIMESCALE
int main(){
    UQFFEvaporationTimescale evap;
    std::cout<<"M_cross_std="<<evap.compute_M_cross_standard()<<" kg\n";
    std::cout<<"M_cross_UQFF="<<evap.compute_M_cross_UQFF()<<" kg\n";
    double M=1e12;
    std::cout<<"tau_std(1e12kg)="<<evap.compute_tau_standard(M)<<" s\n";
    std::cout<<"tau_UQFF(1e12kg)="<<evap.compute_tau_UQFF(M)<<" s\n";
    evap.sensitivity_Um_sweep(0.0,5.0,0.1,M,"evap_Um_sensitivity.csv");
    evap.simulate_timescale_sweep(1e8,1e20,100,"evap_timescale_sweep.csv");
    return 0;
}
#endif
