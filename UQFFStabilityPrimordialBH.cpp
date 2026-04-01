#include "UQFFStabilityPrimordialBH.h"
UQFFStabilityPrimordialBH::UQFFStabilityPrimordialBH(){}
double UQFFStabilityPrimordialBH::compute_tau_std(double M) const{return _tau_std(M);}
double UQFFStabilityPrimordialBH::compute_tau_UQFF(double M) const{
    double tau=_tau_std(M)/(1.0-F_TRZ);
    tau*=(RHO_UA/RHO_SCM);
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=Um/(K_B*T_H); if(ex>700)ex=700;
    tau*=std::exp(ex);
    for(auto& m:mods_) tau*=m(M);
    return tau;
}
std::string UQFFStabilityPrimordialBH::classify(double M) const{
    double tu=compute_tau_UQFF(M);
    if(tu>=T_HUBBLE) return "stable_DM";
    if(tu>=0.1*T_HUBBLE) return "marginal";
    return "evaporating";
}
double UQFFStabilityPrimordialBH::pbh_dm_window_min_mass_uqff() const{
    // Binary search for minimum M such that tau_UQFF >= T_HUBBLE
    double lo=1e6, hi=1e25;
    for(int i=0;i<80;i++){
        double mid=(lo+hi)/2.0;
        if(compute_tau_UQFF(mid)>=T_HUBBLE) hi=mid; else lo=mid;
    }
    return (lo+hi)/2.0;
}
void UQFFStabilityPrimordialBH::simulate_mass_stability_map(double M0,double M1,int n,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"M_kg,tau_std_s,tau_UQFF_s,class\n";}}
    double dM=(M1-M0)/std::max(n-1,1);
    for(int i=0;i<n;i++){
        double M=M0+i*dM;
        *os<<M<<","<<compute_tau_std(M)<<","<<compute_tau_UQFF(M)<<","<<classify(M)<<"\n";
    }
    if(f.is_open())f.close();
}
void UQFFStabilityPrimordialBH::add_modifier(std::function<double(double)> m){mods_.push_back(std::move(m));}
void UQFFStabilityPrimordialBH::update_from_file(const std::string& fn){std::ifstream fi(fn);}
#ifdef STANDALONE_UQFFSTABILITYPRIMORDIALBH
int main(){
    UQFFStabilityPrimordialBH pbhs;
    double M=1e12;
    std::cout<<"tau_std="<<pbhs.compute_tau_std(M)<<" s\n";
    std::cout<<"tau_UQFF="<<pbhs.compute_tau_UQFF(M)<<" s\n";
    std::cout<<"class="<<pbhs.classify(M)<<"\n";
    std::cout<<"Min DM mass UQFF="<<pbhs.pbh_dm_window_min_mass_uqff()<<" kg\n";
    pbhs.simulate_mass_stability_map(1e8,1e18,100,"pbh_stability_map.csv");
    return 0;
}
#endif
