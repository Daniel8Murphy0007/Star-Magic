#include "UQFFPBHDarkMatter.h"
UQFFPBHDarkMatter::UQFFPBHDarkMatter(){}
double UQFFPBHDarkMatter::compute_tau_standard(double M) const { return _tau_std(M); }
double UQFFPBHDarkMatter::compute_T_H(double M) const { return _T_H(M); }
double UQFFPBHDarkMatter::compute_tau_prime(double tau_std) const { return tau_std/(1.0-F_TRZ); }
double UQFFPBHDarkMatter::compute_tau_double_prime(double tau_prime) const { return tau_prime*(RHO_UA/RHO_SCM); }
double UQFFPBHDarkMatter::compute_tau_UQFF(double M) const {
    double tau=compute_tau_double_prime(compute_tau_prime(compute_tau_standard(M)));
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=Um/(K_B*T_H); if(ex>700)ex=700;
    tau*=std::exp(ex);
    for(auto& m:mods_) tau*=m(M);
    return tau;
}
double UQFFPBHDarkMatter::compute_f_PBH_enhancement(double M) const {
    // Ratio tau_UQFF/tau_std ~ elevation factor applied to f_PBH
    double tau_s=compute_tau_standard(M);
    if(tau_s<1e-300) return 1.0;
    return compute_tau_UQFF(M)/tau_s;
}
bool UQFFPBHDarkMatter::is_DM_candidate(double M) const {
    // Standard: only M > ~1e12 kg survive; UQFF: window M > ~1e11 kg
    return compute_tau_UQFF(M) > 4.35e17; // older than Hubble time
}
void UQFFPBHDarkMatter::simulate_lifetime_sweep(double M0,double M1,double dM,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"M_kg,tau_std_s,tau_UQFF_s,DM_candidate\n";}}
    for(double M=M0;M<=M1;M+=dM)
        *os<<M<<","<<compute_tau_standard(M)<<","<<compute_tau_UQFF(M)<<","<<(is_DM_candidate(M)?1:0)<<"\n";
    if(f.is_open())f.close();
}
void UQFFPBHDarkMatter::add_modifier(std::function<double(double)> mod){mods_.push_back(std::move(mod));}
void UQFFPBHDarkMatter::update_from_file(const std::string& fn){std::ifstream fi(fn);}
#ifdef STANDALONE_UQFFPBHDARKMATTER
int main(){
    UQFFPBHDarkMatter pbh;
    double M=1e12; // Small PBH 1e12 kg
    std::cout<<"tau_std="<<pbh.compute_tau_standard(M)<<" s\n";
    std::cout<<"tau_UQFF="<<pbh.compute_tau_UQFF(M)<<" s\n";
    std::cout<<"DM candidate: "<<(pbh.is_DM_candidate(M)?"yes":"no")<<"\n";
    pbh.simulate_lifetime_sweep(1e9,1e15,1e12,"pbh_dm_sweep.csv");
    return 0;
}
#endif
