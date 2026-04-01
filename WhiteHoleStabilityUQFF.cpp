#include "WhiteHoleStabilityUQFF.h"
WhiteHoleStabilityUQFF::WhiteHoleStabilityUQFF(){}
double WhiteHoleStabilityUQFF::compute_L_WH(double M) const { return _L_H(M); }
double WhiteHoleStabilityUQFF::compute_T_WH(double M) const { return -_T_H(M); } // negative conceptual
double WhiteHoleStabilityUQFF::compute_L_prime(double L) const { return L*(1.0-F_TRZ); }
double WhiteHoleStabilityUQFF::compute_L_double_prime(double L2) const {
    double ratio=std::abs(1.0-RHO_UA/RHO_SCM); // ~9
    if(ratio<1e-10) return L2;
    return L2/ratio;  // L' * 0.1 confinement
}
double WhiteHoleStabilityUQFF::compute_L_UQFF(double Ldp, double T_WH) const {
    double absT=std::abs(T_WH);
    double r=1e10; // representative r
    double Um=_U_m(r,T_N_REF);
    double ex=-Um/(K_B*absT); if(ex<-700)ex=-700;
    return Ldp*std::exp(ex);
}
double WhiteHoleStabilityUQFF::compute_tau_standard(double M) const {
    return _r_s(M)/C;  // instability timescale
}
double WhiteHoleStabilityUQFF::compute_tau_UQFF(double M) const {
    double tau=compute_tau_standard(M);
    double T_WH=std::abs(compute_T_WH(M));
    // Proof 1
    tau/=(1.0-F_TRZ);
    // Proof 2
    double grad=std::abs(1.0-RHO_UA/RHO_SCM);
    tau*=grad;
    // Proof 3
    double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=Um/(K_B*T_WH); if(ex>700)ex=700;
    tau*=std::exp(ex);
    for(auto& m:mods_) tau*=m(M);
    return tau;
}
double WhiteHoleStabilityUQFF::compute_stability_factor(double M) const {
    double tau_s=compute_tau_standard(M);
    if(tau_s<1e-300) return 1.0;
    return compute_tau_UQFF(M)/tau_s;
}
void WhiteHoleStabilityUQFF::simulate_over_mass(double M0,double M1,double dM,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"M_kg,tau_std_s,tau_UQFF_s,factor\n";}}
    for(double M=M0;M<=M1;M+=dM)
        *os<<M<<","<<compute_tau_standard(M)<<","<<compute_tau_UQFF(M)<<","<<compute_stability_factor(M)<<"\n";
    if(f.is_open())f.close();
}
void WhiteHoleStabilityUQFF::add_mod(std::function<double(double)> m){mods_.push_back(std::move(m));}
void WhiteHoleStabilityUQFF::update_from_file(const std::string& fn){std::ifstream fi(fn);}
#ifdef STANDALONE_WHITEHOLESTABILITYUQFF
int main(){
    WhiteHoleStabilityUQFF whs;
    double M=8.55e36;
    std::cout<<"tau_std="<<whs.compute_tau_standard(M)<<" s\n";
    std::cout<<"tau_UQFF="<<whs.compute_tau_UQFF(M)<<" s\n";
    std::cout<<"stability_factor="<<whs.compute_stability_factor(M)<<"\n";
    whs.simulate_over_mass(1e30,1e40,1e38,"wh_stability_sweep.csv");
    return 0;
}
#endif
