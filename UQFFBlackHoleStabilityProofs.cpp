#include "UQFFBlackHoleStabilityProofs.h"
UQFFBlackHoleStabilityProofs::UQFFBlackHoleStabilityProofs(){}
double UQFFBlackHoleStabilityProofs::compute_tau_Hawking(double M) const{return _tau_std(M);}
double UQFFBlackHoleStabilityProofs::compute_E_barrier(double M) const{return K_B*_T_H(M)*(RHO_SCM/RHO_UA);}
double UQFFBlackHoleStabilityProofs::compute_tau_proof1(double M) const{return _tau_std(M)/(1.0-F_TRZ);}
double UQFFBlackHoleStabilityProofs::compute_tau_proof2(double M) const{return compute_tau_proof1(M)*(RHO_UA/RHO_SCM);}
double UQFFBlackHoleStabilityProofs::compute_tau_proof3(double M) const{
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=Um/(K_B*T_H); if(ex>700)ex=700;
    return compute_tau_proof2(M)*std::exp(ex);
}
double UQFFBlackHoleStabilityProofs::compute_tau_UQFF(double M) const{
    double tau=compute_tau_proof3(M);
    for(auto& m:mods_) tau*=m(M);
    return tau;
}
double UQFFBlackHoleStabilityProofs::compute_stability_factor(double M) const{
    double ts=_tau_std(M); if(ts<1e-300) return 1.0;
    return compute_tau_UQFF(M)/ts;
}
void UQFFBlackHoleStabilityProofs::prove_stability(double M) const{
    std::cout<<"=== UQFF BH Stability Proofs for M="<<M<<" kg ===\n"
             <<"  tau_Hawking = "<<compute_tau_Hawking(M)<<" s\n"
             <<"  Proof 1 (f_TRZ): tau' = "<<compute_tau_proof1(M)<<" s (x"<<1.0/(1.0-F_TRZ)<<")\n"
             <<"  Proof 2 (rho):   tau'' = "<<compute_tau_proof2(M)<<" s (x"<<RHO_UA/RHO_SCM<<")\n"
             <<"  Proof 3 (U_m):   tau_UQFF = "<<compute_tau_proof3(M)<<" s\n"
             <<"  Proof 4 factor = "<<compute_stability_factor(M)<<"\n";
}
void UQFFBlackHoleStabilityProofs::simulate_stability_sweep(double M0,double M1,double dM,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"M_kg,tau_Hawking_s,tau_UQFF_s,factor\n";}}
    for(double M=M0;M<=M1;M+=dM)
        *os<<M<<","<<compute_tau_Hawking(M)<<","<<compute_tau_UQFF(M)<<","<<compute_stability_factor(M)<<"\n";
    if(f.is_open())f.close();
}
void UQFFBlackHoleStabilityProofs::add_mod(std::function<double(double)> m){mods_.push_back(std::move(m));}
void UQFFBlackHoleStabilityProofs::update_from_file(const std::string& fn){std::ifstream fi(fn);}
#ifdef STANDALONE_UQFFBLACKHOLESTABILITYPROOFS
int main(){
    UQFFBlackHoleStabilityProofs sp;
    double M=8.55e36;
    sp.prove_stability(M);
    sp.simulate_stability_sweep(1e20,1e40,1e38,"bh_stability_proofs.csv");
    return 0;
}
#endif
