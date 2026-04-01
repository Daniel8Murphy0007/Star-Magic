#include "UQFFBlackHoleInversion.h"
UQFFBlackHoleInversion::UQFFBlackHoleInversion(unsigned int seed):rng_(seed){}
double UQFFBlackHoleInversion::compute_r_s_UQFF(double M) const { return _r_s(M)*(1.0-RHO_SCM/RHO_UA); }
double UQFFBlackHoleInversion::compute_E_inv(double M) const {
    double rs=compute_r_s_UQFF(M); return G*M*M/rs;
}
double UQFFBlackHoleInversion::compute_P_inv(double M) const {
    double T_H=_T_H(M); double E=compute_E_inv(M);
    double ex=-E/(K_B*T_H); if(ex<-700)ex=-700;
    return F_TRZ*std::exp(ex);
}
double UQFFBlackHoleInversion::compute_Phi_inv(double M) const {
    double delta_rho=RHO_SCM/RHO_UA;
    return (1.0/delta_rho)*(G*M/C)*(1.0+F_TRZ);
}
double UQFFBlackHoleInversion::compute_Theta_inv(double M) const {
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=Um/(K_B*T_H); if(ex>700)ex=700;
    double S_Um=std::exp(ex);
    double Theta=compute_P_inv(M)*compute_Phi_inv(M)*S_Um;
    for(auto& m:mods_) Theta+=m(M);
    return Theta;
}
double UQFFBlackHoleInversion::compute_P_invert_MC(double M,int n,double sigma) const{
    std::normal_distribution<double> nd(0.0,1.0);
    int cnt=0;
    for(int i=0;i<n;i++){
        // Perturb vacuum densities
        double pUA =RHO_UA *std::exp(sigma*nd(rng_));
        double pSCM=RHO_SCM*std::exp(sigma*nd(rng_));
        double T_H=HBAR*C*C*C/(8.0*PI*G*M*K_B);
        double rs_u=_r_s(M)*(1.0-pSCM/pUA);
        double E=G*M*M/rs_u;
        double P=0.1*std::exp(std::max(-E/(K_B*T_H),-700.0));
        double Ph=(pUA/pSCM)*(G*M/C)*(1.0+0.1);
        double r=_r_s(M); double t=T_N_REF;
        double tn=t/T_N_REF;
        double Um=(MU_J/r)*(1.0-std::exp(-GAMMA*t*std::cos(PI*tn)));
        double ex2=Um/(K_B*T_H); if(ex2>700)ex2=700;
        if(P*Ph*std::exp(ex2)>1.0) cnt++;
    }
    return (double)cnt/n;
}
void UQFFBlackHoleInversion::simulate_P_invert_sweep(double M0,double M1,double dM,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"M_kg,Theta_inv,P_invert\n";}}
    for(double M=M0;M<=M1;M+=dM)
        *os<<M<<","<<compute_Theta_inv(M)<<","<<compute_P_invert_MC(M,200,0.05)<<"\n";
    if(f.is_open())f.close();
}
void UQFFBlackHoleInversion::add_mod(std::function<double(double)> mod){mods_.push_back(std::move(mod));}
void UQFFBlackHoleInversion::update_from_file(const std::string& fn){std::ifstream fi(fn);}
#ifdef STANDALONE_UQFFBLACKHOLEINVERSION
int main(){
    UQFFBlackHoleInversion inv;
    double M=8.55e36;
    std::cout<<"Theta_inv="<<inv.compute_Theta_inv(M)<<"\n";
    std::cout<<"P_invert_MC="<<inv.compute_P_invert_MC(M,5000,0.05)<<"\n";
    inv.simulate_P_invert_sweep(1e30,1e40,1e38,"bh_inversion_sweep.csv");
    return 0;
}
#endif
