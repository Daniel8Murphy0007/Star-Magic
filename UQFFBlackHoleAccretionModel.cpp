#include "UQFFBlackHoleAccretionModel.h"
UQFFBlackHoleAccretionModel::UQFFBlackHoleAccretionModel(){}
double UQFFBlackHoleAccretionModel::compute_cs(double T_inf) const{
    return std::sqrt(5.0/3.0*K_B*T_inf/M_P);
}
double UQFFBlackHoleAccretionModel::compute_Mdot_Bondi(double M,double rho_inf,double T_inf) const{
    double cs=compute_cs(T_inf);
    return 4.0*PI*0.25*G*G*M*M*rho_inf/(cs*cs*cs);
}
double UQFFBlackHoleAccretionModel::compute_Mdot_UQFF(double M,double rho_inf,double T_inf) const{
    double rho_eff=rho_inf+RHO_UA-RHO_SCM;
    double s_trz=1.0+F_TRZ;
    double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=-Um/(K_B*T_inf); if(ex<-700)ex=-700;
    double s_um=1.0-std::exp(ex);
    double Mdot=compute_Mdot_Bondi(M,rho_inf,T_inf)*(rho_eff/rho_inf)*s_trz*s_um;
    for(auto& m:mods_) Mdot*=m(M,T_inf);
    return Mdot;
}
double UQFFBlackHoleAccretionModel::compute_L_Edd(double M) const{
    return 4.0*PI*G*M*M_P*C/SIGMA_T;
}
double UQFFBlackHoleAccretionModel::compute_Mdot_Edd(double M) const{
    return compute_L_Edd(M)/(ETA*C*C);
}
double UQFFBlackHoleAccretionModel::compute_Eddington_ratio(double M,double rho_inf,double T_inf) const{
    double Md_Edd=compute_Mdot_Edd(M);
    if(Md_Edd<1e-300) return 0.0;
    return compute_Mdot_UQFF(M,rho_inf,T_inf)/Md_Edd;
}
void UQFFBlackHoleAccretionModel::simulate_M_evolution(double M0,double t0,double t1,double dt,double rho_inf,double T_inf,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"t_s,M_kg,Mdot_UQFF,Mdot_Edd\n";}}
    double M=M0;
    for(double t=t0;t<=t1;t+=dt){
        double Md=compute_Mdot_UQFF(M,rho_inf,T_inf);
        double Ml=-_L_H(M)/(C*C);
        *os<<t<<","<<M<<","<<Md<<","<<compute_Mdot_Edd(M)<<"\n";
        M+=(Md+Ml)*dt;
        if(M<0) M=0;
    }
    if(f.is_open())f.close();
}
void UQFFBlackHoleAccretionModel::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void UQFFBlackHoleAccretionModel::update_from_file(const std::string& fn){std::ifstream fi(fn);}
#ifdef STANDALONE_UQFFBLACKHOLEACCRETIONMODEL
int main(){
    UQFFBlackHoleAccretionModel acc;
    double M=8.55e36, rho=1e-19, T=1e7; // typical SMBH environment
    std::cout<<"Mdot_Bondi="<<acc.compute_Mdot_Bondi(M,rho,T)<<" kg/s\n";
    std::cout<<"Mdot_UQFF="<<acc.compute_Mdot_UQFF(M,rho,T)<<" kg/s\n";
    std::cout<<"Edd_ratio="<<acc.compute_Eddington_ratio(M,rho,T)<<"\n";
    acc.simulate_M_evolution(8.55e36,0,3.15e13,3.15e10,rho,T,"accretion_evolution.csv");
    return 0;
}
#endif
