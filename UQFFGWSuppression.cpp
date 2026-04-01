#include "UQFFGWSuppression.h"
UQFFGWSuppression::UQFFGWSuppression(){}
double UQFFGWSuppression::compute_P_GW_standard(double m1,double m2,double r) const{
    return (32.0/5.0)*std::pow(G,4)/std::pow(C,5)*(m1*m1)*(m2*m2)*(m1+m2)/std::pow(r,5);
}
double UQFFGWSuppression::compute_S_UA() const{
    return 1.0-RHO_UA/RHO_CRIT;
}
double UQFFGWSuppression::compute_S_SCm(double M) const{
    double T_H=_T_H(M); double rs=_r_s(M);
    double ex=-RHO_SCM*rs/(K_B*T_H); if(ex<-700)ex=-700;
    return std::exp(ex);
}
double UQFFGWSuppression::compute_S_TRZ() const{ return 1.0-F_TRZ; }
double UQFFGWSuppression::compute_S_Um(double M) const{
    double rs=_r_s(M); double omega_GW=C/rs;
    double Um=_U_m(rs,T_N_REF);
    double ex=-Um/(omega_GW*C); if(ex<-700)ex=-700;
    return std::exp(ex);
}
double UQFFGWSuppression::compute_P_GW_UQFF(double m1,double m2,double r) const{
    double M=m1+m2;
    double P=compute_P_GW_standard(m1,m2,r);
    P*=compute_S_UA()*compute_S_SCm(M)*compute_S_TRZ()*compute_S_Um(M);
    for(auto& m:mods_) P*=m(m1+m2,r);
    return P;
}
double UQFFGWSuppression::compute_h_suppression(double m1,double m2,double r) const{
    double P_std=compute_P_GW_standard(m1,m2,r);
    double P_uqff=compute_P_GW_UQFF(m1,m2,r);
    if(P_std<1e-300) return 1.0;
    return std::sqrt(P_uqff/P_std);
}
void UQFFGWSuppression::simulate_P_sweep(double r0,double r1,double dr,double m1,double m2,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"r_m,P_std_W,P_UQFF_W,h_supp\n";}}
    for(double r=r0;r<=r1;r+=dr)
        *os<<r<<","<<compute_P_GW_standard(m1,m2,r)<<","<<compute_P_GW_UQFF(m1,m2,r)<<","<<compute_h_suppression(m1,m2,r)<<"\n";
    if(f.is_open())f.close();
}
void UQFFGWSuppression::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void UQFFGWSuppression::update_from_file(const std::string& fn){std::ifstream fi(fn);}
#ifdef STANDALONE_UQFFGWSUPPRESSION
int main(){
    UQFFGWSuppression gws;
    // GW150914: m1=36 Msun, m2=29 Msun, r~350 Rsun
    double m1=36*1.989e30, m2=29*1.989e30, r=350*6.96e8;
    std::cout<<"P_GW_std="<<gws.compute_P_GW_standard(m1,m2,r)<<" W\n";
    std::cout<<"P_GW_UQFF="<<gws.compute_P_GW_UQFF(m1,m2,r)<<" W\n";
    std::cout<<"h_suppression="<<gws.compute_h_suppression(m1,m2,r)<<"\n";
    gws.simulate_P_sweep(1e8,1e12,1e8,m1,m2,"gw_suppression_sweep.csv");
    return 0;
}
#endif
