#include "UQFFSuppressionEquationsHawking.h"
UQFFSuppressionEquationsHawking::UQFFSuppressionEquationsHawking(){}
double UQFFSuppressionEquationsHawking::compute_S1() const { return 1.0+F_TRZ; }
double UQFFSuppressionEquationsHawking::compute_S2() const { return 1.0-RHO_SCM/RHO_UA; }
double UQFFSuppressionEquationsHawking::compute_S3(double M) const {
    double T_H=_T_H(M); double r=_r_s(M);
    double Um=_U_m(r,T_N_REF);
    double ex=-Um/(K_B*T_H); if(ex<-700)ex=-700;
    return std::exp(ex);
}
double UQFFSuppressionEquationsHawking::compute_S_total(double M) const { return compute_S1()*compute_S2()*compute_S3(M); }
double UQFFSuppressionEquationsHawking::compute_T_UQFF(double M) const { return _T_H(M)*compute_S1()*compute_S2(); }
double UQFFSuppressionEquationsHawking::compute_L_UQFF(double M) const {
    double L=_L_H(M)*compute_S2()*compute_S3(M);
    for(auto& m:mods_) L*=m(M);
    return L;
}
double UQFFSuppressionEquationsHawking::compute_dT_dM_standard(double M) const { return -_T_H(M)/M; }
double UQFFSuppressionEquationsHawking::compute_dT_dM_UQFF(double M) const {
    return compute_dT_dM_standard(M)*compute_S1()*compute_S2();
}
void UQFFSuppressionEquationsHawking::sensitivity_sweep_rho_ratio(double r0,double r1,double dR,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"rho_ratio,S2,S_total\n";}}
    double M=8.55e36;
    for(double R=r0;R<=r1;R+=dR){
        double s2=1.0-1.0/R;
        double t_h=_T_H(M); double r2=_r_s(M);
        double Um=_U_m(r2,T_N_REF);
        double ex=-Um/(K_B*t_h); if(ex<-700)ex=-700;
        double s3=std::exp(ex);
        *os<<R<<","<<s2<<","<<compute_S1()*s2*s3<<"\n";
    }
    if(f.is_open())f.close();
}
void UQFFSuppressionEquationsHawking::add_mod(std::function<double(double)> m){mods_.push_back(std::move(m));}
void UQFFSuppressionEquationsHawking::update_from_file(const std::string& fn){std::ifstream fi(fn);}
#ifdef STANDALONE_UQFFSUPPRESSIONEQUATIONSHAWKING
int main(){
    UQFFSuppressionEquationsHawking se;
    double M=8.55e36;
    std::cout<<"S1="<<se.compute_S1()<<"\n";
    std::cout<<"S2="<<se.compute_S2()<<"\n";
    std::cout<<"S3="<<se.compute_S3(M)<<"\n";
    std::cout<<"S_total="<<se.compute_S_total(M)<<"\n";
    std::cout<<"T_UQFF="<<se.compute_T_UQFF(M)<<"\n";
    se.sensitivity_sweep_rho_ratio(2.0,20.0,1.0,"suppression_sweep.csv");
    return 0;
}
#endif
