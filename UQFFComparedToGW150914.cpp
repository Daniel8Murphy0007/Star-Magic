#include "UQFFComparedToGW150914.h"
UQFFComparedToGW150914::UQFFComparedToGW150914(){}
double UQFFComparedToGW150914::compute_chirp_mass(double m1,double m2) const{
    return std::pow(m1*m2,3.0/5.0)/std::pow(m1+m2,1.0/5.0);
}
double UQFFComparedToGW150914::compute_h_GR(double Mc,double d,double f) const{
    // h ~ 4/d * (G Mc/c²)^(5/3) * (pi f)^(2/3) / c  (leading order)
    double gmc=G*Mc/(C*C);
    return (4.0/d)*std::pow(gmc,5.0/3.0)*std::pow(PI*f,2.0/3.0)/C;
}
double UQFFComparedToGW150914::compute_S_SCm(double f) const{
    double lam=C/f;
    double ex=-RHO_SCM*lam; if(ex<-700)ex=-700;
    return std::exp(ex);
}
double UQFFComparedToGW150914::compute_h_UQFF(double Mc,double d,double f) const{
    double h=compute_h_GR(Mc,d,f);
    double omega=2.0*PI*f;
    double Um=MU_J/(_r_s(2*Mc));
    double ex=-Um*omega/(C*C); if(ex<-700)ex=-700;
    h*=(1.0-F_TRZ)*compute_S_SCm(f)*std::exp(ex);
    for(auto& m:mods_) h*=m(Mc,f);
    return h;
}
double UQFFComparedToGW150914::compute_df_dt(double Mc,double f) const{
    return (96.0/5.0)*std::pow(PI,8.0/3.0)*std::pow(G*Mc,5.0/3.0)/std::pow(C,5.0)*std::pow(f,11.0/3.0);
}
double UQFFComparedToGW150914::compute_dphi_UQFF(double f,double t) const{
    return 2.0*PI*f*t + KAPPA*F_TRZ*t;
}
void UQFFComparedToGW150914::simulate_inspiral(double f0,double f1,double df,double m1,double m2,double d,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"f_Hz,h_GR,h_UQFF\n";}}
    double Mc=compute_chirp_mass(m1,m2);
    for(double freq=f0;freq<=f1;freq+=df)
        *os<<freq<<","<<compute_h_GR(Mc,d,freq)<<","<<compute_h_UQFF(Mc,d,freq)<<"\n";
    if(f.is_open())f.close();
}
void UQFFComparedToGW150914::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void UQFFComparedToGW150914::update_from_file(const std::string& fn){std::ifstream fi(fn);}
#ifdef STANDALONE_UQFFCOMPAREDTOGW150914
int main(){
    UQFFComparedToGW150914 gw;
    double Mc=gw.compute_chirp_mass(M1_GW150914,M2_GW150914);
    std::cout<<"Chirp mass="<<Mc/1.989e30<<" Msun\n";
    std::cout<<"h_GR(150Hz)="<<gw.compute_h_GR(Mc,D_GW150914,150.0)<<"\n";
    std::cout<<"h_UQFF(150Hz)="<<gw.compute_h_UQFF(Mc,D_GW150914,150.0)<<"\n";
    gw.simulate_inspiral(10.0,300.0,1.0,M1_GW150914,M2_GW150914,D_GW150914,"gw150914_uqff.csv");
    return 0;
}
#endif
