#include "WhiteHoleRadiationUQFF.h"
WhiteHoleRadiationUQFF::WhiteHoleRadiationUQFF(){}

double WhiteHoleRadiationUQFF::compute_L_H(double M) const { return _L_H(M); }
double WhiteHoleRadiationUQFF::compute_T_H(double M) const { return _T_H(M); }
double WhiteHoleRadiationUQFF::compute_r_s_UQFF(double M) const { return _r_s(M)*(1.0-RHO_SCM/RHO_UA); }
double WhiteHoleRadiationUQFF::compute_L_WH_prime(double L_H) const { return L_H*(1.0+F_TRZ); }
double WhiteHoleRadiationUQFF::compute_L_WH_double_prime(double L_WH_prime) const { return L_WH_prime*(RHO_UA/RHO_SCM); }
double WhiteHoleRadiationUQFF::compute_L_WH_UQFF(double M, double r, double t) const {
    double L_H  = compute_L_H(M);
    double T_H  = compute_T_H(M);
    if(r<=0.0) r = _r_s(M);
    double Um   = _U_m(r, t>0?t:T_N_REF);
    double exp_arg = Um/(K_B*T_H); if(exp_arg>700) exp_arg=700;
    double L = L_H*(1.0+F_TRZ)*(RHO_UA/RHO_SCM)*std::exp(exp_arg);
    for(auto& m:mods_) L+=m(M);
    return L;
}
void WhiteHoleRadiationUQFF::simulate_over_M(double M0,double M1,double dM,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"M_kg,L_WH_UQFF_W\n";}}
    for(double M=M0;M<=M1;M+=dM) *os<<M<<","<<compute_L_WH_UQFF(M)<<"\n";
    if(f.is_open())f.close();
}
void WhiteHoleRadiationUQFF::add_mod(std::function<double(double)> mod){mods_.push_back(std::move(mod));}
void WhiteHoleRadiationUQFF::update_from_file(const std::string& fn){
    std::ifstream fi(fn); if(!fi.is_open())return; std::string l;
    while(std::getline(fi,l)){ auto e=l.find('='); if(e==std::string::npos)continue; }
}
#ifdef STANDALONE_WHITEHOLERADIATIONUQFF
int main(){
    WhiteHoleRadiationUQFF whr;
    double M=8.55e36; // Sgr A*
    std::cout<<"L_H="<<whr.compute_L_H(M)<<"\n";
    std::cout<<"T_H="<<whr.compute_T_H(M)<<"\n";
    std::cout<<"L_WH_UQFF="<<whr.compute_L_WH_UQFF(M)<<"\n";
    whr.simulate_over_M(1e30,1e40,1e38,"wh_radiation_sweep.csv");
    return 0;
}
#endif
