#include "LISAVsLIGOComparisons.h"
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _LISAVsLIGOComparisons_T_H(double M){
    return (LISAVsLIGOComparisons::HBAR*LISAVsLIGOComparisons::C*LISAVsLIGOComparisons::C*LISAVsLIGOComparisons::C)/
           (8.0*M_PI*LISAVsLIGOComparisons::G*M*LISAVsLIGOComparisons::K_B);
}
static inline double _LISAVsLIGOComparisons_L_H(double M){
    return (LISAVsLIGOComparisons::HBAR*std::pow(LISAVsLIGOComparisons::C,6.0))/
           (15360.0*M_PI*LISAVsLIGOComparisons::G*LISAVsLIGOComparisons::G*M*M);
}
static inline double _LISAVsLIGOComparisons_r_s(double M){
    return 2.0*LISAVsLIGOComparisons::G*M/(LISAVsLIGOComparisons::C*LISAVsLIGOComparisons::C);
}
static inline double _LISAVsLIGOComparisons_tau(double M){
    return 5120.0*M_PI*LISAVsLIGOComparisons::G*LISAVsLIGOComparisons::G*M*M*M/
           (LISAVsLIGOComparisons::HBAR*std::pow(LISAVsLIGOComparisons::C,4.0));
}
static inline double _LISAVsLIGOComparisons_Um(double r,double t){
    double tn=t/LISAVsLIGOComparisons::T_N_REF;
    return (LISAVsLIGOComparisons::MU_J/r)*(1.0-std::exp(-LISAVsLIGOComparisons::GAMMA*t*std::cos(M_PI*tn)));
}

LISAVsLIGOComparisons::LISAVsLIGOComparisons(){}
double LISAVsLIGOComparisons::compute_R_supp_LIGO(double M,double f,double t) const{
    double T_H=_LISAVsLIGOComparisons_T_H(M); double r=_LISAVsLIGOComparisons_r_s(M);
    double S_SCm=std::exp(-RHO_SCM*r/(K_B*T_H));
    double Um=_LISAVsLIGOComparisons_Um(r,t);
    double S_Um=std::exp(-Um*2.0*PI*f/(C*C));
    return (1.0-F_TRZ)*S_SCm*S_Um;
}
double LISAVsLIGOComparisons::compute_R_supp_LISA(double M,double f,double T_eff) const{
    double T_H=_LISAVsLIGOComparisons_T_H(M); double r=_LISAVsLIGOComparisons_r_s(M);
    double S_SCm=std::exp(-RHO_SCM*r/(K_B*T_H));
    double S_UA=std::max(0.0,1.0-RHO_UA*L_LISA/(K_B*T_eff));
    return (1.0-F_TRZ)*S_UA*S_SCm;
}
double LISAVsLIGOComparisons::compute_crossover_freq(double M,double T_eff) const{
    // Solve S_Um = S_UA => exp(-Um*2pif/c^2) = 1 - rho_UA*L/k_B T_eff
    double T_H=_LISAVsLIGOComparisons_T_H(M); double r=_LISAVsLIGOComparisons_r_s(M); double Um=_LISAVsLIGOComparisons_Um(r,1e8);
    double S_UA=std::max(1e-300,1.0-RHO_UA*L_LISA/(K_B*T_eff));
    if(Um<=0||S_UA<=0) return 0.0;
    return -std::log(S_UA)*C*C/(2.0*PI*Um);
}
double LISAVsLIGOComparisons::compute_Sn_UQFF(double Sn_GR,double f) const{
    return Sn_GR*(1.0+F_TRZ*(RHO_UA/RHO_SCM)*std::pow(f,-2.0/3.0));
}
void LISAVsLIGOComparisons::compare_sweep(double M,double d,
    double fl_lo,double fl_hi,double fa_lo,double fa_hi,
    const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"f_Hz,R_LIGO,R_LISA\n";}}
    // Sweep LIGO band
    for(double f=fl_lo;f<=fl_hi;f+=(fl_hi-fl_lo)/100.0)
        *os<<f<<","<<compute_R_supp_LIGO(M,f)<<","<<compute_R_supp_LISA(M,f)<<"\n";
    // Sweep LISA band
    for(int i=0;i<=100;i++){
        double f=fa_lo+i*(fa_hi-fa_lo)/100.0;
        *os<<f<<","<<compute_R_supp_LIGO(M,f)<<","<<compute_R_supp_LISA(M,f)<<"\n";
    }
    if(ofs.is_open())ofs.close();
}
void LISAVsLIGOComparisons::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void LISAVsLIGOComparisons::update_from_file(const std::string& fn){
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_LISAVSLIGOCOMPARISONS
int main(){
    LISAVsLIGOComparisons obj;
    double M=20*1.989e30, d=1e25;
    std::cout<<"R_LIGO@150Hz="<<obj.compute_R_supp_LIGO(M,150.0)<<"\n";
    std::cout<<"R_LISA@0.01Hz="<<obj.compute_R_supp_LISA(M,0.01)<<"\n";
    std::cout<<"f_cross="<<obj.compute_crossover_freq(M)<<" Hz\n";
    obj.compare_sweep(M,d,10,2000,1e-4,1,"lisa_vs_ligo.csv");
    return 0;
}
#endif
