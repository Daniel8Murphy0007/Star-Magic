#include "UQFFPredictionsForLISA.h"
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _UQFFPredictionsForLISA_T_H(double M){
    return (UQFFPredictionsForLISA::HBAR*UQFFPredictionsForLISA::C*UQFFPredictionsForLISA::C*UQFFPredictionsForLISA::C)/
           (8.0*M_PI*UQFFPredictionsForLISA::G*M*UQFFPredictionsForLISA::K_B);
}
static inline double _UQFFPredictionsForLISA_L_H(double M){
    return (UQFFPredictionsForLISA::HBAR*std::pow(UQFFPredictionsForLISA::C,6.0))/
           (15360.0*M_PI*UQFFPredictionsForLISA::G*UQFFPredictionsForLISA::G*M*M);
}
static inline double _UQFFPredictionsForLISA_r_s(double M){
    return 2.0*UQFFPredictionsForLISA::G*M/(UQFFPredictionsForLISA::C*UQFFPredictionsForLISA::C);
}
static inline double _UQFFPredictionsForLISA_tau(double M){
    return 5120.0*M_PI*UQFFPredictionsForLISA::G*UQFFPredictionsForLISA::G*M*M*M/
           (UQFFPredictionsForLISA::HBAR*std::pow(UQFFPredictionsForLISA::C,4.0));
}
static inline double _UQFFPredictionsForLISA_Um(double r,double t){
    double tn=t/UQFFPredictionsForLISA::T_N_REF;
    return (UQFFPredictionsForLISA::MU_J/r)*(1.0-std::exp(-UQFFPredictionsForLISA::GAMMA*t*std::cos(M_PI*tn)));
}

UQFFPredictionsForLISA::UQFFPredictionsForLISA(){}
double UQFFPredictionsForLISA::compute_h_GR_SMBH(double Mc,double d,double f) const{
    return (4.0/d)*std::pow(G*Mc/(C*C),5.0/3.0)*std::pow(PI*f,2.0/3.0)/(C*C);
}
double UQFFPredictionsForLISA::compute_S_UA_LISA(double T_eff) const{
    double arg=RHO_UA*L_LISA/(K_B*T_eff); if(arg>700)arg=700;
    return std::max(0.0,1.0-arg);
}
double UQFFPredictionsForLISA::compute_h_UQFF_LISA(double Mc,double d,double f,double T_eff) const{
    double h_gr=compute_h_GR_SMBH(Mc,d,f);
    double S_UA=compute_S_UA_LISA(T_eff);
    double M_tot=Mc*std::pow(4.0,0.2); // approx equal mass
    double r=_UQFFPredictionsForLISA_r_s(M_tot); double T_H=_UQFFPredictionsForLISA_T_H(M_tot);
    double S_SCm=std::exp(-RHO_SCM*r/(K_B*T_H));
    double h=(1.0-F_TRZ)*S_UA*S_SCm*h_gr;
    for(auto& m:mods_) h+=m(f,T_eff);
    return h;
}
double UQFFPredictionsForLISA::compute_omega_GW_UQFF(double omega_GW_GR) const{
    return omega_GW_GR*std::pow(RHO_UA/RHO_CRIT,F_TRZ);
}
double UQFFPredictionsForLISA::compute_EMRI_rate_UQFF(double R_GR) const{
    return R_GR*(1.0+F_TRZ*(RHO_UA/RHO_SCM));
}
double UQFFPredictionsForLISA::compute_phase_mod(double tau_RD) const{
    return KAPPA*F_TRZ*tau_RD;
}
void UQFFPredictionsForLISA::simulate_SMBH_sweep(double Mc,double d,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"f_Hz,h_GR,h_UQFF\n";}}
    double T_eff=2.73; // CMB temperature as proxy
    double lf=std::log10(F_MIN_LISA),hf=std::log10(F_MAX_LISA);
    int N=200;
    for(int i=0;i<=N;i++){
        double f=std::pow(10.0,lf+(hf-lf)*i/N);
        *os<<f<<","<<compute_h_GR_SMBH(Mc,d,f)<<","<<compute_h_UQFF_LISA(Mc,d,f,T_eff)<<"\n";
    }
    if(ofs.is_open())ofs.close();
}
void UQFFPredictionsForLISA::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void UQFFPredictionsForLISA::update_from_file(const std::string& fn){
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_UQFFPREDICTIONSFORLISA
int main(){
    UQFFPredictionsForLISA obj;
    // Sgr A* + M87-type SMBH merger at cosmological distance
    double Mc=1e8*1.989e30; double d=1e27; // ~1 Gpc
    std::cout<<"Omega_GW_UQFF/GR ratio="<<obj.compute_omega_GW_UQFF(1.0)<<"\n";
    std::cout<<"EMRI rate boost="<<obj.compute_EMRI_rate_UQFF(100.0)<<" yr-1\n";
    obj.simulate_SMBH_sweep(Mc,d,"lisa_prediction.csv");
    return 0;
}
#endif
