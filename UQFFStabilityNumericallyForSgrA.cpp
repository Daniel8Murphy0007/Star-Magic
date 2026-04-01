#include "UQFFStabilityNumericallyForSgrA.h"
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _UQFFStabilityNumericallyForSgrA_T_H(double M){
    return (UQFFStabilityNumericallyForSgrA::HBAR*UQFFStabilityNumericallyForSgrA::C*UQFFStabilityNumericallyForSgrA::C*UQFFStabilityNumericallyForSgrA::C)/
           (8.0*M_PI*UQFFStabilityNumericallyForSgrA::G*M*UQFFStabilityNumericallyForSgrA::K_B);
}
static inline double _UQFFStabilityNumericallyForSgrA_L_H(double M){
    return (UQFFStabilityNumericallyForSgrA::HBAR*std::pow(UQFFStabilityNumericallyForSgrA::C,6.0))/
           (15360.0*M_PI*UQFFStabilityNumericallyForSgrA::G*UQFFStabilityNumericallyForSgrA::G*M*M);
}
static inline double _UQFFStabilityNumericallyForSgrA_r_s(double M){
    return 2.0*UQFFStabilityNumericallyForSgrA::G*M/(UQFFStabilityNumericallyForSgrA::C*UQFFStabilityNumericallyForSgrA::C);
}
static inline double _UQFFStabilityNumericallyForSgrA_tau(double M){
    return 5120.0*M_PI*UQFFStabilityNumericallyForSgrA::G*UQFFStabilityNumericallyForSgrA::G*M*M*M/
           (UQFFStabilityNumericallyForSgrA::HBAR*std::pow(UQFFStabilityNumericallyForSgrA::C,4.0));
}
static inline double _UQFFStabilityNumericallyForSgrA_Um(double r,double t){
    double tn=t/UQFFStabilityNumericallyForSgrA::T_N_REF;
    return (UQFFStabilityNumericallyForSgrA::MU_J/r)*(1.0-std::exp(-UQFFStabilityNumericallyForSgrA::GAMMA*t*std::cos(M_PI*tn)));
}

UQFFStabilityNumericallyForSgrA::UQFFStabilityNumericallyForSgrA(){}
double UQFFStabilityNumericallyForSgrA::compute_omega_I_UQFF(double M) const{
    // omega_GR ~ -1/(tau_std) (decays at Hawking rate)
    double tau=_UQFFStabilityNumericallyForSgrA_tau(M);
    double T_H=_UQFFStabilityNumericallyForSgrA_T_H(M); double r=_UQFFStabilityNumericallyForSgrA_r_s(M);
    double Um=_UQFFStabilityNumericallyForSgrA_Um(r,T_N_REF);
    double omega_GR=-1.0/std::max(tau,1e-300);
    return omega_GR*(1.0+F_TRZ*(RHO_UA/RHO_SCM)*Um/(K_B*T_H));
}
double UQFFStabilityNumericallyForSgrA::compute_lyapunov_UQFF(double M) const{
    double tau=_UQFFStabilityNumericallyForSgrA_tau(M);
    double T_H=_UQFFStabilityNumericallyForSgrA_T_H(M); double r=_UQFFStabilityNumericallyForSgrA_r_s(M);
    double Um=_UQFFStabilityNumericallyForSgrA_Um(r,T_N_REF);
    double exp_arg=Um/(K_B*T_H); if(exp_arg>700)exp_arg=700;
    return -(RHO_SCM/RHO_UA)*std::exp(-exp_arg)/std::max(tau,1e-300);
}
double UQFFStabilityNumericallyForSgrA::compute_dM_dt_UQFF(double M,double t) const{
    if(M<=0) return 0.0;
    double L_H=_UQFFStabilityNumericallyForSgrA_L_H(M);
    double T_H=_UQFFStabilityNumericallyForSgrA_T_H(M); double r=_UQFFStabilityNumericallyForSgrA_r_s(M);
    double Um=_UQFFStabilityNumericallyForSgrA_Um(r,t);
    double exp_arg=Um/(K_B*T_H); if(exp_arg>700)exp_arg=700;
    double dM=-(1.0-F_TRZ)*L_H/(C*C)*(RHO_SCM/RHO_UA)*std::exp(-exp_arg);
    for(auto& m:mods_) dM+=m(M,t);
    return dM;
}
std::vector<std::pair<double,double>> UQFFStabilityNumericallyForSgrA::rk4_evolution(double M0,double t_end,double dt) const{
    std::vector<std::pair<double,double>> traj;
    double M=M0,t=0;
    while(t<t_end&&M>0){
        traj.push_back({t,M});
        double k1=compute_dM_dt_UQFF(M,t);
        double k2=compute_dM_dt_UQFF(M+0.5*dt*k1,t+0.5*dt);
        double k3=compute_dM_dt_UQFF(M+0.5*dt*k2,t+0.5*dt);
        double k4=compute_dM_dt_UQFF(M+dt*k3,t+dt);
        M+=dt/6.0*(k1+2*k2+2*k3+k4);
        t+=dt;
    }
    return traj;
}
bool UQFFStabilityNumericallyForSgrA::is_stable(double M) const{
    return compute_lyapunov_UQFF(M)<0.0;
}
void UQFFStabilityNumericallyForSgrA::simulate_stability(double t_end,double dt,const std::string& out) const{
    auto traj=rk4_evolution(M_SGRA,t_end,dt);
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"t_s,M_kg,stable\n";}}
    for(auto& p:traj)
        *os<<p.first<<","<<p.second<<","<<(is_stable(p.second)?1:0)<<"\n";
    if(ofs.is_open())ofs.close();
}
void UQFFStabilityNumericallyForSgrA::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void UQFFStabilityNumericallyForSgrA::update_from_file(const std::string& fn){
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_UQFFSTABILITYNUMERICALLYFORSGRA
int main(){
    UQFFStabilityNumericallyForSgrA obj;
    std::cout<<"Lyapunov SgrA*="<<obj.compute_lyapunov_UQFF(obj.M_SGRA)<<"\n";
    std::cout<<"Stable="<<obj.is_stable(obj.M_SGRA)<<"\n";
    obj.simulate_stability(1e60,1e55,"sgra_stability.csv");
    return 0;
}
#endif
