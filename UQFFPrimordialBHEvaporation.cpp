#include "UQFFPrimordialBHEvaporation.h"
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _UQFFPrimordialBHEvaporation_T_H(double M){
    return (UQFFPrimordialBHEvaporation::HBAR*UQFFPrimordialBHEvaporation::C*UQFFPrimordialBHEvaporation::C*UQFFPrimordialBHEvaporation::C)/
           (8.0*M_PI*UQFFPrimordialBHEvaporation::G*M*UQFFPrimordialBHEvaporation::K_B);
}
static inline double _UQFFPrimordialBHEvaporation_L_H(double M){
    return (UQFFPrimordialBHEvaporation::HBAR*std::pow(UQFFPrimordialBHEvaporation::C,6.0))/
           (15360.0*M_PI*UQFFPrimordialBHEvaporation::G*UQFFPrimordialBHEvaporation::G*M*M);
}
static inline double _UQFFPrimordialBHEvaporation_r_s(double M){
    return 2.0*UQFFPrimordialBHEvaporation::G*M/(UQFFPrimordialBHEvaporation::C*UQFFPrimordialBHEvaporation::C);
}
static inline double _UQFFPrimordialBHEvaporation_tau(double M){
    return 5120.0*M_PI*UQFFPrimordialBHEvaporation::G*UQFFPrimordialBHEvaporation::G*M*M*M/
           (UQFFPrimordialBHEvaporation::HBAR*std::pow(UQFFPrimordialBHEvaporation::C,4.0));
}
static inline double _UQFFPrimordialBHEvaporation_Um(double r,double t){
    double tn=t/UQFFPrimordialBHEvaporation::T_N_REF;
    return (UQFFPrimordialBHEvaporation::MU_J/r)*(1.0-std::exp(-UQFFPrimordialBHEvaporation::GAMMA*t*std::cos(M_PI*tn)));
}

UQFFPrimordialBHEvaporation::UQFFPrimordialBHEvaporation(){}
double UQFFPrimordialBHEvaporation::compute_dM_dt_std(double M) const{
    if(M<=0) return 0;
    return -HBAR*std::pow(C,4.0)/(15360.0*PI*G*G*M*M);
}
double UQFFPrimordialBHEvaporation::compute_dM_dt_UQFF(double M,double t) const{
    if(M<=0) return 0;
    double dM_std=compute_dM_dt_std(M);
    double T_H=_UQFFPrimordialBHEvaporation_T_H(M); double r=_UQFFPrimordialBHEvaporation_r_s(M);
    double Um=_UQFFPrimordialBHEvaporation_Um(r,t);
    double exp_arg=Um/(K_B*T_H); if(exp_arg>700)exp_arg=700;
    double dM=dM_std*(1.0-F_TRZ)*(RHO_SCM/RHO_UA)*std::exp(-exp_arg);
    for(auto& m:mods_) dM+=m(M,t);
    return dM;
}
double UQFFPrimordialBHEvaporation::compute_tau_std(double M) const{return _UQFFPrimordialBHEvaporation_tau(M);}
double UQFFPrimordialBHEvaporation::compute_tau_UQFF(double M) const{
    double tau=_UQFFPrimordialBHEvaporation_tau(M);
    double T_H=_UQFFPrimordialBHEvaporation_T_H(M); double r=_UQFFPrimordialBHEvaporation_r_s(M);
    double Um=_UQFFPrimordialBHEvaporation_Um(r,T_N_REF);
    double exp_arg=Um/(K_B*T_H); if(exp_arg>700)exp_arg=700;
    return tau/(1.0-F_TRZ)*(RHO_UA/RHO_SCM)*std::exp(exp_arg);
}
double UQFFPrimordialBHEvaporation::compute_M_formation(double t_form) const{
    // M_f = (4/3)*pi*(c*t_f/2)^3 * rho_rad(t_f)
    // rho_rad(t_f) ~ 3/(32*pi*G*t_f^2) (radiation dominated)
    double r_form=C*t_form/2.0;
    double rho_rad=3.0/(32.0*PI*G*t_form*t_form);
    return (4.0/3.0)*PI*r_form*r_form*r_form*rho_rad;
}
std::vector<std::pair<double,double>> UQFFPrimordialBHEvaporation::evolve_rk4(double M0,double t_end,double dt) const{
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
void UQFFPrimordialBHEvaporation::simulate_evaporation(double M0,double dt,const std::string& out) const{
    double tau_U=compute_tau_UQFF(M0);
    auto traj=evolve_rk4(M0,tau_U*2.0,dt);
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"t_s,M_kg\n";}}
    for(auto& p:traj) *os<<p.first<<","<<p.second<<"\n";
    if(ofs.is_open())ofs.close();
}
void UQFFPrimordialBHEvaporation::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void UQFFPrimordialBHEvaporation::update_from_file(const std::string& fn){
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_UQFFPRIMORDIALBHEVAPORATION
int main(){
    UQFFPrimordialBHEvaporation obj;
    double Mf=obj.compute_M_formation(1e-23);
    std::cout<<"M_formation="<<Mf<<" kg\n";
    std::cout<<"tau_UQFF="<<obj.compute_tau_UQFF(Mf)<<" s\n";
    obj.simulate_evaporation(1e12,1e50,"pbh_evaporation.csv");
    return 0;
}
#endif
