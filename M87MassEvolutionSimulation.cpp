#include "M87MassEvolutionSimulation.h"
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _M87MassEvolutionSimulation_T_H(double M){
    return (M87MassEvolutionSimulation::HBAR*M87MassEvolutionSimulation::C*M87MassEvolutionSimulation::C*M87MassEvolutionSimulation::C)/
           (8.0*M_PI*M87MassEvolutionSimulation::G*M*M87MassEvolutionSimulation::K_B);
}
static inline double _M87MassEvolutionSimulation_L_H(double M){
    return (M87MassEvolutionSimulation::HBAR*std::pow(M87MassEvolutionSimulation::C,6.0))/
           (15360.0*M_PI*M87MassEvolutionSimulation::G*M87MassEvolutionSimulation::G*M*M);
}
static inline double _M87MassEvolutionSimulation_r_s(double M){
    return 2.0*M87MassEvolutionSimulation::G*M/(M87MassEvolutionSimulation::C*M87MassEvolutionSimulation::C);
}
static inline double _M87MassEvolutionSimulation_tau(double M){
    return 5120.0*M_PI*M87MassEvolutionSimulation::G*M87MassEvolutionSimulation::G*M*M*M/
           (M87MassEvolutionSimulation::HBAR*std::pow(M87MassEvolutionSimulation::C,4.0));
}
static inline double _M87MassEvolutionSimulation_Um(double r,double t){
    double tn=t/M87MassEvolutionSimulation::T_N_REF;
    return (M87MassEvolutionSimulation::MU_J/r)*(1.0-std::exp(-M87MassEvolutionSimulation::GAMMA*t*std::cos(M_PI*tn)));
}

M87MassEvolutionSimulation::M87MassEvolutionSimulation(){}
double M87MassEvolutionSimulation::compute_sound_speed_ISM(double T) const{
    // c_s = sqrt(gamma_adiab * k_B T / m_p) for fully ionized H
    return std::sqrt(5.0/3.0*K_B*T/(1.673e-27));
}
double M87MassEvolutionSimulation::compute_Mdot_Bondi_UQFF(double M,double rho_inf,double T_inf) const{
    double cs=compute_sound_speed_ISM(T_inf);
    double Mdot_B=4.0*PI*G*G*M*M*rho_inf/(cs*cs*cs);
    double rho_eff=rho_inf+RHO_UA-RHO_SCM;
    return Mdot_B*(rho_eff/rho_inf)*(1.0+F_TRZ);
}
double M87MassEvolutionSimulation::compute_dM_dt_evap_UQFF(double M,double t) const{
    if(M<=0) return 0;
    double dM_std=-HBAR*std::pow(C,4.0)/(15360.0*PI*G*G*M*M);
    double T_H=_M87MassEvolutionSimulation_T_H(M); double r=_M87MassEvolutionSimulation_r_s(M);
    double Um=_M87MassEvolutionSimulation_Um(r,t);
    double exp_arg=Um/(K_B*T_H); if(exp_arg>700)exp_arg=700;
    return dM_std*(1.0-F_TRZ)*(RHO_SCM/RHO_UA)*std::exp(-exp_arg);
}
double M87MassEvolutionSimulation::compute_jet_power_UQFF(double M) const{
    // P_BZ ~ B^2 r_s^2 c / 4  (simplified), B from equipartition
    double rs=_M87MassEvolutionSimulation_r_s(M);
    double rho_disk=RHO_ISM*1e6; // enhanced near BH
    double B_eq=std::sqrt(8.0*PI*rho_disk*C*C); // crude equipartition B
    if(B_eq>1e8) B_eq=1e8; // cap
    double Omega_H=C/(4.0*G*M/C/C)*C; // c/(4*Rg)*c
    double P_BZ=0.01*B_eq*B_eq*rs*rs*C*Omega_H*Omega_H/(4.0*PI*C);
    return P_BZ*(1.0+F_TRZ)*std::sqrt(RHO_UA/RHO_SCM);
}
double M87MassEvolutionSimulation::compute_dM_dt_total(double M,double t,double rho_inf,double T_inf) const{
    double dM_acc=compute_Mdot_Bondi_UQFF(M,rho_inf,T_inf);
    double dM_evap=compute_dM_dt_evap_UQFF(M,t);
    double P_jet=compute_jet_power_UQFF(M);
    double dM_jet=-P_jet/(C*C);
    double total=dM_acc+dM_evap+dM_jet;
    for(auto& m:mods_) total+=m(M,t);
    return total;
}
std::vector<std::pair<double,double>> M87MassEvolutionSimulation::evolve(double M0,double t_end,double dt,
    double rho_inf,double T_inf) const{
    std::vector<std::pair<double,double>> traj;
    double M=M0,t=0;
    while(t<t_end&&M>0){
        traj.push_back({t,M});
        double k1=compute_dM_dt_total(M,t,rho_inf,T_inf);
        double k2=compute_dM_dt_total(M+0.5*dt*k1,t+0.5*dt,rho_inf,T_inf);
        double k3=compute_dM_dt_total(M+0.5*dt*k2,t+0.5*dt,rho_inf,T_inf);
        double k4=compute_dM_dt_total(M+dt*k3,t+dt,rho_inf,T_inf);
        M+=dt/6.0*(k1+2*k2+2*k3+k4);
        t+=dt;
    }
    return traj;
}
void M87MassEvolutionSimulation::simulate_over_hubble(double dt,const std::string& out) const{
    auto traj=evolve(M_M87,T_HUBBLE,dt);
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"t_Gyr,M_Msun\n";}}
    for(auto& p:traj)
        *os<<p.first/3.156e16<<","<<p.second/M_SUN<<"\n";
    if(ofs.is_open())ofs.close();
}
void M87MassEvolutionSimulation::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void M87MassEvolutionSimulation::update_from_file(const std::string& fn){
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_M87MASSEVOLUTIONSIMULATION
int main(){
    M87MassEvolutionSimulation obj;
    double M=M_M87;
    std::cout<<"dM/dt_total="<<obj.compute_dM_dt_total(M,0)<<" kg/s\n";
    std::cout<<"Mdot_Bondi_UQFF="<<obj.compute_Mdot_Bondi_UQFF(M,RHO_ISM,T_ISM)<<" kg/s\n";
    obj.simulate_over_hubble(1e13,"m87_mass_evolution.csv");
    return 0;
}
#endif
