#include "AetherSuperfluidDynamics.h"
// â”€â”€â”€ UQFF helpers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
static inline double _AetherSuperfluidDynamics_T_H(double M){
    return (AetherSuperfluidDynamics::HBAR*AetherSuperfluidDynamics::C*AetherSuperfluidDynamics::C*AetherSuperfluidDynamics::C)/
           (8.0*M_PI*AetherSuperfluidDynamics::G*M*AetherSuperfluidDynamics::K_B);
}
static inline double _AetherSuperfluidDynamics_L_H(double M){
    return (AetherSuperfluidDynamics::HBAR*std::pow(AetherSuperfluidDynamics::C,6.0))/
           (15360.0*M_PI*AetherSuperfluidDynamics::G*AetherSuperfluidDynamics::G*M*M);
}
static inline double _AetherSuperfluidDynamics_r_s(double M){
    return 2.0*AetherSuperfluidDynamics::G*M/(AetherSuperfluidDynamics::C*AetherSuperfluidDynamics::C);
}
static inline double _AetherSuperfluidDynamics_tau(double M){
    return 5120.0*M_PI*AetherSuperfluidDynamics::G*AetherSuperfluidDynamics::G*M*M*M/
           (AetherSuperfluidDynamics::HBAR*std::pow(AetherSuperfluidDynamics::C,4.0));
}
static inline double _AetherSuperfluidDynamics_Um(double r,double t){
    double tn=t/AetherSuperfluidDynamics::T_N_REF;
    return (AetherSuperfluidDynamics::MU_J/r)*(1.0-std::exp(-AetherSuperfluidDynamics::GAMMA*t*std::cos(M_PI*tn)));
}

AetherSuperfluidDynamics::AetherSuperfluidDynamics(){}
double AetherSuperfluidDynamics::compute_m_UA() const{return M_UA;}
double AetherSuperfluidDynamics::compute_healing_length(double n_UA) const{
    if(n_UA<=0||G_UA<=0) return 1e300;
    return HBAR/std::sqrt(2.0*M_UA*G_UA*n_UA);
}
double AetherSuperfluidDynamics::compute_sound_speed(double n_UA) const{
    return std::sqrt(G_UA*n_UA/M_UA);
}
double AetherSuperfluidDynamics::compute_vortex_circulation() const{
    return (2.0*PI*HBAR)/M_UA; // h/m_UA
}
double AetherSuperfluidDynamics::compute_g_eff(double r,double M) const{
    double B_sys = 1e-4;  // DPM magnetic field [T]
    double mu_s = B_sys * r * r * r;
    double g_N = mu_s * G * M / (r * r);  // DPM-emergent: mu_s * grad(M_s/r)
    double c_UA=compute_sound_speed(N_UA_REF);
    double boost=1.0+(c_UA*c_UA/(C*C))*F_TRZ*(RHO_UA/RHO_SCM);
    double val=g_N*boost;
    for(auto& m:mods_) val+=m(r,M);
    return val;
}
double AetherSuperfluidDynamics::compute_GP_energy(double n_UA,double psi_sq) const{
    // E/V = -mu*|psi|^2 + g_UA/2 * |psi|^4
    double mu_UA=G_UA*n_UA;
    return -mu_UA*psi_sq + 0.5*G_UA*psi_sq*psi_sq;
}
void AetherSuperfluidDynamics::simulate_radial_profile(double r0,double r1,double dr,double M,
    const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"r_m,g_Newton,g_UQFF,n_UA_profile\n";}}
    for(double r=r0;r<=r1;r+=dr){
        double B_sys = 1e-4;  // DPM magnetic field [T]
        double mu_s = B_sys * r * r * r;
        double g_N = mu_s * G * M / (r * r);  // DPM-emergent: mu_s * grad(M_s/r)
        double g_u=compute_g_eff(r,M);
        // UA density: enhanced near BH by gravitational compression
        double n_local=N_UA_REF*(1.0+_AetherSuperfluidDynamics_r_s(M)/r*F_TRZ);
        *os<<r<<","<<g_N<<","<<g_u<<","<<n_local<<"\n";
    }
    if(ofs.is_open())ofs.close();
}
void AetherSuperfluidDynamics::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void AetherSuperfluidDynamics::update_from_file(const std::string& fn){
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_AETHERSUPERFLUIDDYNAMICS
int main(){
    AetherSuperfluidDynamics obj;
    std::cout<<"m_UA="<<obj.compute_m_UA()<<" kg\n";
    std::cout<<"xi_UA="<<obj.compute_healing_length(obj.N_UA_REF)<<" m\n";
    std::cout<<"c_UA="<<obj.compute_sound_speed(obj.N_UA_REF)<<" m/s\n";
    std::cout<<"kappa_v="<<obj.compute_vortex_circulation()<<" m^2/s\n";
    double M=8.55e36; double rs=2*6.6743e-11*M/(9e16);
    obj.simulate_radial_profile(rs,rs*100,rs,"aether_superfluid.csv");
    return 0;
}
#endif
