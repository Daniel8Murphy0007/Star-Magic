#include "VortexQuantization.h"
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _VortexQuantization_T_H(double M){
    return (VortexQuantization::HBAR*VortexQuantization::C*VortexQuantization::C*VortexQuantization::C)/
           (8.0*M_PI*VortexQuantization::G*M*VortexQuantization::K_B);
}
static inline double _VortexQuantization_L_H(double M){
    return (VortexQuantization::HBAR*std::pow(VortexQuantization::C,6.0))/
           (15360.0*M_PI*VortexQuantization::G*VortexQuantization::G*M*M);
}
static inline double _VortexQuantization_r_s(double M){
    return 2.0*VortexQuantization::G*M/(VortexQuantization::C*VortexQuantization::C);
}
static inline double _VortexQuantization_tau(double M){
    return 5120.0*M_PI*VortexQuantization::G*VortexQuantization::G*M*M*M/
           (VortexQuantization::HBAR*std::pow(VortexQuantization::C,4.0));
}
static inline double _VortexQuantization_Um(double r,double t){
    double tn=t/VortexQuantization::T_N_REF;
    return (VortexQuantization::MU_J/r)*(1.0-std::exp(-VortexQuantization::GAMMA*t*std::cos(M_PI*tn)));
}

VortexQuantization::VortexQuantization(){}
double VortexQuantization::compute_circulation(int n) const{
    return n*(2.0*PI*HBAR)/M_UA;
}
double VortexQuantization::compute_vortex_core(int n) const{
    double xi=HBAR/std::sqrt(2.0*M_UA*G_UA*N_UA);
    return xi*std::exp(-n*PI);
}
double VortexQuantization::compute_vortex_energy(int n,double R_outer,double xi) const{
    double kv=compute_circulation(n);
    if(xi<=0||R_outer<=xi) return 0.0;
    return RHO_UA*kv*kv/(4.0*PI)*std::log(R_outer/xi)*(RHO_UA/RHO_SCM);
}
double VortexQuantization::compute_omega_v_UQFF(double omega_v) const{
    double c_UA=std::sqrt(G_UA*N_UA/M_UA);
    return omega_v*(1.0+F_TRZ*c_UA/C);
}
double VortexQuantization::compute_magnus_force_UQFF(double v_s) const{
    // |F_Magnus| = rho_UA * kappa_1 * v_s per unit length
    double kv=compute_circulation(1);
    return RHO_UA*kv*v_s*(1.0+F_TRZ*(RHO_UA/RHO_SCM));
}
void VortexQuantization::simulate_lattice(int n_max,double R_outer,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"n,kv,a_core,E_v_UQFF\n";}}
    double xi=HBAR/std::sqrt(2.0*M_UA*G_UA*N_UA);
    for(int n=1;n<=n_max;n++){
        double kv=compute_circulation(n);
        double av=compute_vortex_core(n);
        double Ev=compute_vortex_energy(n,R_outer,xi);
        *os<<n<<","<<kv<<","<<av<<","<<Ev<<"\n";
    }
    if(ofs.is_open())ofs.close();
}
void VortexQuantization::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void VortexQuantization::update_from_file(const std::string& fn){
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_VORTEXQUANTIZATION
int main(){
    VortexQuantization obj;
    double xi=1e-5, R=1.0;
    for(int n=1;n<=5;n++){
        std::cout<<"n="<<n<<" kv="<<obj.compute_circulation(n)
                 <<" E_v="<<obj.compute_vortex_energy(n,R,xi)<<"\n";
    }
    obj.simulate_lattice(10,1.0,"vortex_lattice.csv");
    return 0;
}
#endif
