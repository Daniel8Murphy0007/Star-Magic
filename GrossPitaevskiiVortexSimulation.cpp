#include "GrossPitaevskiiVortexSimulation.h"
// ─── UQFF helpers ────────────────────────────────────────────────────────────
static inline double _GrossPitaevskiiVortexSimulation_T_H(double M){
    return (GrossPitaevskiiVortexSimulation::HBAR*GrossPitaevskiiVortexSimulation::C*GrossPitaevskiiVortexSimulation::C*GrossPitaevskiiVortexSimulation::C)/
           (8.0*M_PI*GrossPitaevskiiVortexSimulation::G*M*GrossPitaevskiiVortexSimulation::K_B);
}
static inline double _GrossPitaevskiiVortexSimulation_L_H(double M){
    return (GrossPitaevskiiVortexSimulation::HBAR*std::pow(GrossPitaevskiiVortexSimulation::C,6.0))/
           (15360.0*M_PI*GrossPitaevskiiVortexSimulation::G*GrossPitaevskiiVortexSimulation::G*M*M);
}
static inline double _GrossPitaevskiiVortexSimulation_r_s(double M){
    return 2.0*GrossPitaevskiiVortexSimulation::G*M/(GrossPitaevskiiVortexSimulation::C*GrossPitaevskiiVortexSimulation::C);
}
static inline double _GrossPitaevskiiVortexSimulation_tau(double M){
    return 5120.0*M_PI*GrossPitaevskiiVortexSimulation::G*GrossPitaevskiiVortexSimulation::G*M*M*M/
           (GrossPitaevskiiVortexSimulation::HBAR*std::pow(GrossPitaevskiiVortexSimulation::C,4.0));
}
static inline double _GrossPitaevskiiVortexSimulation_Um(double r,double t){
    double tn=t/GrossPitaevskiiVortexSimulation::T_N_REF;
    return (GrossPitaevskiiVortexSimulation::MU_J/r)*(1.0-std::exp(-GrossPitaevskiiVortexSimulation::GAMMA*t*std::cos(M_PI*tn)));
}

GrossPitaevskiiVortexSimulation::GrossPitaevskiiVortexSimulation(){}

std::vector<double> GrossPitaevskiiVortexSimulation::imaginary_time_step(
    const std::vector<double>& psi0,
    const std::vector<double>& r_grid,
    double M_bh, double dt, int n_steps) const{
    int N=(int)r_grid.size();
    std::vector<double> psi=psi0;
    double dr=r_grid.size()>1?r_grid[1]-r_grid[0]:1.0;
    for(int step=0;step<n_steps;step++){
        std::vector<double> dpsi(N,0.0);
        for(int i=1;i<N-1;i++){
            double r=r_grid[i];
            double lap=(psi[i+1]-2*psi[i]+psi[i-1])/(dr*dr)+
                       2.0/r*(psi[i+1]-psi[i-1])/(2*dr);
            double kin=-HBAR*HBAR/(2.0*M_UA)*lap;
            double V_grav=-G*M_bh*M_UA/r;
            double gp_int=G_UA*psi[i]*psi[i]*psi[i];
            double Um_val=_GrossPitaevskiiVortexSimulation_Um(r,1e8+step*dt);
            dpsi[i]=psi[i]-dt/HBAR*(kin/psi[i]+V_grav+gp_int+Um_val*psi[i]/(HBAR));
        }
        // boundary
        dpsi[0]=dpsi[1]; dpsi[N-1]=0.0;
        // normalise
        double norm=0.0;
        for(int i=0;i<N;i++) norm+=dpsi[i]*dpsi[i]*r_grid[i]*r_grid[i]*dr;
        norm=std::sqrt(norm*4*PI*N_UA);
        if(norm>0) for(int i=0;i<N;i++) psi[i]=dpsi[i]/norm*std::sqrt(N_UA);
    }
    return psi;
}
double GrossPitaevskiiVortexSimulation::compute_chemical_potential(
    const std::vector<double>& psi,
    const std::vector<double>& r_grid,
    double M_bh) const{
    double mu=0.0,norm=0.0;
    int N=(int)r_grid.size();
    double dr=N>1?r_grid[1]-r_grid[0]:1.0;
    for(int i=0;i<N;i++){
        double r=r_grid[i];
        double psi2=psi[i]*psi[i];
        double V=-G*M_bh*M_UA/r;
        mu+=(V+G_UA*psi2)*psi2*r*r*dr;
        norm+=psi2*r*r*dr;
    }
    return norm>0?mu/norm:0.0;
}
void GrossPitaevskiiVortexSimulation::simulate(double r_min,double r_max,int nr,double M_bh,
    int n_steps,const std::string& out) const{
    std::vector<double> r_grid(nr);
    double dr=(r_max-r_min)/(nr-1);
    for(int i=0;i<nr;i++) r_grid[i]=r_min+i*dr;
    // Initial guess: Gaussian
    std::vector<double> psi(nr);
    double r_mid=(r_min+r_max)/2.0;
    for(int i=0;i<nr;i++)
        psi[i]=std::exp(-(r_grid[i]-r_mid)*(r_grid[i]-r_mid)/(2.0*r_mid*r_mid));
    // Normalise
    double norm=0.0;
    for(int i=0;i<nr;i++) norm+=psi[i]*psi[i]*r_grid[i]*r_grid[i]*dr;
    if(norm>0){for(auto& p:psi) p/=std::sqrt(norm*4*PI/N_UA);}
    psi=imaginary_time_step(psi,r_grid,M_bh,1e-30,n_steps);
    double mu=compute_chemical_potential(psi,r_grid,M_bh);
    std::ostream* os=&std::cout; std::ofstream ofs;
    if(!out.empty()){ofs.open(out);if(ofs.is_open()){os=&ofs;*os<<"r_m,psi_sq,mu_eV\n";}}
    for(int i=0;i<nr;i++)
        *os<<r_grid[i]<<","<<psi[i]*psi[i]<<","<<mu/1.602e-19<<"\n";
    if(ofs.is_open())ofs.close();
}
void GrossPitaevskiiVortexSimulation::add_mod(std::function<double(double,double)> m){mods_.push_back(std::move(m));}
void GrossPitaevskiiVortexSimulation::update_from_file(const std::string& fn){
    std::ifstream fi(fn);if(!fi.is_open())return;
    std::string l;while(std::getline(fi,l)){auto e=l.find('=');if(e==std::string::npos)continue;}
}
#ifdef STANDALONE_GROSSPITAEVSKIIVORTEXSIMULATION
int main(){
    GrossPitaevskiiVortexSimulation obj;
    double M=8.55e36; double rs=2*6.6743e-11*M/(9e16);
    obj.simulate(rs,rs*50,100,M,50,"gp_vortex.csv");
    std::cout<<"GP simulation complete.\n";
    return 0;
}
#endif
