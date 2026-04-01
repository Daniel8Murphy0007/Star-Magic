#include "UQFFAdvancementsAndTHzHoles.h"
UQFFAdvancementsAndTHzHoles::UQFFAdvancementsAndTHzHoles(){}
double UQFFAdvancementsAndTHzHoles::compute_f_THz(double T_c) const{
    return K_B*T_c/(2.0*PI*HBAR);
}
double UQFFAdvancementsAndTHzHoles::compute_f_Hawking_M(double M) const{
    // f_H = k_B T_H / (2 pi hbar)
    return K_B*_T_H(M)/(2.0*PI*HBAR);
}
double UQFFAdvancementsAndTHzHoles::compute_L_THz_UQFF(double M,double T_c) const{
    double f_THz_val=compute_f_THz(T_c);
    double f_H=compute_f_Hawking_M(M);
    if(f_H<1e-300) f_H=1e-300;
    double ratio_f4=std::pow(f_THz_val/f_H,4.0); if(ratio_f4>1e200) ratio_f4=1e200;
    double L=_L_H(M)*ratio_f4*(RHO_UA/RHO_SCM);
    for(auto& m:mods_) L*=m(M);
    return L;
}
double UQFFAdvancementsAndTHzHoles::compute_tau_RD_standard() const{
    // tau_RD_std ~ 1e13 yr in seconds
    return 1e13*3.156e7;
}
double UQFFAdvancementsAndTHzHoles::compute_tau_RD_UQFF() const{
    return compute_tau_RD_standard()*(RHO_UA/RHO_SCM)*(1.0+F_TRZ);
}
double UQFFAdvancementsAndTHzHoles::compute_fusion_suppression(double rho_core, double T_core) const{
    // Gamma_UQFF = sigma_fus * n^2 * (1 - rho_SCm/rho_UA)
    // sigma_fus ~ 1e-34 m^2 at T_core=1e7 K (pp chain)
    double sigma_fus=1e-34;
    double n=rho_core/M_E; // rough number density
    return sigma_fus*n*n*(1.0-RHO_SCM/RHO_UA);
}
double UQFFAdvancementsAndTHzHoles::compute_FAS(int n_papers) const{
    return n_papers*(1.0+F_TRZ)*std::sqrt(RHO_UA/RHO_SCM);
}
void UQFFAdvancementsAndTHzHoles::print_self_consistent_cycle() const{
    std::cout<<"UQFF Self-Consistent Cycle (PAPER_657-673):\n"
             <<"  KB7 (657) -> BH_Bounce (658) -> BH_WH_Trans (659) -> WH_Rad (660)\n"
             <<"  -> PBH_DM (661) -> Hawking (662) -> BH_Inv (663) -> WH_Stab (664)\n"
             <<"  -> Suppression (665) -> GW_Supp (666) -> BH_Stab (667) -> PBH_Stab (668)\n"
             <<"  -> GW150914 (669) -> Accretion (670) -> dM/dt (671) -> tau_evap (672)\n"
             <<"  -> THz_Holes (673) -> back to KB7\n"
             <<"  FAS = "<<compute_FAS(673)<<"\n";
}
void UQFFAdvancementsAndTHzHoles::simulate_THz_sweep(double M0,double M1,double dM,double T_c,const std::string& out) const{
    std::ostream* os=&std::cout; std::ofstream f;
    if(!out.empty()){f.open(out);if(f.is_open()){os=&f;*os<<"M_kg,f_Hawking_Hz,L_THz_UQFF_W\n";}}
    for(double M=M0;M<=M1;M+=dM)
        *os<<M<<","<<compute_f_Hawking_M(M)<<","<<compute_L_THz_UQFF(M,T_c)<<"\n";
    if(f.is_open())f.close();
}
void UQFFAdvancementsAndTHzHoles::add_mod(std::function<double(double)> m){mods_.push_back(std::move(m));}
void UQFFAdvancementsAndTHzHoles::update_from_file(const std::string& fn){std::ifstream fi(fn);}
#ifdef STANDALONE_UQFFADVANCEMENTSANDTHZHOLES
int main(){
    UQFFAdvancementsAndTHzHoles thz;
    std::cout<<"f_THz(100K)="<<thz.compute_f_THz()<<" Hz\n";
    std::cout<<"tau_RD_std="<<thz.compute_tau_RD_standard()<<"s\n";
    std::cout<<"tau_RD_UQFF="<<thz.compute_tau_RD_UQFF()<<"s\n";
    std::cout<<"FAS="<<thz.compute_FAS(673)<<"\n";
    thz.print_self_consistent_cycle();
    thz.simulate_THz_sweep(1e15,1e42,1e40,100.0,"thz_holes_sweep.csv");
    return 0;
}
#endif
