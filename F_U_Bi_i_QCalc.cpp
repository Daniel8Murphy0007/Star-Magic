// ======================================================================
// F_U_Bi_i_QCalc.cpp  -  Universal Buoyancy Equation Catalogue
// Author  : Daniel T. Murphy
// Project : Star-Magic / UQFF (Unified Quantum Field Framework)
// Version : v4.83   |   Date: 2026-03-21
// ----------------------------------------------------------------------
// SOURCES:
//   grok_share_c020496d9e.txt  (6,955 lines)
//   BB_C_Equations_04Sept2025.pdf
//   UQFF+Framework_Progress_Completion_Calibration_22Sept2025.pdf
//   UQFF_99_9_Complete_14Sept2025.pdf
// ----------------------------------------------------------------------
// PHYSICS:
//   F_U_Bi_i = Universal Buoyancy (outside->inside)
//              Governs ALL mass/non-mass actions at every scale.
//   F_U_Bi   = Universal Buoyancy (inside->outside) [negligibility factor]
//   Together -> suspension and perpetual motion.
//
//   MASTER INTEGRAL (grok_share_c020496d9e.txt, lines 1-60):
//   F_U_Bi_i = Int_0^{x2} [
//     -F0
//     + (me*c^2/r^2)*DPM_momentum*cos(theta)
//     + (GM/r^2)*DPM_gravity
//     + rho_vac[UA]*DPM_stability
//     + k_LENR*(omega_LENR/omega0)^2
//     + k_act*cos(omega_act*t)
//     + k_DE*L_X
//     + 2qB0V*sin(theta)*DPM_resonance*P_pol
//     + k_neutron*sigma_n
//     + k_rel*(Ecm_eff/Ecm)^2
//     + k_UV*L_UV  +  k_mm*L_mm*f_mm
//   ] dx
//
//   KEY NUMERICAL RESULTS:
//     F_U_Bi_i (LENR dominant)  ~= +2.11e208 N  (x2 ~= -1.35e172 m)
//     F_U_Bi_i (F_rel dominant) ~= -8.31e211 N
//     F_rel  =  4.31e33 N   (2024 LEP validation)
//     F_LENR =  1.56e36 N
//
// SECTIONS:
//   A  29 per-system g_UQFF equations (astrophysical systems)
//   B   6 Compressed UQFF backbone + Triadic Master equations
//   C  10 Sub-equations  (Um, [SSq], t_n, f_Ub, Ug2, vacuum series)
//   D  12 F_U_Bi_i master integral component force equations
//   E  79 unique F_UBii variant types  (BB_C_Equations catalogue)
//   F  68 unique Um   variant types
//   G  25 numerical solutions and calibration constants
//   H   7 Lambda-CDM / MOND comparison and validation notes
// ======================================================================

// Requires: C++20 or later
// Build:  cl /EHsc /std:c++20 F_U_Bi_i_QCalc.cpp /Fe:F_U_Bi_i_QCalc.exe
//         g++ -std=c++20 -O2 F_U_Bi_i_QCalc.cpp -o F_U_Bi_i_QCalc

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <numbers>   // C++20: std::numbers::pi etc.
#include <ranges>    // C++20: range views

// -------------------------------------------------------------------
// BuoyancyEquation: core data struct for every catalogued equation
// -------------------------------------------------------------------
struct BuoyancyEquation {
    int         id;               // Catalogue number
    std::string name;             // Short identifier
    std::string expression;       // Mathematical expression (ASCII notation)
    std::string context;          // Physical context / domain
    std::string validation_note;  // Key validation comment or cross-reference
    double      numerical_result; // Computed value (0.0 = not numerically resolved)
    std::string units;            // SI or natural units
    std::string system;           // Astrophysical system or domain
    std::string section;          // Section letter A-H
};

// -------------------------------------------------------------------
// printSection: display all equations in a named section
// -------------------------------------------------------------------
static void printSection(const std::string& title,
                         const std::vector<BuoyancyEquation>& eqs)
{
    const int W = 90;
    std::cout << "\n" << std::string(W, '=') << "\n";
    std::cout << "  SECTION: " << title
              << "  [" << eqs.size() << " equations]\n";
    std::cout << std::string(W, '-') << "\n";
    for (const auto& eq : eqs) {
        std::cout << "\n  [" << std::setw(5) << eq.id << "]  "
                  << eq.name << "  (" << eq.section << ")\n";
        std::cout << "    EXPR : " << eq.expression << "\n";
        std::cout << "    CTX  : " << eq.context    << "\n";
        if (!eq.validation_note.empty())
            std::cout << "    NOTE : " << eq.validation_note << "\n";
        if (eq.numerical_result != 0.0)
            std::cout << "    NUM  : "
                      << std::scientific << std::setprecision(4)
                      << eq.numerical_result
                      << "  " << eq.units << "\n";
    }
    std::cout << std::string(W, '-') << "\n";
}

// Forward declarations
std::vector<BuoyancyEquation> buildSectionA();
std::vector<BuoyancyEquation> buildSectionB();
std::vector<BuoyancyEquation> buildSectionC();
std::vector<BuoyancyEquation> buildSectionD();
std::vector<BuoyancyEquation> buildSectionE();
std::vector<BuoyancyEquation> buildSectionF();
std::vector<BuoyancyEquation> buildSectionG();
std::vector<BuoyancyEquation> buildSectionH();


// ======================================================
// SECTION A:  29 per-system g_UQFF equations
// UQFF backbone (ALL systems):
//   g_UQFF(r,t) = (G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)
//               + (Ug1+Ug2+Ug3+Ug4) + Lambda*c^2/3
//               + (hbar/sqrt(dx*dp))*Int[psi*H*psi dV]*(2pi/t_Hub)
//               + rho_f*V*g + (M_vis+M_DM)*(drho/rho + 3GM/r^3)
//               + SYSTEM_SPECIFIC_TERM
//   H(t,z) = H0*sqrt(0.3*(1+z)^3 + 0.7)
// ======================================================
std::vector<BuoyancyEquation> buildSectionA() {
    return {
        { 1, "g_student", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) [BASE: G*M_sun(t)/r(t)^2*(1+H0*t)]", "Student Guide Universe - cosmological base reference system", "Doc 1: calibration anchor; H0=67.4 km/s/Mpc; D_univ~92.77 Gly", 0.0, "m/s^2", "Student Guide Universe", "A" },
        { 2, "g_magnetar_SGR1745", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) + M_mag + D(t) + G*M_BH/r_BH^2", "Magnetar SGR 1745-2900 near Sgr A* - extreme B-field ~10^15 G", "Doc 2: M_mag=magnetic dipole torque; D(t)=drag; B_crit=4.4e13 T", 0.0, "m/s^2", "SGR 1745-2900 Magnetar", "A" },
        { 3, "g_sgrA_star", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) * sin(30deg) + G*M(t)^2/(c^4*r)*(dOmega/dt)^2", "Sagittarius A* SMBH - 4e6 Msun; galactic center accretion", "Doc 3: Kerr frame-dragging; a_spin~0.9; GWTC-4 compatible", 0.0, "m/s^2", "Sagittarius A*", "A" },
        { 4, "g_tapestry", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) + rho*v_wind^2", "Tapestry Starbirth Region - molecular cloud star formation; stellar winds", "Doc 4: v_wind~2000 km/s; SFR~10 Msun/yr", 0.0, "m/s^2", "Tapestry Starbirth Region", "A" },
        { 5, "g_westerlund1", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) + rho*v_wind^2 [cluster winds; implied]", "Westerlund 1 - massive young star cluster; ~10^5 Msun within 1 pc", "Doc 5: dense cluster; wind-feedback dominated", 0.0, "m/s^2", "Westerlund 1", "A" },
        { 6, "g_westerlund2", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) + rho*v_wind^2", "Westerlund 2 - star-forming giant HII region; eta Carinae neighbourhood", "Doc 6: FU_g1~2.43e-40 N; R(t)~-2.29e-41 N; FU_Bi~6.14e-32 N", 0.0, "m/s^2", "Westerlund 2", "A" },
        { 7, "g_pillars", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) * (1-E(t)) + rho*v_wind^2", "Pillars of Creation (M16 Eagle Nebula) - photoevaporation pillars", "Doc 7: E(t)=evaporation fraction; FU_g1~3.95e-41 N; FU_Bi~9.79e-33 N", 0.0, "m/s^2", "Pillars of Creation", "A" },
        { 8, "g_rings_relativity", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) * (1+L(t))", "Rings of Relativity - gravitational lensing system; Einstein ring", "Doc 8: L(t)=lensing amplification factor; arc geometry", 0.0, "m/s^2", "Rings of Relativity", "A" },
        { 9, "g_ngc2525", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) + G*M_BH/r_BH^2 - M_SN(t)", "NGC 2525 - SN Ia host galaxy; tidal interaction; AGN feedback", "Doc 10: M_SN(t)=supernova mass-loss term; Hubble SN Ia 2018", 0.0, "m/s^2", "NGC 2525", "A" },
        { 10, "g_ngc3603", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) * (1-P(t)) + rho*v_wind^2", "NGC 3603 - most luminous HII region in Milky Way; O/WR star cluster", "Doc 11: P(t)=pressure ionization fraction", 0.0, "m/s^2", "NGC 3603", "A" },
        { 11, "g_bubble_nebula", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) * (1+E(t)) + rho*v_exp^2", "Bubble Nebula (NGC 7635) - stellar wind blown bubble around O-star", "Doc 12: E(t)=expansion enhancement; v_exp~40 km/s", 0.0, "m/s^2", "Bubble Nebula NGC 7635", "A" },
        { 12, "g_antennae", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) * (1-M_coll(t)) + rho*v_sf^2", "Antennae Galaxies (NGC 4038/4039) - ongoing major merger; starburst", "Doc 14: M_coll(t)=collision deceleration; SFR~20 Msun/yr", 0.0, "m/s^2", "Antennae Galaxies", "A" },
        { 13, "g_horsehead", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) * (1-E(t)) + P_rad", "Horsehead Nebula (B33) - dark nebula silhouette against IC 434", "Doc 15: P_rad=radiation pressure from sigma Orionis", 0.0, "m/s^2", "Horsehead Nebula B33", "A" },
        { 14, "g_ngc1275", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) + F_BH + M_fil", "NGC 1275 Perseus Cluster - brightest cluster galaxy; AGN jet filaments", "Doc 16: F_BH=AGN feedback force; M_fil=filament infall momentum", 0.0, "m/s^2", "NGC 1275 Perseus Cluster", "A" },
        { 15, "g_hubble_udf", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) * (1+M_evo(t)) * (1-M_merge(t))", "Hubble Ultra Deep Field - z=0-7 galaxy population; cosmic web", "Doc 18: M_evo(t)=evolution modifier; M_merge(t)=merger fraction", 0.0, "m/s^2", "Hubble Ultra Deep Field", "A" },
        { 16, "g_ngc1792", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) * (1+M_sf(t)) + F_sn", "NGC 1792 - active star-forming spiral; supernovae feedback", "Doc 19: M_sf(t)=SF modulator; F_sn=supernova kick force", 0.0, "m/s^2", "NGC 1792", "A" },
        { 17, "g_sombrero", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) + G*M_BH/r_BH^2 + D_dust", "Sombrero Galaxy (M104) - prominent dust lane; 10^9 Msun SMBH", "Doc 20: D_dust=dust-lane pressure; bulge-disk decomposition", 0.0, "m/s^2", "Sombrero Galaxy M104", "A" },
        { 18, "g_saturn", "G*M_Sun/r_orbit^2 + T_ring + F_wind + (G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3)", "Saturn system - ring dynamics; magnetosphere; Enceladus plumes", "Doc 22: T_ring=ring torque; F_wind=solar wind; Cassini data", 0.0, "m/s^2", "Saturn System", "A" },
        { 19, "g_m16_eagle", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) * (1+M_sf(t)) - E_rad", "M16 Eagle Nebula - young open cluster; Pillars host nebula", "Doc 23: E_rad=EUV erosion rate from cluster stars", 0.0, "m/s^2", "M16 Eagle Nebula", "A" },
        { 20, "g_crab_nebula", "G*M_NS/r(t)^2 + F_wind + M_mag + (G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3)", "Crab Nebula (M1) - pulsar wind nebula; 33 ms pulsar; 1054 CE SN", "Doc 24: M_mag=pulsar magnetic torque; F_wind=pulsar nebula wind; dE/dt~5e38 erg/s", 0.0, "m/s^2", "Crab Nebula M1", "A" },
        { 21, "g_hydrogen_atom", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) * (1+P_term) + F_tech", "Hydrogen Atom - quantum scale; proton-electron Coulomb + nuclear", "Doc 27: P_term=pairing energy; F_tech=technology coupling term", 0.0, "m/s^2", "Hydrogen Atom", "A" },
        { 22, "g_hydrogen_resonance", "H_res=A_res*sin(2pi*f_res*t)+U_dp*SC_m*k_nuc+S_shell; A_res=k_A*Z*(A/A_H)*(1+d_pair); f_res=(E_bind/h)*(A_H/A)*(1+S_shell); U_dp=k*(A1*A2/f_dp^2)*cos(phi_dp); k_nuc=k0*(N/Z)*(1+d_pair); S_shell=0.1*(Z_magic+N_magic)", "Hydrogen Resonance - 7 sub-equations; nuclear + shell resonance", "Doc 28: SC_m~1; 7 coupled sub-equations for resonance frequency", 0.0, "m/s^2", "Hydrogen Resonance", "A" },
        { 23, "g_universe_diameter", "D_univ=2*D_p*(1+H(z)*t0)*(1+Lambda*c^2/(3*H0^2))*(1+(hbar/sqrt(dx*dp))*Int[psi*H*psi]/(G*M_tot))*(1+k*r_c^2)", "Universe Diameter - full UQFF cosmological scale equation", "Doc 26/29: D_univ~92.77 Gly (matches ~93 Gly observed); H0=67.4 km/s/Mpc", 0.0, "m/s^2", "Observable Universe", "A" },
        { 24, "g_ngc2264", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) * (1+M_sf(t)) + delta_SFR", "NGC 2264 - Christmas Tree Cluster + Cone Nebula; pre-main-sequence stars", "Doc 18b / SOURCE115: 26D polynomial node; SFR~0.5 Msun/yr", 0.0, "m/s^2", "NGC 2264", "A" },
        { 25, "g_tadpole", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) + G*dM_tidal/r_tidal^2", "Tadpole Galaxy (UGC 10214) - 280 kpc tidal tail from z~0.3 companion", "SOURCE115: tidal disruption; 26D resonance node", 0.0, "m/s^2", "Tadpole Galaxy UGC 10214", "A" },
        { 26, "g_mice_galaxies", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) * (1-M_coll(t)) + rho*v_tidal^2", "Mice Galaxies (NGC 4676 A/B) - equal-mass merger; counter-rotating tails", "SOURCE115: merger-induced starburst; 26D polynomial master", 0.0, "m/s^2", "Mice Galaxies NGC 4676", "A" },
        { 27, "g_m42_orion", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) + rho_PDR*v_shock^2 + P_OB_rad", "M42 Orion Nebula - photodissociation region; Trapezium OB cluster", "SOURCE115: PDR driven; Chandra X-ray data; 26D node", 0.0, "m/s^2", "M42 Orion Nebula", "A" },
        { 28, "g_carina_nebula", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) + rho*v_wind^2 + F_eta_car", "Carina Nebula - eta Carinae LBV; multiple HII regions; ~60 O-stars", "SOURCE115: F_eta_car=LBV Great Eruption analogue; SFR~3 Msun/yr", 0.0, "m/s^2", "Carina Nebula", "A" },
        { 29, "g_eso137", "(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3) - rho_ICM*v_ram^2", "ESO 137-001 - ram pressure stripping in Norma cluster; x-ray tail", "SOURCE115: v_ram~2000 km/s; ICM rho~1e-27 kg/m^3; Chandra 2025", 0.0, "m/s^2", "ESO 137-001", "A" },
    };
}

// ======================================================
// SECTION B:  Compressed UQFF backbone + Triadic Master Equations
// Source: grok_share_c020496d9e.txt (triadic masters ~lines 150-300)
// ======================================================
std::vector<BuoyancyEquation> buildSectionB() {
    return {
        // B-1  Compressed UQFF backbone (H(t,z) expansion)
        {  101, "H_tz_backbone",
           "H(t,z) = H0*sqrt(0.3*(1+z)^3 + 0.7)  [embedded in all g_UQFF]",
           "Hubble parameter evolution - Lambda-CDM shared backbone in UQFF",
           "H0=67.4 km/s/Mpc (Planck 2018); used by all 29 systems in Section A",
           67.4, "km/s/Mpc", "Cosmological backbone", "B" },

        // B-2  FU_g1 — Compressed gravity triadic master
        {  102, "FU_g1_compressed",
           "FU_g1 = Sum_k [ k_k*(f_UA1*f_SCm1*R_EB1)*(f_UA2*f_SCm2*R_EB2)/r^2*G_k"
           " + k4*rho_vac[SCm]*M_BH/r*exp(-alpha*t)*cos(pi*t_n)"
           "*(1+f_feedback)*exp(-[SSq]*n/26) ]",
           "Triadic compressed gravity - dual vacuum-density coupling across 26 layers",
           "Westerlund 2: FU_g1~2.43e-40 N (collapse driver); Pillars: ~3.95e-41 N",
           2.43e-40, "N", "Universal triadic gravity", "B" },

        // B-3  R(t) — Resonance master (26-layer oscillatory erosion)
        {  103, "R_t_resonance",
           "R(t) = Sum_{i=1}^{26} [R_Ug1i*cos(omega_Ug1i*t)"
           " + R_Ug2i*cos(omega_Ug2i*t)"
           " + R_Ug3i*cos(omega_Ug3i*t) + R_Ug4i*cos(omega_Ug4i*t)]",
           "26-layer resonance master - oscillatory erosion of buoyancy",
           "Westerlund 2: R(t)~-2.29e-41 N (oscillatory erosion)",
           -2.29e-41, "N", "26-layer resonance", "B" },

        // B-4  R_Ug1i — Sub-term for resonance amplitude and frequency
        {  104, "R_Ug1i_subterm",
           "R_Ug1i = F_Ug1i*(1+M_sf(t))*exp(-[SSq]*i/26);"
           " omega_Ug1i = 2*pi/(T_sf/i)*(1+[SSq])",
           "Individual layer amplitude and frequency for 26-layer resonance",
           "[SSq] calibrated: 0.5 (low-n), 5.26 (n=26 cosmic); T_sf=star-formation timescale",
           0.0, "N", "26-layer resonance sub-term", "B" },

        // B-5  FU_Bi — Universal Buoyancy triadic master
        {  105, "FU_Bi_buoyancy",
           "FU_Bi = Sum_k [ k_Ub_k*f_UA*f_SCm*R_EB/r^2"
           "*H_k(nu_THz, U_b, geometry_k)*f_Ub*exp(-(pi-t_n)) ]",
           "Triadic universal buoyancy (inside->outside); THz-frequency geometric coupling",
           "Westerlund 2: FU_Bi~6.14e-32 N (f_Ub*2.20e8); Pillars: ~9.79e-33 N",
           6.14e-32, "N", "Universal triadic buoyancy", "B" },

        // B-6  H_k and f_Ub sub-terms
        {  106, "H_k_f_Ub_subterms",
           "H_k = cos(phi)*f(nu_THz);"
           " f_Ub = k_Ub*Delta_k_eta*(rho_vac[UA]/rho_vac[SCm])*(V_little/V_big)",
           "Geometric THz filter H_k and buoyancy geometric ratio f_Ub",
           "f_Ub encodes little/big volume ratio and vacuum density contrast",
           0.0, "dimensionless", "Buoyancy sub-terms", "B" },
    };
}

// ======================================================
// SECTION C:  Sub-equations
// Source: grok_share_c020496d9e.txt (~lines 150-400)
// ======================================================
std::vector<BuoyancyEquation> buildSectionC() {
    return {
        // C-1  Um — Universal Magnetism (general form)
        {  201, "Um_general",
           "Um = Sum_j [mu_j(t,rho_vac[SCm])/r_j"
           " * (1-exp(-gamma*t)*cos(pi*t_n)) * phi^j]"
           " * P_SCm * E_react"
           " * (1+1e13*f_Heaviside) * (1+f_quasi) * exp(-[SSq])",
           "Universal Magnetism - vacuum-density weighted magnetic moment sum"
           " across all scales; Heaviside step for THz threshold",
           "Um~3.78e-6 J/m^3; gamma~5e-5 day^-1; phi~1.02;"
           " f_Heaviside=1 above 1 THz; f_quasi=quasi-particle correction",
           3.78e-6, "J/m^3", "Universal magnetism", "C" },

        // C-2  delta_n — Pseudo-monopole phase
        {  202, "delta_n_pseudomonopole",
           "delta_n = phi * (2*pi*n/6)",
           "Pseudo-monopole phase angle for n-th vacuum density layer",
           "phi~1.02; n=1..26; each layer contributes 60-degree phase offset",
           0.0, "rad", "Pseudo-monopole phase", "C" },

        // C-3  rho_vac layered series
        {  203, "rho_vac_UA_SCm_series",
           "rho_vac[UA]:[SCm] = rho_vac[UA]"
           " * (rho_vac[SCm]/rho_vac[UA])^n"
           " * exp(-[SSq]*n/26) * exp(-(pi-t_n))",
           "Layered vacuum density ratio across n=1..26 UQFF layers",
           "Bridges [UA] (universal awareness) and [SCm] (superconducting mass) states",
           0.0, "kg/m^3", "Layered vacuum density", "C" },

        // C-4  E_neutrino — neutrino energy from vacuum density
        {  204, "E_neutrino_vacuum",
           "E_neutrino ~ rho_vac[UA]:[SCm]"
           " * exp(-[SSq]*n/26 * exp(-(pi-t_n)))"
           " * (U_m / rho_vac[UA])",
           "Neutrino energy derived from layered vacuum density ratio",
           "E_neutrino~1.05e5 eV; vacuum-mediated neutrino mass generation",
           1.05e5, "eV", "Neutrino energy", "C" },

        // C-5  Decay Rate — vacuum-mediated decay
        {  205, "decay_rate_vacuum",
           "Decay_Rate ~ (rho_vac[SCm]/rho_vac[UA])"
           " * exp(-[SSq]*n/26 * exp(-(pi-t_n)))",
           "Vacuum-mediated particle decay rate - ratio of SCm to UA vacuum densities",
           "Decay_Rate~0.0583; scales with [SSq] layer index",
           0.0583, "s^-1", "Vacuum decay rate", "C" },

        // C-6  [SSq] — Calibrated suppression/squeezing factor
        {  206, "SSq_calibration",
           "[SSq] = log(rho_vac[SCm]/rho_vac[UA]) * n * exp(-(pi-t_n));"
           " n=1..26",
           "UQFF squeezing/suppression factor; calibrated against 96 astrophysical systems",
           "[SSq]=0.5 (low-n); [SSq]=5.26 (n=26 cosmic); Q_26([SSq]=5.26)~6.63e21",
           0.5, "dimensionless", "UQFF calibration constant", "C" },

        // C-7  t_n — Normalized UQFF time
        {  207, "t_n_normalized_time",
           "t_n = t/t_Hubble * (1 + H(z)*t0)",
           "Normalized UQFF time coordinate; scales with Hubble time",
           "t_Hubble~13.8 Gyr; t_n in [0,1] for cosmic epochs",
           0.0, "dimensionless", "Normalized time", "C" },

        // C-8  Vacuum Density Series (infinite sum)
        {  208, "vacuum_density_series",
           "V_density_series = Sum_{n=1}^{inf} (1/n^26) * [SSq]^n",
           "Infinite vacuum density series; converges for [SSq]<1/e",
           "26th-power convergence ensures rapid cutoff for large n",
           0.0, "kg/m^3", "Vacuum density series", "C" },

        // C-9  Buoyancy Harmonics H_m
        {  209, "buoyancy_harmonics_Hm",
           "H_m = Sum_{k=1}^{m} (1/k) * f_Ub",
           "Harmonic series for buoyancy amplitude; logarithmically divergent",
           "f_Ub=buoyancy geometric ratio; H_m~ln(m)*f_Ub for large m",
           0.0, "dimensionless", "Buoyancy harmonics", "C" },

        // C-10  Ug2 — Harmonic gravity series
        {  210, "Ug2_harmonic_gravity",
           "Ug2 = Sum_{m=1}^{inf} H_m * (1-exp(-[SSq]*m))"
           " * cos(omega_Ug2*t_n)",
           "Ug2 gravity component - buoyancy harmonic series with [SSq] suppression",
           "Oscillates at omega_Ug2; damped by [SSq] for large m; part of triadic sum",
           0.0, "m/s^2", "Ug2 harmonic gravity", "C" },
    };
}

// ======================================================
// SECTION D:  F_U_Bi_i Master Integral Component Forces
// Source: grok_share_c020496d9e.txt (lines ~1-80)
// Master integral broken into 12 additive force terms
// F_U_Bi_i(LENR dominant) ~= +2.11e208 N  (x2 ~= -1.35e172 m)
// F_U_Bi_i(F_rel dominant) ~= -8.31e211 N
// ======================================================
std::vector<BuoyancyEquation> buildSectionD() {
    return {
        // D-1  F0 — Reference force (offset)
        {  301, "F0_reference",
           "F0 = -F0  [subtracted; sets buoyancy zero-point]",
           "Reference force offset; subtracted in integral; sets absolute buoyancy baseline",
           "Calibrated per system; typically ~1 N for terrestrial reference",
           0.0, "N", "Reference force", "D" },

        // D-2  F_DPM_momentum — DPM momentum coupling
        {  302, "F_DPM_momentum",
           "F_mom = (me*c^2/r^2)*DPM_momentum*cos(theta)",
           "Electron rest-mass energy / r^2 * DPM momentum * cos(theta)",
           "DPM=Dynamic Phase Modulation; couples electron inertia to angular geometry",
           0.0, "N", "DPM momentum force", "D" },

        // D-3  F_DPM_gravity — DPM gravitational coupling
        {  303, "F_DPM_gravity",
           "F_grav = (G*M/r^2)*DPM_gravity",
           "Newtonian gravity * DPM stability modifier",
           "DPM_gravity bridges UQFF buoyancy to Newtonian; unity for low-density",
           0.0, "N", "DPM gravity force", "D" },

        // D-4  F_DPM_stability — Vacuum stability coupling
        {  304, "F_DPM_stability",
           "F_stab = rho_vac[UA]*DPM_stability",
           "Vacuum density [UA] * DPM stability factor",
           "Couples quantum vacuum to macroscopic buoyancy; dominant in low-density voids",
           0.0, "N", "Vacuum stability force", "D" },

        // D-5  F_LENR — Low Energy Nuclear Reaction coupling (DOMINANT)
        {  305, "F_LENR",
           "F_LENR = k_LENR * (omega_LENR/omega0)^2",
           "LENR resonance force - nuclear frequency ratio squared",
           "F_LENR~1.56e36 N; DOMINANT term; k_LENR from Widom-Larsen LENR coupling",
           1.56e36, "N", "LENR resonance", "D" },

        // D-6  F_act — Active modulation force
        {  306, "F_act",
           "F_act = k_act * cos(omega_act*t)",
           "Oscillatory activation force at frequency omega_act",
           "F_act~1e-6 N; small modulation; omega_act=activation resonance frequency",
           1.0e-6, "N", "Active modulation force", "D" },

        // D-7  F_DE — Dark Energy coupling via X-ray luminosity
        {  307, "F_DE",
           "F_DE = k_DE * L_X",
           "Dark energy coupling via X-ray luminosity proxy",
           "F_DE~1e-3 N; k_DE=dark energy-X-ray coupling constant",
           1.0e-3, "N", "Dark energy coupling", "D" },

        // D-8  F_res — Magnetic resonance force
        {  308, "F_res",
           "F_res = 2*q*B0*V*sin(theta)*DPM_resonance*P_pol",
           "Lorentz-resonance: charge*field*velocity*DPM_resonance*polarisation",
           "F_res~1.12e-25 N; DPM_resonance~1.67e7; P_pol=polarisation fraction",
           1.12e-25, "N", "Magnetic resonance force", "D" },

        // D-9  F_neutron — Neutron cross-section force
        {  309, "F_neutron",
           "F_neutron = k_neutron * sigma_n",
           "Neutron interaction force via total cross-section sigma_n",
           "F_neutron~1e7 N; k_neutron=neutron coupling constant",
           1.0e7, "N", "Neutron cross-section force", "D" },

        // D-10  F_rel — Relativistic particle collision force (LEP 2024)
        {  310, "F_rel",
           "F_rel = k_rel * (Ecm_eff/Ecm)^2",
           "Relativistic CM-energy ratio squared; 2024 LEP re-analysis validation",
           "F_rel=4.31e33 N; key cross-validation anchor for all F_UBii variants",
           4.31e33, "N", "Relativistic LEP force", "D" },

        // D-11  F_UV and F_mm — Multi-wavelength coupling forces
        {  311, "F_UV_F_mm",
           "F_UV = k_UV*L_UV  +  F_mm = k_mm*L_mm*f_mm;"
           " k_UV=k_mm=1e-30 N/W; f_mm=1.05",
           "UV (GALEX/Spitzer) + millimetre (ALMA) luminosity coupling forces",
           "k_UV=k_mm=1e-30 N/W; f_mm=1.05 (protoplanetary disk correction)",
           0.0, "N", "UV+mm luminosity coupling", "D" },

        // D-12  Derived F_hyb, F_hier, DeltaF, f_z_CGM
        {  312, "F_derived_composite",
           "F_hyb=P_pol*f_mm*omega0^-1;"
           " F_hier=Sum(v_i/c)^2*omega0^-1;"
           " DeltaF=Int(F_rel*exp(-t/tau)dt);"
           " f_z_CGM~1.46e-73",
           "Hybrid, hierarchical, time-integrated, and CGM redshift correction forces",
           "f_z_CGM~1.46e-73; DeltaF integrates F_rel over age tau;"
           " F_hier accounts for velocity hierarchy",
           1.46e-73, "dimensionless", "Derived composite forces", "D" },
    };
}

// ======================================================
// SECTION E:  79 unique F_UBii variant types
// Source: BB_C_Equations_04Sept2025.pdf
//         (extracted from grok_share_c020496d9e.txt ~lines 692-6955)
// General form:
//   F_UBii,TYPE = +/-F_rel * (physics_factor/E_LEP) * Q_wave * scaling
//   F_rel = 4.31e33 N  (2024 LEP)  Q_wave_mean = 3.97e4 J/m^3
// ======================================================
std::vector<BuoyancyEquation> buildSectionE() {
    return {
        { 710, "F_UBii_anyons", "F_UBii = +/-F_rel*(E_anyons/E_LEP)*Q_wave*g(r,t)*exp(-delta*c^2/(2*sigma^2))", "Anyon systems - fractional statistics in 2+1D; topological phases", "frac. quantum Hall analogue; [SSq] controls topological suppression; F_rel=4.31e33 N", 0.0, "N", "Topological condensed matter", "E" },
        { 721, "F_UBii_DE", "F_UBii = +/-F_rel*(rho_DE*c^2/E_LEP)*Q_wave*(8*pi*G*rho_tot/3)*(1+w(z))", "Dark energy density - late-universe acceleration; Lambda-CDM / quintessence", "w(z)=-1 for Lambda-CDM; measured w=-1.03+/-0.03 (Planck 2018)", 0.0, "N", "Cosmological dark energy", "E" },
        { 724, "F_UBii_inf", "F_UBii = +/-F_rel*(V(phi)/E_LEP)*Q_wave*3*H^2*exp(N)/(1+eps)", "Inflation - slow-roll scalar field V(phi); N e-folds; eps=slow-roll parameter", "Starobinsky / Higgs inflation compatible; eps=1-n_s/2~0.02", 0.0, "N", "Inflationary cosmology", "E" },
        { 727, "F_UBii_GW", "F_UBii = +/-F_rel*(rho_GW/E_LEP)*Q_wave*(32*pi*G*hdot^2/c^2)*exp(-t/tau)", "Gravitational waves - energy density rho_GW; quadrupole radiation", "tau=GW decay timescale; LIGO GWTC-4.0 compatible", 0.0, "N", "Gravitational wave astronomy", "E" },
        { 730, "F_UBii_merger", "F_UBii = +/-F_rel*(L_GW/E_LEP)*Q_wave*(32/5)*(G^(5/3)/c^5)*M^(5/3)*(pi*f)^(10/3)", "GW merger luminosity - Peters formula; inspiral power radiated", "Chirp mass M_c=(m1*m2)^(3/5)/(m1+m2)^(1/5); f=GW frequency", 0.0, "N", "Compact binary merger", "E" },
        { 739, "F_UBii_CR", "F_UBii = +/-F_rel*(E_max/E_LEP)*Q_wave*(4/3)*(v/c)^2*(c^2*E/lambda)*Z", "Cosmic rays - maximum energy Hillas criterion; DSA acceleration", "E_max=Z*e*B*R (Hillas); Z=charge; v/c=shock velocity ratio", 0.0, "N", "Cosmic ray acceleration", "E" },
        { 747, "F_UBii_MHD", "F_UBii = +/-F_rel*((E_M/t_eddy)/E_LEP)*Q_wave*(3/2)*(Re_m^(1/2)-1)", "MHD dynamo - magnetic energy over eddy turnover time; Re_m=magnetic Reynolds", "Re_m>1 for dynamo action; turbulent induction; Alfven wave spectrum", 0.0, "N", "MHD dynamo", "E" },
        { 753, "F_UBii_star", "F_UBii = +/-F_rel*((0.34*sigma_v^3/(G^2*rho*m*lnLambda))/E_LEP)*Q_wave*exp(-t/t_re)", "Star cluster relaxation timescale - velocity dispersion sigma_v; Coulomb log", "t_re=relaxation time; Spitzer 1987; lnLambda~10-15 for globular clusters", 0.0, "N", "Star cluster dynamics", "E" },
        { 770, "F_UBii_clus", "F_UBii = +/-F_rel*((gamma+1)/(2*gamma)*M^2/(gamma-1+2/M^2)/E_LEP)*Q_wave*(rho2/rho1)", "Cluster shock - Rankine-Hugoniot jump conditions; M=Mach number", "Bow shocks in galaxy clusters; rho2/rho1=(gamma+1)/(gamma-1) for strong shocks", 0.0, "N", "Galaxy cluster shock", "E" },
        { 772, "F_UBii_void", "F_UBii = +/-F_rel*(H*r/E_LEP)*Q_wave*(-3/5)*(Omega_m*a+Omega_Lambda)^(-3/2)*delta_v0", "Void expansion - linear theory; Hubble flow in underdense region", "delta_v0=initial void density contrast; voids ~10-30 Mpc; Omega_Lambda drives late", 0.0, "N", "Cosmic voids", "E" },
        { 774, "F_UBii_reion", "F_UBii = +/-F_rel*(ndot_gamma*eps_esc*f_star/E_LEP)*Q_wave*Int(ndot_gamma*dt)/(4*pi*n_H)^(1/3)", "Reionization - ionising photon budget; f_star=fraction that escape", "Epoch z~6-10; JWST early universe; eps_esc~0.1-0.2", 0.0, "N", "Cosmic reionization", "E" },
        { 776, "F_UBii_ISM", "F_UBii = +/-F_rel*(pi*c_s^2/(G*rho)/E_LEP)*Q_wave*(pi*c_s^2/(G*rho))", "ISM Jeans stability - thermal sound speed c_s; critical wavelength", "lambda_J=sqrt(pi*c_s^2/(G*rho)); c_s~0.2 km/s cold ISM; ~1 km/s warm", 0.0, "N", "Interstellar medium", "E" },
        { 782, "F_UBii_merg2", "F_UBii = +/-F_rel*((5/256*c^5/G^(5/3)*(1/(pi*f_i))^(8/3)*M^(5/3))/E_LEP)*Q_wave*Int(df/fdot)", "GW merger coalescence time - integrating inspiral phase from f_i", "Peters (1964) coalescence time; f_i=initial GW frequency", 0.0, "N", "GW coalescence time", "E" },
        { 800, "F_UBii_chirp", "F_UBii = +/-F_rel*((c^3/G*(5/96*pi^(-8/3)*f^(-11/3)*fdot)^(3/5))/E_LEP)*Q_wave*Int(df/fdot)", "Chirp mass from GW inspiral - frequency derivative fdot measurement", "M_chirp=(c^3/G)*(5/96*pi^{-8/3}*f^{-11/3}*fdot)^{3/5}; model-independent", 0.0, "N", "GW chirp mass", "E" },
        { 802, "F_UBii_qnm", "F_UBii = +/-F_rel*((c^3/(2*pi*G*M)*(0.3737+0.088*a_f))/E_LEP)*Q_wave*exp(-t/tau_qnm)", "BH quasi-normal modes - ringdown frequency; Kerr parameter a_f", "tau_qnm~M; l=2,m=2 dominant mode; Berti+2009 fits; GW150914 ringdown", 0.0, "N", "BH ringdown QNM", "E" },
        { 804, "F_UBii_bz", "F_UBii = +/-F_rel*((1/32*B^2*R_H^4*Omega_H^2/c)/E_LEP)*Q_wave*(a*c/(2*R_H))*(1+t_var)", "Blandford-Znajek jet power - frame-dragging; BH spin a; horizon radius R_H", "P_BZ=(kappa/16pi)*Phi_BH^2*Omega_H^2/c; AGN/microquasar jets", 0.0, "N", "BZ jet power", "E" },
        { 807, "F_UBii_reljet", "F_UBii = +/-F_rel*((Gamma(z)*c)/E_LEP)*Q_wave*(1/theta^2)*(z/z_acc)", "Relativistic jet - bulk Lorentz factor Gamma(z); opening angle theta", "z_acc=acceleration zone; VLBI proper motion; Gamma~5-20 AGN jets", 0.0, "N", "Relativistic jet", "E" },
        { 809, "F_UBii_puls", "F_UBii = +/-F_rel*((P/(2*Pdot))/E_LEP)*Q_wave*(2/3*B^2*R^6*Omega^4*sin^2(alpha)/c^3)", "Pulsar spin-down luminosity - period P; period derivative Pdot", "Magnetic dipole radiation; tau_c=P/(2*Pdot) characteristic age; Crab tau~1250 yr", 0.0, "N", "Pulsar spin-down", "E" },
        { 811, "F_UBii_glitch", "F_UBii = +/-F_rel*((I_s/I*nu0*(1-exp(-t/tau_q)))/E_LEP)*Q_wave*(Delta_Omega*exp(-t/tau_q))", "Pulsar glitch recovery - crustal superfluid vortex unpinning", "I_s/I=fraction of superfluid moment of inertia; tau_q=recovery timescale", 0.0, "N", "Pulsar glitch", "E" },
        { 813, "F_UBii_fire", "F_UBii = +/-F_rel*((eta_E/(M*c^2))/E_LEP)*Q_wave*(r/R0)", "GRB fireball - Lorentz factor eta_E=E/(Mc^2); above photosphere radius", "r>R_s for optically thin; Gamma_0~100-1000; GRB afterglow onset", 0.0, "N", "GRB fireball", "E" },
        { 815, "F_UBii_after", "F_UBii = +/-F_rel*(F_nu/E_LEP)*Q_wave*Gamma_B*eps_e^(1/2)*t^(-3/8); F_nu ~ nu^(-(p-1)/2) * t^(-3*(p-1)/4) for nu_m<nu<nu_c", "GRB afterglow synchrotron flux - electron spectrum index p~2.2", "nu_m=injection freq; nu_c=cooling freq; ISM or wind CBM", 0.0, "N", "GRB afterglow", "E" },
        { 817, "F_UBii_cmb", "F_UBii = +/-F_rel*((C_l)/E_LEP)*Q_wave; C_l = 2*pi*Int(k^2*dk*P(k)*|Delta_l^T(k)|^2); l~200", "CMB power spectrum C_l - acoustic peak structure; transfer function", "l~200 first peak; n_s~0.96; WMAP/Planck; k*eta_*~pi", 0.0, "N", "CMB power spectrum", "E" },
        { 819, "F_UBii_recomb", "F_UBii = +/-F_rel*((tau_recomb(z))/E_LEP)*Q_wave; tau = Int(n_e*sigma_T*c*(dt/dz)*dz); x_e ~ (T/n_b)^(3/4)*exp(-B_H/T)", "Recombination optical depth - Thomson scattering; ionisation fraction x_e", "z_rec~1100; delta_z~200; Saha equation; last scattering surface", 0.0, "dimensionless", "Recombination", "E" },
        { 821, "F_UBii_agn", "F_UBii = +/-F_rel*((mdot_out*v_out)/E_LEP)*Q_wave; mdot_out*v_out = f(v_out)*L_AGN/c", "AGN wind momentum flux - radiation momentum driving", "momentum boost factor f~1-20; quasar feedback; CfA survey", 0.0, "N", "AGN wind feedback", "E" },
        { 824, "F_UBii_bz2", "F_UBii = +/-F_rel*((kappa/(16*pi)*Phi_BH^2*Omega_BH^2/c)/E_LEP)*Q_wave*(a*c^3/(2*G*M))*[0.05 to 1]", "BZ variant with explicit spin-dependent efficiency factor", "kappa~0.044; Phi_BH=BH magnetic flux; a=spin parameter; eff=0.05-1", 0.0, "N", "BZ spin variant", "E" },
        { 826, "F_UBii_duty", "F_UBii = +/-F_rel*((1-exp(-t/tau_cool))/E_LEP)*(1+mdot_acc/mdot_Edd)^(-1)*Q_wave*(r^3/(G*M))", "AGN duty cycle - cooling time vs accretion Eddington ratio", "tau_cool=cooling time; duty~0.01-0.1 for L>0.1 L_Edd", 0.0, "N", "AGN duty cycle", "E" },
        { 829, "F_UBii_photo", "F_UBii = +/-F_rel*((eps*L_X*R_p^3/(G*M_p^2)*K(xi))/E_LEP)*Q_wave*(F_X/(G*M_p/R_p))", "Photoevaporation - X-ray flux on planet atmosphere; ionisation parameter xi", "K(xi)~correction; F_X=X-ray flux; eps=efficiency~0.1-0.3", 0.0, "N", "Planetary photoevaporation", "E" },
        { 831, "F_UBii_mig", "F_UBii = +/-F_rel*((f*(G*M_p)^2*Sigma*Omega^2*r_p^4/(H/r)^3)/E_LEP)*Q_wave; tau_mig=M_p/mdot_acc~1e6 yr", "Planet migration - Type I Lindblad-Corotation torque in disk", "Tanaka+2002; f~correction; tau_mig~1 Myr; disk-planet interaction", 0.0, "N", "Planet migration", "E" },
        { 853, "F_UBii_termv", "F_UBii = +/-F_rel*((sqrt(2*G*M*(1-Gamma)/r_launch))/E_LEP)*Q_wave*(tau*L/c)", "Terminal velocity - radiation-driven wind; Eddington factor Gamma", "Gamma=L/L_Edd; r_launch=launch radius; mass-loss by radiation pressure", 0.0, "m/s", "Stellar wind terminal velocity", "E" },
        { 855, "F_UBii_upar", "F_UBii = +/-F_rel*((Q_H/(4*pi*r^2*n_H*c))/E_LEP)*Q_wave*(Gamma/(n_H*c))", "Ionisation parameter U - ionising photon flux Q_H over n_H*c", "xi=L_X/(n_H*r^2) alternative form; U~1e-2 AGN NLR", 0.0, "dimensionless", "Ionisation parameter", "E" },
        { 857, "F_UBii_coup", "F_UBii = +/-F_rel*((0.5*mdot_w*v_w^2/(mdot_acc*c^2/10))/E_LEP)*Q_wave; coupling eff=[0.05 to 0.1]", "Wind coupling efficiency - kinetic luminosity to accretion power ratio", "eps_k~0.05-0.1; feedback models King+2015; AGN-ISM momentum coupling", 0.0, "dimensionless", "AGN coupling efficiency", "E" },
        { 877, "F_UBii_ps", "F_UBii = +/-F_rel*((dn/dM_PS)/E_LEP)*Q_wave; dn/dM = sqrt(2/pi)*rho0/M*(delta_c/sigma(M))*|dln(sigma)/dln(M)|*exp(-delta_c^2/(2*sigma^2))", "Press-Schechter mass function - collapse fraction vs halo mass M", "delta_c~1.686; sigma(M)=rms density fluctuation; Sheth-Tormen improved", 0.0, "Msun^-1 Mpc^-3", "Halo mass function", "E" },
        { 879, "F_UBii_sfe", "F_UBii = +/-F_rel*((f_b*mdot_halo/M_halo*H(z))/E_LEP)*(1+M_halo/M_crit)^(-1)*Q_wave; M_crit~1e12 Msun", "Star formation efficiency - baryon fraction; halo mass relative to M_crit", "SFR/M_halo peaks at M_crit~10^11.7 Msun; Moster+2013", 0.0, "yr^-1", "Star formation efficiency", "E" },
        { 881, "F_UBii_fb", "F_UBii = +/-F_rel*((eta*mdot_star*c^2)/E_LEP)*Q_wave; eta=radiative efficiency; mdot_out/SFR~1-10", "Stellar/AGN feedback - energy deposition rate; mass loading", "eta~0.1 for BH; eta~0.0015 for stellar; mdot_out~SFR for momentum-driven", 0.0, "W", "AGN/stellar feedback", "E" },
        { 899, "F_UBii_grow", "F_UBii = +/-F_rel*((delta)/E_LEP)*Q_wave; delta_ddot + 2*H*delta_dot = (3/2)*Omega_m*H^2*delta/a^3", "Density perturbation growth - Jeans instability in expanding universe", "D(a)=growth factor; D~a matter domination; D~const Lambda domination", 0.0, "dimensionless", "Structure growth", "E" },
        { 901, "F_UBii_haw", "F_UBii = +/-F_rel*((T_H)/E_LEP)*Q_wave*(kappa/(2*pi)); T_H = hbar*c^3/(8*pi*G*M*k_B)", "Hawking temperature - quantum BH radiation; surface gravity kappa", "T_H~6e-8 K (Msun BH); evaporation power P=hbar*c^6/(15360*pi*G^2*M^2)", 0.0, "K", "Hawking radiation temperature", "E" },
        { 903, "F_UBii_bhent", "F_UBii = +/-F_rel*((S_BH)/E_LEP)*Q_wave; S_BH = 4*pi*k_B*G*M^2/(hbar*c) = A/(4*l_pl^2)", "Bekenstein-Hawking entropy - horizon area in Planck units", "S_BH~1e77 k_B (Msun BH); information paradox; holographic bound", 0.0, "J/K", "BH entropy", "E" },
        { 905, "F_UBii_evapl", "F_UBii = +/-F_rel*((t_evap)/E_LEP)*Q_wave; t_evap = 5120*pi*G^2*M^3/(hbar*c^4); tau ~ M^3", "Hawking evaporation lifetime - cubic mass dependence", "t_evap~2e67 yr (Msun); primordial BH M~5e14 g evaporating now", 0.0, "s", "BH evaporation lifetime", "E" },
        { 907, "F_UBii_lqcf", "F_UBii = +/-F_rel*((H_LQC^2)/E_LEP)*Q_wave; H^2 = (8*pi*G*rho/3)*(1 - rho/rho_crit); rho_crit = 0.41*rho_Planck", "LQC Friedmann - loop quantum gravity bounce; quantum geometry correction", "Singularity replaced by bounce; rho_crit=0.41*rho_Pl~5e96 kg/m^3", 0.0, "s^-2", "LQC modified Friedmann", "E" },
        { 909, "F_UBii_lqcp", "F_UBii = +/-F_rel*((P_LQC(k))/E_LEP)*Q_wave; P(k) ~ k^(n_s-1)*(1+k/k_star)^(-alpha); k_star=bounce scale; alpha~4", "LQC power spectrum - suppressed at bounce scale k*; tilted at small k", "LQG imprint on CMB; n_s~0.96; alpha~4 from polymer quantisation", 0.0, "dimensionless", "LQC power spectrum", "E" },
        { 911, "F_UBii_bounc", "F_UBii = +/-F_rel*((t_bounce)/E_LEP)*Q_wave; t_bounce = sqrt(3/(16*pi*G*rho_crit)); rho_crit ~ rho_Planck", "LQC bounce time - timescale of quantum-gravity-driven contraction reversal", "t_bounce~Planck time at maximum density rho_crit~rho_Pl", 0.0, "s", "LQC bounce timescale", "E" },
        { 913, "F_UBii_ent", "F_UBii = +/-F_rel*((S_ent)/E_LEP)*Q_wave*exp(-t/tau_dec); S_ent = -Tr(rho_A*ln(rho_A))", "Quantum entanglement entropy - von Neumann entropy of subsystem A", "AdS/CFT: S_ent = A_RT/(4G); Page curve; decoherence time tau_dec", 0.0, "bits", "Quantum entanglement entropy", "E" },
        { 915, "F_UBii_ang", "F_UBii = +/-F_rel*((Jdot)/E_LEP)*Q_wave*T_B*exp(-t/tau_disk); Jdot = mdot*r^2*Omega", "Angular momentum transport - disk-star coupling; magnetic braking", "T_B=braking torque; tau_disk=disk lifetime; applies to accretion disks", 0.0, "N*m", "Angular momentum transport", "E" },
        { 917, "F_UBii_jetvel", "F_UBii = +/-F_rel*((v_j)/E_LEP)*Q_wave*(B/sqrt(4*pi*rho))*(1+t/tau_A); v_j ~= v_K*(r_A/r0)^(1/2)", "Jet velocity - Alfven speed at launch radius; Keplerian scaling", "r_A=Alfven radius; v_j~0.1-0.9c for AGN jets; magneto-centrifugal", 0.0, "m/s", "Jet launch velocity", "E" },
        { 919, "F_UBii_jshock", "F_UBii = +/-F_rel*((rho1*v1^2+P1)/E_LEP)*Q_wave*(v2/v1)*(gamma+1)/(gamma-1+2/M^2)", "J-type (jump) shock - Rankine-Hugoniot; discontinuous velocity jump", "strong shock: rho2/rho1=(gamma+1)/(gamma-1)~4; v2=v1/4; M>>1", 0.0, "Pa", "J-shock Rankine-Hugoniot", "E" },
        { 921, "F_UBii_cshock", "F_UBii = +/-F_rel*((rho_n*v_n*nu_ni*(v_i-v_n))/E_LEP)*Q_wave*v_s*exp(-z/L_d)", "C-type (continuous) shock - ion-neutral friction; magnetic support", "nu_ni=ion-neutral collision rate; L_d=dissipation scale; MHD shocks", 0.0, "N/m^3", "C-shock ion-neutral", "E" },
        { 923, "F_UBii_halo", "F_UBii = +/-F_rel*((dN_halo/dMdz)/E_LEP)*Q_wave*(1+z)^m; dN/dM~2*pi*sigma_M*sigma_m*|d delta_c/dz|*exp(-delta_c^2/(2*(sigma_m^2-sigma_M^2)))*dSigma_M/dM", "Halo formation rate - extended Press-Schechter; progenitor mass function", "m=merger number; Lacey-Cole 1993; conditional probability", 0.0, "Msun^-1 Mpc^-3 yr^-1", "Halo formation", "E" },
        { 925, "F_UBii_disk", "F_UBii = +/-F_rel*((mdot_star_disk)/E_LEP)*Q_wave; mdot_star = eps*M_gas/t_dyn; t_dyn = sqrt(3*pi/(32*G*rho)) for Q<1", "Disk star formation - gravitational instability Q<1; Kennicutt-Schmidt", "eps~0.01-0.03; t_dyn=dynamical time; Toomre Q criterion", 0.0, "Msun/yr", "Disk star formation", "E" },
        { 927, "F_UBii_burst", "F_UBii = +/-F_rel*((mdot_burst)/E_LEP)*Q_wave*(1+t/tau_merge); mdot_burst = mdot_gas_inflow*eps_burst", "Starburst inflow rate - merger-triggered gas funnelling", "eps_burst~0.5-1; tau_merge=merger timescale; Arp220 analogue", 0.0, "Msun/yr", "Starburst", "E" },
        { 929, "F_UBii_sedov", "F_UBii = +/-F_rel*((R_sedov)/E_LEP)*Q_wave; R(t) = (E*t^2/rho)^(1/5); d/dt(M*v) = 0; R ~ t^(2/5)", "Sedov-Taylor blast wave - adiabatic supernova remnant expansion", "Sedov 1946; E=explosion energy; rho=ISM density; valid phase 2", 0.0, "m", "Sedov-Taylor remnant", "E" },
        { 931, "F_UBii_dsa", "F_UBii = +/-F_rel*((s_index)/E_LEP)*Q_wave; N(E) ~ E^(-s); s = (4/3)*u_s^2/(r*dp); N(E) ~ E^(-2)*(u_s/c)", "Diffusive shock acceleration - power-law spectrum from first-order Fermi", "s~2 for strong shocks; u_s=shock speed; dp=momentum diffusion", 0.0, "GeV^-1 cm^-3", "DSA cosmic ray spectrum", "E" },
        { 941, "F_UBii_tov", "F_UBii = +/-F_rel*((dP_dr_TOV)/E_LEP)*Q_wave; dP/dr = -G*m(r)*rho(r)/r^2*(1+P/(rho*c^2))*(1+4*pi*r^3*P/(m*c^2))*(1-2*G*m/(r*c^2))^(-1)", "TOV equation - relativistic stellar structure for neutron stars", "Maximum mass ~2-2.5 Msun; GW170817 EOS constraints; stiffness", 0.0, "Pa/m", "Neutron star TOV", "E" },
        { 965, "F_UBii_nfwrot", "F_UBii = +/-F_rel*((v_NFW^2)/E_LEP)*Q_wave; v(r)^2 = 4*pi*G*rho_s*r_s^2*[ln(1+x)-x/(1+x)]/r; x=r/r_s", "NFW dark matter halo rotation curve - Navarro-Frenk-White profile", "rho_s=characteristic density; r_s=scale radius; flat rotation confirmed", 0.0, "m/s", "NFW rotation curve", "E" },
        { 967, "F_UBii_sidm", "F_UBii = +/-F_rel*((t_cross_SIDM)/E_LEP)*Q_wave; t_cross = 1e10*(rho/1e8Msun_kpc3)^(-1)*(sigma_m/1cm2_g)^(-1) yr; Gamma*t~1", "Self-interacting dark matter - core formation timescale; cross-section sigma/m", "sigma/m~1 cm^2/g for core-halo; Bullet Cluster limit <2 cm^2/g", 0.0, "yr", "SIDM core crossing", "E" },
        { 981, "F_UBii_IGM", "F_UBii = +/-F_rel*((E_IGM)/E_LEP)*Q_wave; E_IGM = k_B*T/m_p; E_total = (3/2)*n*k_B*T", "IGM thermal energy - post-reionization gas temperature T~10^4 K", "Ly-alpha forest; T~1e4 K at z~3; heating by photoionization", 0.0, "J/kg", "IGM thermal energy", "E" },
        { 983, "F_UBii_gal", "F_UBii = +/-F_rel*((dn_gal/dMdz)/E_LEP)*Q_wave; dn/dM ~ M_halo*H(z)^2*(2*delta_c/sigma(M))*exp(-delta_c^2/(2*sigma^2))/M^2", "Galaxy formation halo rate - abundance matching framework", "Moster+2013; stellar-to-halo mass relation; M_halo peak~1e12 Msun", 0.0, "Msun^-1 Mpc^-3 yr^-1", "Galaxy formation rate", "E" },
        { 985, "F_UBii_quant", "F_UBii = +/-F_rel*((P_quant(k))/E_LEP)*Q_wave*(1+f_NL/5); P(k) = H^2/(8*pi^2*eps*M_pl^2)", "Quantum inflation power spectrum - Harrison-Zel'dovich with non-Gaussianity", "f_NL=non-Gaussianity amplitude; f_NL<10 (Planck 2018); eps=slow-roll", 0.0, "dimensionless", "Inflationary power spectrum", "E" },
        { 989, "F_UBii_w", "F_UBii = +/-F_rel*((rho_DE*c^2*(1+w))/E_LEP)*Q_wave; d ln(rho_DE)/d ln(a) = -3*(1+w)", "Dark energy equation of state w = p/(rho*c^2)", "w=-1: Lambda; w=-1+/-0.04 DESI 2024; rhodot+3H(rho+p)=0", 0.0, "Pa", "Dark energy EOS", "E" },
        { 991, "F_UBii_BBN", "F_UBii = +/-F_rel*((eta_BBN)/E_LEP)*Q_wave; eta = n_b/n_gamma = 6e-10; reaction rates = H(T)", "Big Bang Nucleosynthesis baryon-to-photon ratio; He-4 abundance", "Y_p~0.247; D/H~2.5e-5; eta measured via CMB anisotropies", 6e-10, "dimensionless", "BBN", "E" },
        { 993, "F_UBii_relax", "F_UBii = +/-F_rel*((t_relax_star)/E_LEP)*Q_wave*exp(-t/t_relax); t_relax = 0.34*sigma_v^3/(G^2*rho*m*lnLambda)", "Dynamical relaxation - stellar cluster two-body relaxation time", "lnLambda~Coulomb logarithm ~10-15; cluster evaporation after t~t_relax", 0.0, "yr", "Stellar dynamical relaxation", "E" },
        { 995, "F_UBii_acc", "F_UBii = +/-F_rel*((L_Edd)/E_LEP)*Q_wave*exp(t/t_Sal); L_Edd = 4*pi*G*M*m_p*c/(sigma_T); mdot_Edd = L_Edd/(eps*c^2)", "Eddington accretion luminosity - radiation pressure limit", "t_Sal=45 Myr e-folding; L_Edd~1.26e38*(M/Msun) erg/s", 0.0, "W", "Eddington accretion", "E" },
        { 1046, "F_UBii_diff", "F_UBii = +/-F_rel*((D_CR)/E_LEP)*Q_wave; D_CR = 1e28*(E/10GeV)^(0.3 to 0.6) cm^2/s; D ~ (E/(B*deltaB))^(1/2)", "Cosmic ray diffusion coefficient - Bohm-like energy scaling", "Leaky-box model; grammage; B~3 muG ISM; measured via B/C ratio", 1e+28, "cm^2/s", "Cosmic ray diffusion", "E" },
        { 1066, "F_UBii_ng", "F_UBii = +/-F_rel*((f_NL)/E_LEP)*Q_wave; f_NL = (5/6)*(Gamma3-3*Gamma*Gammadot^2+2*Gammadot^3)/Gamma^4; B(k1,k2,k3) = (6/5)*f_NL*(P(k1)*P(k2)+perms)", "Non-Gaussianity bispectrum fNL - inflationary self-interaction", "Squeezed limit f_NL=5*(n_s-1)/12; f_NL<10 Planck 2018; CMB+LSS", 0.0, "dimensionless", "Primordial non-Gaussianity", "E" },
        { 1086, "F_UBii_wde", "F_UBii = +/-F_rel*((w_de)/E_LEP)*Q_wave; w = p/(rho*c^2) = -1 + (1/3)*d ln(rho_DE)/d ln(a); rhodot + 3*H*(rho+p) = 0", "w-dark energy density evolution - continuity equation solution", "rho ~ a^(-3*(1+w)); w=-1.03+/-0.03 DESI 2024", -1.03, "dimensionless", "w-dark energy", "E" },
        { 1106, "F_UBii_disk2", "F_UBii = +/-F_rel*((mdot_star_disk2)/E_LEP)*Q_wave; mdot_star = eps*M_gas/t_dyn; Q = sigma_v*kappa/(pi*G*Sigma) < 1; t_dyn = sqrt(3*pi/(32*G*rho))", "Disk SF Toomre Q variant - gravitational instability with epicyclic freq kappa", "Q<1 unstable; kappa=epicycle frequency; Sigma=surface density; Kennicutt 1998", 0.0, "Msun/yr", "Disk SF Toomre variant", "E" },
        { 1746, "F_UBii_reh", "F_UBii = +/-F_rel*((T_reh)/E_LEP)*Q_wave; T_reh = (30*V_end/(pi^2*g_star))^(1/4)*exp(-3*N_reh/4)", "Reheating temperature after inflation - inflaton decay; N_reh e-folds", "g_star~100-200 relativistic dof; T_reh~1e14 GeV high; ~MeV low", 0.0, "GeV", "Reheating temperature", "E" },
        { 1805, "F_UBii_eta", "F_UBii = +/-F_rel*((eta_baryons)/E_LEP)*Q_wave; eta = n_b/n_gamma = 6e-10; BBN rates = H(T)", "Baryon-to-photon ratio - baryon asymmetry; measured via D/H and CMB", "Baryogenesis; Sakharov conditions; eta=6.1e-10 (Planck 2018)", 6.1e-10, "dimensionless", "Baryon asymmetry", "E" },
        { 1825, "F_UBii_termv2", "F_UBii = +/-F_rel*((v_inf)/E_LEP)*Q_wave*(tau*L/c); v_inf = sqrt(2*G*M*(1-Gamma)/r_launch)", "Line-driven wind terminal velocity - OB star mass loss; Eddington-modified", "v_inf~1-3*v_esc; Gamma=L/L_Edd; CAK 1975; CMFGEN models", 0.0, "m/s", "Line-driven wind", "E" },
        { 1845, "F_UBii_metal", "F_UBii = +/-F_rel*((Z_metal)/E_LEP)*Q_wave; dZ/dt = Y*SFR/mdot_gas - Z*mdot_out/M_gas; Z ~ 0.1 Z_sun", "Metallicity evolution - closed box + outflow; yield Y; loading factor", "Effective yield Y*(1+eta_out)^(-1); mass-metallicity relation", 0.0, "Z_sun", "Metallicity evolution", "E" },
        { 1865, "F_UBii_rev", "F_UBii = +/-F_rel*((l_rev)/E_LEP)*Q_wave; l_rev = (alpha_dyn*eta)^(1/2)*l_force; dB/dt = curl(alpha*B - eta*curl*B)", "Magnetic field reversal / turbulent dynamo - mean-field alpha-Omega dynamo", "alpha=helicity; eta=resistivity; reversal scale l_rev; solar ~11 yr cycle", 0.0, "m", "Magnetic dynamo reversal", "E" },
        { 1885, "F_UBii_ent2", "F_UBii = +/-F_rel*((S_ent_page)/E_LEP)*Q_wave; S_ent = -Tr(rho_A*ln(rho_A)); Page curve: S increases then decreases", "Entanglement entropy Page curve variant - BH information paradox", "Penington+2019; island rule; replica wormholes; unitarity restored", 0.0, "bits", "BH Page curve", "E" },
        { 1945, "F_UBii_voidden", "F_UBii = +/-F_rel*((delta_v)/E_LEP)*Q_wave; delta_v(a) = -3/5*(Omega_m*a+Omega_Lambda)^(-3/2)*delta_v0; ~ a^(-1)", "Void density contrast evolution - linear theory; void grows as a^(-1)", "delta_v0=initial contrast; underdense regions evacuate faster", 0.0, "dimensionless", "Void density", "E" },
        { 1965, "F_UBii_jeans", "F_UBii = +/-F_rel*((lambda_J)/E_LEP)*Q_wave; lambda_J = sqrt(pi*c_s^2/(G*rho)); dP = -rho*grad(Phi)", "Jeans instability length - thermal support vs self-gravity", "lambda_J~0.2 pc cold HI; ~1 pc warm; M_J=collapse threshold mass", 0.0, "m", "Jeans instability", "E" },
        { 1985, "F_UBii_relax2", "F_UBii = +/-F_rel*((t_relax2)/E_LEP)*Q_wave*exp(-t/t_relax); t_relax = 0.34*sigma_v^3/(G^2*rho*m*lnLambda)", "Dynamical relaxation II - identical to Spitzer but applied to galaxy clusters", "lnLambda~30 for clusters; t_relax~Gyr; drives isothermal cores", 0.0, "yr", "Cluster dynamical relaxation", "E" },
        { 2005, "F_UBii_conv", "F_UBii = +/-F_rel*((t_conv)/E_LEP)*Q_wave; t_conv = H_p/v_conv; v_conv ~= (alpha^2*g*deltaT*H_p/(4*T))^(1/3)", "Convective turnover time - mixing-length theory; alpha=MLT parameter", "H_p=pressure scale height; alpha~1.8 solar calibration; stellar interiors", 0.0, "s", "Convective turnover", "E" },
        { 933, "F_UBii_chirp2", "F_UBii = +/-F_rel*((M_chirp2)/E_LEP)*Q_wave; M_chirp = (m1*m2)^(3/5)/(m1+m2)^(1/5); GW strain h ~ M_chirp^(5/3)/d", "Chirp mass variant - strain amplitude scaling; source distance d", "h~1e-21 at 100 Mpc for 30+30 Msun binary; LIGO sensitivity limit", 0.0, "Msun", "GW chirp mass variant", "E" },
        { 1058, "F_UBii_dyn", "F_UBii = +/-F_rel*((B_dyn)/E_LEP)*Q_wave; dB/dt = curl(v x B - eta*curl*B); growth rate ~ v_turb/l_turb - eta/l_turb^2", "Turbulent dynamo field growth - kinematic vs resistive competition", "Minimum turbulent Rm > Rm_crit~100; field saturates at B^2~rho*v_turb^2", 0.0, "T", "Turbulent dynamo B-field", "E" },
        { 1706, "F_UBii_deb", "F_UBii = +/-F_rel*((rho_rad)/E_LEP)*Q_wave; rho_rad = pi^2*k_B^4*T^4/(30*hbar^3*c^5)*g_star", "Radiation energy density - Stefan-Boltzmann with relativistic dof count g*", "g_star=106.75 SM at T>200 GeV; rho_rad~a^(-4); reheating endpoint", 0.0, "J/m^3", "Radiation energy density", "E" },
        { 1786, "F_UBii_vir", "F_UBii = +/-F_rel*((sigma_v_vir)/E_LEP)*Q_wave; sigma_v^2 = G*M/(3*r); (virial theorem for gas)", "Virial velocity dispersion - kinetic-gravitational energy equipartition", "2*KE + PE = 0; sigma_v~200 km/s galaxy clusters; ~5 km/s globulars", 0.0, "m/s", "Virial velocity dispersion", "E" },
        { 992, "F_UBii_pec", "F_UBii = +/-F_rel*((v_pec)/E_LEP)*Q_wave; v_pec = -(f*H/(4*pi))*grad^-1*(delta); f(Omega_m) ~ Omega_m^(0.55)", "Peculiar velocity - linear growth rate; density-velocity correspondence", "RSD; beta=f*sigma8; 2M++ peculiar velocity survey; Tully-Fisher", 0.0, "m/s", "Peculiar velocity", "E" },
        { 788, "F_UBii_merg3", "F_UBii = +/-F_rel*((L_GW_merg)/E_LEP)*Q_wave; L_GW = (32/5)*(G^4/c^5)*(m1*m2*(m1+m2)/r^5); tau_merge = 12/85*(a^4/G^3*m1*m2*(m1+m2))", "GW merger power and coalescence time (quadrupole formula full form)", "Peters 1964; a=semi-major axis; tau shrinks as a^4; GW inspiral", 0.0, "W", "GW merger luminosity", "E" },
    };
}

// ======================================================
// SECTION F:  68 unique Um (Universal Magnetism) variant types
// Source: BB_C_Equations_04Sept2025.pdf
//         (grok_share_c020496d9e.txt ~lines 693-2100)
// General form (base):
//   Um = Sum_j[mu_j(t,rho_vac[SCm])/r_j*(1-exp(-gamma*t)*cos(pi*t_n))*phi^j]
//       * P_SCm * E_react * (1+1e13*f_Heaviside) * (1+f_quasi) * exp(-[SSq])
//       * TYPE_SPECIFIC_FACTOR
// ======================================================
std::vector<BuoyancyEquation> buildSectionF() {
    return {
        { 693, "Um_general", "Um = Sum_j[mu_j(t,rho_vac[SCm])/r_j*(1-exp(-gamma*t)*cos(pi*t_n))*phi^j]*P_SCm*E_react*(1+1e13*f_Heaviside)*(1+f_quasi)*exp(-[SSq])", "Universal Magnetism general form - all-scale vacuum magnetic moment summation", "gamma~5e-5 day^-1; phi~1.02; Um~3.78e-6 J/m^3; base of all 68 variants", 3.78e-06, "J/m^3", "Universal magnetism", "F" },
        { 724, "Um_cosmo", "Um_cosmo = Um_general * (k*c^2/a^2) / H^2", "Cosmological Um - scaled by curvature/Hubble ratio; FRW correction", "k=curvature; a=scale factor; H=Hubble parameter; applies to cosmic epochs", 0.0, "J/m^3", "Cosmological Um", "F" },
        { 746, "Um_reh", "Um_reh = Um_general * (30*V_end/(pi^2*g_star))^(1/4) * exp(-3*N_reh/4)", "Reheating Um - inflaton decay epoch magnetic moment; reheat temperature factor", "V_end=inflaton potential at end; g_star~106; N_reh e-folds of reheating", 0.0, "J/m^3", "Reheating Um", "F" },
        { 748, "Um_MHD", "Um_MHD = Um_general * 4*pi*rho*v_turb * M_A", "MHD turbulent Um - magnetic energy from turbulent velocity * Alfven Mach M_A", "M_A=v_turb/v_A; v_A=Alfven speed; sub-Alfvenic turbulence M_A<1", 0.0, "J/m^3", "MHD turbulent Um", "F" },
        { 770, "Um_BBN", "Um_BBN = Um_general * (3/(32*pi*G*rho_rad)); t_BBN ~ 180 s", "BBN epoch Um - cosmic time at nucleosynthesis; Friedmann time-density", "t_BBN~180 s at T~0.1 MeV; primordial B-field constraint from BBN", 0.0, "J/m^3", "BBN Um", "F" },
        { 802, "Um_qnm", "Um_qnm = Um_general * a*c^3/(2*G*M); l=2, m=2", "Quasi-normal mode Um - BH ringdown; Kerr spin a; dominant harmonic", "omega_QNM = c^3/(2piGM)*(0.3737+0.088*a); tau_QNM ~ M", 0.0, "J/m^3", "QNM Um", "F" },
        { 804, "Um_bz", "Um_bz = Um_general * B^2*Omega_H^2*R_H^4; ~ a^2", "Blandford-Znajek Um - jet power magnetic factor; horizon angular velocity", "Omega_H=a*c^3/(2*G*M*r_H); P_BZ~a^2 for small spin", 0.0, "J/m^3", "BZ Um", "F" },
        { 830, "Um_photo", "Um_photo = Um_general * eps_photo; eps_photo = 0.1 to 0.3 * L_X/F_X", "Photoevaporation Um - X-ray efficiency fraction; EUV/FUV heating", "xi=L_X/(n*r^2); eps~0.1-0.3; planetary atmosphere evaporation", 0.0, "J/m^3", "Photoevaporation Um", "F" },
        { 832, "Um_mig", "Um_mig = Um_general * Sigma*(H/r)^(-3)", "Migration Um - disk surface density / aspect ratio; Type I torque scaling", "Sigma=surface density; H/r~0.05; migration rate ~ Sigma*(H/r)^(-3)", 0.0, "J/m^3", "Migration Um", "F" },
        { 854, "Um_termv", "Um_termv = Um_general * (L/L_Edd); ~ 1", "Terminal velocity Um - Eddington luminosity ratio; stellar wind driving", "L/L_Edd~1 for OB stars; v_inf~3*v_esc; mass-loss rate critical", 0.0, "J/m^3", "Terminal velocity Um", "F" },
        { 878, "Um_ps", "Um_ps = Um_general * sigma_M_rms", "Press-Schechter Um - rms density fluctuation on scale M; mass variance", "sigma_M = rms of delta smoothed on sphere of mass M; sigma~1 at M*", 0.0, "J/m^3", "Press-Schechter Um", "F" },
        { 880, "Um_sfe", "Um_sfe = Um_general * mdot_halo; mdot ~ M^(1.1)*(1+z)^(2.25)", "Star formation efficiency Um - halo accretion rate scaling", "Dekel+2013; cold stream mass flux; SFR~mdot_halo*f_b*eps_sf", 0.0, "J/m^3", "SFE Um", "F" },
        { 904, "Um_bhent", "Um_bhent = Um_general; dM = T_H*dS", "BH entropy Um - thermodynamic 1st law; Hawking temperature and entropy", "dS_BH=dA/(4*l_pl^2); T_H*dS=dM*c^2; 1st law of BH mechanics", 0.0, "J/m^3", "BH entropy Um", "F" },
        { 906, "Um_evapl", "Um_evapl = Um_general; dM/dt = -P_evap/c^2; P~hbar*c^6/(G^2*M^2)", "Evaporation lifetime Um - BH Hawking radiation power rate", "P_Hawking~hbar*c^6/(15360*pi*G^2*M^2); tau~M^3; Planck mass remnant?", 0.0, "J/m^3", "BH evaporation Um", "F" },
        { 924, "Um_ang", "Um_ang = Um_general * G*M/r^3 * r_A", "Angular momentum Um - specific angular momentum at Alfven radius r_A", "r_A=Alfven radius in stellar wind; angular momentum shed per unit mass", 0.0, "J/m^3", "Angular momentum Um", "F" },
        { 926, "Um_jetvel", "Um_jetvel = Um_general * Omega*r_A*(r_A/r0)^(1/2)", "Jet velocity Um - magneto-centrifugal launch; Alfven lever arm", "Blandford-Payne 1982; r_A/r0=lever arm ratio; efficiency eta_MHD", 0.0, "J/m^3", "Jet velocity Um", "F" },
        { 944, "Um_chirp", "Um_chirp = Um_general * (32/5)*(G*M_c^(5/3)/c^5)*(pi*f)^(10/3)", "Chirp mass GW Um - inspiral power spectral density", "M_c=chirp mass; GW luminosity P_GW~M_c^(10/3)*f^(10/3)", 0.0, "J/m^3", "Chirp mass Um", "F" },
        { 992, "Um_pec", "Um_pec = Um_general * (-(f*H/(4*pi))*grad^-1*(delta))", "Peculiar velocity Um - linear density-velocity relation", "f=d ln D/d ln a~Omega_m^0.55; redshift space distortion", 0.0, "J/m^3", "Peculiar velocity Um", "F" },
        { 994, "Um_ion", "Um_ion = Um_general * Int(alpha_B*n_e^2*C*dt)", "Ionisation Um - recombination rate integral; clumping factor C", "alpha_B=case B recombination; C=<n^2>/<n>^2 ~3; reionization budget", 0.0, "J/m^3", "Reionisation Um", "F" },
        { 1012, "Um_deb", "Um_deb = Um_general * rho_rad; rho_rad=pi^2*k_B^4*T^4/(30*hbar^3*c^5)*g_star", "Debye/radiation epoch Um - radiation energy density at temperature T", "g_star=106.75 SM; scales as T^4; dominates for T>Teq~3000 K", 0.0, "J/m^3", "Radiation epoch Um", "F" },
        { 1014, "Um_fried1", "Um_fried1 = Um_general * sqrt(Omega_m*(1+z)^3+Omega_k*(1+z)^2+Omega_Lambda)", "Friedmann Um - H(z) normalised Hubble parameter; all energy components", "H(z)=H0*E(z); E^2=Omega_m*(1+z)^3+Omega_k*(1+z)^2+Omega_Lambda", 0.0, "J/m^3", "Friedmann Um", "F" },
        { 1034, "Um_damp", "Um_damp = Um_general * Q_QNM; Q~10 for l=2, m=2", "QNM damping Um - quality factor Q of ringdown oscillation", "Q=pi*f*tau; l=2,m=2 most excited; Q~2-12 for a=0-0.99", 0.0, "J/m^3", "QNM damping Um", "F" },
        { 1036, "Um_arnett", "Um_arnett = Um_general * t_d^2; t_d^2 = 3*kappa*M_ej/(4*pi*c*v_ej^2)", "Arnett diffusion time Um - SN Ia light curve peak timescale", "Arnett 1982; M_max ~ mdot_Ni; t_rise~20 days for typical SN Ia", 0.0, "J/m^3", "Arnett SN diffusion Um", "F" },
        { 1058, "Um_dyn", "Um_dyn = Um_general * eta_resistivity", "Dynamo resistivity Um - Ohmic diffusivity eta; reconnection coupling", "eta=c^2/(4*pi*sigma); resistive MHD; fast reconnection for large Rm", 0.0, "J/m^3", "Dynamo resistivity Um", "F" },
        { 1080, "Um_metal", "Um_metal = Um_general; d(Z*M_gas)/dt = Y*SFR - Z*mdot_out", "Metallicity Um - metal production minus outflow; yield Y; loading", "Effective yield Y_eff=Y/(1+eta_load); mass-metallicity scatter", 0.0, "J/m^3", "Metallicity Um", "F" },
        { 1082, "Um_cool", "Um_cool = Um_general * n_e*n_i*Lambda_cool(T)", "Cooling Um - bremsstrahlung + line cooling; Lambda(T) cooling function", "Lambda~10^-23 erg cm^3/s at T~10^7 K; tcool~3kT/(n*Lambda)", 0.0, "J/m^3", "Cooling function Um", "F" },
        { 1100, "Um_rev", "Um_rev = Um_general * eta_resistivity", "Reversal Um - same resistivity coupling as Um_dyn but in reversal context", "Magnetic polarity reversal; alpha-Omega dynamo; l_rev scale", 0.0, "J/m^3", "Magnetic reversal Um", "F" },
        { 1102, "Um_wde", "Um_wde = Um_general * rho_DE; rho_DE ~ a^(-3*(1+w))", "w-dark energy Um - quintessence density evolution with w", "w=-1: constant; w>-1: decreasing; tracker solutions phi^2/2~V", 0.0, "J/m^3", "w-dark energy Um", "F" },
        { 1252, "Um_vir", "Um_vir = Um_general; sigma_v^2 = G*M/(3*r) for gas; 2*KE+PE=0", "Virial Um - equilibrium cluster velocity dispersion", "T_vir=m*G*M/(3*k_B*r); X-ray temperature diagnostic", 0.0, "J/m^3", "Virial Um", "F" },
        { 1272, "Um_ms", "Um_ms = Um_general * L_MS; L_MS ~ mu^4*M^3 [Eddington main seq.]", "Main sequence luminosity Um - mass-luminosity scaling", "mu=mean molecular weight; L~M^4 for stars; opacity-dominated", 0.0, "J/m^3", "Main sequence Um", "F" },
        { 1274, "Um_ml", "Um_ml = Um_general; L ~ M^3.5 from HR diagram empirical fit", "Mass-luminosity Um - observational fit from HR diagram", "Empirical; L/Lsun~(M/Msun)^3.5; valid for 0.1-50 Msun", 0.0, "J/m^3", "Mass-luminosity Um", "F" },
        { 753, "Um_star", "Um_star = Um_general * (sigma_v^3/(G^2*rho*m*lnLambda))", "Star cluster Um - relaxation timescale coupling; stellar encounters", "t_re=0.34*sigma_v^3/(G^2*rho*m*lnLambda); globular cluster dynamics", 0.0, "J/m^3", "Star cluster Um", "F" },
        { 881, "Um_fb", "Um_fb = Um_general * eta_fb*mdot_star*c^2", "Feedback Um - efficiency-weighted energy injection rate", "eta_fb~0.001 stellar winds; ~0.1 AGN; quenches SF in massive halos", 0.0, "J/m^3", "Feedback Um", "F" },
        { 899, "Um_grow", "Um_grow = Um_general * delta; delta_ddot+2*H*delta_dot=(3/2)*Omega_m*H^2*delta", "Growth Um - density perturbation amplitude coupling", "D(a)~a matter dom; growth rate f=d ln D/d ln a~Omega_m^0.55", 0.0, "J/m^3", "Growth Um", "F" },
        { 901, "Um_haw", "Um_haw = Um_general * T_H; T_H=hbar*c^3/(8*pi*G*M*k_B)", "Hawking temperature Um - thermal radiation coupling from BH horizon", "T_H~6e-8 K (Msun); BH magnetic moment analogue in thermal field", 0.0, "J/m^3", "Hawking Um", "F" },
        { 907, "Um_lqcf", "Um_lqcf = Um_general * H_LQC; H^2=(8*pi*G*rho/3)*(1-rho/rho_crit)", "LQC Friedmann Um - loop quantum bounce Hubble rate coupling", "Polymerisation correction; maximum rho=rho_crit=0.41*rho_Pl", 0.0, "J/m^3", "LQC Friedmann Um", "F" },
        { 911, "Um_bounc", "Um_bounc = Um_general * t_bounce; t_bounce=sqrt(3/(16*pi*G*rho_crit))", "Bounce time Um - LQC quantum bounce timescale", "t_bounce~Planck time; connects pre/post-bounce epochs", 0.0, "J/m^3", "LQC bounce Um", "F" },
        { 913, "Um_ent", "Um_ent = Um_general * S_ent*exp(-t/tau_dec); S=-Tr(rho_A*ln rho_A)", "Entanglement entropy Um - decoherence-weighted von Neumann entropy", "tau_dec=decoherence time; AdS/CFT Ryu-Takayanagi formula", 0.0, "J/m^3", "Entanglement Um", "F" },
        { 919, "Um_jshock", "Um_jshock = Um_general * (rho1*v1^2+P1)", "J-shock Um - total momentum flux across discontinuous shock", "Rankine-Hugoniot; adiabatic index gamma=5/3 monatomic gas", 0.0, "J/m^3", "J-shock Um", "F" },
        { 921, "Um_cshock", "Um_cshock = Um_general * rho_n*v_n*nu_ni", "C-shock Um - ion-neutral friction rate; continuous shock structure", "nu_ni~1.5e-9 cm^3/s * n_i; L_d=v_di/nu_ni", 0.0, "J/m^3", "C-shock Um", "F" },
        { 923, "Um_halo", "Um_halo = Um_general * (dN_halo/dM_dz)", "Halo formation Um - merger tree rate coupling", "Lacey-Cole 1993; conditional probability; EPS framework", 0.0, "J/m^3", "Halo Um", "F" },
        { 925, "Um_disk", "Um_disk = Um_general * eps*M_gas/t_dyn", "Disk SF Um - gravitational instability star formation rate", "Q<1 for instability; Kennicutt-Schmidt; Sigma_SFR~Sigma_gas^1.4", 0.0, "J/m^3", "Disk SF Um", "F" },
        { 927, "Um_burst", "Um_burst = Um_general * mdot_gas_inflow*eps_burst", "Starburst Um - merger-driven gas inflow rate", "eps_burst~0.5-1; triggered by tidal compression; Arp 220 analogue", 0.0, "J/m^3", "Starburst Um", "F" },
        { 929, "Um_sedov", "Um_sedov = Um_general * (E*t^2/rho)^(1/5)", "Sedov-Taylor Um - blast wave radius coupling to magnetism", "R_sedov~t^(2/5); transitions to radiative at R~30 pc", 0.0, "J/m^3", "Sedov-Taylor Um", "F" },
        { 931, "Um_dsa", "Um_dsa = Um_general * u_s^2/c^2", "DSA Um - shock velocity squared coupling; particle acceleration", "N(E)~E^(-2) for strong shock; maximum E_max=Z*e*B*R_Hillas", 0.0, "J/m^3", "DSA Um", "F" },
        { 941, "Um_tov", "Um_tov = Um_general * (dP/dr_TOV)", "TOV Um - neutron star pressure gradient coupling", "Stiff EOS: R~13 km; soft EOS: R~10 km; GW170817 R>10.7 km", 0.0, "J/m^3", "TOV Um", "F" },
        { 995, "Um_acc", "Um_acc = Um_general * (L_Edd/(eps*c^2))*exp(t/t_Sal)", "Accretion Um - Eddington rate exponential growth; Salpeter time", "t_Sal=45 Myr; e-folding to SMBH from seed; quasar duty cycle", 0.0, "J/m^3", "Eddington accretion Um", "F" },
        { 981, "Um_IGM", "Um_IGM = Um_general * k_B*T_IGM/m_p", "IGM thermal Um - post-reionisation warm-hot gas coupling", "T_IGM~1e4-1e7 K; WHIM 30-40% of baryons; Ly-alpha forest", 0.0, "J/m^3", "IGM Um", "F" },
        { 983, "Um_gal", "Um_gal = Um_general * M_halo*H(z)^2", "Galaxy formation Um - halo mass times Hubble growth factor", "SFR~M_halo*f_b*eps_sf*H(z); cosmic SFR peaks z~2", 0.0, "J/m^3", "Galaxy formation Um", "F" },
        { 985, "Um_quant", "Um_quant = Um_general * H^2/(8*pi^2*eps*M_pl^2)", "Quantum inflation Um - inflationary power spectrum amplitude", "A_s~2.2e-9; eps~0.004 Starobinsky; pivot k_0=0.05 Mpc^-1", 0.0, "J/m^3", "Quantum inflation Um", "F" },
        { 989, "Um_w", "Um_w = Um_general * rho_DE*(1+w)", "DE equation-of-state Um - quintessence pressure coupling", "w+1 tracks field kinetic energy; w~-0.95 DESI 2024", 0.0, "J/m^3", "DE EOS Um", "F" },
        { 991, "Um_ng", "Um_ng = Um_general * f_NL; f_NL<10 (Planck 2018)", "Non-Gaussianity Um - primordial bispectrum coupling", "Local f_NL; equilateral f_NL; orthogonal f_NL; multi-field inflation", 0.0, "J/m^3", "Non-Gaussianity Um", "F" },
        { 965, "Um_nfwrot", "Um_nfwrot = Um_general * v_NFW^2; v^2=4*pi*G*rho_s*r_s^2*(ln(1+x)-x/(1+x))/r", "NFW rotation curve Um - halo profile velocity coupling", "rho_s*r_s^3=const; concentration c=r_vir/r_s; c~10-15 Milky Way", 0.0, "J/m^3", "NFW Um", "F" },
        { 967, "Um_sidm", "Um_sidm = Um_general * (sigma_m/m)", "SIDM cross-section Um - self-interaction velocity-dependent coupling", "sigma/m~1 cm^2/g for core formation; v^(-4) for light mediator", 0.0, "J/m^3", "SIDM Um", "F" },
        { 772, "Um_void", "Um_void = Um_general * (-3/5)*(Omega_m*a+Omega_Lambda)^(-3/2)*delta_v0", "Void density Um - linear underdense expansion", "Voids evacuate at rate ~ a^(-1); supervoids ~100 Mpc empty", 0.0, "J/m^3", "Void Um", "F" },
        { 774, "Um_reion", "Um_reion = Um_general * ndot_gamma*eps_esc", "Reionization Um - ionising photon emission rate coupling", "eps_esc~0.1-0.2; z_reion~6-10; completed by z~6 (JWST 2024)", 0.0, "J/m^3", "Reionization Um", "F" },
        { 776, "Um_ISM", "Um_ISM = Um_general * pi*c_s^2/(G*rho)", "ISM Jeans Um - thermal support scale coupling to magnetism", "lambda_J=sqrt(pi*c_s^2/(G*rho)); B-field adds magnetic Jeans scale", 0.0, "J/m^3", "ISM Um", "F" },
        { 782, "Um_merg2", "Um_merg2 = Um_general * (5/(256*pi^(8/3)))*c^5*M_c^(-5/3)*f^(-11/3)", "Merger coalescence Um - GW inspiral frequency-mass coupling", "Time to coalescence from f_i; dE/dt=-P_GW; orbit shrinks", 0.0, "J/m^3", "GW coalescence Um", "F" },
        { 739, "Um_CR", "Um_CR = Um_general * E_max * Z", "Cosmic ray Um - maximum energy times charge; Hillas criterion", "E_max=Z*e*B*R; Auger UHECRs E>1e19 eV; ankle E~3e18 eV", 0.0, "J/m^3", "Cosmic ray Um", "F" },
        { 821, "Um_agn", "Um_agn = Um_general * L_AGN/c", "AGN radiation Um - momentum flux from quasar luminosity", "L_AGN up to 10^48 erg/s; momentum coupling to host ISM", 0.0, "J/m^3", "AGN Um", "F" },
        { 1965, "Um_jeans", "Um_jeans = Um_general * sqrt(pi*c_s^2/(G*rho))", "Jeans scale Um - instability wavelength magnetic coupling", "M_J=(4*pi/3)*(lambda_J/2)^3*rho; cloud collapse criterion", 0.0, "J/m^3", "Jeans Um", "F" },
        { 1985, "Um_relax", "Um_relax = Um_general * t_relax; t_relax=0.34*sigma_v^3/(G^2*rho*m*lnLambda)", "Relaxation Um - two-body relaxation timescale coupling", "Phase mixing; ergodic distribution; violent relaxation Lynden-Bell", 0.0, "J/m^3", "Relaxation Um", "F" },
        { 1945, "Um_voidden", "Um_voidden = Um_general * delta_v(a); delta_v~a^(-1)", "Void density contrast Um - underdense region evolution", "Linear density contrast in void; compensated void profile", 0.0, "J/m^3", "Void density Um", "F" },
        { 2005, "Um_conv", "Um_conv = Um_general * v_conv; v_conv~(alpha^2*g*deltaT*H_p/(4*T))^(1/3)", "Convective Um - mixing-length velocity magnetic coupling", "alpha=MLT parameter; H_p=pressure scale height; stellar convection zones", 0.0, "J/m^3", "Convective Um", "F" },
        { 1805, "Um_eta", "Um_eta = Um_general * eta_b; eta_b=n_b/n_gamma=6e-10", "Baryon-photon ratio Um - BBN era baryon asymmetry coupling", "Measured via D/H primordial; consistent with CMB Omega_b", 6e-10, "J/m^3", "Baryon-photon Um", "F" },
        { 1846, "Um_GW", "Um_GW = Um_general * rho_GW; rho_GW=(c^2/(32*pi*G))*<hdot_ij^2>", "Gravitational wave Um - GW energy density coupling to magnetism", "PTA nHz background; LIGO HF; UQFF connects GW to Um", 0.0, "J/m^3", "Gravitational wave Um", "F" },
        { 1847, "Um_inf", "Um_inf = Um_general * V(phi); V~(3*H^2*M_pl^2/8*pi)*exp(N)", "Inflation potential Um - slow-roll potential energy coupling", "Chaotic/Starobinsky/natural inflation; r<0.06 Planck", 0.0, "J/m^3", "Inflation potential Um", "F" },
        { 1848, "Um_glitch", "Um_glitch = Um_general * Delta_Omega/Omega; Delta_Omega~1e-6 to 1e-4", "Pulsar glitch Um - fractional spin-up coupling to magnetism", "Vortex unpinning; crustal superfluid; Vela glitch Delta_nu/nu~1e-6", 0.0, "J/m^3", "Pulsar glitch Um", "F" },
        { 1849, "Um_after", "Um_after = Um_general * F_nu_afterglow; F_nu ~ nu^(-(p-1)/2)*t^(-alpha_t)", "GRB afterglow Um - synchrotron flux spectrum coupling", "p~2.2; alpha_t=3(p-1)/4; jet break steepens decay", 0.0, "J/m^3", "GRB afterglow Um", "F" },
        { 1850, "Um_cmb_pol", "Um_cmb_pol = Um_general * C_l_BB; CMB B-mode polarisation power", "CMB B-mode Um - primordial GW through lensing B-mode coupling", "r<0.06 Planck; BICEP/Keck; dusty CMB foreground subtraction challenge", 0.0, "J/m^3", "CMB B-mode Um", "F" },
        { 1851, "Um_duty", "Um_duty = Um_general * (1-exp(-t/tau_cool))*(1+mdot_acc/mdot_Edd)^(-1)", "AGN duty cycle Um - time-averaged accretion activity coupling", "Duty_cycle~0.01-0.1; merger-triggered AGN episodes; calorimetry", 0.0, "J/m^3", "AGN duty Um", "F" },
    };
}

// ======================================================
// SECTION G:  Numerical Solutions and Calibration Constants
// Source: grok_share_c020496d9e.txt; UQFF Calibration 22Sept2025 PDF
// ======================================================
std::vector<BuoyancyEquation> buildSectionG() {
    return {
        // G-1  F_U_Bi_i total (LENR dominant)
        {  701, "F_UBii_total_LENR",
           "F_U_Bi_i (LENR dominant) ~= +2.11e208 N; x2 ~= -1.35e172 m",
           "Master integral evaluated with LENR term dominating; extremely large buoyancy",
           "Key result: LENR coupling drives macroscopic buoyancy across hierarchy",
           2.11e208, "N", "Universal buoyancy integral", "G" },

        // G-2  F_U_Bi_i total (F_rel dominant — negative)
        {  702, "F_UBii_total_Frel",
           "F_U_Bi_i (F_rel dominant) ~= -8.31e211 N",
           "Master integral with F_rel dominant; negative indicates inward suspension",
           "Negative buoyancy: inside-looking suppression dominates at relativistic scales",
           -8.31e211, "N", "Universal buoyancy integral (F_rel)", "G" },

        // G-3  F_rel — Relativistic LEP anchor
        {  703, "F_rel_LEP",
           "F_rel = 4.31e33 N  [2024 LEP re-analysis cross-validation]",
           "Primary validation anchor; 2024 LEP collision energy ratio squared",
           "All 79 F_UBii types normalised to F_rel; LEP CM energy~208 GeV",
           4.31e33, "N", "2024 LEP validation", "G" },

        // G-4  F_LENR — Low Energy Nuclear Reaction
        {  704, "F_LENR_value",
           "F_LENR = k_LENR * (omega_LENR/omega0)^2 ~= 1.56e36 N",
           "LENR resonance force; dominant term in F_U_Bi_i master integral",
           "Widom-Larsen LENR coupling; k_LENR calibrated to cold fusion data",
           1.56e36, "N", "LENR resonance force", "G" },

        // G-5  FU_g1 Westerlund 2
        {  705, "FUg1_Westerlund2",
           "FU_g1 (Westerlund 2) ~= 2.43e-40 N",
           "Triadic compressed gravity for Westerlund 2 star cluster",
           "Collapse drive; numerically validated against JWST NIRCam photometry",
           2.43e-40, "N", "Westerlund 2", "G" },

        // G-6  R(t) Westerlund 2
        {  706, "Rt_Westerlund2",
           "R(t) (Westerlund 2) ~= -2.29e-41 N",
           "26-layer oscillatory erosion term for Westerlund 2",
           "Negative: oscillatory erosion of FU_g1; net buoyancy = 2.43e-40 - 2.29e-41",
           -2.29e-41, "N", "Westerlund 2 resonance", "G" },

        // G-7  FU_Bi Westerlund 2
        {  707, "FUBi_Westerlund2",
           "FU_Bi (Westerlund 2) ~= 6.14e-32 N  [f_Ub * 2.20e8]",
           "Universal Buoyancy (inside->outside) for Westerlund 2",
           "Dominant over FU_g1 by 8 orders of magnitude; buoyancy > gravity",
           6.14e-32, "N", "Westerlund 2 buoyancy", "G" },

        // G-8  FU_g1 Pillars of Creation
        {  708, "FUg1_Pillars",
           "FU_g1 (Pillars of Creation) ~= 3.95e-41 N",
           "Triadic gravity for Pillars of Creation M16",
           "Photodissociation region; wind-driven collapse suppressed",
           3.95e-41, "N", "Pillars of Creation", "G" },

        // G-9  R(t) Pillars of Creation
        {  709, "Rt_Pillars",
           "R(t) (Pillars of Creation) ~= -1.12e-42 N",
           "Resonance erosion for Pillars; 26-layer oscillation",
           "Smaller than Westerlund 2; diffuse nebular environment",
           -1.12e-42, "N", "Pillars resonance", "G" },

        // G-10  FU_Bi Pillars of Creation
        {  710, "FUBi_Pillars",
           "FU_Bi (Pillars of Creation) ~= 9.79e-33 N  [f_Ub * 2.20e7]",
           "Universal Buoyancy for Pillars; factor f_Ub*2.20e7",
           "Buoyancy > gravity by 8 orders; consistent with Westerlund 2 ratio",
           9.79e-33, "N", "Pillars buoyancy", "G" },

        // G-11  U_m energy density
        {  711, "Um_density_value",
           "U_m ~= 3.78e-6 J/m^3",
           "Universal Magnetism energy density; evaluated at calibration point",
           "Calibrated against 48-system ensemble; Q_wave_mean=3.97e4 J/m^3",
           3.78e-6, "J/m^3", "Universal magnetism", "G" },

        // G-12  E_neutrino
        {  712, "E_neutrino_value",
           "E_neutrino ~= 1.05e5 eV",
           "Neutrino energy from layered vacuum density ratio; UQFF prediction",
           "Consistent with KATRIN upper bound; vacuum-mediated mass generation",
           1.05e5, "eV", "Neutrino energy", "G" },

        // G-13  Decay Rate
        {  713, "decay_rate_value",
           "Decay_Rate ~= 0.0583",
           "Vacuum-mediated particle decay rate; ratio rho_SCm/rho_UA * exp(-[SSq])",
           "Dimensionless; calibrated at n=1; scales with layer index",
           0.0583, "s^-1", "Vacuum decay rate", "G" },

        // G-14  [SSq] calibrations
        {  714, "SSq_calibration_values",
           "[SSq]=0.5 (low-n layers); [SSq]=5.26 (n=26, cosmic scale)",
           "UQFF squeezing factor calibrated against 96 astrophysical systems",
           "Q_26([SSq]=5.26) ~= 6.63e21 (26th Ramanujan polynomial at cosmic scale)",
           0.5, "dimensionless", "SSq calibration", "G" },

        // G-15  gamma (decay constant)
        {  715, "gamma_temporal",
           "gamma ~= 5e-5 day^-1",
           "Temporal decay constant in Um general form; slow cosmic evolution",
           "gamma^-1 ~= 20000 days ~= 55 years; stellar cluster timescale",
           5.0e-5, "day^-1", "Um temporal decay", "G" },

        // G-16  k_UV and k_mm coupling constants
        {  716, "k_UV_k_mm",
           "k_UV = k_mm = 1e-30 N/W;  f_mm = 1.05",
           "UV (GALEX/Spitzer) and mm (ALMA) luminosity coupling constants",
           "f_mm=1.05 protoplanetary disk correction; k validated to ALMA data",
           1.0e-30, "N/W", "UV+mm coupling", "G" },

        // G-17  DPM_resonance
        {  717, "DPM_resonance_value",
           "DPM_resonance ~= 1.67e7",
           "Dynamic Phase Modulation resonance amplitude; dimensionless",
           "Calibrated to Sgr A* THz emission cycle f_TRZ=5.95e-4 Hz",
           1.67e7, "dimensionless", "DPM resonance", "G" },

        // G-18  Q_wave ensemble statistics
        {  718, "Q_wave_statistics",
           "Q_wave mean = 3.97e4 J/m^3;  Q_wave std = 6.35e4 J/m^3  [48 systems]",
           "Wave energy density statistics across 48-system calibration ensemble",
           "Log-normal distribution; std > mean indicates wide dynamic range",
           3.97e4, "J/m^3", "Wave energy density ensemble", "G" },

        // G-19  phi calibration
        {  719, "phi_calibration",
           "phi = sin(pi*t_n) + 0.01*cos(2*pi*f_flare*t) ~= 1.02",
           "Calibrated phi for pseudo-monopole phase; includes flare modulation",
           "f_flare=stellar flare frequency; 0.01 coefficient small perturbation",
           1.02, "dimensionless", "Phi calibration", "G" },

        // G-20  f_TRZ — Sgr A* 28-minute cycle
        {  720, "f_TRZ_SgrA",
           "f_TRZ = 5.95e-4 Hz  [Sgr A* 28-minute near-IR flare cycles]",
           "THz resonance zero-point frequency; Sgr A* quasi-periodic oscillation",
           "28.0 min period; Chandra X-ray confirmed; VLBI 2025 corroborated",
           5.95e-4, "Hz", "Sgr A* QPO frequency", "G" },

        // G-21  H2O-H2 cross-section
        {  721, "CS_H2O_H2",
           "CS sigma(Dj=2, E=400 cm^-1) = 11.65 Angstrom^2  [H2O-H2 refit]",
           "Rotationally inelastic cross-section; H2O vibrational pumping in ISM",
           "Refit to JWST-era interstellar water line ratios; UQFF coupling anchor",
           11.65, "Angstrom^2", "H2O-H2 collision cross-section", "G" },

        // G-22  D_universe prediction
        {  722, "D_universe_prediction",
           "D_universe ~= 92.77 Gly  [z=1100, H0=67.4; matches ~93 Gly observed]",
           "UQFF prediction of observable universe diameter; Planck + buoyancy terms",
           "Lambda-CDM predicts 93.1 Gly; UQFF prediction 92.77 Gly: 0.4% accuracy",
           92.77, "Gly", "Observable universe diameter", "G" },

        // G-23  Q_26 Ramanujan 26th polynomial
        {  723, "Q26_Ramanujan",
           "Q_26 ([SSq]=5.26) ~= 6.63e21",
           "26th Ramanujan polynomial evaluated at cosmic [SSq]=5.26",
           "26D polynomial from SOURCE115; encodes 19-system master gravity",
           6.63e21, "dimensionless", "Ramanujan 26th polynomial", "G" },

        // G-24  UQFF solvability
        {  724, "UQFF_solvability",
           "UQFF solvability = 99.9% (Grok 4 analysis, Sept 14-21 2025)",
           "Completeness measure of UQFF equation set vs observational data",
           "96% match to JWST/Chandra 2025 data; 6688+ physics terms registered",
           99.9, "percent", "UQFF solvability", "G" },

        // G-25  H_SCm and U_UA calibrated values
        {  725, "H_SCm_U_UA_kappa",
           "H_SCm~0.99; U_UA~0.0001; k_eta=1e-113; beta_i~0.603; kappa=0.0005/day",
           "Core UQFF calibration constants; kappa=daily decay; beta_i=coupling",
           "Calibrated Sept 2025; [SSq]=0.57 central value; consistent across 29 systems",
           0.0005, "day^-1", "UQFF core calibration", "G" },
    };
}

// ======================================================
// SECTION H:  Lambda-CDM / MOND Comparison and Validation Notes
// Source: grok_share_c020496d9e.txt (~lines 600-800)
// ======================================================
std::vector<BuoyancyEquation> buildSectionH() {
    return {
        // H-1  Lambda-CDM backbone shared with UQFF
        {  801, "LambdaCDM_H_z",
           "H(z) = H0*sqrt(Omega_m*(1+z)^3 + Omega_Lambda);"
           " H0=67.4 km/s/Mpc; Omega_m=0.3; Omega_Lambda=0.7",
           "Lambda-CDM Hubble evolution - shared backbone embedded in all UQFF g(r,t)",
           "Planck 2018; UQFF extends via: vacuum buoyancy+[SSq]+entanglement terms",
           67.4, "km/s/Mpc", "Lambda-CDM cosmology", "H" },

        // H-2  MOND deep-field limit
        {  802, "MOND_deep_acceleration",
           "a_MOND = sqrt(g_N * a0);  a0 ~= 1.2e-10 m/s^2",
           "MOND deep-MOND limit - geometric mean of Newtonian and a0 threshold",
           "UQFF extends MOND via vacuum density ratios rho_vac[UA]/rho_vac[SCm]",
           1.2e-10, "m/s^2", "MOND deep-field", "H" },

        // H-3  UQFF advantage over Lambda-CDM
        {  803, "UQFF_vs_LambdaCDM",
           "g_UQFF = g_LambdaCDM + FU_g1 + R(t) + FU_Bi + [SSq]_terms + ent_terms",
           "UQFF unifies Lambda-CDM and MOND via buoyancy and entanglement extensions",
           "96% match to JWST/Chandra 2025; explains rotation curves without new particles",
           0.0, "m/s^2", "UQFF unification", "H" },

        // H-4  UQFF vs MOND at all scales
        {  804, "UQFF_vs_MOND",
           "UQFF -> MOND as rho_vac[SCm]/rho_vac[UA] -> 1 (low density limit);"
           " UQFF -> Newtonian as [SSq] -> 0 (dense limit)",
           "UQFF reduces to MOND in voids and to GR in dense regions",
           "Smooth interpolation; no free parameter beyond calibrated [SSq]",
           0.0, "m/s^2", "UQFF-MOND limit", "H" },

        // H-5  UQFF vs Observation summary
        {  805, "UQFF_96pct_match",
           "96% match to: JWST NIRCam (2023-2025); Chandra X-ray; ALMA mm;"
           " LIGO GWTC-4.0; Planck 2018; DESI DR1; Gaia DR4",
           "Cross-dataset validation at 96% solvability across all catalogued systems",
           "Remaining 4%: unresolved dark matter substructure and precision baryon models",
           96.0, "percent", "UQFF observational match", "H" },

        // H-6  Buoyancy-gravity ratio summary
        {  806, "buoyancy_over_gravity_ratio",
           "FU_Bi / FU_g1 ~= 6.14e-32 / 2.43e-40 ~= 2.53e8  [Westerlund 2]",
           "Buoyancy exceeds triadic gravity by 8 orders of magnitude at all tested systems",
           "Buoyancy dominance implies perpetual suspension without dark energy tuning",
           2.53e8, "dimensionless", "Buoyancy-gravity ratio", "H" },

        // H-7  F_U_Bi_i vs standard model forces
        {  807, "FUBii_vs_SM_forces",
           "F_U_Bi_i(LENR)~2.11e208 N >> F_strong~1e5 N >> F_EM~1e2 N >> F_grav~1e-38 N",
           "F_U_Bi_i at LENR-dominant integration dwarfs all SM forces",
           "Hierarchy: F_UBii >> F_LENR > F_rel > F_neutron >> F_res; LENR dominant",
           2.11e208, "N", "Force hierarchy", "H" },
    };
}

// ======================================================
// main():  Print complete F_U_Bi_i equation catalogue
// ======================================================
int main()
{
    std::cout << "\n";
    std::cout << "=====================================================\n";
    std::cout << "  F_U_Bi_i_QCalc.cpp  -  Universal Buoyancy Catalogue\n";
    std::cout << "  Author: Daniel T. Murphy   Version: v4.83\n";
    std::cout << "  UQFF (Unified Quantum Field Framework) v4.83\n";
    std::cout << "=====================================================\n";

    auto secA = buildSectionA();
    auto secB = buildSectionB();
    auto secC = buildSectionC();
    auto secD = buildSectionD();
    auto secE = buildSectionE();
    auto secF = buildSectionF();
    auto secG = buildSectionG();
    auto secH = buildSectionH();

    printSection("A — 29 Per-System g_UQFF Equations (Astrophysical Systems)", secA);
    printSection("B — Compressed UQFF Backbone + Triadic Master Equations",    secB);
    printSection("C — Sub-Equations (Um, [SSq], t_n, f_Ub, Ug2, Vacuum Series)", secC);
    printSection("D — F_U_Bi_i Master Integral Component Force Equations",     secD);
    printSection("E — 79 Unique F_UBii Variant Types (BB_C_Equations)",        secE);
    printSection("F — 68 Unique Um (Universal Magnetism) Variant Types",       secF);
    printSection("G — Numerical Solutions and Calibration Constants",          secG);
    printSection("H — Lambda-CDM / MOND Comparison and Validation Notes",      secH);

    // Totals summary
    std::size_t total = secA.size()+secB.size()+secC.size()+secD.size()
                       +secE.size()+secF.size()+secG.size()+secH.size();
    std::cout << "\n\n=====================================================\n";
    std::cout << "  CATALOGUE TOTALS:\n";
    std::cout << "    Section A (g_UQFF systems):   " << secA.size() << " equations\n";
    std::cout << "    Section B (Triadic masters):  " << secB.size() << " equations\n";
    std::cout << "    Section C (Sub-equations):    " << secC.size() << " equations\n";
    std::cout << "    Section D (Force components): " << secD.size() << " equations\n";
    std::cout << "    Section E (F_UBii variants):  " << secE.size() << " equations\n";
    std::cout << "    Section F (Um variants):      " << secF.size() << " equations\n";
    std::cout << "    Section G (Constants):        " << secG.size() << " equations\n";
    std::cout << "    Section H (Validation notes): " << secH.size() << " equations\n";
    std::cout << "    ------------------------------------------\n";
    std::cout << "    TOTAL CATALOGUED:             " << total << " equations\n";
    std::cout << "=====================================================\n";
    std::cout << "\n  KEY RESULTS:\n";
    std::cout << "    F_U_Bi_i (LENR dominant)  ~= +2.11e208 N\n";
    std::cout << "    F_U_Bi_i (F_rel dominant) ~= -8.31e211 N\n";
    std::cout << "    F_rel  = 4.31e33 N  (2024 LEP)\n";
    std::cout << "    F_LENR = 1.56e36 N\n";
    std::cout << "    [SSq] calibrated: 0.5 (low-n) to 5.26 (n=26 cosmic)\n";
    std::cout << "    UQFF solvability: 99.9% (Grok 4, Sept 2025)\n";
    std::cout << "    96% match: JWST/Chandra/ALMA/LIGO/Planck/DESI/Gaia\n";
    std::cout << "\n";

    return 0;
}
