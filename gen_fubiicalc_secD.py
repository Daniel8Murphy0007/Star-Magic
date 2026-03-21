"""
gen_fubiicalc_secD.py
Returns C++ buildSectionD() — 12 F_U_Bi_i component force equations
from the master integral. Source: grok_share_c020496d9e.txt lines ~1-80.
"""


def get_section_D() -> str:
    return r"""
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
"""
