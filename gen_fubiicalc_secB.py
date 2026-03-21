"""
gen_fubiicalc_secB.py
Returns C++ buildSectionB() — Compressed UQFF backbone + 6 Triadic Master equations.
Source: grok_share_c020496d9e.txt lines ~1-200; UQFF Calibration PDF.
"""


def get_section_B() -> str:
    return r"""
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
"""
