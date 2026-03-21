"""
gen_fubiicalc_secC.py
Returns C++ buildSectionC() — 10 Sub-equations:
  Um (Universal Magnetism), delta_n, rho_vac series,
  E_neutrino, Decay Rate, [SSq], t_n, vacuum density series,
  buoyancy harmonics H_m, Ug2 harmonic series.
Source: grok_share_c020496d9e.txt lines ~150-400.
"""


def get_section_C() -> str:
    return r"""
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
"""
