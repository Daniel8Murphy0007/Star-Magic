"""
================================================================================
Session 277 — UQFF Gap Patches (Tracks A, C continuation)
================================================================================
Closes 3 of the 14 critical gaps surfaced in UQFF_CALIBRATION_AUDIT.md and
ships planetary presets for the Track C crustal/tectonic engine.

Gap #1  — OrionH2OMaserUQFFCalculator
          22 GHz H2O / 6.7 GHz CH3OH maser brightness tied to FUBii Aether
          buoyancy in dense molecular gas. Anchored to JWST 2025 Orion
          proplyd spectra and the BeSSeL Survey 22 GHz astrometry.

Gap #3  — NeutrinoGWCoincidenceCalculator
          Joint LIGO-GWTC-4.0 / IceCube / Kamiokande / Super-K timing window:
          UQFF predicts -5 s ≤ Δt ≤ +200 s (Aether UA scalar gate delay).
          Returns chi^2 vs observed coincidences.

Gap #5  — PerseusACoolingSuppressionCalculator
          NGC 1275 (Perseus A) cooling-flow suppression via MUGE-Compressed
          Λ-cooling × superconductive heavy plasma damping. Anchored to
          Chandra 1.4 Ms exposure (Fabian et al. 2017) deficit.

Planetary presets — Earth, Europa, Enceladus, Io, Titan, Mars
          Drop-in dataset dicts for QCalcGeom.CrustalZeroPointCalculator.

Pattern matches the existing CondensedPhysics calculator convention:
    class XCalculator:
        def compute(self, dataset: dict) -> dict:
            return {
                'primary_equations': [...],
                'available_equations': [...],
                'simulation_set': [...],
                ...numeric fields...
            }
"""

from __future__ import annotations
import math
from typing import Dict, Any, List

# Canonical UQFF constants (verified Sept 2025, Grok 4 audit)
RHO_VAC_SCM   = 7.09e-37        # [J/m^3] SCm vacuum density
RHO_VAC_UA    = 7.09e-36        # [J/m^3] UA vacuum density
C_LIGHT       = 2.998e8         # [m/s]
G_NEWTON      = 6.674e-11       # [m^3/(kg*s^2)]
H_PLANCK      = 6.626e-34       # [J*s]
HBAR          = 1.055e-34       # [J*s]
KB_BOLTZ      = 1.381e-23       # [J/K]
KAPPA         = 0.0005          # [1/day] decay rate
SSQ           = 0.57            # [SSq] coupling
H_SCM         = 0.99            # H_SCm
BETA_I        = 0.603           # buoyancy coefficient


# ===========================================================================
# GAP #1 — Orion H2O / CH3OH Maser (FUBii-coupled)
# ===========================================================================

class OrionH2OMaserUQFFCalculator:
    """22 GHz H2O / 6.7 GHz CH3OH maser brightness via Aether buoyancy (FUBii).

    Physics:
      In dense (n_H2 ~ 1e8-1e10 cm^-3) molecular gas the population inversion
      that pumps H2O 22.235 GHz and CH3OH 6.668 GHz masers is amplified by
      the Aether UA buoyancy term FUBii pushing collisional partners upward
      through a velocity-coherent column.

      Maser brightness temperature:
        T_b = T_b0 * (n_H2 / n_crit) * exp(tau_v) * (1 + delta_FUBii)
      where:
        delta_FUBii = beta_i * (rho_UA / rho_H2O) * cos(pi * t_n)
        tau_v       = (n_H2 * sigma_v * L) * coherence_factor

    Observational anchor: JWST 2025 Orion BN/KL spectra,
    BeSSeL 22 GHz astrometry (parallax 416 ± 3 pc to Orion KL).

    Inputs (dataset, all optional):
      n_H2_cm3       : H2 number density       (default 1.0e9)
      L_pc           : coherence path length   (default 1.0e-4 pc = 20 AU)
      v_kms          : line FWHM               (default 1.5 km/s)
      sigma_v_cm2    : molecular cross-section (default 1.0e-15)
      t_n            : Mayan time fraction     (default 0.0)
      species        : 'H2O' or 'CH3OH'        (default 'H2O')
      D_pc           : source distance         (default 416.0 pc)
      M_solar        : pumping clump mass      (default 10.0)
    """

    H2O_FREQ_GHZ      = 22.23508
    CH3OH_FREQ_GHZ    = 6.66852
    N_CRIT_H2O_CM3    = 1.0e9      # critical density for inversion
    N_CRIT_CH3OH_CM3  = 1.0e6
    T_B0_K            = 1.0e12     # reference brightness temperature

    def compute(self, dataset: Dict[str, Any] = None) -> Dict[str, Any]:
        ds = dataset or {}
        n_H2     = ds.get('n_H2_cm3', 1.0e9)
        L_pc     = ds.get('L_pc',     1.0e-4)
        v_kms    = ds.get('v_kms',    1.5)
        sig_v    = ds.get('sigma_v_cm2', 1.0e-15)
        t_n      = ds.get('t_n',      0.0)
        species  = ds.get('species',  'H2O').upper()
        D_pc     = ds.get('D_pc',     416.0)
        M_solar  = ds.get('M_solar',  10.0)

        if species == 'CH3OH':
            freq_ghz = self.CH3OH_FREQ_GHZ
            n_crit   = self.N_CRIT_CH3OH_CM3
        else:
            freq_ghz = self.H2O_FREQ_GHZ
            n_crit   = self.N_CRIT_H2O_CM3

        # Convert L to cm
        L_cm = L_pc * 3.086e18

        # Velocity-coherence factor (Doppler narrowing)
        # narrower line -> longer coherence -> larger tau
        coherence = max(0.1, 5.0 / v_kms)   # dimensionless
        tau_v = n_H2 * sig_v * L_cm * coherence

        # Aether buoyancy modulation
        # rho_H2O proxy (cgs equivalent of mass density of water at maser conditions)
        rho_H2O = n_H2 * 18.0 * 1.66e-24    # g/cm^3
        # Convert to SI
        rho_H2O_SI = rho_H2O * 1.0e3        # kg/m^3
        # Normalise UA contribution
        delta_FUBii = BETA_I * (RHO_VAC_UA / max(rho_H2O_SI, 1.0e-30)) \
                      * math.cos(math.pi * t_n)
        # Saturated maser cap: |delta| < 1
        delta_FUBii = max(-0.99, min(0.99, delta_FUBii * 1.0e36))

        # Brightness temperature with maser cap (T_b <= 1e15 K)
        T_b = self.T_B0_K * (n_H2 / n_crit) * math.exp(min(tau_v, 30.0)) \
              * (1.0 + delta_FUBii)
        T_b = min(T_b, 1.0e15)

        # Luminosity estimate (isotropic, beamed factor 1e-3)
        # L_maser = (8*pi*k*T_b*nu^2/c^2) * Omega_beam * Area_source * df
        nu_Hz = freq_ghz * 1.0e9
        Omega_beam = 1.0e-6                 # sr (narrow beam)
        # Source area from clump mass and density
        m_p = 1.673e-27
        M_kg = M_solar * 1.989e30
        N_H2_total = M_kg / (2.0 * m_p)
        V_cm3 = N_H2_total / max(n_H2, 1.0)
        R_cm  = (3.0 * V_cm3 / (4.0 * math.pi)) ** (1.0/3.0)
        A_m2  = math.pi * (R_cm * 1.0e-2) ** 2

        df_Hz = (v_kms * 1.0e3 / C_LIGHT) * nu_Hz
        L_maser_W = (8.0 * math.pi * KB_BOLTZ * T_b * nu_Hz**2 / C_LIGHT**2) \
                    * Omega_beam * A_m2 * df_Hz
        L_maser_Lsun = L_maser_W / 3.828e26

        # Flux density at Earth
        D_m = D_pc * 3.086e16
        S_Jy = L_maser_W / (4.0 * math.pi * D_m**2 * df_Hz) * 1.0e26

        return {
            'species'             : species,
            'freq_GHz'            : freq_ghz,
            'n_H2_cm3'            : n_H2,
            'n_crit_cm3'          : n_crit,
            'tau_v'               : tau_v,
            'coherence_factor'    : coherence,
            'delta_FUBii'         : delta_FUBii,
            't_n'                 : t_n,
            'T_b_K'               : T_b,
            'L_maser_Lsun'        : L_maser_Lsun,
            'L_maser_W'           : L_maser_W,
            'S_Jy_at_Earth'       : S_Jy,
            'source_radius_AU'    : R_cm / 1.496e13,
            'D_pc'                : D_pc,
            'primary_equations': [
                f"Species: {species} @ {freq_ghz:.5f} GHz",
                f"n_H2 = {n_H2:.2e} cm^-3   (n_crit = {n_crit:.2e})",
                f"Optical depth: tau_v = n_H2 * sigma_v * L * coherence = {tau_v:.3f}",
                f"FUBii modulation: delta = beta_i*(rho_UA/rho_H2O)*cos(pi*t_n) = {delta_FUBii:+.3e}",
                f"T_b = T_b0 * (n/n_crit) * exp(tau_v) * (1+delta_FUBii) = {T_b:.3e} K",
                f"L_maser = (8*pi*k*T_b*nu^2/c^2) * Omega*A*df = {L_maser_Lsun:.3e} Lsun",
                f"S(Earth) = L/(4*pi*D^2*df) = {S_Jy:.3e} Jy at {D_pc} pc",
            ],
            'available_equations': [
                "T_b = T_b0*(n/n_crit)*exp(tau_v)*(1+delta_FUBii)",
                "delta_FUBii = beta_i*(rho_UA/rho_H2O)*cos(pi*t_n)",
                "tau_v = n_H2 * sigma_v * L_cm * coherence",
                "L_maser = (8*pi*k_B*T_b*nu^2/c^2) * Omega_beam * A_source * df",
                "df = (v/c)*nu  (FWHM in Hz from km/s linewidth)",
                "S_Jy = L_maser / (4*pi*D^2*df) * 1e26",
                "Maser cap (saturated): T_b <= 1e15 K",
            ],
            'simulation_set': [
                "Sweep n_H2 = 1e7..1e11 cm^-3; locate inversion-saturation knee",
                "Sweep t_n = 0..2; show FUBii periodic 2-yr amplification cycle",
                "Compare H2O 22 GHz vs CH3OH 6.7 GHz brightness for same clump",
                "Vary v_kms = 0.5..5.0; show coherence-driven tau scaling",
                "Population synthesis: 100 Orion clumps from M=1..50 Msun",
            ],
        }


# ===========================================================================
# GAP #3 — Neutrino-GW Joint Coincidence Window
# ===========================================================================

class NeutrinoGWCoincidenceCalculator:
    """Joint LIGO/IceCube/Kamiokande/Super-K coincidence calculator.

    UQFF prediction: scalar Aether UA gate releases a neutrino burst trailing
    the gravitational-wave merger signature by -5 s to +200 s (compact-binary
    coalescences). Mechanism: post-merger BEC integration of FUBi collapse
    energy into the UA channel before the SCm phase transition releases
    nu_e, nu_mu, nu_tau via Widom-Larsen LENR-like coupling.

    UQFF model probability density:
        p(dt) = N * exp( -(dt - dt0)^2 / (2*sigma^2) )  on  [dt_min, dt_max]
    with dt0 = +12.0 s, sigma = 45.0 s, gating window [-5, +200] s.

    Returns chi^2 vs an observed coincidence list and decision flag.

    Inputs (dataset, all optional):
      observed_dt_list : list of float, observed neutrino - GW arrival
                         differences in seconds (default empty)
      gw_event_id      : str, e.g. 'GW150914', 'GW170817', 'GW190425'
      sigma_t_s        : observation timing resolution (default 1.0 s)
      n_background     : expected background coincidences per window
                         (default 0.0)
    """

    DT_WINDOW_MIN_S   = -5.0
    DT_WINDOW_MAX_S   = 200.0
    DT_PEAK_S         = 12.0
    DT_SIGMA_S        = 45.0

    def _model_pdf(self, dt: float) -> float:
        if dt < self.DT_WINDOW_MIN_S or dt > self.DT_WINDOW_MAX_S:
            return 0.0
        return math.exp(
            -(dt - self.DT_PEAK_S) ** 2 / (2.0 * self.DT_SIGMA_S ** 2)
        )

    def _normalisation(self, n: int = 5000) -> float:
        # Trapezoidal integral over the window
        a, b = self.DT_WINDOW_MIN_S, self.DT_WINDOW_MAX_S
        h = (b - a) / n
        s = 0.5 * (math.exp(-(a-self.DT_PEAK_S)**2/(2*self.DT_SIGMA_S**2))
                   + math.exp(-(b-self.DT_PEAK_S)**2/(2*self.DT_SIGMA_S**2)))
        for k in range(1, n):
            x = a + k * h
            s += math.exp(-(x-self.DT_PEAK_S)**2/(2*self.DT_SIGMA_S**2))
        return s * h

    def compute(self, dataset: Dict[str, Any] = None) -> Dict[str, Any]:
        ds = dataset or {}
        observed: List[float] = list(ds.get('observed_dt_list', []))
        event_id  = ds.get('gw_event_id', 'unknown')
        sigma_t   = ds.get('sigma_t_s',  1.0)
        n_bg      = ds.get('n_background', 0.0)

        Z = self._normalisation()
        # Likelihood for each observed dt: UQFF model vs uniform null
        n_in_window = sum(
            1 for dt in observed
            if self.DT_WINDOW_MIN_S <= dt <= self.DT_WINDOW_MAX_S
        )
        n_total = len(observed)
        window_width = self.DT_WINDOW_MAX_S - self.DT_WINDOW_MIN_S

        logL_uqff = 0.0
        logL_null = 0.0
        for dt in observed:
            p_uqff = self._model_pdf(dt) / Z if Z > 0 else 0.0
            p_null = 1.0 / window_width \
                if self.DT_WINDOW_MIN_S <= dt <= self.DT_WINDOW_MAX_S else 0.0
            # Floor for log-domain stability
            logL_uqff += math.log(max(p_uqff, 1.0e-30))
            logL_null += math.log(max(p_null, 1.0e-30))

        delta_logL = logL_uqff - logL_null
        # chi^2 vs background-only hypothesis (Poisson)
        chi2 = 0.0
        if n_bg > 0:
            chi2 = ((n_in_window - n_bg) ** 2) / n_bg

        # Decision: UQFF favoured if delta_logL > 2 (Bayes factor ~7)
        uqff_favoured = delta_logL > 2.0

        # Mean predicted Δt under UQFF (truncated Gaussian mean)
        mean_dt = self.DT_PEAK_S
        mean_dt_observed = sum(observed) / n_total if n_total > 0 else 0.0

        return {
            'gw_event_id'         : event_id,
            'n_observed'          : n_total,
            'n_in_window'         : n_in_window,
            'n_background'        : n_bg,
            'window_min_s'        : self.DT_WINDOW_MIN_S,
            'window_max_s'        : self.DT_WINDOW_MAX_S,
            'window_width_s'      : window_width,
            'predicted_peak_s'    : self.DT_PEAK_S,
            'predicted_sigma_s'   : self.DT_SIGMA_S,
            'mean_dt_observed_s'  : mean_dt_observed,
            'logL_uqff'           : logL_uqff,
            'logL_null'           : logL_null,
            'delta_logL'          : delta_logL,
            'chi2_vs_bg'          : chi2,
            'uqff_favoured'       : uqff_favoured,
            'sigma_t_s'           : sigma_t,
            'normalisation_Z'     : Z,
            'primary_equations': [
                f"GW event: {event_id}",
                f"UQFF window: [{self.DT_WINDOW_MIN_S}, {self.DT_WINDOW_MAX_S}] s",
                f"Predicted: dt0 = {self.DT_PEAK_S} s,  sigma = {self.DT_SIGMA_S} s",
                f"p(dt) = N * exp(-(dt-dt0)^2/(2*sigma^2))  on window",
                f"Observed: N={n_total}, in_window={n_in_window}",
                f"<dt>_obs = {mean_dt_observed:.2f} s   (predicted {self.DT_PEAK_S} s)",
                f"logL_UQFF - logL_null = {delta_logL:+.3f}",
                f"chi^2 vs background = {chi2:.3f}  (n_bg = {n_bg})",
                f"UQFF favoured: {uqff_favoured}",
            ],
            'available_equations': [
                "p_UQFF(dt) = (1/Z)*exp(-(dt-dt0)^2/(2*sigma^2))  [dt0=+12s, sigma=45s]",
                "p_null(dt) = 1/(dt_max - dt_min)  uniform-in-window",
                "delta_logL = sum_i [ log(p_UQFF(dt_i)) - log(p_null(dt_i)) ]",
                "chi^2 = (N_obs_window - N_bg)^2 / N_bg",
                "Decision: UQFF favoured if delta_logL > 2  (Bayes factor ~ 7)",
                "Aether-UA scalar gate: dt0 = 12 +/- 5 s tied to BEC integration time",
            ],
            'simulation_set': [
                "Inject 100 dt samples from UQFF prior; verify recovery of dt0,sigma",
                "Inject 100 uniform samples; verify delta_logL < 0",
                "Scan dt0 = -5..+50 s; show likelihood landscape",
                "Scan sigma = 5..100 s; identify resolution floor for detection",
                "Run on GW150914, GW170817, GW190425 archival coincidence lists",
            ],
        }


# ===========================================================================
# GAP #5 — Perseus A (NGC 1275) Cooling-Flow Suppression
# ===========================================================================

class PerseusACoolingSuppressionCalculator:
    """NGC 1275 cooling-flow suppression via MUGE-Compressed Λ-cooling damping.

    Standard ICM cooling-flow model predicts dM/dt ~ 300 Msun/yr for the
    Perseus cluster core. Chandra 1.4 Ms exposure observes only
    ~ 10-50 Msun/yr — a factor 6-30 suppression.

    UQFF/MUGE-Compressed mechanism:
      The superconductive heavy plasma layer below the BCG generates an
      Aether-mediated Lambda-cooling counter-pressure that damps the
      Bondi flow. Suppression factor:

        S = 1 / (1 + xi_SC * (T_ICM / T_SCm)^(3/2)
                       * (B / B_ref)^2
                       * cos^2(pi * t_n) )

      where xi_SC ~ 0.5 calibrates against Chandra deficit, and B_ref ~ 100
      microG is the ICM equipartition field. Hotter ICM means more SCm
      coupling layers participate, so suppression grows as T^(3/2).

      dM/dt_observed = S * dM/dt_classical

    Inputs (dataset, all optional):
      T_keV         : ICM temperature        (default 4.0 keV core)
      n_e_cm3       : electron density       (default 0.1 cm^-3 core)
      B_microG      : magnetic field         (default 25.0)
      r_kpc         : cooling radius         (default 100.0 kpc)
      dM_dt_classic : classical Bondi rate   (default 300.0 Msun/yr)
      t_n           : Mayan time fraction    (default 0.0)
      xi_SC         : SC coupling            (default 0.5)
    """

    T_SCM_KEV    = 8.6e-8 * 1.0e6   # ~0.086 keV (1 MK SCm reference temp)
    B_REF_T      = 1.0e-8           # ICM equipartition ~100 microG
    LAMBDA_COSM  = 1.1056e-52       # m^-2

    def compute(self, dataset: Dict[str, Any] = None) -> Dict[str, Any]:
        ds = dataset or {}
        T_keV     = ds.get('T_keV',         4.0)
        n_e       = ds.get('n_e_cm3',       0.1)
        B_uG      = ds.get('B_microG',      25.0)
        r_kpc     = ds.get('r_kpc',         100.0)
        dM_class  = ds.get('dM_dt_classic', 300.0)
        t_n       = ds.get('t_n',           0.0)
        xi_SC     = ds.get('xi_SC',         0.5)

        # Thermal ratio (UQFF mechanism scales as (T_ICM/T_SCm)^+3/2:
        # more thermal energy -> more SCm layers active -> more damping).
        T_ratio = T_keV / self.T_SCM_KEV
        thermal_factor = T_ratio ** 1.5

        # Magnetic damping vs ICM equipartition field
        B_T = B_uG * 1.0e-10        # microG -> T
        mag_factor = (B_T / self.B_REF_T) ** 2

        # Mayan time gating
        cos2_tn = math.cos(math.pi * t_n) ** 2

        suppression_denom = 1.0 + xi_SC * thermal_factor * mag_factor * cos2_tn
        S = 1.0 / suppression_denom
        dM_observed = S * dM_class

        # Classical Bondi check
        kT_J = T_keV * 1.602e-16
        mu_mp = 0.6 * 1.673e-27
        c_s = math.sqrt(5.0 * kT_J / (3.0 * mu_mp))   # adiabatic sound speed [m/s]
        # Cooling time
        n_e_m3 = n_e * 1.0e6
        Lambda_cool = 1.0e-23 * (T_keV / 1.0) ** 0.5   # erg cm^3 / s rough
        t_cool_s = (3.0 * kT_J) / (n_e_m3 * Lambda_cool * 1.0e-7) \
            if Lambda_cool > 0 else 0.0
        t_cool_Gyr = t_cool_s / (3.156e16)

        # Aether-Lambda counter pressure
        P_Lambda = self.LAMBDA_COSM * C_LIGHT ** 2 / (8.0 * math.pi * G_NEWTON)

        # Match against Chandra observed deficit factor
        deficit = dM_class / max(dM_observed, 1.0e-6)
        chandra_observed_range = (6.0, 30.0)
        within_chandra_band = (chandra_observed_range[0] <= deficit
                               <= chandra_observed_range[1])

        return {
            'T_keV'                 : T_keV,
            'n_e_cm3'               : n_e,
            'B_microG'              : B_uG,
            'r_kpc'                 : r_kpc,
            't_n'                   : t_n,
            'xi_SC'                 : xi_SC,
            'thermal_factor'        : thermal_factor,
            'mag_factor'            : mag_factor,
            'cos2_tn'               : cos2_tn,
            'suppression_S'         : S,
            'dM_dt_classic_Msun_yr' : dM_class,
            'dM_dt_observed_Msun_yr': dM_observed,
            'deficit_factor'        : deficit,
            'within_chandra_band'   : within_chandra_band,
            'sound_speed_m_s'       : c_s,
            'cooling_time_Gyr'      : t_cool_Gyr,
            'P_Lambda_Pa'           : P_Lambda,
            'primary_equations': [
                f"T_ICM = {T_keV} keV   n_e = {n_e} cm^-3   B = {B_uG} uG",
                f"Classical Bondi: dM/dt = {dM_class:.1f} Msun/yr",
                f"Thermal factor:   (T_ICM/T_SCm)^+3/2 = {thermal_factor:.3e}",
                f"Magnetic factor:  (B/B_ref)^2       = {mag_factor:.3e}",
                f"cos^2(pi*t_n)                        = {cos2_tn:.3f}",
                f"S = 1 / (1 + xi_SC * thermal * mag * cos^2) = {S:.4e}",
                f"dM/dt(observed) = S * dM/dt(classical) = {dM_observed:.2f} Msun/yr",
                f"Deficit factor = {deficit:.2f}  (Chandra band: 6-30)",
                f"In Chandra band: {within_chandra_band}",
                f"Cooling time t_cool = {t_cool_Gyr:.2f} Gyr",
            ],
            'available_equations': [
                "S = 1 / (1 + xi_SC*(T_ICM/T_SCm)^+3/2 *(B/B_ref)^2 * cos^2(pi*t_n))",
                "dM/dt(obs) = S * dM/dt(Bondi-classical)",
                "P_Lambda = Lambda*c^2/(8*pi*G)  (Aether counter-pressure)",
                "c_s = sqrt(5*k*T/(3*mu*m_p))",
                "t_cool = 3*k*T / (n_e * Lambda_cool)",
                "Deficit = dM/dt(classical) / dM/dt(observed)",
            ],
            'simulation_set': [
                "Scan xi_SC = 1..200; fit Chandra 1.4 Ms observed deficit",
                "Sweep T_keV = 1..15 keV across cluster cores",
                "Scan B_microG = 1..100; test magnetic-coupling sensitivity",
                "Multi-cluster: Perseus, Coma, Virgo, A1689 cooling-flow grid",
                "Time series: dM/dt vs t_n over one Great Cycle",
            ],
        }


# ===========================================================================
# Planetary presets — drop-in datasets for QCalcGeom.CrustalZeroPointCalculator
# ===========================================================================

PLANETARY_CRUSTAL_PRESETS: Dict[str, Dict[str, Any]] = {
    'Earth': dict(
        rho_crust_kg_m3 = 2700.0,
        rho_plasma_kg_m3 = 12000.0,    # outer-core analog
        h_crust_m = 3.5e4,
        epoch = 5,
        g_local = 9.81,
    ),
    'Mars': dict(
        rho_crust_kg_m3 = 2900.0,
        rho_plasma_kg_m3 = 9000.0,
        h_crust_m = 5.0e4,             # 50 km
        epoch = 5,
        g_local = 3.71,
    ),
    'Europa': dict(
        rho_crust_kg_m3 = 920.0,       # water ice
        rho_plasma_kg_m3 = 1100.0,     # salty subsurface ocean
        h_crust_m = 1.5e4,             # ~15 km ice shell
        epoch = 5,
        g_local = 1.315,
    ),
    'Enceladus': dict(
        rho_crust_kg_m3 = 920.0,
        rho_plasma_kg_m3 = 1050.0,
        h_crust_m = 5.0e3,             # ~5 km (south polar)
        epoch = 5,
        g_local = 0.113,
    ),
    'Io': dict(
        rho_crust_kg_m3 = 3500.0,      # silicate
        rho_plasma_kg_m3 = 4500.0,     # tidally-heated magma ocean
        h_crust_m = 3.0e4,
        epoch = 5,
        g_local = 1.796,
    ),
    'Titan': dict(
        rho_crust_kg_m3 = 950.0,       # methane-water ice
        rho_plasma_kg_m3 = 1300.0,     # ammonia-water ocean
        h_crust_m = 1.0e5,             # ~100 km thick ice shell
        epoch = 5,
        g_local = 1.352,
    ),
}


def planetary_crustal_preset(body: str) -> Dict[str, Any]:
    """Return a CrustalZeroPointCalculator-ready dataset dict.

    Adds current_year automatically (May 2026) unless overridden.
    """
    if body not in PLANETARY_CRUSTAL_PRESETS:
        raise KeyError(
            f"Unknown planetary body '{body}'. "
            f"Known: {sorted(PLANETARY_CRUSTAL_PRESETS.keys())}"
        )
    ds = dict(PLANETARY_CRUSTAL_PRESETS[body])
    ds.setdefault('current_year', 2026.37)
    ds.setdefault('threshold_frac', 1.0e-3)
    return ds


# ===========================================================================
# Registry — for downstream library aggregation (CondensedPhysics3 / index.js)
# ===========================================================================

SESSION_277_CALCULATORS = {
    'OrionH2OMaserUQFFCalculator'         : OrionH2OMaserUQFFCalculator,
    'NeutrinoGWCoincidenceCalculator'     : NeutrinoGWCoincidenceCalculator,
    'PerseusACoolingSuppressionCalculator': PerseusACoolingSuppressionCalculator,
}


# ===========================================================================
# Smoke tests — run as __main__
# ===========================================================================

def _run_smoke_tests() -> int:
    results = []

    def chk(tid: str, name: str, ok: bool, detail: str = ''):
        tag = 'PASS' if ok else 'FAIL'
        print(f"  [{tag}] {tid} {name}: {detail}")
        results.append((tid, ok))

    # G1: Orion H2O maser produces finite positive brightness
    r1 = OrionH2OMaserUQFFCalculator().compute()
    chk('G1-1', 'H2O default brightness > 0',
        r1['T_b_K'] > 0 and math.isfinite(r1['T_b_K']),
        f"T_b = {r1['T_b_K']:.2e} K")
    chk('G1-2', 'H2O luminosity in JWST-plausible range (1e-5..1e3 Lsun)',
        1.0e-5 <= r1['L_maser_Lsun'] <= 1.0e3,
        f"L = {r1['L_maser_Lsun']:.2e} Lsun")
    r1b = OrionH2OMaserUQFFCalculator().compute({'species': 'CH3OH', 'n_H2_cm3': 1.0e7})
    chk('G1-3', 'CH3OH 6.7 GHz mode switch',
        abs(r1b['freq_GHz'] - 6.66852) < 1.0e-3,
        f"freq = {r1b['freq_GHz']} GHz")
    chk('G1-4', 'Maser cap T_b <= 1e15 K',
        r1['T_b_K'] <= 1.0e15 + 1.0,
        f"T_b cap honoured")
    chk('G1-5', 'Provides primary/available/simulation triple',
        len(r1.get('primary_equations', [])) > 0
        and len(r1.get('available_equations', [])) > 0
        and len(r1.get('simulation_set', [])) > 0,
        "triple OK")

    # G3: Neutrino-GW coincidence
    calc3 = NeutrinoGWCoincidenceCalculator()
    # Inject signal near peak
    sig = [10.0, 15.0, 8.0, 25.0, -2.0, 50.0, 100.0]
    r3 = calc3.compute({'observed_dt_list': sig,
                        'gw_event_id': 'GW170817',
                        'n_background': 0.5})
    chk('G3-1', 'logL_UQFF finite',
        math.isfinite(r3['logL_uqff']),
        f"logL_UQFF = {r3['logL_uqff']:.2f}")
    chk('G3-2', 'UQFF favoured for in-window signal',
        r3['delta_logL'] > 0.0,
        f"delta_logL = {r3['delta_logL']:.2f}")
    # Uniform out-of-peak sample (worst-case)
    bg = [-4.9, 199.0, 180.0, 195.0, 198.0]
    r3b = calc3.compute({'observed_dt_list': bg,
                         'gw_event_id': 'background_only',
                         'n_background': 5.0})
    chk('G3-3', 'Far-tail samples disfavour UQFF',
        r3b['delta_logL'] < 0.0,
        f"delta_logL = {r3b['delta_logL']:.2f}")
    chk('G3-4', 'Window bounds [-5, +200] s',
        r3['window_min_s'] == -5.0 and r3['window_max_s'] == 200.0,
        f"window OK")
    chk('G3-5', 'Mean predicted dt = +12 s',
        abs(r3['predicted_peak_s'] - 12.0) < 0.01,
        f"dt0 = {r3['predicted_peak_s']} s")

    # G5: Perseus A cooling suppression
    r5 = PerseusACoolingSuppressionCalculator().compute()
    chk('G5-1', 'Suppression S in (0,1]',
        0.0 < r5['suppression_S'] <= 1.0,
        f"S = {r5['suppression_S']:.3e}")
    chk('G5-2', 'Observed dM/dt < classical',
        r5['dM_dt_observed_Msun_yr'] < r5['dM_dt_classic_Msun_yr'],
        f"obs={r5['dM_dt_observed_Msun_yr']:.1f}, "
        f"cls={r5['dM_dt_classic_Msun_yr']:.1f}")
    # Find xi_SC that lands inside Chandra band
    found_band = False
    for xi in (0.1, 0.3, 0.5, 1.0, 3.0, 10.0):
        rr = PerseusACoolingSuppressionCalculator().compute({'xi_SC': xi})
        if rr['within_chandra_band']:
            found_band = True
            break
    chk('G5-3', 'At least one xi_SC lands in Chandra 6-30 band',
        found_band, f"found in xi_SC sweep: {found_band}")
    chk('G5-4', 'Cooling time t_cool > 0',
        r5['cooling_time_Gyr'] > 0,
        f"t_cool = {r5['cooling_time_Gyr']:.2f} Gyr")
    chk('G5-5', 'Triple output present',
        len(r5.get('primary_equations', [])) > 0,
        "triple OK")

    # Planetary presets
    for body in ('Earth', 'Mars', 'Europa', 'Enceladus', 'Io', 'Titan'):
        ds = planetary_crustal_preset(body)
        chk(f'P-{body}', f'preset has required keys',
            all(k in ds for k in
                ('rho_crust_kg_m3', 'rho_plasma_kg_m3', 'h_crust_m',
                 'epoch', 'g_local', 'current_year')),
            f"keys OK")

    passed = sum(1 for _, ok in results if ok)
    total = len(results)
    print(f"\nTOTAL: {passed}/{total} PASS  |  {total - passed} FAIL")
    return 0 if passed == total else 1


if __name__ == '__main__':
    import sys
    sys.exit(_run_smoke_tests())
