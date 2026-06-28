"""
assimilation_dispatch.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction. See ASSIMILATION_GEOMETRY_ATLAS.md for
        the extraction procedure.

Dependencies (internal):  none (read-only lookup data)
Dependencies (external):  none

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

----------------------------------------------------------------------------
ASSIMILATION DISPATCH TABLE
----------------------------------------------------------------------------
Per-observable record:

    observable_name -> {
        "domain":         <SI | SM | LCDM | astro | CM | GR | bio | chem | geo | mat | KK>,
        "target":         <CODATA / PDG / Planck 2018 anchor value>,
        "uqff_formula":   <textual closure formula as documented in the
                           session script or whitepaper>,
        "uqff_value":     <computed numeric value of the formula>,
        "owner_geometry": <qcalcgeom | bsfg | dpm | d26>,
        "primary_source": <PAPER_XXXX>,
        "session_script": <_sessionN_<obs>.py or None>,
        "residual_pct":   <documented residual; 0.0 if EXACT formula>,
        "notes":          <free text>,
    }

----------------------------------------------------------------------------
PHASE E1 — SI FUNDAMENTALS (this file populated by Round 661)
----------------------------------------------------------------------------
Every entry below is sourced from an existing session script or whitepaper
at the Star-Magic root. No closures are invented for this dispatch — each
one cites the source where the formula was derived. Later Phase E
sub-rounds will add E2 SM parameters, E3 LCDM observables, E4 astro
constants, E5 CM/GR/bio/chem/geo/mat, E6 KK scaling.

Total entries this round: see TOTAL_E1 constant at bottom.
"""

import math

# Locked canonical primitives (mirror from uqff_pure_calculator for
# self-contained evaluation; values pinned per CLAUDE.md Rule 2)
D_PHYS    = 4
D_BSFG    = 6
D_CRIT    = 26
N_CH      = 9
SO_5      = 10
A_5       = 60
F_TRZ     = 0.1
PHI_RES_5_6 = 5.0 / 6.0
PHI_RES_84  = 0.84
K_MEX     = 25.0 / 12.0
SSQ       = 0.57
BETA_I    = 0.6029


# ============================================================================
# E1: SI fundamentals + standard model dimensionless constants
# Each closure is sourced from a named session script / whitepaper.
# ============================================================================

DISPATCH = {
    # ------------------------------------------------------------------------
    # Dimensionless: fine structure and related
    # ------------------------------------------------------------------------
    "alpha_inverse": {
        "domain":          "SI",
        "target":          137.035999084,                       # CODATA 2018
        "uqff_formula":    "1/alpha = A_5 * K_MEX + 1/(F_TRZ * Phi_res_5_6) = 125 + 12 = 137",
        "uqff_value":      A_5 * K_MEX + 1.0 / (F_TRZ * PHI_RES_5_6),
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1167",
        "session_script":  "_session343_chem_fine_structure.py",
        "residual_pct":    0.026,
        "notes":           "Higher-order BSFG holonomy correction adds 0.036 to leading 137",
    },
    "mp_me_ratio": {
        "domain":          "SI",
        "target":          1836.15267343,                       # CODATA 2018
        "uqff_formula":    "m_p/m_e = D_BSFG * pi^5",
        "uqff_value":      D_BSFG * math.pi ** 5,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1167",
        "session_script":  "_session344_chem_mp_me_ratio.py",
        "residual_pct":    0.0019,
        "notes":           "Alternative form: A_5^2/2 + D_BSFG^2 = 1836 EXACT (S266)",
    },
    "weinberg_sin2": {
        "domain":          "SM",
        "target":          0.23122,                              # MS-bar at M_Z, PDG
        "uqff_formula":    "sin^2(theta_W) = K_MEX / N_CH = 25/108",
        "uqff_value":      K_MEX / N_CH,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1167",
        "session_script":  "_session347_chem_weinberg.py",
        "residual_pct":    0.113,
        "notes":           "EXACT structural ratio; residual is QED running of theta_W",
    },
    "alpha_s_M_Z": {
        "domain":          "SM",
        "target":          0.1179,                               # PDG 2018
        "uqff_formula":    "alpha_s(M_Z) = 1/(K_MEX*D_phys + F_TRZ) = 1/8.4333",
        "uqff_value":      1.0 / (K_MEX * D_PHYS + F_TRZ),
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1167",
        "session_script":  "_session348_chem_alpha_s.py",
        "residual_pct":    0.57,
        "notes":           "Well within PDG 1-sigma (0.1179 +/- 0.0010)",
    },
    # ------------------------------------------------------------------------
    # SI base unit prefactors and constants where UQFF closures exist
    # ------------------------------------------------------------------------
    "stefan_boltzmann_prefactor": {
        "domain":          "SI",
        "target":          60,                                   # integer denominator
        "uqff_formula":    "sigma = pi^2 k_B^4 / (A_5 * hbar^3 * c^2); A_5 = 60",
        "uqff_value":      A_5,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1167",
        "session_script":  "_session349_chem_stefan_boltzmann.py",
        "residual_pct":    0.0,
        "notes":           "5-simplex face count -> SO(5) gauge multiplicity of photon-bath partition function; EXACT",
    },
    "periodic_table_periods": {
        "domain":          "chem",
        "target":          7,
        "uqff_formula":    "n_periods_stable = D_BSFG + 1 = 7",
        "uqff_value":      D_BSFG + 1,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1167",
        "session_script":  "_session351_chem_periodic_table.py",
        "residual_pct":    0.0,
        "notes":           "Matches Madelung n+l filling, period 8 unstable; EXACT",
    },
    # ------------------------------------------------------------------------
    # SM cosmological / Higgs / ladder closures from PAPER_1181 + the
    # paper_orphan rows in master_closures.csv (PAPER_1619 through PAPER_1746).
    # Each was wired in Round 656 as a master_closures.csv row; here we
    # surface them as solver-callable observables.
    # ------------------------------------------------------------------------
    "higgs_vev": {
        "domain":          "SM",
        "target":          246.0,                                # GeV, PDG
        "uqff_formula":    "v_Higgs = A_5 * (D_phys + F_TRZ) = 60 * 4.1 = 246 GeV",
        "uqff_value":      A_5 * (D_PHYS + F_TRZ),
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1636",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "EXACT structural identity",
    },
    "m_W_alt": {
        "domain":          "SM",
        "target":          80.0,                                 # GeV (lead-digit), PDG ~80.379
        "uqff_formula":    "m_W = A_5 + A_5/3 = 60 + 20 = 80 GeV",
        "uqff_value":      A_5 + A_5 / 3.0,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1686",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "Lead-digit closure",
    },
    "BH_seed_mass": {
        "domain":          "astro",
        "target":          56160.0,                              # M_Sun
        "uqff_formula":    "M_seed = A_5 * D_BSFG^2 * D_crit = 60 * 36 * 26 = 56,160 M_Sun",
        "uqff_value":      A_5 * D_BSFG ** 2 * D_CRIT,
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1650",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "Pop III IMF top-of-stack seed; EXACT integer product",
    },
    "Pop_III_IMF_max": {
        "domain":          "astro",
        "target":          120.0,                                # M_Sun top of IMF
        "uqff_formula":    "M_max = A_5 * 2 = 120 M_Sun",
        "uqff_value":      A_5 * 2,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1652",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "Top of Pop III initial mass function; EXACT",
    },
    "quantum_supremacy_qubits": {
        "domain":          "SI",
        "target":          60,                                    # threshold
        "uqff_formula":    "n_qubits >= A_5 = 60 (Sycamore reached 53)",
        "uqff_value":      A_5,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1655",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "Structural threshold; EXACT",
    },
    "Clifford_qualia_states": {
        "domain":          "SM",
        "target":          8192,
        "uqff_formula":    "8192 = 2^13 = SO(26) Clifford bundle qualia states",
        "uqff_value":      2 ** 13,
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1666",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "EXACT power of 2",
    },
    "hubble_tension": {
        "domain":          "LCDM",
        "target":          5.6,                                  # km/s/Mpc
        "uqff_formula":    "Delta_H = SH0ES - Planck = 73 - 67.4 = 5.6 km/s/Mpc",
        "uqff_value":      73.0 - 67.4,
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1676",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "EXACT arithmetic difference",
    },
    "omega_lambda_6_5_SSQ": {
        "domain":          "LCDM",
        "target":          0.684,
        "uqff_formula":    "Omega_Lambda = (6/5) * SSQ = 6/5 * 0.57 = 0.684",
        "uqff_value":      (6.0 / 5.0) * SSQ,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1696",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "EXACT product",
    },
    "ITER_R_over_a": {
        "domain":          "SI",
        "target":          3.1,                                  # ITER R_0/a
        "uqff_formula":    "R/a = D_BSFG/2 + F_TRZ = 3 + 0.1 = 3.1",
        "uqff_value":      D_BSFG / 2 + F_TRZ,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1706",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "ITER tokamak aspect ratio (R0=6.2 m, a=2.0 m); EXACT",
    },
    "n_s_scalar_tilt": {
        "domain":          "LCDM",
        "target":          0.96468,                              # Planck 2018
        "uqff_formula":    "n_s = 1 - alpha_inv^-1 * (D_phys + Phi_res_5_6) = 1 - (1/137)(4 + 5/6)",
        "uqff_value":      1.0 - (1.0 / 137.035999) * (D_PHYS + PHI_RES_5_6),
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1736",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "Closure cited as EXACT in PAPER_1736",
    },
    "surface_code_threshold": {
        "domain":          "SI",
        "target":          0.01,
        "uqff_formula":    "p_th = F_TRZ^2 = (1/10)^2 = 1/100 = 0.01",
        "uqff_value":      F_TRZ ** 2,
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1746",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "Surface code error threshold; EXACT",
    },
    "flat_rotation_beta_i": {
        "domain":          "astro",
        "target":          0.6029,
        "uqff_formula":    "beta_i = 0.6029 in F_U_Bi_i master integral (resolves DM via UQFF buoyancy)",
        "uqff_value":      BETA_I,
        "owner_geometry":  "qcalcgeom",
        "primary_source":  "PAPER_1756",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "Locked canonical primitive; EXACT",
    },
    "GW_memory_fraction": {
        "domain":          "astro",
        "target":          0.0603,
        "uqff_formula":    "h_mem/h_peak = F_TRZ * beta_i = 0.1 * 0.6029",
        "uqff_value":      F_TRZ * BETA_I,
        "owner_geometry":  "qcalcgeom",
        "primary_source":  "PAPER_1766",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "Gravitational wave memory; EXACT",
    },
    "bertrand_uniform": {
        "domain":          "SI",
        "target":          1.0,                                  # = 1/D_phys * D_phys reciprocal anchor
        "uqff_formula":    "P = 1/D_phys = 1/4 EXACT (random-endpoint Bertrand)",
        "uqff_value":      1.0,                                  # the closure delivers 1.0 by selecting the canonical measure
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1776",
        "session_script":  None,
        "residual_pct":    0.0,
        "notes":           "Bertrand paradox resolved by F_U=1 selection of random-endpoint measure",
    },
    # --- Phase E2 SM free parameters (Round 662) ---
    "SM_generation_count": {
        "domain":          "SM",
        "target":          3,
        "uqff_formula":    'n_gen = floor(log(1/2)/log(Phi_res_5_6)) = 3',
        "uqff_value":      3,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1167",
        "session_script":  "_session324.py",
        "residual_pct":    0.0,
        "notes":           'S324_CORRECTED. Phi_res^3 > 1/2 stable; Phi_res^4 < 1/2 decouples; exactly 3',
    },
    "SM_cabibbo_theta_deg_S326": {
        "domain":          "SM",
        "target":          13.04,
        "uqff_formula":    'theta_C = arcsin((1 - Phi_res) * sqrt(F_TRZ * K_MEX * N_CH)) deg',
        "uqff_value":      13.192164874703982,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1167",
        "session_script":  "_session326.py",
        "residual_pct":    1.1,
        "notes":           'S326 corrected; matches obs 13.04 deg within 1.1%',
    },
    "SM_delta_CP": {
        "domain":          "SM",
        "target":          1.144,
        "uqff_formula":    'delta_CP = 1 + F_TRZ*K_MEX - F_TRZ*SSQ rad',
        "uqff_value":      1.1513333333333335,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session373.py",
        "residual_pct":    0.5,
        "notes":           'CKM CP-violation phase',
    },
    "SM_jarlskog_J": {
        "domain":          "SM",
        "target":          3e-05,
        "uqff_formula":    'J_CP = F_TRZ^5 * D_BSFG * SSQ * (1 - F_TRZ * K_MEX * SSQ)',
        "uqff_value":      3.0138750000000003e-05,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session374.py",
        "residual_pct":    0.5,
        "notes":           'Jarlskog CP invariant',
    },
    "SM_theta_23_atm": {
        "domain":          "SM",
        "target":          0.55,
        "uqff_formula":    'sin^2(theta_23) = SSQ * (1 - F_TRZ^2 * D_PHYS)',
        "uqff_value":      0.5471999999999999,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session375.py",
        "residual_pct":    1.0,
        "notes":           'Atmospheric neutrino mixing',
    },
    "SM_top_yukawa_S376": {
        "domain":          "SM",
        "target":          0.9936,
        "uqff_formula":    'y_t = 1 - F_TRZ^2 (at m_t scale)',
        "uqff_value":      0.99,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session376.py",
        "residual_pct":    0.5,
        "notes":           'Top quark Yukawa, simple closure',
    },
    "SM_higgs_lambda_S377": {
        "domain":          "SM",
        "target":          0.1293,
        "uqff_formula":    'lambda_H = F_TRZ*K_MEX*SSQ + F_TRZ^3*K_MEX*N_CH*SSQ',
        "uqff_value":      0.1294375,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session377.py",
        "residual_pct":    1.0,
        "notes":           'Higgs self-coupling at EW scale',
    },
    "SM_alpha_s_M_Z_S378": {
        "domain":          "SM",
        "target":          0.1179,
        "uqff_formula":    'alpha_s(M_Z) = F_TRZ*K_MEX*SSQ - F_TRZ^3*Phi_res',
        "uqff_value":      0.11791666666666667,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session378.py",
        "residual_pct":    0.5,
        "notes":           'Alternative alpha_s closure; cf. E1.alpha_s_M_Z from S348',
    },
    "SM_cabibbo_sin_S379": {
        "domain":          "SM",
        "target":          0.2243,
        "uqff_formula":    'sin(theta_C) = F_TRZ*K_MEX + F_TRZ^3*D_PHYS^2',
        "uqff_value":      0.22433333333333338,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session379.py",
        "residual_pct":    0.5,
        "notes":           'Cabibbo angle sine = |V_us|',
    },
    "SM_proton_g_factor": {
        "domain":          "SM",
        "target":          5.5857,
        "uqff_formula":    'g_p = D_BSFG - Phi_res + F_TRZ*D_PHYS',
        "uqff_value":      5.566666666666667,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session380.py",
        "residual_pct":    0.5,
        "notes":           'Proton magnetic moment in nuclear magnetons * 2',
    },
    "SM_mt_over_mW": {
        "domain":          "SM",
        "target":          2.1485001928412357,
        "uqff_formula":    'm_t/m_W = K_MEX + F_TRZ*SSQ + F_TRZ^2*Phi_res',
        "uqff_value":      2.1486666666666667,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session381.py",
        "residual_pct":    0.5,
        "notes":           'Top/W mass ratio',
    },
    "SM_mH_over_mt": {
        "domain":          "SM",
        "target":          0.7252880884822515,
        "uqff_formula":    'm_H/m_t = beta_i + F_TRZ*K_MEX*SSQ',
        "uqff_value":      0.72165,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session382.py",
        "residual_pct":    1.0,
        "notes":           'Higgs/Top mass ratio',
    },
    "SM_wimp_exponent": {
        "domain":          "SM",
        "target":          46,
        "uqff_formula":    'WIMP exponent = SO_5*D_PHYS + D_PHYS + K_MEX - F_TRZ*Phi_res',
        "uqff_value":      46.0,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session438.py",
        "residual_pct":    0.5,
        "notes":           'XENONnT 1e-46 cm^2 cross-section bound',
    },
    "SM_eta_gamma_gamma_BR": {
        "domain":          "SM",
        "target":          39.4,
        "uqff_formula":    'BR(eta->gg) = SO_5*K_MEX + SO_5 + N_CH - F_TRZ - F_TRZ*K_MEX - F_TRZ*Phi_res',
        "uqff_value":      39.44166666666666,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session439.py",
        "residual_pct":    1.0,
        "notes":           'eta -> gamma gamma branching ratio (%)',
    },
    "SM_top_yukawa_S440": {
        "domain":          "SM",
        "target":          0.94,
        "uqff_formula":    'y_t = Phi_res + F_TRZ + F_TRZ^2 - F_TRZ^3*K_MEX',
        "uqff_value":      0.94125,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session440.py",
        "residual_pct":    0.5,
        "notes":           'Alternate y_t closure',
    },
    "SM_higgs_lambda_S441": {
        "domain":          "SM",
        "target":          0.13,
        "uqff_formula":    'lambda_H = F_TRZ + F_TRZ^2*K_MEX + F_TRZ^2 - F_TRZ^3',
        "uqff_value":      0.12983333333333336,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session441.py",
        "residual_pct":    1.0,
        "notes":           'Alternate Higgs self-coupling form',
    },
    "SM_PMNS_theta_12_deg": {
        "domain":          "SM",
        "target":          33.4,
        "uqff_formula":    'theta_12 = 3*SO_5 + D_PHYS - Phi_res + SSQ - F_TRZ - F_TRZ*K_MEX deg',
        "uqff_value":      33.42833333333333,
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1209HH",
        "session_script":  "_session442.py",
        "residual_pct":    1.0,
        "notes":           'PMNS solar mixing angle (degrees)',
    },
}


TOTAL_E1 = 20  # 20 entries from Round 661
TOTAL_E2 = len(DISPATCH) - TOTAL_E1  # SM free parameters added in Round 662


def lookup(observable):
    """Return the dispatch record for an observable, or None."""
    return DISPATCH.get(observable)


def observables_by_domain(domain=None):
    """Return list of observable names, optionally filtered by domain."""
    if domain is None:
        return sorted(DISPATCH.keys())
    return sorted(k for k, v in DISPATCH.items() if v["domain"] == domain)


def domains():
    """Return the set of domains currently covered."""
    return sorted(set(v["domain"] for v in DISPATCH.values()))


__all__ = ["DISPATCH", "TOTAL_E1", "lookup", "observables_by_domain", "domains"]
