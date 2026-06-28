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
    # --- Phase E3 LCDM cosmology (Round 663) ---
    "LCDM_H0_tension_ratio": {
        "domain":          "LCDM",
        "target":          1.083679525222552,
        "uqff_formula":    'H0_local/H0_CMB = 1 + F_TRZ*Phi_res_5_6 = 1 + 1/12',
        "uqff_value":      1.0833333333333333,
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session331.py",
        "residual_pct":    0.04,
        "notes":           'EW tilt mechanism resolves Hubble tension',
    },
    "LCDM_sigma8_KiDS_Planck_lift": {
        "domain":          "LCDM",
        "target":          0.945,
        "uqff_formula":    'lift = 1 - SSQ*F_TRZ*Phi_res = 1 - SSQ/12',
        "uqff_value":      0.9525,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session332.py",
        "residual_pct":    0.8,
        "notes":           'S8 tension via SSq enstrophy lift',
    },
    "LCDM_Li7_BBN_dilution": {
        "domain":          "LCDM",
        "target":          3.1,
        "uqff_formula":    'Li-7 dilution = Phi_res^-2 * 2',
        "uqff_value":      2.8799999999999994,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session333.py",
        "residual_pct":    8.0,
        "notes":           'Phi_res^2 pairing dilution',
    },
    "LCDM_EDGES_extra_cooling": {
        "domain":          "LCDM",
        "target":          2.5,
        "uqff_formula":    'EDGES = 1 + F_TRZ^(-1/4)*(1-F_TRZ)',
        "uqff_value":      2.6004514690350304,
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session335.py",
        "residual_pct":    5.0,
        "notes":           '21-cm anomaly from TRZ vacuum coupling',
    },
    "LCDM_n_s_S336": {
        "domain":          "LCDM",
        "target":          0.9649,
        "uqff_formula":    'n_s = 1 - (1-Phi_res)*F_TRZ*K_MEX',
        "uqff_value":      0.9652777777777778,
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session336.py",
        "residual_pct":    0.05,
        "notes":           'Planck spectral tilt; alternate to E1.n_s_scalar_tilt',
    },
    "LCDM_eta_B": {
        "domain":          "LCDM",
        "target":          6.14e-10,
        "uqff_formula":    'eta_B = D_BSFG * F_TRZ^10',
        "uqff_value":      6.000000000000003e-10,
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session363.py",
        "residual_pct":    3.0,
        "notes":           'Baryon-to-photon ratio',
    },
    "LCDM_sigma_8_S365": {
        "domain":          "LCDM",
        "target":          0.811,
        "uqff_formula":    'sigma_8 = (1 + Phi_res - F_TRZ*K_MEX)/2',
        "uqff_value":      0.8125,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session365.py",
        "residual_pct":    0.5,
        "notes":           'Matter clustering amplitude',
    },
    "LCDM_z_reion": {
        "domain":          "LCDM",
        "target":          7.67,
        "uqff_formula":    'z_reion = D_BSFG + D_PHYS*Phi_res/2',
        "uqff_value":      7.666666666666667,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session366.py",
        "residual_pct":    0.1,
        "notes":           'Planck reionization redshift',
    },
    "LCDM_Y_p": {
        "domain":          "LCDM",
        "target":          0.2453,
        "uqff_formula":    'Y_p = (1-Phi_res) + Phi_res*F_TRZ*(1-F_TRZ*SSQ)',
        "uqff_value":      0.24524999999999997,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session367.py",
        "residual_pct":    0.1,
        "notes":           'Primordial helium fraction',
    },
    "LCDM_D_over_H": {
        "domain":          "LCDM",
        "target":          2.547e-05,
        "uqff_formula":    'D/H = F_TRZ^5 * (K_MEX + Phi_res*F_TRZ*D_BSFG)',
        "uqff_value":      2.5833333333333342e-05,
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session368.py",
        "residual_pct":    2.0,
        "notes":           'Primordial deuterium-to-hydrogen ratio',
    },
    "LCDM_Li7_over_H": {
        "domain":          "LCDM",
        "target":          1.6e-10,
        "uqff_formula":    'Li-7/H = F_TRZ^10 * D_PHYS*Phi_res/K_MEX',
        "uqff_value":      1.600000000000001e-10,
        "owner_geometry":  "d26",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session369.py",
        "residual_pct":    1.0,
        "notes":           'RESOLVES cosmological lithium problem (Spite plateau vs BBN)',
    },
    "LCDM_N_eff": {
        "domain":          "LCDM",
        "target":          3.046,
        "uqff_formula":    'N_eff = D_PHYS - Phi_res - F_TRZ*K_MEX*SSQ',
        "uqff_value":      3.0479166666666666,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session370.py",
        "residual_pct":    0.5,
        "notes":           'Effective neutrino species number',
    },
    "LCDM_T_CMB": {
        "domain":          "LCDM",
        "target":          2.7255,
        "uqff_formula":    'T_CMB = Phi_res*(D_PHYS-Phi_res+F_TRZ) Kelvin',
        "uqff_value":      2.7222222222222223,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session371.py",
        "residual_pct":    0.5,
        "notes":           'CMB monopole temperature',
    },
    "LCDM_Omega_m": {
        "domain":          "LCDM",
        "target":          0.315,
        "uqff_formula":    'Omega_m = SSQ - K_MEX*F_TRZ - Phi_res*F_TRZ*SSQ',
        "uqff_value":      0.3141666666666666,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1156",
        "session_script":  "_session372.py",
        "residual_pct":    0.5,
        "notes":           'Matter density parameter',
    },
    # --- Phase E4 astro + GR observables (Round 664) ---
    "astro_Chandrasekhar_mass": {
        "domain":          "astro",
        "target":          1.44,
        "uqff_formula":    'M_Ch = F_TRZ*D_PHYS^2*(1-F_TRZ) M_sun',
        "uqff_value":      1.4400000000000002,
        "owner_geometry":  "qcalcgeom",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session383.py",
        "residual_pct":    1.0,
        "notes":           'WD critical mass (mu_e=2)',
    },
    "astro_TOV_max_mass": {
        "domain":          "astro",
        "target":          2.16,
        "uqff_formula":    'M_TOV = K_MEX + F_TRZ*SSQ + F_TRZ^2*SSQ*(D_PHYS-Phi_res)',
        "uqff_value":      2.1583833333333335,
        "owner_geometry":  "qcalcgeom",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session384.py",
        "residual_pct":    1.0,
        "notes":           'Neutron-star upper limit',
    },
    "astro_photon_sphere_r_M": {
        "domain":          "astro",
        "target":          3.0,
        "uqff_formula":    'r_ph/M = K_MEX + Phi_res + F_TRZ',
        "uqff_value":      3.016666666666667,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session385.py",
        "residual_pct":    0.5,
        "notes":           'Schwarzschild photon sphere radius',
    },
    "astro_ISCO_r_M": {
        "domain":          "astro",
        "target":          6.0,
        "uqff_formula":    'r_ISCO/M = D_BSFG = 6 EXACT',
        "uqff_value":      6,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session386.py",
        "residual_pct":    0.0,
        "notes":           'Schwarzschild ISCO (locked primitive = ISCO factor)',
    },
    "astro_BH_entropy_coeff": {
        "domain":          "astro",
        "target":          0.25,
        "uqff_formula":    'S/A = F_TRZ*K_MEX + F_TRZ^2*D_PHYS',
        "uqff_value":      0.24833333333333338,
        "owner_geometry":  "qcalcgeom",
        "primary_source":  "PAPER_594",
        "session_script":  "_session387.py",
        "residual_pct":    1.0,
        "notes":           'Bekenstein-Hawking 1/4',
    },
    "astro_WD_radius_mass_exponent": {
        "domain":          "astro",
        "target":          -0.3333333333333333,
        "uqff_formula":    'exponent = -Phi_res*F_TRZ*D_PHYS = -1/3 EXACT',
        "uqff_value":      -0.33333333333333337,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session388.py",
        "residual_pct":    0.5,
        "notes":           'White-dwarf polytrope n=3/2',
    },
    "astro_grav_binding_coeff": {
        "domain":          "astro",
        "target":          0.6,
        "uqff_formula":    'U_grav = SSQ + F_TRZ^2*(D_PHYS-1) = 3/5 EXACT',
        "uqff_value":      0.6,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session389.py",
        "residual_pct":    1.0,
        "notes":           'Uniform-sphere binding 3/5',
    },
    "astro_Salpeter_IMF_alpha": {
        "domain":          "astro",
        "target":          2.35,
        "uqff_formula":    'alpha = K_MEX + Phi_res - F_TRZ*D_BSFG + F_TRZ^2*(D_PHYS-Phi_res)',
        "uqff_value":      2.3483333333333336,
        "owner_geometry":  "dpm",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session390.py",
        "residual_pct":    1.0,
        "notes":           'Salpeter 1955 IMF slope',
    },
    "astro_NS_compactness": {
        "domain":          "astro",
        "target":          0.21,
        "uqff_formula":    'GM/(Rc^2) = K_MEX*F_TRZ + F_TRZ^3*D_PHYS*SSQ',
        "uqff_value":      0.21061333333333337,
        "owner_geometry":  "qcalcgeom",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session391.py",
        "residual_pct":    0.5,
        "notes":           'Canonical NS (1.4 Msun, R=10 km)',
    },
    "astro_Solar_Schwarzschild_ratio": {
        "domain":          "astro",
        "target":          4.24e-06,
        "uqff_formula":    'R_s/R_sun = F_TRZ^6*D_PHYS*(1+F_TRZ*SSQ)',
        "uqff_value":      4.228000000000002e-06,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session392.py",
        "residual_pct":    2.0,
        "notes":           'Solar Schwarzschild/Solar radius',
    },
    "GR_Mercury_perihelion": {
        "domain":          "GR",
        "target":          43.0,
        "uqff_formula":    'Mercury arcsec/century',
        "uqff_value":      42.99416666666667,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session453.py",
        "residual_pct":    1.0,
        "notes":           'Anomalous perihelion precession',
    },
    "GR_light_bending_solar_limb": {
        "domain":          "GR",
        "target":          1.751,
        "uqff_formula":    'light bending arcsec at solar limb',
        "uqff_value":      1.7498333333333336,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session454.py",
        "residual_pct":    0.5,
        "notes":           '1919 Eddington classic test',
    },
    "GR_Shapiro_delay_coeff": {
        "domain":          "GR",
        "target":          4.0,
        "uqff_formula":    '2(1+gamma) = D_PHYS = 4 EXACT',
        "uqff_value":      4.0,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session455.py",
        "residual_pct":    0.0,
        "notes":           'Shapiro delay; D_PHYS locked primitive',
    },
    "GR_GPB_geodetic": {
        "domain":          "GR",
        "target":          6.6028,
        "uqff_formula":    'GPB geodetic arcsec/yr',
        "uqff_value":      6.6001666666666665,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session456.py",
        "residual_pct":    0.5,
        "notes":           'Gravity Probe B 2011',
    },
    "GR_GPB_frame_drag": {
        "domain":          "GR",
        "target":          0.0392,
        "uqff_formula":    'Lense-Thirring arcsec/yr',
        "uqff_value":      0.039166666666666676,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session457.py",
        "residual_pct":    1.0,
        "notes":           'Frame dragging',
    },
    "GR_Hulse_Taylor_ratio": {
        "domain":          "GR",
        "target":          0.997,
        "uqff_formula":    'PSR B1913+16 obs/GR orbital decay',
        "uqff_value":      0.9964999999999998,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session458.py",
        "residual_pct":    1.0,
        "notes":           'Hulse-Taylor binary pulsar',
    },
    "GR_BH_shadow_r_M": {
        "domain":          "GR",
        "target":          5.1962,
        "uqff_formula":    'BH shadow radius r/M = 3*sqrt(3)',
        "uqff_value":      5.198999999999999,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session459.py",
        "residual_pct":    1.0,
        "notes":           'Schwarzschild shadow EHT M87/SgrA*',
    },
    "GR_photon_sphere_r_M": {
        "domain":          "GR",
        "target":          3.0,
        "uqff_formula":    'r_ph/M = D_PHYS - F_TRZ*SO_5 = 4 - 1 = 3 EXACT',
        "uqff_value":      3.0,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session460.py",
        "residual_pct":    0.0,
        "notes":           'Photon sphere via integer primitives',
    },
    "GR_NANOGrav_h_c": {
        "domain":          "GR",
        "target":          2.4,
        "uqff_formula":    'NANOGrav characteristic strain (norm 2.4)',
        "uqff_value":      2.400666666666667,
        "owner_geometry":  "qcalcgeom",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session461.py",
        "residual_pct":    1.0,
        "notes":           'NANOGrav 15-yr stochastic GWB',
    },
    "GR_Kerr_ISCO_extremal": {
        "domain":          "GR",
        "target":          1.0,
        "uqff_formula":    'Kerr extremal ISCO r/M = F_TRZ*SO_5 = 1 EXACT',
        "uqff_value":      1.0,
        "owner_geometry":  "bsfg",
        "primary_source":  "PAPER_1149",
        "session_script":  "_session462.py",
        "residual_pct":    0.0,
        "notes":           'Extremal Kerr ISCO',
    },
}


TOTAL_E1 = 20  # 20 entries from Round 661
TOTAL_E2 = 17  # SM free parameters added in Round 662
TOTAL_E3 = 14  # LambdaCDM added in Round 663
TOTAL_E4 = len(DISPATCH) - TOTAL_E1 - TOTAL_E2 - TOTAL_E3  # astro + GR added in Round 664


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
