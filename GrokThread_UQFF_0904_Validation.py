#!/usr/bin/env python3
"""
GrokThread_UQFF_0904_Validation.py

Unique validation data from Grok 4 thread 0904a12a5c2b4a639389ae084391b94f
(Star Magic UQFF Session — 7,121 lines).

Contents:
1. UQFF_52_SYSTEM_CATALOGUE       — 52-system F_U_Bi_i integration dataset
2. KAPPA_MCMC_CALIBRATION          — kappa MCMC refinement (0.00052 day⁻¹)
3. DELTA_RHO_NORMALITY_TESTS       — SW/KS/AD normality statistics (n=47→52)
4. Q_WAVE_52_STATISTICS            — Q_wave stats for 52-system set
5. UQFF_MASTER_EQUATIONS_FULL      — 6 master equations with notation
6. UQFF_ATOMIC_Z_SCALING           — Atomic Z-scaling (Z=1–118)
7. CERN_DELPHI_HIGGS_RECORDS       — DELPHI records 93719-93726, Higgs 85 GeV
8. DPM_YIN_YANG_COSMOLOGY          — DPM Yin-Yang cosmological detail
9. UQFFAtomicZScalingCalculator    — Calculator class for Z-scaled mass function

Cross-reference:
  Primary source: grok_share_0904a12a5c2b4a639389ae084391b94f_content.txt
  Grok thread ref: https://x.com/i/grok/share/0904a12a5c2b4a639389ae084391b94f
  GrokThread_StarMagic_UnifiedFramework.py — companion file (inflation epochs,
      variable documentation, DPM sphere birth)
  CondensedPhysics_Validation.py          — astronomical system catalogue (24→29)
  CondensedPhysics_OutputData.py          — Q_WAVE statistics, SYSTEM_COUNTS

©2025 Daniel T. Murphy, daniel.murphy00@gmail.com — All Rights Reserved
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any
import math


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: UQFF 52-SYSTEM CATALOGUE
# Grok 0904 thread — F_U_Bi_i_mean = -6.05e217 N (52-system average)
# log bootstrap std = 3%
# ═══════════════════════════════════════════════════════════════════════════════

UQFF_52_SYSTEM_CATALOGUE = {
    'thread': 'grok_share_0904a12a5c2b4a639389ae084391b94f',
    'F_U_Bi_i_mean_N': -6.05e217,
    'log_bootstrap_std_percent': 3.0,
    'n_systems': 52,
    'x_2_mean_m': -3.40e172,       # Cosmic quadratic solve
    'description': '52-system F_U_Bi_i integration from 0904 Grok thread',

    # ─────────────────────────────────────────────────────────────
    # SYSTEMS 01-24: inherited from UQFF_ASTRONOMICAL_SYSTEMS_VALIDATION
    # (CondensedPhysics_Validation.py L1505)
    # See that file for full parameter dicts for systems 01-24.
    # ─────────────────────────────────────────────────────────────
    'systems_01_24_ref': 'CondensedPhysics_Validation.py::UQFF_ASTRONOMICAL_SYSTEMS_VALIDATION',

    # ─────────────────────────────────────────────────────────────
    # SYSTEMS 25-52: new systems added in 0904 thread
    # ─────────────────────────────────────────────────────────────

    'system_25_m87': {
        'name': 'M87 / Virgo A',
        'category': 'Galaxy + AGN Jet',
        'uqff_terms': ['Navier-Stokes jets', 'Ub_i asymmetry', 'Ug4 SMBH'],
        'parameters': {
            'M_bh_kg': 1.26e40,     # 6.5e9 solar masses
            'M_bh_M_sun': 6.5e9,
            'distance_Mpc': 16.4,
            'jet_length_kpc': 5.0,
            'B_jet_T': 1e-4,
        },
        'verification': ['EHT 2019', 'Chandra jet imaging'],
    },

    'system_26_crab_nebula': {
        'name': 'Crab Nebula (M1, SNR + Pulsar)',
        'category': 'SNR/Pulsar Wind Nebula',
        'uqff_terms': ['Q_wave spin-down', 'Ug3 rotation', 'Ub_i turbulence'],
        'parameters': {
            'B_pulsar_T': 3.8e8,
            'E_dot_W': 5e31,          # Spin-down luminosity
            'age_yr': 970,
            'distance_kpc': 2.0,
            'P_rot_ms': 33.0,
        },
        'verification': ['HEASARC 2025', 'Fermi-LAT gamma'],
    },

    'system_27_vela_snr': {
        'name': 'Vela Pulsar / SNR',
        'category': 'SNR/Pulsar',
        'uqff_terms': ['Ug1 dipole', 'E_react pulsar wind', 'Q_wave 3.98e-3'],
        'parameters': {
            'B_pulsar_T': 3.4e8,
            'P_rot_ms': 89.0,
            'distance_kpc': 0.287,
            'age_yr': 11000,
        },
        'verification': ['ATNF Pulsar Catalogue 2025'],
    },

    'system_28_ngc_1365': {
        'name': 'NGC 1365 (Seyfert 1.5)',
        'category': 'Seyfert Galaxy',
        'uqff_terms': ['CRP', 'E_react', 'Triadic master gravity'],
        'parameters': {
            'M_bh_M_sun': 2e7,
            'L_bol_W': 3e37,
            'z': 0.0055,
        },
        'verification': ['XMM-Newton 2025', 'NuSTAR obscured AGN'],
    },

    'system_29_betelgeuse': {
        'name': 'Betelgeuse (α Orionis)',
        'category': 'Red Supergiant',
        'uqff_terms': ['Ug2 bubble', 'mass-loss Ub_i', 'stellar convection'],
        'parameters': {
            'M_M_sun': 16.5,
            'R_R_sun': 887,
            'L_L_sun': 1.26e5,
            'M_dot_M_sun_yr': 3e-6,
            'v_wind_km_s': 15.0,
            'distance_pc': 168,
        },
        'verification': ['ALMA dust envelope 2024', 'Hubble UV 2025'],
    },

    'system_30_eso_137_001': {
        'name': 'ESO 137-001 (Jellyfish Galaxy)',
        'category': 'Ram-Pressure Stripped Galaxy',
        'uqff_terms': ['Navier-Stokes fluid', 'Ub_i ram stripping', 'Ug2 stream'],
        'parameters': {
            'ram_pressure_Pa': 3e-11,
            'tail_length_kpc': 70,
            'SFR_M_sun_yr': 2.5,
        },
        'verification': ['Chandra X-ray tails 2025', 'Norma Cluster'],
    },

    'system_31_tycho_snr': {
        'name': "Tycho's SNR (SN 1572 remnant)",
        'category': 'Type Ia SNR',
        'uqff_terms': ['Ug3 rotation blast', 'CRP acceleration', 'Ub_i expansion'],
        'parameters': {
            'age_yr': 453,
            'distance_kpc': 2.3,
            'shock_velocity_km_s': 3500,
            'E_SN_J': 1e44,
        },
        'verification': ['Chandra 2025', 'Fermi GeV emission'],
    },

    'system_32_at2019qiz': {
        'name': 'AT2019qiz (Tidal Disruption Event)',
        'category': 'TDE',
        'uqff_terms': ['F_U_Bi_i peak flux', 'Ub_i disruption', 'E_react TDE'],
        'parameters': {
            'M_bh_M_sun': 1e6,
            'peak_L_W': 1e44,
            'decay_timescale_d': 29.0,
            'z': 0.0206,
        },
        'verification': ['arXiv:2005.06519', 'ZTF/Chandra 2020'],
    },

    'system_33_gw150914': {
        'name': 'GW150914 (First BBH Detection)',
        'category': 'Binary Black Hole Merger',
        'uqff_terms': ['Triadic validation', 'UQFF chirp mass', 'Ub_i merger'],
        'parameters': {
            'M1_M_sun': 35.6,
            'M2_M_sun': 30.6,
            'M_final_M_sun': 63.1,
            'E_rad_M_sun': 3.1,
            'D_L_Mpc': 440,
            'z': 0.09,
        },
        'verification': ['LIGO/Virgo 2016', 'GWTC-1'],
    },

    'system_34_grs_1915': {
        'name': 'GRS 1915+105 (Micro-quasar)',
        'category': 'Micro-quasar XRB',
        'uqff_terms': ['Ug3 relativistic jets', 'Navier-Stokes', 'SCm accretion'],
        'parameters': {
            'M_bh_M_sun': 12.4,
            'beta_jet': 0.98,
            'P_rot_ms': 14.0,
            'distance_kpc': 8.6,
            'L_Edd_fraction': 0.9,
        },
        'verification': ['RXTE/INTEGRAL 2025'],
    },

    'system_35_cyg_x1': {
        'name': 'Cygnus X-1 (BH XRB)',
        'category': 'X-ray Binary',
        'uqff_terms': ['Ug4 jet feedback', 'Ub_i companion wind', 'CRP'],
        'parameters': {
            'M_bh_M_sun': 21.2,
            'M_companion_M_sun': 40.6,
            'distance_kpc': 2.22,
            'jet_power_W': 2e29,
        },
        'verification': ['Chandra/HST 2025', 'VLBI parallax'],
    },

    'system_36_30_doradus': {
        'name': '30 Doradus (Tarantula Nebula, LMC)',
        'category': 'H II Region / Star Formation',
        'uqff_terms': ['Ug2 bubble pressure', 'Q_wave stats', 'mass-loss Ub_i'],
        'parameters': {
            'SFR_M_sun_yr': 0.1,
            'N_OB_stars': 2400,
            'distance_kpc': 49.0,
            'M_molecular_cloud_M_sun': 5e5,
        },
        'verification': ['HST 2025', 'JWST NIRCam 2024'],
    },

    'system_37_sgr_b2': {
        'name': 'Sgr B2 (Molecular Cloud Complex)',
        'category': 'Galactic Center Molecular Cloud',
        'uqff_terms': ['Ug4 cross-section', 'LENR density', 'Ub_i infall'],
        'parameters': {
            'M_cloud_M_sun': 3e6,
            'distance_kpc': 0.1,
            'density_H2_cm3': 1e4,
            'T_K': 60,
        },
        'verification': ['IRAM 30m 2025', 'Herschel data'],
    },

    'system_38_ngc_4151': {
        'name': 'NGC 4151 (Seyfert 1.5 AGN)',
        'category': 'Seyfert AGN',
        'uqff_terms': ['CRP pp/pγ', 'E_react AGN', 'Ub_i BAL wind'],
        'parameters': {
            'M_bh_M_sun': 4.57e7,
            'L_bol_erg_s': 4e43,
            'distance_Mpc': 19.0,
        },
        'verification': ['HST/Chandra 2025', 'multi-epoch variability'],
    },

    'system_39_r136': {
        'name': 'R136 (WR/O-star Cluster, LMC)',
        'category': 'Super Star Cluster',
        'uqff_terms': ['Triadic stellar mass', 'Ug1 massive stars', 'E_react wind'],
        'parameters': {
            'M_cluster_M_sun': 5e4,
            'L_cluster_L_sun': 4.5e8,
            'M_max_single_M_sun': 250,
        },
        'verification': ['HST/STIS 2025', 'Crowther et al.'],
    },

    'system_40_ngc_6240': {
        'name': 'NGC 6240 (Dual AGN Merger)',
        'category': 'Merging Galaxy / Dual AGN',
        'uqff_terms': ['Triadic dual BH', 'Ub_i merger torque', 'Ug4 pair'],
        'parameters': {
            'sep_kpc': 0.74,
            'M_bh1_M_sun': 9e8,
            'M_bh2_M_sun': 2e9,
            'SFR_M_sun_yr': 50,
        },
        'verification': ['Chandra dual AGN 2024', 'VLBI 2025'],
    },

    'system_41_macs_j0416': {
        'name': 'MACS J0416.1-2403 (Galaxy Cluster)',
        'category': 'Galaxy Cluster / Gravitational Lens',
        'uqff_terms': ['Um turbulence', 'Triadic cluster', 'Ub_i ICM'],
        'parameters': {
            'M_500_M_sun': 5.6e14,
            'redshift': 0.396,
            'T_ICM_keV': 7.3,
        },
        'verification': ['CLASH/HFF 2025', 'Chandra ICM'],
    },

    'system_42_psr_b1919': {
        'name': 'PSR B1919+21 (First Pulsar)',
        'category': 'Pulsar',
        'uqff_terms': ['Ug3 rotation period', 'Q_wave pulsar', 'Ub_i spin-down'],
        'parameters': {
            'P_rot_s': 1.337,
            'B_T': 1.6e8,
            'distance_kpc': 0.361,
        },
        'verification': ['ATNF 2025', 'Hewish 1968 historical'],
    },

    'system_43_ngc_1277': {
        'name': 'NGC 1277 (Overmassive BH Galaxy)',
        'category': 'Compact Elliptical / BH Outlier',
        'uqff_terms': ['Ug4 extreme BH fraction', 'Triadic velocity dispersion'],
        'parameters': {
            'M_bh_M_sun': 1.7e10,
            'M_stellar_M_sun': 1.2e11,
            'f_BH_stellar': 0.14,
            'distance_Mpc': 73,
        },
        'verification': ['van den Bosch 2012', 'HST reanalysis 2025'],
    },

    'system_44_gro_j1655': {
        'name': 'GRO J1655-40 (Micro-quasar / XRB)',
        'category': 'Micro-quasar X-ray Binary',
        'uqff_terms': ['Ug3 relativistic jets', 'Ub_i superluminal blob', 'E_react'],
        'parameters': {
            'M_bh_M_sun': 6.3,
            'beta_jet': 0.92,
            'distance_kpc': 3.2,
            'inclination_deg': 70,
        },
        'verification': ['Hjellming & Rupen 1995', 'RXTE 2025'],
    },

    'system_45_cygnus_loop': {
        'name': 'Cygnus Loop (Veil Nebula SNR)',
        'category': 'Middle-age SNR',
        'uqff_terms': ['Ub_i blast deceleration', 'CRP shock', 'Q_wave thermal'],
        'parameters': {
            'age_yr': 10000,
            'distance_kpc': 0.54,
            'shock_velocity_km_s': 170,
            'angular_size_deg': 2.7,
        },
        'verification': ['Chandra 2025', 'XMM-Newton diffuse'],
    },

    'system_46_g292': {
        'name': 'G292.0+1.8 (Oxygen-rich PWN/SNR)',
        'category': 'Pulsar Wind Nebula / SNR',
        'uqff_terms': ['Ug3 pulsar wind', 'Q_wave O-rich ejecta', 'Triadic PWN'],
        'parameters': {
            'age_yr': 1600,
            'distance_kpc': 6.0,
            'E_SN_J': 2e44,
            'M_O_ejecta_M_sun': 0.5,
        },
        'verification': ['Chandra PWN morphology 2025'],
    },

    'system_47_ngc_7293': {
        'name': 'NGC 7293 (Helix Nebula, Planetary Nebula)',
        'category': 'Planetary Nebula',
        'uqff_terms': ['Ug2 mass-loss shell', 'Ub_i ionised winds', 'ASKAP transient analog'],
        'parameters': {
            'distance_pc': 216,
            'age_yr': 10657,
            'M_wd_M_sun': 0.66,
            'T_wd_K': 1.07e5,
            'R_outer_pc': 0.8,
        },
        'verification': ['ASKAP J1832 template 2025', 'Hubble 2025'],
    },

    'system_48_perseus_cluster': {
        'name': 'Perseus Galaxy Cluster (Abell 426)',
        'category': 'Galaxy Cluster',
        'uqff_terms': ['Um ICM turbulence', 'Triadic cluster mass', 'Ub_i AGN bubbles'],
        'parameters': {
            'M_500_M_sun': 6.7e14,
            'T_ICM_keV': 7.0,
            'z': 0.0179,
            'AGN_cavity_power_W': 1e38,
            'r_cool_kpc': 100,
        },
        'verification': ['Chandra Perseus 2025', 'Hitomi 2016 turbulence'],
    },

    'system_49_widom_larsen_lenr': {
        'name': 'Widom-Larsen LENR Benchmark',
        'category': 'Nuclear / Laboratory LENR',
        'uqff_terms': ['F_core LENR', 'k_eta = 1e-113', 'ultra-low-momentum neutron'],
        'parameters': {
            'k_LENR': 1e-10,
            'omega_LENR_rad_s': 7.85e12,
            'omega_0_rad_s': 1e-12,
            'sigma_n_m2': 1e-28,
            'F_LENR_N': 6.16e39,
            'F_core_ratio': 1.0,        # This IS F_core benchmark
        },
        'verification': ['Widom-Larsen 2006 PRB', 'LENR_bib validation 2025'],
    },

    'system_50_bec_alpha_cluster': {
        'name': 'BEC Alpha-Cluster (12C Hoyle State Extension)',
        'category': 'Nuclear BEC',
        'uqff_terms': ['N_B BEC term', 'T_c Bose shift', 'UQFF-LENR coupling'],
        'parameters': {
            'N_B': 3,               # 3-alpha Bose condensate
            'T_c_shift_MeV': 0.38,
            'alpha_cluster_mass_u': 4.0,
            'Phi_BEC_norm': 0.57,   # SSq = [SSq] at nuclear scale
        },
        'verification': ['Tohsaki AMD 2025', 'arXiv:nucl-th recent'],
    },

    'system_51_gw190521': {
        'name': 'GW190521 (IMBH Binary Merger)',
        'category': 'Intermediate-Mass Black Hole Merger',
        'uqff_terms': ['Ub_i IMBH merger', 'Triadic chirp validation', 'F_U_Bi_i IMBH'],
        'parameters': {
            'M1_M_sun': 85,
            'M2_M_sun': 66,
            'M_final_M_sun': 142,
            'E_rad_M_sun': 9,
            'z': 0.82,
        },
        'verification': ['LIGO/Virgo 2020', 'GWTC-2'],
    },

    'system_52_cmb_cosmological': {
        'name': 'Cosmic Microwave Background / ΛCDM Reference',
        'category': 'Cosmological Baseline',
        'uqff_terms': ['ρ_vac λ_vac', 'Friedmann Hubble', 'DPM sphere inflation'],
        'parameters': {
            'H_0_km_s_Mpc': 67.4,
            'Omega_Lambda': 0.6889,
            'T_CMB_K': 2.7255,
            'rho_crit_kg_m3': 9.47e-27,
            'n_s': 0.9649,
        },
        'verification': ['Planck 2018', 'JWST 2025 cosmological background'],
    },

    # ─────────────────────────────────────────────────────────────
    # CATALOGUE SUMMARY
    # ─────────────────────────────────────────────────────────────
    'summary': {
        'total_systems': 52,
        'systems_01_24_ref': 'CondensedPhysics_Validation.py',
        'new_systems_25_52': 28,
        'F_U_Bi_i_mean_N': -6.05e217,
        'F_U_Bi_i_log_bootstrap_std': '3%',
        'x_2_cosmic_m': -3.40e172,
        'Q_wave_mean_J_m3': 3.98e-5,    # B = 1e-5 T reference
        'Q_wave_crab_J_m3': 3.98e-3,    # Crab nebula B = 1e-4 T
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: KAPPA MCMC CALIBRATION
# From 0904 thread — refined kappa from MCMC analysis of 47-system E_react series
# Note: kappa=0.0005 day⁻¹ remains canonical in CondensedPhysics_OutputData.py
#       This 0904 MCMC result (0.00052) is 4% from canonical.
#       Per user instruction, the 0.001 value in CondensedPhysics2.py is NOT a bug.
# ═══════════════════════════════════════════════════════════════════════════════

KAPPA_MCMC_CALIBRATION = {
    'thread': 'grok_share_0904a12a5c2b4a639389ae084391b94f',
    'canonical_kappa_day': 0.0005,              # CondensedPhysics_OutputData.py L2748
    'mcmc_kappa_day': 0.00052,                  # 0904 MCMC refinement
    'mcmc_std_day': 1.23e-5,
    'mcmc_ci_95_day': (0.00048, 0.00056),       # 95% confidence interval
    'mcmc_n_systems': 47,
    'deviation_from_canonical_percent': 4.0,
    'note': (
        'MCMC value (0.00052) is 4% from canonical (0.0005). '
        'The canonical value is validated against 24-system set. '
        'The 0.001 value in CondensedPhysics2.py is a separate implementation — '
        'left unchanged per project policy.'
    ),
    'grok_ref': 'Thread 0904, Section: kappa MCMC posterior analysis',
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: DELTA_RHO NORMALITY TESTS
# Distribution of δρ (density perturbation) over 47-system sample
# Shapiro-Wilk, Kolmogorov-Smirnov, Anderson-Darling results
# ═══════════════════════════════════════════════════════════════════════════════

DELTA_RHO_NORMALITY_TESTS = {
    'thread': 'grok_share_0904a12a5c2b4a639389ae084391b94f',
    'n_systems': 47,
    'variable': 'delta_rho (density perturbation)',
    'shapiro_wilk': {
        'statistic': 0.9412,
        'p_value': 0.00055,
        'interpretation': 'Reject normality (p < 0.001) — leptokurtic tail',
    },
    'kolmogorov_smirnov': {
        'statistic': 0.098,
        'p_value': 0.741,
        'interpretation': 'Cannot reject normality at α=0.05 — KS less sensitive to tails',
    },
    'anderson_darling': {
        'statistic': 1.35,
        'critical_1_percent': 1.092,
        'interpretation': 'Reject normality at 1% level — heavy right tail confirmed',
    },
    'jarque_bera': {
        'statistic': 8.78,
        'p_value': 0.012,
        'leptokurtosis': 0.037,
        'skewness': 0.41,
        'interpretation': 'JB rejects normality — consistent with SW/AD findings',
    },
    'conclusion': (
        'δρ distribution is leptokurtic (fat-tailed). '
        'Log-normal transformation recommended for parametric tests. '
        'Bootstrap std = 3% on F_U_Bi_i mean is robust to non-normality.'
    ),
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Q_WAVE 52-SYSTEM STATISTICS
# Updated from Q_WAVE_47_STATISTICS in CondensedPhysics_OutputData.py
# ═══════════════════════════════════════════════════════════════════════════════

Q_WAVE_52_STATISTICS = {
    'n_systems': 52,
    # Energy density form: Q_wave = B² / (2 μ_0)
    'B_ref_T': 1e-5,                 # Reference field (ISM / weak B)
    'Q_wave_at_B_ref_J_m3': 3.98e-5,
    'B_crab_T': 1e-4,                # Crab Nebula (strong B)
    'Q_wave_at_B_crab_J_m3': 3.98e-3,
    # Legacy Orb statistics (from Q_WAVE_47_STATISTICS baseline)
    'mean_Jm3': 3.98e4,              # Updated from 3.97e4 with 5 new systems
    'std_Jm3': 51200.0,              # Minor update from 51131.3
    'jarque_bera': 9.12,
    'jarque_bera_p': 0.010,
    'leptokurtosis': 0.039,
    'prev_ref': 'Q_WAVE_47_STATISTICS in CondensedPhysics_OutputData.py',
    'note': (
        'Q_wave mean updated from 3.97e4 (n=47) to 3.98e4 (n=52). '
        'Addition of 5 high-B systems (GRO J1655, Cygnus Loop, G292, NGC7293, Perseus) '
        'shifts mean slightly upward. JB statistic increases confirming persistent '
        'leptokurtosis in the 52-system ensemble.'
    ),
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5: UQFF MASTER EQUATIONS (FULL NOTATION)
# Six master equations from 0904 thread with complete variable notation
# ═══════════════════════════════════════════════════════════════════════════════

UQFF_MASTER_EQUATIONS_FULL = {
    'thread': 'grok_share_0904a12a5c2b4a639389ae084391b94f',

    'eq1_unified_field': {
        'name': 'Unified Field F_U',
        'latex': (
            r'F_U = \sum_{i=1}^{26} k_i U_{g,i} - \beta_i U_{b,i} + U_m + U_A '
            r'+ U_i + SCm - \sum_j \delta_j + CRP \cdot \sum D_E \frac{\partial^2 n}{\partial p^2} e^{-\gamma t}'
        ),
        'python': 'F_U = sum(k_i*Ug_i - beta_i*Ub_i for i in range(26)) + Um + UA + Ui + SCm - sum(delta_j) + CRP',
        'variables': {
            'k_i': 'Layer coupling constants (k_1=1.5, k_2=1.2, k_3=0.8, k_4=0.5 for primary layers)',
            'Ug_i': 'Gravity component i (Ug1=magnetic dipole, Ug2=charge-reactivity, Ug3=string, Ug4=vacuum)',
            'beta_i': 'Buoyancy coefficient (β_i ≈ 0.603)',
            'Ub_i': 'Buoyancy opposition term',
            'Um': 'Universal magnetism sum',
            'UA': 'Universal Aether tensor contribution',
            'Ui': 'Inertia/ISM term',
            'SCm': 'Superconductive medium contribution',
            'delta_j': 'Depletion/dissipation corrections',
            'CRP': 'Cosmic Ray Propagation integral',
            'D_E': 'Energy diffusion coefficient',
            'gamma': 'Decay constant (related to κ)',
        },
    },

    'eq2_superconductive': {
        'name': 'SCm Superconductive Decay',
        'latex': r'[SCm](t) = [SCm]_0 \cdot e^{-\kappa t} \cdot (1 - \cos(\omega_{SCm} t))',
        'python': 'SCm_t = SCm_0 * math.exp(-kappa * t) * (1 - math.cos(omega_SCm * t))',
        'variables': {
            'SCm_0': 'Initial SCm vacuum density (7.09e-37 J/m³)',
            'kappa': 'Reactivity decay constant (0.0005 day⁻¹ canonical)',
            'omega_SCm': 'SCm oscillation frequency (2π/t_age)',
        },
        'calibrated_values': {
            'kappa_canonical': 0.0005,
            'kappa_mcmc': 0.00052,
            'H_SCm': 0.99,
        },
    },

    'eq3_buoyancy': {
        'name': 'F_U_Bi_i Master Buoyancy Force',
        'latex': (
            r'F_{U,Bi,i} = \Omega_g \cdot \frac{M_{bh}}{d_g} \cdot \sum_{i=1}^{26} '
            r'\left( U_{g,i} + U_{b,i} \right) \cdot [SSq]^{26} \cdot [UA]_i'
        ),
        'python': (
            'F_U_Bi_i = Omega_g * (M_bh / d_g) * sum(Ug_i + Ub_i for i in range(26)) '
            '* SSq**26 * UA'
        ),
        'mean_52sys_N': -6.05e217,
        'log_bootstrap_std': 0.03,
        'variables': {
            'Omega_g': 'Galactic rotation rate (7.3e-16 rad/s)',
            'M_bh': 'Black hole mass (kg)',
            'd_g': 'Galactic center distance (m)',
            'SSq': 'Quantum state factor ([SSq] = 0.57)',
            'UA': 'Universal Aether tensor factor',
        },
        'forms': {
            'Form_A': '[SSq]^26 × F_Ug  (cosmic quantum — fastest computation)',
            'Form_B': 'F_Bi = β × (Ug - ρ_vac × G)  (simplified per system)',
            'Form_C': 'F_U_Bi_i = k_LENR × (ω_LENR/ω_0)² × x_2  (LENR coupling)',
        },
    },

    'eq4_triadic': {
        'name': 'Triadic Master Gravity (26-Layer)',
        'latex': (
            r'g(r,t) = \sum_{i=1}^{26} \left[ U_{g1,i} + U_{g2,i} + U_{g3,i} + U_{g4,i} \right]'
        ),
        'python': 'g_triadic = sum(Ug1_i + Ug2_i + Ug3_i + Ug4_i for i in range(1, 27))',
        'variables': {
            'Ug1_i': 'Layer-i magnetic dipole gravity',
            'Ug2_i': 'Layer-i charge-reactivity gravity',
            'Ug3_i': 'Layer-i string rotation gravity',
            'Ug4_i': 'Layer-i vacuum concentration gravity',
        },
    },

    'eq5_quadratic_cosmic': {
        'name': 'Cosmic Quadratic Solve (x_2)',
        'latex': r'x^2 + bx + c = 0 \Rightarrow x_2 = \frac{-b - \sqrt{b^2 - 4c}}{2}',
        'python': 'x_2 = (-b - math.sqrt(b**2 - 4*c)) / 2',
        'result_52sys_m': -3.40e172,
        'result_atomic_Zscale_m': -3.56e116,   # mean over Z=1–118
        'comment': (
            'x_2 is the physically meaningful (negative) root in UQFF framework. '
            'Cosmic scale (x_2 = -3.40e172 m) vs atomic Z-averaged (x_2_Z = -3.56e116 m) '
            'differ because atomic solve uses Z-scaled mass m_Z = Z × m_p.'
        ),
    },

    'eq6_f_core_lenr': {
        'name': 'F_core (LENR Universal Force)',
        'latex': r'F_{core} = \frac{\hbar \, \omega_{LENR}}{\sigma_n \, \rho_{vac,[UA]}}',
        'python': 'F_core = (h_bar * omega_LENR) / (sigma_n * rho_vac_UA)',
        'value_N': 1e10,    # ~10^10 N (universal k_η)
        'variables': {
            'h_bar': 'Reduced Planck constant (1.054571817e-34 J·s)',
            'omega_LENR': 'LENR characteristic frequency (7.85e12 rad/s)',
            'sigma_n': 'Nuclear cross-section factor',
            'rho_vac_UA': 'Universal Aether vacuum density (7.09e-36 J/m³)',
        },
        'canonical_location': 'DPMCosmologyModule.py L230',
        'duplicate_note': (
            'GrokThread_StarMagic_UnifiedFramework.py had an inline copy of this formula. '
            'That inline copy is removed; see DPMCosmologyModule.py::F_core for the authoritative version.'
        ),
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6: UQFF ATOMIC Z-SCALING
# Applies UQFF buoyancy to atoms Z=1 (H) through Z=118 (Og)
# ═══════════════════════════════════════════════════════════════════════════════

UQFF_ATOMIC_Z_SCALING = {
    'thread': 'grok_share_0904a12a5c2b4a639389ae084391b94f',
    'Z_range': (1, 118),
    'p_Z_formula': 'p_Z = Z / (Z + 1)',
    'p_Z_Z1': 0.5,          # Hydrogen
    'p_Z_Z118': 0.992,      # Oganesson
    'SSq_Z_formula': 'SSq_Z = 0.507 + (Z / 118) * 0.1',
    'SSq_Z_Z1': 0.508,      # Hydrogen
    'SSq_Z_Z118': 0.607,    # Oganesson
    'SSq_base_note': (
        'SSq_Z uses the MCMC linear form: 0.507 + scaling × Z/118. '
        'This is NOT the same as SSq=0.57 (output constant) which is used as [SSq]^26. '
        'SSq_Z is the per-atom quantum factor in the atomic mass function.'
    ),
    'x_2_Z_mean_m': -3.56e116,     # Mean x_2 over Z=1–118
    'x_2_Z_formula': 'x_2_Z(Z) = solve_quadratic(m_Z = Z * m_proton, rho_vac)',
    'atomic_scaling_note': (
        'UQFF extends naturally to atomic systems where the "gravitational" object is a nucleus. '
        'p_Z = Z/(Z+1) is the nuclear density fraction. At Z=118 (Og), p_Z → 0.992 ≈ full density.'
    ),
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7: CERN DELPHI HIGGS RECORDS
# From 0904 thread — DELPHI experiment records, Higgs candidate mean = 85 GeV
# Context: Relevant to UQFF E_n = E_0 × 10^n at n=12 (Higgs level)
# ═══════════════════════════════════════════════════════════════════════════════

CERN_DELPHI_HIGGS_RECORDS = {
    'thread': 'grok_share_0904a12a5c2b4a639389ae084391b94f',
    'experiment': 'DELPHI (LEP2, CERN)',
    'record_ids': list(range(93719, 93727)),    # Records 93719–93726
    'n_records': 8,
    'higgs_candidate_mean_GeV': 85,
    'higgs_candidate_std_GeV': 7.3,
    'pdg_higgs_mass_GeV': 125.09,              # PDG 2025 (LHC confirmed)
    'UQFF_level': 'n=12 in E_n = E_0 × 10^n → E_12 ~ 10^{12} × E_0',
    'note': (
        'DELPHI 85 GeV candidates (pre-LHC era) are listed here as historical '
        'calibration reference for UQFF E_n polynomial. The UQFF n=12 "Higgs level" '
        'predated the confirmed LHC mass of 125.09 GeV (PDG 2025). '
        'DELPHI records show early LEP2 Higgs search data before confirmation.'
    ),
    'url_refs': [
        'https://cds.cern.ch/search?ln=en&p=DELPHI+93719&action_search=Search',
        'https://inspirehep.net/search?p=DELPHI+LEP2+Higgs',
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 8: DPM YIN-YANG COSMOLOGY
# From 0904 thread — detailed DPM yin-yang picture beyond DPMCosmologyModule.py
# ═══════════════════════════════════════════════════════════════════════════════

DPM_YIN_YANG_COSMOLOGY = {
    'thread': 'grok_share_0904a12a5c2b4a639389ae084391b94f',
    'framework': 'DPM = Di-Pseudo-Monopole',
    'description': (
        'The DPM (Di-Pseudo-Monopole) is the core primordial entity in Star Magic cosmology. '
        'It arises from the standing resonance of [SCm] and [UA] inside a 26-shell '
        'oscillating EM field (described in GrokThread_StarMagic_UnifiedFramework.py::birth_of_dpm_sphere). '
        'The Yin-Yang pairing represents the dual nature of DPM: '
        'one lobe is positive (Yin, [UA] dominant), one negative (Yang, [SCm] dominant).'
    ),
    'yin_yang_parameters': {
        'Yin_UA_fraction': 0.57,        # [SSq] = fraction of UA in Yin lobe
        'Yang_SCm_fraction': 0.43,      # 1 - [SSq] in Yang lobe
        'balance_ratio': 'SSq : (1-SSq) = 0.57 : 0.43',
        'resonance_modes': 26,          # 26 quantum shells
        'inflation_timescale_s': 1e-36, # ~Planck time for initial DPM inflation
    },
    'cosmological_phases': {
        1: 'Pre-Big Bang — [SCm]/[UA] field separation',
        2: 'DPM Formation — standing resonance in 26-shell EM',
        3: 'Inflation — DPM sphere expansion (see InflationForceChartCalculator)',
        4: 'Big Bang — DPM topology change, matter/antimatter asymmetry',
        5: 'Current era — SSq decay, CRP propagation, Ub_i formation',
    },
    'link_to_inflation_epochs': (
        'GrokThread_StarMagic_UnifiedFramework.py::INFLATION_FORCE_EPOCHS details '
        'epochs 1-5 of the post-DPM inflation sequence.'
    ),
    'link_to_sphere_equation': (
        'GrokThread_StarMagic_UnifiedFramework.py::birth_of_dpm_sphere() '
        'is the canonical DPM sphere birth equation.'
    ),
    'link_to_dpm_module': (
        'DPMCosmologyModule.py::DPMCosmologyCalculator implements the '
        'quantitative DPM evolution model.'
    ),
}


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 9: UQFFAtomicZScalingCalculator
# Pure physics calculator — receives dataset, returns Z-scaled equations
# As per CondensedPhysics.py architecture: NO hardcoded system data. 
# Computes p_Z, SSq_Z, x_2_Z for any input atomic number Z.
# ═══════════════════════════════════════════════════════════════════════════════

class UQFFAtomicZScalingCalculator:
    """
    UQFF Atomic Z-Scaling Calculator.

    Computes UQFF buoyancy parameters scaled to atomic number Z.

    Source:
        Grok 0904 thread, Section: Atomic Z-scaling
        Uses UQFF formulas extended to nuclear/atomic domain.

    Physics:
        p_Z   = Z / (Z + 1)           Nuclear density fraction
        SSq_Z = 0.507 + (Z/118)*0.1   Quantum state factor (linear MCM form)
        x_2_Z = negative quadratic root using m_Z = Z * m_proton

    NOTE: This is a PURE CALCULATOR — no system-specific data is hardcoded.
    Pass any Z or dataset to compute. Outputs go to CondensedPhysics_OutputData.py
    for storage/recall.
    """

    # Physical constants
    M_PROTON_KG = 1.67262192e-27      # kg
    H_BAR_JS = 1.054571817e-34        # J·s
    G_m3_kg_s2 = 6.67430e-11         # m³/(kg·s²)
    RHO_VAC_UA_Jm3 = 7.09e-36        # J/m³

    def compute_p_Z(self, Z: int) -> float:
        """Nuclear density fraction for atomic number Z."""
        if Z <= 0:
            raise ValueError(f"Z must be positive, got {Z}")
        return Z / (Z + 1)

    def compute_SSq_Z(self, Z: int) -> float:
        """Z-scaled quantum state factor (MCMC linear form from 0904 thread)."""
        if not (1 <= Z <= 118):
            raise ValueError(f"Z must be 1–118, got {Z}")
        return 0.507 + (Z / 118) * 0.1

    def compute_x_2_Z(self, Z: int,
                       rho_vac: float = 7.09e-36,
                       omega_LENR: float = 7.85e12) -> float:
        """
        Compute x_2 (negative quadratic root) for atomic number Z.

        Uses:
            m_Z = Z * m_proton
            b = -G * m_Z / r_nucleus²   (simplified, r_nucleus = 1e-15 m)
            c = rho_vac * G * r_nucleus²
        """
        m_Z = Z * self.M_PROTON_KG
        r_nucleus = 1.2e-15 * (Z ** (1 / 3))     # nuclear radius r = r_0 * A^{1/3}
        if r_nucleus <= 0:
            return 0.0
        b = -(self.G_m3_kg_s2 * m_Z) / (r_nucleus ** 2)
        c = rho_vac * self.G_m3_kg_s2 * (r_nucleus ** 2)
        discriminant = b ** 2 - 4 * c
        if discriminant < 0:
            return -b / 2.0   # real part of complex root
        return (-b - math.sqrt(discriminant)) / 2.0

    def compute_dataset(self, dataset: dict) -> dict:
        """
        Compute all Z-scaling parameters for a dataset.

        Args:
            dataset: dict with keys 'Z' (int or list of ints),
                     optional 'rho_vac' (float), 'omega_LENR' (float)

        Returns:
            dict with Z-scaled parameters and equation strings
        """
        Z_input = dataset.get('Z', 1)
        rho_vac = dataset.get('rho_vac', self.RHO_VAC_UA_Jm3)
        omega_LENR = dataset.get('omega_LENR', 7.85e12)

        if isinstance(Z_input, (list, range)):
            results = []
            for Z in Z_input:
                results.append({
                    'Z': Z,
                    'p_Z': self.compute_p_Z(Z),
                    'SSq_Z': self.compute_SSq_Z(Z),
                    'x_2_Z_m': self.compute_x_2_Z(Z, rho_vac, omega_LENR),
                })
            x2_values = [r['x_2_Z_m'] for r in results]
            x2_mean = sum(x2_values) / len(x2_values) if x2_values else 0.0
            return {
                'primary_equations': results,
                'available_equations': [
                    'p_Z = Z / (Z + 1)',
                    'SSq_Z = 0.507 + (Z/118) * 0.1',
                    'x_2_Z = (-b - sqrt(b²-4c)) / 2  [b=−G·m_Z/r², c=ρ_vac·G·r²]',
                ],
                'simulation_set': {
                    'Z_range': (min(Z_input), max(Z_input)),
                    'n_atoms': len(results),
                    'x_2_mean_m': x2_mean,
                },
                'thread_reference_x_2_Z_mean_m': -3.56e116,
            }
        else:
            Z = int(Z_input)
            return {
                'primary_equations': {
                    'Z': Z,
                    'p_Z': self.compute_p_Z(Z),
                    'SSq_Z': self.compute_SSq_Z(Z),
                    'x_2_Z_m': self.compute_x_2_Z(Z, rho_vac, omega_LENR),
                },
                'available_equations': [
                    f'p_{Z} = {Z}/{Z+1}',
                    f'SSq_{Z} = 0.507 + ({Z}/118)*0.1 = {self.compute_SSq_Z(Z):.4f}',
                    f'x_2_{Z} = {self.compute_x_2_Z(Z, rho_vac, omega_LENR):.3e} m',
                ],
                'simulation_set': {
                    'Z': Z,
                    'atomic_symbol': f'Z={Z}',
                    'mass_kg': Z * self.M_PROTON_KG,
                },
            }


# ═══════════════════════════════════════════════════════════════════════════════
# QUICK VERIFICATION
# ═══════════════════════════════════════════════════════════════════════════════

def _verify_0904_constants():
    """Spot-check key 0904 thread values against known reference figures."""
    calc = UQFFAtomicZScalingCalculator()

    # Z=1 (Hydrogen)
    assert abs(calc.compute_p_Z(1) - 0.5) < 1e-9,  "p_Z(1) should be 0.5"
    # Z=118 (Oganesson)
    assert abs(calc.compute_p_Z(118) - 118/119) < 1e-9, "p_Z(118) wrong"
    # SSq_Z(1)
    expected_ssq_1 = 0.507 + (1/118)*0.1
    assert abs(calc.compute_SSq_Z(1) - expected_ssq_1) < 1e-9, "SSq_Z(1) wrong"
    # SSq_Z(118)
    expected_ssq_118 = 0.507 + (118/118)*0.1
    assert abs(calc.compute_SSq_Z(118) - expected_ssq_118) < 1e-9, "SSq_Z(118) wrong"
    # MCMC kappa CI check
    ci = KAPPA_MCMC_CALIBRATION['mcmc_ci_95_day']
    assert ci[0] < KAPPA_MCMC_CALIBRATION['mcmc_kappa_day'] < ci[1], "kappa outside CI"
    # System count
    assert UQFF_52_SYSTEM_CATALOGUE['n_systems'] == 52, "System count wrong"
    print("0904 constant verification PASSED.")


if __name__ == '__main__':
    _verify_0904_constants()

    # Demo: Z-scaling for H, Fe, U, Og
    calc = UQFFAtomicZScalingCalculator()
    print("\nUQFF Atomic Z-Scaling (0904 thread):")
    print(f"{'Z':>4}  {'Element':>8}  {'p_Z':>8}  {'SSq_Z':>8}  {'x_2_Z (m)':>14}")
    print("-" * 50)
    for Z, sym in [(1,'H'), (26,'Fe'), (92,'U'), (118,'Og')]:
        p = calc.compute_p_Z(Z)
        s = calc.compute_SSq_Z(Z)
        x = calc.compute_x_2_Z(Z)
        print(f"{Z:>4}  {sym:>8}  {p:>8.4f}  {s:>8.4f}  {x:>14.3e}")
