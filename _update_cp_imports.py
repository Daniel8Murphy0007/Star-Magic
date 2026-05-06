#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
_update_cp_imports.py  -- one-time import consolidation for CP1-CP4

Replaces all try/except import blocks for scm_vacuum_manifold and
ua_vacuum_manifold with direct imports from the consolidated
dpm_vacuum_manifold.  The DPM block is also cleaned up.

Run once, then delete this script.
"""

import os, shutil

CWD = os.path.dirname(os.path.abspath(__file__))

def read_text(path):
    with open(path, 'r', encoding='utf-8') as f:
        return f.read()

def write_text(path, content):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(content)


# ===========================================================================
# COMMON BLOCKS (identical in CP1, CP2, CP3, CP4)
# ===========================================================================

TOP_OLD = '''\
try:
    from scm_vacuum_manifold import derive_from_quantum_chain as _derive_qc
    _RHO_VAC_SCM, _ = _derive_qc(n_levels=26, f_SCm=0.57)   # J/m\xb3 SCm energy density
    _RHO_VAC_UA,  _ = _derive_qc(n_levels=26, f_SCm=5.7)    # J/m\xb3 UA  energy density (10x)
except ImportError:
    # Fallback: canonical numeric values if scm_vacuum_manifold not on path
    _RHO_VAC_SCM = 6.333333e+05   # J/m\xb3 \u2014 SCm vacuum energy density (Quantum Chain)
    _RHO_VAC_UA  = 6.333333e+06   # J/m\xb3 \u2014 UA  vacuum energy density (Quantum Chain)'''

TOP_NEW = '''\
from dpm_vacuum_manifold import derive_from_quantum_chain as _derive_qc
_RHO_VAC_SCM, _ = _derive_qc(n_levels=26, f_SCm=0.57)   # J/m\xb3 SCm energy density
_RHO_VAC_UA,  _ = _derive_qc(n_levels=26, f_SCm=5.7)    # J/m\xb3 UA  energy density (10x)'''

# ---------------------------------------------------------------------------
UA_OLD = '''\
# ==================== UA VACUUM MANIFOLD IMPORTS [ua_vacuum_manifold.py] ====================
try:
    from ua_vacuum_manifold import (
        ua_layer_density             as _ua_layer_density,
        ua_dpm_total_density         as _ua_dpm_total_density,
        ua_dpm_buoyancy_factor       as _ua_dpm_buoyancy_factor,
        ua_calibration_ratio         as _ua_calibration_ratio,
        DPM_DENSITY_RATIO            as _UA_DPM_DENSITY_RATIO,
        E_PHONON                     as _UA_E_PHONON,
        S26_3                        as _UA_S26_3,
        PHI_RES                      as _UA_PHI_RES,
        DELTA_UA_FOURTH              as _UA_DELTA_UA_FOURTH,
        F_U_Bi_i_DPM                 as _ua_F_U_Bi_i_DPM,
        ua_lenr_comparison           as _ua_lenr_comparison,
        ua_casimir_comparison        as _ua_casimir_comparison,
        ua_cosmological_acceleration as _ua_cosmo_accel,
        ua_rotation_curve_flat       as _ua_rotation_curve,
        ua_hubble_tension_modulation as _ua_hubble_tension,
        ua_dark_energy_substitute    as _ua_dark_energy,
    )
    _UA_MANIFOLD_LOADED = True
except ImportError:
    _UA_MANIFOLD_LOADED = False
    def _ua_layer_density(layer=1, t_n_val=0.0): return _RHO_VAC_SCM
    def _ua_dpm_total_density(t_n_val=0.0): return 4 * _RHO_VAC_SCM
    def _ua_dpm_buoyancy_factor(t_n_val=0.0): return 4.0
    def _ua_calibration_ratio(): return 10.0
    _UA_E_PHONON = 6.62607015e-34 * 1.25e12
    _UA_S26_3 = 1.4531e26
    _UA_PHI_RES = 0.84
    _UA_DELTA_UA_FOURTH = 0.0
    _UA_DPM_DENSITY_RATIO = 10.0
    def _ua_F_U_Bi_i_DPM(t_n_val=0.0, **kw): return 0.0
    def _ua_lenr_comparison(): return {}
    def _ua_casimir_comparison(): return {}
    def _ua_cosmo_accel(z=0.0): return 0.0
    def _ua_rotation_curve(r=1e20, v0=220e3): return v0
    def _ua_hubble_tension(t=0.0): return 0.0
    def _ua_dark_energy(t_n_val=0.5): return 0.0'''

UA_NEW = '''\
# ==================== UA VACUUM MANIFOLD IMPORTS (from dpm_vacuum_manifold) ====================
from dpm_vacuum_manifold import (
    ua_layer_density             as _ua_layer_density,
    ua_dpm_total_density         as _ua_dpm_total_density,
    ua_dpm_buoyancy_factor       as _ua_dpm_buoyancy_factor,
    ua_calibration_ratio         as _ua_calibration_ratio,
    DPM_DENSITY_RATIO            as _UA_DPM_DENSITY_RATIO,
    E_PHONON                     as _UA_E_PHONON,
    S26_3                        as _UA_S26_3,
    PHI_RES                      as _UA_PHI_RES,
    DELTA_UA_FOURTH              as _UA_DELTA_UA_FOURTH,
    F_U_Bi_i_DPM                 as _ua_F_U_Bi_i_DPM,
    ua_lenr_comparison           as _ua_lenr_comparison,
    ua_casimir_comparison        as _ua_casimir_comparison,
    ua_cosmological_acceleration as _ua_cosmo_accel,
    ua_rotation_curve_flat       as _ua_rotation_curve,
    ua_hubble_tension_modulation as _ua_hubble_tension,
    ua_dark_energy_substitute    as _ua_dark_energy,
)
_UA_MANIFOLD_LOADED = True'''

# ---------------------------------------------------------------------------
DPM_OLD = '''\
# ==================== DPM VACUUM MANIFOLD IMPORTS [dpm_vacuum_manifold.py] ====================
try:
    from dpm_vacuum_manifold import (
        PERIODIC_TABLE    as _DPM_PERIODIC_TABLE,
        ELEMENT           as _DPM_ELEMENT,
        E_CRACK           as _DPM_E_CRACK,
        M_0_DPM           as _DPM_M0,
        DPMBody           as _DPMBody,
        DPM_DENSITY_RATIO as _DPM_DENSITY_RATIO,
    )
    _DPM_MANIFOLD_LOADED = True
except ImportError:
    _DPM_MANIFOLD_LOADED = False
    _DPM_PERIODIC_TABLE = []
    _DPM_ELEMENT = {}
    _DPM_E_CRACK = float(_RHO_VAC_SCM * (2.99792458e8) ** 2) / 0.57
    _DPM_M0 = _DPM_E_CRACK / (2.99792458e8) ** 2
    _DPMBody = None
    _DPM_DENSITY_RATIO = 10.0'''

DPM_NEW = '''\
# ==================== DPM VACUUM MANIFOLD IMPORTS ====================
from dpm_vacuum_manifold import (
    PERIODIC_TABLE    as _DPM_PERIODIC_TABLE,
    ELEMENT           as _DPM_ELEMENT,
    E_CRACK           as _DPM_E_CRACK,
    M_0_DPM           as _DPM_M0,
    DPMBody           as _DPMBody,
    DPM_DENSITY_RATIO as _DPM_DENSITY_RATIO,
)
_DPM_MANIFOLD_LOADED = True'''

# ===========================================================================
# SCm blocks — CP1 / CP2 (same text, private aliases with leading _)
# ===========================================================================
SCM_CP12_OLD = '''\
# --- import updated module-level constants from canonical [pdf/scm_vacuum_manifold.py] ---
try:
    from scm_vacuum_manifold import (
        E_phonon                       as _E_PHONON_SCM,   # overrides hardcoded value if canonical available
        S26_3                          as _S26_3,
        Phi_resonance                  as _PHI_RESONANCE,
        KER_SCm                        as _KER_SCM,
        scaling_factor                 as _SCALING_SCM,    # exact 630 eV normalizer
        KAPPA_FLOAT                    as _KAPPA_FLOAT_SCM, # float(KAPPA) = 0.0005
        F_TRZ                          as _F_TRZ_SCM,      # canonical 0.1
        coleman_guillespie_scm         as _scm_coleman_guillespie,  # beta->1.25THz phonon
        neutrino_oscillation_prob_lenr as _scm_neutrino_osc,        # P_osc via S26_3*cos(pi*t_n)
        quark_production_prob_ui       as _scm_quark_prod,          # |Phi|^2*cos(pi*t_n)*Ui
        mckubre_lenr                   as _scm_mckubre,             # Pd-D sub-barrier via F_U_Bi_i
        s26_3_from_vds                 as _scm_s26_3_from_vds,      # canonical 1.4531e26
        qgp_energy_density_scm         as _scm_qgp_energy_density,  # QGP tokamak energy [J]
        strange_quark_matter_density   as _scm_sqm_density,         # SQM (density, B_eff)
        mit_bag_scm                    as _scm_mit_bag,             # MIT bag B_eff [J/m^3]
        ads_cft_scm_dual               as _scm_ads_cft_dual,        # SCm holographic mapping dict
        scm_gw_metric_perturbation     as _scm_gw_metric_pert,      # h(f_gw, r) GW strain
    )
except ImportError:
    pass  # fallback values already set above
    _F_TRZ_SCM = 0.1
    def _scm_coleman_guillespie(decay_rate=1.0e6, t_n=-100.0, Gamma=1.0e12): return 0.0
    def _scm_neutrino_osc(t_n=-100.0): return 0.0
    def _scm_quark_prod(t_n=-100.0, Gamma=1.0e12): return 0.0
    def _scm_mckubre(PdD_loading=0.9, volume=1.0e-6, t_n=-100.0): return 0.0
    _scm_s26_3_from_vds = lambda: 1.4531e26
    def _scm_qgp_energy_density(T_plasma=1.0e11): return 0.0
    def _scm_sqm_density(): return (1.0e18, 0.0)
    def _scm_mit_bag(): return 0.0
    def _scm_ads_cft_dual(): return {}
    def _scm_gw_metric_pert(f_gw=100.0, r_detector=3.086e22): return 0.0'''

SCM_CP12_NEW = '''\
# --- SCm constants from dpm_vacuum_manifold (consolidated) ---
from dpm_vacuum_manifold import (
    E_phonon                       as _E_PHONON_SCM,
    S26_3                          as _S26_3,
    Phi_resonance                  as _PHI_RESONANCE,
    KER_SCm                        as _KER_SCM,
    scaling_factor                 as _SCALING_SCM,
    KAPPA_FLOAT                    as _KAPPA_FLOAT_SCM,
    F_TRZ                          as _F_TRZ_SCM,
    coleman_guillespie_scm         as _scm_coleman_guillespie,
    neutrino_oscillation_prob_lenr as _scm_neutrino_osc,
    quark_production_prob_ui       as _scm_quark_prod,
    mckubre_lenr                   as _scm_mckubre,
    s26_3_from_vds                 as _scm_s26_3_from_vds,
    qgp_energy_density_scm         as _scm_qgp_energy_density,
    strange_quark_matter_density   as _scm_sqm_density,
    mit_bag_scm                    as _scm_mit_bag,
    ads_cft_scm_dual               as _scm_ads_cft_dual,
    scm_gw_metric_perturbation     as _scm_gw_metric_pert,
)'''

# ===========================================================================
# SCm block — CP3 (no leading _, bare public names)
# ===========================================================================
SCM_CP3_OLD = '''\
# --- import from canonical source when available [pdf/scm_vacuum_manifold.py] ---
try:
    from scm_vacuum_manifold import (
        E_phonon       as E_PHONON_SCM,
        S26_3          as S26_3,
        Phi_resonance  as PHI_RESONANCE,
        KER_SCm        as KER_SCM,
        scaling_factor as SCALING_SCM,    # exact 630 eV normalizer
        KAPPA_FLOAT    as KAPPA_FLOAT,    # float(KAPPA) = 0.0005
        F_TRZ                       as F_TRZ,
        coleman_guillespie_scm      as coleman_guillespie_scm,
        neutrino_oscillation_prob_lenr as neutrino_oscillation_prob_lenr,
        quark_production_prob_ui    as quark_production_prob_ui,
        mckubre_lenr                as mckubre_lenr,
        s26_3_from_vds              as s26_3_from_vds,
        qgp_energy_density_scm      as qgp_energy_density_scm,
        strange_quark_matter_density as strange_quark_matter_density,
        mit_bag_scm                 as mit_bag_scm,
        ads_cft_scm_dual            as ads_cft_scm_dual,
        scm_gw_metric_perturbation  as scm_gw_metric_perturbation,
    )
except ImportError:
    pass  # fallback values already set above
    F_TRZ = 0.1
    def coleman_guillespie_scm(decay_rate=1.0e6, t_n=-100.0, Gamma=1.0e12): return 0.0
    def neutrino_oscillation_prob_lenr(t_n=-100.0): return 0.0
    def quark_production_prob_ui(t_n=-100.0, Gamma=1.0e12): return 0.0
    def mckubre_lenr(PdD_loading=0.9, volume=1.0e-6, t_n=-100.0): return 0.0
    s26_3_from_vds = lambda: 1.4531e26
    def qgp_energy_density_scm(T_plasma=1.0e11): return 0.0
    def strange_quark_matter_density(): return (1.0e18, 0.0)
    def mit_bag_scm(): return 0.0
    ads_cft_scm_dual = lambda: {}
    def scm_gw_metric_perturbation(f_gw=100.0, r_detector=3.086e22): return 0.0'''

SCM_CP3_NEW = '''\
# --- SCm constants from dpm_vacuum_manifold (consolidated) ---
from dpm_vacuum_manifold import (
    E_phonon       as E_PHONON_SCM,
    S26_3          as S26_3,
    Phi_resonance  as PHI_RESONANCE,
    KER_SCm        as KER_SCM,
    scaling_factor as SCALING_SCM,
    KAPPA_FLOAT    as KAPPA_FLOAT,
    F_TRZ          as F_TRZ,
    coleman_guillespie_scm         as coleman_guillespie_scm,
    neutrino_oscillation_prob_lenr as neutrino_oscillation_prob_lenr,
    quark_production_prob_ui       as quark_production_prob_ui,
    mckubre_lenr                   as mckubre_lenr,
    s26_3_from_vds                 as s26_3_from_vds,
    qgp_energy_density_scm         as qgp_energy_density_scm,
    strange_quark_matter_density   as strange_quark_matter_density,
    mit_bag_scm                    as mit_bag_scm,
    ads_cft_scm_dual               as ads_cft_scm_dual,
    scm_gw_metric_perturbation     as scm_gw_metric_perturbation,
)'''

# ===========================================================================
# SCm block — CP4 (_SCM_ prefix, extra functions)
# ===========================================================================
SCM_CP4_OLD = '''\
# ---------------------------------------------------------------------------
# SCm VACUUM MANIFOLD MODULE (scm_vacuum_manifold.py \u2014 clean 27FEB2026_A.docx thread)
# Provides: SSQ, KAPPA, RHO_VAC_SCM, THZ_PHONON, NEG_TIME_RANGE,
#           compute_F_U_Bi_i_numerical(), vds_numerical(), export_all_to_latex()
# ---------------------------------------------------------------------------
try:
    from scm_vacuum_manifold import (
        SSQ          as _SCM_SSQ,           # [SSq] = 0.57
        KAPPA        as _SCM_KAPPA,         # \u03ba = 5.0 \u00d7 10^{-4} day^{-1}
        RHO_VAC_SCM  as _SCM_RHO_VAC,      # _RHO_VAC_SCM kg/m\xb3 vacuum manifold baseline
        THZ_PHONON   as _SCM_THz,          # 1.25 THz Gaussian phonon activation
        E_phonon     as _SCM_E_PHONON,     # h * f_THz [J]  \u2014 new module-level const
        S26_3        as _SCM_S26_3,        # 1.4531e26 Ramanujan amplification
        Phi_resonance as _SCM_PHI_RES,     # 0.84 on-resonance Gaussian factor
        KER_SCm      as _SCM_KER_SCm,      # E_phonon * S26_3 * Phi_resonance [J]
        scaling_factor as _SCM_SCALING,     # exact 630 eV normalizer
        KAPPA_FLOAT  as _SCM_KAPPA_FLOAT,   # float(KAPPA) = 0.0005
        compute_F_U_Bi_i_numerical as _scm_F_U_Bi_i_num,
        monte_carlo_fubi_i         as _scm_monte_carlo_fubi_i,
        vds_numerical              as _scm_vds_num,
        export_all_to_latex        as _scm_export_latex,
        parkhomov_excess_heat       as _scm_parkhomov,      # Ni-H excess heat [kW]
        pons_fleischmann_excess_heat as _scm_pons_fleischmann, # Pd-D excess heat [kW]
        F_TRZ                       as _SCM_F_TRZ,
        coleman_guillespie_scm      as _scm_coleman_guillespie,
        neutrino_oscillation_prob_lenr as _scm_neutrino_osc,
        quark_production_prob_ui    as _scm_quark_prod,
        mckubre_lenr                as _scm_mckubre,
        s26_3_from_vds              as _scm_s26_3_from_vds,
        qgp_energy_density_scm      as _scm_qgp_energy_density,
        strange_quark_matter_density as _scm_sqm_density,
        mit_bag_scm                 as _scm_mit_bag,
        ads_cft_scm_dual            as _scm_ads_cft_dual,
        scm_gw_metric_perturbation  as _scm_gw_metric_pert,
    )
    _SCM_MANIFOLD_LOADED = True
except ImportError:
    _SCM_MANIFOLD_LOADED = False
    _SCM_SSQ      = 0.57
    _SCM_KAPPA    = 5.0e-4
    _SCM_RHO_VAC  = _RHO_VAC_SCM
    _SCM_THz      = 1.25e12
    _SCM_E_PHONON = 6.62607015e-34 * 1.25e12
    _SCM_S26_3    = 1.4531e26
    _SCM_PHI_RES  = 0.84
    _SCM_KER_SCm  = _SCM_E_PHONON * _SCM_S26_3 * _SCM_PHI_RES
    _SCM_SCALING  = 630 * 1.60217662e-19 / (_SCM_E_PHONON * _SCM_S26_3 * _SCM_PHI_RES)  # exact 630 eV normalizer
    _SCM_KAPPA_FLOAT = 0.0005  # float(KAPPA)
    def _scm_F_U_Bi_i_num(**kw): return 0.0
    def _scm_monte_carlo_fubi_i(n_samples=10000): return 0.0, 0.0, [0.0, 0.0]
    def _scm_vds_num(terms=1000): return 0.0
    def _scm_export_latex(): return {}
    def _scm_parkhomov(N_clusters=2.0e18, t_hours=1.0):
        import math as _m_pk
        _energy_per_cluster_j = 630 * 1.60217662e-19
        return N_clusters * _energy_per_cluster_j * _m_pk.exp(-5e-4 * t_hours * 24) / 1e3
    def _scm_pons_fleischmann(PdD_loading=0.9, volume=1e-6):
        return PdD_loading * volume * _SCM_E_PHONON * _SCM_S26_3 * _SCM_PHI_RES * 0.001 * 1e6 / 1e3
    _SCM_F_TRZ = 0.1
    def _scm_coleman_guillespie(decay_rate=1.0e6, t_n=-100.0, Gamma=1.0e12): return 0.0
    def _scm_neutrino_osc(t_n=-100.0): return 0.0
    def _scm_quark_prod(t_n=-100.0, Gamma=1.0e12): return 0.0
    def _scm_mckubre(PdD_loading=0.9, volume=1.0e-6, t_n=-100.0): return 0.0
    _scm_s26_3_from_vds = lambda: 1.4531e26
    def _scm_qgp_energy_density(T_plasma=1.0e11): return 0.0
    def _scm_sqm_density(): return (1.0e18, 0.0)
    def _scm_mit_bag(): return 0.0
    def _scm_ads_cft_dual(): return {}
    def _scm_gw_metric_pert(f_gw=100.0, r_detector=3.086e22): return 0.0'''

SCM_CP4_NEW = '''\
# ---------------------------------------------------------------------------
# SCm VACUUM MANIFOLD IMPORTS (from dpm_vacuum_manifold, consolidated)
# ---------------------------------------------------------------------------
from dpm_vacuum_manifold import (
    SSQ          as _SCM_SSQ,
    KAPPA        as _SCM_KAPPA,
    RHO_VAC_SCM  as _SCM_RHO_VAC,
    THZ_PHONON   as _SCM_THz,
    E_phonon     as _SCM_E_PHONON,
    S26_3        as _SCM_S26_3,
    Phi_resonance as _SCM_PHI_RES,
    KER_SCm      as _SCM_KER_SCm,
    scaling_factor as _SCM_SCALING,
    KAPPA_FLOAT  as _SCM_KAPPA_FLOAT,
    compute_F_U_Bi_i_numerical as _scm_F_U_Bi_i_num,
    monte_carlo_fubi_i         as _scm_monte_carlo_fubi_i,
    vds_numerical              as _scm_vds_num,
    export_all_to_latex        as _scm_export_latex,
    parkhomov_excess_heat       as _scm_parkhomov,
    pons_fleischmann_excess_heat as _scm_pons_fleischmann,
    F_TRZ                       as _SCM_F_TRZ,
    coleman_guillespie_scm      as _scm_coleman_guillespie,
    neutrino_oscillation_prob_lenr as _scm_neutrino_osc,
    quark_production_prob_ui    as _scm_quark_prod,
    mckubre_lenr                as _scm_mckubre,
    s26_3_from_vds              as _scm_s26_3_from_vds,
    qgp_energy_density_scm      as _scm_qgp_energy_density,
    strange_quark_matter_density as _scm_sqm_density,
    mit_bag_scm                 as _scm_mit_bag,
    ads_cft_scm_dual            as _scm_ads_cft_dual,
    scm_gw_metric_perturbation  as _scm_gw_metric_pert,
)
_SCM_MANIFOLD_LOADED = True'''


def do_replace(content, old, new, label):
    if old in content:
        return content.replace(old, new, 1), True
    else:
        print(f"  WARNING: '{label}' not found — check for text mismatch")
        return content, False


# ===========================================================================
# Process each file
# ===========================================================================

files = {
    'CondensedPhysics.py':  ('cp12', 'CP1'),
    'CondensedPhysics2.py': ('cp12', 'CP2'),
    'CondensedPhysics3.py': ('cp3',  'CP3'),
    'CondensedPhysics4.py': ('cp4',  'CP4'),
}

for filename, (style, label) in files.items():
    path = os.path.join(CWD, filename)
    orig = read_text(path)
    content = orig

    # Backup
    shutil.copy2(path, path + '.bak')

    ok_count = 0

    # 1. Top try/except block
    content, ok = do_replace(content, TOP_OLD, TOP_NEW, f'{label} TOP')
    if ok: ok_count += 1

    # 2. SCm try/except block (style-dependent)
    if style == 'cp12':
        content, ok = do_replace(content, SCM_CP12_OLD, SCM_CP12_NEW, f'{label} SCm-CP12')
    elif style == 'cp3':
        content, ok = do_replace(content, SCM_CP3_OLD, SCM_CP3_NEW, f'{label} SCm-CP3')
    elif style == 'cp4':
        content, ok = do_replace(content, SCM_CP4_OLD, SCM_CP4_NEW, f'{label} SCm-CP4')
    if ok: ok_count += 1

    # 3. UA try/except block
    content, ok = do_replace(content, UA_OLD, UA_NEW, f'{label} UA')
    if ok: ok_count += 1

    # 4. DPM try/except block
    content, ok = do_replace(content, DPM_OLD, DPM_NEW, f'{label} DPM')
    if ok: ok_count += 1

    write_text(path, content)
    print(f"{label}: {ok_count}/4 replacements made — {filename}")

print("DONE.")
