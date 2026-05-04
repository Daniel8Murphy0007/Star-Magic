#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
_fix_remaining_imports.py
Fix stale scm_vacuum_manifold / ua_vacuum_manifold imports in:
  QCalc.py, qcalcgeom_helpers.py, 99system_master_equation.py,
  99system_wstp_gamma.py, _chain_trace_SSq.py
"""
import re, shutil

def read_text(path):
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        return f.read()

def write_text(path, content):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(content)

def do_replace(content, old, new, label):
    if old in content:
        return content.replace(old, new, 1), True
    print(f"  WARNING: '{label}' not found")
    return content, False

def do_regex(content, pattern, replacement, label):
    # replacement may contain backslashes — use a lambda to avoid regex substitution escaping
    if hasattr(pattern, 'pattern'):
        compiled = pattern
    else:
        compiled = re.compile(pattern, re.DOTALL)
    result = compiled.search(content)
    if result:
        new_content = content[:result.start()] + replacement + content[result.end():]
        return new_content, True
    print(f"  WARNING: regex '{label}' not matched")
    return content, False

# ─────────────────────────────────────────────────────────────────────────────
# COMMON BLOCKS
# ─────────────────────────────────────────────────────────────────────────────

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

# CP3-style block (bare names, no underscore prefix) — used by 99system*.py
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

# ─────────────────────────────────────────────────────────────────────────────
# QCalc.py — second block (_QC suffix aliases)
# ─────────────────────────────────────────────────────────────────────────────
SCM_QC_OLD = '''\
# --- import from canonical source when available [pdf/scm_vacuum_manifold.py] ---
try:
    from scm_vacuum_manifold import (
        E_phonon                       as _E_PHONON_SCM_QC,
        S26_3                          as _S26_3_QC,
        Phi_resonance                  as _PHI_RESONANCE_QC,
        KER_SCm                        as _KER_SCM_QC,
        scaling_factor                 as _SCALING_SCM_QC,    # exact 630 eV normalizer
        KAPPA_FLOAT                    as _KAPPA_FLOAT_QC,    # float(KAPPA) = 0.0005
        F_TRZ                          as _F_TRZ_QC,          # canonical 0.1
        coleman_guillespie_scm         as _scm_coleman_guillespie_qc,  # beta->1.25THz phonon
        neutrino_oscillation_prob_lenr as _scm_neutrino_osc_qc,        # P_osc via S26_3*cos(pi*t_n)
        quark_production_prob_ui       as _scm_quark_prod_qc,          # |Phi|^2*cos(pi*t_n)*Ui
        mckubre_lenr                   as _scm_mckubre_qc,             # Pd-D sub-barrier via F_U_Bi_i
        s26_3_from_vds                 as _scm_s26_3_from_vds_qc,      # canonical 1.4531e26
        qgp_energy_density_scm         as _scm_qgp_energy_density_qc,  # QGP tokamak energy [J]
        strange_quark_matter_density   as _scm_sqm_density_qc,         # SQM (density, B_eff)
        mit_bag_scm                    as _scm_mit_bag_qc,             # MIT bag B_eff [J/m^3]
        ads_cft_scm_dual               as _scm_ads_cft_dual_qc,        # SCm holographic mapping dict
        scm_gw_metric_perturbation     as _scm_gw_metric_pert_qc,      # h(f_gw, r) GW strain
    )
except ImportError:
    pass  # fallback values already set above
    _F_TRZ_QC = 0.1
    def _scm_coleman_guillespie_qc(decay_rate=1.0e6, t_n=-100.0, Gamma=1.0e12): return 0.0
    def _scm_neutrino_osc_qc(t_n=-100.0): return 0.0
    def _scm_quark_prod_qc(t_n=-100.0, Gamma=1.0e12): return 0.0
    def _scm_mckubre_qc(PdD_loading=0.9, volume=1.0e-6, t_n=-100.0): return 0.0
    _scm_s26_3_from_vds_qc = lambda: 1.4531e26
    def _scm_qgp_energy_density_qc(T_plasma=1.0e11): return 0.0
    def _scm_sqm_density_qc(): return (1.0e18, 0.0)
    def _scm_mit_bag_qc(): return 0.0
    def _scm_ads_cft_dual_qc(): return {}
    def _scm_gw_metric_pert_qc(f_gw=100.0, r_detector=3.086e22): return 0.0'''

SCM_QC_NEW = '''\
# --- SCm constants from dpm_vacuum_manifold (consolidated, _QC aliases) ---
from dpm_vacuum_manifold import (
    E_phonon                       as _E_PHONON_SCM_QC,
    S26_3                          as _S26_3_QC,
    Phi_resonance                  as _PHI_RESONANCE_QC,
    KER_SCm                        as _KER_SCM_QC,
    scaling_factor                 as _SCALING_SCM_QC,
    KAPPA_FLOAT                    as _KAPPA_FLOAT_QC,
    F_TRZ                          as _F_TRZ_QC,
    coleman_guillespie_scm         as _scm_coleman_guillespie_qc,
    neutrino_oscillation_prob_lenr as _scm_neutrino_osc_qc,
    quark_production_prob_ui       as _scm_quark_prod_qc,
    mckubre_lenr                   as _scm_mckubre_qc,
    s26_3_from_vds                 as _scm_s26_3_from_vds_qc,
    qgp_energy_density_scm         as _scm_qgp_energy_density_qc,
    strange_quark_matter_density   as _scm_sqm_density_qc,
    mit_bag_scm                    as _scm_mit_bag_qc,
    ads_cft_scm_dual               as _scm_ads_cft_dual_qc,
    scm_gw_metric_perturbation     as _scm_gw_metric_pert_qc,
)'''

# ─────────────────────────────────────────────────────────────────────────────
# qcalcgeom_helpers.py — _SCM_ prefix with _SCM_LOADED flag (CP4 style)
# Uses regex because header comment may have mojibake
# ─────────────────────────────────────────────────────────────────────────────
QCALCGEOM_PAT = re.compile(
    r'# \u2500+ SCm Vacuum Manifold module \u2500+\n'
    r'# Source: scm_vacuum_manifold\.py.*?'
    r'def _scm_gw_metric_pert\(f_gw=100\.0, r_detector=3\.086e22\): return 0\.0',
    re.DOTALL
)
QCALCGEOM_NEW = '''\
# ── SCm Vacuum Manifold (from dpm_vacuum_manifold, consolidated) ────────────
from dpm_vacuum_manifold import (
    SSQ          as _SCM_SSQ,
    KAPPA        as _SCM_KAPPA,
    RHO_VAC_SCM  as _SCM_RHO_VAC,
    THZ_PHONON   as _SCM_THz,
    compute_F_U_Bi_i_numerical as _scm_F_U_Bi_i_num,
    vds_numerical              as _scm_vds_num,
    E_phonon                   as _SCM_E_PHONON,
    S26_3                      as _SCM_S26_3,
    Phi_resonance              as _SCM_PHI_RES,
    KER_SCm                    as _SCM_KER_SCm,
    scaling_factor             as _SCM_SCALING,
    KAPPA_FLOAT                as _SCM_KAPPA_FLOAT,
    F_TRZ                      as _SCM_F_TRZ,
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
)
_SCM_LOADED = True'''

# Same pattern (without header variation) for QCalc.py top SCm block
QCALC_SCM_TOP_PAT = re.compile(
    r'# \u2500+ SCm Vacuum Manifold module \u2500+\n'
    r'# Source: scm_vacuum_manifold\.py.*?'
    r'def _scm_gw_metric_pert\(f_gw=100\.0, r_detector=3\.086e22\): return 0\.0',
    re.DOTALL
)

# ─────────────────────────────────────────────────────────────────────────────
# PROCESS FILES
# ─────────────────────────────────────────────────────────────────────────────

results = {}

# ── QCalc.py ────────────────────────────────────────────────────────────────
path = 'QCalc.py'
shutil.copy2(path, path + '.bak')
c = read_text(path)
ok = 0
c, b = do_replace(c, TOP_OLD, TOP_NEW, 'QCalc TOP'); ok += b
c, b = do_regex(c, QCALC_SCM_TOP_PAT, QCALCGEOM_NEW, 'QCalc SCm-top-block'); ok += b
c, b = do_replace(c, SCM_QC_OLD, SCM_QC_NEW, 'QCalc SCm-QC'); ok += b
write_text(path, c)
results[path] = ok
print(f"QCalc.py: {ok}/3 replacements")

# ── qcalcgeom_helpers.py ────────────────────────────────────────────────────
path = 'qcalcgeom_helpers.py'
shutil.copy2(path, path + '.bak')
c = read_text(path)
ok = 0
c, b = do_regex(c, QCALCGEOM_PAT, QCALCGEOM_NEW, 'qcalcgeom SCm block'); ok += b
write_text(path, c)
results[path] = ok
print(f"qcalcgeom_helpers.py: {ok}/1 replacements")

# ── 99system_master_equation.py ─────────────────────────────────────────────
path = '99system_master_equation.py'
shutil.copy2(path, path + '.bak')
c = read_text(path)
ok = 0
c, b = do_replace(c, TOP_OLD, TOP_NEW, '99master TOP'); ok += b
c, b = do_replace(c, SCM_CP3_OLD, SCM_CP3_NEW, '99master SCm'); ok += b
write_text(path, c)
results[path] = ok
print(f"99system_master_equation.py: {ok}/2 replacements")

# ── 99system_wstp_gamma.py ──────────────────────────────────────────────────
path = '99system_wstp_gamma.py'
shutil.copy2(path, path + '.bak')
c = read_text(path)
ok = 0
c, b = do_replace(c, TOP_OLD, TOP_NEW, '99wstp TOP'); ok += b
c, b = do_replace(c, SCM_CP3_OLD, SCM_CP3_NEW, '99wstp SCm'); ok += b
write_text(path, c)
results[path] = ok
print(f"99system_wstp_gamma.py: {ok}/2 replacements")

# ── _chain_trace_SSq.py ─────────────────────────────────────────────────────
path = '_chain_trace_SSq.py'
shutil.copy2(path, path + '.bak')
c = read_text(path)
old = 'from scm_vacuum_manifold import SSQ'
new = 'from dpm_vacuum_manifold import SSQ'
c, b = do_replace(c, old, new, '_chain_trace SSQ'); ok = b
write_text(path, c)
results[path] = ok
print(f"_chain_trace_SSq.py: {int(ok)}/1 replacements")

print("\nDONE.")
