#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fix CP4 SCm import block using a start/end marker approach."""
import re, shutil

path = 'CondensedPhysics4.py'
shutil.copy2(path, path + '.bak2')

with open(path, 'r', encoding='utf-8') as f:
    content = f.read()

# Match the CP4-specific section: starts at the header comment block
# and ends after the last fallback function in the except clause.
# Use regex with DOTALL so . matches newlines.
pattern = re.compile(
    r'# ---------------------------------------------------------------------------\n'
    r'# SCm VACUUM MANIFOLD MODULE \(scm_vacuum_manifold\.py.*?'   # header (with mojibake — any chars)
    r'def _scm_gw_metric_pert\(f_gw=100\.0, r_detector=3\.086e22\): return 0\.0',
    re.DOTALL
)

replacement = '''\
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

new_content, count = pattern.subn(replacement, content, count=1)
if count == 1:
    with open(path, 'w', encoding='utf-8') as f:
        f.write(new_content)
    print(f"CP4 SCm block replaced ({count}).")
else:
    print("ERROR: pattern not found in CP4. No changes made.")
