# -*- coding: utf-8 -*-
"""
_chain_trace_26layer.py
26-Layer DPM Amplification -- First-Principles Mass Derivation

Shows:
  1. Each layer's physical contribution (SCm_i, UA_i, B0_i, w_i = i^6)
  2. Cumulative sum toward A_26 = 1,307,798,101
  3. AMU derived from M_0_DPM × A_26 (no PDG lookup for the derivation)
  4. Proton, neutron, C-12, Fe-56 from vacuum constants only
  5. Scale_factor residual per element vs ACP chain
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dpm_vacuum_manifold import (
    M_0_DPM, AMU, C_LIGHT, HBAR, OMEGA_CW, OMEGA_CCW, SSQ, N_LAYERS,
    chain_26layer_weights,
    chain_26layer_amplification,
    chain_derive_nucleon_mass,
    chain_derive_particle_masses,
    chain_acp_M_proto,
    RHO_VAC_SCM,
)

SEP  = "=" * 72
SEP2 = "-" * 72


# ---------------------------------------------------------------------------
# 0. Constants summary
# ---------------------------------------------------------------------------
print(SEP)
print("26-LAYER DPM AMPLIFICATION  --  First-Principles Mass Derivation")
print(SEP)
print(f"\n  M_0_DPM          = {M_0_DPM:.6e} kg   (= rho_SCm / [SSq])")
print(f"  rho_SCm          = {float(RHO_VAC_SCM):.4e} kg/m^3  (vacuum SCm density)")
print(f"  [SSq]            = {float(SSQ):.4f}          (E_crack gate)")
print(f"  AMU (PDG)        = {AMU:.6e} kg")
print(f"  AMU / M_0_DPM    = {AMU/M_0_DPM:.6e}  (target for 26-layer sum)")
print(f"  N_LAYERS         = {N_LAYERS}")
print()

# ---------------------------------------------------------------------------
# 1. Per-layer breakdown
# ---------------------------------------------------------------------------
print(SEP)
print("LAYER-BY-LAYER BREAKDOWN  (w_i = [SCm]_i × [UA]_i × B0_i = i^6)")
print(SEP)
print(f"  {'i':>3}  {'[SCm]_i=i^2':>12}  {'[UA]_i=i':>10}  "
      f"{'B0_i=i^3':>10}  {'w_i=i^6':>14}  cumulative_sum")
print(SEP2)

layers   = chain_26layer_weights()
cum_sum  = 0
for lyr in layers:
    i      = lyr["i"]
    cum_sum += lyr["w_i"]
    m_derived = M_0_DPM * cum_sum
    print(f"  {i:>3}  {lyr['SCm_i']:>12,}  {lyr['UA_i']:>10,}  "
          f"{lyr['B0_i']:>10,}  {lyr['w_i']:>14,}  "
          f"Σ={cum_sum:>14,}  M={m_derived:.4e} kg")

print()

# ---------------------------------------------------------------------------
# 2. Amplification factor summary
# ---------------------------------------------------------------------------
print(SEP)
print("26-LAYER AMPLIFICATION FACTOR  (A_26)")
print(SEP)
amp = chain_26layer_amplification()
print(f"\n  A_26 (exact integer)      = {amp['A_26_exact']:,}")
print(f"  AMU_derived = M_0×A_26    = {amp['AMU_derived_kg']:.6e} kg")
print(f"  AMU_observed (PDG)        = {amp['AMU_observed_kg']:.6e} kg")
print(f"  Error                     = {amp['error_pct']:+.2f}%")
print(f"  f_SSq_gate                = {amp['f_SSq_gate']:.6f}  "
      f"  (≈1: residual is [SSq]=0.57 E_crack gate)")
print()
print(f"  PREDICTED rho_SCm         = {amp['rho_SCm_predicted_kg_m3']:.4e} kg/m^3")
print(f"  CANONICAL  rho_SCm        = {amp['rho_SCm_canonical_kg_m3']:.4e} kg/m^3")
print(f"  Prediction error          = {amp['rho_prediction_error_pct']:+.2f}%")
print()
print(f"  {amp['derivation_note']}")
print()

# ---------------------------------------------------------------------------
# 3. Nucleon mass for A = 1, 12, 56
# ---------------------------------------------------------------------------
print(SEP)
print("NUCLEON MASS DERIVATION  (M = A × M_0_DPM × A_26, no PDG input)")
print(SEP)
for A_val in [1, 12, 56]:
    res = chain_derive_nucleon_mass(A_val)
    print(f"\n  A={A_val:>3}: derived = {res['M_derived_kg']:.6e} kg"
          f"  PDG_ref = {res['M_PDG_ref_kg']:.6e} kg"
          f"  error = {res['error_pct']:+.2f}%")
    print(f"         formula: A × M_0_DPM × A_26  (A_26={res['A_26']:,})")
    print(f"         inputs:  {res['inputs_used']}")
print()

# ---------------------------------------------------------------------------
# 4. Full particle mass derivation
# ---------------------------------------------------------------------------
print(SEP)
print("PARTICLE MASSES FROM VACUUM CONSTANTS  (no PDG mass input)")
print(SEP)
pm = chain_derive_particle_masses()
A_26 = pm["A_26"]

print(f"\n  Shared:  A_26 = {A_26:,}  |  M_0_DPM = {pm['M_0_DPM_kg']:.6e} kg")
print()

p = pm["proton"]
print(f"  PROTON:")
print(f"    Formula:  {p['formula']}")
print(f"    Derived:  {p['derived_kg']:.6e} kg")
print(f"    Observed: {p['observed_kg']:.6e} kg")
print(f"    Error:    {p['error_pct']:+.2f}%")
print()

n_d = pm["neutron"]
print(f"  NEUTRON:")
print(f"    Formula:  {n_d['formula']}")
print(f"    Leading-order derived:  {n_d['derived_kg']:.6e} kg")
print(f"    Observed:               {n_d['observed_kg']:.6e} kg")
print(f"    Error (leading order):  {n_d['error_pct']:+.2f}%")
print(f"    Observed n-p split:     {n_d['delta_M_np_observed']:.4e} kg  (1.293 MeV/c^2)")
print(f"    hbar*Dw/2c^2 estimate:  {n_d['delta_M_np_hbar_est']:.4e} kg  (too small by ~11 orders)")
print(f"    Mechanism: {n_d['mechanism']}")
print()

e = pm["electron"]
print(f"  ELECTRON:")
print(f"    Observed: {e['observed_kg']:.4e} kg")
print(f"    mp/me ratio (obs): {e['mp_me_ratio_obs']:.4f}")
print(f"    Note: {e['note']}")
print()

c12 = pm["carbon_12"]
print(f"  CARBON-12:")
print(f"    Formula:  {c12['formula']}")
print(f"    Derived:  {c12['derived_kg']:.6e} kg")
print(f"    Observed: {c12['observed_kg']:.6e} kg")
print(f"    Error:    {c12['error_pct']:+.2f}%")
print()

fe = pm["iron_56"]
print(f"  IRON-56:")
print(f"    Formula:  {fe['formula']}")
print(f"    Derived:  {fe['derived_kg']:.6e} kg")
print(f"    Observed: {fe['observed_kg']:.6e} kg")
print(f"    Error:    {fe['error_pct']:+.2f}%")
print(f"    Note:     {fe['note']}")
print()

print(f"  SUMMARY: {pm['summary']}")
print()

# ---------------------------------------------------------------------------
# 5. Scale_factor residual vs ACP chain
# ---------------------------------------------------------------------------
print(SEP)
print("ACP CHAIN vs 26-LAYER DERIVATION  (scale_factor audit)")
print(SEP)
print()
print(f"  The current ACP chain gives M_proto(Z) = M_0_DPM × (1-exp(-Z/10)) × Z")
print(f"  The 26-layer derivation gives M_obs(A) = A × M_0_DPM × A_26")
print()
print(f"  {'Z':>3}  {'M_proto (ACP)':>16}  {'M_0×A_26 (26L)':>16}  {'ratio':>12}  {'A_26/ACP':>12}")
print(SEP2)
for Z_check in [1, 2, 6, 8, 20, 26]:
    m_proto = chain_acp_M_proto(Z_check)
    m_26L   = M_0_DPM * A_26
    ratio   = m_26L / m_proto if m_proto > 0 else float('nan')
    acp_val = (1.0 - (2.718281828 ** (-Z_check / 10.0))) * Z_check
    print(f"  {Z_check:>3}  {m_proto:>16.4e}  {m_26L:>16.4e}  {ratio:>12.4e}  {A_26/acp_val:>12.4e}")
print()
print(f"  Note: ACP(Z) = (1-exp(-Z/10)) × Z gives the SEED mass at the chain crossing.")
print(f"  The 26-layer formula M_0×A_26 gives the FINAL nucleon mass (one bundle = 1 AMU).")
print(f"  They are complementary views: ACP traces the resonance ACTIVATION; ")
print(f"  the 26-layer sum gives the AMPLITUDE of each resonance unit.")
print()

print(SEP)
print("DONE -- exit 0")
print(SEP)
